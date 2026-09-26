#!/usr/bin/env python3
"""
Render a run's forecast maps and its HTML viewer.

    python src/make_maps.py --run-dir DATA/tensors_3d/obs_<stamp>
    python src/make_maps.py --run-dir ... --forecast forecast_w15.npz --out maps_w15

Writes <run-dir>/maps/<product>_fHHH.png, <run-dir>/maps/index.html, and
<run-dir>/maps/data/ (hHHH.js: every map field per hour, for the hover
readout; sHHH.js: every model column per hour, for the click-for-sounding
panel). The data files are JavaScript so the page works from file://.
Open index.html in a browser (copy the maps/ folder, or tunnel to the
server). Every product is drawn for every whole hour the forecast reached;
hour 0 is the analysis. Analysis-with-observations maps come from the run's
own archived observations, and error maps from its verification archive
entry, once `verify` has run.

Nothing is fetched except, once, the Natural Earth state and coast lines
(cached under <data>/static; see src/maps/geography.py).
"""

import argparse
import gzip
import json
import sys
import time
from datetime import datetime
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for p in (ROOT, ROOT / "src", ROOT / "src" / "analysis", ROOT / "src" / "dynamics",
          ROOT / "src" / "verification"):
    sys.path.insert(0, str(p))

import config                                                     # noqa: E402
from maps import derive, geography, render, viewer                # noqa: E402
import matplotlib.pyplot as plt                                   # noqa: E402


def run_time_of(run_dir, f=None):
    rt = str(f.get("run_time", "")) if f is not None else ""
    if rt:
        try:
            return datetime.fromisoformat(rt[:19])
        except ValueError:
            pass
    stamp = "_".join(Path(run_dir).name.split("_")[-2:])
    return datetime.strptime(stamp, "%Y%m%d_%H")


def load_raw(raw_dir):
    raw = {}
    for path in sorted(Path(raw_dir).glob("*.gz")):
        name = path.name[:-3]
        if name.startswith("mrms"):
            continue
        with gzip.open(path, "rb") as fh:
            raw[name] = fh.read().decode("utf-8", errors="replace")
    return raw


def station_reports(run_dir, cycle, domain):
    """Surface reports and 850/500 hPa sounding values that passed QC."""
    raw_dir = Path(run_dir) / "observations"
    if not raw_dir.exists():
        return [], {}, set()
    import sources
    from observations import run_qc
    obs = sources.reparse(load_raw(raw_dir), cycle, sources.analysis_box(domain))
    kept, _, _ = run_qc(obs)
    box = geography.domain_box(domain, margin=0.3)
    inside = lambda o: box[0] <= o.lon <= box[1] and box[2] <= o.lat <= box[3]
    sfc, ua = {}, {850: {}, 500: {}}
    for o in kept:
        if not inside(o):
            continue
        if o.pressure is None:
            r = sfc.setdefault(o.station, {"station": o.station, "lat": o.lat, "lon": o.lon})
            if o.variable == "TMP":
                r["tf"] = (o.value - 273.15) * 9 / 5 + 32
            elif o.variable == "UGRD":
                r["u"] = o.value
            elif o.variable == "VGRD":
                r["v"] = o.value
            elif o.variable == "PMSL":
                r["pmsl"] = o.value / 100.0
        else:
            for L in ua:
                if abs(o.pressure - L * 100.0) < 50.0:
                    r = ua[L].setdefault(o.station, {"station": o.station, "lat": o.lat, "lon": o.lon})
                    if o.variable == "TMP":
                        r["tc"] = o.value - 273.15
                    elif o.variable == "UGRD":
                        r["u"] = o.value
                    elif o.variable == "VGRD":
                        r["v"] = o.value
    withheld = set()
    try:
        prov = json.loads(str(np.load(Path(run_dir) / "obs_analysis_f00.npz",
                                      allow_pickle=True)["provenance"]))
        withheld = set(prov.get("withheld", []))
    except Exception:                                   # noqa: BLE001
        pass
    return list(sfc.values()), {L: list(v.values()) for L, v in ua.items()}, withheld


def verification_pairs(run_dir):
    """{lead_hour: {"t": [...], "w": [...]}} from the run's verification archive."""
    vj = Path(run_dir) / "verified.json"
    if not vj.exists():
        return {}
    mpath = Path(json.load(open(vj))["archive"]) / "matches.jsonl"
    if not mpath.exists():
        return {}
    rows = [json.loads(l) for l in open(mpath) if l.strip()]
    out, wind = {}, {}
    for r in rows:
        if r.get("pressure") not in (None, "") or r.get("lead_hours") is None:
            continue
        h = int(round(float(r["lead_hours"])))
        if r["variable"] == "TMP":
            out.setdefault(h, {"t": [], "w": []})["t"].append(
                {"lat": r["lat"], "lon": r["lon"], "err": r["error"]})
        elif r["variable"] in ("UGRD", "VGRD"):
            wind.setdefault((h, r["station"], r["valid_time"]), {}).update(
                {r["variable"]: (r["forecast"], r["observation"]), "lat": r["lat"], "lon": r["lon"]})
    for (h, _, _), w in wind.items():
        if "UGRD" in w and "VGRD" in w:
            fs = np.hypot(w["UGRD"][0], w["VGRD"][0]); os_ = np.hypot(w["UGRD"][1], w["VGRD"][1])
            out.setdefault(h, {"t": [], "w": []})["w"].append(
                {"lat": w["lat"], "lon": w["lon"], "err": (fs - os_) * render.MS_TO_KT})
    return out


# ---------------------------------------------------------------------------
# Data for the viewer: hover values and soundings
# ---------------------------------------------------------------------------

def pack(channels):
    """
    [(name, 2-D array), ...] -> JSON-able dict of int16 base64 with a scale
    and offset per channel; NaN -> -32768. Quantisation error is below
    1/65000 of each channel's range.
    """
    import base64
    names, scale, offset, q = [], [], [], []
    for name, a in channels:
        a = np.asarray(a, dtype=float)
        ok = np.isfinite(a)
        lo = float(a[ok].min()) if ok.any() else 0.0
        hi = float(a[ok].max()) if ok.any() else 1.0
        sc = (hi - lo) / 65000.0 or 1e-6
        qi = np.where(ok, np.round((np.where(ok, a, lo) - lo) / sc) - 32500, -32768)
        names.append(name); scale.append(sc); offset.append(lo + 32500 * sc)
        q.append(qi.astype("<i2"))
    blob = base64.b64encode(np.stack(q).tobytes()).decode("ascii")
    return {"names": names, "scale": scale, "offset": offset, "b64": blob}


def write_js(path, var, hour, obj):
    with open(path, "w", encoding="utf-8") as fh:
        fh.write(f"window.{var} = window.{var} || {{}};\n"
                 f"window.{var}[{hour}] = {json.dumps(obj)};\n")


def hover_channels(d, dt1=None):
    """Every map field at one hour, in display units (see HOVER_META)."""
    kt = render.MS_TO_KT
    tf = lambda K: (K - 273.15) * 9 / 5 + 32
    ch = [("terrain", d["terrain"]), ("mslp", d["mslp"] / 100.0),
          ("t_low_f", tf(d["T_low"])),
          ("wind_low_u", d["u_low"] * kt), ("wind_low_v", d["v_low"] * kt),
          ("omega700", -d["omega700"] * 10.0), ("avort500", d["avort500"] * 1e5)]
    if "thick_1000_500" in d:
        ch.append(("thick", d["thick_1000_500"] / 10.0))
    if dt1 is not None:
        ch.append(("dt1_f", dt1 * 9 / 5))
    if "T_2m" in d:
        ch.append(("t2m_f", tf(d["T_2m"])))
    for L in (925, 850, 700, 500, 300, 250):
        below = d[f"below{L}"]
        ch += [(f"t{L}", np.where(below, np.nan, d[f"T{L}"] - 273.15)),
               (f"z{L}", np.where(below, np.nan, d[f"Z{L}"] / 10.0)),
               (f"wind{L}_u", np.where(below, np.nan, d[f"u{L}"] * kt)),
               (f"wind{L}_v", np.where(below, np.nan, d[f"v{L}"] * kt))]
    return ch


HOVER_META = {   # name: (label, unit, decimals); winds are shown as dir/speed
    "terrain": ("Terrain", "m", 0), "mslp": ("MSLP", "mb", 1), "thick": ("1000-500 thickness", "dam", 0),
    "t_low_f": ("Temperature (lowest level)", "\u00b0F", 0), "dt1_f": ("1-hr change", "\u00b0F", 1),
    "t2m_f": ("2 m temperature", "\u00b0F", 0), "wind_low": ("Wind (lowest level)", "kt", 0),
    "omega700": ("700 mb vertical velocity", "-\u00b5b/s", 1), "avort500": ("500 mb abs. vorticity", "10\u207b\u2075 s\u207b\u00b9", 1),
    **{f"t{L}": (f"{L} mb temperature", "\u00b0C", 1) for L in (925, 850, 700, 500, 300, 250)},
    **{f"z{L}": (f"{L} mb height", "dam", 0) for L in (925, 850, 700, 500, 300, 250)},
    **{f"wind{L}": (f"{L} mb wind", "kt", 0) for L in (925, 850, 700, 500, 300, 250)},
}


SOUNDING_STRIDE = 2      # soundings on every 2nd grid point (24 km): ~0.5 MB an hour


def sounding_channels_sigma(f, i, stride=SOUNDING_STRIDE):
    p, T, u, v, Z, ps = derive.column_for_sounding(f, i)
    s = (slice(None, None, stride), slice(None, None, stride))
    ch = [("ps", ps[s])]
    for k in range(T.shape[0]):
        ch += [(f"T{k}", T[k][s]), (f"u{k}", u[k][s]), (f"v{k}", v[k][s])]
    return ch


def sounding_channels_plev(an, stride=SOUNDING_STRIDE):
    feats = np.asarray(an["features"], dtype=float)
    c = an["channels"]
    s = (slice(None, None, stride), slice(None, None, stride))
    ch = []
    for k in range(feats.shape[1]):
        for name in ("TMP", "UGRD", "VGRD", "RH", "HGT"):
            ch.append((f"{name}{k}", feats[c.index(name), k][s]))
    return ch


def match_hours(hours, every=1, tol=0.1):
    """
    {whole hour: snapshot index}. Snapshots land a few seconds after each
    output time (the model takes whole steps: 1.0023 h, 2.0045 h, ... at
    dt = 17.1 s), so each hour takes the nearest snapshot within tol hours.
    Requiring exact hours dropped every forecast hour (prompt 123).
    """
    hrs = np.asarray(hours, dtype=float)
    out = {}
    if not len(hrs):
        return out
    for H in range(every, int(np.floor(hrs.max() + tol)) + 1, every):
        j = int(np.argmin(np.abs(hrs - H)))
        if abs(hrs[j] - H) < tol:
            out[H] = j
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--forecast", default="forecast.npz",
                    help="forecast file inside the run dir (default forecast.npz)")
    ap.add_argument("--out", default="maps", help="output folder inside the run dir")
    ap.add_argument("--every", type=int, default=1, help="hours between maps (default 1)")
    ap.add_argument("--products", default="", help="comma list to render only some products")
    ap.add_argument("--no-data", action="store_true",
                    help="skip the hover and sounding data files")
    a = ap.parse_args()

    t0 = time.time()
    run_dir = Path(a.run_dir)
    out = run_dir / a.out
    (out / "data").mkdir(parents=True, exist_ok=True)
    domain = config.DOMAIN
    only = {s_.strip() for s_ in a.products.split(",") if s_.strip()}
    want = lambda k: not only or k in only

    fpath = run_dir / a.forecast
    f = derive.load_forecast(fpath, domain) if fpath.exists() else None
    apath = run_dir / "obs_analysis_f00.npz"
    an = derive.load_analysis(apath) if apath.exists() else None
    if f is None and an is None:
        print(f"nothing to map in {run_dir}")
        return 1
    from geo import cell_centres, spacing_m
    if f is not None:
        lat, lon, dx, dy = f["lat"], f["lon"], f["dx"], f["dy"]
        terrain = np.asarray(f["terrain"], float)
    else:
        ny, nx = an["features"].shape[-2:]
        lat, lon = cell_centres(domain, ny, nx)
        dy, dx = spacing_m(domain, ny, nx)
        terrain = np.asarray(np.load(run_dir / "terrain.npz")["terrain"], float)
    cycle = run_time_of(run_dir, f)
    m = render.Mapper(lat, lon, geography.load_boundaries(domain))

    stopped = str(f.get("stopped", "")) if f is not None else ""
    status = ("forecast " + stopped) if stopped else ("forecast complete" if f is not None else "analysis only")
    first_guess = ""
    av = run_dir / "availability.json"
    if av.exists():
        first_guess = json.load(open(av)).get("first_guess", "")
    left = f"first guess: {first_guess}" if first_guess else ""

    # Hour -> source. Snapshots land a few seconds after each output time
    # (a whole number of model steps), so match the nearest within 0.1 h.
    snaps = {}
    if an is not None:
        snaps[0] = ("a", None)
    z_low = None
    if f is not None:
        for H, j in match_hours(f["hours"], a.every).items():
            snaps[H] = ("f", j)
        z_low = derive.snapshot(f, 0)["z_low_agl"]

    n_img = 0
    t_low = {}
    for hour in sorted(snaps):
        kind, i = snaps[hour]
        d = (derive.analysis_snapshot(an, terrain, dx, dy, lat, z_low) if kind == "a"
             else derive.snapshot(f, i))
        t_low[hour] = d["T_low"]
        prev = hour - a.every
        dt1 = t_low[hour] - t_low[prev] if prev in t_low and a.every == 1 else None
        if dt1 is not None:
            d["dT1"] = dt1
        sub = render.time_labels(cycle, hour)
        right = "hour 0 = the analysis" if kind == "a" else (
            f"run stopped: {stopped}" if stopped else "")
        for key, fn in render.FORECAST_PRODUCTS.items():
            if not want(key) or (key == "tchg" and dt1 is None):
                continue
            fig = fn(m, d, sub, (left, right))
            fig.savefig(out / f"{key}_f{hour:03d}.png", dpi=render.DPI); plt.close(fig); n_img += 1
        if kind == "a":
            sfc, ua, withheld = station_reports(run_dir, cycle, domain)
            note = (left, "analysis and the reports it was built from")
            figs = []
            if want("an_sfc"):
                figs.append(("an_sfc", render.product_an_sfc(m, d, sub, note, sfc, withheld)))
            if want("an_mslp"):
                figs.append(("an_mslp", render.product_an_mslp(m, d, sub, note, sfc, withheld)))
            for L, tlv, hs in ((850, np.arange(-30, 32, 1), 3), (500, np.arange(-45, 2, 1), 6)):
                if want(f"an_{L}"):
                    figs.append((f"an_{L}", render._an_level(m, d, sub, note, L, ua.get(L, []), tlv, hs)))
            for key, fig in figs:
                fig.savefig(out / f"{key}_f000.png", dpi=render.DPI); plt.close(fig); n_img += 1
        if not a.no_data:
            write_js(out / "data" / f"h{hour:03d}.js", "NWPH", hour, pack(hover_channels(d, dt1)))
            snd = ({"stride": SOUNDING_STRIDE, "nx": int(np.ceil(lat.shape[1] / SOUNDING_STRIDE)),
                    "kind": "plev", "levels_hPa": [float(x) for x in an["levels_hPa"]],
                    **pack(sounding_channels_plev(an))} if kind == "a" else
                   {"stride": SOUNDING_STRIDE, "nx": int(np.ceil(lat.shape[1] / SOUNDING_STRIDE)),
                    "kind": "sigma", "sigma": [float(x) for x in f["sigma"]],
                    "p_top": float(f["p_top"]), **pack(sounding_channels_sigma(f, i))})
            write_js(out / "data" / f"s{hour:03d}.js", "NWPS", hour, snd)

    ver = verification_pairs(run_dir)
    for hour, pr in sorted(ver.items()):
        if hour not in snaps:
            continue
        d = {"terrain": terrain}
        sub = render.time_labels(cycle, hour)
        note = (left, "surface stations; forecast at the station's height")
        for key, k in (("err_t", "t"), ("err_w", "w")):
            if want(key) and pr[k]:
                fig = render.product_error(m, d, sub, note, pr[k], k)
                fig.savefig(out / f"{key}_f{hour:03d}.png", dpi=render.DPI); plt.close(fig); n_img += 1

    # The viewer lists every image in the folder, not only this invocation's,
    # so a later partial render (the error maps after `verify`) adds to it.
    products = {}
    for key, (g, lab, hover) in render.PRODUCTS.items():
        hrs_ = sorted(int(q.stem.rsplit("_f", 1)[1]) for q in out.glob(f"{key}_f[0-9][0-9][0-9].png"))
        if hrs_:
            products[key] = {"group": g, "label": lab, "hours": hrs_, "hover": hover}
    data_hours = sorted(int(q.stem[1:]) for q in (out / "data").glob("h[0-9][0-9][0-9].js"))
    grid = {"ny": int(lat.shape[0]), "nx": int(lat.shape[1]),
            "lat0": float(lat[0, 0]), "lon0": float(lon[0, 0]),
            "dlat": float(lat[1, 0] - lat[0, 0]), "dlon": float(lon[0, 1] - lon[0, 0])}
    viewer.write_viewer(out / "index.html", f"{run_dir.name}  (init {cycle:%Y-%m-%d %H}Z)",
                        products, status,
                        "Hover for values; click for the model sounding. Near-surface fields are "
                        "the lowest model level. The model is dry: no precipitation or moisture "
                        "products, and forecast soundings have no dewpoint.",
                        extra={"frame": m.pixel_frame(), "grid": grid, "hover_meta": HOVER_META,
                               "data_hours": data_hours, "init": cycle.isoformat()})
    print(f"maps: {n_img} images drawn, {len(products)} products in the viewer, hours "
          f"{min(snaps)}-{max(snaps)} ({len(snaps)} times) -> {out / 'index.html'}  "
          f"({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

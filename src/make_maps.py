#!/usr/bin/env python3
"""
Render a run's forecast maps and its HTML viewer.

    python src/make_maps.py --run-dir DATA/tensors_3d/obs_<stamp>
    python src/make_maps.py --run-dir ... --forecast forecast_w15.npz --out maps_w15

Writes <run-dir>/maps/<product>_fHHH.png and <run-dir>/maps/index.html.
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


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--run-dir", required=True)
    ap.add_argument("--forecast", default="forecast.npz",
                    help="forecast file inside the run dir (default forecast.npz)")
    ap.add_argument("--out", default="maps", help="output folder inside the run dir")
    ap.add_argument("--every", type=float, default=1.0,
                    help="hours between maps (default 1; snapshots must exist)")
    ap.add_argument("--products", default="",
                    help="comma list to render only some products")
    a = ap.parse_args()

    t0 = time.time()
    run_dir = Path(a.run_dir)
    out = run_dir / a.out
    out.mkdir(parents=True, exist_ok=True)
    domain = config.DOMAIN
    only = {s.strip() for s in a.products.split(",") if s.strip()}
    want = lambda k: not only or k in only

    fpath = run_dir / a.forecast
    f = derive.load_forecast(fpath, domain) if fpath.exists() else None
    apath = run_dir / "obs_analysis_f00.npz"
    an = derive.load_analysis(apath) if apath.exists() else None
    if f is None and an is None:
        print(f"nothing to map in {run_dir}")
        return 1
    if f is not None:
        lat, lon, dx, dy, terrain = f["lat"], f["lon"], f["dx"], f["dy"], np.asarray(f["terrain"], float)
    else:
        from geo import cell_centres, spacing_m
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

    snaps = {}                                      # hour -> diagnostics
    z_low = None
    if f is not None:
        for i, h in enumerate(f["hours"]):
            if abs(h / a.every - round(h / a.every)) < 1e-6:
                snaps[int(round(h))] = ("f", i)
        z_low = derive.snapshot(f, 0)["z_low_agl"]
    if an is not None:
        snaps[0] = ("a", None)

    products = {}
    def add(key, hour, fig):
        fig.savefig(out / f"{key}_f{hour:03d}.png", dpi=render.DPI)
        plt.close(fig)
        g, lab = render.PRODUCTS[key]
        products.setdefault(key, {"group": g, "label": lab, "hours": []})["hours"].append(hour)

    n_img = 0
    for hour in sorted(snaps):
        kind, i = snaps[hour]
        d = (derive.analysis_snapshot(an, terrain, dx, dy, lat, z_low) if kind == "a"
             else derive.snapshot(f, i))
        sub = render.time_labels(cycle, hour)
        right = "hour 0 = the analysis" if kind == "a" else (
            f"run stopped: {stopped}" if stopped else "")
        for key, fn in render.FORECAST_PRODUCTS.items():
            if want(key):
                add(key, hour, fn(m, d, sub, (left, right))); n_img += 1
        if kind == "a":
            sfc, ua, withheld = station_reports(run_dir, cycle, domain)
            note = (left, "analysis and the reports it was built from")
            if want("an_sfc"):
                add("an_sfc", 0, render.product_an_sfc(m, d, sub, note, sfc, withheld)); n_img += 1
            if want("an_mslp"):
                add("an_mslp", 0, render.product_an_mslp(m, d, sub, note, sfc, withheld)); n_img += 1
            for L, tlv, hs in ((850, np.arange(-30, 32, 1), 3), (500, np.arange(-45, 2, 1), 6)):
                if want(f"an_{L}"):
                    add(f"an_{L}", 0, render._an_level(m, d, sub, note, L, ua.get(L, []), tlv, hs)); n_img += 1

    ver = verification_pairs(run_dir)
    for hour, pr in sorted(ver.items()):
        if hour not in snaps:
            continue
        kind, i = snaps[hour]
        d = {"terrain": terrain}
        sub = render.time_labels(cycle, hour)
        note = (left, "surface stations; forecast at the station's height")
        for key, k in (("err_t", "t"), ("err_w", "w")):
            if want(key) and pr[k]:
                add(key, hour, render.product_error(m, d, sub, note, pr[k], k)); n_img += 1

    # The viewer lists every image in the folder, not only this invocation's,
    # so a later partial render (the error maps after `verify`) adds to it.
    products = {}
    for key, (g, lab) in render.PRODUCTS.items():
        hrs = sorted(int(q.stem.rsplit("_f", 1)[1]) for q in out.glob(f"{key}_f[0-9][0-9][0-9].png"))
        if hrs:
            products[key] = {"group": g, "label": lab, "hours": hrs}
    viewer.write_viewer(out / "index.html", f"{run_dir.name}  (init {cycle:%Y-%m-%d %H}Z)",
                        products, status,
                        "Near-surface fields are the lowest model level. The model is dry: "
                        "no precipitation or moisture products. Dashed outline: the model domain.")
    print(f"maps: {n_img} images drawn, {len(products)} products in the viewer, hours "
          f"{min(snaps)}-{max(snaps)} -> {out / 'index.html'}  ({time.time() - t0:.0f} s)")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

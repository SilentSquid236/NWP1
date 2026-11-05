"""
GFS lateral-boundary frames for an observation-only run (prompt 266, test GE).

The INITIAL STATE stays our own analysis (obs_analysis_f00.npz, observations
only). GFS supplies nothing but the frames the edge relaxation reads at hours
1..H, written as gfs_edge_f01..fHH.npz next to the analysis. Until now those
edges held the analysis fixed for 24 h (P-53, frozen edges).

Source: NOAA GFS 0.5-degree pressure-level files on the AWS open-data bucket
(noaa-gfs-bdp-pds). Only the 90 messages needed (TMP, RH, UGRD, VGRD, HGT on 18
levels) are fetched by byte range from the .idx inventory, about 20 MB per
forecast hour. The GFS cycle used is the one LAG hours before ours (default 6):
at our cycle's start time (+1:45) the GFS of the same hour is not yet out.

Processing, per frame:
  * horizontal: bilinear from the 0.5-degree lat/lon grid to our cell centres;
  * vertical: our 875 and 825 hPa levels (not in GFS) linear in ln p between
    the neighbours;
  * time: GFS 0.5-degree is 3-hourly; hourly frames are linear in time.

Three modes (RH clipped 0-100 in all):
  abs   frame(h) = GFS(t0 + h)
  tend  frame(h) = analysis + [GFS(t0 + h) - GFS(t0)]
        The edges start exactly at our analysis and take only the GFS CHANGE,
        so there is no jump at hour 0 and GFS's own state never enters.
  ramp  frame(h) = tend + min(h / RAMP_H, 1) [GFS(t0) - analysis]
        Starts as tend, becomes abs after RAMP_H hours (default 6).

Usage (production, tools/daily.sh with NWP_GFS_EDGES=<mode>):
  python tools/gfs_edges.py --run-dir data/tensors_3d/obs_20261028_00 \
         --out-dir data/tensors_3d/obs_20261028_00/gfs_edges --mode abs
  python src/forecast.py --run-dir ... --edge-frames .../gfs_edges ...
Only a run given --edge-frames reads the frames. The dry run that becomes the
next cycle's first guess is not given them, so no GFS value reaches an analysis.
--copy-analysis also copies obs_analysis_f00.npz and terrain.npz to --out-dir
(test GE ran each arm from such a copy).
"""
import argparse
import os
import re
import shutil
import sys
import time
import urllib.request
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent))
import config  # noqa: E402

BUCKET = "https://noaa-gfs-bdp-pds.s3.amazonaws.com"
VARS = {"TMP": "TMP", "RH": "RH", "UGRD": "UGRD", "VGRD": "VGRD", "HGT": "HGT"}
# our level -> GFS levels it comes from (two entries: ln-p interpolation)
MISSING = {875: (900, 850), 825: (850, 800)}


def gfs_levels():
    need = set()
    for p in config.PRESSURE_LEVELS:
        need.update(MISSING.get(p, (p,)))
    return sorted(need, reverse=True)


def url_for(cycle, fxx):
    d, h = cycle.strftime("%Y%m%d"), cycle.strftime("%H")
    return f"{BUCKET}/gfs.{d}/{h}/atmos/gfs.t{h}z.pgrb2.0p50.f{fxx:03d}"


def fetch(url, rng=None, tries=4):
    for k in range(tries):
        try:
            req = urllib.request.Request(url)
            if rng:
                req.add_header("Range", f"bytes={rng[0]}-{'' if rng[1] is None else rng[1]}")
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:  # network hiccup: retry, then give up loudly
            if k == tries - 1:
                raise RuntimeError(f"{url} {rng}: {e}")
            time.sleep(5 * (k + 1))


def subset_grib(cycle, fxx, cache):
    """The needed messages of one GFS file, cached as a small GRIB2 file."""
    out = cache / f"gfs.{cycle:%Y%m%d}" / f"{cycle:%H}" / f"sub.0p50.f{fxx:03d}.grib2"
    if out.exists() and out.stat().st_size > 0:
        return out
    url = url_for(cycle, fxx)
    idx = fetch(url + ".idx").decode().strip().splitlines()
    rows = []
    for ln in idx:
        f = ln.split(":")
        rows.append((int(f[1]), f[3], f[4]))
    levs = {f"{p} mb" for p in gfs_levels()}
    want = []
    for i, (off, var, lev) in enumerate(rows):
        if var in VARS and lev in levs:
            end = rows[i + 1][0] - 1 if i + 1 < len(rows) else None
            want.append((off, end))
    n_exp = len(VARS) * len(levs)
    if len(want) != n_exp:
        raise RuntimeError(f"{url}: {len(want)} messages matched, expected {n_exp}")
    # merge contiguous byte ranges
    merged = []
    for a, b in sorted(want):
        if merged and merged[-1][1] is not None and merged[-1][1] + 1 == a:
            merged[-1] = (merged[-1][0], b)
        else:
            merged.append((a, b))
    out.parent.mkdir(parents=True, exist_ok=True)
    tmp = out.with_suffix(".part")
    with open(tmp, "wb") as fh:
        for rg in merged:
            fh.write(fetch(url, rg))
    os.replace(tmp, out)
    return out


def decode(path):
    """{(var, level_hPa): 2-D field}, plus the 1-D lat (N->S) and lon (0..360)."""
    import eccodes
    fields, lat, lon = {}, None, None
    with open(path, "rb") as fh:
        while True:
            g = eccodes.codes_grib_new_from_file(fh)
            if g is None:
                break
            try:
                name = eccodes.codes_get(g, "shortName")
                lev = int(eccodes.codes_get(g, "level"))
                ni, nj = eccodes.codes_get(g, "Ni"), eccodes.codes_get(g, "Nj")
                vals = eccodes.codes_get_values(g).reshape(nj, ni)
                if lat is None:
                    la0, la1 = eccodes.codes_get(g, "latitudeOfFirstGridPointInDegrees"), \
                        eccodes.codes_get(g, "latitudeOfLastGridPointInDegrees")
                    lo0, lo1 = eccodes.codes_get(g, "longitudeOfFirstGridPointInDegrees"), \
                        eccodes.codes_get(g, "longitudeOfLastGridPointInDegrees")
                    lat, lon = np.linspace(la0, la1, nj), np.linspace(lo0, lo1, ni)
            finally:
                eccodes.codes_release(g)
            key = {"t": "TMP", "r": "RH", "u": "UGRD", "v": "VGRD", "gh": "HGT"}.get(name)
            if key:
                fields[(key, lev)] = vals
    return fields, lat, lon


def bilinear_weights(lat1d, lon1d, la, lo):
    """Weights to sample a regular lat/lon grid at points (la, lo)."""
    lo = np.mod(lo, 360.0)
    if lat1d[0] > lat1d[-1]:                  # N->S: work in flipped index
        fy = (lat1d[0] - la) / (lat1d[0] - lat1d[1])
    else:
        fy = (la - lat1d[0]) / (lat1d[1] - lat1d[0])
    fx = (lo - lon1d[0]) / (lon1d[1] - lon1d[0])
    j0, i0 = np.floor(fy).astype(int), np.floor(fx).astype(int)
    wy, wx = fy - j0, fx - i0
    return j0, i0, wy, wx


def sample(field, w):
    j0, i0, wy, wx = w
    f00, f01 = field[j0, i0], field[j0, i0 + 1]
    f10, f11 = field[j0 + 1, i0], field[j0 + 1, i0 + 1]
    return (1 - wy) * ((1 - wx) * f00 + wx * f01) + wy * ((1 - wx) * f10 + wx * f11)


def gfs_frame(cycle, fxx, cache, la, lo):
    """One GFS time as a [C, L, Y, X] array on our grid and levels."""
    flds, lat1d, lon1d = decode(subset_grib(cycle, fxx, cache))
    w = bilinear_weights(lat1d, lon1d, la, lo)
    out = np.empty((len(config.CHANNELS), len(config.PRESSURE_LEVELS)) + la.shape)
    for c, ch in enumerate(config.CHANNELS):
        for k, p in enumerate(config.PRESSURE_LEVELS):
            if p in MISSING:
                pa, pb = MISSING[p]
                a = sample(flds[(ch, pa)], w)
                b = sample(flds[(ch, pb)], w)
                t = (np.log(p) - np.log(pa)) / (np.log(pb) - np.log(pa))
                out[c, k] = (1 - t) * a + t * b
            else:
                out[c, k] = sample(flds[(ch, p)], w)
    return out


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__.split("\n\n")[0])
    p.add_argument("--run-dir", required=True, help="dir holding obs_analysis_f00.npz")
    p.add_argument("--out-dir", default=None)
    p.add_argument("--hours", type=int, default=24)
    p.add_argument("--mode", choices=("abs", "tend", "ramp"), default="tend")
    p.add_argument("--ramp-h", type=float, default=6.0)
    p.add_argument("--lag", type=int, default=6, help="GFS cycle = ours minus LAG hours")
    p.add_argument("--step", type=int, default=3, help="GFS output interval (0.5-degree: 3 h)")
    p.add_argument("--cache", default=None, help="GRIB subset cache (default data/gfs_cache)")
    p.add_argument("--keep-days", type=float, default=3.0,
                   help="delete cached GFS cycles older than this many days (0 = keep all)")
    p.add_argument("--copy-analysis", action="store_true")
    a = p.parse_args(argv)

    run = Path(a.run_dir)
    out = Path(a.out_dir) if a.out_dir else run
    cache = Path(a.cache) if a.cache else Path(os.environ.get("NWP_DATA_ROOT", HERE.parent / "data")) / "gfs_cache"
    z = np.load(run / "obs_analysis_f00.npz", allow_pickle=False)
    A = z["features"].astype(float)
    meta = {k: z[k] for k in z.files if k != "features"}
    assert list(meta["channels"]) == list(config.CHANNELS)
    assert list(meta["levels_hPa"]) == list(config.PRESSURE_LEVELS)
    la, lo = meta["lat"].astype(float), meta["lon"].astype(float)
    t0 = datetime.strptime(str(meta["run_time"])[:13], "%Y-%m-%dT%H")
    gcy = t0 - timedelta(hours=a.lag)
    print(f"GFS edges: analysis {t0:%Y-%m-%d %HZ}, GFS cycle {gcy:%Y-%m-%d %HZ} "
          f"(f{a.lag:03d}..f{a.lag + a.hours:03d}, every {a.step} h), mode {a.mode}")

    # the GFS times bracketing hours 0..H, fetched once each
    need = sorted({a.step * ((a.lag + h) // a.step) for h in range(a.hours + 1)}
                  | {a.step * -(-(a.lag + h) // a.step) for h in range(a.hours + 1)})
    t_fetch = time.time()
    G = {f: gfs_frame(gcy, f, cache, la, lo) for f in need}
    print(f"  fetched/decoded {len(need)} GFS times in {time.time() - t_fetch:.0f} s")

    def at(fx):
        lo_, hi_ = a.step * (fx // a.step), a.step * -(-fx // a.step)
        if lo_ == hi_:
            return G[lo_]
        w = (fx - lo_) / (hi_ - lo_)
        return (1 - w) * G[lo_] + w * G[hi_]

    out.mkdir(parents=True, exist_ok=True)
    if a.copy_analysis and out.resolve() != run.resolve():
        for f in ("obs_analysis_f00.npz", "terrain.npz"):
            if (run / f).exists():
                shutil.copy2(run / f, out / f)
    for f in out.glob("gfs_edge_f*.npz"):
        f.unlink()
    g0 = at(a.lag)
    ir = config.CHANNELS.index("RH")
    d0 = g0 - A
    print("  GFS(t0) minus analysis, domain mean / RMS, by channel at 850 hPa: "
          + ", ".join(f"{ch} {d0[c, 6].mean():+.2f}/{np.sqrt((d0[c, 6] ** 2).mean()):.2f}"
                      for c, ch in enumerate(config.CHANNELS)))
    for h in range(1, a.hours + 1):
        g = at(a.lag + h)
        if a.mode == "abs":
            f = g.copy()
        elif a.mode == "tend":
            f = A + (g - g0)
        else:
            f = A + (g - g0) + min(h / a.ramp_h, 1.0) * (g0 - A)
        f[ir] = np.clip(f[ir], 0.0, 100.0)
        m = dict(meta)
        m["valid_time"] = np.array((t0 + timedelta(hours=h)).strftime("%Y-%m-%dT%H:%M:%S"))
        m["source"] = np.array(f"GFS 0.5deg {gcy:%Y%m%d%H} f{a.lag + h:03d} edges ({a.mode}{'' if a.mode != 'ramp' else f' {a.ramp_h:g} h'})")
        m.pop("provenance", None); m.pop("first_guess", None); m.pop("p_msl", None)
        np.savez_compressed(out / f"gfs_edge_f{h:02d}.npz", features=f.astype(np.float32), **m)
    print(f"  wrote gfs_edge_f01..f{a.hours:02d} to {out}")
    if a.keep_days > 0:
        cut = time.time() - a.keep_days * 86400.0
        old = [d for d in cache.glob("gfs.*/*") if d.is_dir() and d.stat().st_mtime < cut]
        for d in old:
            shutil.rmtree(d, ignore_errors=True)
        for d in cache.glob("gfs.*"):
            if d.is_dir() and not any(d.iterdir()):
                d.rmdir()
        if old:
            print(f"  cache: removed {len(old)} GFS cycles older than {a.keep_days:g} days")
    return 0


if __name__ == "__main__":
    sys.exit(main())

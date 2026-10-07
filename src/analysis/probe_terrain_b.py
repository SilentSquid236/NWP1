"""
P-56 test B: is it the TERRAIN FIELD? Same observations, HRRR's terrain.

    python src/analysis/probe_terrain_b.py <obs run dir> <YYYY-mm-ddTHH>

DIAGNOSTIC ONLY. HRRR's surface height is used here as a second terrain
field -- static geography, not weather -- to answer one question, and the
output directory is named so it can never be mistaken for a run. Nothing
here is archived or verified.

Rebuilds the analysis from the run's ARCHIVED observations (no new fetch
except HRRR's two surface fields), on HRRR terrain block-averaged onto the
same 12 km cells ETOPO was, and writes <run dir>_testB_hrrrterrain/ ready for
`python src/forecast.py --run-dir ... --hours 24`. It also prints both
terrain fields' height and slope statistics, since a steeper or rougher
field is the most direct way terrain could differ.

Prediction (research log, 2026-09-22): if the terrain field is the cause, the
rebuilt run passes the hour its ETOPO twin died at.
"""

import sys
from datetime import datetime
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent
for p in (HERE, SRC, SRC.parent, SRC / "verification", SRC / "dynamics"):
    sys.path.insert(0, str(p))

import numpy as np

import config
import geo
import build
import sources
import ingest_obs
from observations import run_qc


def hrrr_terrain_on_grid(when, ny, nx):
    import ingest_hrrr as IH
    from herbie import Herbie
    H = Herbie(when.strftime("%Y-%m-%d %H:%M"), model="hrrr", product="prs",
               fxx=0, verbose=False, save_dir=IH.herbie_save_dir())
    ds = IH._as_single_dataset(IH._open_hrrr(H, IH.HRRR_SURFACE_SEARCH,
                                             allow_full=False, verbose=True))
    ysl, xsl = IH.domain_slice(ds)
    z = np.asarray(ds[IH.resolve_surface(ds, "HGT")].values)[ysl, xsl]
    lat = np.asarray(ds.latitude.values)[ysl, xsl]
    lon = np.asarray(ds.longitude.values)[ysl, xsl]
    lon = np.where(lon > 180, lon - 360, lon)
    return np.clip(geo.block_average(lat.ravel(), lon.ravel(), z.ravel(),
                                     config.DOMAIN, ny, nx), 0.0, None)


def slope_stats(t, dy, dx):
    gy, gx = np.gradient(t, dy, dx)
    s = np.hypot(gx, gy)
    return f"max {t.max():6.0f} m, mean {t.mean():5.0f} m, slope max {s.max():.4f}, p99 {np.percentile(s, 99):.4f}"


def main():
    run = Path(sys.argv[1])
    cycle = datetime.strptime(sys.argv[2], "%Y-%m-%dT%H")
    etopo = np.load(run / "terrain.npz")["terrain"].astype(float)
    ny, nx = etopo.shape
    dy, dx = geo.spacing_m(config.DOMAIN, ny, nx)
    hrrr = hrrr_terrain_on_grid(cycle, ny, nx)
    print(f"  ETOPO : {slope_stats(etopo, dy, dx)}")
    print(f"  HRRR  : {slope_stats(hrrr, dy, dx)}")
    d = hrrr - etopo
    print(f"  HRRR - ETOPO: mean {d.mean():+.0f} m, RMS {np.sqrt((d**2).mean()):.0f} m, "
          f"range {d.min():+.0f} .. {d.max():+.0f} m")

    raw = ingest_obs.load_raw(run / "observations")
    obs = sources.reparse(raw, cycle, sources.analysis_box(config.DOMAIN))
    kept, _, _ = run_qc(obs)
    feats, meta = build.build_analysis(cycle, kept, hrrr, config.DOMAIN,
                                       config.PRESSURE_LEVELS)
    meta["source"] = "observations (TEST B: HRRR terrain, diagnostic only)"
    out = run.parent / f"{run.name}_testB_hrrrterrain"
    out.mkdir(exist_ok=True)
    np.savez_compressed(out / "obs_analysis_f00.npz", features=feats, **meta)
    np.savez_compressed(out / "terrain.npz", terrain=hrrr.astype(np.float32))
    print(f"\nwrote {out}\nnext: python -u src/forecast.py --run-dir {out} --hours 24")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

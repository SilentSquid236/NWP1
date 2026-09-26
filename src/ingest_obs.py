"""
Build one cycle's initial state from observations.

    python src/ingest_obs.py --cycle 2026-09-21T12
    python src/ingest_obs.py --cycle 2026-09-21T12 --skip ndbc      # what does a source add?
    python src/ingest_obs.py --cycle 2026-09-21T12 --from-raw       # rebuild, no network

Writes, into <data>/tensors_3d/obs_<YYYYmmdd_HH>/:

    observations/*.gz       every payload VERBATIM, written first
    obs_analysis_f00.npz    the analysis, in the ingest format forecast.py reads
    terrain.npz             ETOPO terrain on the same grid
    availability.json       which sources answered, QC, first guess, and the
                            analysis scored against the WITHHELD stations

WHAT GOES IN (2026-09-22, prompts 94-101)

Every reliable observation at or before the cycle time -- never after it. A
missing source is skipped and logged. Where soundings are missing the upper
air comes from the previous run's 6-hour forecast. No model output of anyone
else's enters at any point.
"""

import argparse
import gzip
import json
import sys
import time
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "analysis"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "verification"))

import resources
RESOURCE_PLAN = resources.apply()

import numpy as np

import config
import sources
import geo
import build
from observations import run_qc

CYCLE_HOURS = (0, 6, 12, 18)


def run_dir_for(cycle):
    return config.TENSOR_DIR / f"obs_{cycle:%Y%m%d_%H}"


def store_raw(raw_dir, raw):
    """Verbatim, compressed, BEFORE anything is parsed (see src/verify.py)."""
    raw_dir.mkdir(parents=True, exist_ok=True)
    for name, payload in raw.items():
        path = raw_dir / f"{name}.gz"
        data = payload if isinstance(payload, bytes) else payload.encode("utf-8")
        with gzip.open(path, "wb") as f:
            f.write(data)


def load_raw(raw_dir):
    raw = {}
    for path in sorted(Path(raw_dir).glob("*.gz")):
        name = path.name[:-3]
        if name.startswith("mrms"):
            continue
        with gzip.open(path, "rb") as f:
            raw[name] = f.read().decode("utf-8", errors="replace")
    return raw


def previous_run(cycle):
    """
    What the run 6 h earlier left: (forecast spec, its RH, its analysis, note).

    The forecast is offered only if the file exists; whether it reaches +6 h
    is decided by build.first_guess_from_forecast, which refuses a lead time
    it does not have. The analysis is the fallback when it does not.
    """
    prev = cycle - timedelta(hours=6)
    d = run_dir_for(prev)
    fc = d / "forecast.npz"
    a = d / "obs_analysis_f00.npz"
    rh = analysis = None
    if a.exists():
        z = np.load(a, allow_pickle=False)
        feats = z["features"].astype(float)
        rh = feats[build.CHANNELS.index("RH")]
        analysis = (str(a), feats)
    spec = (str(fc), 6.0) if fc.exists() else None
    note = (str(fc) if spec else f"no forecast in {d}") + \
        ("" if analysis else "; no analysis either")
    return spec, rh, analysis, note


def score_withheld(feats, terrain, obs, domain):
    """The analysis at the withheld ASOS stations: its only honest score."""
    Tsfc = build.value_at_height(feats[0].astype(float), feats[4].astype(float),
                                 terrain + 2.0, "TMP")
    err = []
    for o in obs:
        if (o.variable != "TMP" or not build.is_withheld(o.station, o.source)
                or o.elevation is None or not sources.in_box(o.lat, o.lon, domain)):
            continue
        tg = float(geo.bilinear(Tsfc, o.lat, o.lon, domain))
        zg = float(geo.bilinear(terrain, o.lat, o.lon, domain))
        if np.isfinite(tg):
            err.append(tg - build.LAPSE * (o.elevation - zg) - o.value)
    err = np.asarray(err)
    if err.size == 0:
        return {"n": 0}
    return {"n": int(err.size), "rmse_K": float(np.sqrt(np.mean(err ** 2))),
            "bias_K": float(err.mean())}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--cycle", required=True, help="e.g. 2026-09-21T12 (00/06/12/18Z)")
    p.add_argument("--skip", default="", help="comma-separated sources to leave out")
    p.add_argument("--from-raw", action="store_true",
                   help="rebuild from this cycle's archived payloads; no network")
    p.add_argument("--spacing-km", type=float, default=12.0)
    p.add_argument("--no-previous", action="store_true",
                   help="ignore the previous run (cold start)")
    p.add_argument("--max-slope", type=float, default=0.0086,
                   help="smooth terrain to at most this slope (P-56); "
                        "0 disables")
    args = p.parse_args()

    cycle = None
    for fmt in ("%Y-%m-%dT%H", "%Y-%m-%d %H", "%Y-%m-%dT%H:%M"):
        try:
            cycle = datetime.strptime(args.cycle, fmt)
            break
        except ValueError:
            pass
    if cycle is None:
        p.error(f"could not parse --cycle {args.cycle!r}")
    if cycle.hour not in CYCLE_HOURS or cycle.minute:
        p.error(f"cycles are 00, 06, 12, 18Z; got {cycle:%H:%M}Z")
    now = datetime.now(timezone.utc).replace(tzinfo=None)
    if cycle > now:
        p.error(f"{cycle:%Y-%m-%d %H}Z is in the future")

    t0 = time.time()
    out = run_dir_for(cycle)
    out.mkdir(parents=True, exist_ok=True)
    domain = config.DOMAIN
    box = sources.analysis_box(domain)
    print("Observation analysis")
    print(config.describe())
    print(resources.describe(RESOURCE_PLAN))
    print(f"  cycle          : {cycle:%Y-%m-%d %H}Z  (observations at or before it only)")
    print(f"  output         : {out}\n")

    raw_dir = out / "observations"
    if args.from_raw:
        raw = load_raw(raw_dir)
        if not raw:
            print(f"--from-raw, but nothing archived in {raw_dir}")
            return 1
        obs = sources.reparse(raw, cycle, box)
        table = [{"name": "archive", "status": "ok", "n": len(obs)}]
        print(f"  rebuilt from {len(raw)} archived payloads: {len(obs)} obs")
    else:
        skip = {s.strip() for s in args.skip.split(",") if s.strip()}
        results = sources.collect(cycle, domain, skip=skip)
        raw = {}
        for r in results:
            raw.update(r.raw)
        store_raw(raw_dir, raw)                         # BEFORE any parsing
        obs = sources.all_observations(results)
        table = [{"name": r.name, "status": r.status, "n": len(r.obs),
                  "platforms": r.platforms, "seconds": round(r.seconds, 1),
                  "message": r.message} for r in results]

    kept, rejected, qc = run_qc(obs)
    print(f"\n  QC             : {qc}")

    ny, nx = geo.grid_shape(domain, args.spacing_km * 1000.0)
    terrain, tsrc = geo.load_terrain(domain, ny, nx, config.DATA_ROOT / "static")
    raw_max = terrain.max()
    if args.max_slope > 0:
        terrain, n_pass, slope = geo.limit_slope(terrain, domain, args.max_slope)
        smooth_msg = (f"; slope limited to {slope:.4f} in {n_pass} passes, "
                      f"peak {raw_max:.0f} -> {terrain.max():.0f} m")
    else:
        smooth_msg = "; NOT slope-limited (--max-slope 0)"
    print(f"  terrain        : {terrain.min():.0f}-{terrain.max():.0f} m on "
          f"{ny}x{nx} ({tsrc}){smooth_msg}")

    prev, prev_rh, prev_an, prev_msg = (None, None, None, "cold start requested") \
        if args.no_previous else previous_run(cycle)
    print(f"  previous run   : {prev_msg}")

    try:
        feats, meta = build.build_analysis(cycle, kept, terrain, domain,
                                           config.PRESSURE_LEVELS,
                                           previous_forecast=prev,
                                           previous_rh=prev_rh,
                                           previous_analysis=prev_an)
    except ValueError as e:
        print(f"\n  ANALYSIS FAILED: {e}")
        json.dump({"cycle": cycle.isoformat(), "sources": table, "qc": qc,
                   "error": str(e)}, open(out / "availability.json", "w"), indent=2)
        return 1

    np.savez_compressed(out / "obs_analysis_f00.npz", features=feats, **meta)
    np.savez_compressed(out / "terrain.npz", terrain=terrain.astype(np.float32))

    wh = score_withheld(feats, terrain, kept, domain)
    prov = json.loads(meta["provenance"])
    json.dump({"cycle": cycle.isoformat(), "sources": table, "qc": qc,
               "first_guess": meta["first_guess"], "withheld_score": wh,
               "soundings_per_level": prov["soundings_per_level"],
               "surface_obs_used": prov["surface_obs_used"],
               "seconds": round(time.time() - t0, 1)},
              open(out / "availability.json", "w"), indent=2)

    print(f"  first guess    : {meta['first_guess']}")
    print(f"  withheld score : {wh}")
    print(f"\nWrote obs_analysis_f00.npz {feats.shape} in {time.time()-t0:.0f} s -> {out}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

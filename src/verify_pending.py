"""
Verify every archived forecast whose window has closed, once.

    python src/verify_pending.py            # verify what is due
    python src/verify_pending.py --dry      # list what is due, touch nothing

Prompt 101 (2026-09-22): after the forecast time is over, the forecast is
checked against surface observations. A run cannot do that itself -- the
observations it would be scored against do not exist yet -- so this is a
separate job with its own cron entry, outside the 1.5 h run budget.

A run is DUE when its last output hour plus `--latency-h` has passed. It is
verified at most once: a `verified.json` beside its forecast records the
result, and a run that raised is retried next time rather than marked done.
Nothing in data/ is overwritten; `verify.verify` appends to the archive and
skips pairs it already holds.
"""

import argparse
import os
import json
import sys
import traceback
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import numpy as np

import config
from verify import verify


def run_time_of(run_dir):
    """obs_20260921_12 / analysis_20260921_12 -> datetime."""
    stamp = "_".join(Path(run_dir).name.split("_")[-2:])
    return datetime.strptime(stamp, "%Y%m%d_%H")


def pending(root, now, latency_h, max_age_h=None, product="forecast.npz", marker="verified.json"):
    """Runs under root whose window has closed and that carry no `marker` yet.

    product (2026-10-20, NWP_12KM_PRODUCT): the file to verify in each run dir,
    e.g. forecast_nhm.npz for the moist 12 km product; a run dir without it
    falls back to forecast.npz. Returns (due, waiting) of (dir, run time,
    hours, ready time, forecast path).
    """
    due, waiting = [], []
    for fc0 in sorted(Path(root).glob("*_????????_??/forecast.npz")):
        d = fc0.parent
        fc = d / product if (d / product).exists() else fc0
        if (d / marker).exists():
            continue
        try:
            rt = run_time_of(d)
        except ValueError:
            continue
        if max_age_h is not None and now - rt > timedelta(hours=max_age_h):
            continue
        z = np.load(fc, allow_pickle=False)
        last = float(np.asarray(z["times_s"]).max()) / 3600.0
        ready = rt + timedelta(hours=last + latency_h)
        (due if now >= ready else waiting).append((d, rt, last, ready, fc))
    return due, waiting


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--latency-h", type=float, default=2.0,
                   help="hours after the window closes before observations "
                        "are complete in the archive (default 2)")
    p.add_argument("--dry", action="store_true")
    p.add_argument("--max-age-3km-h", type=float, default=72.0,
                   help="verify 3 km runs (tensors_3km) only if their run time is within "
                        "this many hours (default 72: the September campaign reruns, "
                        "verified by their own jobs, are left alone)")
    args = p.parse_args()

    now = datetime.now(timezone.utc).replace(tzinfo=None)
    # The 12 km runs, and (P-75, CAM stage S7) the 3 km runs of the last
    # --max-age-3km-h hours. The 3 km pairs go to their own archive: the
    # archive de-duplicates by station and time, so one archive for both grids
    # would keep only the first grid's pairs. Their forecasts are not copied.
    # NWP_12KM_PRODUCT (2026-10-20): the 12 km file that is verified as THE 12 km
    # forecast (default forecast.npz). When it is the moist run, the dry run is also
    # verified, for runs of the last --max-age-3km-h hours, into verification_12dry
    # (marker verified_dry.json), so the two can be compared over the overlap.
    product = os.environ.get("NWP_12KM_PRODUCT", "forecast.npz")
    sets = [(config.TENSOR_DIR, config.DATA_ROOT / "verification", None, True, product, "verified.json"),
            (config.DATA_ROOT / "tensors_3km", config.DATA_ROOT / "verification_3km",
             args.max_age_3km_h, False, "forecast.npz", "verified.json")]
    if product != "forecast.npz":
        sets.append((config.TENSOR_DIR, config.DATA_ROOT / "verification_12dry",
                     args.max_age_3km_h, False, "forecast.npz", "verified_dry.json"))
    status = 0
    print(f"Deferred verification, {now:%Y-%m-%d %H:%M}Z")
    for root, archive, max_age, copy, prod, marker in sets:
        if not Path(root).exists():
            continue
        due, waiting = pending(root, now, args.latency_h, max_age, prod, marker)
        print(f"  due {len(due)}, waiting {len(waiting)}  (runs under {root})")
        for d, rt, last, ready, _ in waiting:
            print(f"    waiting  {d.name}  +{last:.0f} h, due {ready:%Y-%m-%d %H:%M}Z")
        if args.dry:
            for d, rt, last, _, fc in due:
                print(f"    due      {d.name}  +{last:.0f} h  {fc.name}")
            continue
        for d, rt, last, _, fc in due:
            print(f"\n  {d.name}: run {rt:%Y-%m-%d %H}Z, {last:.0f} h, {fc.name}")
            try:
                matches, paths = verify(fc, archive, run_time=rt,
                                        copy_forecast=copy)
            except Exception as e:                          # noqa: BLE001
                traceback.print_exc()
                print(f"    FAILED ({type(e).__name__}); will retry next run")
                status = 1
                continue
            errs = {}
            for m in matches:
                errs.setdefault(m["variable"], []).append(m["error"])
            summary = {v: {"n": len(e), "bias": float(np.mean(e)),
                           "rmse": float(np.sqrt(np.mean(np.square(e))))}
                       for v, e in errs.items()}
            json.dump({"verified_at": now.isoformat(), "run_time": rt.isoformat(),
                       "hours": last, "n_matches": len(matches), "by_variable": summary,
                       "archive": str(paths["dir"]), "forecast": fc.name},
                      open(d / marker, "w"), indent=2)
            for v, s in summary.items():
                print(f"    {v:5s} n={s['n']:5d} bias {s['bias']:+.2f} rmse {s['rmse']:.2f}")
    return status

if __name__ == "__main__":
    raise SystemExit(main())

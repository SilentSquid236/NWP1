"""
Fetch HRRR pressure-level data, subset to the configured domain, and write
[C, L, Y, X] arrays as compressed .npz.

HRRR supplies INITIAL and BOUNDARY conditions for the physics core. It is not
training data and it is never used as verification truth -- see
docs/DATA_ASSIMILATION.md.

Two modes:

  analysis  (default) -- successive hourly F00 analyses. Best estimate of the
                         real atmosphere at each hour; use these to drive the
                         lateral boundaries of a limited-area run.

  forecast            -- one model run's F00..FNN sequence, for comparing our
                         forecast against HRRR's own at matching lead times.

Examples:
    python src/ingest_hrrr.py --start 2026-08-01 --hours 24
    python src/ingest_hrrr.py --start 2026-08-01T00 --hours 6 --dry-run
    python src/ingest_hrrr.py --start 2026-08-01T12 --mode forecast --hours 18
"""

import argparse
import sys
import time
import traceback
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import resources
RESOURCE_PLAN = resources.apply()

import numpy as np
import config
from netpolicy import estimate_ingest_mb, max_mbps


def _open_hrrr(H, search):
    """
    Herbie renamed searchString -> search around 2024. Support both so this
    doesn't break on whichever version the server happens to carry.
    """
    try:
        return H.xarray(search)
    except TypeError:
        return H.xarray(searchString=search)


def _as_single_dataset(ds):
    """cfgrib sometimes returns a list of datasets; merge them."""
    if isinstance(ds, (list, tuple)):
        import xarray as xr
        if len(ds) == 1:
            return ds[0]
        return xr.merge(ds, compat="override", combine_attrs="override")
    return ds


def domain_slice(ds):
    """
    Bounding-box indices for the configured domain.

    HRRR is on a Lambert conformal grid, so latitude and longitude are 2D
    fields, not axes -- you cannot .sel() on them. We build a mask and take
    the enclosing rectangle so the result stays a regular array.
    """
    d = config.DOMAIN
    lat = ds.latitude.values
    lon = ds.longitude.values
    lon = np.where(lon > 180, lon - 360, lon)      # 0..360 -> -180..180

    mask = ((lat >= d["lat_min"]) & (lat <= d["lat_max"]) &
            (lon >= d["lon_min"]) & (lon <= d["lon_max"]))
    if not mask.any():
        raise ValueError(
            f"Domain {d} does not intersect the HRRR grid. "
            f"grid lat {lat.min():.1f}..{lat.max():.1f}, "
            f"lon {lon.min():.1f}..{lon.max():.1f}")

    ys, xs = np.where(mask)
    return slice(ys.min(), ys.max() + 1), slice(xs.min(), xs.max() + 1)


def extract_state(ds, ysl, xsl):
    """Build the [C, L, Y, X] tensor in config.CHANNELS / PRESSURE_LEVELS order."""
    levels = config.PRESSURE_LEVELS

    have = set(int(v) for v in np.atleast_1d(ds.isobaricInhPa.values))
    missing = [lv for lv in levels if lv not in have]
    if missing:
        raise ValueError(f"HRRR is missing requested levels: {missing}. "
                         f"Available: {sorted(have, reverse=True)}")

    # Explicit .sel puts levels in OUR order rather than trusting file order.
    sub = ds.sel(isobaricInhPa=levels)

    planes = []
    for name in config.CHANNELS:
        if name not in sub:
            raise KeyError(f"Variable {name} not in downloaded fields: "
                           f"{list(sub.data_vars)}")
        arr = sub[name].values[:, ysl, xsl]      # [L, Y, X]
        planes.append(arr)

    return np.ascontiguousarray(np.stack(planes, axis=0), dtype=np.float32)


def fetch_hour(when, fxx=0, verbose=True):
    """Download one HRRR field set and return (tensor, metadata)."""
    from herbie import Herbie

    # One regex for all five variables across all isobaric levels -- Herbie
    # byte-range downloads only the matching GRIB messages, so we never pull
    # the full ~130 MB file.
    search = r"^(?:TMP|RH|UGRD|VGRD|HGT):\d+ mb:"

    H = Herbie(when.strftime("%Y-%m-%d %H:%M"), model="hrrr",
               product="prs", fxx=fxx, verbose=False)

    ds = _as_single_dataset(_open_hrrr(H, search))
    ysl, xsl = domain_slice(ds)
    state = extract_state(ds, ysl, xsl)

    lat = ds.latitude.values[ysl, xsl]
    lon = ds.longitude.values[ysl, xsl]
    lon = np.where(lon > 180, lon - 360, lon)

    meta = {
        "valid_time": (when + timedelta(hours=fxx)).isoformat(),
        "run_time": when.isoformat(),
        "fxx": fxx,
        "channels": np.array(config.CHANNELS),
        "levels_hPa": np.array(config.PRESSURE_LEVELS),
        "lat": np.ascontiguousarray(lat, dtype=np.float32),
        "lon": np.ascontiguousarray(lon, dtype=np.float32),
        "source": f"HRRR prs f{fxx:02d}",
    }
    if verbose:
        print(f"    shape {state.shape}  {state.nbytes / 1e6:.1f} MB  "
              f"T[0] {state[0, 0].mean():.1f}K")
    return state, meta


def write_state(path: Path, state, meta):
    """
    Save as compressed .npz -- no torch dependency anywhere in the physics
    pipeline. The whole chain (ingest -> dynamics -> verification) is numpy.
    """
    path.parent.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, hrrr_features=state, **meta)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--start", required=True,
                   help="First valid time, e.g. 2026-08-01 or 2026-08-01T12")
    p.add_argument("--hours", type=int, default=24,
                   help="How many hourly steps to fetch (default 24)")
    p.add_argument("--mode", choices=["analysis", "forecast"], default="analysis")
    p.add_argument("--run-name", default=None,
                   help="Output subfolder name (default: derived from --start)")
    p.add_argument("--overwrite", action="store_true")
    p.add_argument("--dry-run", action="store_true",
                   help="Fetch one hour, report shape, write nothing.")
    p.add_argument("--pause", type=float, default=2.0,
                   help="Seconds to pause between hourly downloads (default 2). "
                        "Herbie manages its own transfers, so pacing between "
                        "files is how we stay polite on a shared link.")
    args = p.parse_args()

    for fmt in ("%Y-%m-%dT%H", "%Y-%m-%d %H", "%Y-%m-%d"):
        try:
            start = datetime.strptime(args.start, fmt)
            break
        except ValueError:
            continue
    else:
        p.error(f"Could not parse --start {args.start!r}")

    run_name = args.run_name or f"{args.mode}_{start:%Y%m%d_%H}"
    out_dir = config.TENSOR_DIR / run_name

    print("HRRR ingestion")
    print(config.describe())
    print(f"  mode           : {args.mode}")
    print(f"  window         : {start:%Y-%m-%d %H}Z + {args.hours}h")
    print(f"  output         : {out_dir}")
    est = estimate_ingest_mb(args.hours, config.N_LEVELS, config.N_CHANNELS)
    print(f"  est. download  : ~{est['total_download_MB']:.0f} MB "
          f"({est['per_hour_download_MB']:.0f} MB/hour), "
          f"pausing {args.pause:.1f}s between files")
    print(f"  bandwidth      : shared link -- cap {max_mbps():.0f} MB/s "
          f"(NWP_MAX_MBPS to override)\n")

    ok = failed = skipped = 0
    for i in range(args.hours):
        if args.mode == "analysis":
            when, fxx = start + timedelta(hours=i), 0
        else:
            when, fxx = start, i

        label = f"f{i:02d}"
        out_path = out_dir / f"live_hrrr_{label}.npz"

        if out_path.exists() and not args.overwrite and not args.dry_run:
            print(f"  {label}  skip (exists)")
            skipped += 1
            continue

        print(f"  {label}  {when:%Y-%m-%d %H}Z fxx={fxx}")
        try:
            state, meta = fetch_hour(when, fxx)
        except Exception as e:
            print(f"    FAILED: {type(e).__name__}: {e}")
            if "--debug" in sys.argv:
                traceback.print_exc()
            failed += 1
            continue

        if args.dry_run:
            print("\nDry run complete -- nothing written.")
            print(f"  state {state.shape} matches config "
                  f"[{config.N_CHANNELS}, {config.N_LEVELS}, Y, X]: "
                  f"{state.shape[:2] == (config.N_CHANNELS, config.N_LEVELS)}")
            return 0

        write_state(out_path, state, meta)
        ok += 1

        # Pace the downloads. This machine is shared with ~30 users behind one
        # connection; back-to-back GRIB fetches are noticeable to all of them.
        if args.pause > 0 and i < args.hours - 1:
            time.sleep(args.pause)

    print(f"\nDone. {ok} written, {skipped} skipped, {failed} failed -> {out_dir}")
    if failed:
        print("Failures are usually missing archive hours or a transient S3 error; "
              "re-run to fill gaps (existing files are skipped).")
    return 0 if ok or skipped else 1


if __name__ == "__main__":
    raise SystemExit(main())

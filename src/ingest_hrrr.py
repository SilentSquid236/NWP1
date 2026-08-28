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
import os
import shutil
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


def _download_subset(H, search, allow_full=True, verbose=True):
    """
    Download the matching GRIB messages and return the local path.

    Herbie's byte-range subsetting shells out to `curl`. When curl is missing
    -- or the .idx index cannot be fetched, so byte ranges cannot be computed
    -- the download quietly produces NO FILE, and Herbie then hands the
    non-existent path to cfgrib. The result is a FileNotFoundError deep in
    cfgrib that says nothing about the actual cause.

    So: download explicitly, verify the file exists and is non-empty, and
    report the real problem here. Fall back to the full file only if asked,
    since it is roughly 8x the transfer on a shared link.
    """
    path = None
    try:
        path = H.download(search)
    except TypeError:
        path = H.download(searchString=search)
    except Exception as e:
        if verbose:
            print(f"    subset download raised {type(e).__name__}: {e}")

    if path is not None and Path(path).exists() and Path(path).stat().st_size > 0:
        return Path(path)

    # Work out WHY, so the message is actionable.
    reasons = []
    if shutil.which("curl") is None:
        reasons.append("curl not found on PATH (Herbie uses it for byte-range "
                       "subsetting; install it, or use --full-file)")
    try:
        n = len(H.inventory(search))
        if n == 0:
            reasons.append(f"the search regex matched 0 of "
                           f"{len(H.inventory())} GRIB messages -- the "
                           f"variable naming does not fit this file")
        else:
            reasons.append(f"{n} messages matched, so the index is fine; "
                           f"the transfer itself produced no file")
    except Exception as e:
        reasons.append(f"could not read the GRIB index ({type(e).__name__}: "
                       f"{e}) -- without it byte ranges cannot be computed")

    if allow_full:
        if verbose:
            print(f"    subset produced no file; falling back to the FULL "
                  f"file (~8x the transfer)")
            for r in reasons:
                print(f"      cause: {r}")
        try:
            full = H.download()
            if full is not None and Path(full).exists() and Path(full).stat().st_size > 0:
                return Path(full)
        except Exception as e:
            reasons.append(f"full-file download also failed: "
                           f"{type(e).__name__}: {e}")

    raise RuntimeError("HRRR download produced no file.\n        "
                       + "\n        ".join(f"- {r}" for r in reasons))


def _open_hrrr(H, search, allow_full=True, verbose=True):
    """Download, verify, then open. Never hand cfgrib a path that may not exist."""
    import xarray as xr
    path = _download_subset(H, search, allow_full=allow_full, verbose=verbose)

    import cfgrib
    ds = cfgrib.open_datasets(str(path), backend_kwargs={"indexpath": ""})
    return ds


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
        var = resolve_variable(sub, name)
        arr = sub[var].values[:, ysl, xsl]       # [L, Y, X]
        planes.append(arr)

    return np.ascontiguousarray(np.stack(planes, axis=0), dtype=np.float32)


# Matches ":TMP:850 mb:anl", ":UGRD:500 mb:1 hour fcst", etc.
HRRR_SEARCH = r":(?:TMP|RH|UGRD|VGRD|HGT):\d+ mb:"

# cfgrib renames GRIB variables to CF short names on the way in, so the names
# we ASK Herbie for are not the names we get back:
#
#     GRIB   cfgrib
#     TMP -> t          RH   -> r
#     UGRD -> u         VGRD -> v
#     HGT -> gh
#
# Candidates are tried in order, so this also survives a cfgrib version that
# uses the long CF name or passes the GRIB name through unchanged.
CF_ALIASES = {
    "TMP":  ("t", "TMP", "temperature", "air_temperature"),
    "RH":   ("r", "RH", "relative_humidity", "r2"),
    "UGRD": ("u", "UGRD", "eastward_wind", "u_component_of_wind"),
    "VGRD": ("v", "VGRD", "northward_wind", "v_component_of_wind"),
    "HGT":  ("gh", "HGT", "geopotential_height", "z"),
}


def resolve_variable(ds, name):
    """
    Find the dataset variable holding our channel, trying known aliases.

    Returns the resolved name, or raises with everything that was tried and
    everything that was available -- the two facts needed to fix it.
    """
    for cand in CF_ALIASES.get(name, (name,)):
        if cand in ds:
            return cand
    raise KeyError(
        f"channel {name!r} not found. Tried {list(CF_ALIASES.get(name, (name,)))}; "
        f"dataset has {sorted(ds.data_vars)}")


def herbie_save_dir():
    """
    Where Herbie caches downloaded GRIB.

    Herbie defaults to ~/data, which on a shared machine is usually a small
    quota'd home directory. A failed write there does not raise -- the file
    simply never appears, and the failure surfaces later as a confusing
    FileNotFoundError from xarray. Keep the cache beside our own data.
    """
    d = Path(os.environ.get("NWP_HERBIE_DIR", config.DATA_ROOT / "herbie_cache"))
    d.mkdir(parents=True, exist_ok=True)
    return d


def fetch_hour(when, fxx=0, verbose=True, stride=1, allow_full=True):
    """Download one HRRR field set and return (array, metadata)."""
    from herbie import Herbie

    # One regex for all five variables across all isobaric levels. Herbie
    # byte-range downloads only the matching GRIB messages, so we never pull
    # the full ~130 MB file.
    #
    # NOTE THE LEADING COLON, and the absence of "^". HRRR index entries are
    # of the form ":HGT:50 mb:anl" -- they BEGIN with a colon, so anchoring
    # with ^ matches nothing at all. A regex that matches zero messages makes
    # Herbie download no file, which surfaces much later as a FileNotFoundError
    # from cfgrib about a path that was never created.
    search = HRRR_SEARCH

    H = Herbie(when.strftime("%Y-%m-%d %H:%M"), model="hrrr",
               product="prs", fxx=fxx, verbose=False,
               save_dir=herbie_save_dir())

    ds = _as_single_dataset(_open_hrrr(H, search, allow_full=allow_full,
                                       verbose=verbose))
    ysl, xsl = domain_slice(ds)
    state = extract_state(ds, ysl, xsl)

    lat = ds.latitude.values[ysl, xsl]
    lon = ds.longitude.values[ysl, xsl]
    lon = np.where(lon > 180, lon - 360, lon)

    if stride > 1:
        # Subsample AFTER the domain cut, so the geographic extent is
        # unchanged and only the spacing coarsens. build_grid() derives dx
        # from the domain extent and array shape, so the dynamics picks up
        # the coarser spacing automatically with no other change.
        state = np.ascontiguousarray(state[:, :, ::stride, ::stride])
        lat = np.ascontiguousarray(lat[::stride, ::stride])
        lon = np.ascontiguousarray(lon[::stride, ::stride])

    meta = {
        "valid_time": (when + timedelta(hours=fxx)).isoformat(),
        "run_time": when.isoformat(),
        "fxx": fxx,
        "channels": np.array(config.CHANNELS),
        "levels_hPa": np.array(config.PRESSURE_LEVELS),
        "lat": np.ascontiguousarray(lat, dtype=np.float32),
        "lon": np.ascontiguousarray(lon, dtype=np.float32),
        "source": f"HRRR prs f{fxx:02d}",
        "stride": stride,
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
    p.add_argument("--stride", type=int, default=1,
                   help="Subsample the grid by this factor. 1 = native 3 km, "
                        "2 = 6 km, 4 = 12 km. The hydrostatic core resolves "
                        "nothing extra at 3 km, and cost scales as stride^-3, "
                        "so 4 is the sensible default for iteration.")
    p.add_argument("--no-full-fallback", action="store_true",
                   help="Do not fall back to downloading the full GRIB file "
                        "when subsetting fails. The full file is ~8x the "
                        "transfer; on a shared link that matters.")
    p.add_argument("--debug", action="store_true",
                   help="Full traceback on every failure.")
    p.add_argument("--max-failures", type=int, default=3,
                   help="Give up after this many failures with no successes.")
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
    print(f"  stride         : {args.stride} "
          f"(~{3 * args.stride} km effective spacing)")
    print(f"  window         : {start:%Y-%m-%d %H}Z + {args.hours}h")
    print(f"  output         : {out_dir}")
    est = estimate_ingest_mb(args.hours, config.N_LEVELS, config.N_CHANNELS)
    print(f"  est. download  : ~{est['total_download_MB']:.0f} MB "
          f"({est['per_hour_download_MB']:.0f} MB/hour), "
          f"pausing {args.pause:.1f}s between files")
    print(f"  bandwidth      : shared link -- cap {max_mbps():.0f} MB/s "
          f"(NWP_MAX_MBPS to override)")
    print(f"  herbie cache   : {herbie_save_dir()}\n")

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
            state, meta = fetch_hour(when, fxx, stride=args.stride,
                                     allow_full=not args.no_full_fallback)
        except Exception as e:
            print(f"    FAILED: {type(e).__name__}: {e}")
            # Print the full chain on the FIRST failure. The surface error is
            # often FileNotFoundError from xarray, which says nothing about
            # why the download did not happen.
            if failed == 0 or args.debug:
                traceback.print_exc()
            failed += 1

            if args.dry_run:
                print("\nDry run stopped at the first failure. "
                      "Run diagnose_herbie.py for the cause.")
                return 1

            # Do not hammer a source that is clearly not working: 13 identical
            # failures tell us nothing 12 more times than 1 does.
            if failed >= args.max_failures and ok == 0:
                print(f"\nAborting: {failed} consecutive failures, none "
                      f"succeeded. Fix the cause before retrying "
                      f"(try: python diagnose_herbie.py).")
                return 1
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

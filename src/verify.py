"""
Verify a forecast against observations, and archive the pairing.

    python src/verify.py --forecast runs/analysis_20260904_00/forecast.npz
    python src/verify.py --forecast ... --report-only      # no network

WHY THIS IS THE URGENT ONE (P-07 in docs/PROBLEMS.md)

Every other open problem can be worked on tomorrow with no loss. This one
cannot: verification needs forecasts and observations PAIRED IN TIME, and an
observation not captured on the day is not obtainable afterwards in the form
that matters -- the archive of what this model predicted and what actually
happened. Learned post-processing needs a history that cannot be
reconstructed later. A day not archived is gone.

WHAT GETS STORED, AND IN WHAT ORDER

The order is deliberate and it is the main design decision here.

  1. RAW OBSERVATIONS, verbatim, compressed, written FIRST -- before any
     parsing, QC or matching is attempted. This is the irreplaceable part.
     Everything else can be recomputed from it.
  2. The forecast file, copied beside them.
  3. The matched pairs, as JSONL.

Matched pairs are DERIVED data. If the observation operator changes -- and it
will; the elevation correction is a standard lapse rate that is wrong on
exactly the nights it matters most -- every match can be recomputed from (1)
and (2). So a failure in parsing, QC or matching must never cost the raw
observations, which is why they are written before anything can throw.

VERIFICATION USES OBSERVATIONS ONLY

Never HRRR, never any model output. HRRR may seed a forecast; it may never
score one. That is a standing constraint of this project and the reason
`src/verification/fetchers.py` talks to IEM and MRMS rather than to an
archive of analyses.
"""

import argparse
import gzip
import json
import shutil
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
sys.path.insert(0, str(Path(__file__).resolve().parent / "dynamics"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "verification"))

import faulthandler
import signal

import numpy as np

import config

# Where is it stuck? `kill -USR1 <pid>` prints a traceback of every thread to
# stderr and the process carries on. Nothing to install, which matters on a
# server where nothing CAN be installed -- py-spy and friends are not an
# option here. A long IEM query and a genuine hang look identical from
# outside; this is how to tell them apart.
faulthandler.enable()
if hasattr(signal, "SIGUSR1"):
    faulthandler.register(signal.SIGUSR1)
from grid import CGrid
from sigma_operator import SigmaInterpolator
from observations import run_qc
from scoring import ForecastArchive, scores, report as score_report
import fetchers


# Networks covering the Northeast domain. ASOS only: automated, hourly,
# quality-controlled at source, and dense enough to score a 12 km grid.
NORTHEAST_NETWORKS = ["ME_ASOS", "NH_ASOS", "VT_ASOS", "MA_ASOS", "RI_ASOS",
                      "CT_ASOS", "NY_ASOS", "NJ_ASOS", "PA_ASOS"]

# What we verify, and against which model field.
#
# The names are the project's GRIB-style channel names -- the same ones
# `config.CHANNELS`, `observations.RANGE_LIMITS` and the fetchers already use.
# Inventing a second vocabulary here ("temperature", "u_wind") silently
# matched nothing at all: every observation fell through to "variable not
# verified" and the archive came out empty while every other check passed.
VERIFIED = {"TMP", "UGRD", "VGRD"}


# ---------------------------------------------------------------------------
# Archive layout
# ---------------------------------------------------------------------------

def archive_paths(root, run_time):
    """
    One directory per forecast run, so a run is self-describing and a partial
    day is obvious from the file listing rather than from a database.
    """
    root = Path(root)
    stamp = run_time.strftime("%Y%m%d_%H")
    d = root / stamp
    return {
        "dir": d,
        "raw": d / "observations",
        "forecast": d / "forecast.npz",
        "matches": d / "matches.jsonl",
        "meta": d / "run.json",
    }


def store_raw(raw_dir, name, text):
    """
    Write a fetched observation payload verbatim, compressed.

    Verbatim matters: a parser bug found in six months should be fixable
    against what the service actually sent, not against what we thought it
    meant at the time.
    """
    raw_dir.mkdir(parents=True, exist_ok=True)
    path = raw_dir / f"{name}.csv.gz"
    with gzip.open(path, "wt", encoding="utf-8") as f:
        f.write(text)
    return path


def load_raw(raw_dir, name):
    path = Path(raw_dir) / f"{name}.csv.gz"
    if not path.exists():
        return None
    with gzip.open(path, "rt", encoding="utf-8") as f:
        return f.read()


# ---------------------------------------------------------------------------
# Forecast loading
# ---------------------------------------------------------------------------

def load_forecast(path):
    """Load a sigma-coordinate forecast written by forecast.py."""
    z = np.load(path, allow_pickle=False)
    need = ("times_s", "u", "v", "theta", "pi", "sigma", "p_top", "terrain")
    missing = [k for k in need if k not in z.files]
    if missing:
        raise SystemExit(
            f"{path} is missing {missing}.\n"
            f"It has {sorted(z.files)}. A forecast written before the sigma "
            f"port has no pi/sigma/terrain and cannot be verified against "
            f"station elevations.")
    return {k: z[k] for k in z.files}


def build_interpolator(fc, snapshot):
    """Grid and operator for one snapshot of a forecast."""
    ny, nx = fc["theta"].shape[-2:]
    lat0 = 0.5 * (config.DOMAIN["lat_min"] + config.DOMAIN["lat_max"])
    dy = (config.DOMAIN["lat_max"] - config.DOMAIN["lat_min"]) * 111_132.0 / ny
    dx = (config.DOMAIN["lon_max"] - config.DOMAIN["lon_min"]) * 111_320.0 * \
        np.cos(np.radians(lat0)) / nx
    gr = CGrid(nx, ny, dx, dy, edge_mode="replicate")
    return SigmaInterpolator(gr, config.DOMAIN, fc["sigma"], float(fc["p_top"]),
                             fc["pi"][snapshot], fc["terrain"])


# ---------------------------------------------------------------------------
# Matching
# ---------------------------------------------------------------------------

def match_snapshot(fc, snapshot, obs, valid_time, lead_hours,
                   window_min=30):
    """
    Pair one forecast snapshot with the observations valid near it.

    Skips are counted by reason rather than silently dropped: a sudden fall in
    the match count should be visible as a change in one reason, not as a
    smaller number.
    """
    op = build_interpolator(fc, snapshot)
    theta = fc["theta"][snapshot]
    u = fc["u"][snapshot]
    v = fc["v"][snapshot]

    matches = []
    skipped = {}

    def skip(why):
        skipped[why] = skipped.get(why, 0) + 1

    for o in obs:
        if not getattr(o, "passed", True):
            skip("failed qc")
            continue
        dt_min = abs((o.time - valid_time).total_seconds()) / 60.0
        if dt_min > window_min:
            skip("outside time window")
            continue

        info = {}
        if o.variable == "TMP":
            # Surface temperature goes through the elevation correction; an
            # upper-air sounding does not, because it already carries its own
            # pressure.
            if o.pressure is None:
                value, info = op.station_temperature(
                    theta, o.lat, o.lon, getattr(o, "elevation", None))
            else:
                value = op.temperature(theta, o.lat, o.lon, o.pressure)
        elif o.variable == "UGRD":
            value = op.at_observation(u, o.lat, o.lon, o.pressure)
        elif o.variable == "VGRD":
            value = op.at_observation(v, o.lat, o.lon, o.pressure)
        else:
            skip(f"variable not verified: {o.variable}")
            continue

        if value is None:
            skip("outside domain or column")
            continue

        rec = {
            "time": o.time.isoformat(),
            "valid_time": valid_time.isoformat(),
            "station": o.station,
            "source": o.source,
            "variable": o.variable,
            "lat": o.lat, "lon": o.lon,
            "pressure": o.pressure,
            "station_elev_m": getattr(o, "elevation", None),
            "lead_hours": lead_hours,
            "forecast": float(value),
            "observation": float(o.value),
            "error": float(value - o.value),
            "error_std": getattr(o, "error_std", None),
            "time_offset_min": round(dt_min, 1),
        }
        rec.update(info)
        matches.append(rec)

    return matches, skipped


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def verify(forecast_path, archive_root, run_time=None, window_min=30,
           report_only=False, networks=None, verbose=True):
    fc = load_forecast(forecast_path)
    times_s = np.asarray(fc["times_s"], dtype=float)

    run_time = run_time or datetime.now(timezone.utc).replace(
        minute=0, second=0, microsecond=0, tzinfo=None)
    paths = archive_paths(archive_root, run_time)
    paths["dir"].mkdir(parents=True, exist_ok=True)

    valid_times = [run_time + timedelta(seconds=float(t)) for t in times_s]
    first, last = valid_times[0], valid_times[-1]

    if verbose:
        print(f"  forecast       : {forecast_path}")
        print(f"  run time       : {run_time:%Y-%m-%d %H}Z")
        print(f"  valid times    : {len(valid_times)}, "
              f"{first:%H}Z to {last:%H}Z")
        print(f"  archive        : {paths['dir']}")

    # ---- 1. RAW OBSERVATIONS FIRST -------------------------------------
    #
    # Before parsing, before QC, before matching. Everything downstream is
    # recomputable; this is not.
    name = f"asos_{first:%Y%m%d%H}_{last:%Y%m%d%H}"
    text = load_raw(paths["raw"], name)
    if text is not None:
        if verbose:
            print(f"  observations   : cached ({len(text)/1e3:.0f} kB)")
    elif report_only:
        raise SystemExit("--report-only, but no cached observations for this "
                         "window. Run once without it to fetch them.")
    else:
        if verbose:
            print(f"  observations   : fetching ASOS, "
                  f"{len(networks or NORTHEAST_NETWORKS)} networks...")
        text = fetchers.fetch_asos_text(
            networks or NORTHEAST_NETWORKS,
            first - timedelta(minutes=window_min),
            last + timedelta(minutes=window_min))
        store_raw(paths["raw"], name, text)
        if verbose:
            print(f"                   stored {len(text)/1e3:.0f} kB verbatim")

    # ---- 2. THE FORECAST, BESIDE THEM ----------------------------------
    if not paths["forecast"].exists():
        shutil.copy2(forecast_path, paths["forecast"])

    # ---- 3. DERIVED: parse, QC, match ----------------------------------
    obs = fetchers.parse_asos_csv(text)
    kept, rejected, qc_info = run_qc(obs)
    if verbose:
        print(f"  parsed         : {len(obs)} observations, "
              f"{len(kept)} passed QC, {len(rejected)} rejected")

    archive = ForecastArchive(paths["matches"])
    seen = {(r.get("valid_time"), r.get("station"), r.get("variable"))
            for r in archive.load()}

    all_matches, all_skips = [], {}
    for i, vt in enumerate(valid_times):
        lead = float(times_s[i]) / 3600.0
        m, sk = match_snapshot(fc, i, kept, vt, lead, window_min)
        m = [r for r in m
             if (r["valid_time"], r["station"], r["variable"]) not in seen]
        all_matches.extend(m)
        for k, n in sk.items():
            all_skips[k] = all_skips.get(k, 0) + n

    if all_matches:
        archive.append(all_matches, run_time=run_time,
                       extra={"model": "primitive_sigma"})

    with open(paths["meta"], "w") as f:
        json.dump({
            "run_time": run_time.isoformat(),
            "forecast": str(forecast_path),
            "valid_times": [t.isoformat() for t in valid_times],
            "observations_raw": name,
            "n_observations": len(obs),
            "n_passed_qc": len(kept),
            "qc": qc_info,
            "n_matches": len(all_matches),
            "skipped": all_skips,
            "written": datetime.now(timezone.utc).isoformat(),
        }, f, indent=2)

    if verbose:
        print(f"  matched        : {len(all_matches)} new pairs "
              f"({len(seen)} already archived)")
        if all_skips:
            for k, n in sorted(all_skips.items(), key=lambda kv: -kv[1]):
                print(f"                   skipped {n:5d}  {k}")

    return all_matches, paths


def main():
    p = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--forecast", required=True,
                   help="forecast.npz written by src/forecast.py")
    p.add_argument("--archive", default=None,
                   help="Archive root (default: <data>/verification)")
    p.add_argument("--run-time", default=None,
                   help="Forecast run time, e.g. 2026-09-04T00. Defaults to "
                        "the current hour, which is only right if you are "
                        "verifying a forecast you just made.")
    p.add_argument("--window-min", type=float, default=30.0,
                   help="Match observations within this many minutes")
    p.add_argument("--report-only", action="store_true",
                   help="Use cached observations; touch no network.")
    p.add_argument("--summary", action="store_true",
                   help="Print scores for the whole archive and exit.")
    args = p.parse_args()

    root = Path(args.archive or (config.DATA_ROOT / "verification"))

    if args.summary:
        recs = []
        for m in sorted(root.glob("*/matches.jsonl")):
            recs.extend(ForecastArchive(m).load())
        if not recs:
            print(f"No archived matches under {root}.")
            return 1
        print(f"\nVerification archive: {root}")
        print(f"{len(recs)} pairs from "
              f"{len(set(r['run_time'] for r in recs))} runs\n")
        for var in sorted(set(r["variable"] for r in recs)):
            sub = [r for r in recs if r["variable"] == var]
            print(f"{var}:")
            print(score_report(sub))
        return 0

    run_time = None
    if args.run_time:
        for fmt in ("%Y-%m-%dT%H", "%Y-%m-%d %H", "%Y-%m-%d"):
            try:
                run_time = datetime.strptime(args.run_time, fmt)
                break
            except ValueError:
                continue
        if run_time is None:
            p.error(f"could not parse --run-time {args.run_time!r}")

    print("Forecast verification")
    print(config.describe())
    matches, paths = verify(args.forecast, root, run_time=run_time,
                            window_min=args.window_min,
                            report_only=args.report_only)

    if matches:
        print()
        print(score_report(matches))
    print(f"\nArchived -> {paths['dir']}")
    print("The raw observations are stored verbatim: every match here can be "
          "recomputed if the observation operator changes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

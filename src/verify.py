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
#
# Widened 2026-09-22 to every state and province the domain touches: the old
# list stopped at Pennsylvania, so the southern and western thirds of the
# domain (DE, MD, VA, WV, OH) and all of Canada were never scored. Stations
# outside the domain are dropped by the operator, not by this list.
NORTHEAST_NETWORKS = ["ME_ASOS", "NH_ASOS", "VT_ASOS", "MA_ASOS", "RI_ASOS",
                      "CT_ASOS", "NY_ASOS", "NJ_ASOS", "PA_ASOS", "DE_ASOS",
                      "MD_ASOS", "VA_ASOS", "WV_ASOS", "OH_ASOS",
                      "CA_ON_ASOS", "CA_QC_ASOS", "CA_NB_ASOS", "CA_NS_ASOS"]

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


def forecast_mslp(fc, snapshot):
    """Sea-level pressure (Pa) of one snapshot, as maps/derive.py computes it."""
    import sys as _sys
    from pathlib import Path as _P
    _sys.path.insert(0, str(_P(__file__).resolve().parent / "maps"))
    from derive import column_state, mslp
    ter = np.asarray(fc["terrain"], float)
    p, T, Z, ps = column_state(np.asarray(fc["theta"][snapshot], float), np.asarray(fc["pi"][snapshot], float),
                               fc["sigma"], float(fc["p_top"]), ter)
    return mslp(ps, T, Z, ter)


def hourly_rain(fc, snapshot):
    """
    Rain (mm) in the hour ending at `snapshot`, from the accumulated "rain"
    field: the difference from the snapshot one hour earlier, or the
    accumulation itself when the snapshot is at +1 h. False when neither
    exists (output not hourly).
    """
    t = np.asarray(fc["times_s"], float)
    r = np.asarray(fc["rain"], float)
    tk = t[snapshot]
    if abs(tk - 3600.0) < 1.0:
        return np.maximum(r[snapshot], 0.0)
    j = np.where(np.abs(t - (tk - 3600.0)) < 1.0)[0]
    if len(j) == 0:
        return False
    return np.maximum(r[snapshot] - r[j[0]], 0.0)


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

Z0_FIXED = 0.1                 # the operators' roughness before 2026-10-10
Z0_LAND_DEFAULT, Z0_SEA_DEFAULT = 1.0, 0.0002     # forecast.py defaults since test AM


def model_z0_map(fc):
    """
    The roughness map (m) the forecast's drag used, and where it came from.

    Forecasts written since 2026-10-10 carry it ("z0"). For older ones the
    default land/sea map is rebuilt as forecast.py builds it: land = raw ETOPO
    terrain above 0 m, 1.0 m over land and 0.0002 m over water (test AM). A
    run made with --z0 or other --z0-land/--z0-sea values is NOT recoverable
    this way; those runs should be verified with --z0-operator fixed.
    """
    if "z0" in fc:
        return np.asarray(fc["z0"], float), "forecast file"
    ter = np.asarray(fc["terrain"], float)
    raw = Path(config.DATA_ROOT) / "static" / f"terrain_etopo_{ter.shape[0]}x{ter.shape[1]}.npz"
    land = (np.load(raw)["terrain"] > 0.0) if raw.exists() else (ter > 0.0)
    return np.where(land, Z0_LAND_DEFAULT, Z0_SEA_DEFAULT), f"default land/sea map ({raw.name if raw.exists() else 'run terrain'})"


def match_snapshot(fc, snapshot, obs, valid_time, lead_hours,
                   window_min=30, wind_operator="log10m", surface_operator="standard",
                   z0_operator="fixed", z0_map=None):
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
    tg = None
    if surface_operator == "similarity":
        if "tg" not in fc:
            raise ValueError("--surface-operator similarity needs a forecast run with "
                             "--land-surface (no ground temperature 'tg' in the file)")
        tg = fc["tg"][snapshot]

    matches = []
    skipped = {}
    mslp2d = None            # sea-level pressure field, built on first use (P-76)
    acc1h = None             # one-hour rain field, built on first use (S7)
    logz0 = None if z0_map is None else np.log(np.maximum(np.asarray(z0_map, float), 1e-6))

    def z0_at(lat, lon):
        # the model's roughness at the station: log-z0 interpolated, so a
        # coastal station between 1 m and 0.0002 m gets a geometric mean
        if z0_operator != "model" or logz0 is None:
            return Z0_FIXED
        g = op.horizontal(logz0, lat, lon)
        return Z0_FIXED if g is None else float(np.exp(g))

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
            if o.pressure is None and tg is not None:
                value, info = op.station_temperature_similarity(
                    theta, u, v, tg, o.lat, o.lon, getattr(o, "elevation", None),
                    z0=z0_at(o.lat, o.lon))
            elif o.pressure is None:
                value, info = op.station_temperature(
                    theta, o.lat, o.lon, getattr(o, "elevation", None))
            else:
                value = op.temperature(theta, o.lat, o.lon, o.pressure)
        elif o.variable in ("UGRD", "VGRD"):
            value = op.at_observation(u if o.variable == "UGRD" else v,
                                      o.lat, o.lon, o.pressure)
            # P-68: a 10 m anemometer against the lowest level (~300 m above
            # the ground) needs a reduction, or every wind score is mostly
            # the speed ratio. "lowest" keeps the old operator for comparison.
            if value is not None and o.pressure is None and tg is not None:
                d = op.surface_similarity(theta, u, v, tg, o.lat, o.lon, z0=z0_at(o.lat, o.lon))
                if d is not None:
                    info = {"wind_operator": "similarity", "wind_factor": d["wind_factor"],
                            "z0": z0_at(o.lat, o.lon),
                            "zeta": d["zeta"], "Ri_bulk": d["Ri_bulk"],
                            "model_level_agl_m": d["model_level_agl_m"],
                            "forecast_lowest_level": float(value)}
                    value = value * d["wind_factor"]
            elif (value is not None and o.pressure is None
                    and wind_operator == "log10m"):
                f, agl = op.wind_10m_factor(theta, o.lat, o.lon, z0=z0_at(o.lat, o.lon))
                if f is not None:
                    info = {"wind_operator": "log10m", "wind_factor": f,
                            "z0": z0_at(o.lat, o.lon),
                            "model_level_agl_m": agl,
                            "forecast_lowest_level": float(value)}
                    value = value * f
        elif o.variable == "PMSL" and o.pressure is None:
            # P-76: sea-level pressure, reduced exactly as the maps reduce it
            # (maps/derive.py), then interpolated to the station
            if mslp2d is None:
                mslp2d = forecast_mslp(fc, snapshot)
            value = op.horizontal(mslp2d, o.lat, o.lon)
        elif o.variable == "APCP1H" and o.pressure is None:
            # S7: one-hour precipitation. The forecast's "rain" is accumulated
            # surface rain (mm) since the start; the hour ending at this
            # snapshot is its difference from the snapshot one hour earlier
            # (or from zero at the start). METARs report the hour ending at
            # about :53, within the 30 min window of the snapshot.
            if "rain" not in fc:
                skip("precipitation: forecast has no rain")
                continue
            if acc1h is None:
                acc1h = hourly_rain(fc, snapshot)
            if acc1h is False:
                skip("precipitation: no snapshot one hour earlier")
                continue
            value = op.horizontal(acc1h, o.lat, o.lon)
        elif o.variable == "DPT" and o.pressure is None:
            # CAM stage S7: 2 m dewpoint from a moist forecast's lowest-level qv
            if "qv" not in fc:
                skip("dewpoint: forecast has no moisture")
                continue
            value, info = op.station_dewpoint(fc["qv"][snapshot], theta, o.lat, o.lon,
                                              getattr(o, "elevation", None))
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

def resolve_surface_operator(fc):
    """
    The 'auto' surface operator (default since 2026-10-03).

    Similarity needs the model's ground temperature, so a forecast run with the
    land surface is scored with it. Persistence is NOT, even when it carries a
    held ground temperature: holding the cycle-time stability all day reduces
    its daytime wind as if it were still night (the 'psim' reference, test AK),
    which is not a fair do-nothing forecast. Persistence keeps the standard
    operator, as in every comparison since test AD.
    """
    stopped = str(np.asarray(fc["stopped"]).item()) if "stopped" in fc else ""
    return "similarity" if ("tg" in fc and stopped != "persistence") else "standard"


def verify(forecast_path, archive_root, run_time=None, window_min=30,
           report_only=False, networks=None, verbose=True,
           wind_operator="log10m", surface_operator="auto", copy_forecast=True,
           z0_operator="model"):
    fc = load_forecast(forecast_path)
    z0_map, z0_src = (model_z0_map(fc) if z0_operator == "model" else (None, f"fixed {Z0_FIXED:g} m"))
    if verbose:
        print(f"  operator z0    : {z0_src}")
    times_s = np.asarray(fc["times_s"], dtype=float)
    if surface_operator == "auto":
        surface_operator = resolve_surface_operator(fc)
        if verbose:
            print(f"  surface op.    : {surface_operator} (auto)")

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
    # (copy_forecast=False for the 3 km runs: 1.6 GB each, and they stay in
    # tensors_3km beside their verified.json; P-75)
    if copy_forecast and not paths["forecast"].exists():
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
        m, sk = match_snapshot(fc, i, kept, vt, lead, window_min,
                               wind_operator=wind_operator,
                               surface_operator=surface_operator,
                               z0_operator=z0_operator, z0_map=z0_map)
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
            "wind_operator": wind_operator,
            "surface_operator": surface_operator,
            "z0_operator": z0_src,
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
    p.add_argument("--wind-operator", choices=("lowest", "log10m"), default="log10m",
                   help="Surface wind: reduced to 10 m with the model's "
                        "neutral log law (the default since 2026-10-02, "
                        "P-68), or the lowest model level as it is (the "
                        "old default; scores before 2026-10-02 used it)")
    # DEFAULT CHANGED 2026-10-10 (P-78, test P2): the drag's own roughness.
    # Every wind score reported before 2026-10-10 used "fixed" (0.1 m).
    p.add_argument("--z0-operator", choices=("fixed", "model"), default="model",
                   help="roughness in the 10 m wind and 2 m temperature operators: model = the "
                        "drag's own map (the forecast's z0 field, or the default land 1.0 m / "
                        "water 0.0002 m map rebuilt for older files; default since 2026-10-10, "
                        "test P2), fixed = 0.1 m everywhere (all scores before)")
    p.add_argument("--surface-operator", choices=("auto", "standard", "similarity"),
                   default="auto",
                   help="standard: 2 m temperature from the lowest level by a "
                        "standard lapse rate, wind by --wind-operator. "
                        "similarity: both from Monin-Obukhov similarity between "
                        "the model's ground temperature and the lowest level "
                        "(needs a --land-surface forecast; test AK). auto (default "
                        "since 2026-10-03): similarity when the forecast carries a "
                        "ground temperature and is not persistence, else standard")
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
                            report_only=args.report_only,
                            wind_operator=args.wind_operator,
                            surface_operator=args.surface_operator,
                            z0_operator=args.z0_operator)

    if matches:
        print()
        print(score_report(matches))
    print(f"\nArchived -> {paths['dir']}")
    print("The raw observations are stored verbatim: every match here can be "
          "recomputed if the observation operator changes.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

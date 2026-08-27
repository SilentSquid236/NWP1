"""
Verification scoring and the forecast-observation archive.

Two jobs:

1. SCORE a forecast against observations, broken down by whatever dimension
   matters -- station, lead time, level, variable.
2. ARCHIVE every matched pair, because post-processing cannot be trained on
   data that was never stored, and none of it is recoverable later.

The archive is deliberately a plain append-only JSONL file. It has to survive
crashes, be readable by anything, and be trivially concatenated across runs.
A database would be faster to query and much easier to lose.
"""

import json
from collections import defaultdict
from datetime import datetime
from pathlib import Path

import numpy as np


# ---------------------------------------------------------------------------
# Scores
# ---------------------------------------------------------------------------

def scores(pairs):
    """
    pairs : iterable of (forecast, observation)

    bias : mean error -- the systematic part, what post-processing removes
    rmse : total error
    mae  : robust to outliers
    """
    f = np.array([p[0] for p in pairs], dtype=float)
    o = np.array([p[1] for p in pairs], dtype=float)
    m = np.isfinite(f) & np.isfinite(o)
    f, o = f[m], o[m]

    if len(f) == 0:
        return {"n": 0, "bias": np.nan, "rmse": np.nan, "mae": np.nan}

    e = f - o
    out = {
        "n": int(len(f)),
        "bias": float(e.mean()),
        "rmse": float(np.sqrt((e**2).mean())),
        "mae": float(np.abs(e).mean()),
    }
    if len(f) > 1 and o.std() > 0 and f.std() > 0:
        out["corr"] = float(np.corrcoef(f, o)[0, 1])
    return out


def scores_by(matches, key):
    """
    Break scores down by an attribute of the match dicts.

    Aggregate scores hide the failures that matter. A model can look fine
    overall while being badly wrong at one station, one level, or one lead
    time -- and that is usually where the physics bug is.
    """
    groups = defaultdict(list)
    for m in matches:
        groups[m[key]].append((m["forecast"], m["observation"]))
    return {k: scores(v) for k, v in sorted(groups.items(), key=lambda kv: str(kv[0]))}


def skill_vs_reference(matches, reference_matches):
    """
    Fractional RMSE improvement over a reference forecast.

    The reference should be something trivial -- persistence, climatology, or
    HRRR. Beating nothing is not a result. Negative means the reference is
    better, which is information worth having early.
    """
    a = scores([(m["forecast"], m["observation"]) for m in matches])
    b = scores([(m["forecast"], m["observation"]) for m in reference_matches])
    if not a["n"] or not b["n"] or not np.isfinite(b["rmse"]) or b["rmse"] == 0:
        return np.nan
    return float(1.0 - a["rmse"] / b["rmse"])


# ---------------------------------------------------------------------------
# Matching model to observations
# ---------------------------------------------------------------------------

def match_forecast_to_obs(field3d, interpolator, observations, valid_time,
                          variable, lead_hours=None, time_window_min=30):
    """
    Interpolate the model to each observation and pair them up.

    Observations outside the domain, outside the model's vertical range, or
    outside the time window are SKIPPED rather than matched approximately.
    Every skip is counted, so a sudden drop in match count is visible instead
    of silently shrinking the sample.
    """
    matches = []
    skipped = defaultdict(int)

    for obs in observations:
        if obs.variable != variable:
            skipped["wrong variable"] += 1
            continue
        if not obs.passed:
            skipped["failed qc"] += 1
            continue
        dt_min = abs((obs.time - valid_time).total_seconds()) / 60.0
        if dt_min > time_window_min:
            skipped["outside time window"] += 1
            continue

        model_value = interpolator.at_observation(field3d, obs.lat, obs.lon,
                                                  obs.pressure)
        if model_value is None:
            skipped["outside domain or levels"] += 1
            continue

        matches.append({
            "time": obs.time.isoformat(),
            "station": obs.station,
            "source": obs.source,
            "variable": variable,
            "lat": obs.lat, "lon": obs.lon,
            "pressure": obs.pressure,
            "lead_hours": lead_hours,
            "forecast": model_value,
            "observation": obs.value,
            "error": model_value - obs.value,
            "error_std": obs.error_std,
        })

    return matches, dict(skipped)


# ---------------------------------------------------------------------------
# Archive
# ---------------------------------------------------------------------------

class ForecastArchive:
    """
    Append-only JSONL store of matched forecast/observation pairs.

    This exists because learned post-processing needs a history that cannot be
    reconstructed after the fact. Archiving has to start with the very first
    forecast; anything not written is gone.
    """

    def __init__(self, path):
        self.path = Path(path)
        self.path.parent.mkdir(parents=True, exist_ok=True)

    def append(self, matches, run_time=None, extra=None):
        run = (run_time or datetime.utcnow()).isoformat()
        with open(self.path, "a") as f:
            for m in matches:
                rec = dict(m)
                rec["run_time"] = run
                if extra:
                    rec.update(extra)
                f.write(json.dumps(rec) + "\n")
        return len(matches)

    def load(self, variable=None, station=None, since=None):
        if not self.path.exists():
            return []
        out = []
        with open(self.path) as f:
            for line in f:
                line = line.strip()
                if not line:
                    continue
                try:
                    r = json.loads(line)
                except json.JSONDecodeError:
                    continue        # a truncated final line after a crash
                if variable and r.get("variable") != variable:
                    continue
                if station and r.get("station") != station:
                    continue
                if since and r.get("time", "") < since:
                    continue
                out.append(r)
        return out

    def __len__(self):
        return len(self.load())

    def __repr__(self):
        return f"ForecastArchive({self.path}, {len(self)} records)"


def report(matches, by=("station", "lead_hours")):
    """Human-readable verification summary."""
    lines = []
    overall = scores([(m["forecast"], m["observation"]) for m in matches])
    lines.append(f"Overall: n={overall['n']}  bias={overall['bias']:+.3f}  "
                 f"rmse={overall['rmse']:.3f}  mae={overall['mae']:.3f}")
    for key in by:
        if not matches or key not in matches[0]:
            continue
        lines.append(f"\nBy {key}:")
        for k, s in scores_by(matches, key).items():
            if s["n"]:
                lines.append(f"  {str(k):>12}  n={s['n']:4d}  "
                             f"bias={s['bias']:+7.3f}  rmse={s['rmse']:7.3f}")
    return "\n".join(lines)

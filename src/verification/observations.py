"""
Observation records and quality control.

One flat record type for every source -- radiosonde, METAR, buoy -- so that QC
and verification stay source-agnostic. Anything source-specific belongs in the
fetcher that produces these, not downstream.

QUALITY CONTROL PHILOSOPHY

A bad observation is worse than no observation. A single mis-decoded value
injects a large false innovation into assimilation, or a large false error
into verification, and either one is invisible in aggregate statistics.

Every rejection is recorded with a reason. Silent dropping makes "no obs
available" indistinguishable from "all obs rejected", and those need
completely different fixes.
"""

from dataclasses import dataclass, field
from datetime import datetime, timedelta
from typing import Optional

import numpy as np


@dataclass
class Observation:
    time: datetime
    lat: float
    lon: float
    variable: str                    # TMP, RH, UGRD, VGRD, HGT
    value: float                     # SI units, converted at ingest
    source: str                      # raob | metar | buoy | ...
    station: str = ""
    elevation: Optional[float] = None    # m MSL, station height
    pressure: Optional[float] = None     # Pa; None for surface obs
    error_std: float = 1.0               # obs + representativeness error

    # Filled in by QC / verification
    qc_flag: str = "unchecked"
    qc_reason: str = ""

    @property
    def is_surface(self):
        return self.pressure is None

    @property
    def passed(self):
        return self.qc_flag == "good"


# Physical bounds. Values outside these are decode errors, not weather.
RANGE_LIMITS = {
    "TMP": (170.0, 340.0),        # K
    "RH": (0.0, 105.0),           # %, slight slack for supersaturation
    "UGRD": (-150.0, 150.0),      # m/s
    "VGRD": (-150.0, 150.0),
    "HGT": (-500.0, 60000.0),     # gpm
}

# Default observation + representativeness error by source and variable.
# Representativeness usually DOMINATES instrument error: we are comparing a
# point measurement with a grid-cell average several km across.
DEFAULT_ERROR_STD = {
    ("raob", "TMP"): 1.0,
    ("raob", "RH"): 10.0,
    ("raob", "UGRD"): 2.5,
    ("raob", "VGRD"): 2.5,
    ("raob", "HGT"): 10.0,
    ("metar", "TMP"): 1.5,
    ("metar", "RH"): 12.0,
    ("metar", "UGRD"): 2.0,
    ("metar", "VGRD"): 2.0,
}


def default_error_std(source, variable, fallback=2.0):
    return DEFAULT_ERROR_STD.get((source, variable), fallback)


# ---------------------------------------------------------------------------
# Quality control
# ---------------------------------------------------------------------------

def range_check(obs):
    """Reject physically impossible values."""
    lim = RANGE_LIMITS.get(obs.variable)
    if lim is None:
        return True, ""
    lo, hi = lim
    if not np.isfinite(obs.value):
        return False, "non-finite value"
    if not (lo <= obs.value <= hi):
        return False, f"out of range [{lo}, {hi}]: {obs.value:g}"
    return True, ""


def gross_error_check(obs, background, n_sigma=5.0):
    """
    Reject when the observation disagrees wildly with the model background.

    This is the workhorse: it catches decode errors, misplaced stations, and
    unit mistakes. The threshold must be generous -- at 3 sigma you start
    rejecting real extreme weather, which is precisely what you most want to
    capture.
    """
    if background is None or not np.isfinite(background):
        return True, ""
    d = abs(obs.value - background)
    if d > n_sigma * obs.error_std:
        return False, (f"gross error: |{obs.value:g} - {background:g}| = "
                       f"{d:.2f} > {n_sigma}*{obs.error_std:g}")
    return True, ""


def buddy_check(obs, neighbours, n_sigma=4.0, min_buddies=2):
    """
    Reject an observation that disagrees sharply with nearby observations of
    the same variable. Catches a broken instrument that reports plausible but
    wrong values -- which the range and gross-error checks both miss.

    Skipped when there are too few neighbours: in sparse regions there is no
    consensus to compare against, and rejecting on one buddy is worse than
    not checking.
    """
    vals = [o.value for o in neighbours
            if o.variable == obs.variable and o.station != obs.station
            and np.isfinite(o.value)]
    if len(vals) < min_buddies:
        return True, "too few buddies to check"

    med = float(np.median(vals))
    spread = float(np.median(np.abs(np.asarray(vals) - med))) * 1.4826  # robust sigma
    spread = max(spread, 0.5 * obs.error_std)      # never divide by ~0

    if abs(obs.value - med) > n_sigma * spread:
        return False, (f"buddy check: {obs.value:g} vs neighbour median "
                       f"{med:g} (spread {spread:.2f})")
    return True, ""


def run_qc(observations, backgrounds=None, blacklist=None,
           gross_sigma=5.0, buddy_sigma=4.0, buddy_radius_km=150.0):
    """
    Apply all gates in order. Returns (kept, rejected, summary).

    backgrounds : optional dict station-> model value at the ob location, for
                  the gross-error check. Without it that gate is skipped.
    """
    blacklist = set(blacklist or ())
    summary = {"total": len(observations), "good": 0, "blacklist": 0,
               "range": 0, "gross": 0, "buddy": 0}

    # Group by (variable, time) for the buddy check.
    kept, rejected = [], []

    for obs in observations:
        if obs.station in blacklist:
            obs.qc_flag, obs.qc_reason = "rejected", "blacklisted station"
            summary["blacklist"] += 1
            rejected.append(obs)
            continue

        ok, why = range_check(obs)
        if not ok:
            obs.qc_flag, obs.qc_reason = "rejected", why
            summary["range"] += 1
            rejected.append(obs)
            continue

        bg = (backgrounds or {}).get(obs.station)
        ok, why = gross_error_check(obs, bg, gross_sigma)
        if not ok:
            obs.qc_flag, obs.qc_reason = "rejected", why
            summary["gross"] += 1
            rejected.append(obs)
            continue

        kept.append(obs)

    # Buddy check on survivors, using geographic neighbours.
    final = []
    for obs in kept:
        nb = [o for o in kept
              if o is not obs
              and o.variable == obs.variable
              and abs((o.time - obs.time).total_seconds()) < 1800
              and _haversine_km(obs.lat, obs.lon, o.lat, o.lon) < buddy_radius_km]
        ok, why = buddy_check(obs, nb, buddy_sigma)
        if not ok:
            obs.qc_flag, obs.qc_reason = "rejected", why
            summary["buddy"] += 1
            rejected.append(obs)
        else:
            obs.qc_flag, obs.qc_reason = "good", why
            summary["good"] += 1
            final.append(obs)

    return final, rejected, summary


def _haversine_km(lat1, lon1, lat2, lon2):
    R = 6371.0
    p1, p2 = np.radians(lat1), np.radians(lat2)
    dp = p2 - p1
    dl = np.radians(lon2 - lon1)
    a = np.sin(dp / 2)**2 + np.cos(p1) * np.cos(p2) * np.sin(dl / 2)**2
    return float(2 * R * np.arcsin(np.sqrt(a)))

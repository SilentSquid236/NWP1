"""
Build the cycle-time analysis from observations: the model's initial state.

    fields, meta = build_analysis(cycle, observations, terrain, domain,
                                  previous_forecast=path_or_None)

ORDER, AND WHY

  1. FIRST GUESS on the pressure levels. The previous run's forecast valid at
     this time if one exists (prompt 100); otherwise the mean of the
     soundings available, level by level; otherwise a standard atmosphere.
  2. UPPER AIR. Soundings are interpolated in log(p) onto the analysis levels
     (never extrapolated beyond the profile), and their innovations against
     the first guess are analysed level by level. Far from every sounding the
     first guess stands -- which is what "use the previous run where
     soundings can't be found" means at a missing site.
  3. SEA-LEVEL PRESSURE from ASOS and buoys: the densest, best-measured field
     there is, and the anchor for every height.
  4. HEIGHTS, hydrostatically: 1000 hPa from sea-level pressure, then upward
     through the analysed virtual temperature. Never analysed independently:
     heights perturbed separately from temperature were one of the four worst
     test-setup defects in this project (a state that is not hydrostatic is a
     pressure-gradient force nothing balances).
  5. SURFACE BLEND. Surface temperature, humidity and wind innovations (at
     the station, against the first guess at the grid terrain height) are
     analysed in 2-D and applied to the lowest levels with a weight that
     falls from 1 at the ground to 0 at `blend_depth_m` above it.
  6. HEIGHTS AGAIN, because step 5 changed the temperatures.

WITHHELD STATIONS

One ASOS station in five (by a hash of its identifier, so the set is the same
every run) is never used. The analysis is scored against those alone; any
other score of a cycle-time analysis is a score against its own inputs.
"""

import json
import zlib
from collections import defaultdict

import numpy as np

import barnes
import geo

RD = 287.05
G0 = 9.80665
EPS = 0.622
LAPSE = 0.0065
CHANNELS = ("TMP", "RH", "UGRD", "VGRD", "HGT")


# ---------------------------------------------------------------------------
# Small physics
# ---------------------------------------------------------------------------

def saturation_vapour_pressure(T):
    """Pa, Bolton (1980)."""
    return 611.2 * np.exp(17.67 * (T - 273.15) / (T - 29.65))


def virtual_temperature(T, RH, p):
    e = np.clip(RH, 0.0, 100.0) / 100.0 * saturation_vapour_pressure(T)
    q = EPS * e / (p - (1.0 - EPS) * e)
    return T * (1.0 + 0.608 * q)


def standard_temperature(p_pa):
    """US standard atmosphere, troposphere + isothermal lower stratosphere."""
    p = np.asarray(p_pa, float)
    t = 288.15 * (p / 101_325.0) ** (RD * LAPSE / G0)
    return np.maximum(t, 216.65)


def is_withheld(station, source):
    """Fixed one-in-five of ASOS stations, the same set in every run."""
    return source == "asos" and zlib.crc32(station.encode()) % 5 == 0


# ---------------------------------------------------------------------------
# Soundings onto the analysis levels
# ---------------------------------------------------------------------------

def superob_soundings(obs, levels_pa, max_gap_pa=15_000.0):
    """
    {var: [(level_index, station, lat, lon, value, error_std), ...]}

    Each sounding is interpolated in log(p) to each analysis level it
    brackets. A level is skipped when the nearest reported levels either
    side are more than `max_gap_pa` apart -- a profile with a hole is not
    evidence about what is in the hole.
    """
    prof = defaultdict(list)
    for o in obs:
        if o.source == "raob" and o.pressure:
            prof[(o.station, o.variable)].append(o)
    out = defaultdict(list)
    for (stn, var), rows in prof.items():
        rows = sorted(rows, key=lambda o: o.pressure)
        p = np.array([o.pressure for o in rows])
        v = np.array([o.value for o in rows])
        p, keep = np.unique(p, return_index=True)
        v = v[keep]
        if p.size < 2:
            continue
        lp = np.log(p)
        for k, pk in enumerate(levels_pa):
            if pk < p[0] or pk > p[-1]:
                continue
            j = int(np.searchsorted(p, pk))
            if j < p.size and p[j] == pk:
                val = v[j]
            else:
                if p[j] - p[j - 1] > max_gap_pa:
                    continue
                w = (np.log(pk) - lp[j - 1]) / (lp[j] - lp[j - 1])
                val = (1 - w) * v[j - 1] + w * v[j]
            o = rows[0]
            out[var].append((k, stn, o.lat, o.lon, float(val), o.error_std))
    return out


# ---------------------------------------------------------------------------
# First guess
# ---------------------------------------------------------------------------

def first_guess_from_soundings(so, levels_pa, shape):
    """Horizontally uniform: the mean sounding, level by level."""
    L = len(levels_pa)
    fallback = {"TMP": standard_temperature(np.asarray(levels_pa)),
                "RH": np.full(L, 50.0), "UGRD": np.zeros(L), "VGRD": np.zeros(L)}
    bg, used = {}, {}
    for var in ("TMP", "RH", "UGRD", "VGRD"):
        prof = fallback[var].astype(float).copy()
        n = np.zeros(L, int)
        for k in range(L):
            vals = [r[4] for r in so.get(var, []) if r[0] == k]
            if vals:
                prof[k] = float(np.mean(vals))
                n[k] = len(vals)
        bg[var] = np.broadcast_to(prof[:, None, None], (L,) + shape).copy()
        used[var] = n
    label = ("sounding_mean" if any(u.any() for u in used.values())
             else "standard_atmosphere")
    return bg, label


def sigma_to_pressure_levels(field, p3d, levels_pa, kind):
    """
    Column interpolation of a sigma-level field onto pressure levels, in
    log(p). Below the ground (p > surface) temperature follows the standard
    lapse rate from the lowest level and winds hold; above the top model
    level everything holds.
    """
    L = len(levels_pa)
    nz, ny, nx = field.shape
    out = np.empty((L, ny, nx))
    lp = np.log(p3d)
    for k, pk in enumerate(levels_pa):
        lk = np.log(pk)
        j = (p3d < pk).sum(axis=0)                  # levels above pk
        below = j >= nz
        above = j == 0
        jj = np.clip(j, 1, nz - 1)
        a = np.take_along_axis(lp, (jj - 1)[None], 0)[0]
        b = np.take_along_axis(lp, jj[None], 0)[0]
        fa = np.take_along_axis(field, (jj - 1)[None], 0)[0]
        fb = np.take_along_axis(field, jj[None], 0)[0]
        w = np.clip((lk - a) / (b - a), 0.0, 1.0)
        v = (1 - w) * fa + w * fb
        v = np.where(above, field[0], v)
        if kind == "TMP":
            ext = field[-1] * (pk / p3d[-1]) ** (RD * LAPSE / G0)
            v = np.where(below, ext, v)
        else:
            v = np.where(below, field[-1], v)
        out[k] = v
    return out


def first_guess_from_forecast(path, lead_hours, levels_pa, rh_default=None):
    """
    The previous run's forecast at `lead_hours`, on pressure levels.

    Returns None -- never a substitute -- if the file lacks that lead time
    (the run was cut short by its deadline, or died): the caller then falls
    back and records that it did.
    """
    z = np.load(path, allow_pickle=False)
    t = np.asarray(z["times_s"], float) / 3600.0
    hit = np.where(np.abs(t - lead_hours) < 1e-3)[0]
    if hit.size == 0:
        return None
    i = int(hit[0])
    sigma, p_top = np.asarray(z["sigma"], float), float(z["p_top"])
    pi = z["pi"][i].astype(float)
    p3d = p_top + sigma[:, None, None] * pi[None]
    theta = z["theta"][i].astype(float)
    T = theta * (p3d / 100_000.0) ** (RD / 1004.6)
    bg = {"TMP": sigma_to_pressure_levels(T, p3d, levels_pa, "TMP"),
          "UGRD": sigma_to_pressure_levels(z["u"][i].astype(float), p3d, levels_pa, "U"),
          "VGRD": sigma_to_pressure_levels(z["v"][i].astype(float), p3d, levels_pa, "V")}
    shape = bg["TMP"].shape
    # The model is dry: humidity is not carried, so its first guess is the
    # previous ANALYSIS's humidity if given, else a flat 50 %. It enters only
    # through virtual temperature, a <1 % effect on heights.
    bg["RH"] = (np.full(shape, 50.0) if rh_default is None
                else np.asarray(rh_default, float).reshape(shape))
    return bg


# ---------------------------------------------------------------------------
# Analysis steps
# ---------------------------------------------------------------------------

def analyse_upper(bg, so, lat, lon, domain, L_km=450.0, lam=0.25):
    """Level-by-level increments from the sounding superobs."""
    info = {}
    for var in ("TMP", "RH", "UGRD", "VGRD"):
        rows = so.get(var, [])
        per_level = {}
        for k in range(bg[var].shape[0]):
            r = [x for x in rows if x[0] == k]
            if not r:
                continue
            la = np.array([x[2] for x in r]); lo = np.array([x[3] for x in r])
            y = np.array([x[4] for x in r]); s = np.array([x[5] for x in r])
            d = y - geo.bilinear(bg[var][k], la, lo, domain, clip=True)
            bg[var][k] += barnes.barnes_increments(lat, lon, la, lo, d, s, L_km,
                                                   lam=lam, passes=2)
            per_level[k] = len(r)
        info[var] = per_level
    bg["RH"] = np.clip(bg["RH"], 1.0, 100.0)
    return bg, info


def analyse_pmsl(obs, lat, lon, L_km=300.0):
    ob = [o for o in obs if o.variable == "PMSL"]
    if not ob:
        return None, 0
    la = np.array([o.lat for o in ob]); lo = np.array([o.lon for o in ob])
    y = np.array([o.value for o in ob]); s = np.array([o.error_std for o in ob])
    b = float(np.median(y))
    inc = barnes.barnes_increments(lat, lon, la, lo, y - b, s, L_km,
                                   lam=0.1, passes=3)
    return b + inc, len(ob)


def hydrostatic_heights(T, RH, pmsl, levels_pa):
    """Heights (m) on every level, from sea-level pressure upward."""
    p = np.asarray(levels_pa, float)[:, None, None]
    Tv = virtual_temperature(T, RH, p)
    z = np.empty_like(T)
    # Sea level to the first level: mean Tv of the layer, one fixed-point step.
    tv_mean = Tv[0].copy()
    for _ in range(2):
        z0 = RD * tv_mean / G0 * np.log(pmsl / p[0, 0, 0])
        tv_mean = Tv[0] + 0.5 * LAPSE * z0
    z[0] = z0
    for k in range(1, T.shape[0]):
        z[k] = z[k - 1] + RD / G0 * 0.5 * (Tv[k - 1] + Tv[k]) * np.log(p[k - 1, 0, 0] / p[k, 0, 0])
    return z


def value_at_height(F, Z, zt, kind):
    """Column value at height zt (2-D); z increases with level index."""
    L = F.shape[0]
    j = (Z < zt[None]).sum(axis=0)
    low = j == 0
    high = j >= L
    jj = np.clip(j, 1, L - 1)
    za = np.take_along_axis(Z, (jj - 1)[None], 0)[0]
    zb = np.take_along_axis(Z, jj[None], 0)[0]
    fa = np.take_along_axis(F, (jj - 1)[None], 0)[0]
    fb = np.take_along_axis(F, jj[None], 0)[0]
    w = np.clip((zt - za) / (zb - za), 0.0, 1.0)
    v = (1 - w) * fa + w * fb
    if kind == "TMP":
        v = np.where(low, F[0] + LAPSE * (Z[0] - zt), v)
    else:
        v = np.where(low, F[0], v)
    return np.where(high, F[-1], v)


def blend_surface(fields, Z, terrain, obs, lat, lon, domain,
                  blend_depth_m=1000.0, L_km=120.0, max_elev_mismatch_m=600.0):
    """Analyse surface innovations and apply them to the lowest levels."""
    info = {}
    w = np.clip(1.0 - (Z - terrain[None]) / blend_depth_m, 0.0, 1.0)
    for var, dz in (("TMP", 2.0), ("RH", 2.0), ("UGRD", 10.0), ("VGRD", 10.0)):
        ob = [o for o in obs if o.variable == var and o.pressure is None]
        if not ob:
            info[var] = 0
            continue
        la = np.array([o.lat for o in ob]); lo = np.array([o.lon for o in ob])
        y = np.array([o.value for o in ob]); s = np.array([o.error_std for o in ob])
        elev = np.array([o.elevation if o.elevation is not None else np.nan
                         for o in ob])
        z_grid = geo.bilinear(terrain, la, lo, domain, clip=True)
        mismatch = np.where(np.isfinite(elev), elev - z_grid, 0.0)
        if var == "TMP":
            # Bring the station's reading to the grid terrain height.
            y = y + LAPSE * mismatch
        use = np.abs(mismatch) <= max_elev_mismatch_m
        sfc = value_at_height(fields[var], Z, terrain + dz, var)
        d = y - geo.bilinear(sfc, la, lo, domain, clip=True)
        d = np.where(use, d, np.nan)
        inc = barnes.barnes_increments(lat, lon, la, lo, d, s, L_km,
                                       lam=0.25, passes=2)
        fields[var] = fields[var] + w * inc[None]
        info[var] = int(np.isfinite(d).sum())
    fields["RH"] = np.clip(fields["RH"], 1.0, 100.0)
    return fields, info


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def first_guess_from_analysis(features):
    """A previous ANALYSIS [C, L, Y, X] as a first guess (winds, T, RH)."""
    f = np.asarray(features, float)
    return {c: f[CHANNELS.index(c)].copy() for c in ("TMP", "RH", "UGRD", "VGRD")}


def build_analysis(cycle, observations, terrain, domain, levels_hpa,
                   previous_forecast=None, previous_rh=None,
                   previous_analysis=None):
    """
    Returns (features [C, L, Y, X] float32, meta dict).

    `observations` must already have passed QC; this function only removes
    the withheld stations and anything valid after the cycle time.
    """
    levels_pa = np.asarray(levels_hpa, float) * 100.0
    ny, nx = terrain.shape
    lat, lon = geo.cell_centres(domain, ny, nx)

    late = [o for o in observations if o.time > cycle]
    if late:
        raise ValueError(f"{len(late)} observations are valid after the cycle "
                         f"time; nothing observed later may enter the run")
    withheld = sorted({o.station for o in observations
                       if is_withheld(o.station, o.source)})
    used = [o for o in observations if not is_withheld(o.station, o.source)]

    so = superob_soundings(used, levels_pa)

    bg, bg_label = None, None
    if previous_forecast is not None:
        path, lead = previous_forecast
        bg = first_guess_from_forecast(path, lead, levels_pa, previous_rh)
        if bg is not None:
            bg_label = f"previous_forecast:{path}:+{lead:g}h"
    # No usable previous forecast (the last run died or was cut short before
    # +6 h), and no soundings at this cycle (06Z, 18Z): the previous run's
    # ANALYSIS is 6 h old but built from real soundings, which is far closer
    # to the atmosphere than a standard one. Found 2026-09-21 18Z, when the
    # 12Z run had died at 3.75 h and the upper air fell to a standard
    # atmosphere -- no jet, no gradients (P-56).
    if bg is None and not so.get("TMP") and previous_analysis is not None:
        path, feats_prev = previous_analysis
        bg = first_guess_from_analysis(feats_prev)
        bg_label = f"previous_analysis:{path} (6 h old; previous forecast unusable)"
    if bg is None:
        bg, bg_label = first_guess_from_soundings(so, levels_pa, (ny, nx))
        if previous_forecast is not None:
            bg_label += " (previous forecast lacked the lead time)"

    fields, upper_info = analyse_upper(bg, so, lat, lon, domain)

    pmsl, n_p = analyse_pmsl(used, lat, lon)
    if pmsl is None:
        raise ValueError("no sea-level pressure observations: heights cannot "
                         "be anchored, and a run without an anchor is not "
                         "an analysis")
    Z = hydrostatic_heights(fields["TMP"], fields["RH"], pmsl, levels_pa)
    fields, sfc_info = blend_surface(fields, Z, terrain, used, lat, lon, domain)
    fields["HGT"] = hydrostatic_heights(fields["TMP"], fields["RH"], pmsl,
                                        levels_pa)

    feats = np.stack([fields[c] for c in CHANNELS]).astype(np.float32)
    stations = defaultdict(set)
    for o in used:
        stations[o.source].add(o.station)
    meta = {
        "valid_time": cycle.isoformat(),
        "run_time": cycle.isoformat(),
        "channels": np.array(CHANNELS),
        "levels_hPa": np.asarray(levels_hpa),
        "lat": lat.astype(np.float32),
        "lon": lon.astype(np.float32),
        "p_msl": pmsl.astype(np.float32),
        "source": "observations",
        "first_guess": bg_label,
        "provenance": json.dumps({
            "first_guess": bg_label,
            "soundings_per_level": {v: {str(k): n for k, n in d.items()}
                                    for v, d in upper_info.items()},
            "surface_obs_used": sfc_info,
            "pmsl_obs": n_p,
            "platforms_used": {s: sorted(v) for s, v in stations.items()},
            "withheld": withheld,
        }),
    }
    return feats, meta

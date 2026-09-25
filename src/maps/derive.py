"""
Map diagnostics from a saved forecast (u, v, theta, pi on sigma levels).

Everything a map product needs is derived here from the model's own
variables with the model's own constants (dynamics/sigma.py). No field is
taken from any other source.

  * pressure      p = p_top + sigma * pi, with pi = p_s - p_top
  * temperature   T = theta (p / P0)^kappa
  * height        hydrostatic integration upward from the terrain, exactly
                  as sigma.hydrostatic_geopotential does it
  * pressure levels  linear in ln p between model levels. Below the ground,
                  T follows the standard lapse rate down from the lowest
                  level and heights follow hypsometrically, so 1000 hPa
                  exists (for thickness) over high ground, as on every
                  operational chart. Such points are flagged in `below`.
  * MSLP          the surface pressure reduced along a standard-lapse-rate
                  column from the ground temperature
  * winds         de-staggered from the C-grid faces to the cell centres
                  (u sits on western faces, v on southern faces)
  * vorticity     relative vorticity from centred differences on the
                  model's dx, dy, plus f = 2 Omega sin(lat)

"Near-surface" fields are the LOWEST MODEL LEVEL, a few hundred metres
above the ground. The model has no 2 m or 10 m diagnosis, and its maps say
so instead of pretending.
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "src", ROOT / "src" / "dynamics", ROOT / "src" / "analysis"):
    if str(p) not in sys.path:
        sys.path.insert(0, str(p))

from sigma import RD, KAPPA, P0, G0            # noqa: E402

LAPSE = 0.0065                                 # K/m, standard atmosphere
OMEGA = 7.292e-5
PRESSURE_LEVELS_HPA = (1000, 850, 700, 500, 250)


def load_forecast(path, domain=None):
    """A forecast.npz as a plain dict, with lat/lon of the cell centres."""
    z = np.load(path, allow_pickle=True)
    f = {k: z[k] for k in z.files}
    f["path"] = str(path)
    ny, nx = f["theta"].shape[-2:]
    if domain is None:
        import config
        domain = config.DOMAIN
    from geo import cell_centres, spacing_m
    f["lat"], f["lon"] = cell_centres(domain, ny, nx)
    f["dy"], f["dx"] = spacing_m(domain, ny, nx)
    f["hours"] = np.asarray(f["times_s"], dtype=float) / 3600.0
    f["domain"] = domain
    return f


def destagger(u, v):
    """C-grid face winds -> cell-centre winds (edges replicate)."""
    uc = 0.5 * (u + np.concatenate([u[..., 1:], u[..., -1:]], axis=-1))
    vc = 0.5 * (v + np.concatenate([v[..., 1:, :], v[..., -1:, :]], axis=-2))
    return uc, vc


def column_state(theta, pi, sigma, p_top, terrain):
    """p, T and height (m) on model levels. Index 0 is the lid, -1 the ground."""
    sig = np.asarray(sigma, dtype=float).reshape(-1, 1, 1)
    p = float(p_top) + sig * pi[None]
    T = theta * (p / P0) ** KAPPA
    ps = float(p_top) + pi
    phi = np.empty_like(theta, dtype=float)
    phi[-1] = G0 * terrain + RD * T[-1] * np.log(ps / p[-1])
    for k in range(theta.shape[0] - 2, -1, -1):
        phi[k] = phi[k + 1] + RD * 0.5 * (T[k] + T[k + 1]) * np.log(p[k + 1] / p[k])
    return p, T, phi / G0, ps


def ground_temperature(T, Z, terrain):
    """T at the ground, extrapolated from the lowest level at the standard lapse."""
    return T[-1] + LAPSE * (Z[-1] - terrain)


def mslp(ps, T, Z, terrain):
    """Surface pressure reduced to sea level along a standard-lapse column."""
    Tg = ground_temperature(T, Z, terrain)
    return ps * ((Tg + LAPSE * terrain) / Tg) ** (G0 / (RD * LAPSE))


def to_pressure(field, p, pk, kind, ps=None, T=None, Z=None, terrain=None):
    """
    Interpolate a model-level field to pressure pk (Pa), linear in ln p.

    Above the lid the top value holds. Below the lowest model level:
    kind "T" extrapolates at the standard lapse rate (in height), kind "Z"
    hypsometrically from the ground, and anything else (winds) holds the
    lowest level. Returns (values, below_ground mask).
    """
    nz = p.shape[0]
    lp, lk = np.log(p), np.log(pk)
    j = (p < pk).sum(axis=0)                 # model levels above pk
    jj = np.clip(j, 1, nz - 1)
    a = np.take_along_axis(lp, (jj - 1)[None], 0)[0]
    b = np.take_along_axis(lp, jj[None], 0)[0]
    fa = np.take_along_axis(field, (jj - 1)[None], 0)[0]
    fb = np.take_along_axis(field, jj[None], 0)[0]
    w = np.clip((lk - a) / (b - a), 0.0, 1.0)
    out = (1 - w) * fa + w * fb
    out = np.where(j == 0, field[0], out)
    low = j >= nz                             # pk below the lowest level
    below = pk > ps if ps is not None else low
    if low.any():
        if kind == "T":
            # T at pk along the standard lapse from the lowest level:
            # T = T_k (pk / p_k)^(R LAPSE / g)
            ext = field[-1] * (pk / p[-1]) ** (RD * LAPSE / G0)
            out = np.where(low, ext, out)
        elif kind == "Z":
            Tk = T[-1]
            ext = Z[-1] + Tk / LAPSE * (1.0 - (pk / p[-1]) ** (RD * LAPSE / G0))
            out = np.where(low, ext, out)
        else:
            out = np.where(low, field[-1], out)
    return out, below


def relative_vorticity(u, v, dx, dy):
    """dv/dx - du/dy on cell centres (centred, one-sided at the edges)."""
    return np.gradient(v, dx, axis=-1) - np.gradient(u, dy, axis=-2)


def coriolis(lat):
    return 2.0 * OMEGA * np.sin(np.radians(lat))


def snapshot(f, i, levels_hpa=PRESSURE_LEVELS_HPA):
    """
    Every diagnostic for snapshot i, as a dict of 2-D arrays (SI units:
    K, Pa, m, m/s, 1/s). Keys: 'ps', 'mslp', 'terrain', 'T_low', 'u_low',
    'v_low', 'z_low_agl', and per level L: 'T<L>', 'Z<L>', 'u<L>', 'v<L>',
    'below<L>', 'avort<L>'; plus 'thick_1000_500'.
    """
    theta = np.asarray(f["theta"][i], dtype=float)
    pi = np.asarray(f["pi"][i], dtype=float)
    terrain = np.asarray(f["terrain"], dtype=float)
    uc, vc = destagger(np.asarray(f["u"][i], float), np.asarray(f["v"][i], float))
    p, T, Z, ps = column_state(theta, pi, f["sigma"], f["p_top"], terrain)
    d = {"ps": ps, "terrain": terrain, "mslp": mslp(ps, T, Z, terrain),
         "T_low": T[-1], "u_low": uc[-1], "v_low": vc[-1],
         "z_low_agl": Z[-1] - terrain, "hour": float(f["hours"][i])}
    fcor = coriolis(f["lat"])
    for L in levels_hpa:
        pk = L * 100.0
        d[f"T{L}"], d[f"below{L}"] = to_pressure(T, p, pk, "T", ps=ps)
        d[f"Z{L}"], _ = to_pressure(Z, p, pk, "Z", ps=ps, T=T, Z=Z,
                                    terrain=terrain)
        d[f"u{L}"], _ = to_pressure(uc, p, pk, "wind", ps=ps)
        d[f"v{L}"], _ = to_pressure(vc, p, pk, "wind", ps=ps)
        d[f"avort{L}"] = relative_vorticity(d[f"u{L}"], d[f"v{L}"],
                                            f["dx"], f["dy"]) + fcor
    if 1000 in levels_hpa and 500 in levels_hpa:
        d["thick_1000_500"] = d["Z500"] - d["Z1000"]
    return d


def load_analysis(path):
    """obs_analysis_f00.npz as a dict (features [C, L, Y, X], levels, p_msl)."""
    z = np.load(path, allow_pickle=True)
    a = {k: z[k] for k in z.files}
    a["channels"] = [str(c) for c in a["channels"]]
    a["levels_hPa"] = np.asarray(a["levels_hPa"], dtype=float)
    return a


def _interp_levels(F, levels_pa, pk):
    """Column values at pk from pressure-level fields F [L, Y, X], linear in ln p."""
    lp = np.log(np.asarray(levels_pa, dtype=float))
    order = np.argsort(lp)
    lp, F = lp[order], F[order]
    lk = np.log(pk)
    j = int(np.clip(np.searchsorted(lp, lk), 1, len(lp) - 1))
    w = np.clip((lk - lp[j - 1]) / (lp[j] - lp[j - 1]), 0.0, 1.0)
    return (1 - w) * F[j - 1] + w * F[j]


def analysis_snapshot(a, terrain, dx, dy, lat, z_low_agl=None,
                      levels_hpa=PRESSURE_LEVELS_HPA):
    """
    The same keys as snapshot(), from the ANALYSIS on pressure levels, for
    hour 0 (a forecast's snapshots start at its first output hour).

    The near-surface fields are taken at `z_low_agl` above the ground (the
    forecast's lowest-level height), so hour 0 is comparable with hour 1.
    """
    import build                                    # analysis/build.py
    feats = np.asarray(a["features"], dtype=float)
    ch = a["channels"]
    T, U, V, Zf = (feats[ch.index(c)] for c in ("TMP", "UGRD", "VGRD", "HGT"))
    pa = a["levels_hPa"] * 100.0
    zt = terrain + (np.nanmean(z_low_agl) if z_low_agl is not None else 2.0)
    d = {"terrain": terrain, "mslp": np.asarray(a["p_msl"], dtype=float),
         "T_low": build.value_at_height(T, Zf, zt, "TMP"),
         "u_low": build.value_at_height(U, Zf, zt, "UGRD"),
         "v_low": build.value_at_height(V, Zf, zt, "VGRD"),
         "z_low_agl": zt - terrain, "hour": 0.0,
         "T_2m": build.value_at_height(T, Zf, terrain + 2.0, "TMP")}
    fcor = coriolis(lat)
    for L in levels_hpa:
        pk = L * 100.0
        for key, F in (("T", T), ("Z", Zf), ("u", U), ("v", V)):
            d[f"{key}{L}"] = _interp_levels(F, pa, pk)
        d[f"below{L}"] = d[f"Z{L}"] < terrain
        d[f"avort{L}"] = relative_vorticity(d[f"u{L}"], d[f"v{L}"], dx, dy) + fcor
    if 1000 in levels_hpa and 500 in levels_hpa:
        d["thick_1000_500"] = d["Z500"] - d["Z1000"]
    return d

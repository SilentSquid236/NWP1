"""
Sounding indices for every model column at once (2026-10-20, prompt 232).

The viewer's click sounding shows these numbers; they are computed here, once,
for every column and hour, and checked against SHARPpy (Blumberg et al. 2017)
on sampled model columns (tools/check_indices.py). One implementation, vectorised
over columns: arrays are (levels, columns) with level 0 at the SURFACE and the
levels rising.

Parcels (as SHARPpy defines them):
  SB  the surface parcel (2 m temperature and dewpoint, surface pressure);
  ML  the mean potential temperature and mixing ratio of the lowest 100 hPa,
      lifted from the surface;
  MU  the level of highest equivalent potential temperature in the lowest 300 hPa.
Lifting follows SHARPpy (sharptab/thermo.py, BSD-3; from NSHARP, Hart and Korotky
1991): the LCL temperature from its `lcltemp` polynomial and the LCL pressure from
Poisson's equation, then the saturated ascent from the Wobus function: one
conserved `thetam` per parcel, and `satlift` solved by secant iteration at each
level (ported here to arrays). Buoyancy uses virtual temperature for both parcel
and environment (Doswell and Rasmussen 1994). As in SHARPpy's parcelx: CAPE sums
the positive layer areas above the LCL up to the EL; CIN sums the negative areas
from the parcel's origin to the LFC that lie below 500 hPa; CAPE 0 gives CIN 0.

THE MODEL LID IS AT 200 hPa. A parcel still buoyant there has its EL set to the
lid, its CAPE truncated, and el_at_lid = 1.

Kinematics: bulk shear from the 10 m wind; the Bunkers et al. (2000) right-moving
storm motion by SHARPpy's non-parcel method (winds.non_parcel_bunkers_motion: the
0-6 km non-pressure-weighted mean wind, 7.5 m/s to the right of the 0-6 km shear
vector); storm-relative helicity for that motion over 0-1 and 0-3 km. Lapse rates
use virtual temperature, as SHARPpy's params.lapse_rate does.
"""
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "dynamics"))
from moisture import es_bolton, dewpoint  # noqa: E402

RD, CP, LV, G, EPS = 287.05, 1005.7, 2.501e6, 9.80665, 0.622
KAP = RD / CP
KT = 1.0 / 0.514444


def mixr_sat(T, p):
    es = es_bolton(T)
    return EPS * es / np.maximum(p - es, 1.0)


def virt(T, r):
    return T * (1.0 + r / EPS) / (1.0 + r)


def lcl(T, Td, p):
    """(T_LCL, p_LCL), K and Pa: SHARPpy thermo.lcltemp and thalvl."""
    Td = np.minimum(Td, T)
    tc, sd = T - 273.15, T - Td
    dlt = sd * (1.2185 + 0.001278 * tc + sd * (-0.00219 + 1.173e-5 * sd - 0.0000052 * tc))
    TL = T - dlt
    return TL, p * (TL / T) ** (1.0 / KAP_SHARPPY)


def theta_e(T, Td, p):
    """Equivalent potential temperature, Bolton (1980) eq. 43."""
    r = EPS * es_bolton(Td) / np.maximum(p - es_bolton(Td), 1.0)
    TL, _ = lcl(T, Td, p)
    th = T * (1e5 / p) ** (0.2854 * (1.0 - 0.28 * r))
    return th * np.exp((3376.0 / TL - 2.54) * r * (1.0 + 0.81 * r))


def wobf(t):
    """The Wobus function (deg C in, deg C out), SHARPpy thermo.wobf, for arrays."""
    t = np.asarray(t, float) - 20.0
    npol = 1.0 + t * (-8.841660499999999e-3 + t * (1.4714143e-4 + t * (-9.671989000000001e-7
                       + t * (-3.2607217e-8 + t * (-3.8598073e-10)))))
    npol = 15.13 / np.power(npol, 4)
    ppol = t * (4.9618922e-07 + t * (-6.1059365e-09 + t * (3.9401551e-11 + t * (-1.2588129e-13
                                                                             + t * (1.6688280e-16)))))
    ppol = 1.0 + t * (3.6182989e-03 + t * (-1.3603273e-05 + ppol))
    ppol = (29.93 / np.power(ppol, 4)) + (0.96 * t) - 14.8
    return np.where(t <= 0, npol, ppol)


def satlift(p_hpa, thetam, conv=0.1, it_max=50):
    """Temperature (deg C) on the saturated adiabat thetam (deg C) at p_hpa: SHARPpy
    thermo.satlift, the secant iteration done for all points at once."""
    p_hpa, thetam = np.broadcast_arrays(np.asarray(p_hpa, float), np.asarray(thetam, float))
    pwrp = np.power(p_hpa / 1000.0, KAP_SHARPPY)
    t1 = (thetam + 273.15) * pwrp - 273.15
    e1 = wobf(t1) - wobf(thetam)
    rate = np.ones_like(t1)
    t2 = t1 - e1 * rate
    e2 = (t2 + 273.15) / pwrp - 273.15
    e2 = e2 + wobf(t2) - wobf(e2) - thetam
    eor = e2 * rate
    for _ in range(it_max):
        if np.nanmax(np.abs(eor)) <= conv:
            break
        rate = (t2 - t1) / np.where(e2 - e1 == 0, 1e-12, e2 - e1)
        t1, e1 = t2, e2
        t2 = t1 - e1 * rate
        e2 = (t2 + 273.15) / pwrp - 273.15
        e2 = e2 + wobf(t2) - wobf(e2) - thetam
        eor = e2 * rate
    return np.where(np.abs(p_hpa - 1000.0) <= 0.001, thetam, t2 - eor)


KAP_SHARPPY = 0.28571426      # SHARPpy's ROCP


def sharppy_vappres(t_c):
    """Vapour pressure (hPa) over water at t_c (deg C): SHARPpy thermo.vappres."""
    t = np.asarray(t_c, float)
    pol = t * (1.1112018e-17 + (t * -3.0994571e-20))
    pol = t * (2.1874425e-13 + (t * (-1.789232e-15 + pol)))
    pol = t * (4.3884180e-09 + (t * (-2.988388e-11 + pol)))
    pol = t * (7.8736169e-05 + (t * (-6.111796e-07 + pol)))
    pol = 0.99999683 + (t * (-9.082695e-03 + pol))
    return 6.1078 / pol ** 8


def sharppy_mixratio(p_hpa, t_c):
    """Mixing ratio (g/kg) at dewpoint t_c and p_hpa, with SHARPpy's enhancement factor."""
    x = 0.02 * (t_c - 12.5 + (7500.0 / p_hpa))
    wfw = 1.0 + (0.0000045 * p_hpa) + (0.0014 * x * x)
    fwesw = wfw * sharppy_vappres(t_c)
    return 621.97 * (fwesw / (p_hpa - fwesw))


def sharppy_temp_at_mixrat(w, p_hpa):
    """Temperature (deg C) at mixing ratio w (g/kg) and p_hpa: SHARPpy thermo.temp_at_mixrat."""
    c1, c2, c3, c4, c5, c6 = 0.0498646455, 2.4082965, 7.07475, 38.9114, 0.0915, 1.2035
    x = np.log10(w * p_hpa / (622.0 + w))
    return (np.power(10.0, (c1 * x) + c2) - c3 + (c4 * np.power(np.power(10.0, c5 * x) - c6, 2))) - 273.15


def lift(pe, T0, r0, p0, Td0=None):
    """Parcel virtual temperature at every level pe (L, N) from (T0, r0, p0) (N,).
    Levels below p0 (higher pressure) are NaN. Also returns (T_LCL, p_LCL).
    Td0 (K) overrides the dewpoint computed from r0 (the ML parcel, as SHARPpy)."""
    Td0 = dewpoint(r0, p0) if Td0 is None else Td0
    TL, pL = lcl(T0, Td0, p0)
    thta = (TL - 273.15 + 273.15) * (1e5 / pL) ** KAP_SHARPPY - 273.15      # theta at 1000 hPa, deg C
    thetam = thta - wobf(thta) + wobf(TL - 273.15)
    above = pe <= p0[None] + 1e-6
    dry = above & (pe >= pL[None])
    wet = above & (pe < pL[None])
    Tdry = T0[None] * (pe / p0[None]) ** KAP_SHARPPY
    Twet = satlift(pe / 100.0, np.broadcast_to(thetam, pe.shape)) + 273.15
    Tv = np.where(dry, virt(Tdry, r0[None]), np.nan)
    Tv = np.where(wet, virt(Twet, mixr_sat(Twet, pe)), Tv)
    return Tv, TL, pL


def _interp_lnp(pe, a, pq):
    """a (L, N) at pressure pq (N,), linear in ln p; NaN outside."""
    lp, lq = np.log(pe), np.log(pq)
    out = np.full(pq.shape, np.nan)
    for k in range(pe.shape[0] - 1):
        m = (lp[k] >= lq) & (lq >= lp[k + 1]) & np.isnan(out)
        w = (lp[k] - lq) / np.where(lp[k] - lp[k + 1] == 0, 1, lp[k] - lp[k + 1])
        out = np.where(m, a[k] + w * (a[k + 1] - a[k]), out)
    return out


def buoyancy_areas(pe, dTv, p_origin, p_lcl, ze):
    """CAPE, CIN (J/kg), LFC and EL pressure (Pa), LFC height (m, same datum as ze),
    el_at_lid (0/1). dTv = parcel minus environment virtual temperature (L, N)."""
    L, N = pe.shape
    cape, cin = np.zeros(N), np.zeros(N)
    lfc, el = np.full(N, np.nan), np.full(N, np.nan)
    # LFC: lowest level at or above the LCL where dTv > 0 and the layer below is <= 0
    pos = (dTv > 0) & (pe <= p_lcl[None] + 1e-6)
    for k in range(L):
        m = np.isnan(lfc) & pos[k]
        if not m.any():
            continue
        if k == 0:
            lfc = np.where(m, pe[0], lfc)
            continue
        a, b = dTv[k - 1], dTv[k]
        # crossing between k-1 and k (linear in ln p); if level k-1 is below the LCL, not before the LCL
        w = np.where(np.isfinite(a) & (a <= 0) & (b > a), -a / np.where(b - a == 0, 1, b - a), 0.0)
        pc = np.exp(np.log(pe[k - 1]) + w * (np.log(pe[k]) - np.log(pe[k - 1])))
        # a > 0: level k-1 is below the LCL (else it would have been chosen) and already
        # buoyant, so free convection starts at the LCL itself
        pc = np.where(np.isfinite(a) & (a <= 0), pc, p_lcl)
        lfc = np.where(m, np.minimum(pc, p_lcl), lfc)
    # EL: the highest crossing from positive to non-positive above the LFC; buoyant at the lid -> lid
    top = pe[-1]
    el_lid = np.isfinite(lfc) & (dTv[-1] > 0)
    for k in range(L - 1, 0, -1):
        a, b = dTv[k - 1], dTv[k]
        m = np.isnan(el) & np.isfinite(lfc) & (a > 0) & (b <= 0) & (pe[k] < lfc)
        w = a / np.where(a - b == 0, 1, a - b)
        pc = np.exp(np.log(pe[k - 1]) + w * (np.log(pe[k]) - np.log(pe[k - 1])))
        el = np.where(m, pc, el)
    el = np.where(el_lid, top, el)
    el = np.where(np.isfinite(lfc) & np.isnan(el), lfc, el)
    # areas, layer by layer, with each layer split at its zero crossing
    for k in range(L - 1):
        a, b = dTv[k], dTv[k + 1]
        p1, p2 = pe[k], pe[k + 1]
        ok = np.isfinite(a) & np.isfinite(b)
        lp1, lp2 = np.log(p1), np.log(p2)
        cross = ok & (a * b < 0)
        w = np.where(cross, a / np.where(a - b == 0, 1, a - b), 0.0)
        lpc = lp1 + w * (lp2 - lp1)
        # piecewise parts: [p1, pc] with mean a/2, [pc, p2] with mean b/2; no crossing: (a+b)/2 over all
        part1 = np.where(cross, RD * 0.5 * a * (lp1 - lpc), RD * 0.5 * (a + b) * (lp1 - lp2))
        part2 = np.where(cross, RD * 0.5 * b * (lpc - lp2), 0.0)
        mid1 = np.exp(0.5 * (lp1 + np.where(cross, lpc, lp2)))
        mid2 = np.exp(0.5 * (lpc + lp2))
        for part, mid in ((part1, mid1), (part2, mid2)):
            inside = ok & np.isfinite(lfc) & (mid <= p_lcl) & (mid >= el)
            below = ok & (mid > np.where(np.isfinite(lfc), lfc, -1)) & (mid <= p_origin) & (mid > 50000.0)
            cape += np.where(inside & (part > 0), part, 0.0)
            cin += np.where(below & (part < 0), part, 0.0)
    cin = np.where(np.isfinite(lfc) & (np.floor(cape) > 0), cin, 0.0)
    z_lfc = _interp_lnp(pe, ze, np.where(np.isfinite(lfc), lfc, pe[0]))
    z_lfc = np.where(np.isfinite(lfc), z_lfc, np.nan)
    return cape, cin, lfc, el, z_lfc, el_lid.astype(float)


def _interp_lnp_z(z, p, zq):
    """Pressure at height zq (m AGL, scalar), interpolated in ln p."""
    return np.exp(_at_height(z, np.log(p), zq))


def _at_height(z, a, zq):
    """a (L, N) at height zq (scalar or (N,)), linear; NaN above the top."""
    zq = np.broadcast_to(np.asarray(zq, float), z.shape[1:])
    out = np.full(z.shape[1:], np.nan)
    for k in range(z.shape[0] - 1):
        m = (z[k] <= zq) & (zq <= z[k + 1]) & np.isnan(out)
        w = (zq - z[k]) / np.where(z[k + 1] - z[k] == 0, 1, z[k + 1] - z[k])
        out = np.where(m, a[k] + w * (a[k + 1] - a[k]), out)
    return out


def column_indices(p, T, r, z, u, v, terrain):
    """Every index for N columns. Inputs (L, N), level 0 = the surface level
    (2 m T, 2 m dewpoint as r, surface pressure, 10 m wind, z = terrain).
    Returns a dict of (N,) arrays. Units: J/kg, m AGL, hPa, K, mm, C/km, kt, m2/s2."""
    agl = z - terrain[None]
    tve = virt(T, r)
    out = {}
    # SB
    T0, r0, p0 = T[0], r[0], p[0]
    tvp, TL, pL = lift(p, T0, r0, p0)
    sb = buoyancy_areas(p, tvp - tve, p0, pL, agl)
    out.update(sbcape=sb[0], sbcin=sb[1], sblfc=sb[4], sbel=sb[3] / 100, el_at_lid=sb[5],
               sblcl=_interp_lnp(p, agl, pL))
    li = _interp_lnp(p, tve, np.full_like(p0, 50000.0)) - _interp_lnp(p, tvp, np.full_like(p0, 50000.0))
    out["li"] = li
    # ML: the lowest 100 hPa as SHARPpy defines it (params.mean_theta / mean_mixratio,
    # exact=True): plain means over [bottom value, every level inside the layer twice,
    # top value] of theta, and of dewpoint and pressure, the mixing ratio being that of
    # the mean dewpoint at the mean pressure. Lifted from the surface.
    top = p0 - 10000.0
    th = T * (1e5 / p) ** KAP_SHARPPY
    Tdp = dewpoint(r, p)
    inside = (p < p0[None]) & (p > top[None])
    inside[0] = False
    n_in = inside.sum(0)
    th_top = _interp_lnp(p, th, top); td_top = _interp_lnp(p, Tdp, top)
    thm = (th[0] + 2.0 * np.where(inside, th, 0.0).sum(0) + th_top) / (2.0 * n_in + 2.0)
    tdm = (Tdp[0] + 2.0 * np.where(inside, Tdp, 0.0).sum(0) + td_top) / (2.0 * n_in + 2.0)
    pm = (p0 + 2.0 * np.where(inside, p, 0.0).sum(0) + top) / (2.0 * n_in + 2.0)
    w_ml = sharppy_mixratio(pm / 100.0, tdm - 273.15)                  # g/kg
    rm_ = w_ml / 1000.0
    td_ml = sharppy_temp_at_mixrat(w_ml, p0 / 100.0) + 273.15
    tvp, TLm, pLm = lift(p, thm * (p0 / 1e5) ** KAP_SHARPPY, rm_, p0, Td0=td_ml)
    ml = buoyancy_areas(p, tvp - tve, p0, pLm, agl)
    out.update(mlcape=ml[0], mlcin=ml[1], mllcl=_interp_lnp(p, agl, pLm))
    # MU: highest theta-e in the lowest 300 hPa
    Td = dewpoint(r, p)
    the = np.where(p >= p0 - 30000.0, theta_e(T, Td, p), -np.inf)
    k0 = np.argmax(the, axis=0)
    cols = np.arange(p.shape[1])
    pmu, Tmu, rmu = p[k0, cols], T[k0, cols], r[k0, cols]
    tvp, TLu, pLu = lift(p, Tmu, rmu, pmu)
    mu = buoyancy_areas(p, tvp - tve, pmu, pLu, agl)
    out.update(mucape=mu[0], mucin=mu[1], mup=pmu / 100)
    # moisture and lapse rates
    dp = p[:-1] - p[1:]
    out["pw"] = np.sum(0.5 * (r[:-1] / (1 + r[:-1]) + r[1:] / (1 + r[1:])) * dp, axis=0) / G   # mm
    T3 = _at_height(agl, tve, 3000.0)
    out["lr03"] = (tve[0] - T3) / 3.0
    t7, t5 = _interp_lnp(p, tve, np.full_like(p0, 70000.0)), _interp_lnp(p, tve, np.full_like(p0, 50000.0))
    z7, z5 = _interp_lnp(p, agl, np.full_like(p0, 70000.0)), _interp_lnp(p, agl, np.full_like(p0, 50000.0))
    out["lr75"] = (t7 - t5) / ((z5 - z7) / 1000.0)
    # kinematics on heights
    hz = np.arange(0.0, 6001.0, 100.0)
    U = np.array([_at_height(agl, u, h) for h in hz]); V = np.array([_at_height(agl, v, h) for h in hz])
    U[0], V[0] = u[0], v[0]
    for name, h in (("shr1", 1000), ("shr6", 6000)):
        i = int(h / 100)
        out[name] = np.hypot(U[i] - U[0], V[i] - V[0]) * KT
    # SHARPpy non-parcel Bunkers: mean wind sampled evenly in pressure (every hPa in
    # SHARPpy; 120 samples here) from the surface to 6 km, shear = V(6 km) - V(sfc)
    p6 = _interp_lnp_z(agl, p, 6000.0)
    fr = np.linspace(0.0, 1.0, 121)[:, None]
    pq = p0[None] + fr * (p6 - p0)[None]          # evenly in p, as SHARPpy's dp = -1 hPa
    us = np.array([_interp_lnp(p, u, q) for q in pq]); vs = np.array([_interp_lnp(p, v, q) for q in pq])
    us[0], vs[0] = u[0], v[0]
    mu_, mv_ = np.nanmean(us, 0), np.nanmean(vs, 0)
    su, sv = U[60] - U[0], V[60] - V[0]
    mag = np.maximum(np.hypot(su, sv), 1e-6)
    cu, cv = mu_ + 7.5 * sv / mag, mv_ - 7.5 * su / mag
    out["rm_u"], out["rm_v"] = cu * KT, cv * KT
    for name, h in (("srh1", 1000), ("srh3", 3000)):
        n = int(h / 100)
        a = (U[1:n + 1] - cu) * (V[:n] - cv) - (U[:n] - cu) * (V[1:n + 1] - cv)
        out[name] = np.nansum(a, axis=0)
    return out

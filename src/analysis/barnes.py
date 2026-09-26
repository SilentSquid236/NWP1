"""
Successive-correction (Barnes) analysis of increments against a first guess.

    inc = barnes_increments(lat_g, lon_g, lat_o, lon_o, d, sigma_o, L_km)
    analysis = background + inc

WHY INCREMENTS AND A BACKGROUND WEIGHT

Classical Barnes analyses the observations themselves and normalises the
weights, so a point 800 km from the nearest sounding gets that sounding's
value at full strength. With a first guess (the previous run's forecast, or
the mean sounding) the right behaviour far from every observation is to KEEP
the first guess. So this analyses innovations d = y - H(x_b), and adds a
constant `lam` to the weight sum:

    inc(x) = sum_i w_i(x) d_i / (sum_i w_i(x) + lam)
    w_i(x) = exp(-r_i(x)^2 / L^2) * (s_ref / s_i)^2

Near a good observation the sum of weights is >> lam and the increment is the
observation's; far from all of them it goes to zero and the first guess
stands. `lam` is the ratio of observation to first-guess error variance in
disguise -- 0.25 says the first guess is trusted about half as much as an
observation of reference quality. The (s_ref/s_i)^2 factor weights each
observation by its error: a buoy temperature (1 K) counts more than an ASOS
one (1.5 K).

Each further pass analyses what the previous pass left unexplained AT THE
OBSERVATIONS, with a length scale multiplied by `gamma` -- the standard Barnes
refinement, which recovers smaller scales where the data are dense.

Distances are on a local tangent plane, which is accurate to well under 1 %
across a 2000 km box -- far below the representativeness error.
"""

import numpy as np

KM_PER_DEG_LAT = 111.132
KM_PER_DEG_LON = 111.320


def sq_distance_km(lat_a, lon_a, lat_b, lon_b):
    """Squared distance (km^2), every a against every b: shape (Na, Nb)."""
    lat_a = np.asarray(lat_a, float).ravel()[:, None]
    lon_a = np.asarray(lon_a, float).ravel()[:, None]
    lat_b = np.asarray(lat_b, float).ravel()[None, :]
    lon_b = np.asarray(lon_b, float).ravel()[None, :]
    c = np.cos(np.radians(0.5 * (lat_a + lat_b)))
    dy = (lat_a - lat_b) * KM_PER_DEG_LAT
    dx = (lon_a - lon_b) * KM_PER_DEG_LON * c
    return dx * dx + dy * dy


def _weights(r2, L_km, sigma_o, s_ref):
    return np.exp(-r2 / (L_km * L_km)) * (s_ref / np.asarray(sigma_o, float))[None, :] ** 2


def barnes_increments(lat_g, lon_g, lat_o, lon_o, d, sigma_o, L_km,
                      lam=0.25, passes=2, gamma=0.35, s_ref=None,
                      return_fit=False):
    """
    Increment field at grid points (shape of lat_g) from innovations d.

    With return_fit=True also returns the analysed increment AT the
    observations after the last pass -- what leave-one-out and the residual
    checks compare against.
    """
    lat_g = np.asarray(lat_g, float)
    d = np.asarray(d, float).ravel()
    sigma_o = np.broadcast_to(np.asarray(sigma_o, float), d.shape)
    ok = np.isfinite(d)
    shape = lat_g.shape
    if not ok.any():
        z = np.zeros(shape)
        return (z, np.zeros_like(d)) if return_fit else z
    lat_o = np.asarray(lat_o, float).ravel()[ok]
    lon_o = np.asarray(lon_o, float).ravel()[ok]
    d_ok, s_ok = d[ok], sigma_o[ok]
    s_ref = float(np.median(s_ok)) if s_ref is None else float(s_ref)

    r2_go = sq_distance_km(lat_g, lon_g, lat_o, lon_o)
    r2_oo = sq_distance_km(lat_o, lon_o, lat_o, lon_o)

    inc_g = np.zeros(r2_go.shape[0])
    inc_o = np.zeros(d_ok.size)
    L = float(L_km)
    for _ in range(max(1, int(passes))):
        resid = d_ok - inc_o
        wg = _weights(r2_go, L, s_ok, s_ref)
        wo = _weights(r2_oo, L, s_ok, s_ref)
        inc_g += (wg @ resid) / (wg.sum(axis=1) + lam)
        inc_o += (wo @ resid) / (wo.sum(axis=1) + lam)
        L *= gamma ** 0.5
    inc_g = inc_g.reshape(shape)
    if return_fit:
        fit = np.full(d.shape, np.nan)
        fit[ok] = inc_o
        return inc_g, fit
    return inc_g


def mean_spacing_km(lat_o, lon_o):
    """Mean nearest-neighbour distance between observations."""
    if len(lat_o) < 2:
        return np.inf
    r2 = sq_distance_km(lat_o, lon_o, lat_o, lon_o)
    np.fill_diagonal(r2, np.inf)
    return float(np.mean(np.sqrt(r2.min(axis=1))))

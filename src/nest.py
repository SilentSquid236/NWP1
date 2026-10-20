"""
Lateral boundaries for the 3 km run from the 12 km forecast (CAM stage S5).

The 3 km model starts from its own analysis, but its edges follow the 12 km
forecast of the same cycle, hour by hour, instead of being held at the start
(the user's choice, prompt 172). This module turns the 12 km forecast file
(forecast.py output: u, v, theta, pi on sigma levels, hourly) into boundary
frames on the 3 km grid and its own sigma levels:

1. Horizontal: bilinear, separately for each C-grid staggering. Both grids
   are regular in latitude and longitude over the same domain (geo.py), so
   a 3 km point's position in the 12 km grid is a fractional index:
   centres at (i + 1/2) cells, u points on the west faces (i cells),
   v points on the south faces (j cells). Outside the outermost 12 km points
   the edge value is used.
2. Terrain: the 12 km column sits over 12 km terrain h12; the 3 km column
   over h3. Surface pressure is moved hydrostatically from h12 to h3 with the
   temperature of the lowest 12 km level and a standard lapse rate.
3. Vertical: the 3 km sigma levels are placed in pressure with that surface
   pressure, and u, v and ln theta are interpolated linearly in ln p from the
   12 km column (ln theta is exactly linear in ln p for an isothermal layer). Below the lowest 12 km level and above the highest, the nearest
   level's value is kept (both grids have the same 200 hPa lid).

pi in the frames is p_s - p_top, as everywhere in NWP1.
"""
import numpy as np

from constants import RD, CP, KAPPA, P0, G0 as G   # one definition (was 287.04 / 1004.5)
LAPSE = 0.0065


def fractional_index(n_src, n_dst, dst_stagger, src_stagger=None, span=None):
    """Positions of destination points in SOURCE index units, one axis.

    Staggering "c": cell centres (i + 1/2 cells); "f": the low face of each
    cell (i cells). Both grids span the same interval of the domain, unless
    span = (a, b) gives the destination interval as fractions of the source
    one (test N1: a 3 km grid inside a larger 12 km grid; None = (0, 1)).
    src_stagger defaults to dst_stagger.
    """
    src_stagger = dst_stagger if src_stagger is None else src_stagger
    x = (np.arange(n_dst) + (0.5 if dst_stagger == "c" else 0.0)) / n_dst
    if span is not None:
        x = span[0] + x * (span[1] - span[0])
    return x * n_src - (0.5 if src_stagger == "c" else 0.0)


def _weights(fi, n):
    fi = np.clip(fi, 0.0, n - 1.0)
    i0 = np.minimum(np.floor(fi).astype(int), n - 2) if n > 1 else np.zeros(fi.shape, int)
    w = fi - i0
    return i0, np.minimum(i0 + 1, n - 1), w


def regrid(a, dst_shape, stagger_y="c", stagger_x="c", src_y=None, src_x=None, spans=None):
    """Bilinear regrid of a (..., ny, nx) field to (..., NY, NX), same domain.

    stagger_y/x: where the destination points are; src_y/x: where the source
    points are (default: the same staggering). spans = ((ay, by), (ax, bx)):
    the destination box as fractions of the source box (None: the same box).
    """
    ny, nx = a.shape[-2:]
    NY, NX = dst_shape
    sy, sx = (None, None) if spans is None else spans
    j0, j1, wy = _weights(fractional_index(ny, NY, stagger_y, src_y, sy), ny)
    i0, i1, wx = _weights(fractional_index(nx, NX, stagger_x, src_x, sx), nx)
    a0 = a[..., j0, :]
    a1 = a[..., j1, :]
    row = (1.0 - wy)[:, None] * a0 + wy[:, None] * a1
    return (1.0 - wx) * row[..., i0] + wx * row[..., i1]


def surface_pressure_to(ps_src, h_src, h_dst, t_low):
    """p_s moved hydrostatically from h_src to h_dst (layer-mean temperature)."""
    dh = h_dst - h_src
    t_mean = t_low - 0.5 * LAPSE * dh
    return ps_src * np.exp(-G * dh / (RD * t_mean))


def remap_column(q, p_src, p_dst, extrapolate_top=False):
    """q given at p_src (k, ...) increasing downward; values at p_dst (m, ...).

    Linear in ln p; constant beyond the end levels, except above the top
    level when extrapolate_top is set (theta: the 3 km top level can sit above
    the 12 km one, and a constant theta there would be a neutral layer under
    the lid).
    """
    lps, lpd = np.log(p_src), np.log(p_dst)
    n = q.shape[0]
    out = np.empty(p_dst.shape)
    for m in range(p_dst.shape[0]):
        x = lpd[m]
        j = (lps < x[None]).sum(axis=0)          # levels above the target
        jj = np.clip(j, 1, n - 1)
        xa = np.take_along_axis(lps, (jj - 1)[None], 0)[0]
        xb = np.take_along_axis(lps, jj[None], 0)[0]
        qa = np.take_along_axis(q, (jj - 1)[None], 0)[0]
        qb = np.take_along_axis(q, jj[None], 0)[0]
        w = (x - xa) / (xb - xa)
        w = np.clip(w, -np.inf if extrapolate_top else 0.0, 1.0)
        w = np.where((j > 0) | (w >= 0.0), np.clip(w, 0.0, 1.0), w)
        out[m] = (1.0 - w) * qa + w * qb
    return out


def frame_from_coarse(u, v, theta, pi, sigma_src, p_top, h_src, h_dst, sigma_dst, spans=None):
    """
    One 12 km state (u, v, theta on sigma_src levels; pi = p_s - p_top) as a
    3 km state on sigma_dst over h_dst. Returns (u, v, theta, pi) on the
    destination grid with its C-grid staggering.
    """
    NY, NX = h_dst.shape
    sig_s = np.asarray(sigma_src, float)[:, None, None]
    sig_d = np.asarray(sigma_dst, float)[:, None, None]
    ps_src = pi + p_top
    t_low = theta[-1] * ((p_top + sig_s[-1, 0, 0] * pi) / P0) ** KAPPA

    def column(stagger_y, stagger_x, h_d):
        # ps, h and the low-level temperature are cell-centred on the source
        ps = regrid(ps_src, (NY, NX), stagger_y, stagger_x, "c", "c", spans=spans)
        hs = regrid(h_src, (NY, NX), stagger_y, stagger_x, "c", "c", spans=spans)
        tl = regrid(t_low, (NY, NX), stagger_y, stagger_x, "c", "c", spans=spans)
        ps_d = surface_pressure_to(ps, hs, h_d, tl)
        p_s = p_top + sig_s * (ps - p_top)[None]
        p_d = p_top + sig_d * (ps_d - p_top)[None]
        return ps_d, p_s, p_d

    # centres: theta and pi
    ps_c, p_s, p_d = column("c", "c", h_dst)
    # ln theta linear in ln p: exact for an isothermal or constant-lapse-in-ln-p
    # column, where theta itself is an exponential of ln p.
    th = np.exp(remap_column(np.log(regrid(theta, (NY, NX), spans=spans)), p_s, p_d, extrapolate_top=True))
    # u points: west faces; terrain there is the mean of the two cells
    h_u = 0.5 * (h_dst + np.concatenate([h_dst[:, :1], h_dst[:, :-1]], axis=1))
    _, p_su, p_du = column("c", "f", h_u)
    uu = remap_column(regrid(u, (NY, NX), "c", "f", spans=spans), p_su, p_du)
    # v points: south faces
    h_v = 0.5 * (h_dst + np.concatenate([h_dst[:1], h_dst[:-1]], axis=0))
    _, p_sv, p_dv = column("f", "c", h_v)
    vv = remap_column(regrid(v, (NY, NX), "f", "c", spans=spans), p_sv, p_dv)
    return uu, vv, th, ps_c - p_top


def scalar_from_coarse(q, pi, sigma_src, p_top, h_src, h_dst, sigma_dst, theta=None, spans=None):
    """A cell-centred 12 km scalar (a mixing ratio) on the 3 km grid and
    sigma_dst: regridded and remapped linear in ln p, held beyond the end
    levels, so a non-negative field stays non-negative (CAM stage S7)."""
    NY, NX = h_dst.shape
    sig_s = np.asarray(sigma_src, float)[:, None, None]
    sig_d = np.asarray(sigma_dst, float)[:, None, None]
    ps_src = pi + p_top
    t_low = theta[-1] * ((p_top + sig_s[-1, 0, 0] * pi) / P0) ** KAPPA
    ps = regrid(ps_src, (NY, NX), spans=spans)
    hs = regrid(h_src, (NY, NX), spans=spans)
    tl = regrid(t_low, (NY, NX), spans=spans)
    ps_d = surface_pressure_to(ps, hs, h_dst, tl)
    p_s = p_top + sig_s * (ps - p_top)[None]
    p_d = p_top + sig_d * (ps_d - p_top)[None]
    return np.maximum(remap_column(regrid(q, (NY, NX), spans=spans), p_s, p_d), 0.0)


def spans_for(src, dst, tol=1e-4):
    """((ay, by), (ax, bx)): box dst as fractions of box src, or None when the
    two boxes are the same (to tol degrees) -- then nothing changes (test N1)."""
    keys = ("lat_min", "lat_max", "lon_min", "lon_max")
    if src is None or dst is None or all(abs(float(src[k]) - float(dst[k])) <= tol for k in keys):
        return None
    ly, lx = src["lat_max"] - src["lat_min"], src["lon_max"] - src["lon_min"]
    sp = (((dst["lat_min"] - src["lat_min"]) / ly, (dst["lat_max"] - src["lat_min"]) / ly),
          ((dst["lon_min"] - src["lon_min"]) / lx, (dst["lon_max"] - src["lon_min"]) / lx))
    if min(sp[0][0], sp[1][0]) < -1e-6 or max(sp[0][1], sp[1][1]) > 1 + 1e-6:
        raise ValueError(f"the 3 km box {dst} is not inside the 12 km box {src}")
    return sp


def domain_of_file(z):
    """The box a forecast file covers, from its cell-centre lat/lon, or None."""
    if not {"lat", "lon"} <= set(z.files):
        return None
    try:
        la, lo = np.asarray(z["lat"], float), np.asarray(z["lon"], float)
    except ValueError:          # stored as None (object array)
        return None
    if la.ndim == 0 or lo.ndim == 0:
        return None
    la = la[:, 0] if la.ndim == 2 else la
    lo = lo[0, :] if lo.ndim == 2 else lo
    dla, dlo = (la[-1] - la[0]) / (la.size - 1), (lo[-1] - lo[0]) / (lo.size - 1)
    return {"lat_min": la[0] - 0.5 * dla, "lat_max": la[-1] + 0.5 * dla,
            "lon_min": lo[0] - 0.5 * dlo, "lon_max": lo[-1] + 0.5 * dlo}


def frames_from_forecast(path, h_dst, sigma_dst, hours, moist=False, dst_domain=None):
    """[(t_seconds, u, v, theta, pi)] from a 12 km forecast file, up to `hours`.

    moist (CAM stage S7): each tuple gets a sixth item, a dict of qv, qc, qr on
    the destination grid when the file has them (empty otherwise).
    dst_domain (test N1): the 3 km box; when the file covers a larger box the
    frames are taken from the part of it the 3 km grid covers.
    """
    z = np.load(path, allow_pickle=False)
    has_q = moist and all(k in z.files for k in ("qv", "qc", "qr"))
    h_src = np.asarray(z["terrain"], float)
    spans = spans_for(domain_of_file(z), dst_domain) if dst_domain is not None else None
    out = []
    for n, t in enumerate(np.asarray(z["times_s"], float)):
        if t > hours * 3600.0 + 1e-6:
            break
        fr = (float(t),) + frame_from_coarse(
            z["u"][n], z["v"][n], z["theta"][n], z["pi"][n], z["sigma"],
            float(z["p_top"]), h_src, h_dst, sigma_dst, spans=spans)
        if moist:
            extra = {}
            if has_q:
                th_n, pi_n = np.asarray(z["theta"][n], float), np.asarray(z["pi"][n], float)
                for k in ("qv", "qc", "qr"):
                    extra[k] = scalar_from_coarse(np.asarray(z[k][n], float), pi_n, z["sigma"],
                                                  float(z["p_top"]), h_src, h_dst, sigma_dst,
                                                  theta=th_n, spans=spans)
            fr = fr + (extra,)
        out.append(fr)
    return out

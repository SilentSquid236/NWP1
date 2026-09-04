"""
Pressure levels to sigma levels.

WHY THIS IS THE BLOCKING PIECE

The sigma core is the model. Everything measured since the coordinate change
-- terrain, boundary layer, convective adjustment -- lives there. But the
end-to-end driver still built a `Primitive3D` on pressure levels, because HRRR
delivers isobaric data and nothing converted it. So the only core reachable
from a real forecast was the one that diverges in two to three hours
(P-14 in docs/PROBLEMS.md), while the one that reaches 12/12 could only be run
on idealised states. This module is what connects them.

WHAT IT HAS TO GET RIGHT

    sigma = (p - p_top) / (p_s - p_top)

so the target pressure of every model level depends on the SURFACE PRESSURE of
that column, which depends on the terrain. Three steps, in order:

  1. terrain height h  ->  surface pressure p_s, by finding the pressure at
     which the analysis geopotential height equals h. HRRR already gives
     geopotential height on every pressure level, so this is an interpolation
     rather than a hydrostatic guess.
  2. p_s  ->  the target pressure of each sigma level.
  3. analysis columns  ->  those pressures, interpolated in log(p).

INTERPOLATION IN LOG(p), NOT p

Between two pressure levels a field varies far more nearly linearly in log(p)
than in p, because log(p) is close to linear in height. At 12 km with 20
levels from 1000 to 200 hPa the gaps are large (1000 -> 975 -> 950 near the
ground, but 300 -> 250 -> 200 aloft), so the difference is not cosmetic.

EXTRAPOLATION, WHICH IS WHERE THIS GETS DANGEROUS

Sigma levels near the ground over LOW terrain sit at pressures higher than the
analysis's lowest level -- 1000 hPa is not the surface at sea level, it is
about 100 m up. Something has to be said about the layer below, and the choice
matters:

  * theta is extrapolated along the LAPSE RATE of the lowest two analysis
    levels, not held constant. Holding theta constant makes the near-surface
    layer exactly neutral, which the convective adjustment then reads as being
    on the edge of instability everywhere on the first step.
  * wind is held constant below the lowest level. Extrapolating a shear
    downward produces surface winds that the drag scheme then fights.

Above the top analysis level the fields are held constant. That region is
inside the sponge, so its detail is absorbed rather than forecast.
"""

import numpy as np

from sigma import RD, G0, P0, KAPPA


def interp_log_p(field, p_src, p_dst, extrapolate="hold"):
    """
    Interpolate columns from source pressures to target pressures in log(p).

    `field` is (nsrc, ny, nx), `p_src` is 1D (nsrc,) in Pa, `p_dst` is
    (nz, ny, nx) in Pa. Source levels may be given in either order.

    `extrapolate` controls what happens outside the source range:
      "hold"    -- nearest source value (used for wind)
      "linear"  -- continue the gradient of the two nearest source levels,
                   in log(p) (used for theta)
    """
    field = np.asarray(field, dtype=float)
    p_src = np.asarray(p_src, dtype=float).ravel()
    p_dst = np.asarray(p_dst, dtype=float)

    if field.shape[0] != p_src.size:
        raise ValueError(f"field has {field.shape[0]} levels but p_src has "
                         f"{p_src.size}")

    x_src = np.log(p_src)
    order = np.argsort(x_src)              # searchsorted needs increasing x
    x_src = x_src[order]
    f_src = field[order]

    x = np.log(p_dst)
    n = x_src.size

    # Bracketing indices, clipped so the arithmetic is always defined; the
    # out-of-range cases are then overwritten explicitly below.
    hi = np.clip(np.searchsorted(x_src, x), 1, n - 1)
    lo = hi - 1

    x0 = x_src[lo]
    x1 = x_src[hi]
    f0 = np.take_along_axis(f_src, lo, axis=0) if f_src.ndim == 3 else f_src[lo]
    f1 = np.take_along_axis(f_src, hi, axis=0) if f_src.ndim == 3 else f_src[hi]

    w = (x - x0) / (x1 - x0)
    out = f0 + w * (f1 - f0)

    if extrapolate == "hold":
        below = x < x_src[0]               # lower pressure than any source
        above = x > x_src[-1]              # higher pressure than any source
        out = np.where(below, f_src[0], out)
        out = np.where(above, f_src[-1], out)
    elif extrapolate != "linear":
        raise ValueError("extrapolate must be 'hold' or 'linear'")

    return out


def surface_pressure_from_heights(z, p_levels_pa, terrain):
    """
    Surface pressure implied by the analysis, at the model's terrain height.

    `z` is geopotential height (nsrc, ny, nx) in metres, `p_levels_pa` the
    matching pressures, `terrain` the model orography (ny, nx) in metres.

    Found by interpolating log(p) against height rather than assuming a lapse
    rate, so it inherits the analysis's own stratification. Where the terrain
    lies below the lowest analysis level -- which is most of a coastal domain,
    since 1000 hPa sits about 100 m above sea level -- the hydrostatic
    relation extends it downward using the temperature implied by the lowest
    two levels.
    """
    z = np.asarray(z, dtype=float)
    p = np.asarray(p_levels_pa, dtype=float).ravel()
    terrain = np.asarray(terrain, dtype=float)

    order = np.argsort(-p)                 # highest pressure (lowest) first
    p = p[order]
    z = z[order]

    lnp = np.log(p).reshape(-1, 1, 1) * np.ones_like(z)

    # Walk down the levels, blending between the pair that brackets the
    # terrain height. Heights DECREASE with index here (index 0 is the lowest
    # level), so the bracket is z[k] >= h > z[k+1].
    # After sorting by descending pressure, index 0 is the HIGHEST pressure
    # and therefore the LOWEST height: heights increase with index. Getting
    # this backwards leaves every column above the lowest analysis level
    # pinned at that level's pressure, which showed up as a 253 hPa error over
    # 2500 m terrain -- caught by comparing against a standard atmosphere,
    # where the right answer is known in closed form.
    ln_ps = np.full(terrain.shape, lnp[0])
    found = np.zeros(terrain.shape, dtype=bool)
    for k in range(z.shape[0] - 1):
        bracket = (~found) & (terrain >= z[k]) & (terrain < z[k + 1])
        if bracket.any():
            w = (terrain - z[k]) / np.maximum(z[k + 1] - z[k], 1e-6)
            ln_ps = np.where(bracket, lnp[k] + w * (lnp[k + 1] - lnp[k]), ln_ps)
            found |= bracket

    # Terrain above the highest analysis level would mean a mountain above
    # 200 hPa; refuse rather than silently clamp.
    too_high = (~found) & (terrain >= z[-1])
    if too_high.any():
        raise ValueError(f"terrain reaches {terrain.max():.0f} m, above the "
                         f"top analysis level at {z[-1].max():.0f} m")

    # Terrain below the lowest analysis level: hydrostatic extension.
    # dln(p)/dz = -g/(R T), with T from the lowest layer.
    below = terrain < z[0]
    if below.any():
        # A representative temperature for the layer being extended.
        dz = np.maximum(z[1] - z[0], 1.0)
        T_layer = G0 * dz / (RD * np.log(p[0] / p[1]))
        ln_ps = np.where(below,
                         lnp[0] + G0 * (z[0] - terrain) / (RD * T_layer),
                         ln_ps)

    return np.exp(ln_ps)


def theta_from_T(T, p):
    return T * (P0 / p) ** KAPPA


def pressure_to_sigma(u, v, T, z, p_levels_pa, terrain, lev,
                      p_surface=None):
    """
    Convert an isobaric analysis to a sigma-coordinate model state.

    Returns (pi, u_sigma, v_sigma, theta_sigma) ready to assign to a
    `PrimitiveSigma`. `p_surface` may be supplied directly when the analysis
    carries it; otherwise it is derived from the geopotential heights.
    """
    p_levels_pa = np.asarray(p_levels_pa, dtype=float).ravel()

    p_s = (surface_pressure_from_heights(z, p_levels_pa, terrain)
           if p_surface is None else np.asarray(p_surface, dtype=float))
    pi = p_s - lev.p_top
    if np.any(pi <= 0):
        raise ValueError(
            f"surface pressure at or below the model lid: min p_s "
            f"{p_s.min():.0f} Pa, p_top {lev.p_top:.0f} Pa. Terrain too high "
            f"for this lid, or the heights are in the wrong units.")

    p_dst = lev.pressure(pi)

    u_s = interp_log_p(u, p_levels_pa, p_dst, extrapolate="hold")
    v_s = interp_log_p(v, p_levels_pa, p_dst, extrapolate="hold")

    # Theta is interpolated, not temperature: theta is the model's variable
    # and is conserved by dry adiabatic motion, so interpolating it does not
    # invent heating the way interpolating T through a deep layer does.
    theta_src = theta_from_T(np.asarray(T, dtype=float),
                             p_levels_pa.reshape(-1, 1, 1))
    th_s = interp_log_p(theta_src, p_levels_pa, p_dst, extrapolate="linear")

    return pi, u_s, v_s, th_s


def hydrostatic_theta(phi_target, pi, lev, phi_surface):
    """
    Choose theta on the model's levels so the model's OWN discrete hydrostatic
    integral reproduces `phi_target`.

    WHY THIS EXISTS (P-46 in docs/PROBLEMS.md)

    Interpolating theta reproduces the analysis's stratification but not the
    model's discrete integral over sigma layers. Measured over a 1500 m
    mountain: the model's geopotential differs from the analysis's by only
    0.06-2.7 m in the horizontal MEAN -- the integration is fine -- but the
    error SPREADS horizontally, from 3 m at the ground to 140 m at the lid,
    because a column over the mountain samples the analysis at different
    pressures than a column beside it. A horizontally varying geopotential
    error is a pressure-gradient force with nothing balancing it: 2.7 m/s per
    hour of acceleration on an atmosphere at rest, against 1e-10 m/s per hour
    for a state the model built itself.

    HOW

    `hydrostatic_geopotential` is triangular -- phi[k] depends only on T[k]
    and the levels below it -- so it inverts exactly by the reverse recursion:

        T[-1] = (phi[-1] - phi_s) / (R ln(p_s / p[-1]))
        T[k]  = 2 (phi[k] - phi[k+1]) / (R ln(p[k+1]/p[k])) - T[k+1]

    THE CATCH, AND WHY THE CALLER GETS TOLD

    That `- T[k+1]` makes the inverse an alternating recursion: an error at one
    level flips sign and persists upward as a (-1)^k sawtooth. The inversion is
    exact, but exactness is not the same as usefulness -- a temperature profile
    with a 2dx vertical wiggle is worse than a small geopotential error. The
    returned info reports the sawtooth amplitude and whether the result is
    still statically stable, so the caller can refuse it rather than integrate
    a profile that only looks right through the operator that produced it.
    """
    phi_target = np.asarray(phi_target, dtype=float)
    p = lev.pressure(pi)
    p_s = lev.p_top + pi

    T = np.empty_like(phi_target)
    T[-1] = (phi_target[-1] - phi_surface) / (RD * np.log(p_s / p[-1]))
    for k in range(lev.nz - 2, -1, -1):
        T[k] = (2.0 * (phi_target[k] - phi_target[k + 1])
                / (RD * np.log(p[k + 1] / p[k])) - T[k + 1])

    theta = T * (P0 / p) ** KAPPA

    # Sawtooth diagnostic: the amplitude of the (-1)^k component of T, which
    # is the mode the inverse recursion amplifies.
    sign = (-1.0) ** np.arange(lev.nz).reshape(-1, 1, 1)
    saw = float(np.abs((T * sign).mean(axis=0)).max())
    stable = bool(np.all(theta[:-1] > theta[1:]))

    info = {"sawtooth_K": saw,
            "statically_stable": stable,
            "T_min": float(np.nanmin(T)),
            "T_max": float(np.nanmax(T))}
    return theta, info

"""
Vertical machinery for a hydrostatic pressure-coordinate model.

Pressure as the vertical coordinate is the simplest choice that matches the
data we have: HRRR gives us fields on isobaric surfaces already, so no
interpolation is needed on the way in or out.

Its limitation is real and worth stating up front: **isobaric surfaces
intersect terrain**. A pressure-coordinate model therefore assumes a flat
lower boundary. For the Northeast domain, with the Appalachians running
through it, that is a genuine approximation -- terrain drives a great deal of
the real weather there. The fix is a terrain-following (sigma or hybrid)
coordinate, which is the natural next step once the dry core is proven.

Conventions used throughout:

  * Level index 0 is the BOTTOM (highest pressure), index nz-1 is the TOP.
  * p decreases with increasing index.
  * omega = dp/dt is the vertical velocity in pressure coordinates. It is
    NEGATIVE for rising motion, because rising air moves toward lower p.
"""

import numpy as np

RD = 287.058        # J/(kg K)  gas constant for dry air
CP = 1004.6         # J/(kg K)  specific heat at constant pressure
KAPPA = RD / CP     # ~0.2857
P0 = 100_000.0      # Pa, reference pressure for potential temperature
G0 = 9.80665


def exner(p):
    """Exner function pi = (p/p0)^kappa. T = theta * pi."""
    return (p / P0) ** KAPPA


def theta_from_T(T, p):
    return T / exner(p)


def T_from_theta(theta, p):
    return theta * exner(p)


class PressureLevels:
    """
    Immutable description of the vertical grid.

    levels_hPa : list, ordered bottom -> top (decreasing pressure)
    """

    def __init__(self, levels_hPa):
        p = np.asarray(levels_hPa, dtype=float) * 100.0     # -> Pa
        if np.any(np.diff(p) >= 0):
            raise ValueError("levels must be ordered bottom -> top "
                             "(strictly decreasing pressure)")
        self.p = p
        self.nz = len(p)
        self.pi = exner(p)

        # Layer thicknesses between levels: dp[k] = p[k] - p[k+1] > 0
        self.dp = p[:-1] - p[1:]

    @property
    def p_bottom(self):
        return self.p[0]

    @property
    def p_top(self):
        return self.p[-1]

    def column(self, shape):
        """Broadcast p to (nz, ny, nx) for elementwise use."""
        return self.p.reshape(-1, *([1] * len(shape)))

    def __repr__(self):
        return (f"PressureLevels({self.nz} levels, "
                f"{self.p[0]/100:.0f} -> {self.p[-1]/100:.0f} hPa)")


def hydrostatic_geopotential(theta, lev, phi_surface=0.0):
    """
    Integrate the hydrostatic relation upward to get geopotential.

        dPhi/dp = -R T / p        =>       dPhi = -R T dln(p)

    Between adjacent levels this is the hypsometric equation. T is taken as
    the layer mean, which makes the discrete integral second-order accurate.

    Returns Phi with the same shape as theta.
    """
    T = T_from_theta(theta, lev.p.reshape(-1, 1, 1))
    phi = np.empty_like(theta)
    phi[0] = phi_surface

    for k in range(lev.nz - 1):
        T_layer = 0.5 * (T[k] + T[k + 1])
        dlnp = np.log(lev.p[k] / lev.p[k + 1])          # > 0 going up
        phi[k + 1] = phi[k] + RD * T_layer * dlnp

    return phi


def diagnose_omega(div, lev):
    """
    Vertical velocity from mass continuity.

        du/dx + dv/dy + d(omega)/dp = 0     =>     omega(p) = -int_top^p D dp'

    With a rigid lid and a flat lower boundary, omega must vanish at BOTH
    ends. Integrating from the top generally leaves a non-zero residual at the
    bottom, because the column-integrated divergence is not exactly zero on a
    discrete grid. We remove that residual with a correction linear in p --
    the standard treatment. Physically it is the adjustment that keeps total
    column mass fixed.

    div : (nz, ny, nx) horizontal divergence
    """
    nz = lev.nz
    omega = np.zeros_like(div)

    # Integrate downward from the top (index nz-1) toward the bottom.
    for k in range(nz - 2, -1, -1):
        d_layer = 0.5 * (div[k] + div[k + 1])
        dp = lev.p[k] - lev.p[k + 1]                    # > 0
        omega[k] = omega[k + 1] + d_layer * dp

    # omega[0] should be zero; distribute the residual linearly in p.
    residual = omega[0]
    frac = ((lev.p - lev.p_top) /
            (lev.p_bottom - lev.p_top)).reshape(-1, 1, 1)
    omega -= residual * frac

    return -omega


def static_stability(theta, lev):
    """
    dtheta/dp. Negative means theta increases upward = statically stable.

    A dry adiabatic atmosphere has dtheta/dp = 0. Positive values mean theta
    DEcreases with height, which is convectively unstable -- a dry hydrostatic
    model has no way to remove that instability, so it will simply amplify.
    """
    dth = np.gradient(theta, lev.p, axis=0)
    return dth

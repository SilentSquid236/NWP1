"""
Surface layer: drag on the lowest model level.

WHY THIS MATTERS HERE

The Richardson-number mixing scheme redistributes shear within a column but
cannot remove momentum from it. Real momentum leaves the atmosphere through
the ground. Without that sink, the near-surface wind has nothing holding it
back, shear against the layers above keeps building, and the mixing scheme is
left fighting a source it cannot switch off.

Adding drag also produces a physical signature that is impossible to fake: the
EKMAN SPIRAL. Friction slows the near-surface wind, Coriolis no longer
balances the pressure gradient, and the flow turns across the isobars toward
low pressure -- typically 10-30 degrees over land. Any implementation that
gets the sign or the magnitude wrong fails that test visibly.

FORMULATION

Bulk aerodynamic, with a log-law drag coefficient:

    tau = rho * Cd * |V| * V                (surface stress, N/m^2)
    Cd  = (kappa / ln(z1/z0))^2 * F(Ri_b)   (neutral value, times stability)

kappa = 0.4 is von Karman's constant, z1 the height of the lowest model level,
z0 the roughness length. F(Ri_b) is the Louis-type stability correction:
drag increases when the surface layer is unstable (convective gusts reach
down) and falls toward zero when it is strongly stable (a decoupled nocturnal
inversion). Neglecting F entirely would apply daytime drag to a calm, clear
night -- a known way to destroy nocturnal boundary-layer structure.

ROUGHNESS LENGTHS (m), conventional values:

    open sea            0.0002
    grassland           0.03
    cropland            0.1
    forest              1.0
    urban               1.5

The Northeast domain is mostly forest and cropland, so 0.1-1.0; the coastal
and offshore part of the domain is four orders of magnitude smoother, which
is why a single domain-wide value is a poor approximation there.
"""

import numpy as np

from sigma import RD, G0, P0, KAPPA

VON_KARMAN = 0.4

ROUGHNESS = {
    "sea": 0.0002,
    "grass": 0.03,
    "cropland": 0.1,
    "forest": 1.0,
    "urban": 1.5,
}


def lowest_level_height(theta, pi, lev):
    """
    Height of the lowest full level above ground, hydrostatically.

    This is z1 in the log law, and the drag coefficient depends on it
    logarithmically -- so it must come from the model's own layer thickness,
    not a hardcoded constant, or the drag silently changes with the vertical
    grid.
    """
    p = lev.pressure(pi)
    p_half = lev.pressure_half(pi)
    T = theta * (p / P0) ** KAPPA
    # From the ground (p_half[-1]) up to the lowest full level (p[-1]).
    return RD * T[-1] / G0 * np.log(p_half[-1] / p[-1])


def neutral_drag_coefficient(z1, z0):
    """Cd for a neutrally stratified surface layer, from the log law."""
    z1 = np.maximum(z1, 2.0 * z0)          # level must sit above the roughness
    return (VON_KARMAN / np.log(z1 / z0)) ** 2


def bulk_richardson(u1, v1, theta1, theta_s, z1):
    """
    Bulk Richardson number across the surface layer.

    Positive => stable (surface colder than the air above, e.g. a clear
    night). Negative => unstable (surface warmer, convective).
    """
    speed2 = np.maximum(u1 ** 2 + v1 ** 2, 0.01)
    return (G0 * z1 * (theta1 - theta_s)) / (theta1 * speed2)


def stability_function(Ri_b, ri_crit=0.2):
    """
    Louis-type correction to the neutral drag coefficient.

    Unstable (Ri_b < 0): enhanced, capped -- convective plumes couple the
    surface to the flow above.
    Stable (0 <= Ri_b < ri_crit): suppressed quadratically.
    Strongly stable (Ri_b >= ri_crit): decoupled, essentially no drag. This is
    what lets a nocturnal inversion form instead of being mixed away.
    """
    F = np.ones_like(Ri_b)

    unstable = Ri_b < 0
    F[unstable] = np.minimum(1.0 + 10.0 * np.abs(Ri_b[unstable]), 4.0)

    stable = (Ri_b >= 0) & (Ri_b < ri_crit)
    F[stable] = (1.0 - Ri_b[stable] / ri_crit) ** 2

    F[Ri_b >= ri_crit] = 0.0
    return F


def surface_drag(u, v, theta, pi, lev, z0=0.1, theta_s=None, ri_crit=0.2):
    """
    Momentum tendency from surface stress, applied to the lowest level.

        du/dt = -tau_x / (rho * dz1)  =  -Cd |V| u / dz1

    Returns (du, dv, info). Only the lowest level is affected; the mixing
    scheme spreads the effect upward from there, which is how a boundary layer
    actually deepens.
    """
    p = lev.pressure(pi)
    p_half = lev.pressure_half(pi)
    T = theta * (p / P0) ** KAPPA

    z1 = lowest_level_height(theta, pi, lev)
    dz1 = RD * T[-1] / G0 * np.log(p_half[-1] / p_half[-2])

    u1, v1 = u[-1], v[-1]
    speed = np.sqrt(u1 ** 2 + v1 ** 2)

    cd_n = neutral_drag_coefficient(z1, z0)
    if theta_s is None:
        F = np.ones_like(u1)
        Ri_b = np.zeros_like(u1)
    else:
        Ri_b = bulk_richardson(u1, v1, theta[-1], theta_s, z1)
        F = stability_function(Ri_b, ri_crit)
    cd = cd_n * F

    du = np.zeros_like(u)
    dv = np.zeros_like(v)
    du[-1] = -cd * speed * u1 / dz1
    dv[-1] = -cd * speed * v1 / dz1

    info = {"z1": z1, "dz1": dz1, "cd": cd, "cd_neutral": cd_n,
            "Ri_bulk": Ri_b, "u_star": np.sqrt(cd * speed ** 2)}
    return du, dv, info


def drag_stability_dt(info, safety=0.4):
    """
    Explicit drag is a relaxation with timescale dz1 / (Cd |V|); the timestep
    must stay well inside it or the lowest level oscillates.
    """
    cd = np.asarray(info["cd"])
    us = np.asarray(info["u_star"])
    speed = np.where(us > 0, us ** 2 / np.maximum(cd, 1e-12), 0.0)
    rate = cd * np.sqrt(np.maximum(speed, 0.0)) / np.asarray(info["dz1"])
    rmax = float(np.max(rate))
    return np.inf if rmax <= 0 else safety / rmax

"""
Richardson-number dependent vertical mixing.

WHY THIS IS PHYSICS, NOT A NUMERICAL PATCH

The model develops shear instability and has nothing to dissipate it.
Measured on flat ground with a realistic jet, survival tracks the vertical
shear between adjacent levels:

    3.2 m/s per level  ->  12/12 forecast hours
    4.2 m/s per level  ->   3/12
    8.4 m/s per level  ->   2/12

A shear of ~8 m/s across a ~400 m layer with N ~ 0.015 gives a Richardson
number near 0.5, approaching the Ri = 0.25 threshold below which shear
instability is expected in the real atmosphere. The instability is REAL. What
is missing is the turbulence that would mix it away.

Every operational model carries a scheme like this. Without one, any flow that
develops locally low Ri -- which real analyses do, in thin layers -- has no
sink and grows without bound.

    Ri = N^2 / S^2,     N^2 = (g/theta) dtheta/dz,   S^2 = (du/dz)^2 + (dv/dz)^2

    Ri >= Ri_c        no mixing (stable, laminar)
    0 < Ri < Ri_c     mixing increases as Ri falls
    Ri <= 0           statically unstable: maximum mixing

The Louis (1979) family of stability functions is the classic formulation;
this is the same shape with the constants exposed.
"""

import numpy as np

from sigma import RD, G0, P0, KAPPA

RI_CRIT = 0.25          # below this, shear overcomes stratification
K_MAX = 100.0           # m^2/s ceiling on the eddy diffusivity
MIXING_LENGTH = 150.0   # m


def richardson(u, v, theta, pi, lev):
    """
    Gradient Richardson number at layer interfaces, plus N^2 and shear^2.

    Everything is computed on HALF levels, between the full levels where u, v
    and theta live -- that is where the shear and the stratification are both
    naturally defined.
    """
    p = lev.pressure(pi)
    T = theta * (p / P0) ** KAPPA

    # Layer thickness in metres, hydrostatic.
    T_half = 0.5 * (T[:-1] + T[1:])
    dz = RD * T_half / G0 * np.log(p[1:] / p[:-1])       # >0, index 0 = top
    dz = np.maximum(np.abs(dz), 1.0)

    dth = theta[1:] - theta[:-1]
    th_half = 0.5 * (theta[:-1] + theta[1:])

    # z increases as index decreases (0 is the lid), so dtheta/dz flips sign.
    N2 = -(G0 / th_half) * dth / dz

    du = (u[1:] - u[:-1]) / dz
    dv = (v[1:] - v[:-1]) / dz
    S2 = du ** 2 + dv ** 2

    with np.errstate(divide="ignore", invalid="ignore"):
        Ri = np.where(S2 > 1e-12, N2 / S2, np.inf)
    return Ri, N2, S2, dz


def eddy_diffusivity(Ri, S2, ri_crit=RI_CRIT, k_max=K_MAX,
                     mixing_length=MIXING_LENGTH):
    """
    K = l^2 * |S| * f(Ri), with f falling to zero at Ri_c.

    f(Ri) = (1 - Ri/Ri_c)^2 for 0 <= Ri < Ri_c, 1 for Ri <= 0 (static
    instability mixes at full strength), 0 above Ri_c.
    """
    S = np.sqrt(np.maximum(S2, 0.0))

    f = np.zeros_like(Ri)
    unstable = Ri <= 0
    marginal = (Ri > 0) & (Ri < ri_crit)
    f[unstable] = 1.0
    f[marginal] = (1.0 - Ri[marginal] / ri_crit) ** 2

    K = mixing_length ** 2 * S * f
    return np.clip(K, 0.0, k_max)


def vertical_mixing(u, v, theta, pi, lev, ri_crit=RI_CRIT, k_max=K_MAX,
                    mixing_length=MIXING_LENGTH):
    """
    Tendencies from turbulent vertical mixing: d/dz ( K du/dz ), etc.

    Returns (du, dv, dtheta), each shaped like the input. Fluxes vanish at the
    lid and the ground, so this redistributes momentum and heat within a
    column without creating or destroying either.
    """
    Ri, N2, S2, dz = richardson(u, v, theta, pi, lev)
    K = eddy_diffusivity(Ri, S2, ri_crit, k_max, mixing_length)

    def mix(a):
        flux = K * (a[1:] - a[:-1]) / dz            # at interfaces
        out = np.zeros_like(a)
        # Divergence of the flux; zero flux through top and bottom boundaries.
        out[1:-1] = (flux[1:] - flux[:-1]) / (0.5 * (dz[1:] + dz[:-1]))
        out[0] = flux[0] / dz[0]
        out[-1] = -flux[-1] / dz[-1]
        return out

    return mix(u), mix(v), mix(theta), K


def mixing_stability_dt(K, dz, safety=0.4):
    """Explicit diffusion limit: dt <= safety * dz^2 / K."""
    Kmax = float(np.max(K))
    if Kmax <= 0:
        return np.inf
    return float(safety * np.min(dz) ** 2 / Kmax)

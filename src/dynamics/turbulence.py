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

from backend import xp_of

from sigma import RD, G0, P0, KAPPA

RI_CRIT = 0.25          # below this, shear overcomes stratification
K_MAX = 200.0           # m^2/s ceiling on the eddy diffusivity
MIXING_LENGTH = 150.0   # m

# WHY K_MAX IS 200 AND NOT 100 (P-40, measured 2026-09-11/12)
#
# 4000 m terrain, 8-level sponge, clean and filtered, 12-hour ceiling:
#
#   K_MAX      100   110   125   150   200   250   300   1000
#   survived  6/12  6/12  7/12  8/12  8/12  8/12  8/12   8/12
#
# A ramp between 100 and 150, flat above it. At 100 the ceiling was
# truncating the diffusivity the scheme itself asked for on 0.11% of
# interfaces -- about one in a thousand, in the breaking region -- and that
# truncation cost two forecast hours. Above ~600 the formula never asks for
# more, so the parameter is inert there: at K_MAX 1000 the realized maximum
# is 605 and the clip fraction is 0.00%.
#
# It is NOT suppression, which is the failure mode that matters here (L2,
# and how the first sponge failed in P-16). Compared at a common hour,
# max|u| is 54.7 m/s at every setting and the jet is marginally STRONGER
# with more mixing, 48.1 -> 48.5; mid-level N^2 is 2.214e-04 at every
# setting to four figures. Nothing is being flattened.
#
# The production case pays nothing: 2500 m with a 5-level sponge is 12/12
# at 100, 150 and 300, with max|u| 44.0 / 43.9 / 43.9 and identical jet and
# stratification. 200 is chosen over 150 for margin, and is still at the
# bottom of the 10^2-10^3 m^2/s observed in breaking mountain waves.
#
# What this does NOT explain: a higher ceiling does not reduce the
# overturning fraction (0.369% -> 0.379%, slightly the wrong way), so the
# mechanism behind the two hours is open. See P-40.

# THESE THREE ARE DEFAULTS, NOT KNOBS. Do not set them at runtime.
#
# The values below are bound into the function signatures at IMPORT time, so
# `import turbulence; turbulence.K_MAX = 400` changes nothing that has already
# been imported -- the callers keep the value from the moment the module was
# read. That is not a hypothetical: a K_MAX ladder run that way returned peak
# |v| of 44.1 at both 100 and 400, identical to one decimal place, and before
# that it produced a RECORDED NEGATIVE RESULT (P-40: 6/12, 6/12, 6/12 at 100 /
# 300 / 1000) that was read as a clean elimination rather than as a broken
# experiment. Re-run properly, that ladder separates -- see the table
# above.
#
# To vary any of them, pass the value down: PrimitiveSigma(..., k_max=...)
# carries it as instance state and hands it to vertical_mixing on every call.
# test_primitive_sigma.test_mixing_knobs_are_connected asserts that the path
# actually works, so a future refactor that re-freezes it fails a suite instead
# of quietly returning identical numbers.


def richardson(u, v, theta, pi, lev):
    """
    Gradient Richardson number at layer interfaces, plus N^2 and shear^2.

    Everything is computed on HALF levels, between the full levels where u, v
    and theta live -- that is where the shear and the stratification are both
    naturally defined.
    """
    xp = xp_of(u, theta, pi)
    p = lev.pressure(pi)
    T = theta * (p / P0) ** KAPPA

    # Layer thickness in metres, hydrostatic.
    T_half = 0.5 * (T[:-1] + T[1:])
    dz = RD * T_half / G0 * xp.log(p[1:] / p[:-1])       # >0, index 0 = top
    dz = xp.maximum(xp.abs(dz), 1.0)

    dth = theta[1:] - theta[:-1]
    th_half = 0.5 * (theta[:-1] + theta[1:])

    # z increases as index decreases (0 is the lid), so dtheta/dz flips sign.
    N2 = -(G0 / th_half) * dth / dz

    du = (u[1:] - u[:-1]) / dz
    dv = (v[1:] - v[:-1]) / dz
    S2 = du ** 2 + dv ** 2

    with xp.errstate(divide="ignore", invalid="ignore"):
        Ri = xp.where(S2 > 1e-12, N2 / S2, xp.inf)
    return Ri, N2, S2, dz


def eddy_diffusivity(Ri, S2, ri_crit=RI_CRIT, k_max=K_MAX,
                     mixing_length=MIXING_LENGTH):
    """
    K = l^2 * |S| * f(Ri), with f falling to zero at Ri_c.

    f(Ri) = (1 - Ri/Ri_c)^2 for 0 <= Ri < Ri_c, 1 for Ri <= 0 (static
    instability mixes at full strength), 0 above Ri_c.
    """
    xp = xp_of(Ri, S2)
    S = xp.sqrt(xp.maximum(S2, 0.0))

    f = xp.zeros_like(Ri)
    unstable = Ri <= 0
    marginal = (Ri > 0) & (Ri < ri_crit)
    f[unstable] = 1.0
    f[marginal] = (1.0 - Ri[marginal] / ri_crit) ** 2

    K = mixing_length ** 2 * S * f
    return xp.clip(K, 0.0, k_max)


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

    xp = xp_of(u)

    def mix(a):
        flux = K * (a[1:] - a[:-1]) / dz            # at interfaces
        out = xp.zeros_like(a)
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

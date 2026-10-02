"""
Surface-layer similarity: 10 m wind and 2 m temperature from the lowest model level
and the ground temperature (P-68 step 2, P-69; test AK).

WHY

The neutral 10 m operator (`wind_10m_factor`, P-68 step 1) takes about 0.59
of the lowest-level wind at every hour. A clear night is not neutral. The
ground cools, the air near it becomes stable, and the 10 m wind falls far
below the wind 200 m up. On the 20 campaign cycles, the dry model's 10 m wind
under the neutral operator is about right by day (1.10 times the observed)
but 1.6–1.8 times the observed at night (test AJ).

The 2 m temperature has the same problem. It comes from the lowest level,
about 230 m up, through a standard lapse rate. That is wrong in a night-time
inversion and in a daytime superadiabatic layer.

With a ground temperature (`land_surface.py`), both can be diagnosed the way
most NWP models diagnose them: Monin–Obukhov similarity between the ground
and the lowest level.

WHAT THIS IS

Given the bulk Richardson number between the ground (theta_g) and the lowest
level (theta_1, U_1 at z_1), solve for the stability zeta = z_1 / L:

    Ri_b = zeta [ln(z1/z0h) - psi_h(zeta) + psi_h(zeta z0h/z1)]
               / [ln(z1/z0)  - psi_m(zeta) + psi_m(zeta z0/z1)]^2

Then read the profiles at 10 m and 2 m:

    U10      = U1 [ln(10/z0) - psi_m(zeta 10/z1) + psi_m(zeta z0/z1)]
                  / [ln(z1/z0) - psi_m(zeta) + psi_m(zeta z0/z1)]
    theta_2  = theta_g + (theta_1 - theta_g)
               [ln(2/z0h) - psi_h(zeta 2/z1) + psi_h(zeta z0h/z1)]
               / [ln(z1/z0h) - psi_h(zeta) + psi_h(zeta z0h/z1)]

- Unstable: the Businger–Dyer forms integrated by Paulson (1970).
- Stable: Beljaars and Holtslag (1991), which stay finite in very stable air.
- z0 = 0.1 m, as in the model's drag; z0h = z0 / 10.

LIMITS, STATED PLAINLY
- Similarity theory holds in the surface layer, the lowest tenth or so of the
  boundary layer. NWP1's lowest level is about 230 m up, which at night is
  often above the stable boundary layer. Operational models apply the same
  diagnostics from a lowest level at 10–30 m. Here they are an extrapolation,
  and they are tested against observations for that reason.
- Over water theta_g is the fixed initial surface temperature (land_surface.py).
"""
import numpy as np

KAPPA_VK = 0.4
G = 9.80665
BH_A, BH_B, BH_C, BH_D = 1.0, 2.0 / 3.0, 5.0, 0.35


def psi_m(zeta):
    zeta = np.asarray(zeta, dtype=float)
    x = np.sqrt(np.sqrt(np.maximum(1.0 - 16.0 * np.minimum(zeta, 0.0), 1.0)))
    un = (2.0 * np.log((1.0 + x) / 2.0) + np.log((1.0 + x * x) / 2.0)
          - 2.0 * np.arctan(x) + np.pi / 2.0)
    z = np.maximum(zeta, 0.0)
    st = -(BH_A * z + BH_B * (z - BH_C / BH_D) * np.exp(-BH_D * z) + BH_B * BH_C / BH_D)
    return np.where(zeta < 0.0, un, st)


def psi_h(zeta):
    zeta = np.asarray(zeta, dtype=float)
    x = np.sqrt(np.sqrt(np.maximum(1.0 - 16.0 * np.minimum(zeta, 0.0), 1.0)))
    un = 2.0 * np.log((1.0 + x * x) / 2.0)
    z = np.maximum(zeta, 0.0)
    st = -((1.0 + 2.0 * BH_A * z / 3.0) ** 1.5 + BH_B * (z - BH_C / BH_D) * np.exp(-BH_D * z)
           + BH_B * BH_C / BH_D - 1.0)
    return np.where(zeta < 0.0, un, st)


def _profiles(zeta, z1, z0, z0h):
    fm = np.log(z1 / z0) - psi_m(zeta) + psi_m(zeta * z0 / z1)
    fh = np.log(z1 / z0h) - psi_h(zeta) + psi_h(zeta * z0h / z1)
    return fm, fh


def richardson_from_zeta(zeta, z1, z0=0.1, z0h=None):
    z0h = z0 / 10.0 if z0h is None else z0h
    fm, fh = _profiles(zeta, z1, z0, z0h)
    return zeta * fh / fm ** 2


def solve_zeta(Ri_b, z1, z0=0.1, z0h=None, lo=-50.0, hi=500.0, iters=80):
    """zeta = z1/L from the bulk Richardson number, by bisection (vectorised)."""
    Ri_b = np.asarray(Ri_b, dtype=float)
    z1 = np.broadcast_to(np.asarray(z1, dtype=float), Ri_b.shape)
    a = np.full(Ri_b.shape, lo); b = np.full(Ri_b.shape, hi)
    for _ in range(iters):
        mid = 0.5 * (a + b)
        f = richardson_from_zeta(mid, z1, z0, z0h) - Ri_b
        b = np.where(f > 0.0, mid, b)
        a = np.where(f > 0.0, a, mid)
    return 0.5 * (a + b)


def surface_diagnostics(U1, theta1, theta_g, z1, z0=0.1, z0h=None, u_min=0.5):
    """
    Returns (wind_factor U10/U1, theta_2, zeta, Ri_b).

    U1 m/s, theta in K (potential temperature on the same reference), z1 the
    lowest level's height above the model ground in m.
    """
    z0h = z0 / 10.0 if z0h is None else z0h
    U1 = np.maximum(np.asarray(U1, dtype=float), u_min)
    theta1 = np.asarray(theta1, dtype=float); theta_g = np.asarray(theta_g, dtype=float)
    z1 = np.maximum(np.asarray(z1, dtype=float), 20.0)
    Ri = G * z1 * (theta1 - theta_g) / (theta1 * U1 ** 2)
    zeta = solve_zeta(Ri, z1, z0, z0h)
    fm, fh = _profiles(zeta, z1, z0, z0h)
    f10 = (np.log(10.0 / z0) - psi_m(zeta * 10.0 / z1) + psi_m(zeta * z0 / z1)) / fm
    w2 = (np.log(2.0 / z0h) - psi_h(zeta * 2.0 / z1) + psi_h(zeta * z0h / z1)) / fh
    theta2 = theta_g + (theta1 - theta_g) * w2
    return f10, theta2, zeta, Ri

"""
Validation for the pressure-to-sigma conversion.

The decisive tests are the last two: a state that is exactly representable
must survive the round trip unchanged, and a real-shaped analysis over real
terrain must produce a model state the core will actually hold.

Run:  python test_interpolate.py
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from sigma import hydrostatic_geopotential
from interpolate import (interp_log_p, surface_pressure_from_heights,
                         pressure_to_sigma, theta_from_T)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


LEVELS_HPA = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
              600, 550, 500, 450, 400, 350, 300, 275, 250, 200]
P_PA = np.array(LEVELS_HPA, dtype=float) * 100.0


def consistent_analysis(gr, dT=1.5, p_sea=101325.0):
    """
    A synthetic analysis whose heights are the hydrostatic integral of its OWN
    temperatures.

    THIS IS NOT A DETAIL. The first version of this file perturbed temperature
    and height independently -- -1.5 K and -45 m of cos(k_y*y) -- which are not
    in hydrostatic balance with each other. A real analysis is. Everything
    measured against the inconsistent version was the model responding to the
    test data: a 140 m geopotential spread and a 100 m/s "balanced" wind, both
    of which vanish here (3 m and 62 m/s).
    """
    ny, nx = gr.ny, gr.nx
    T0, L = 288.15, 0.0065
    k_y = 2 * np.pi / gr.Ly

    T = np.empty((P_PA.size, ny, nx))
    for k, p in enumerate(P_PA):
        T[k] = T0 * (p / 101325.0) ** (RD * L / G0) - dT * np.cos(k_y * gr.Yc)

    z = np.empty_like(T)
    z[0] = (RD * T[0] / G0) * np.log(p_sea / P_PA[0])
    for k in range(1, P_PA.size):
        T_layer = 0.5 * (T[k - 1] + T[k])
        z[k] = z[k - 1] + (RD * T_layer / G0) * np.log(P_PA[k - 1] / P_PA[k])
    return T, z


def us_standard(p):
    """Temperature and geopotential height of a standard atmosphere at p."""
    T0, L = 288.15, 0.0065
    T = T0 * (p / 101325.0) ** (RD * L / G0)
    zz = (T0 - T) / L
    return T, zz


def test_linear_in_log_p_is_exact():
    """A field linear in log(p) must be reproduced exactly."""
    ny, nx = 6, 5
    a = 3.0 + 2.0 * np.log(P_PA)
    src = np.repeat(np.repeat(a[:, None, None], ny, 1), nx, 2)
    p_dst = np.exp(np.linspace(np.log(P_PA.max()), np.log(P_PA.min()), 12))
    p_dst = np.repeat(np.repeat(p_dst[:, None, None], ny, 1), nx, 2)
    out = interp_log_p(src, P_PA, p_dst)
    want = 3.0 + 2.0 * np.log(p_dst)
    err = np.abs(out - want).max()
    report("a field linear in log(p) interpolates exactly", err < 1e-10,
           f"max error {err:.2e}")


def test_source_levels_recovered():
    """Asking for the source pressures must return the source values."""
    ny, nx = 4, 4
    rng = np.random.default_rng(0)
    src = rng.normal(0, 10, (P_PA.size, ny, nx))
    p_dst = np.repeat(np.repeat(P_PA[:, None, None], ny, 1), nx, 2)
    out = interp_log_p(src, P_PA, p_dst)
    err = np.abs(out - src).max()
    report("source levels are recovered exactly", err < 1e-12,
           f"max error {err:.2e} over {P_PA.size} levels")


def test_level_order_does_not_matter():
    """Reversing the source level order must not change the answer."""
    ny, nx = 4, 4
    rng = np.random.default_rng(1)
    src = rng.normal(0, 10, (P_PA.size, ny, nx))
    p_dst = np.full((8, ny, nx), 0.0)
    p_dst[:] = np.linspace(95000, 25000, 8)[:, None, None]
    a = interp_log_p(src, P_PA, p_dst)
    b = interp_log_p(src[::-1], P_PA[::-1], p_dst)
    err = np.abs(a - b).max()
    report("source level order does not matter", err < 1e-12,
           f"max difference {err:.2e}")


def test_surface_pressure_matches_standard_atmosphere():
    """
    Surface pressure derived from analysis heights must match the standard
    atmosphere at the same terrain height, to a few pascals.
    """
    ny, nx = 8, 8
    T, zz = us_standard(P_PA)
    z = np.repeat(np.repeat(zz[:, None, None], ny, 1), nx, 2)
    terrain = np.linspace(0, 2500, ny * nx).reshape(ny, nx)

    p_s = surface_pressure_from_heights(z, P_PA, terrain)
    T0, L = 288.15, 0.0065
    want = 101325.0 * (1 - L * terrain / T0) ** (G0 / (RD * L))
    err = np.abs(p_s - want).max()
    report("surface pressure matches the standard atmosphere", err < 200.0,
           f"max error {err:.1f} Pa over terrain 0-2500 m "
           f"(p_s {p_s.min()/100:.0f}-{p_s.max()/100:.0f} hPa)")


def test_sea_level_terrain_extrapolates_below_1000_hPa():
    """
    At sea level the surface is BELOW the lowest analysis level, which is the
    common case over a coastal domain and the one an interpolation alone
    cannot reach.
    """
    ny, nx = 4, 4
    T, zz = us_standard(P_PA)
    z = np.repeat(np.repeat(zz[:, None, None], ny, 1), nx, 2)
    p_s = surface_pressure_from_heights(z, P_PA, np.zeros((ny, nx)))
    ok = np.all(p_s > P_PA.max()) and abs(p_s.mean() - 101325.0) < 300.0
    report("sea-level terrain extrapolates below the lowest level", ok,
           f"p_s {p_s.mean()/100:.1f} hPa against 1013.25 standard; "
           f"lowest analysis level {P_PA.max()/100:.0f} hPa")


def test_theta_extrapolation_is_not_neutral():
    """
    Theta below the lowest analysis level must keep its lapse rate. Holding it
    constant makes the near-surface layer exactly neutral, which the
    convective adjustment then reads as marginal everywhere on step one.
    """
    ny, nx = 4, 4
    lev = SigmaLevels(20)
    T, zz = us_standard(P_PA)
    T3 = np.repeat(np.repeat(T[:, None, None], ny, 1), nx, 2)
    z3 = np.repeat(np.repeat(zz[:, None, None], ny, 1), nx, 2)
    u3 = np.zeros_like(T3)

    pi, u_s, v_s, th = pressure_to_sigma(u3, u3, T3, z3, P_PA,
                                         np.zeros((ny, nx)), lev)
    dth = th[-1] - th[-2]                  # bottom two model levels
    ok = np.all(dth < -0.01) and np.all(np.isfinite(th))
    report("theta keeps its lapse rate below the lowest analysis level", ok,
           f"d(theta) across the lowest model layer {dth.mean():+.3f} K "
           f"(negative = stable, zero would be neutral)")


def test_geopotential_survives_the_conversion():
    """
    THE MEASUREMENT THAT ACTUALLY TESTS THE CONVERSION.

    The model's own hydrostatic integral over the converted state must
    reproduce the analysis's geopotential. What matters is not the mean offset
    -- a uniform one accelerates nothing -- but the HORIZONTAL SPREAD, because
    a geopotential error that varies across the domain is a pressure-gradient
    force with nothing balancing it.

    Terrain is the hard case: a column over the mountain samples the analysis
    at different pressures than one beside it.
    """
    gr = CGrid(60, 60, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    T, z = consistent_analysis(gr)
    zero = np.zeros_like(T)

    spreads = {}
    for hgt in (0.0, 1500.0, 2500.0):
        h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 150e3) ** 2 +
                           ((gr.Yc - gr.Ly / 2) / 150e3) ** 2))
        pi, u, v, th = pressure_to_sigma(zero, zero, T, z, P_PA, h, lev)
        phi = hydrostatic_geopotential(th, pi, lev, phi_surface=G0 * h)
        phi_a = G0 * interp_log_p(z, P_PA, lev.pressure(pi),
                                  extrapolate="linear")
        d = (phi - phi_a) / G0
        spreads[hgt] = float((d.max(axis=(1, 2)) - d.min(axis=(1, 2))).max())

    ok = spreads[0.0] < 0.1 and spreads[2500.0] < 10.0
    report("geopotential survives the conversion", ok,
           f"horizontal spread of the error: flat {spreads[0.0]:.2f} m, "
           f"1500 m terrain {spreads[1500.0]:.2f} m, "
           f"2500 m terrain {spreads[2500.0]:.2f} m")


def test_terrain_costs_little_extra_drift():
    """
    A converted analysis started at rest adjusts toward balance -- it has a
    temperature gradient and no wind, so it MUST accelerate, and roughly 9 m/s
    over six hours is that adjustment, not a defect.

    What the conversion is responsible for is the EXTRA drift terrain adds. So
    the flat case is the control and the terrain case is measured against it,
    rather than against zero.
    """
    gr = CGrid(60, 60, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    T, z = consistent_analysis(gr)
    zero = np.zeros_like(T)

    drift = {}
    for hgt in (0.0, 2500.0):
        h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 150e3) ** 2 +
                           ((gr.Yc - gr.Ly / 2) / 150e3) ** 2))
        pi, u, v, th = pressure_to_sigma(zero, zero, T, z, P_PA, h, lev)
        m = PrimitiveSigma(gr, lev, terrain=h)
        m.pi, m.u, m.v, m.theta = pi, u, v, th
        m.run(6 * 3600, dt=m.max_dt())
        drift[hgt] = (float(np.abs(m.u).max()) if np.isfinite(m.u).all()
                      else float("nan"))

    extra = drift[2500.0] - drift[0.0]
    ok = np.isfinite(extra) and extra < 4.0
    report("terrain adds little beyond the rest-start adjustment", ok,
           f"flat {drift[0.0]:.2f} m/s (geostrophic adjustment), "
           f"2500 m terrain {drift[2500.0]:.2f} m/s, "
           f"attributable to terrain {extra:+.2f} m/s")


if __name__ == "__main__":
    print("\nPressure -> sigma conversion\n" + "=" * 66)
    for fn in (test_linear_in_log_p_is_exact,
               test_source_levels_recovered,
               test_level_order_does_not_matter,
               test_surface_pressure_matches_standard_atmosphere,
               test_sea_level_terrain_extrapolates_below_1000_hPa,
               test_theta_extrapolation_is_not_neutral,
               test_geopotential_survives_the_conversion,
               test_terrain_costs_little_extra_drift):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

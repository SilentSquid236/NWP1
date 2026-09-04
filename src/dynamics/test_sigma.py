"""
Validation for the sigma vertical coordinate.

The critical tests are the two that pressure coordinates could not pass:
sigma_dot vanishing at both boundaries WITHOUT a correction, and a resting
atmosphere staying at rest over sloping terrain (where the sigma
pressure-gradient force splits into two large terms that must cancel).

Run:  python test_sigma.py
"""

import numpy as np

from grid import CGrid
from sigma import (SigmaLevels, hydrostatic_geopotential, continuity,
                   vertical_advection, pressure_gradient_force,
                   RD, G0, P0, KAPPA)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# ---------------------------------------------------------------------------
def test_level_structure():
    """Sigma must run 0 -> 1 monotonically, with resolution packed near ground."""
    lev = SigmaLevels(20)
    top_dz = lev.dsigma[0]
    bot_dz = lev.dsigma[-1]
    ok = (lev.sigma_half[0] == 0.0 and abs(lev.sigma_half[-1] - 1.0) < 1e-12
          and np.all(lev.dsigma > 0) and bot_dz > top_dz)
    report("sigma levels span 0-1, packed toward the surface", ok,
           f"{lev.nz} levels; d(sigma) {top_dz:.4f} at lid -> {bot_dz:.4f} at "
           f"ground ({bot_dz/top_dz:.1f}x finer below)")


# ---------------------------------------------------------------------------
def test_hydrostatic_exact_isothermal():
    """
    Isothermal atmosphere has the closed form Phi = R T ln(p_s / p).
    Tests the integration independently of everything else.
    """
    lev = SigmaLevels(20)
    pi = np.full((3, 3), 95000.0)          # p_s = 100000, p_top = 5000
    T0 = 250.0
    p = lev.pressure(pi)
    theta = T0 / (p / P0) ** KAPPA

    phi = hydrostatic_geopotential(theta, pi, lev)
    exact = RD * T0 * np.log((lev.p_top + pi) / p)
    err = np.abs(phi - exact).max()

    ok = err < 1e-6
    report("hydrostatic integration exact for isothermal atmosphere", ok,
           f"max error {err:.2e} m2/s2 ({err/G0:.2e} m); "
           f"model top at {phi[0,1,1]/G0:.0f} m")


# ---------------------------------------------------------------------------
def test_sigma_dot_zero_at_boundaries():
    """
    THE test pressure coordinates failed.

    For ARBITRARY winds, sigma_dot must vanish at sigma = 0 and sigma = 1 with
    NO correction applied. In the pressure-coordinate version this required a
    linear correction whose residual fed a feedback loop and destroyed the
    forecast. Here it falls out of the formulation.
    """
    gr = CGrid(32, 30, 12e3, 12e3)
    lev = SigmaLevels(20)
    rng = np.random.default_rng(0)

    pi = 95000.0 + rng.normal(0, 500, (gr.ny, gr.nx))
    u = rng.normal(0, 15, (lev.nz, gr.ny, gr.nx))
    v = rng.normal(0, 15, (lev.nz, gr.ny, gr.nx))

    dpi, sd = continuity(u, v, pi, lev, gr)

    top = np.abs(sd[0]).max()
    bot = np.abs(sd[-1]).max()
    interior = np.abs(sd).max()

    ok = top < 1e-15 and bot < 1e-15
    report("sigma_dot vanishes at lid and ground with NO correction", ok,
           f"|sigma_dot| top {top:.2e}, ground {bot:.2e}, "
           f"interior max {interior:.2e} 1/s")


# ---------------------------------------------------------------------------
def test_mass_conservation():
    """
    Column mass is pi integrated over the domain. With no flux through the
    lateral boundaries (periodic), the total must not change.
    """
    gr = CGrid(32, 32, 12e3, 12e3, edge_mode="periodic")
    lev = SigmaLevels(20)
    rng = np.random.default_rng(1)

    pi = 95000.0 + 300 * np.sin(2 * np.pi * gr.Xc / gr.Lx)
    u = 20.0 + rng.normal(0, 2, (lev.nz, gr.ny, gr.nx))
    v = rng.normal(0, 2, (lev.nz, gr.ny, gr.nx))

    dpi, _ = continuity(u, v, pi, lev, gr)
    total_tendency = dpi.sum()

    ok = abs(total_tendency) < 1e-6 * abs(pi.sum())
    report("total column mass tendency sums to zero (periodic)", ok,
           f"sum(d(pi)/dt) = {total_tendency:.3e} vs total mass "
           f"{pi.sum():.3e} Pa")


# ---------------------------------------------------------------------------
def test_vertical_advection_of_constant():
    """
    Transporting a constant field must produce exactly zero tendency, whatever
    sigma_dot is doing. Flux form guarantees this; advective form does not,
    and the residual is worst near the boundaries.
    """
    gr = CGrid(16, 16, 12e3, 12e3)
    lev = SigmaLevels(20)
    rng = np.random.default_rng(2)

    pi = np.full((gr.ny, gr.nx), 95000.0)
    u = rng.normal(0, 20, (lev.nz, gr.ny, gr.nx))
    v = rng.normal(0, 20, (lev.nz, gr.ny, gr.nx))
    _, sd = continuity(u, v, pi, lev, gr)

    const = np.full((lev.nz, gr.ny, gr.nx), 300.0)
    tend = vertical_advection(const, sd, lev)

    ok = np.abs(tend).max() < 1e-9
    report("vertical advection of a constant field is exactly zero", ok,
           f"max tendency {np.abs(tend).max():.2e} K/s with "
           f"max|sigma_dot| {np.abs(sd).max():.2e} 1/s")


# ---------------------------------------------------------------------------
def test_pressure_gradient_over_terrain():
    """
    THE sigma-coordinate weakness, tested directly.

    An isothermal atmosphere at rest over a SLOPE must stay at rest. On a
    sloping sigma surface the pressure-gradient force splits into two large
    terms -- grad(Phi) along sigma, and R T grad(ln p_s) -- which nearly
    cancel. Their difference is the physical force. Over steep terrain the
    cancellation is the dominant error source in sigma models.

    Here we check the two terms cancel to a small residual acceleration.
    """
    gr = CGrid(48, 48, 12e3, 12e3)
    lev = SigmaLevels(20)

    # Terrain: a ridge rising 1500 m across the domain.
    h = 1500.0 * np.exp(-((gr.Yc - gr.Ly / 2) / 200e3) ** 2)
    phi_s = G0 * h

    # Isothermal atmosphere in exact hydrostatic balance above that terrain.
    T0 = 260.0
    p_sea = 101325.0
    p_s = p_sea * np.exp(-G0 * h / (RD * T0))
    pi = p_s - lev.p_top

    p = lev.pressure(pi)
    theta = T0 / (p / P0) ** KAPPA
    phi = hydrostatic_geopotential(theta, pi, lev, phi_surface=phi_s)

    # Hydrostatic consistency: for T = T(p) the true force is EXACTLY zero.
    T = theta * (p / P0) ** KAPPA
    term1 = -gr.dy_backward(phi)
    term2 = -gr.h_to_v(RD * T) * gr.dy_backward(np.log(p))
    _, net = pressure_gradient_force(phi, theta, pi, lev, gr)

    big = max(np.abs(term1).max(), np.abs(term2).max())
    resid = float(np.abs(net).max())
    ok = resid < 1e-12
    report("PGF is hydrostatically consistent over a 1500 m ridge", ok,
           f"individual terms up to {big:.3f} m/s2, net residual "
           f"{resid:.2e} m/s2 -- must be machine zero, not merely small")


# ---------------------------------------------------------------------------
def test_terrain_following_surface():
    """
    The ground is sigma = 1 regardless of terrain height, and the pressure
    there equals surface pressure exactly. In pressure coordinates this is
    the condition that could not be represented at all.
    """
    lev = SigmaLevels(20)
    h = np.array([[0.0, 500.0], [1500.0, 3000.0]])
    p_s = 101325.0 * np.exp(-G0 * h / (RD * 260.0))
    pi = p_s - lev.p_top

    p_half = lev.pressure_half(pi)
    ground = p_half[-1]
    lid = p_half[0]

    ok = (np.abs(ground - p_s).max() < 1e-9
          and np.abs(lid - lev.p_top).max() < 1e-9)
    report("ground is sigma=1 at every terrain height", ok,
           f"terrain 0-3000 m -> p_s {p_s.min()/100:.0f}-{p_s.max()/100:.0f} hPa; "
           f"p at sigma=1 matches p_s exactly, lid fixed at "
           f"{lev.p_top/100:.0f} hPa")


if __name__ == "__main__":
    print("\nSigma coordinate validation\n" + "=" * 64)
    for fn in (test_level_structure,
               test_hydrostatic_exact_isothermal,
               test_sigma_dot_zero_at_boundaries,
               test_mass_conservation,
               test_vertical_advection_of_constant,
               test_pressure_gradient_over_terrain,
               test_terrain_following_surface):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 64)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

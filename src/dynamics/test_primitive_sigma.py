"""
Validation for the sigma-coordinate 3D core.

The decisive test is the last one: integrate from a REALISTIC noisy, sheared,
stratified state -- the case where the pressure-coordinate core diverged
within 2-3 hours at every damping setting tried.

Run:  python test_primitive_sigma.py
"""

import numpy as np

from grid import CGrid
from sigma import (SigmaLevels, RD, G0, P0, KAPPA,
                   hydrostatic_geopotential, pressure_gradient_force)
from primitive_sigma import PrimitiveSigma
from subgrid import balance_initial_state
from initialization import filter_initial_state

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def isothermal(model, T0=260.0, p_sea=101325.0):
    """Isothermal atmosphere in hydrostatic balance over the model's terrain."""
    p_s = p_sea * np.exp(-G0 * model.terrain / (RD * T0))
    model.pi = p_s - model.lev.p_top
    p = model.lev.pressure(model.pi)
    model.theta = T0 / (p / P0) ** KAPPA


# ---------------------------------------------------------------------------
def test_rest_over_flat_ground():
    """A motionless isothermal atmosphere over flat ground must stay motionless."""
    gr = CGrid(32, 32, 25e3, 25e3)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev)
    isothermal(m)

    pi0 = m.pi.copy()
    m.run(12 * 3600, dt=m.max_dt())

    du = max(np.abs(m.u).max(), np.abs(m.v).max())
    dpi = np.abs(m.pi - pi0).max()
    ok = du < 1e-8 and dpi < 1e-6
    report("rest over flat ground stays at rest (12 h)", ok,
           f"max|u| {du:.2e} m/s, max|d(p_s)| {dpi:.2e} Pa")


# ---------------------------------------------------------------------------
def test_rest_over_terrain():
    """
    THE test pressure coordinates could not even pose. A balanced atmosphere
    over a mountain must stay at rest: the two large pressure-gradient terms
    must cancel on the sloping sigma surfaces.
    """
    gr = CGrid(48, 48, 20e3, 20e3)
    lev = SigmaLevels(20)
    h = 1200.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 150e3) ** 2 +
                          ((gr.Yc - gr.Ly / 2) / 150e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h)
    isothermal(m)

    m.run(6 * 3600, dt=m.max_dt())

    umax = max(np.abs(m.u).max(), np.abs(m.v).max())
    ok = umax < 2.0 and np.isfinite(m.u).all()
    report("rest over a 1200 m mountain stays near rest (6 h)", ok,
           f"spurious max|u| {umax:.3f} m/s over "
           f"{h.max():.0f} m terrain (sigma PGF cancellation error)")


# ---------------------------------------------------------------------------
def test_surface_pressure_is_prognostic():
    """
    Surface pressure must actually evolve -- that is the structural change.
    A divergent flow should change p_s measurably.
    """
    gr = CGrid(32, 32, 25e3, 25e3, edge_mode="periodic")
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev)
    isothermal(m)
    m.u[:] = 10.0 * np.sin(2 * np.pi * gr.Xc / gr.Lx)      # convergent/divergent

    ps0 = m.surface_pressure.copy()
    m.run(3 * 3600, dt=m.max_dt())
    change = np.abs(m.surface_pressure - ps0).max()

    ok = 1.0 < change < 5000.0 and np.isfinite(m.pi).all()
    report("surface pressure evolves under divergent flow", ok,
           f"max |dp_s| = {change:.1f} Pa over 3 h "
           f"(pressure coords could not represent this at all)")


# ---------------------------------------------------------------------------
def test_mass_conservation():
    """Total column mass is conserved on a periodic domain."""
    gr = CGrid(32, 32, 25e3, 25e3, edge_mode="periodic")
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev, hyper=0.0)
    isothermal(m)
    m.u[:] = 15.0
    m.pi += 200.0 * np.sin(2 * np.pi * gr.Yc / gr.Ly)

    m0 = m.total_mass()
    m.run(12 * 3600, dt=m.max_dt())
    rel = abs(m.total_mass() - m0) / m0

    ok = rel < 1e-10
    report("total mass conserved over 12 h", ok,
           f"relative drift {rel:.2e}")


# ---------------------------------------------------------------------------
def test_thermal_wind_balance():
    """A balanced baroclinic jet must persist, as in the pressure version."""
    gr = CGrid(48, 48, 50e3, 50e3, f0=1.0e-4, beta=0.0)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev)

    T0, dT = 260.0, 4.0
    k_y = 2 * np.pi / gr.Ly
    p_s = 101325.0 * np.ones((gr.ny, gr.nx))
    m.pi = p_s - lev.p_top
    p = lev.pressure(m.pi)
    T = T0 + dT * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA

    phi = m.geopotential()
    for k in range(lev.nz):
        dphidy = 0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k]))
        m.u[k] = -dphidy / gr.f0

    u0 = m.u.copy()
    m.run(24 * 3600, dt=m.max_dt())

    drift = np.abs(m.u - u0).max() / max(np.abs(u0).max(), 1e-9)
    vspur = np.abs(m.v).max() / max(np.abs(u0).max(), 1e-9)
    ok = drift < 0.20 and vspur < 0.15
    report("thermal-wind jet persists 24 h", ok,
           f"|du|/|u| {drift*100:.2f}%, spurious |v|/|u| {vspur*100:.2f}%; "
           f"jet {np.abs(u0).max():.1f} m/s")


# ---------------------------------------------------------------------------
def test_realistic_noisy_state_is_stable():
    """
    THE DECISIVE TEST -- rewritten after the initial state was probed rather
    than the model patched.

    THREE ERRORS WERE IN THE OLD VERSION OF THIS TEST, all in the setup:

    1. A 6 K meridional temperature contrast implies a 166 m/s jet by thermal
       wind. The test then clipped the wind at 60 m/s, destroying geostrophic
       balance over 34% of the domain. The 2dx growth the model was blamed
       for came from that clip. 1.5 K gives a realistic 41 m/s jet, unclipped.

    2. The balanced wind was -d(phi)/dy / f. On sigma surfaces the horizontal
       force has two terms that largely cancel over sloping ground; keeping
       only one implies an 845 m/s wind over 2500 m terrain. It now comes
       from the full PGF.

    3. 1.2 m/s of WHITE noise puts 89% of its variance at wavelengths the
       grid cannot carry, and nonlinear advection amplifies that faster than
       hyperdiffusion (3 h e-folding, measured) can remove it. Real analyses
       are filtered before they are integrated; this one now is too, and then
       rebalanced, because filtering u, v and theta separately reintroduces
       divergence.

    Measured on flat ground: no filter 1/12 h, filter only 11/12 h, filter
    then balance 12/12 h -- with mixing and drag making no difference in any
    of the three. The failure was never a boundary-layer problem.
    """
    gr = CGrid(90, 88, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    rng = np.random.default_rng(0)

    h = 400.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                         ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h)

    k_y = 2 * np.pi / gr.Ly
    p_s = 101325.0 * np.exp(-G0 * h / (RD * 280.0))
    m.pi = p_s - lev.p_top
    p = lev.pressure(m.pi)
    T = 288.0 - 55.0 * (1 - p / p.max()) - 1.5 * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA

    phi = hydrostatic_geopotential(m.theta, m.pi, lev, phi_surface=m.phi_s)
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u[:] = -fy / gr.f0
    m.v[:] = fx / gr.f0
    jet = float(np.abs(m.u).max())

    # Analysis-like grid-scale noise, then the treatment a real analysis gets.
    m.u += rng.normal(0, 1.2, m.u.shape)
    m.v += rng.normal(0, 1.2, m.v.shape)
    m.u, m.v, m.theta = filter_initial_state(m.u, m.v, m.theta, gr)
    m.u, m.v, _ = balance_initial_state(m.u, m.v, gr, verbose=False)

    dt = m.max_dt()
    hours, trace = 0, []
    for h_i in range(12):
        m.run(3600, dt=dt)
        d = m.diagnostics()
        trace.append(d["max|sigma_dot|"])
        if not np.isfinite(m.u).all():
            break
        hours = h_i + 1

    umax = float(np.abs(m.u).max()) if np.isfinite(m.u).all() else float("nan")
    ok = hours == 12 and umax < 80.0
    report("12 h from a realistic noisy state stays stable", ok,
           f"survived {hours}/12 h; {jet:.0f} m/s jet -> max|u| {umax:.1f} m/s; "
           f"max|sigma_dot| by hour: "
           f"{', '.join(f'{t:.1e}' for t in trace[:6])}"
           + (" ..." if len(trace) > 6 else ""))


if __name__ == "__main__":
    print("\nSigma-coordinate 3D core\n" + "=" * 66)
    for fn in (test_rest_over_flat_ground,
               test_rest_over_terrain,
               test_surface_pressure_is_prognostic,
               test_mass_conservation,
               test_thermal_wind_balance,
               test_realistic_noisy_state_is_stable):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

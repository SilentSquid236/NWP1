"""
Analytic validation for the dry hydrostatic primitive-equation core.

Run:  python test_primitive3d.py
"""

import numpy as np

from grid import CGrid
from vertical import (PressureLevels, hydrostatic_geopotential, diagnose_omega,
                      theta_from_T, T_from_theta, RD, G0)
from primitive3d import Primitive3D

LEVELS = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 750,
          700, 650, 600, 550, 500, 450, 400, 300, 250, 200]

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def isothermal_state(model, T0=250.0):
    for k in range(model.lev.nz):
        model.theta[k] = theta_from_T(T0, model.lev.p[k])


# ---------------------------------------------------------------------------
def test_hydrostatic_exact_isothermal():
    """
    For an isothermal atmosphere the hypsometric relation has a closed form:

        Phi(p) = R T ln(p_s / p)

    The discrete integral must reproduce it exactly, because layer-mean T is
    exact when T is constant.
    """
    lev = PressureLevels(LEVELS)
    T0 = 250.0
    theta = np.empty((lev.nz, 3, 3))
    for k in range(lev.nz):
        theta[k] = theta_from_T(T0, lev.p[k])

    phi = hydrostatic_geopotential(theta, lev)[:, 1, 1]
    exact = RD * T0 * np.log(lev.p_bottom / lev.p)
    err = np.abs(phi - exact).max()

    ok = err < 1e-6
    report("hydrostatic integration exact for isothermal atmosphere", ok,
           f"max error {err:.2e} m2/s2 ({err/G0:.2e} m); "
           f"top at {phi[-1]/G0:.0f} m")


# ---------------------------------------------------------------------------
def test_rest_stays_at_rest():
    """
    A horizontally uniform, stably stratified, motionless atmosphere is an
    exact steady solution. Any drift means a spurious pressure-gradient or
    advection term.
    """
    gr = CGrid(32, 32, 50e3, 50e3)
    lev = PressureLevels(LEVELS)
    m = Primitive3D(gr, lev)
    isothermal_state(m)

    th0 = m.theta.copy()
    m.run(12 * 3600, dt=m.max_dt())

    du = max(np.abs(m.u).max(), np.abs(m.v).max())
    dth = np.abs(m.theta - th0).max()
    ok = du < 1e-9 and dth < 1e-9
    report("stratified atmosphere at rest stays at rest (12 h)", ok,
           f"max|u| = {du:.2e} m/s, max|dtheta| = {dth:.2e} K")


# ---------------------------------------------------------------------------
def test_omega_zero_for_nondivergent_flow():
    """
    omega is diagnosed from horizontal divergence. A purely zonal flow with no
    x-variation is non-divergent, so omega must vanish -- and must vanish at
    both the lid and the ground by construction.
    """
    gr = CGrid(32, 32, 50e3, 50e3)
    lev = PressureLevels(LEVELS)
    m = Primitive3D(gr, lev)
    isothermal_state(m)
    m.u[:] = 20.0                       # uniform zonal wind

    om = m.omega()
    ok = np.abs(om).max() < 1e-12
    report("omega vanishes for non-divergent flow", ok,
           f"max|omega| = {np.abs(om).max():.2e} Pa/s "
           f"(top {np.abs(om[-1]).max():.1e}, bottom {np.abs(om[0]).max():.1e})")


# ---------------------------------------------------------------------------
def test_omega_boundary_conditions():
    """
    For an ARBITRARY divergent flow, omega must still be zero at the top and
    the bottom -- rigid lid, flat ground. This is what the linear correction
    in diagnose_omega() enforces.
    """
    gr = CGrid(32, 32, 50e3, 50e3)
    lev = PressureLevels(LEVELS)
    rng = np.random.default_rng(0)
    div = rng.normal(0, 1e-5, (lev.nz, gr.ny, gr.nx))

    om = diagnose_omega(div, lev)
    top = np.abs(om[-1]).max()
    bot = np.abs(om[0]).max()

    ok = top < 1e-12 and bot < 1e-12
    report("omega = 0 at lid and ground for arbitrary divergence", ok,
           f"|omega| top = {top:.2e}, bottom = {bot:.2e}, "
           f"interior max = {np.abs(om).max():.2e} Pa/s")


# ---------------------------------------------------------------------------
def test_thermal_wind_balance():
    """
    THE 3D balance test. A meridional temperature gradient requires vertical
    wind shear -- thermal wind:

        f du/dp = (R/p) dT/dy

    Build a baroclinic state satisfying it, and the flow should persist. This
    couples the hydrostatic integration, the pressure gradient, and Coriolis;
    an error in any of the three destroys the balance within hours.
    """
    gr = CGrid(48, 48, 50e3, 50e3, f0=1.0e-4, beta=0.0)
    lev = PressureLevels(LEVELS)
    m = Primitive3D(gr, lev)

    # Temperature decreasing poleward, periodic in y.
    T0, dT = 260.0, 8.0
    k_y = 2 * np.pi / gr.Ly
    for k in range(lev.nz):
        T = T0 + dT * np.cos(k_y * gr.Yc)
        m.theta[k] = theta_from_T(T, lev.p[k])

    # u in geostrophic balance with the resulting geopotential:
    #   f u = -dPhi/dy
    phi = m.geopotential()
    for k in range(lev.nz):
        dphidy = 0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k]))
        m.u[k] = -dphidy / gr.f0

    u0 = m.u.copy()
    shear = np.abs(u0[-1] - u0[0]).max()
    m.run(24 * 3600, dt=m.max_dt())

    drift = np.abs(m.u - u0).max() / max(np.abs(u0).max(), 1e-9)
    v_spur = np.abs(m.v).max() / max(np.abs(u0).max(), 1e-9)
    ok = drift < 0.10 and v_spur < 0.05
    report("thermal-wind balanced jet persists 24 h", ok,
           f"|du|/|u| = {drift*100:.2f}%, spurious |v|/|u| = {v_spur*100:.2f}%; "
           f"jet {np.abs(u0).max():.1f} m/s, shear {shear:.1f} m/s over depth")


# ---------------------------------------------------------------------------
def test_geostrophic_imbalance_converges():
    """
    The initial geostrophic state is balanced only to discretisation accuracy,
    so a small spurious dv/dt is expected. What matters is that it CONVERGES
    at second order under grid refinement -- that is the difference between
    truncation error and a bug. A bug would not converge.
    """
    def imbalance(n, dx):
        gr = CGrid(n, n, dx, dx, f0=1.0e-4, beta=0.0)
        lev = PressureLevels(LEVELS)
        m = Primitive3D(gr, lev)
        k_y = 2 * np.pi / gr.Ly
        for k in range(lev.nz):
            m.theta[k] = theta_from_T(260.0 + 8.0 * np.cos(k_y * gr.Yc), lev.p[k])
        phi = m.geopotential()
        for k in range(lev.nz):
            m.u[k] = -0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k])) / gr.f0
        _, dv, _ = m.tendencies(m.u, m.v, m.theta)
        return np.abs(dv).max()

    # Same physical domain, doubling resolution.
    a = imbalance(24, 25e3)
    b = imbalance(48, 12.5e3)
    c = imbalance(96, 6.25e3)
    r1, r2 = a / b, b / c

    ok = 3.5 < r1 < 4.5 and 3.5 < r2 < 4.5
    report("geostrophic imbalance converges at 2nd order", ok,
           f"|dv/dt| {a:.2e} -> {b:.2e} -> {c:.2e}, "
           f"ratios {r1:.2f}, {r2:.2f} (2nd order = 4.0)")


# ---------------------------------------------------------------------------
def test_theta_conservation():
    """
    Dry dynamics has no heat source, so mass-weighted theta is conserved.
    Drift means the advection scheme is creating or destroying heat.
    """
    gr = CGrid(48, 48, 50e3, 50e3)
    lev = PressureLevels(LEVELS)
    m = Primitive3D(gr, lev)
    isothermal_state(m, 260.0)

    # A warm bubble plus flow to move it around.
    blob = 3.0 * np.exp(-(((gr.Xc - gr.Lx/2) / 300e3)**2 +
                          ((gr.Yc - gr.Ly/2) / 300e3)**2))
    for k in range(lev.nz):
        m.theta[k] += blob * np.exp(-k / 6)
    m.u[:] = 15.0

    t0 = m.total_theta()
    m.run(12 * 3600, dt=m.max_dt())
    rel = abs(m.total_theta() - t0) / abs(t0)

    ok = rel < 1e-4 and np.isfinite(m.theta).all()
    report("mass-weighted theta nearly conserved (advective form)", ok,
           f"relative drift = {rel:.2e} over 12 h "
           f"(advective form is not exactly conservative -- see README)")


# ---------------------------------------------------------------------------
def test_static_stability_sign():
    """
    A standard atmosphere must come out statically STABLE: theta increasing
    with height, i.e. dtheta/dp < 0. A sign error here would make the model
    spontaneously convect.
    """
    lev = PressureLevels(LEVELS)
    from vertical import static_stability
    theta = np.empty((lev.nz, 2, 2))
    for k in range(lev.nz):
        T = 288.0 - 6.5e-3 * (RD * 250.0 * np.log(lev.p_bottom/lev.p[k]) / G0)
        theta[k] = theta_from_T(T, lev.p[k])

    dthdp = static_stability(theta, lev)[:, 0, 0]
    ok = np.all(dthdp < 0)
    report("standard atmosphere is statically stable", ok,
           f"dtheta/dp in [{dthdp.min():.2e}, {dthdp.max():.2e}] K/Pa "
           f"(all negative = stable); theta {theta[0,0,0]:.1f} -> "
           f"{theta[-1,0,0]:.1f} K")


if __name__ == "__main__":
    print("\nDry primitive-equation validation\n" + "=" * 62)
    for fn in (test_hydrostatic_exact_isothermal,
               test_rest_stays_at_rest,
               test_omega_zero_for_nondivergent_flow,
               test_omega_boundary_conditions,
               test_thermal_wind_balance,
               test_geostrophic_imbalance_converges,
               test_theta_conservation,
               test_static_stability_sign):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

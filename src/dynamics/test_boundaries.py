"""
Validation for limited-area boundary treatment.

The question a boundary scheme must answer: when a wave reaches the edge, does
it leave, or does it bounce back and pollute the forecast? These tests measure
that directly.

Run:  python test_boundaries.py
"""

import numpy as np

from grid import CGrid
from shallow_water import ShallowWaterModel, G
from boundaries import (relaxation_weights, DaviesRelaxation, BoundaryDriver,
                        run_limited_area)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# ---------------------------------------------------------------------------
def test_weight_profile():
    """alpha must be ~1 at the perimeter, exactly 0 in the interior, smooth between."""
    w = relaxation_weights(80, 80, width=10)
    edge = w[0, 40]
    interior = w[40, 40]
    inner_edge = w[10, 40]
    monotonic = np.all(np.diff(w[:11, 40]) <= 1e-12)

    ok = edge > 0.99 and interior == 0.0 and inner_edge < 1e-12 and monotonic
    report("relaxation weights taper correctly", ok,
           f"edge {edge:.3f}, zone-edge {inner_edge:.2e}, interior {interior:.1f}, "
           f"monotonic {monotonic}")


# ---------------------------------------------------------------------------
def test_steady_state_preserved():
    """
    If the external state equals the model state and both are steady, nothing
    should change. Catches sign errors in the relaxation itself.
    """
    gr = CGrid(64, 64, 20e3, 20e3, edge_mode="replicate")
    m = ShallowWaterModel(gr, H=10_000.0)
    relax = DaviesRelaxation(gr, width=8)
    ext = {"u": m.u.copy(), "v": m.v.copy(), "h": m.h.copy()}
    driver = BoundaryDriver([0.0], [ext])

    run_limited_area(m, driver, relax, 6 * 3600)

    err = max(np.abs(m.u).max(), np.abs(m.v).max(), np.abs(m.h - m.H).max())
    ok = err < 1e-9
    report("steady state preserved under relaxation", ok,
           f"max deviation = {err:.2e}")


# ---------------------------------------------------------------------------
def test_relaxation_pulls_toward_external():
    """The boundary zone must actually track the external state."""
    gr = CGrid(64, 64, 20e3, 20e3, edge_mode="replicate")
    m = ShallowWaterModel(gr, H=10_000.0)
    relax = DaviesRelaxation(gr, width=8)

    target = {"u": np.full_like(m.u, 5.0),
              "v": np.zeros_like(m.v),
              "h": np.full_like(m.h, 10_000.0)}
    driver = BoundaryDriver([0.0], [target])

    run_limited_area(m, driver, relax, 3 * 3600)

    edge_u = m.u[0, :].mean()
    interior_u = m.u[32, 32]
    ok = abs(edge_u - 5.0) < 0.05
    report("boundary zone tracks external state", ok,
           f"edge u = {edge_u:.3f} m/s (target 5.0), interior u = {interior_u:.3f}")


# ---------------------------------------------------------------------------
def test_wave_exits_without_reflecting():
    """
    THE test. A gravity-wave pulse propagates outward and hits the boundary.

    With relaxation, the wave should leave the domain. Without it (a rigid
    replicate edge), it reflects and the energy stays inside. We compare
    interior energy after the wave has had time to reach the edge and exit.
    """
    def interior_energy(use_relaxation):
        gr = CGrid(120, 120, 10e3, 10e3, f0=0.0, beta=0.0,
                   edge_mode="replicate")
        H = 10_000.0
        m = ShallowWaterModel(gr, H=H)
        m.h += 20.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 60e3) ** 2 +
                               ((gr.Yc - gr.Ly / 2) / 60e3) ** 2))

        # Time for the pulse to cross half the domain and exit.
        c = m.gravity_wave_speed
        t_run = 1.6 * (gr.Lx / 2) / c

        relax = DaviesRelaxation(gr, width=12)
        quiet = {"u": np.zeros_like(m.u), "v": np.zeros_like(m.v),
                 "h": np.full_like(m.h, H)}
        driver = BoundaryDriver([0.0], [quiet])

        if use_relaxation:
            run_limited_area(m, driver, relax, t_run)
        else:
            m.run(t_run)

        core = slice(30, 90)
        dh = m.h[core, core] - H
        return float((dh ** 2).sum())

    with_relax = interior_energy(True)
    without = interior_energy(False)
    ratio = with_relax / without if without > 0 else np.inf

    ok = ratio < 0.2
    report("outgoing wave exits instead of reflecting", ok,
           f"residual interior signal: relaxed {with_relax:.3e}, "
           f"rigid {without:.3e} -> {ratio:.1%} of the reflecting case")


# ---------------------------------------------------------------------------
def test_time_interpolation():
    """
    Boundary data arrives hourly but the model steps every few seconds.
    Interpolation must be linear in time and exact at the frame times.
    """
    a = {"h": np.zeros((4, 4))}
    b = {"h": np.full((4, 4), 100.0)}
    d = BoundaryDriver([0.0, 3600.0], [a, b])

    mid = d.at(1800.0)["h"][0, 0]
    at0 = d.at(0.0)["h"][0, 0]
    at1 = d.at(3600.0)["h"][0, 0]
    before = d.at(-500.0)["h"][0, 0]
    after = d.at(9999.0)["h"][0, 0]

    ok = (abs(mid - 50.0) < 1e-12 and at0 == 0.0 and at1 == 100.0
          and before == 0.0 and after == 100.0)
    report("boundary data interpolates linearly in time", ok,
           f"t=0 {at0:.1f}, t=1800 {mid:.1f} (want 50), t=3600 {at1:.1f}, "
           f"clamped outside range")


# ---------------------------------------------------------------------------
def test_interior_still_conserves_mass():
    """
    Relaxation deliberately adds and removes mass at the boundary -- that is
    what an open boundary means. But with a quiescent external state the
    interior should stay close to its initial mass rather than drifting.
    """
    gr = CGrid(80, 80, 20e3, 20e3, edge_mode="replicate")
    H = 10_000.0
    m = ShallowWaterModel(gr, H=H)
    m.h += 10.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 100e3) ** 2 +
                           ((gr.Yc - gr.Ly / 2) / 100e3) ** 2))
    relax = DaviesRelaxation(gr, width=10)
    quiet = {"u": np.zeros_like(m.u), "v": np.zeros_like(m.v),
             "h": np.full_like(m.h, H)}
    driver = BoundaryDriver([0.0], [quiet])

    core = slice(20, 60)
    m0 = float(m.h[core, core].sum())
    run_limited_area(m, driver, relax, 12 * 3600)
    m1 = float(m.h[core, core].sum())
    rel = abs(m1 - m0) / m0

    ok = rel < 1e-3 and np.isfinite(m.h).all()
    report("interior mass stable over 12 h with open boundaries", ok,
           f"relative change = {rel:.2e}")


if __name__ == "__main__":
    print("\nLimited-area boundary validation\n" + "=" * 60)
    for fn in (test_weight_profile,
               test_steady_state_preserved,
               test_relaxation_pulls_toward_external,
               test_wave_exits_without_reflecting,
               test_time_interpolation,
               test_interior_still_conserves_mass):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 60)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

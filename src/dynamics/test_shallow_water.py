"""
Analytic validation for the shallow-water core.

These are not unit tests of code paths -- they are physics tests with known
answers. A dynamical core that passes them has correct advection, Coriolis,
pressure gradient and mass conservation. One that fails them is wrong in a way
that will be invisible and unfixable once moisture and 3D are added.

Run:  python test_shallow_water.py
"""

import numpy as np

from grid import CGrid
from shallow_water import ShallowWaterModel, G

PASS, FAIL = "PASS", "FAIL"
results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{PASS if ok else FAIL}] {name}\n        {detail}")


# ---------------------------------------------------------------------------
def test_rest_stays_at_rest():
    """
    A motionless atmosphere of uniform depth must remain motionless forever.

    This catches sign errors and spurious pressure-gradient forcing. It should
    hold to machine precision -- every term is identically zero.
    """
    gr = CGrid(64, 64, 10e3, 10e3)
    m = ShallowWaterModel(gr, H=10_000.0)
    m.run(6 * 3600, dt=m.max_dt())

    max_u = max(np.abs(m.u).max(), np.abs(m.v).max())
    dh = np.abs(m.h - m.H).max()
    ok = max_u < 1e-12 and dh < 1e-9
    report("resting atmosphere stays at rest (6 h)", ok,
           f"max|u| = {max_u:.2e} m/s, max|h - H| = {dh:.2e} m")


# ---------------------------------------------------------------------------
def test_mass_conservation():
    """
    Flux-form continuity on a periodic domain conserves mass exactly: the
    divergence sums telescope to zero. Drift here means the flux form is wrong.
    """
    gr = CGrid(64, 64, 10e3, 10e3)
    m = ShallowWaterModel(gr, H=10_000.0)

    # A blob, plus a sheared flow to actually move mass around.
    m.h += 200.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 60e3) ** 2 +
                            ((gr.Yc - gr.Ly / 2) / 60e3) ** 2))
    m.u += 20.0 * np.sin(2 * np.pi * gr.Yc / gr.Ly)

    m0 = m.total_mass()
    m.run(12 * 3600, dt=m.max_dt())
    rel = abs(m.total_mass() - m0) / m0

    ok = rel < 1e-12
    report("mass conserved over 12 h with active flow", ok,
           f"relative drift = {rel:.2e}")


# ---------------------------------------------------------------------------
def test_gravity_wave_speed():
    """
    A depth perturbation with no rotation spreads as a gravity wave at
    c = sqrt(gH). Measuring the front position tests the pressure gradient and
    divergence terms against an exact analytic speed.
    """
    H = 1000.0
    c_exact = np.sqrt(G * H)

    gr = CGrid(400, 4, 2e3, 2e3, f0=0.0, beta=0.0)     # no rotation
    m = ShallowWaterModel(gr, H=H)

    x0 = gr.Lx / 2
    m.h += 1.0 * np.exp(-((gr.Xc - x0) / 8e3) ** 2)    # small => linear

    t_run = 600.0
    m.run(t_run, dt=m.max_dt())

    # The initial bump splits into two pulses travelling in opposite
    # directions at +-c. Track the RIGHT-GOING PEAK: a threshold-based "front"
    # is not well defined for a Gaussian, but the peak is.
    row = m.h[gr.ny // 2] - H
    right = row[gr.xc > x0]
    x_right = gr.xc[gr.xc > x0]
    peak_x = x_right[np.argmax(right)]
    c_measured = (peak_x - x0) / t_run

    err = abs(c_measured - c_exact) / c_exact
    ok = err < 0.05
    report("gravity wave speed matches sqrt(gH)", ok,
           f"exact {c_exact:.1f} m/s, measured {c_measured:.1f} m/s, "
           f"error {err * 100:.1f}%")


# ---------------------------------------------------------------------------
def test_geostrophic_balance():
    """
    A jet in geostrophic balance -- Coriolis exactly opposing the pressure
    gradient -- is a steady solution and must persist.

        f u = -g dh/dy

    This is the single most important test for an atmospheric model: nearly all
    large-scale flow is close to this balance. A sign error in Coriolis, or
    mismatched staggering between the two terms, destroys it within hours.
    """
    gr = CGrid(64, 64, 20e3, 20e3, f0=1.0e-4, beta=0.0)
    H = 10_000.0
    m = ShallowWaterModel(gr, H=H)

    # The balanced state must itself be periodic in y, or the wrap-around seam
    # introduces a huge spurious gradient. A tanh jet is NOT periodic; a
    # sinusoid is. This is a property of the domain, not of the physics.
    amp = 100.0
    k = 2 * np.pi / gr.Ly
    m.h = H + amp * np.cos(k * gr.Yc)

    # u = -(g/f) dh/dy, exactly, at the u points (same rows as h).
    m.u = -(G / gr.f0) * (-amp * k * np.sin(k * gr.Yc))
    m.v[:] = 0.0

    u0 = m.u.copy()
    m.run(24 * 3600, dt=m.max_dt())

    drift = np.abs(m.u - u0).max() / np.abs(u0).max()
    v_spur = np.abs(m.v).max()
    ok = drift < 0.05 and v_spur < 1.0
    report("geostrophic jet persists for 24 h", ok,
           f"|du|/|u| = {drift * 100:.2f}%, spurious max|v| = {v_spur:.3f} m/s "
           f"(jet peak {np.abs(u0).max():.1f} m/s)")


# ---------------------------------------------------------------------------
def test_cfl_limit_is_real():
    """
    Stepping past the CFL limit must blow up. If it does not, max_dt() is
    reporting a limit the scheme does not actually have, and long runs would
    fail unpredictably instead of immediately.
    """
    gr = CGrid(64, 64, 10e3, 10e3)
    m = ShallowWaterModel(gr, H=10_000.0)
    m.h += 50.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 50e3) ** 2 +
                           ((gr.Yc - gr.Ly / 2) / 50e3) ** 2))

    dt_max = m.max_dt()

    safe = ShallowWaterModel(gr, H=10_000.0); safe.h = m.h.copy()
    safe.run(3600, dt=dt_max * 0.5)
    safe_ok = np.isfinite(safe.h).all() and np.abs(safe.u).max() < 1e3

    wild = ShallowWaterModel(gr, H=10_000.0); wild.h = m.h.copy()
    with np.errstate(all="ignore"):
        try:
            wild.run(3600, dt=dt_max * 4.0)
            blew_up = (not np.isfinite(wild.h).all()) or np.abs(wild.u).max() > 1e4
        except (FloatingPointError, ValueError):
            blew_up = True

    ok = safe_ok and blew_up
    report("CFL limit is real (stable below, unstable above)", ok,
           f"dt_max = {dt_max:.1f}s; at 0.5x stable = {safe_ok}, "
           f"at 4x unstable = {blew_up}")


# ---------------------------------------------------------------------------
def test_energy_drift():
    """
    Advective-form momentum does not conserve energy exactly, but drift over a
    day should be small. Large drift signals an unstable or badly staggered
    scheme; growth means energy is being created from nothing.
    """
    gr = CGrid(64, 64, 20e3, 20e3)
    m = ShallowWaterModel(gr, H=10_000.0)
    rng = np.random.default_rng(0)
    m.h += 30.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 150e3) ** 2 +
                           ((gr.Yc - gr.Ly / 2) / 150e3) ** 2))

    e0 = m.total_energy()
    m.run(24 * 3600, dt=m.max_dt())
    rel = (m.total_energy() - e0) / e0

    ok = abs(rel) < 0.05
    report("energy drift bounded over 24 h", ok,
           f"relative change = {rel * 100:+.3f}%")


if __name__ == "__main__":
    print("\nShallow-water analytic validation\n" + "=" * 60)
    for fn in (test_rest_stays_at_rest,
               test_mass_conservation,
               test_gravity_wave_speed,
               test_geostrophic_balance,
               test_cfl_limit_is_real,
               test_energy_drift):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 60)
    n_ok = sum(results)
    print(f"{n_ok}/{len(results)} passed\n")
    raise SystemExit(0 if n_ok == len(results) else 1)

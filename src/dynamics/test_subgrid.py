"""
Tests for dissipation and stochastic variance -- the machinery that lets the
model run on real, evolving atmospheres instead of analytic balanced states.

Run:  python test_subgrid.py
"""

import numpy as np

from grid import CGrid
from vertical import PressureLevels, theta_from_T
from primitive3d import Primitive3D
from subgrid import (hyperdiffusion, recommended_hyper_coeff,
                     StochasticPerturbation, perturb_initial_state)

LEVELS = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 750,
          700, 650, 600, 550, 500, 450, 400, 300, 250, 200]

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# ---------------------------------------------------------------------------
def test_hyperdiffusion_is_scale_selective():
    """
    The whole point of biharmonic damping: it must annihilate grid-scale noise
    while barely touching resolved features. Damping rate goes as k^4, so a
    wave 8x longer should be damped ~4096x more slowly.
    """
    gr = CGrid(64, 64, 10e3, 10e3)
    K = recommended_hyper_coeff(gr)

    def checkerboard():
        return (np.indices((gr.ny, gr.nx)).sum(axis=0) % 2) * 2.0 - 1.0

    rates = {}
    for field, label in [(checkerboard(), "2dx (grid scale)"),
                         (np.sin(2 * np.pi * 8 * np.arange(gr.nx) / gr.nx)[None, :]
                          * np.ones((gr.ny, 1)), "8dx"),
                         (np.sin(2 * np.pi * 4 * np.arange(gr.nx) / gr.nx)[None, :]
                          * np.ones((gr.ny, 1)), "16dx")]:
        tend = hyperdiffusion(field, gr, K)
        rates[label] = np.abs(tend).max() / np.abs(field).max()

    ratio = rates["2dx (grid scale)"] / rates["16dx"]
    ok = ratio > 100
    report("hyperdiffusion is scale-selective", ok,
           f"damping rate 2dx {rates['2dx (grid scale)']:.2e}/s vs "
           f"16dx {rates['16dx']:.2e}/s -> {ratio:.0f}x more damping at grid scale")


# ---------------------------------------------------------------------------
def test_hyperdiffusion_damping_time():
    """The auto coefficient should damp the 2dx wave on roughly the requested time."""
    gr = CGrid(64, 64, 10e3, 10e3)
    target = 3 * 3600.0
    K = recommended_hyper_coeff(gr, damping_time=target)

    # True 2D checkerboard -- the mode the coefficient is defined against.
    a = (np.indices((gr.ny, gr.nx)).sum(axis=0) % 2) * 2.0 - 1.0
    rate = np.abs(hyperdiffusion(a, gr, K)).max() / np.abs(a).max()
    e_time = 1.0 / rate

    ok = 0.8 * target < e_time < 1.25 * target
    report("grid-scale damping time near target", ok,
           f"e-folding {e_time/3600:.2f} h (target {target/3600:.1f} h)")


# ---------------------------------------------------------------------------
def test_noise_is_removed_but_signal_kept():
    """
    Practical check: a smooth wave plus grid-scale noise. After diffusion the
    noise should be largely gone and the wave largely intact.
    """
    gr = CGrid(64, 64, 10e3, 10e3)
    K = recommended_hyper_coeff(gr)

    x = np.arange(gr.nx)
    signal = np.sin(2 * np.pi * 3 * x / gr.nx)[None, :] * np.ones((gr.ny, 1))
    noise = 0.3 * ((np.indices((gr.ny, gr.nx)).sum(axis=0) % 2) * 2.0 - 1.0)
    a = signal + noise

    dt = 60.0
    for _ in range(int(6 * 3600 / dt)):
        a = a + dt * hyperdiffusion(a, gr, K)

    # Project onto each component.
    sig_amp = np.abs((a * signal).mean() / (signal * signal).mean())
    noise_amp = np.abs((a * noise).mean() / (noise * noise).mean())

    ok = sig_amp > 0.95 and noise_amp < 0.2
    report("noise removed, resolved signal preserved (6 h)", ok,
           f"signal retained {sig_amp*100:.1f}%, noise retained {noise_amp*100:.1f}%")


# ---------------------------------------------------------------------------
def test_stochastic_field_statistics():
    """Perturbation field must have the requested amplitude and be smooth."""
    gr = CGrid(96, 96, 20e3, 20e3)
    p = StochasticPerturbation(gr, amplitude=0.3, length_scale=400e3, seed=1)

    std = p.r.std()
    # Smoothness: neighbour correlation should be high for a smooth field.
    corr = np.corrcoef(p.r[:, :-1].ravel(), p.r[:, 1:].ravel())[0, 1]

    ok = 0.15 < std < 0.45 and corr > 0.95
    report("stochastic field has correct amplitude and is smooth", ok,
           f"std = {std:.3f} (target 0.30), neighbour correlation = {corr:.4f}")


# ---------------------------------------------------------------------------
def test_stochastic_temporal_correlation():
    """
    AR(1) evolution: correlation after one step should be exp(-dt/tau), and
    the variance must stay stationary rather than decaying to zero.
    """
    gr = CGrid(64, 64, 20e3, 20e3)
    tau = 6 * 3600.0
    p = StochasticPerturbation(gr, amplitude=0.3, tau=tau, seed=2)

    dt = 3600.0
    r0 = p.r.copy()
    p.advance(dt)
    corr = np.corrcoef(r0.ravel(), p.r.ravel())[0, 1]
    expected = np.exp(-dt / tau)

    # Long run: variance must not drift.
    stds = []
    for _ in range(60):
        p.advance(dt)
        stds.append(p.r.std())
    drift = abs(stds[-1] - stds[0]) / stds[0]

    ok = abs(corr - expected) < 0.15 and drift < 0.5
    report("stochastic field decorrelates on tau, variance stationary", ok,
           f"1 h correlation {corr:.3f} (expect {expected:.3f}); "
           f"std {stds[0]:.3f} -> {stds[-1]:.3f} over 60 h")


# ---------------------------------------------------------------------------
def test_ensemble_members_diverge():
    """
    THE point of stochastic variance. Two runs from slightly different initial
    states must diverge -- that is atmospheric predictability, and a model
    that keeps them identical is not representing the real system.
    """
    gr = CGrid(40, 40, 50e3, 50e3, f0=1.0e-4, beta=1.6e-11)
    lev = PressureLevels(LEVELS)
    rng = np.random.default_rng(7)

    def member(perturb):
        m = Primitive3D(gr, lev)
        k_y = 2 * np.pi / gr.Ly
        for k in range(lev.nz):
            m.theta[k] = theta_from_T(260.0 + 8.0 * np.cos(k_y * gr.Yc), lev.p[k])
        phi = m.geopotential()
        for k in range(lev.nz):
            m.u[k] = -0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k])) / gr.f0
        if perturb:
            m.theta = perturb_initial_state(m.theta, rng, 0.1, gr)
        m.run(12 * 3600, dt=m.max_dt())
        return m

    a = member(False)
    b = member(True)

    spread = np.abs(a.u - b.u).max()
    ok = spread > 0.1 and np.isfinite(spread)
    report("perturbed ensemble members diverge over 12 h", ok,
           f"max|u_a - u_b| = {spread:.3f} m/s from a 0.1 K initial perturbation")


# ---------------------------------------------------------------------------
def test_model_actually_evolves():
    """
    A model that only preserves balanced states is not forecasting anything.
    Start from a baroclinically unstable jet with a small perturbation and
    confirm the perturbation GROWS -- weather developing from a smooth state,
    which is what baroclinic instability does in the real atmosphere.
    """
    gr = CGrid(48, 48, 60e3, 60e3, f0=1.0e-4, beta=1.6e-11)
    lev = PressureLevels(LEVELS)
    m = Primitive3D(gr, lev)

    # Strong meridional temperature contrast = strong vertical shear.
    k_y = 2 * np.pi / gr.Ly
    for k in range(lev.nz):
        m.theta[k] = theta_from_T(258.0 + 6.0 * np.cos(k_y * gr.Yc), lev.p[k])
    phi = m.geopotential()
    for k in range(lev.nz):
        m.u[k] = -0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k])) / gr.f0

    # Small wave perturbation to seed the instability.
    seed = 0.5 * np.sin(4 * np.pi * gr.Xc / gr.Lx) * np.sin(k_y * gr.Yc)
    for k in range(lev.nz):
        m.theta[k] += seed

    def eddy_energy(mm):
        up = mm.u - mm.u.mean(axis=2, keepdims=True)
        vp = mm.v - mm.v.mean(axis=2, keepdims=True)
        return float((up**2 + vp**2).sum())

    # The seed is in theta, so eddy WIND energy starts at exactly zero --
    # measuring growth from t=0 divides by nothing and reports a meaningless
    # ratio. Measure the growth rate between two later times instead, which
    # is how a baroclinic growth rate is actually diagnosed.
    m.run(24 * 3600, dt=m.max_dt() * 0.5)
    e1 = eddy_energy(m)
    m.run(24 * 3600, dt=m.max_dt() * 0.5)
    e2 = eddy_energy(m)

    growth = e2 / e1 if e1 > 0 else np.inf
    # Amplitude e-folding time from energy growth over 24 h.
    e_fold_h = 24.0 / (0.5 * np.log(growth)) if growth > 1 else np.inf
    vmax = np.abs(m.v).max()

    # Real baroclinic waves e-fold in roughly 1-3 days and reach tens of m/s.
    ok = (growth > 1.2 and np.isfinite(m.u).all() and 5.0 < vmax < 100.0
          and 6.0 < e_fold_h < 400.0)
    report("baroclinic perturbation grows (model produces weather)", ok,
           f"eddy energy x{growth:.2f} over 24 h (day 1 -> day 2), "
           f"amplitude e-folding {e_fold_h:.0f} h; max|v| {vmax:.1f} m/s")


if __name__ == "__main__":
    print("\nDissipation and stochastic variance\n" + "=" * 62)
    for fn in (test_hyperdiffusion_is_scale_selective,
               test_hyperdiffusion_damping_time,
               test_noise_is_removed_but_signal_kept,
               test_stochastic_field_statistics,
               test_stochastic_temporal_correlation,
               test_ensemble_members_diverge,
               test_model_actually_evolves):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

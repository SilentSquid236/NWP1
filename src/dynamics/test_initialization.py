"""
Validation for the initialization filter.

Run:  python test_initialization.py
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from initialization import (spectral_lowpass, grid_scale_energy,
                            filter_initial_state)
import probe_failure as P

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def test_smooth_field_survives():
    """A field with no grid-scale content must pass through nearly unchanged."""
    gr = CGrid(90, 88, 12e3, 12e3)
    a = np.sin(2 * np.pi * gr.Xc / gr.Lx) * np.cos(2 * np.pi * gr.Yc / gr.Ly)
    b = spectral_lowpass(a, gr)
    err = np.abs(b - a).max() / np.abs(a).max()
    report("a domain-scale wave passes the filter unchanged", err < 1e-3,
           f"max relative change {err:.2e}")


def test_grid_scale_removed():
    """A 2dx checkerboard must be removed essentially completely."""
    gr = CGrid(90, 88, 12e3, 12e3)
    j, i = np.arange(gr.ny)[:, None], np.arange(gr.nx)[None, :]
    a = (-1.0) ** (i + j) * np.ones((gr.ny, gr.nx))
    b = spectral_lowpass(a, gr)
    report("a 2dx checkerboard is removed", np.abs(b).max() < 1e-10,
           f"residual amplitude {np.abs(b).max():.2e} from 1.0")


def test_energy_diagnostic():
    """White noise is mostly grid-scale; a smooth wave is not."""
    gr = CGrid(90, 88, 12e3, 12e3)
    rng = np.random.default_rng(0)
    noise = rng.normal(0, 1, (gr.ny, gr.nx))
    wave = np.sin(4 * np.pi * gr.Xc / gr.Lx)
    en, ew = grid_scale_energy(noise, gr), grid_scale_energy(wave, gr)
    report("grid-scale energy separates noise from structure",
           en > 0.5 and ew < 0.01,
           f"white noise {en*100:.0f}% at <=4dx, smooth wave {ew*100:.3f}%")


def test_filter_preserves_mean_profile():
    """Filtering theta must not change the stratification."""
    gr = CGrid(48, 48, 25e3, 25e3)
    rng = np.random.default_rng(1)
    th = np.linspace(290, 400, 20)[:, None, None] + rng.normal(0, 2, (20, 48, 48))
    u = rng.normal(0, 1, (20, 48, 48))
    _, _, tf = filter_initial_state(u, u, th, gr)
    d = np.abs(tf.mean(axis=(1, 2)) - th.mean(axis=(1, 2))).max()
    report("filtering theta preserves the mean profile", d < 1e-9,
           f"max change in level-mean theta {d:.2e} K")


def test_noisy_run_survives_after_filtering():
    """
    THE POINT OF THE MODULE.

    1.2 m/s of white noise kills the run in 1 hour with every combination of
    mixing and drag. Filtering the same state, then rebalancing it, should
    complete 12 hours -- which would show the failure was unresolved initial
    variance rather than a defect in the dynamics.

    ORDER MATTERS AND WAS MEASURED. Filtering changes u, v and theta
    separately, so it puts divergence back into a balanced state; balancing
    afterwards removes it. Measured: no filter 1/12 h, filter only 11/12 h,
    filter then balance 12/12 h.
    """
    from subgrid import balance_initial_state

    m = P.build(0.0, 1.2, dT=1.5, clip=None)
    noise_before = float(np.sqrt(((m.u - spectral_lowpass(m.u, m.grid)) ** 2).mean()))

    m.u, m.v, m.theta = filter_initial_state(m.u, m.v, m.theta, m.grid)
    m.u, m.v, _ = balance_initial_state(m.u, m.v, m.grid, verbose=False)
    noise_after = float(np.sqrt(((m.u - spectral_lowpass(m.u, m.grid)) ** 2).mean()))

    dt = m.max_dt()
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    umax = float(np.abs(m.u).max()) if np.isfinite(m.u).all() else float("nan")
    report("12 h from a FILTERED, REBALANCED noisy state stays stable",
           done == 12 and umax < 60.0,
           f"survived {done}/12 h; max|u| {umax:.1f} m/s; sub-4dx wind rms "
           f"{noise_before:.3f} -> {noise_after:.3f} m/s")


if __name__ == "__main__":
    print("\nInitialization filter\n" + "=" * 66)
    for fn in (test_smooth_field_survives, test_grid_scale_removed,
               test_energy_diagnostic, test_filter_preserves_mean_profile,
               test_noisy_run_survives_after_filtering):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

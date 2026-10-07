"""
Validation for dry convective adjustment.

Run:  python test_convection.py
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels, P0, KAPPA, RD, G0
from convection import dry_convective_adjustment, unstable_fraction

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def column(theta_profile, ny=4, nx=4, p_s=101325.0):
    lev = SigmaLevels(len(theta_profile))
    pi = np.full((ny, nx), p_s - lev.p_top)
    th = np.repeat(np.asarray(theta_profile, float)[:, None, None], ny, 1)
    th = np.repeat(th, nx, 2)
    return lev, pi, th


def test_stable_column_untouched():
    """A stably stratified column must come back bit-identical."""
    lev, pi, th = column(np.linspace(400, 290, 20))   # decreasing with index
    u = np.zeros_like(th) + 10.0
    t2, u2, v2, info = dry_convective_adjustment(th, u, u, pi, lev)
    ok = np.array_equal(t2, th) and info["sweeps"] == 0
    report("a stable column is untouched", ok,
           f"max|dtheta| {np.abs(t2-th).max():.2e} K, {info['sweeps']} sweeps")


def test_inversion_removed():
    """An overturned layer must come back neutral."""
    prof = np.linspace(400, 290, 20)
    prof[10], prof[11] = prof[11], prof[10]      # invert one pair
    lev, pi, th = column(prof)
    u = np.zeros_like(th)
    t2, _, _, info = dry_convective_adjustment(th, u, u, pi, lev)
    ok = info["unstable_after"] == 0.0 and info["unstable_before"] > 0
    report("an overturned pair is mixed to neutral", ok,
           f"unstable interfaces {info['unstable_before']*100:.1f}% -> "
           f"{info['unstable_after']*100:.1f}% in {info['sweeps']} sweeps")


def test_enthalpy_conserved():
    """Mass-weighted theta must be conserved to round-off."""
    rng = np.random.default_rng(0)
    prof = np.linspace(400, 290, 20) + rng.normal(0, 15, 20)   # badly mixed
    lev, pi, th = column(prof)
    u = rng.normal(0, 10, th.shape)
    dm = lev.dsigma[:, None, None] * pi
    h0 = float((dm * th).sum())
    m0 = float((dm * u).sum())
    t2, u2, _, info = dry_convective_adjustment(th, u, u, pi, lev)
    h1 = float((dm * t2).sum())
    m1 = float((dm * u2).sum())
    rel_h = abs(h1 - h0) / abs(h0)
    rel_m = abs(m1 - m0) / max(abs(m0), 1e-12)
    ok = rel_h < 1e-12 and rel_m < 1e-10 and info["unstable_after"] == 0.0
    report("enthalpy and momentum conserved by the adjustment", ok,
           f"relative change: heat {rel_h:.2e}, momentum {rel_m:.2e}; "
           f"{info['sweeps']} sweeps to neutral")


def test_fully_inverted_column_converges():
    """A completely inverted column must reach a single well-mixed value."""
    lev, pi, th = column(np.linspace(290, 400, 20))   # increasing = unstable
    u = np.zeros_like(th)
    t2, _, _, info = dry_convective_adjustment(th, u, u, pi, lev,
                                               max_sweeps=200)
    spread = float(t2.max() - t2.min())
    ok = info["unstable_after"] == 0.0 and spread < 1e-6
    report("a fully inverted column mixes to uniform theta", ok,
           f"spread {t2.min():.3f}-{t2.max():.3f} K = {spread:.2e} after "
           f"{info['sweeps']} sweeps")


def test_momentum_mixing_can_be_disabled():
    """mix_momentum=False must leave the wind alone."""
    prof = np.linspace(400, 290, 20)
    prof[5], prof[6] = prof[6], prof[5]
    lev, pi, th = column(prof)
    u = np.arange(20, dtype=float)[:, None, None] * np.ones_like(th)
    _, u2, _, _ = dry_convective_adjustment(th, u, u, pi, lev,
                                            mix_momentum=False)
    report("momentum mixing can be switched off", np.array_equal(u2, u),
           f"max|du| {np.abs(u2-u).max():.2e} m/s")


# ---------------------------------------------------------------------------
# Pool-adjacent-violators (P-63)
# ---------------------------------------------------------------------------
from convection import dry_convective_adjustment_pav


def _noisy(ny=30, nx=30, amp=3.0, seed=1):
    rng = np.random.default_rng(seed)
    lev = SigmaLevels(20)
    pi = np.full((ny, nx), 101325.0 - lev.p_top) * (1 + 0.05 * rng.random((ny, nx)))
    th = np.linspace(330, 285, 20)[:, None, None] + amp * rng.standard_normal((20, ny, nx))
    u, v = rng.standard_normal((2, 20, ny, nx))
    return lev, pi, th, u, v


def _reference_pav(x, w, tol=1e-10):
    """Scalar PAV on one column, ground (last index) up. Slow and obvious."""
    blocks = []                                   # [sum w*x, sum w, n layers]
    for k in range(len(x) - 1, -1, -1):
        blocks.append([w[k] * x[k], w[k], 1])
        while len(blocks) > 1 and blocks[-1][0] / blocks[-1][1] < blocks[-2][0] / blocks[-2][1] - tol:
            s, m, n = blocks.pop()
            blocks[-1][0] += s; blocks[-1][1] += m; blocks[-1][2] += n
    out = []
    for s, m, n in blocks:
        out += [s / m] * n
    return np.array(out[::-1])


def test_pav_removes_instability_and_conserves():
    lev, pi, th, u, v = _noisy()
    t2, u2, v2, info = dry_convective_adjustment_pav(th, u, v, pi, lev)
    w = np.asarray(lev.dsigma)[:, None, None] * pi[None]
    rel = max(float(np.abs((w * a2).sum(0) - (w * a).sum(0)).max() / np.abs((w * a).sum(0)).max())
              for a, a2 in ((th, t2), (u, u2), (v, v2)))
    ok = info["unstable_before"] > 0.05 and info["unstable_after"] == 0 and rel < 1e-12
    report("PAV removes all instability in one pass and conserves theta, u, v", ok,
           f"unstable {info['unstable_before']:.3f} -> {info['unstable_after']}, "
           f"column drift {rel:.1e}")


def test_pav_matches_scalar_reference():
    lev, pi, th, u, v = _noisy(ny=5, nx=10, seed=2)
    t2, _, _, _ = dry_convective_adjustment_pav(th, u, v, pi, lev)
    w = np.asarray(lev.dsigma)
    err = max(float(np.abs(t2[:, j, i] - _reference_pav(th[:, j, i], w * pi[j, i])).max())
              for j in range(5) for i in range(10))
    report("vectorised PAV equals a scalar reference PAV, column by column", err < 1e-9,
           f"max difference {err:.1e} K over 50 columns")


def test_pav_is_least_change():
    """PAV is the least mixing; the converged sweep scheme mixes more."""
    lev, pi, th, u, v = _noisy(seed=3)
    tp = dry_convective_adjustment_pav(th, u, v, pi, lev)[0]
    ts, _, _, si = dry_convective_adjustment(th, u, v, pi, lev, max_sweeps=10000)
    w = np.asarray(lev.dsigma)[:, None, None] * pi[None]
    cp, cs = (w * (tp - th) ** 2).sum(0), (w * (ts - th) ** 2).sum(0)
    # Relative 1e-9: where both schemes form the same blocks the two sums agree
    # only to rounding (2e-11 relative was seen), not to 1e-12.
    ok = si["unstable_after"] == 0 and bool((cp <= cs * (1 + 1e-9)).all()) and cp.sum() < cs.sum()
    report("PAV changes each column no more than the converged sweep scheme", ok,
           f"mass-weighted squared change: PAV {cp.sum():.3e}, sweep {cs.sum():.3e} "
           f"({si['sweeps']} sweeps); max theta difference {np.abs(tp - ts).max():.2f} K")


def test_pav_clears_neutral_chain():
    """The P-63 state: a neutral block nudged by 1e-4 K defeats the 20-sweep cap."""
    lev, pi, th, u, v = _noisy(ny=40, nx=40, amp=0.0, seed=4)
    th[-10:] = th[-10:].mean(0)
    th += 1e-4 * np.random.default_rng(5).standard_normal(th.shape)
    s20 = dry_convective_adjustment(th, u, v, pi, lev)[3]
    p = dry_convective_adjustment_pav(th, u, v, pi, lev)[3]
    ok = s20["sweeps"] == 20 and s20["unstable_after"] > 0 and p["unstable_after"] == 0
    report("PAV clears the near-neutral chain the 20-sweep cap leaves behind", ok,
           f"before {s20['unstable_before']:.3f}; sweep cap 20 leaves {s20['unstable_after']:.3f}; "
           f"PAV leaves {p['unstable_after']}")


def test_pav_stable_column_untouched():
    lev, pi, th = column(np.linspace(400, 290, 20))
    u = np.zeros_like(th) + 10.0
    t2, u2, v2, info = dry_convective_adjustment_pav(th, u, u, pi, lev)
    ok = np.array_equal(t2, th) and np.array_equal(u2, u) and info["sweeps"] == 0
    report("PAV leaves a stable column bit-identical", ok,
           f"max|dtheta| {np.abs(t2 - th).max():.1e} K, {info['columns']} columns touched")


def test_pav_single_inversion_matches_sweep():
    prof = np.linspace(400, 290, 20)
    prof[10], prof[11] = prof[11], prof[10]
    lev, pi, th = column(prof)
    u = np.arange(20, dtype=float)[:, None, None] * np.ones_like(th)
    a = dry_convective_adjustment(th, u, u, pi, lev)
    b = dry_convective_adjustment_pav(th, u, u, pi, lev)
    err = max(float(np.abs(x - y).max()) for x, y in zip(a[:3], b[:3]))
    report("one inverted pair: PAV and sweep give the same answer", err < 1e-12,
           f"max difference {err:.1e}")


def test_pav_torch_matches_numpy():
    try:
        import torch
    except ImportError:
        print("  [SKIP] PAV on torch tensors\n        torch is not installed here")
        return
    lev, pi, th, u, v = _noisy(seed=6)
    a = dry_convective_adjustment_pav(th, u, v, pi, lev)
    b = dry_convective_adjustment_pav(*(torch.from_numpy(x) for x in (th, u, v, pi)), lev)
    ok = isinstance(b[0], torch.Tensor) and b[0].dtype == torch.float64 and all(
        np.array_equal(x, y.numpy()) for x, y in zip(a[:3], b[:3]))
    report("PAV on torch tensors equals NumPy and returns float64 tensors", ok,
           f"types {type(b[0]).__name__}, {b[0].dtype}")


if __name__ == "__main__":
    print("\nDry convective adjustment\n" + "=" * 66)
    for fn in (test_stable_column_untouched, test_inversion_removed,
               test_enthalpy_conserved, test_fully_inverted_column_converges,
               test_momentum_mixing_can_be_disabled,
               test_pav_removes_instability_and_conserves,
               test_pav_matches_scalar_reference,
               test_pav_is_least_change,
               test_pav_clears_neutral_chain,
               test_pav_stable_column_untouched,
               test_pav_single_inversion_matches_sweep,
               test_pav_torch_matches_numpy):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

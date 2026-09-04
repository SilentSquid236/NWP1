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


if __name__ == "__main__":
    print("\nDry convective adjustment\n" + "=" * 66)
    for fn in (test_stable_column_untouched, test_inversion_removed,
               test_enthalpy_conserved, test_fully_inverted_column_converges,
               test_momentum_mixing_can_be_disabled):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

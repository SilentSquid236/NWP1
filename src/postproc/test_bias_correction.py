"""
Tests for adaptive bias correction.

The important ones are the negative tests: a corrector that only ever helps in
its own test suite is not being tested honestly. These check that it declines
to act without evidence, that it cannot run away, and that it does NOT improve
things when there is no systematic bias to remove.

Run:  python test_bias_correction.py
"""

import numpy as np

from bias_correction import KalmanBiasCorrector, verify, skill_score

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# ---------------------------------------------------------------------------
def test_learns_constant_bias():
    """A station running persistently warm should be corrected toward truth."""
    rng = np.random.default_rng(0)
    kbc = KalmanBiasCorrector(obs_var=1.0, min_samples=5)
    true_bias = 2.5

    raw, corr, obs = [], [], []
    for i in range(60):
        truth = 273.0 + 10 * np.sin(i / 8) + rng.normal(0, 1.0)
        fcst = truth + true_bias
        c = kbc.apply("KALB|f06|T", fcst)
        if i >= 10:
            raw.append(fcst); corr.append(c); obs.append(truth)
        kbc.update("KALB|f06|T", fcst, truth)

    est = kbc.bias["KALB|f06|T"]
    ss = skill_score(corr, raw, obs)
    ok = abs(est - true_bias) < 0.5 and ss > 0.3
    report("learns a constant warm bias", ok,
           f"estimated {est:.2f} K (true {true_bias}), "
           f"RMSE improved {ss*100:.0f}%")


# ---------------------------------------------------------------------------
def test_no_correction_without_evidence():
    """Must refuse to act before min_samples -- a single outlier is not a bias."""
    kbc = KalmanBiasCorrector(min_samples=10)
    for i in range(5):
        kbc.update("X", 10.0, 5.0)          # huge apparent bias
    early = kbc.correction("X")
    for i in range(10):
        kbc.update("X", 10.0, 5.0)
    later = kbc.correction("X")

    ok = early == 0.0 and later > 1.0
    report("no correction applied before enough samples", ok,
           f"after 5 updates: {early:.2f} (must be 0), "
           f"after 15: {later:.2f}")


# ---------------------------------------------------------------------------
def test_correction_is_capped():
    """A pathological error stream must not produce an unbounded correction."""
    kbc = KalmanBiasCorrector(max_correction=3.0, min_samples=2)
    for _ in range(100):
        kbc.update("Y", 100.0, 0.0)         # absurd 100-unit error
    c = kbc.correction("Y")

    ok = abs(c) <= 3.0
    report("correction respects the hard cap", ok,
           f"raw bias estimate {kbc.bias['Y']:.1f}, applied {c:.1f} (cap 3.0)")


# ---------------------------------------------------------------------------
def test_does_not_help_when_no_bias():
    """
    NEGATIVE TEST. With unbiased random error there is nothing systematic to
    remove, so the corrector should be roughly neutral. If it showed large
    gains here it would be fitting noise -- which would mean it also destroys
    signal on real data.
    """
    rng = np.random.default_rng(1)
    kbc = KalmanBiasCorrector(obs_var=1.0, min_samples=5)

    raw, corr, obs = [], [], []
    for i in range(200):
        truth = 280.0 + rng.normal(0, 3.0)
        fcst = truth + rng.normal(0, 1.5)      # unbiased noise
        c = kbc.apply("Z", fcst)
        if i >= 20:
            raw.append(fcst); corr.append(c); obs.append(truth)
        kbc.update("Z", fcst, truth)

    ss = skill_score(corr, raw, obs)
    ok = abs(ss) < 0.10
    report("neutral when there is no systematic bias (does not fit noise)", ok,
           f"skill score {ss*100:+.1f}% (want near zero, |ss| < 10%)")


# ---------------------------------------------------------------------------
def test_adapts_to_regime_change():
    """
    A bias that shifts -- seasonal change, instrument swap, model update --
    must be tracked. This is what the Kalman filter buys over a fixed mean.
    """
    rng = np.random.default_rng(2)
    kbc = KalmanBiasCorrector(process_var=0.05, obs_var=1.0, min_samples=5)

    for i in range(80):
        truth = 275.0 + rng.normal(0, 0.5)
        kbc.update("W", truth + 3.0, truth)      # warm bias
    before = kbc.bias["W"]

    for i in range(80):
        truth = 275.0 + rng.normal(0, 0.5)
        kbc.update("W", truth - 2.0, truth)      # flips cold
    after = kbc.bias["W"]

    ok = before > 2.0 and after < -1.0
    report("tracks a bias that changes sign", ok,
           f"bias estimate {before:+.2f} -> {after:+.2f} K after the regime flip")


# ---------------------------------------------------------------------------
def test_keys_are_independent():
    """Stations and lead times must not contaminate each other."""
    kbc = KalmanBiasCorrector(min_samples=3)
    for _ in range(20):
        kbc.update("A|f06|T", 5.0, 0.0)      # +5
        kbc.update("B|f06|T", -3.0, 0.0)     # -3
        kbc.update("A|f24|T", 0.0, 0.0)      # unbiased

    a6, b6, a24 = (kbc.correction("A|f06|T"), kbc.correction("B|f06|T"),
                   kbc.correction("A|f24|T"))
    ok = a6 > 3.0 and b6 < -2.0 and abs(a24) < 0.5
    report("per-station, per-lead-time keys stay independent", ok,
           f"A/f06 {a6:+.2f}, B/f06 {b6:+.2f}, A/f24 {a24:+.2f}")


# ---------------------------------------------------------------------------
def test_gain_decreases_with_evidence():
    """
    The Kalman gain should start high (learn fast from nothing) and settle as
    evidence accumulates -- that is the adaptivity, and it should be visible.
    """
    kbc = KalmanBiasCorrector(process_var=0.001, obs_var=1.0)
    gains = []
    for i in range(30):
        gains.append(kbc.gain("G"))
        kbc.update("G", 1.0, 0.0)

    ok = gains[0] > gains[5] > gains[-1] and gains[-1] < 0.2
    report("Kalman gain falls as evidence accumulates", ok,
           f"gain {gains[0]:.3f} -> {gains[5]:.3f} -> {gains[-1]:.3f}")


if __name__ == "__main__":
    print("\nAdaptive bias correction\n" + "=" * 62)
    for fn in (test_learns_constant_bias,
               test_no_correction_without_evidence,
               test_correction_is_capped,
               test_does_not_help_when_no_bias,
               test_adapts_to_regime_change,
               test_keys_are_independent,
               test_gain_decreases_with_evidence):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

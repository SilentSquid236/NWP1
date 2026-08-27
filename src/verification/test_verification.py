"""
Tests for the observation operator, QC, and verification scoring.

Run:  python test_verification.py
"""

import sys, tempfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "dynamics"))

from grid import CGrid
from vertical import PressureLevels
from observations import (Observation, run_qc, range_check, buddy_check,
                          default_error_std)
from obs_operator import GridInterpolator, elevation_correct_temperature
from scoring import (scores, scores_by, match_forecast_to_obs, ForecastArchive,
                     skill_vs_reference)

DOMAIN = {"name": "test", "lat_min": 37.0, "lat_max": 47.5,
          "lon_min": -82.0, "lon_max": -66.0}
LEVELS = [1000, 975, 950, 925, 900, 875, 850, 825, 800, 750,
          700, 650, 600, 550, 500, 450, 400, 300, 250, 200]

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def make_interp(n=64):
    dx = (DOMAIN["lon_max"] - DOMAIN["lon_min"]) * 111320.0 * \
         np.cos(np.radians(42.25)) / n
    dy = (DOMAIN["lat_max"] - DOMAIN["lat_min"]) * 111132.0 / n
    gr = CGrid(n, n, dx, dy)
    return GridInterpolator(gr, DOMAIN, PressureLevels(LEVELS)), gr


# ---------------------------------------------------------------------------
def test_bilinear_exact_for_linear_field():
    """
    Bilinear interpolation reproduces a linear field EXACTLY. This is the
    strongest available check on the horizontal operator: any indexing or
    coordinate-mapping error breaks exactness immediately.
    """
    interp, gr = make_interp()
    # Field varying linearly in x and y.
    field = 2.0 * gr.Xc / gr.Lx + 3.0 * gr.Yc / gr.Ly

    errs = []
    rng = np.random.default_rng(0)
    for _ in range(50):
        lat = rng.uniform(38.0, 46.5)
        lon = rng.uniform(-80.5, -67.5)
        got = interp.horizontal(field, lat, lon)
        x, y = interp.lonlat_to_xy(lat, lon)
        want = 2.0 * x / gr.Lx + 3.0 * y / gr.Ly
        if got is not None:
            errs.append(abs(got - want))

    worst = max(errs)
    ok = worst < 1e-10 and len(errs) > 40
    report("bilinear interpolation exact for a linear field", ok,
           f"max error {worst:.2e} over {len(errs)} random points")


# ---------------------------------------------------------------------------
def test_out_of_domain_returns_none():
    """
    Points outside the domain must return None, not an extrapolated number.
    A fabricated match scored as if real is worse than a missing match.
    """
    interp, gr = make_interp()
    field = np.ones((gr.ny, gr.nx))

    outside = [(30.0, -74.0), (55.0, -74.0), (42.0, -95.0), (42.0, -50.0)]
    got = [interp.horizontal(field, la, lo) for la, lo in outside]
    inside = interp.horizontal(field, 42.0, -74.0)

    ok = all(g is None for g in got) and inside is not None
    report("outside-domain points return None, not extrapolation", ok,
           f"4 exterior points -> {got}, interior -> {inside}")


# ---------------------------------------------------------------------------
def test_vertical_interp_linear_in_logp():
    """
    Vertical interpolation must be linear in log(p). A profile constructed
    linear in log(p) is reproduced exactly; the same profile interpolated
    linearly in p would show visible error in the mid-troposphere.
    """
    interp, _ = make_interp()
    lev = interp.levels

    col = 3.0 + 2.0 * np.log(lev.p)         # exactly linear in log p
    errs = []
    for p_test in [95000.0, 82500.0, 62500.0, 47500.0, 22500.0]:
        got = interp.vertical(col, p_test)
        want = 3.0 + 2.0 * np.log(p_test)
        errs.append(abs(got - want))

    # How wrong would linear-in-p be, for comparison?
    p_mid = 47500.0
    k = np.searchsorted(-lev.p, -p_mid)
    w = (p_mid - lev.p[k]) / (lev.p[k - 1] - lev.p[k])
    lin_p = col[k] + w * (col[k - 1] - col[k])
    lin_p_err = abs(lin_p - (3.0 + 2.0 * np.log(p_mid)))

    ok = max(errs) < 1e-10
    report("vertical interpolation is exact in log(p)", ok,
           f"max error {max(errs):.2e}; linear-in-p would err "
           f"{lin_p_err:.2e} at 475 hPa")


# ---------------------------------------------------------------------------
def test_qc_rejects_bad_values():
    """Range and gross-error gates must fire, with reasons recorded."""
    t = datetime(2026, 8, 1, 12)
    obs = [
        Observation(t, 42.0, -74.0, "TMP", 288.0, "metar", "GOOD", error_std=1.5),
        Observation(t, 42.1, -74.1, "TMP", 999.0, "metar", "BADRANGE", error_std=1.5),
        Observation(t, 42.2, -74.2, "TMP", 250.0, "metar", "BADGROSS", error_std=1.5),
        Observation(t, 42.3, -74.3, "TMP", float("nan"), "metar", "NAN", error_std=1.5),
    ]
    bg = {"GOOD": 288.5, "BADGROSS": 288.0, "BADRANGE": 288.0, "NAN": 288.0}
    kept, rejected, summary = run_qc(obs, backgrounds=bg)

    names = {o.station for o in kept}
    reasons = {o.station: o.qc_reason for o in rejected}
    ok = names == {"GOOD"} and summary["range"] == 2 and summary["gross"] == 1
    report("QC rejects out-of-range, non-finite, and gross errors", ok,
           f"kept {sorted(names)}; summary {summary}; "
           f"example reason: {list(reasons.values())[0][:45]}")


# ---------------------------------------------------------------------------
def test_buddy_check_catches_plausible_liar():
    """
    A station reporting a physically plausible but locally wrong value passes
    range and gross checks. Only the buddy check catches it.
    """
    t = datetime(2026, 8, 1, 12)
    good = [Observation(t, 42.0 + 0.1 * i, -74.0, "TMP", 288.0 + 0.2 * i,
                        "metar", f"S{i}", error_std=1.5) for i in range(5)]
    liar = Observation(t, 42.2, -74.05, "TMP", 300.0, "metar", "LIAR",
                       error_std=1.5)

    ok_range, _ = range_check(liar)
    passed, why = buddy_check(liar, good)

    ok = ok_range and not passed
    report("buddy check catches a plausible but wrong station", ok,
           f"range check passed ({ok_range}), buddy rejected: {why[:60]}")


# ---------------------------------------------------------------------------
def test_scores_match_hand_calculation():
    """Verify the arithmetic against values computed by hand."""
    pairs = [(10.0, 8.0), (12.0, 11.0), (9.0, 10.0), (15.0, 12.0)]
    # errors: +2, +1, -1, +3 -> bias 1.25, rmse sqrt((4+1+1+9)/4)=sqrt(3.75)
    s = scores(pairs)
    ok = (abs(s["bias"] - 1.25) < 1e-12 and
          abs(s["rmse"] - np.sqrt(3.75)) < 1e-12 and
          abs(s["mae"] - 1.75) < 1e-12 and s["n"] == 4)
    report("scores match hand calculation", ok,
           f"bias {s['bias']:.4f} (1.25), rmse {s['rmse']:.4f} "
           f"({np.sqrt(3.75):.4f}), mae {s['mae']:.4f} (1.75)")


# ---------------------------------------------------------------------------
def test_matching_skips_and_counts():
    """
    Unmatched observations must be counted by reason, not silently dropped --
    otherwise a shrinking sample looks like a stable one.
    """
    interp, gr = make_interp()
    lev = interp.levels
    field = np.ones((lev.nz, gr.ny, gr.nx)) * 285.0

    t = datetime(2026, 8, 1, 12)
    obs = [
        Observation(t, 42.0, -74.0, "TMP", 284.0, "metar", "IN"),
        Observation(t, 20.0, -74.0, "TMP", 284.0, "metar", "OUTSIDE"),
        Observation(t + timedelta(hours=5), 42.0, -74.0, "TMP", 284.0, "metar", "LATE"),
        Observation(t, 42.0, -74.0, "RH", 50.0, "metar", "WRONGVAR"),
    ]
    for o in obs:
        o.qc_flag = "good"

    matches, skipped = match_forecast_to_obs(field, interp, obs, t, "TMP")

    ok = (len(matches) == 1 and skipped.get("outside domain or levels") == 1
          and skipped.get("outside time window") == 1
          and skipped.get("wrong variable") == 1)
    report("unmatched observations are counted by reason", ok,
           f"{len(matches)} matched; skipped {skipped}")


# ---------------------------------------------------------------------------
def test_elevation_correction_sign():
    """
    A station BELOW its grid cell should be warmer than the model value.
    A sign error here produces a systematic bias that looks like real forecast
    error, and post-processing would then learn to 'fix' a bug.
    """
    model_T, model_elev = 280.0, 500.0
    valley = elevation_correct_temperature(model_T, model_elev, 100.0)
    peak = elevation_correct_temperature(model_T, model_elev, 900.0)

    ok = valley > model_T and peak < model_T
    report("elevation correction has the right sign", ok,
           f"model {model_T:.1f} K at 500 m -> valley(100 m) {valley:.2f} K, "
           f"peak(900 m) {peak:.2f} K")


# ---------------------------------------------------------------------------
def test_archive_roundtrip_and_recovery():
    """
    The archive must round-trip, and must survive a truncated final line --
    which is exactly what a crash mid-write leaves behind.
    """
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "archive.jsonl"
        arc = ForecastArchive(path)
        recs = [{"time": "2026-08-01T12:00:00", "station": f"S{i}",
                 "variable": "TMP", "forecast": 280.0 + i,
                 "observation": 279.0 + i, "error": 1.0, "lead_hours": 6}
                for i in range(5)]
        arc.append(recs)

        # Simulate a crash mid-write.
        with open(path, "a") as f:
            f.write('{"time": "2026-08-01T13:00:00", "stat')

        loaded = arc.load()
        tmp_only = arc.load(variable="TMP")
        s1 = arc.load(station="S1")

        ok = len(loaded) == 5 and len(tmp_only) == 5 and len(s1) == 1
        report("archive round-trips and survives a truncated line", ok,
               f"{len(loaded)} records recovered after a partial write; "
               f"filter by station -> {len(s1)}")


if __name__ == "__main__":
    print("\nVerification harness\n" + "=" * 62)
    for fn in (test_bilinear_exact_for_linear_field,
               test_out_of_domain_returns_none,
               test_vertical_interp_linear_in_logp,
               test_qc_rejects_bad_values,
               test_buddy_check_catches_plausible_liar,
               test_scores_match_hand_calculation,
               test_matching_skips_and_counts,
               test_elevation_correction_sign,
               test_archive_roundtrip_and_recovery):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

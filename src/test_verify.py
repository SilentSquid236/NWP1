"""
Tests for the verification archiver.

Everything here runs offline: observations are a saved IEM-format payload and
the forecast is fabricated, so the whole archive path is exercised without a
network. The one thing this CANNOT test is the fetch itself (P-06).

The tests are written around the design decision that matters: raw
observations are written before anything can throw, because they are the only
irreplaceable part of the archive.

Run:  python src/test_verify.py
"""
import gzip
import json
import sys
import tempfile
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np
np.seterr(all="ignore")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE / "dynamics"))
sys.path.insert(0, str(HERE / "verification"))

import config
from sigma import SigmaLevels, P0, KAPPA, RD, G0
import verify as V
from scoring import ForecastArchive

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


RUN = datetime(2026, 9, 4, 0)
NY, NX = 40, 44

# Two stations inside the domain, at different elevations.
STATIONS = [("KBOS", 42.36, -71.01, 6.0), ("KBTV", 44.47, -73.15, 100.0)]


def fake_asos_csv(hours=3):
    """An IEM-format ASOS payload: the shape fetch_asos returns."""
    lines = ["station,valid,lon,lat,elevation,tmpf,drct,sknt"]
    for h in range(hours + 1):
        t = RUN + timedelta(hours=h)
        for st, lat, lon, elev in STATIONS:
            lines.append(f"{st},{t:%Y-%m-%d %H:%M},{lon},{lat},{elev},"
                         f"{59.0 + h:.1f},270,10")
    return "\n".join(lines) + "\n"


def fake_forecast(path, hours=3, terrain_height=200.0):
    """A sigma forecast file with the fields verify.py requires."""
    lev = SigmaLevels(20)
    nt = hours
    terrain = np.full((NY, NX), terrain_height)
    pi = np.full((NY, NX), 101325.0 * np.exp(
        -G0 * terrain_height / (RD * 280.0)) - lev.p_top)
    p = lev.pressure(pi)
    T0, L = 288.15, 0.0065
    T = T0 * (p / 101325.0) ** (RD * L / G0)
    theta = T / (p / P0) ** KAPPA

    np.savez_compressed(
        path,
        times_s=np.arange(1, nt + 1) * 3600.0,
        u=np.stack([np.full((20, NY, NX), 8.0)] * nt),
        v=np.stack([np.zeros((20, NY, NX))] * nt),
        theta=np.stack([theta] * nt),
        pi=np.stack([pi] * nt),
        sigma=lev.sigma, p_top=lev.p_top, terrain=terrain,
    )
    return path


def prepared(d, hours=3):
    """A run directory with the observations already cached."""
    d = Path(d)
    fc = fake_forecast(d / "forecast.npz", hours)
    paths = V.archive_paths(d / "archive", RUN)
    first = RUN + timedelta(hours=1)
    last = RUN + timedelta(hours=hours)
    V.store_raw(paths["raw"],
                f"asos_{first:%Y%m%d%H}_{last:%Y%m%d%H}", fake_asos_csv(hours))
    return fc, d / "archive", paths


# ---------------------------------------------------------------------------
def test_raw_observations_are_stored_verbatim():
    """
    The irreplaceable part. Stored bytes must come back identical -- not
    reparsed, not normalised.
    """
    with tempfile.TemporaryDirectory() as d:
        text = fake_asos_csv()
        p = V.store_raw(Path(d) / "obs", "sample", text)
        back = V.load_raw(Path(d) / "obs", "sample")
        ok = back == text and p.suffixes[-2:] == [".csv", ".gz"]
        report("raw observations round-trip verbatim", ok,
               f"{len(text)} chars in, {len(back or '')} back, "
               f"compressed to {p.stat().st_size} bytes")


# ---------------------------------------------------------------------------
def test_archive_written_and_recomputable():
    """
    A full offline pass: raw observations, the forecast copied beside them,
    matches derived. The forecast copy is what makes the archive
    self-contained -- matches can be recomputed from it later.
    """
    with tempfile.TemporaryDirectory() as d:
        fc, root, paths = prepared(d)
        matches, out = V.verify(fc, root, run_time=RUN, report_only=True,
                                verbose=False)
        meta = json.loads(paths["meta"].read_text())
        ok = (len(matches) > 0 and paths["forecast"].exists()
              and paths["matches"].exists() and meta["n_matches"] == len(matches))
        report("archive holds raw obs, the forecast, and the matches", ok,
               f"{len(matches)} pairs, forecast copied "
               f"({paths['forecast'].stat().st_size/1e3:.0f} kB), "
               f"{meta['n_observations']} observations parsed")


# ---------------------------------------------------------------------------
def test_rerun_does_not_duplicate():
    """
    A daily job WILL be run twice -- a retry, a cron overlap, a manual check.
    The archive must not double-count, or every score is silently weighted by
    how often someone re-ran the script.
    """
    with tempfile.TemporaryDirectory() as d:
        fc, root, paths = prepared(d)
        first, _ = V.verify(fc, root, run_time=RUN, report_only=True,
                            verbose=False)
        second, _ = V.verify(fc, root, run_time=RUN, report_only=True,
                             verbose=False)
        total = len(ForecastArchive(paths["matches"]).load())
        ok = len(first) > 0 and len(second) == 0 and total == len(first)
        report("re-running the same day adds nothing", ok,
               f"first run {len(first)} pairs, second run {len(second)}, "
               f"archive holds {total}")


# ---------------------------------------------------------------------------
def test_lead_hours_recorded():
    """
    Skill decays with lead time; an archive that does not record it cannot
    show that, and lead time is not recoverable from the pair alone.
    """
    with tempfile.TemporaryDirectory() as d:
        fc, root, paths = prepared(d)
        matches, _ = V.verify(fc, root, run_time=RUN, report_only=True,
                              verbose=False)
        leads = sorted({m["lead_hours"] for m in matches})
        ok = leads == [1.0, 2.0, 3.0]
        report("lead time is recorded on every pair", ok,
               f"leads present: {leads}")


# ---------------------------------------------------------------------------
def test_elevation_correction_is_recorded():
    """
    The operator's own correction has to be visible in the archive, or a bias
    it causes cannot be told apart from a bias in the model.
    """
    with tempfile.TemporaryDirectory() as d:
        fc, root, paths = prepared(d, hours=1)
        matches, _ = V.verify(fc, root, run_time=RUN, report_only=True,
                              verbose=False)
        temps = [m for m in matches if m["variable"] == "TMP"]
        ok = temps and all("elev_correction_m" in m for m in temps)
        c = [m["elev_correction_m"] for m in temps]
        report("elevation correction is recorded with each temperature pair",
               bool(ok),
               f"{len(temps)} temperature pairs, corrections "
               f"{min(c):+.0f} to {max(c):+.0f} m" if temps else "none")


# ---------------------------------------------------------------------------
def test_report_only_refuses_without_cache():
    """
    --report-only must never quietly reach the network, and must say so when
    it has nothing cached rather than producing an empty archive.
    """
    with tempfile.TemporaryDirectory() as d:
        fc = fake_forecast(Path(d) / "forecast.npz")
        try:
            V.verify(fc, Path(d) / "archive", run_time=RUN, report_only=True,
                     verbose=False)
            ok, why = False, "ran with no cached observations"
        except SystemExit as e:
            ok = "cached" in str(e)
            why = "refused and said to run once without --report-only"
    report("--report-only refuses rather than fetching", ok, why)


# ---------------------------------------------------------------------------
def test_old_forecast_format_is_refused():
    """
    A forecast written before the sigma port has no pi/sigma/terrain. Scoring
    it against station elevations is not possible, and guessing is worse than
    stopping.
    """
    with tempfile.TemporaryDirectory() as d:
        old = Path(d) / "old.npz"
        np.savez_compressed(old, times_s=np.array([3600.0]),
                            u=np.zeros((1, 20, NY, NX)),
                            v=np.zeros((1, 20, NY, NX)),
                            theta=np.zeros((1, 20, NY, NX)))
        try:
            V.load_forecast(old)
            ok, why = False, "accepted a pressure-coordinate forecast"
        except SystemExit as e:
            ok = "pi" in str(e) or "sigma" in str(e)
            why = "refused and named the missing fields"
    report("a pre-sigma forecast file is refused", ok, why)


if __name__ == "__main__":
    print("\nVerification archiver\n" + "=" * 66)
    for fn in (test_raw_observations_are_stored_verbatim,
               test_archive_written_and_recomputable,
               test_rerun_does_not_duplicate,
               test_lead_hours_recorded,
               test_elevation_correction_is_recorded,
               test_report_only_refuses_without_cache,
               test_old_forecast_format_is_refused):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

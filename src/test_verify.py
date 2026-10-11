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


def test_auto_surface_operator():
    """
    Default since 2026-10-03: a land-surface forecast is scored with the
    similarity operator, persistence and forecasts without a ground
    temperature with the standard one (persistence with a held ground
    temperature is not a fair reference; test AK).
    """
    tg = np.zeros((1, NY, NX))
    cases = [({"tg": tg, "stopped": np.array("completed")}, "similarity"),
             ({"tg": tg, "stopped": np.array("persistence")}, "standard"),
             ({"stopped": np.array("completed")}, "standard"),
             ({"tg": tg}, "similarity")]
    got = [V.resolve_surface_operator(fc) for fc, _ in cases]
    ok = got == [want for _, want in cases]
    report("auto operator: similarity for land-surface runs, standard for persistence",
           ok, f"got {got}")


# ---------------------------------------------------------------------------
def test_dewpoint_scored_for_moist_forecasts():
    """
    CAM stage S7: ASOS dwpf becomes DPT observations; a forecast with qv is
    scored against them (lowest-level qv at the station pressure), a dry one
    skips them with a named reason. A uniform qv chosen to give Td = 10 C at
    the model surface must score within 0.5 K of 10 C at a station near the
    model's elevation.
    """
    import moisture as mo
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fc = fake_forecast(d / "forecast.npz", 3, terrain_height=6.0)
        z = dict(np.load(fc))
        ps = z["p_top"] + z["pi"][0]
        qv = float(mo.qsat(283.15, float(ps.mean())))          # Td = 10 C at p_s
        z["qv"] = np.full(z["theta"].shape, qv)
        np.savez_compressed(d / "moist.npz", **z)
        lines = ["station,valid,lon,lat,elevation,tmpf,dwpf,drct,sknt"]
        for h in range(4):
            t = RUN + timedelta(hours=h)
            lines.append(f"KBOS,{t:%Y-%m-%d %H:%M},-71.01,42.36,6.0,59.0,50.0,270,10")
        first, last = RUN + timedelta(hours=1), RUN + timedelta(hours=3)
        for name in ("dry", "moist"):
            paths = V.archive_paths(d / name, RUN)
            V.store_raw(paths["raw"], f"asos_{first:%Y%m%d%H}_{last:%Y%m%d%H}", "\n".join(lines) + "\n")
        m_dry, _ = V.verify(fc, d / "dry", run_time=RUN, report_only=True, verbose=False)
        m_wet, _ = V.verify(d / "moist.npz", d / "moist", run_time=RUN, report_only=True, verbose=False)
        dpt = [m for m in m_wet if m["variable"] == "DPT"]
        td_obs = (50.0 - 32.0) * 5 / 9 + 273.15
        ok = (len(dpt) == 3 and not any(m["variable"] == "DPT" for m in m_dry)
              and all(abs(m["forecast"] - 283.15) < 0.5 for m in dpt)
              and all(abs(m["observation"] - td_obs) < 1e-6 for m in dpt))
        report("dewpoint: DPT pairs for a moist forecast (Td within 0.5 K), none for a dry one", ok,
               f"{len(dpt)} DPT pairs, forecast Td {[round(m['forecast'] - 273.15, 2) for m in dpt]} C, "
               f"obs {td_obs - 273.15:.2f} C; dry run DPT pairs "
               f"{sum(m['variable'] == 'DPT' for m in m_dry)}")


# ---------------------------------------------------------------------------
def test_sea_level_pressure_scored():
    """
    P-76: ASOS mslp becomes PMSL observations, scored against the forecast's
    sea-level pressure reduced as the maps reduce it (maps/derive.py). With the
    station near the model's terrain height the forecast MSLP must equal the
    map MSLP field at that point to 0.01 hPa.
    """
    sys.path.insert(0, str(HERE / "maps"))
    from derive import column_state, mslp
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fc = fake_forecast(d / "forecast.npz", 3, terrain_height=6.0)
        z = dict(np.load(fc))
        p, T, Z, ps = column_state(z["theta"][0], z["pi"][0], z["sigma"], float(z["p_top"]), z["terrain"])
        m = mslp(ps, T, Z, z["terrain"])
        lines = ["station,valid,lon,lat,elevation,tmpf,mslp,drct,sknt"]
        for h in range(4):
            t = RUN + timedelta(hours=h)
            lines.append(f"KBOS,{t:%Y-%m-%d %H:%M},-71.01,42.36,6.0,59.0,1012.3,270,10")
        first, last = RUN + timedelta(hours=1), RUN + timedelta(hours=3)
        paths = V.archive_paths(d / "a", RUN)
        V.store_raw(paths["raw"], f"asos_{first:%Y%m%d%H}_{last:%Y%m%d%H}", "\n".join(lines) + "\n")
        mt, _ = V.verify(fc, d / "a", run_time=RUN, report_only=True, verbose=False)
        pm = [x for x in mt if x["variable"] == "PMSL"]
        ok = (len(pm) == 3 and all(abs(x["observation"] - 101230.0) < 1e-6 for x in pm)
              and all(abs(x["forecast"] - float(m.mean())) < 1.0 for x in pm))
        report("sea-level pressure: PMSL pairs scored against the map MSLP (within 0.01 hPa)", ok,
               f"{len(pm)} PMSL pairs, forecast {[round(x['forecast'] / 100, 2) for x in pm]} hPa, "
               f"map MSLP {float(m.mean()) / 100:.2f} hPa, obs 1012.30 hPa")


# ---------------------------------------------------------------------------
def test_hourly_precipitation_scored():
    """
    S7: ASOS p01i becomes APCP1H observations (inches -> mm; a trace "T" is
    0 mm, "M" is no observation). A moist forecast is scored with the rain of
    the hour ending at each snapshot: the difference of the accumulated "rain"
    field from the snapshot an hour earlier. A forecast without rain skips them
    with a named reason. Two reports 1 km apart -- one wet, one dry -- must
    both pass QC (no buddy check on precipitation).
    """
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fc = fake_forecast(d / "forecast.npz", 3, terrain_height=6.0)
        z = dict(np.load(fc))
        acc = np.array([0.5, 2.5, 2.5])                      # mm since the start at +1, +2, +3 h
        z["rain"] = np.stack([np.full(z["pi"][0].shape, a) for a in acc]).astype(np.float32)
        np.savez_compressed(d / "moist.npz", **z)
        lines = ["station,valid,lon,lat,elevation,tmpf,p01i,drct,sknt"]
        p01 = {1: "0.04", 2: "T", 3: "M"}
        for h in range(1, 4):
            t = RUN + timedelta(hours=h)
            lines.append(f"KBOS,{t:%Y-%m-%d %H:%M},-71.01,42.36,6.0,59.0,{p01[h]},270,10")
            lines.append(f"KXYZ,{t:%Y-%m-%d %H:%M},-71.00,42.36,6.0,59.0,0.00,270,10")
        first, last = RUN + timedelta(hours=1), RUN + timedelta(hours=3)
        for name in ("dry", "moist"):
            paths = V.archive_paths(d / name, RUN)
            V.store_raw(paths["raw"], f"asos_{first:%Y%m%d%H}_{last:%Y%m%d%H}", "\n".join(lines) + "\n")
        m_dry, _ = V.verify(fc, d / "dry", run_time=RUN, report_only=True, verbose=False)
        m_wet, _ = V.verify(d / "moist.npz", d / "moist", run_time=RUN, report_only=True, verbose=False)
        pr = sorted([m for m in m_wet if m["variable"] == "APCP1H"], key=lambda m: (m["station"], m["lead_hours"]))
        bos = [m for m in pr if m["station"] == "KBOS"]
        ok = (len(pr) == 5 and len(bos) == 2
              and abs(bos[0]["observation"] - 1.016) < 1e-9 and abs(bos[0]["forecast"] - 0.5) < 1e-4
              and bos[1]["observation"] == 0.0 and abs(bos[1]["forecast"] - 2.0) < 1e-4
              and not any(m["variable"] == "APCP1H" for m in m_dry))
        report("hourly precipitation: APCP1H pairs (0.04 in -> 1.016 mm, trace -> 0, M dropped), "
               "forecast = hour's rain from the accumulation; wet next to dry both kept; none for a dry run", ok,
               f"{len(pr)} pairs; KBOS obs {[round(m['observation'], 3) for m in bos]} mm, "
               f"forecast {[round(m['forecast'], 3) for m in bos]} mm; dry-run pairs "
               f"{sum(m['variable'] == 'APCP1H' for m in m_dry)}")


# ---------------------------------------------------------------------------
def test_model_z0_operator():
    """
    Test P2 (P-78; model is the default since 2026-10-10): --z0-operator model
    reduces the lowest-level wind with the
    forecast's own roughness (its z0 field). With z0 = 1 m everywhere the
    10 m factor must be ln(10/1)/ln(z1/1), against ln(10/0.1)/ln(z1/0.1)
    with the fixed operator, at the same z1.
    """
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        fc = fake_forecast(d / "forecast.npz", 3, terrain_height=6.0)
        z = dict(np.load(fc))
        z["z0"] = np.full(z["terrain"].shape, 1.0, np.float32)
        np.savez_compressed(d / "z0.npz", **z)
        lines = ["station,valid,lon,lat,elevation,tmpf,drct,sknt"]
        for h in range(4):
            t = RUN + timedelta(hours=h)
            lines.append(f"KBOS,{t:%Y-%m-%d %H:%M},-71.01,42.36,6.0,59.0,270,10")
        first, last = RUN + timedelta(hours=1), RUN + timedelta(hours=3)
        for name in ("fixed", "model"):
            paths = V.archive_paths(d / name, RUN)
            V.store_raw(paths["raw"], f"asos_{first:%Y%m%d%H}_{last:%Y%m%d%H}", "\n".join(lines) + "\n")
        mf, _ = V.verify(d / "z0.npz", d / "fixed", run_time=RUN, report_only=True, verbose=False,
                         z0_operator="fixed")
        mm, _ = V.verify(d / "z0.npz", d / "model", run_time=RUN, report_only=True, verbose=False,
                         z0_operator="model")
        uf = [m for m in mf if m["variable"] == "UGRD"][0]
        um = [m for m in mm if m["variable"] == "UGRD"][0]
        z1 = um["model_level_agl_m"]
        want_f, want_m = np.log(10 / 0.1) / np.log(z1 / 0.1), np.log(10 / 1.0) / np.log(z1 / 1.0)
        ok = (abs(uf["wind_factor"] - want_f) < 1e-9 and abs(um["wind_factor"] - want_m) < 1e-6
              and abs(um["z0"] - 1.0) < 1e-6 and abs(uf["z0"] - 0.1) < 1e-12)
        report("z0 operator: model uses the forecast's z0 (factor ln(10/z0)/ln(z1/z0) at z0 = 1 m); fixed keeps 0.1 m",
               ok, f"z1 {z1:.1f} m; factor fixed {uf['wind_factor']:.4f} (want {want_f:.4f}), "
                   f"model {um['wind_factor']:.4f} (want {want_m:.4f})")


if __name__ == "__main__":
    print("\nVerification archiver\n" + "=" * 66)
    for fn in (test_raw_observations_are_stored_verbatim,
               test_archive_written_and_recomputable,
               test_rerun_does_not_duplicate,
               test_lead_hours_recorded,
               test_elevation_correction_is_recorded,
               test_report_only_refuses_without_cache,
               test_old_forecast_format_is_refused,
               test_auto_surface_operator,
               test_dewpoint_scored_for_moist_forecasts,
               test_sea_level_pressure_scored,
               test_hourly_precipitation_scored,
               test_model_z0_operator):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

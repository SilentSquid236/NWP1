"""
Tests for the observation analysis: sources, Barnes, heights, the builder,
and the two rules that make a run a forecast (nothing after the cycle time;
withheld stations never enter).

The parsers are tested against payloads captured from the live services on
2026-09-22 (testdata/), because every interface defect in this project's
history passed an offline suite built on what the AI assumed a service sent.

Run:  python src/analysis/test_analysis.py
"""

import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (HERE, ROOT, ROOT / "src", ROOT / "src" / "verification",
          ROOT / "src" / "dynamics"):
    sys.path.insert(0, str(p))

import config
import sources
import barnes
import geo
import build
from observations import Observation, run_qc

TD = HERE / "testdata"
CYCLE = datetime(2026, 9, 21, 12)
BOX = sources.analysis_box(config.DOMAIN)
results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def text(name):
    return (TD / name).read_text(encoding="utf-8")


# ---------------------------------------------------------------------------
# Sources
# ---------------------------------------------------------------------------

def test_asos_parse_window_and_one_report():
    obs = sources.parse_asos(text("asos_2026092112_sample.csv"), BOX, CYCLE)
    late = [o for o in obs if o.time > CYCLE]
    early = [o for o in obs if o.time < CYCLE - timedelta(minutes=60)]
    per = {}
    for o in obs:
        per[(o.station, o.variable)] = per.get((o.station, o.variable), 0) + 1
    report("ASOS: only reports in (cycle-60 min, cycle], one per station",
           obs and not late and not early and max(per.values()) == 1,
           f"{len(obs)} obs, {len(late)} after the cycle, max {max(per.values())} per station/variable")


def test_altimeter_to_pressure():
    p = sources.station_pressure_from_altimeter(29.92, 0.0)
    p_hi = sources.station_pressure_from_altimeter(29.92, 1500.0)
    report("altimeter 29.92 inHg at sea level -> 1013.2 hPa; lower aloft",
           abs(p - 101_320) < 30 and 84_000 < p_hi < 86_000,
           f"{p/100:.1f} hPa at 0 m, {p_hi/100:.1f} hPa at 1500 m")


def test_bad_pressure_is_range_rejected():
    o = Observation(time=CYCLE, lat=40.6, lon=-74.7, variable="PMSL",
                    value=10_137.0, source="asos", station="SMQ", error_std=150.0)
    kept, rej, _ = run_qc([o])
    report("a 101 hPa sea-level pressure (altimeter '3.00') is rejected",
           not kept and rej and "range" in rej[0].qc_reason, rej[0].qc_reason if rej else "kept")


def test_raob_table_skips_composites_and_closed():
    sites = sources.parse_raob_table(text("raob_network.geojson"), BOX)
    report("RAOB table: active sites only; RNK included, CHH and _ALY not",
           "KRNK" in sites and "KCHH" not in sites and not any(s.startswith("_") for s in sites),
           f"{len(sites)} sites: {', '.join(sorted(sites))}")


def test_raob_parse_live_sample():
    sites = sources.parse_raob_table(text("raob_network.geojson"), BOX)
    obs = sources.parse_raob(text("raob_KOKX_2026092112.csv"), sites, CYCLE)
    t500 = [o.value for o in obs if o.variable == "TMP" and o.pressure == 50_000]
    report("RAOB: KOKX 2026-09-21 12Z parses; 500 hPa T = 264.65 K",
           obs and t500 and abs(t500[0] - 264.65) < 1e-6 and all(o.pressure for o in obs),
           f"{len(obs)} values, T500 {t500}")


def test_raob_request_is_one_station():
    u = sources.raob_url("KOKX", CYCLE)
    report("RAOB request uses sts/ets and one station (P-54)",
           "sts=" in u and "ets=" in u and u.count("station=") == 1 and "ts1" not in u, u[-90:])


def test_ndbc_latest_not_after_cycle():
    meta = {"lat": 35.0, "lon": -75.4, "elev": 0.0}
    body = text("ndbc_41025_5day_sample.txt")
    cyc = datetime(2026, 9, 22, 12)
    obs = sources.parse_ndbc_5day(body, "41025", meta, cyc)
    ok = obs and all(0 <= (cyc - o.time).total_seconds() <= 3600 for o in obs)
    report("NDBC: the report nearest the cycle, at or before it",
           ok, f"{[(o.variable, round(o.value, 1)) for o in obs]} at {obs[0].time if obs else None}")


def test_collect_survives_a_failing_source():
    def boom(cycle, box, fetcher):
        raise ConnectionError("simulated outage")

    def empty(cycle, box, fetcher):
        return sources.SourceResult("empty", "unavailable", message="nothing")
    res = sources.collect(CYCLE, sources=[("boom", boom), ("empty", empty),
                                          ("skipme", boom)],
                          fetcher=object(), skip={"skipme"}, verbose=False)
    st = [r.status for r in res]
    report("one source failing does not stop the others; statuses kept apart",
           st == ["failed", "unavailable", "disabled"], f"{st}")


def test_no_model_sources():
    banned = ("hrrr", "gfs", "rap", "nam", "meteo", "mesoanalysis", "ecmwf")
    names = [n for n, _ in sources.SOURCES]
    urls = " ".join([sources.IEM, sources.NDBC_ACTIVE, sources.MRMS_BUCKET, geo.ERDDAP]).lower()
    report("no model-derived source in the registry",
           not any(b in n for n in names for b in banned) and not any(b in urls for b in banned),
           f"sources {names}")


# ---------------------------------------------------------------------------
# Analysis
# ---------------------------------------------------------------------------

def _grid(ny=40, nx=44):
    lat, lon = geo.cell_centres(config.DOMAIN, ny, nx)
    return lat, lon


def test_barnes_recovers_analytic_field():
    lat, lon = _grid()
    rng = np.random.default_rng(1)
    la, lo = rng.uniform(34, 50.5, 400), rng.uniform(-86, -62, 400)
    f = lambda a, b: 10 * np.sin(np.radians(b) * 8) + 5 * np.cos(np.radians(a) * 10)
    inc = barnes.barnes_increments(lat, lon, la, lo, f(la, lo), 1.0,
                                   2.5 * barnes.mean_spacing_km(la, lo), lam=0.05, passes=3)
    rel = np.sqrt(np.mean((inc - f(lat, lon)) ** 2)) / np.sqrt(np.mean(f(lat, lon) ** 2))
    report("Barnes recovers a smooth analytic field from 400 points", rel < 0.05,
           f"relative RMS error {rel:.3f}")


def test_barnes_keeps_first_guess_far_away():
    lat, lon = _grid()
    inc = barnes.barnes_increments(lat, lon, [38.0], [-80.0], [5.0], 1.0, 300.0)
    far = inc[geo.cell_centres(config.DOMAIN, *lat.shape)[0] > 45.0]
    near = float(geo.bilinear(inc, 38.0, -80.0, config.DOMAIN))
    report("far from every observation the increment vanishes",
           np.abs(far).max() < 0.05 and near > 3.0,
           f"near {near:.2f}, max beyond 45N {np.abs(far).max():.4f}")


def test_heights_isothermal_exact():
    p = np.asarray(config.PRESSURE_LEVELS, float) * 100
    T = np.full((p.size, 2, 2), 250.0)
    RH = np.zeros_like(T)
    z = build.hydrostatic_heights(T, RH, np.full((2, 2), 101_325.0), p)
    exact = build.RD * 250.0 / build.G0 * np.log(101_325.0 / p)
    err = np.abs(z[:, 0, 0] - exact).max()
    report("hypsometric heights exact for an isothermal dry atmosphere",
           err < 0.5, f"max error {err:.3f} m")


def test_superob_never_extrapolates():
    rows = [Observation(time=CYCLE, lat=40, lon=-75, variable="TMP", value=v,
                        source="raob", station="KXXX", pressure=pp, error_std=1.0)
            for pp, v in ((85_000, 280.0), (70_000, 272.0))]
    so = build.superob_soundings(rows, np.asarray(config.PRESSURE_LEVELS, float) * 100)
    lv = sorted(config.PRESSURE_LEVELS[k] for k, *_ in so["TMP"])
    report("sounding superobs stay inside the reported profile",
           lv and min(lv) >= 700 and max(lv) <= 850, f"levels {lv}")


def _synthetic_obs(cycle):
    rng = np.random.default_rng(2)
    obs = []
    for i in range(150):
        la, lo = rng.uniform(37, 47.5), rng.uniform(-82, -66)
        stn = f"S{i:03d}"
        obs += [Observation(time=cycle, lat=la, lon=lo, variable="TMP", value=290 - (la - 37),
                            source="asos", station=stn, elevation=100.0, error_std=1.5),
                Observation(time=cycle, lat=la, lon=lo, variable="PMSL", value=101_500.0,
                            source="asos", station=stn, elevation=100.0, error_std=150.0)]
    p = np.asarray(config.PRESSURE_LEVELS, float) * 100
    for j, (la, lo) in enumerate(((41, -74), (43, -77), (39, -79), (45, -70))):
        for pp in p:
            obs.append(Observation(time=cycle, lat=la, lon=lo, variable="TMP",
                                   value=float(build.standard_temperature(pp)), source="raob",
                                   station=f"R{j}", pressure=float(pp), error_std=1.0))
    return obs


def test_build_rejects_observations_after_the_cycle():
    obs = _synthetic_obs(CYCLE)
    obs.append(Observation(time=CYCLE + timedelta(minutes=5), lat=40, lon=-75, variable="TMP",
                           value=290.0, source="asos", station="LATE", elevation=10.0))
    ter = np.zeros(geo.grid_shape(config.DOMAIN, 40_000))
    try:
        build.build_analysis(CYCLE, obs, ter, config.DOMAIN, config.PRESSURE_LEVELS)
        report("an observation after the cycle time is refused", False, "accepted")
    except ValueError as e:
        report("an observation after the cycle time is refused", True, str(e)[:80])


def test_withheld_never_enter_and_state_is_sane():
    import json
    obs = _synthetic_obs(CYCLE)
    ter = np.zeros(geo.grid_shape(config.DOMAIN, 40_000))
    feats, meta = build.build_analysis(CYCLE, obs, ter, config.DOMAIN, config.PRESSURE_LEVELS)
    prov = json.loads(meta["provenance"])
    used = set(prov["platforms_used"].get("asos", []))
    wh = set(prov["withheld"])
    z = feats[4]
    report("withheld stations never enter; heights increase upward; shape [C,L,Y,X]",
           wh and not (used & wh) and (np.diff(z, axis=0) > 0).all()
           and feats.shape[:2] == (5, config.N_LEVELS),
           f"{len(wh)} withheld, {len(used)} used, shape {feats.shape}, "
           f"first guess {meta['first_guess']}")


def test_previous_forecast_first_guess():
    import tempfile
    from sigma import SigmaLevels
    lev = SigmaLevels(config.N_LEVELS)
    ny, nx = 6, 7
    ps = np.full((ny, nx), 100_000.0)
    pi = ps - lev.p_top
    p3 = lev.pressure(pi)
    T = build.standard_temperature(p3)
    theta = T * (100_000.0 / p3) ** (build.RD / 1004.6)
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "forecast.npz"
        np.savez(path, times_s=np.array([3600.0 * h for h in range(1, 7)]),
                 u=np.zeros((6, lev.nz, ny, nx)), v=np.zeros((6, lev.nz, ny, nx)),
                 theta=np.repeat(theta[None], 6, 0), pi=np.repeat(pi[None], 6, 0),
                 sigma=lev.sigma, p_top=lev.p_top)
        levels = np.asarray(config.PRESSURE_LEVELS, float) * 100
        bg = build.first_guess_from_forecast(path, 6.0, levels)
        missing = build.first_guess_from_forecast(path, 9.0, levels)
    k = [i for i, p in enumerate(levels) if p <= 95_000]
    err = np.abs(bg["TMP"][k, 0, 0] - build.standard_temperature(levels[k])).max()
    report("previous forecast -> pressure levels; a missing lead gives None, not a substitute",
           err < 0.6 and missing is None, f"max T error {err:.2f} K above ground")


# ---------------------------------------------------------------------------
# Forecast driver
# ---------------------------------------------------------------------------

class _FakeModel:
    def __init__(self, blow_at=None):
        self.u = np.zeros((2, 3, 3)); self.v = self.u.copy()
        self.theta = np.full_like(self.u, 300.0); self.pi = np.full((3, 3), 8e4)
        self.time = 0.0; self.blow_at = blow_at

    def max_dt(self):
        return 60.0

    def step(self, dt):
        self.time += dt
        if self.blow_at is not None and self.time >= self.blow_at:
            self.u += 500.0

    @property
    def surface_pressure(self):
        return self.pi

    def sigma_dot(self):
        return np.zeros(1)


class _NoRelax:
    def apply(self, model, ext):
        pass


class _NoDriver:
    def at(self, t):
        return {}


def test_forecast_stops_inside_the_hour_on_blowup():
    from forecast import run_forecast
    info = {}
    m = _FakeModel(blow_at=5400.0)
    snaps = run_forecast(m, _NoDriver(), _NoRelax(), 6 * 3600, dt=60.0,
                         output_every=3600, progress=False, info=info)
    report("a blow-up is caught within the hour, not at its end (P-52)",
           info["stopped"].startswith("diverged") and m.time < 7200 and len(snaps) == 1,
           f"{info['stopped']}; model time {m.time/3600:.2f} h; {len(snaps)} snapshot kept")


def test_forecast_deadline_keeps_hours_reached():
    from forecast import run_forecast
    info = {}
    snaps = run_forecast(_FakeModel(), _NoDriver(), _NoRelax(), 24 * 3600, dt=60.0,
                         output_every=3600, progress=False, deadline_s=0.0, info=info)
    report("the deadline stops the run and says so",
           info["stopped"].startswith("deadline"), f"{info['stopped']}, {len(snaps)} snapshots")


def test_driving_frames_refuses_mixed_sources():
    import tempfile
    from forecast import driving_frames
    with tempfile.TemporaryDirectory() as d:
        (Path(d) / "obs_analysis_f00.npz").write_bytes(b"")
        files, src = driving_frames(d)
        (Path(d) / "live_hrrr_f00.npz").write_bytes(b"")
        try:
            driving_frames(d)
            mixed = False
        except SystemExit:
            mixed = True
    report("an observation run has one frame; mixed obs/HRRR frames are refused",
           src == "observations" and len(files) == 1 and mixed, f"{src}, refused mixed: {mixed}")


def test_verification_waits_for_the_window():
    import tempfile
    from verify_pending import pending
    with tempfile.TemporaryDirectory() as d:
        run = Path(d) / "obs_20260921_12"
        run.mkdir()
        np.savez(run / "forecast.npz", times_s=np.array([3600.0 * h for h in range(1, 7)]))
        rt = datetime(2026, 9, 21, 12)
        early = pending(d, rt + timedelta(hours=7), 2.0)
        late = pending(d, rt + timedelta(hours=9), 2.0)
        (run / "verified.json").write_text("{}")
        done = pending(d, rt + timedelta(hours=30), 2.0)
    report("a run is verified only after its window + latency, and only once",
           (len(early[0]), len(early[1])) == (0, 1) and len(late[0]) == 1
           and done == ([], []), f"at +7 h {len(early[0])} due; at +9 h {len(late[0])} due; "
                                 f"after verified.json {done}")


if __name__ == "__main__":
    print("=" * 62)
    print("Observation analysis")
    print("=" * 62)
    for fn in (test_asos_parse_window_and_one_report, test_altimeter_to_pressure,
               test_bad_pressure_is_range_rejected, test_raob_table_skips_composites_and_closed,
               test_raob_parse_live_sample, test_raob_request_is_one_station,
               test_ndbc_latest_not_after_cycle, test_collect_survives_a_failing_source,
               test_no_model_sources, test_barnes_recovers_analytic_field,
               test_barnes_keeps_first_guess_far_away, test_heights_isothermal_exact,
               test_superob_never_extrapolates, test_build_rejects_observations_after_the_cycle,
               test_withheld_never_enter_and_state_is_sane, test_previous_forecast_first_guess,
               test_forecast_stops_inside_the_hour_on_blowup,
               test_forecast_deadline_keeps_hours_reached,
               test_driving_frames_refuses_mixed_sources,
               test_verification_waits_for_the_window):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

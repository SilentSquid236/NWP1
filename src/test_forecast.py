"""
End-to-end test of the forecast driver on SYNTHETIC HRRR-shaped input.

No network and no real HRRR file needed: we fabricate arrays with the same
shape and physical character, which is enough to exercise every join between
ingest, dynamics, boundaries, and output. The joins are where integration
bugs live -- units, channel order, array layout, coordinate conventions.

Run:  python src/test_forecast.py
"""

import sys, tempfile
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE / "dynamics"))

import config
from vertical import PressureLevels, theta_from_T, T_from_theta
from primitive3d import Primitive3D
from boundaries import BoundaryDriver
from forecast import (hrrr_to_model_state, build_grid, Relaxation3D,
                      run_forecast, state_to_boundary, load_state)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


NY, NX = 40, 44
LEVELS = PressureLevels(config.PRESSURE_LEVELS)


def synthetic_hrrr(t_offset=0.0):
    """HRRR-shaped [C, L, Y, X] array: realistic magnitudes, right layout."""
    nz = LEVELS.nz
    fields = np.zeros((len(config.CHANNELS), nz, NY, NX), dtype=np.float32)
    y = np.linspace(0, 1, NY)[:, None] * np.ones((1, NX))

    for k in range(nz):
        # Temperature falling with height and poleward, plus a slow trend.
        T = 288.0 - 55.0 * (1 - LEVELS.p[k] / LEVELS.p[0]) - 8.0 * y + t_offset
        fields[0, k] = T                      # TMP
        fields[1, k] = 50.0 + 20.0 * y        # RH
        fields[2, k] = 5.0 + 20.0 * (1 - LEVELS.p[k] / LEVELS.p[0])   # UGRD
        fields[3, k] = 0.0                    # VGRD
        fields[4, k] = 8000.0 * (1 - LEVELS.p[k] / LEVELS.p[0])       # HGT
    return fields


# ---------------------------------------------------------------------------
def test_temperature_to_theta_roundtrip():
    """
    T -> theta -> T must return the original. A slip here is a systematic
    temperature error of tens of kelvin that would look like broken physics.
    """
    fields = synthetic_hrrr()
    u, v, theta = hrrr_to_model_state(fields, LEVELS)
    T_back = T_from_theta(theta, LEVELS.p.reshape(-1, 1, 1))

    err = np.abs(T_back - fields[0]).max()
    ok = err < 1e-4
    report("temperature <-> potential temperature round-trips", ok,
           f"max error {err:.2e} K; theta at surface "
           f"{theta[0].mean():.1f} K, at top {theta[-1].mean():.1f} K")


# ---------------------------------------------------------------------------
def test_theta_increases_upward():
    """
    A realistic atmosphere is statically stable, so theta must increase with
    height even though temperature decreases. If this came out inverted the
    model would convect everywhere immediately.
    """
    fields = synthetic_hrrr()
    _, _, theta = hrrr_to_model_state(fields, LEVELS)
    profile = theta.mean(axis=(1, 2))

    ok = np.all(np.diff(profile) > 0)
    report("theta increases with height (statically stable)", ok,
           f"theta {profile[0]:.1f} -> {profile[-1]:.1f} K while T "
           f"{fields[0].mean(axis=(1,2))[0]:.1f} -> "
           f"{fields[0].mean(axis=(1,2))[-1]:.1f} K")


# ---------------------------------------------------------------------------
def test_grid_matches_data_shape():
    """Grid dimensions must be derived from the data, never assumed."""
    fields = synthetic_hrrr()
    gr = build_grid(fields, LEVELS)

    ok = (gr.nx == NX and gr.ny == NY and gr.edge_mode == "replicate"
          and 1e3 < gr.dx < 100e3 and gr.f0 > 0)
    report("grid built from data dimensions", ok,
           f"{gr.nx}x{gr.ny} (data {NX}x{NY}), dx {gr.dx/1000:.1f} km, "
           f"dy {gr.dy/1000:.1f} km, f0 {gr.f0:.2e}")


# ---------------------------------------------------------------------------
def test_relaxation_covers_perimeter_only():
    """The 3D relaxation must taper horizontally and apply at every level."""
    fields = synthetic_hrrr()
    gr = build_grid(fields, LEVELS)
    relax = Relaxation3D(gr, width=8)

    a = relax.alpha
    ok = (a.shape[0] == 1 and a[0, 0, NX // 2] > 0.99
          and a[0, NY // 2, NX // 2] == 0.0
          and 0.2 < relax.interior_fraction < 0.9)
    report("3D relaxation applies at the perimeter, not the interior", ok,
           f"edge alpha {a[0,0,NX//2]:.3f}, centre {a[0,NY//2,NX//2]:.3f}, "
           f"free interior {relax.interior_fraction:.0%}")


# ---------------------------------------------------------------------------
def test_forecast_runs_and_stays_finite():
    """
    The integration test that matters: initialise from HRRR-shaped data, run
    with boundary relaxation, and confirm the result is physically sane.
    """
    fields = synthetic_hrrr()
    gr = build_grid(fields, LEVELS)
    u, v, th = hrrr_to_model_state(fields, LEVELS)

    times, states = [], []
    for i in range(4):
        fi = synthetic_hrrr(t_offset=0.5 * i)     # slowly warming boundary
        ui, vi, thi = hrrr_to_model_state(fi, LEVELS)
        times.append(i * 3600.0)
        states.append(state_to_boundary(ui, vi, thi))
    driver = BoundaryDriver(times, states)

    model = Primitive3D(gr, LEVELS)
    model.u, model.v, model.theta = u, v, th
    relax = Relaxation3D(gr, width=8)

    th0 = model.theta.copy()
    snaps = run_forecast(model, driver, relax, 3 * 3600,
                         output_every=3600, progress=False)

    finite = np.isfinite(model.u).all() and np.isfinite(model.theta).all()
    umax = np.abs(model.u).max()
    th_range = (model.theta.min(), model.theta.max())
    changed = np.abs(model.theta - th0).max()

    ok = (len(snaps) == 3 and finite and umax < 200.0
          and 200.0 < th_range[0] and th_range[1] < 900.0 and changed > 1e-6)
    report("3 h forecast from HRRR-shaped input stays physical", ok,
           f"{len(snaps)} snapshots, max|u| {umax:.1f} m/s, "
           f"theta {th_range[0]:.0f}-{th_range[1]:.0f} K, "
           f"max change {changed:.3f} K")


# ---------------------------------------------------------------------------
def test_boundaries_hold_edges_to_driver():
    """
    The perimeter must track the driving data while the interior evolves
    freely -- that is the entire contract of a limited-area run.
    """
    fields = synthetic_hrrr()
    gr = build_grid(fields, LEVELS)
    u, v, th = hrrr_to_model_state(fields, LEVELS)

    warm = synthetic_hrrr(t_offset=5.0)         # boundary 5 K warmer
    uw, vw, thw = hrrr_to_model_state(warm, LEVELS)
    driver = BoundaryDriver([0.0], [state_to_boundary(uw, vw, thw)])

    model = Primitive3D(gr, LEVELS)
    model.u, model.v, model.theta = u, v, th
    relax = Relaxation3D(gr, width=8)
    run_forecast(model, driver, relax, 2 * 3600, output_every=3600,
                 progress=False)

    edge_err = np.abs(model.theta[:, 0, :] - thw[:, 0, :]).max()
    interior_diff = np.abs(model.theta[:, NY // 2, NX // 2]
                           - thw[:, NY // 2, NX // 2]).max()

    ok = edge_err < 0.1 and interior_diff > 0.5
    report("edges follow the driver, interior stays free", ok,
           f"edge differs from driver by {edge_err:.4f} K; "
           f"domain centre still differs by {interior_diff:.2f} K")


# ---------------------------------------------------------------------------
def test_npz_roundtrip():
    """Ingest output must load back with channels and levels intact."""
    fields = synthetic_hrrr()
    with tempfile.TemporaryDirectory() as d:
        path = Path(d) / "live_hrrr_f00.npz"
        np.savez_compressed(path, hrrr_features=fields,
                            channels=np.array(config.CHANNELS),
                            levels_hPa=np.array(config.PRESSURE_LEVELS))
        back, meta = load_state(path)

        ok = (back.shape == fields.shape and np.allclose(back, fields)
              and list(meta["channels"]) == list(config.CHANNELS)
              and len(meta["levels_hPa"]) == config.N_LEVELS)
        report("ingested .npz round-trips with metadata", ok,
               f"shape {back.shape}, channels {list(meta['channels'])}, "
               f"{len(meta['levels_hPa'])} levels")


if __name__ == "__main__":
    print("\nForecast driver integration\n" + "=" * 62)
    for fn in (test_temperature_to_theta_roundtrip,
               test_theta_increases_upward,
               test_grid_matches_data_shape,
               test_relaxation_covers_perimeter_only,
               test_forecast_runs_and_stays_finite,
               test_boundaries_hold_edges_to_driver,
               test_npz_roundtrip):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

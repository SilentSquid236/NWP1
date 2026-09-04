"""
End-to-end test of the forecast driver on SYNTHETIC HRRR-shaped input.

No network and no real HRRR file needed: we fabricate arrays with the same
shape and physical character, which is enough to exercise every join between
ingest, conversion, dynamics, boundaries and output. The joins are where
integration bugs live -- units, channel order, array layout, coordinate
conventions.

WHAT CHANGED. The driver now runs the SIGMA core over real terrain, so these
tests cover the pressure-to-sigma conversion and the terrain contract as well.
The old versions of the first two tests checked a T <-> theta round trip
through pressure levels; the equivalent guarantee now lives in
`dynamics/test_interpolate.py`, which can check it against a standard
atmosphere where the right answer is known in closed form.

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
from vertical import PressureLevels
from sigma import SigmaLevels, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from boundaries import BoundaryDriver
from forecast import (hrrr_channels, hrrr_to_sigma_state, load_terrain,
                      build_grid, Relaxation3D, run_forecast,
                      state_to_boundary, load_state)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


NY, NX = 40, 44
LEVELS = PressureLevels(config.PRESSURE_LEVELS)
LEV = SigmaLevels(config.N_LEVELS)


def synthetic_hrrr(t_offset=0.0):
    """
    HRRR-shaped [C, L, Y, X] array: realistic magnitudes, right layout.

    Heights come from a standard atmosphere rather than a linear ramp, because
    the sigma conversion locates the surface by interpolating pressure against
    HEIGHT -- a height field that is only roughly right would make this test
    pass for the wrong reason.
    """
    nz = LEVELS.nz
    fields = np.zeros((len(config.CHANNELS), nz, NY, NX), dtype=np.float32)
    y = np.linspace(0, 1, NY)[:, None] * np.ones((1, NX))

    T0, L = 288.15, 0.0065
    for k in range(nz):
        p = LEVELS.p[k]
        T_std = T0 * (p / 101325.0) ** (RD * L / G0)
        z_std = (T0 - T_std) / L
        fields[0, k] = T_std - 8.0 * y + t_offset          # TMP
        fields[1, k] = 50.0 + 20.0 * y                     # RH
        fields[2, k] = 5.0 + 20.0 * (1 - p / LEVELS.p[0])  # UGRD
        fields[3, k] = 0.0                                 # VGRD
        fields[4, k] = z_std - 8.0 * y * 30.0              # HGT
    return fields


def synthetic_terrain():
    """A ridge, peaking near the 2 km the domain actually contains."""
    y = np.linspace(-1, 1, NY)[:, None]
    x = np.linspace(-1, 1, NX)[None, :]
    return 1800.0 * np.exp(-(x ** 2 + y ** 2) / 0.25)


def write_run(d, n_frames=4, terrain=True):
    """A run directory the driver would accept."""
    d = Path(d)
    for i in range(n_frames):
        np.savez_compressed(
            d / f"live_hrrr_f{i:02d}.npz",
            hrrr_features=synthetic_hrrr(t_offset=0.5 * i),
            channels=np.array(config.CHANNELS),
            levels_hPa=np.array(config.PRESSURE_LEVELS))
    if terrain:
        np.savez_compressed(d / "terrain.npz",
                            terrain=synthetic_terrain().astype(np.float32))
    return d


# ---------------------------------------------------------------------------
def test_channels_require_height():
    """
    HGT is not optional any more: the conversion locates the surface with it.
    A missing height field must say so, not fail later inside an interpolation.
    """
    fields = synthetic_hrrr()
    try:
        hrrr_channels(fields, channels=["TMP", "RH", "UGRD", "VGRD", "PRES"])
        ok, why = False, "accepted a channel list with no HGT"
    except KeyError as e:
        ok = "HGT" in str(e)
        why = f"refused with: {str(e)[:70]}..."
    report("a channel list without HGT is refused", ok, why)


# ---------------------------------------------------------------------------
def test_conversion_is_statically_stable():
    """
    A realistic atmosphere is statically stable, so theta must increase with
    height even though temperature decreases. Inverted, the model would
    convect everywhere on the first step -- and now that convective adjustment
    exists, it would do so silently rather than blowing up.
    """
    fields = synthetic_hrrr()
    pi, u, v, th = hrrr_to_sigma_state(fields, LEV, synthetic_terrain())
    profile = th.mean(axis=(1, 2))
    T = fields[0].mean(axis=(1, 2))
    # index 0 is the model TOP, so theta must decrease with index
    ok = np.all(np.diff(profile) < 0) and np.all(np.isfinite(th))
    report("the converted state is statically stable", ok,
           f"theta {profile[-1]:.1f} K at the ground -> {profile[0]:.1f} K at "
           f"the lid, while T runs {T[0]:.1f} -> {T[-1]:.1f} K")


# ---------------------------------------------------------------------------
def test_surface_pressure_follows_terrain():
    """
    Surface pressure must be lower over the ridge than over the lowland, by
    roughly the hydrostatic amount. This is the join that would silently break
    if terrain and fields were ingested at different strides.
    """
    terrain = synthetic_terrain()
    pi, u, v, th = hrrr_to_sigma_state(synthetic_hrrr(), LEV, terrain)
    ps = pi + LEV.p_top
    lo = ps[terrain < 50.0].mean()
    hi = ps[terrain > terrain.max() - 50.0].mean()
    # Hydrostatic estimate over the ridge height, ~11 Pa/m near the surface.
    want = lo * np.exp(-G0 * terrain.max() / (RD * 280.0))
    ok = hi < lo and abs(hi - want) / want < 0.03
    report("surface pressure falls over the ridge", ok,
           f"lowland {lo/100:.0f} hPa, ridge {hi/100:.0f} hPa over "
           f"{terrain.max():.0f} m; hydrostatic estimate {want/100:.0f} hPa")


# ---------------------------------------------------------------------------
def test_missing_terrain_is_refused():
    """
    A run directory with no terrain must stop the forecast, not silently
    substitute flat ground. A flat Northeast is a different experiment, and
    running it by accident is how a result gets misread later.
    """
    with tempfile.TemporaryDirectory() as d:
        write_run(d, n_frames=1, terrain=False)
        try:
            load_terrain(d, (NY, NX))
            ok, why = False, "accepted a run directory with no terrain"
        except SystemExit as e:
            ok = "terrain" in str(e).lower()
            why = "refused with a message naming the ingest command"
    report("a run with no terrain is refused, not flattened", ok, why)


# ---------------------------------------------------------------------------
def test_terrain_shape_mismatch_is_caught():
    """Terrain and fields ingested at different strides must not be combined."""
    with tempfile.TemporaryDirectory() as d:
        d = Path(d)
        np.savez_compressed(d / "terrain.npz",
                            terrain=np.zeros((NY // 2, NX // 2), np.float32))
        try:
            load_terrain(d, (NY, NX))
            ok, why = False, "accepted mismatched shapes"
        except SystemExit as e:
            ok = "stride" in str(e)
            why = "refused and named --stride as the likely cause"
    report("a terrain/field shape mismatch is caught", ok, why)


# ---------------------------------------------------------------------------
def test_grid_matches_data_shape():
    """Grid dimensions must be derived from the data, never assumed."""
    gr = build_grid(synthetic_hrrr(), LEVELS)
    ok = (gr.nx == NX and gr.ny == NY and gr.edge_mode == "replicate"
          and 1e3 < gr.dx < 100e3 and gr.f0 > 0)
    report("grid built from data dimensions", ok,
           f"{gr.nx}x{gr.ny} (data {NX}x{NY}), dx {gr.dx/1000:.1f} km, "
           f"dy {gr.dy/1000:.1f} km, f0 {gr.f0:.2e}")


# ---------------------------------------------------------------------------
def test_relaxation_covers_perimeter_only():
    """The 3D relaxation must taper horizontally and apply at every level."""
    gr = build_grid(synthetic_hrrr(), LEVELS)
    relax = Relaxation3D(gr, width=8)
    a = relax.alpha
    ok = (a.shape[0] == 1 and a[0, 0, NX // 2] > 0.99
          and a[0, NY // 2, NX // 2] == 0.0
          and 0.2 < relax.interior_fraction < 0.9)
    report("3D relaxation applies at the perimeter, not the interior", ok,
           f"edge alpha {a[0,0,NX//2]:.3f}, centre {a[0,NY//2,NX//2]:.3f}, "
           f"free interior {relax.interior_fraction:.0%}")


# ---------------------------------------------------------------------------
def test_relaxation_drives_surface_pressure():
    """
    Surface pressure is PROGNOSTIC in the sigma core, so the relaxation has to
    hold it at the edges too. Leaving it free while relaxing the wind drives
    the boundary column toward a mass field the incoming flow does not
    support.
    """
    gr = build_grid(synthetic_hrrr(), LEVELS)
    terrain = synthetic_terrain()
    m = PrimitiveSigma(gr, LEV, terrain=terrain)
    pi, u, v, th = hrrr_to_sigma_state(synthetic_hrrr(), LEV, terrain)
    # COPY. `Relaxation3D.apply` updates in place, so handing the model the
    # same array the test then compares against would make the test measure
    # nothing -- it read "edge moved +0 Pa" until this copy was added. The
    # driver has the same hazard and now assigns copies for the same reason.
    m.pi, m.u, m.v, m.theta = pi.copy(), u.copy(), v.copy(), th.copy()

    target = pi - 500.0
    relax = Relaxation3D(gr, width=8)
    relax.apply(m, state_to_boundary(u, v, th, target))

    edge = m.pi[0, NX // 2] - pi[0, NX // 2]
    mid = m.pi[NY // 2, NX // 2] - pi[NY // 2, NX // 2]
    ok = edge < -400.0 and abs(mid) < 1e-9
    report("relaxation drives surface pressure at the edge only", ok,
           f"edge moved {edge:+.0f} Pa toward the driver, "
           f"centre moved {mid:+.1f} Pa")


# ---------------------------------------------------------------------------
def test_forecast_runs_and_stays_finite():
    """
    The integration test that matters: initialise from HRRR-shaped data over
    terrain, run with boundary relaxation, confirm the result is sane.
    """
    fields = synthetic_hrrr()
    terrain = synthetic_terrain()
    gr = build_grid(fields, LEVELS)
    pi, u, v, th = hrrr_to_sigma_state(fields, LEV, terrain)

    times, states = [], []
    for i in range(4):
        pi_i, ui, vi, thi = hrrr_to_sigma_state(
            synthetic_hrrr(t_offset=0.5 * i), LEV, terrain)
        times.append(i * 3600.0)
        states.append(state_to_boundary(ui, vi, thi, pi_i))
    driver = BoundaryDriver(times, states)

    model = PrimitiveSigma(gr, LEV, terrain=terrain)
    model.pi, model.u, model.v, model.theta = pi, u, v, th
    relax = Relaxation3D(gr, width=8)

    th0 = model.theta.copy()
    snaps = run_forecast(model, driver, relax, 3 * 3600,
                         output_every=3600, progress=False)

    finite = np.isfinite(model.u).all() and np.isfinite(model.theta).all()
    umax = float(np.abs(model.u).max())
    ps = model.surface_pressure
    changed = float(np.abs(model.theta - th0).max())

    ok = (len(snaps) == 3 and finite and umax < 200.0
          and 200.0 < model.theta.min() and model.theta.max() < 900.0
          and changed > 1e-6 and ps.min() > 5e4)
    report("3 h forecast from HRRR-shaped input over terrain stays physical",
           ok,
           f"{len(snaps)} snapshots, max|u| {umax:.1f} m/s, "
           f"p_s {ps.min()/100:.0f}-{ps.max()/100:.0f} hPa, "
           f"max theta change {changed:.3f} K")


# ---------------------------------------------------------------------------
def test_boundaries_hold_edges_to_driver():
    """
    The perimeter must track the driving data while the interior evolves
    freely -- the entire contract of a limited-area run.
    """
    fields = synthetic_hrrr()
    terrain = synthetic_terrain()
    gr = build_grid(fields, LEVELS)
    pi, u, v, th = hrrr_to_sigma_state(fields, LEV, terrain)

    pw, uw, vw, thw = hrrr_to_sigma_state(synthetic_hrrr(t_offset=5.0),
                                          LEV, terrain)
    driver = BoundaryDriver([0.0], [state_to_boundary(uw, vw, thw, pw)])

    model = PrimitiveSigma(gr, LEV, terrain=terrain)
    model.pi, model.u, model.v, model.theta = pi, u, v, th
    relax = Relaxation3D(gr, width=8)
    run_forecast(model, driver, relax, 2 * 3600, output_every=3600,
                 progress=False)

    edge_err = float(np.abs(model.theta[:, 0, :] - thw[:, 0, :]).max())
    interior = float(np.abs(model.theta[:, NY // 2, NX // 2]
                            - thw[:, NY // 2, NX // 2]).max())
    ok = edge_err < 0.1 and interior > 0.5
    report("edges follow the driver, interior stays free", ok,
           f"edge differs from driver by {edge_err:.4f} K; "
           f"domain centre still differs by {interior:.2f} K")


# ---------------------------------------------------------------------------
def test_npz_roundtrip():
    """Ingest output must load back with channels and levels intact."""
    with tempfile.TemporaryDirectory() as d:
        write_run(d, n_frames=1)
        back, meta = load_state(Path(d) / "live_hrrr_f00.npz")
        terrain, p_sfc = load_terrain(d, (NY, NX))
        ok = (back.shape == (len(config.CHANNELS), config.N_LEVELS, NY, NX)
              and list(meta["channels"]) == list(config.CHANNELS)
              and len(meta["levels_hPa"]) == config.N_LEVELS
              and terrain.shape == (NY, NX) and p_sfc is None)
        report("run directory round-trips with metadata and terrain", ok,
               f"fields {back.shape}, terrain {terrain.shape} "
               f"{terrain.min():.0f}-{terrain.max():.0f} m, "
               f"surface pressure absent (derived from heights instead)")


if __name__ == "__main__":
    print("\nForecast driver integration\n" + "=" * 62)
    for fn in (test_channels_require_height,
               test_conversion_is_statically_stable,
               test_surface_pressure_follows_terrain,
               test_missing_terrain_is_refused,
               test_terrain_shape_mismatch_is_caught,
               test_grid_matches_data_shape,
               test_relaxation_covers_perimeter_only,
               test_relaxation_drives_surface_pressure,
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

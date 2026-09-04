"""
Validation for the sigma-coordinate observation operator.

The decisive test is the first one: a surface observation must be compared
against the BOTTOM of the model column. The base-class operator returned the
lid for exactly this call, which is a seventy-kelvin error that looks like a
plausible number.

Run:  python src/verification/test_sigma_operator.py
"""
import sys
from pathlib import Path

import numpy as np
np.seterr(all="ignore")

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE.parent))
sys.path.insert(0, str(HERE.parent.parent))
sys.path.insert(0, str(HERE.parent / "dynamics"))

import config
from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA, RD, G0
from obs_operator import GridInterpolator
from sigma_operator import SigmaInterpolator

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


DOMAIN = config.DOMAIN
NY, NX = 40, 44


def setup(terrain_height=0.0, ny=NY, nx=NX):
    lat0 = 0.5 * (DOMAIN["lat_min"] + DOMAIN["lat_max"])
    dy = (DOMAIN["lat_max"] - DOMAIN["lat_min"]) * 111_132.0 / ny
    dx = (DOMAIN["lon_max"] - DOMAIN["lon_min"]) * 111_320.0 * \
        np.cos(np.radians(lat0)) / nx
    gr = CGrid(nx, ny, dx, dy, f0=9.81e-5, beta=1.69e-11, edge_mode="replicate")
    lev = SigmaLevels(20)
    terrain = np.full((ny, nx), terrain_height)
    pi = 101325.0 * np.exp(-G0 * terrain / (RD * 280.0)) - lev.p_top
    op = SigmaInterpolator(gr, DOMAIN, lev.sigma, lev.p_top, pi, terrain)
    return gr, lev, op, terrain, pi


def centre():
    return (0.5 * (DOMAIN["lat_min"] + DOMAIN["lat_max"]),
            0.5 * (DOMAIN["lon_min"] + DOMAIN["lon_max"]))


def theta_field(lev, pi, T_surface=288.0, lapse=6.5e-3):
    """A statically stable column, as potential temperature."""
    p = lev.pressure(pi)
    T0 = 288.15
    T = T0 * (p / 101325.0) ** (RD * lapse / G0)
    return T / (p / P0) ** KAPPA, T


# ---------------------------------------------------------------------------
def test_surface_observation_uses_the_bottom():
    """
    THE ONE THAT MATTERS. A surface observation has no pressure, and must be
    matched against the lowest model level. The base class returns field3d[0],
    which in this project is the LID.
    """
    gr, lev, op, terrain, pi = setup()
    th, T = theta_field(lev, pi)
    lat, lon = centre()

    ours = op.at_observation(T, lat, lon, None)
    base = GridInterpolator(gr, DOMAIN).at_observation(T, lat, lon, None)

    ok = (ours is not None and abs(ours - T[-1].mean()) < 1.0
          and abs(ours - base) > 50.0)
    report("a surface observation is matched to the bottom of the column", ok,
           f"sigma operator {ours:.1f} K (model lowest level "
           f"{T[-1].mean():.1f} K); base class would return {base:.1f} K "
           f"(the lid) -- a {abs(ours-base):.0f} K error")


# ---------------------------------------------------------------------------
def test_pressure_interpolation_recovers_levels():
    """Asking for a model level's own pressure must return that level."""
    gr, lev, op, terrain, pi = setup()
    th, T = theta_field(lev, pi)
    lat, lon = centre()
    p = op.pressure_column(lat, lon)

    errs = [abs(op.at_pressure(T, lat, lon, p[k]) - T[k].mean())
            for k in (2, 8, 14, 19)]
    ok = max(errs) < 0.05
    report("interpolating to a model level returns that level", ok,
           f"max error {max(errs):.4f} K across four levels")


# ---------------------------------------------------------------------------
def test_column_pressures_follow_terrain():
    """
    The point of sigma: two columns have DIFFERENT level pressures. An
    operator that used one shared pressure array would be wrong over terrain.
    """
    _, lev, op_flat, _, _ = setup(0.0)
    _, _, op_high, _, _ = setup(1500.0)
    lat, lon = centre()
    p_flat = op_flat.pressure_column(lat, lon)
    p_high = op_high.pressure_column(lat, lon)

    ok = np.all(p_high[1:] < p_flat[1:]) and p_high[-1] < p_flat[-1] - 5000.0
    report("level pressures follow the terrain", ok,
           f"lowest level {p_flat[-1]/100:.0f} hPa over flat ground, "
           f"{p_high[-1]/100:.0f} hPa over 1500 m")


# ---------------------------------------------------------------------------
def test_out_of_range_pressure_returns_none():
    """
    A sounding level below the surface or above the lid must return None, not
    an extrapolated number. A fabricated match is worse than no match.
    """
    gr, lev, op, terrain, pi = setup()
    th, T = theta_field(lev, pi)
    lat, lon = centre()
    below = op.at_pressure(T, lat, lon, 105_000.0)
    above = op.at_pressure(T, lat, lon, 10_000.0)
    ok = below is None and above is None
    report("out-of-column pressures return None", ok,
           f"1050 hPa -> {below}, 100 hPa -> {above}")


# ---------------------------------------------------------------------------
def test_outside_domain_returns_none():
    """Observations outside the domain must not be matched."""
    gr, lev, op, terrain, pi = setup()
    th, T = theta_field(lev, pi)
    ok = (op.at_observation(T, 10.0, -75.0, None) is None
          and op.at_observation(T, 42.0, -120.0, None) is None)
    report("observations outside the domain return None", ok,
           "a station in the tropics and one in the Rockies both refused")


# ---------------------------------------------------------------------------
def test_elevation_correction_has_the_right_sign():
    """
    A station BELOW the model level must come out WARMER. Getting this
    backwards doubles the error instead of removing it, and would look like a
    model cold bias at every valley station.
    """
    gr, lev, op, terrain, pi = setup(terrain_height=500.0)
    th, T = theta_field(lev, pi)
    lat, lon = centre()

    z1 = op.lowest_level_height(th, lat, lon)
    warm, iw = op.station_temperature(th, lat, lon, z1 - 300.0)
    cold, ic = op.station_temperature(th, lat, lon, z1 + 300.0)
    plain = op.temperature(th, lat, lon)

    ok = (warm > plain > cold
          and abs((warm - plain) - 300.0 * 6.5e-3) < 1e-6
          and abs(iw["elev_correction_m"] + 300.0) < 1e-6)
    report("elevation correction warms low stations, cools high ones", ok,
           f"model level at {z1:.0f} m gives {plain:.2f} K; "
           f"300 m lower {warm:.2f} K, 300 m higher {cold:.2f} K "
           f"(6.5 K/km)")


# ---------------------------------------------------------------------------
def test_correction_size_is_recorded():
    """
    The archive has to be able to tell a 5 m correction from a 400 m one
    later, or a bias caused by the operator is indistinguishable from a bias
    in the model.
    """
    gr, lev, op, terrain, pi = setup(terrain_height=800.0)
    th, T = theta_field(lev, pi)
    lat, lon = centre()
    val, info = op.station_temperature(th, lat, lon, 100.0)
    ok = ("elev_correction_m" in info and "elev_correction_K" in info
          and info["elev_correction_m"] < -500.0
          and abs(info["elev_correction_K"]) > 3.0)
    report("the size of the elevation correction is recorded", ok,
           f"station 100 m, model level {info.get('model_level_height_m',0):.0f} m "
           f"-> correction {info.get('elev_correction_m',0):+.0f} m = "
           f"{info.get('elev_correction_K',0):+.2f} K")


if __name__ == "__main__":
    print("\nSigma observation operator\n" + "=" * 66)
    for fn in (test_surface_observation_uses_the_bottom,
               test_pressure_interpolation_recovers_levels,
               test_column_pressures_follow_terrain,
               test_out_of_range_pressure_returns_none,
               test_outside_domain_returns_none,
               test_elevation_correction_has_the_right_sign,
               test_correction_size_is_recorded):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

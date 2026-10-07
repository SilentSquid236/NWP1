"""
Validation for the prescribed diurnal surface heat flux (P-59).

Run:  python test_diurnal.py
"""
import math
from datetime import datetime, timedelta

import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels, CP, G0, P0, KAPPA
from grid import CGrid
from primitive_sigma import PrimitiveSigma
from diurnal import (DiurnalHeating, sin_solar_elevation, surface_heat_flux,
                     solar_declination, H_NIGHT, F_SENSIBLE, S0, TAU)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def test_noon_elevation_albany():
    """Albany, 22 Sep 2026: noon elevation 90 - lat + declination, near 16:55Z."""
    day = datetime(2026, 9, 22)
    times = [day + timedelta(minutes=m) for m in range(24 * 60)]
    e = np.degrees(np.arcsin([float(sin_solar_elevation(42.75, -73.80, t)) for t in times]))
    k = int(np.argmax(e))
    dec = math.degrees(solar_declination(day.timetuple().tm_yday))
    expect = 90.0 - 42.75 + dec
    noon = times[k]
    ok = abs(e[k] - expect) < 0.05 and abs((noon - day.replace(hour=16, minute=55)).total_seconds()) <= 120
    report("noon solar elevation at Albany is 90 - latitude + declination", ok,
           f"max {e[k]:.2f} deg at {noon:%H:%MZ}; expected {expect:.2f} deg near 16:55Z "
           f"(declination {dec:+.2f} deg)")


def test_night_flux_land_and_water():
    lat = np.full((2, 2), 42.0); lon = np.full((2, 2), -74.0)
    land = np.array([[True, False], [True, False]])
    s = sin_solar_elevation(lat, lon, datetime(2026, 9, 23, 6))
    H = surface_heat_flux(s, land)
    ok = bool((s < 0).all()) and np.allclose(H[land], H_NIGHT) and np.all(H[~land] == 0.0)
    report("at 06Z: night cooling over land, nothing over water", ok,
           f"sin e {s.max():.2f}; flux land {H[land].mean():+.1f}, water {H[~land].mean():+.1f} W/m2")


def test_layer_energy_closes():
    """cp * dT * dp1 / g must equal H * dt for the lowest layer."""
    lev = SigmaLevels(20)
    pi = np.full((3, 4), 95000.0 - lev.p_top)
    lat = np.full((3, 4), 42.0); lon = np.full((3, 4), -74.0); land = np.ones((3, 4), bool)
    dh = DiurnalHeating(lat, lon, land, datetime(2026, 9, 22, 17))
    dt = 600.0
    dth = dh.theta_tendency(0.0, pi, lev) * dt
    p1 = lev.p_top + lev.sigma[-1] * pi
    dT = dth * (p1 / P0) ** KAPPA
    energy = CP * dT * lev.dsigma[-1] * pi / G0
    H = dh.flux(0.0)
    err = float(np.abs(energy - H * dt).max() / np.abs(H * dt).max())
    report("the lowest-layer heating closes the energy budget", err < 1e-12 and H.min() > 50,
           f"flux {H.mean():.1f} W/m2 at 17Z, {dth.mean()*3600/dt:.2f} K/h in the lowest layer; "
           f"relative error {err:.1e}")


def test_daily_mean_flux_magnitude():
    """Late-September daily mean at 42N: small and positive, as observed."""
    day = datetime(2026, 9, 25)
    land = np.ones((1, 1), bool)
    H = [surface_heat_flux(sin_solar_elevation(42.0, -74.0, day + timedelta(minutes=m)), land).item()
         for m in range(0, 24 * 60, 10)]
    mean, peak = float(np.mean(H)), float(np.max(H))
    ok = 0.0 < mean < 60.0 and 80.0 < peak < 200.0
    report("daily mean flux small and positive, noon peak 80-200 W/m2", ok,
           f"daily mean {mean:+.1f} W/m2, peak {peak:.0f} W/m2, minimum {min(H):+.0f}")


def test_model_column_budget():
    """
    In the model, after the adjustment has mixed the heated layer, the extra
    column heat over land must equal the prescribed flux; water columns must
    be untouched (to the motion the land-sea contrast starts within 10 steps).
    """
    ny, nx = 12, 16
    gr = CGrid(nx, ny, 12e3, 12e3, f0=1e-4, beta=0.0, edge_mode="replicate")
    lev = SigmaLevels(20)
    pi0 = np.full((ny, nx), 100000.0 - lev.p_top)
    theta0 = np.linspace(340.0, 290.0, 20)[:, None, None] * np.ones((1, ny, nx))
    land = np.zeros((ny, nx), bool); land[:, : nx // 2] = True
    lat = np.full((ny, nx), 42.0); lon = np.full((ny, nx), -74.0)
    out = []
    for heat in (False, True):
        m = PrimitiveSigma(gr, lev)
        m.conv_scheme = "pav"
        m.pi, m.u, m.v, m.theta = pi0.copy(), np.zeros_like(theta0), np.zeros_like(theta0), theta0.copy()
        if heat:
            m.surface_heating = DiurnalHeating(lat, lon, land, datetime(2026, 9, 22, 17),
                                               f_sensible=1.0)   # strong, to force mixing
        rates = []
        for _ in range(10):
            if heat:
                rates.append(m.surface_heating.theta_tendency(m.time + 0.5 * 30.0, m.pi, lev))
            m.step(30.0)
        out.append((m.theta.copy(), m.pi.copy(), rates))
    (t0, p0, _), (t1, p1, rates) = out
    dm = lev.dsigma[:, None, None] * p1[None]
    gained = ((t1 - t0) * dm).sum(0)
    expect = sum(r * 30.0 for r in rates) * lev.dsigma[-1] * p1
    c = (slice(2, -2), slice(2, nx // 2 - 2))         # land, away from the coast
    w = (slice(2, -2), slice(nx // 2 + 2, -2))        # water, away from the coast
    rel = float(np.abs(gained[c] - expect[c]).max() / np.abs(expect[c]).max())
    mixed = int((np.abs(t1[-2][c] - t0[-2][c]) > 1e-6).sum())
    water = float(np.abs(t1[:, w[0], w[1]] - t0[:, w[0], w[1]]).max())
    ok = rel < 0.02 and mixed > 0 and water < 0.05
    report("model: land columns gain exactly the prescribed heat, mixed upward; water does not", ok,
           f"column budget error {rel:.1e}; layer above heated in {mixed} land columns; "
           f"max water change {water:.1e} K")


if __name__ == "__main__":
    print("\nDiurnal surface heating\n" + "=" * 66)
    for fn in (test_noon_elevation_albany, test_night_flux_land_and_water,
               test_layer_energy_closes, test_daily_mean_flux_magnitude,
               test_model_column_budget):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

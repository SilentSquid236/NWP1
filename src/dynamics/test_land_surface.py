"""
Tests for the force-restore ground temperature (land_surface.py).

Run:  python test_land_surface.py
"""
import math
from datetime import datetime

import numpy as np

from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from surface import surface_drag, lowest_level_height
from land_surface import (ForceRestoreSurface, force_restore_coefficient,
                          brutsaert_emissivity, louis_momentum, louis_heat,
                          EPS_A)

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# 1. Constants are the stated ones -------------------------------------------
cg = force_restore_coefficient()
eps = brutsaert_emissivity(12.0, 285.0)
report("force-restore coefficient and clear-sky emissivity",
       abs(cg - 8.53e-6) < 0.02e-6 and abs(eps - EPS_A) < 0.005,
       f"C_G = {cg:.3e} K m^2/J (expect 8.53e-6); Brutsaert(12 hPa, 285 K) = {eps:.3f} "
       f"against EPS_A = {EPS_A}")

# 2. Louis functions -----------------------------------------------------------
a2 = (0.4 / math.log(230.0 / 0.1)) ** 2
zz = np.array(2300.0)
ri = np.array([-1.0, -0.1, 0.0, 0.1, 0.5, 1.0, 5.0])
fm = louis_momentum(ri, a2, zz)
fh = louis_heat(ri, a2, zz)
report("Louis (1979): neutral 1, unstable > 1, stable falling with a long tail",
       abs(fm[2] - 1) < 1e-12 and abs(fh[2] - 1) < 1e-12 and fm[0] > fm[1] > 1
       and np.all(np.diff(fm[2:]) < 0) and fm[-1] > 0 and np.all(fh[3:] < fm[3:]),
       "F_m " + " ".join(f"{x:.3f}" for x in fm) + " | F_h " + " ".join(f"{x:.3f}" for x in fh))


# 3. A column at night and at noon ----------------------------------------------
def column(hour_utc, n_hours, u=5.0, T_air=285.0, lon=-74.0, land=True, long_tail=True):
    lat = np.array([[42.0]]); lo = np.array([[lon]])
    s = ForceRestoreSurface(lat, lo, np.array([[land]]), datetime(2026, 9, 27, hour_utc),
                            np.array([[T_air]]), np.array([[230.0]]), long_tail=long_tail)
    ps, p1, dp1, z1 = (np.array([[100000.0]]), np.array([[97300.0]]),
                       np.array([[5400.0]]), np.array([[230.0]]))
    th1 = np.array([[T_air * (P0 / 97300.0) ** KAPPA]])
    Hs, Tgs, dth = [], [], []
    dt = 60.0
    for k in range(int(n_hours * 60)):
        rate, thg = s.step(k * dt, dt, np.array([[u]]), np.array([[0.0]]), th1, ps, p1, dp1, z1)
        Hs.append(float(s.last["H"][0, 0])); Tgs.append(float(s.Tg[0, 0])); dth.append(float(rate[0, 0]))
    return s, np.array(Hs), np.array(Tgs), np.array(dth)


# local midnight at 74 W is about 05 UTC
s, H, Tg, dth = column(5, 6)
Ta_s = 285.0 + 0.0065 * 230.0
report("night: the ground cools below the air; a downward flux of 10-80 W/m^2 cools the air",
       Tg[-1] < Ta_s - 1.0 and Tg[-1] > Ta_s - 15.0 and -80.0 < H[-1] < -10.0 and dth[-1] < 0
       and s.last["Ri"][0, 0] > 0,
       f"Tg {Tg[0]:.2f} -> {Tg[-1]:.2f} K (air at the ground {Ta_s:.2f} K), H {H[-1]:.1f} W/m^2, "
       f"Ri_b {s.last['Ri'][0, 0]:.3f}")

s0, H0, Tg0, _ = column(5, 6, long_tail=False)
report("night, Louis heat function without the long tail: the flux collapses (why long_tail is the default)",
       abs(H0[-1]) < 0.2 * abs(H[-1]),
       f"after 6 h: H {H0[-1]:.1f} W/m^2 and Tg {Tg0[-1]:.2f} K without the tail, against "
       f"H {H[-1]:.1f} W/m^2 and Tg {Tg[-1]:.2f} K with it")

# local noon at 74 W is about 17 UTC; start 3 h before
s, H, Tg, dth = column(14, 4)
report("day: the ground warms above the air and heats it, peak flux 50-300 W/m^2",
       Tg[-1] > Ta_s + 1.0 and 50.0 < H.max() < 300.0 and dth[-1] > 0
       and s.last["Ri"][0, 0] < 0,
       f"Tg {Tg[0]:.2f} -> {Tg[-1]:.2f} K, max H {H.max():.1f} W/m^2, Rn {s.last['Rn'][0, 0]:.0f}, "
       f"LE {s.last['LE'][0, 0]:.0f}")

s, H, Tg, dth = column(17, 3, land=False)
report("water: the surface temperature stays at its initial value",
       np.ptp(Tg) == 0.0, f"Tg range {np.ptp(Tg):.2e} K over 3 h at noon")

# 4. Drag with the Louis stability -------------------------------------------
gr = CGrid(8, 8, 12e3, 12e3, f0=1.0e-4, beta=0.0, edge_mode="replicate")
lev = SigmaLevels(20)
pi = np.full((8, 8), 101325.0 - lev.p_top)
p = lev.pressure(pi)
assert p.shape == (20, 8, 8), p.shape
theta = (288.0 - 55.0 * (1 - p / p.max())) / (p / P0) ** KAPPA
u = np.full((20, 8, 8), 8.0); v = np.zeros((20, 8, 8))
du0, _, i0 = surface_drag(u, v, theta, pi, lev, z0=0.1)
du_w, _, iw = surface_drag(u, v, theta, pi, lev, z0=0.1, theta_s=theta[-1] + 3.0, stability="louis")
du_c, _, ic = surface_drag(u, v, theta, pi, lev, z0=0.1, theta_s=theta[-1] - 6.0, stability="louis")
du_k, _, _ = surface_drag(u, v, theta, pi, lev, z0=0.1, theta_s=theta[-1] - 6.0, stability="cutoff")
r_w = float(du_w[-1].mean() / du0[-1].mean()); r_c = float(du_c[-1].mean() / du0[-1].mean())
report("drag: stronger over a warm ground, weaker but not zero over a cold one (cutoff gives zero)",
       r_w > 1.0 and 0.0 < r_c < 1.0 and float(np.abs(du_k[-1]).max()) == 0.0,
       f"|du| relative to neutral: warm ground {r_w:.2f}, cold ground {r_c:.3f} "
       f"(Ri_b {float(ic['Ri_bulk'].mean()):.2f}); cutoff form {float(np.abs(du_k[-1]).max()):.1f}")

# 4b. A land/sea roughness map gives each region its own scalar drag -----------
zmap = np.where(np.arange(8)[None, :] < 4, 1.0, 0.0002) * np.ones((8, 8))
du_m, _, _ = surface_drag(u, v, theta, pi, lev, z0=zmap)
du_l, _, _ = surface_drag(u, v, theta, pi, lev, z0=1.0)
du_s, _, _ = surface_drag(u, v, theta, pi, lev, z0=0.0002)
report("roughness map: land columns match z0 = 1.0, water columns match z0 = 0.0002",
       np.allclose(du_m[-1][:, :4], du_l[-1][:, :4]) and np.allclose(du_m[-1][:, 4:], du_s[-1][:, 4:])
       and abs(du_l[-1].mean()) > 3 * abs(du_s[-1].mean()),
       f"|du| land {abs(du_l[-1].mean())*3600:.2f}, water {abs(du_s[-1].mean())*3600:.2f} m/s per hour")

# 5. Coupled to the model: a night cools the lowest layer ----------------------
m = PrimitiveSigma(gr, lev)
m.pi[:] = pi
m.theta = theta.copy()
m.u = u.copy(); m.v = v.copy()
lat = np.full((8, 8), 42.0); lon = np.full((8, 8), -74.0)
p1 = lev.p_top + float(np.asarray(lev.sigma)[-1]) * pi
z1 = lowest_level_height(m.theta, m.pi, lev)
m.land_surface = ForceRestoreSurface(lat, lon, np.ones((8, 8), bool), datetime(2026, 9, 27, 5),
                                     m.theta[-1] * (p1 / P0) ** KAPPA, z1)
m.theta_surface = m.land_surface.Tg / ((lev.p_top + pi) / P0) ** KAPPA
th1_0 = m.theta[-1].mean()
dt = 30.0
for _ in range(240):
    m.step(dt)
d1 = float(m.theta[-1].mean() - th1_0)
ok = np.isfinite(m.u).all() and np.isfinite(m.theta).all() and d1 < -0.05
report("model coupling: 2 h of night cool the lowest layer and the drag sees a stable surface",
       ok and float(m._drag_info["Ri_bulk"].mean()) > 0,
       f"lowest-level theta change {d1:+.3f} K; Ri_b {float(m._drag_info['Ri_bulk'].mean()):.3f}; "
       f"Tg - T2 {float((m.land_surface.Tg - m.land_surface.T2).mean()):+.2f} K")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

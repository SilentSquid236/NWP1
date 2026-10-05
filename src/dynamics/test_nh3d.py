"""
Tests for the 3-D non-hydrostatic core (nh3d.py; CAM stage S2).

Run:  python test_nh3d.py      (about 10 s)
"""
import numpy as np

import nh3d
from grid import CGrid
from sigma import SigmaLevels

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


lev = SigmaLevels(20)

# 1. Rest over a 1500 m mountain on NWP1's levels ----------------------------
gr = CGrid(40, 30, 12000.0, 12000.0, f0=1e-4, beta=0.0, edge_mode="periodic")
X, Y = gr.Xc - gr.Lx / 2, gr.Yc - gr.Ly / 2
m = nh3d.NH3D(gr, lev, terrain=1500.0 * np.exp(-(X ** 2 + Y ** 2) / 6.0e4 ** 2))
M0 = m.total_mass()
for _ in range(60):
    m.step(60.0)
u, v, w = m.winds()
report("rest over a 1500 m mountain for 1 h (dt 60 s): |u|, |v|, |w| < 1e-9 m/s; mass unchanged",
       max(np.abs(u).max(), np.abs(v).max(), np.abs(w).max()) < 1e-9 and abs(m.total_mass() / M0 - 1) < 1e-13,
       f"max|u| {np.abs(u).max():.1e}, |v| {np.abs(v).max():.1e}, |w| {np.abs(w).max():.1e}")

# 2. Inertial oscillation on an f-plane ---------------------------------------
g2 = CGrid(8, 8, 12000.0, 12000.0, f0=1e-4, beta=0.0, edge_mode="periodic")
m = nh3d.NH3D(g2, lev)
m.U = 10.0 * m.h_to_u(m.mu)[None] * np.ones((lev.nz, 8, 8))
for _ in range(int(round(2 * np.pi / 1e-4 / 4 / 60.0))):
    m.step(60.0)
u, v, _ = m.winds()
ue, ve = 10 * np.cos(1e-4 * m.time), -10 * np.sin(1e-4 * m.time)
report("inertial oscillation: after a quarter period u and v match 10 cos ft, -10 sin ft within 0.01 m/s",
       abs(u.mean() - ue) < 0.01 and abs(v.mean() - ve) < 0.01,
       f"u {u.mean():.4f} (expect {ue:.4f}), v {v.mean():.4f} (expect {ve:.4f})")


# 3. The x and y directions are the same code ----------------------------------
def bubble(nx, ny, along):
    g = CGrid(nx, ny, 2000.0, 2000.0, f0=0.0, beta=0.0, edge_mode="periodic")
    mm = nh3d.NH3D(g, lev, K=50.0, theta_ref=lambda z: 300.0 + 0 * z)
    zc = 0.5 * (mm.phib[:-1] + mm.phib[1:]) / 9.80665
    c = (g.Xc - g.Lx / 2) if along == "x" else (g.Yc - g.Ly / 2)
    L = np.sqrt((c[None] / 12000.0) ** 2 + ((zc - 3000.0) / 2000.0) ** 2)
    dT = np.where(L <= 1, -10.0 * (np.cos(np.pi * L) + 1) / 2, 0.0)
    mm.Th = mm.mu[None] * (mm.thb + dT / (mm.pb / nh3d.P0) ** nh3d.KAPPA)
    M0 = mm.total_mass()
    for _ in range(60):
        mm.step(20.0)
    return mm, M0

a, Ma = bubble(64, 4, "x")
b, _ = bubble(4, 64, "y")
ta, tb = a.theta() - a.thb, b.theta() - b.thb
d = np.abs(ta[:, 0, :] - tb[:, :, 0]).max()
report("a cold bubble along x and the same bubble along y give transposed fields to 1e-10 K; mass conserved",
       d < 1e-10 and abs(a.total_mass() / Ma - 1) < 1e-13 and np.abs(a.winds()[0]).max() > 5.0,
       f"difference {d:.1e} K; max|u| {np.abs(a.winds()[0]).max():.1f} m/s; mass change {a.total_mass() / Ma - 1:.1e}")

# 4. The forecast adapter relaxes per unit time, not per step ------------------
# forecast.py's Davies weights were set with the hydrostatic core's ~17 s step;
# at 60 s the same per-step weights relax ~3.5x more weakly per hour. On the
# 28 Sep 06Z case this alone made the NH and hydrostatic runs differ by
# 0.38 m/s rms after 1 h (0.025 m/s once rescaled).
class _Hydro:
    def __init__(self, g):
        self.grid, self.lev, self.terrain = g, lev, np.zeros((g.ny, g.nx))
        self._sponge = np.zeros((lev.nz, 1, 1))
    def max_dt(self):
        return 20.0

class _Relax:
    def __init__(self, g):
        self.alpha2d = np.full((g.ny, g.nx), 0.1)

def _relaxed(dts):
    m = nh3d.NHModel(_Hydro(g2))
    shape = (lev.nz, 8, 8)
    m.set_state(np.zeros(shape), np.zeros(shape), np.full(shape, 300.0), np.full((8, 8), 8.0e4))
    for dt in dts:
        m._last_dt = dt
        m.relax_with(_Relax(g2), {"pi": np.full((8, 8), 8.1e4)})
    return float(m.pi.mean())

one_ref, two_half, one_long = _relaxed([20.0]), _relaxed([10.0, 10.0]), _relaxed([60.0])
expect_long = 8.0e4 + 1.0e3 * (1 - 0.9 ** 3)
report("the NH adapter's edge relaxation depends on elapsed time, not on the number of steps",
       abs(one_ref - 8.01e4) < 1e-6 and abs(two_half - one_ref) < 1e-6 and abs(one_long - expect_long) < 1e-6,
       f"one 20 s step {one_ref:.3f} Pa, two 10 s steps {two_half:.3f} Pa, one 60 s step {one_long:.3f} "
       f"(expect {expect_long:.3f})")

# 5. The September terrain test, with a state the reference profile does not match
# The first hydrostatic core diverged at +3 h on this test (an isothermal
# atmosphere at rest over a 2500 m mountain; docs/instability_growth.png).
# Test 1 starts from the core's own reference profile, where rest is exact by
# construction. Here the state is isothermal at 250 K and the reference is the
# default profile, so the terrain-following pressure-gradient error is live.
# Measured 2026-10-02: max|u| 0.15-0.22 m/s and max|w| 2-4 mm/s, no growth in 24 h.
class _HydroNoPhys(_Hydro):
    drag = mixing = convection = False
    sponge_levels = 0

g5 = CGrid(74, 4, 12000.0, 12000.0, f0=1e-4, beta=0.0, edge_mode="periodic")
X5 = g5.Xc - g5.Lx / 2
h5 = _HydroNoPhys(g5); h5.terrain = 2500.0 * np.exp(-(X5 / 150e3) ** 2)
m5 = nh3d.NHModel(h5)
mu5 = 101325.0 * np.exp(-9.80665 * h5.terrain / (287.04 * 250.0)) - lev.p_top
p5 = lev.p_top + np.asarray(lev.sigma)[:, None, None] * mu5[None]
z5 = np.zeros((lev.nz, 4, 74))
m5.set_state(z5, z5, 250.0 * (1.0e5 / p5) ** (287.04 / 1004.5), mu5)
M5 = m5.core.total_mass(); umax = []
for k in range(360):
    m5.core.step(60.0)
    if (k + 1) % 60 == 0:
        umax.append(float(np.abs(m5.core.winds()[0]).max()))
report("an isothermal atmosphere at rest over 2500 m (state != reference) stays below 0.5 m/s for 6 h without growth",
       max(umax) < 0.5 and umax[-1] < 2 * umax[0] and abs(m5.core.total_mass() / M5 - 1) < 1e-13,
       f"max|u| by hour {', '.join(f'{x:.3f}' for x in umax)} m/s")

# 6. Physics once per step: one evaluation per step instead of one per RK stage
def _calls(mode):
    h6 = _HydroNoPhys(g5); h6.terrain = h5.terrain
    m6 = nh3d.NHModel(h6, physics_every=mode)
    m6.set_state(z5, z5, 250.0 * (1.0e5 / p5) ** (287.04 / 1004.5), mu5)
    for _ in range(5):
        m6.step(60.0)
    return m6.physics_calls, m6

n_stage, a6 = _calls("stage")
n_step, b6 = _calls("step")
d6 = float(np.abs(a6.u - b6.u).max())
report("physics_every='step' evaluates the physics once per step (stage: three times)",
       n_stage == 15 and n_step == 5 and d6 < 1e-12 and b6.core.extra_tendency == b6._physics,
       f"calls: stage {n_stage}, step {n_step}; max|du| {d6:.1e} m/s with no physics active; hook restored")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

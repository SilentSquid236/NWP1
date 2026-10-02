"""
Tests for the horizontal advection schemes in PrimitiveSigma._horiz_adv
("centred2" and "upwind3"; P-67 test AL).

Run:  python test_advection.py
"""
import numpy as np

from grid import CGrid
from sigma import SigmaLevels
from primitive_sigma import PrimitiveSigma

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def model(n, scheme, dx=1.0e4):
    gr = CGrid(n, 4, dx, dx, f0=0.0, beta=0.0, edge_mode="periodic")
    m = PrimitiveSigma(gr, SigmaLevels(10), sponge_levels=2)
    m.advection = scheme
    return m, gr


# 1. Order of accuracy on a smooth periodic wave ------------------------------
def deriv_error(n, scheme):
    m, gr = model(n, scheme)
    x = np.arange(n) * gr.dx
    L = n * gr.dx
    a = np.broadcast_to(np.sin(2 * np.pi * x / L), (4, n)).copy()
    exact = (2 * np.pi / L) * np.cos(2 * np.pi * x / L)
    got = m._horiz_adv(a, np.ones_like(a), np.zeros_like(a))[0]
    return np.abs(got - exact).max() / np.abs(exact).max()

orders = {}
for s in ("centred2", "upwind3"):
    e1, e2 = deriv_error(32, s), deriv_error(64, s)
    orders[s] = np.log2(e1 / e2)
report("order of accuracy: centred2 about 2, upwind3 about 3",
       1.8 < orders["centred2"] < 2.2 and 2.7 < orders["upwind3"] < 3.3,
       f"centred2 {orders['centred2']:.2f}, upwind3 {orders['upwind3']:.2f}")

# 2. A uniform field is not changed; a wind of either sign works ---------------
m, gr = model(16, "upwind3")
a = np.full((4, 16), 3.0)
z1 = np.abs(m._horiz_adv(a, np.full_like(a, 5.0), np.full_like(a, -4.0))).max()
x = np.arange(16) * gr.dx
b = np.broadcast_to(np.sin(2 * np.pi * x / (16 * gr.dx)), (4, 16)).copy()
pos = m._horiz_adv(b, np.ones_like(b), np.zeros_like(b))[0]
neg = m._horiz_adv(b, -np.ones_like(b), np.zeros_like(b))[0]
report("uniform field gives zero; reversing the wind reverses the centred part",
       z1 < 1e-12 and np.allclose(pos + neg, 2 * np.abs(1.0) *
                                  (np.roll(b[0], -2) - 4 * np.roll(b[0], -1) + 6 * b[0]
                                   - 4 * np.roll(b[0], 1) + np.roll(b[0], 2)) / (12 * gr.dx)),
       f"max |tendency| of a uniform field {z1:.1e}")


# 3. One revolution of a narrow bump, RK3 at Courant 0.5 ----------------------
def revolve(scheme, n=100, width=3.0):
    m, gr = model(n, scheme)
    x = np.arange(n)
    a = np.broadcast_to(np.exp(-0.5 * ((x - n / 2) / width) ** 2), (4, n)).copy()
    u = np.ones_like(a) * 10.0
    dt = 0.5 * gr.dx / 10.0
    f = lambda q: -m._horiz_adv(q, u, np.zeros_like(q))
    for _ in range(int(round(n / 0.5))):
        q1 = a + dt / 3 * f(a)
        q2 = a + dt / 2 * f(q1)
        a = a + dt * f(q2)
    return a[0]

c2, u3 = revolve("centred2"), revolve("upwind3")
report("one revolution of a 3-cell bump: upwind3 has far smaller spurious extrema",
       u3.min() > 0.5 * c2.min() and u3.max() <= 1.0 + 1e-9 and c2.min() < -0.02,
       f"centred2 min {c2.min():+.3f} max {c2.max():.3f}; upwind3 min {u3.min():+.3f} max {u3.max():.3f}")

# 4. The 2-dx wave: centred2 does not see it, upwind3 damps it ----------------
for s in ("centred2", "upwind3"):
    m, gr = model(16, s)
    w = np.broadcast_to(np.array([1.0, -1.0] * 8), (4, 16)).copy()
    t = m._horiz_adv(w, np.ones_like(w) * 10.0, np.zeros_like(w))[0]
    orders[s + "_2dx"] = float(np.sign(t[0]) * np.abs(t).max() * np.sign(w[0, 0]))
report("2-dx wave: centred2 tendency 0; upwind3 tendency opposes it (damping)",
       abs(orders["centred2_2dx"]) < 1e-12 and orders["upwind3_2dx"] > 0,
       f"centred2 {orders['centred2_2dx']:.2e}; upwind3 -tendency/amplitude sign "
       f"{'damping' if orders['upwind3_2dx'] > 0 else 'growing'} ({orders['upwind3_2dx']:.2e} 1/s)")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

"""
Tests for the 2-D non-hydrostatic core (nh2d.py; CAM roadmap stage S1a).

Run:  python test_nh2d.py      (about 15 s)
"""
import numpy as np

import nh2d
from nh2d import NH2D, _thomas

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# 1. The tridiagonal solver ----------------------------------------------------
rng = np.random.default_rng(1)
n, m = 9, 4
a, c = rng.normal(size=(n, m)), rng.normal(size=(n, m))
b = 4.0 + np.abs(rng.normal(size=(n, m)))
d = rng.normal(size=(n, m))
x = _thomas(a, b, c, d)
err = 0.0
for j in range(m):
    M = np.diag(b[:, j]) + np.diag(a[1:, j], -1) + np.diag(c[:-1, j], 1)
    err = max(err, np.abs(np.linalg.solve(M, d[:, j]) - x[:, j]).max())
report("tridiagonal solve matches a dense solve", err < 1e-12, f"max difference {err:.1e}")

# 2. A resolved sound pulse travels at the speed of sound ---------------------
mm = NH2D(800, 6, 100.0, 100.0, lambda z: 300.0 + 0 * z, ns=6, div_damp=0.0)
xc = (np.arange(800) + 0.5) * 100.0
mm.pp[:] = 1e-5 * np.exp(-((xc - 40000.0) / 2000.0) ** 2)[None, :]
pos = []
for _ in range(8):
    mm.run(10.0, 0.6)
    prof = mm.pp[2]
    i = np.argmax(prof[400:]) + 400
    y0, y1, y2 = prof[i - 1:i + 2]
    pos.append((mm.time, xc[i] + 50.0 * (y0 - y2) / (y0 - 2 * y1 + y2)))
t, xp = np.array(pos).T
c_meas = np.polyfit(t, xp, 1)[0]
c_th = np.sqrt(nh2d.CP / nh2d.CV * nh2d.RD * 300.0 * mm.pi_c[2])
report("sound speed of a resolved pulse within 0.5 %", abs(c_meas / c_th - 1) < 0.005,
       f"measured {c_meas:.1f} m/s, theory {c_th:.1f} m/s")


# 3. The Straka et al. (1993) density current at 200 m -------------------------
def straka(dx, dt, ns=6):
    L, H = 51200.0, 6400.0
    s = NH2D(int(L / dx), int(H / dx), dx, dx, lambda z: 300.0 + 0 * z, K=75.0, ns=ns)
    xx = (np.arange(s.nx) + 0.5) * dx - L / 2
    X, Z = np.meshgrid(xx, s.z_c)
    Ld = np.sqrt((X / 4000.0) ** 2 + ((Z - 3000.0) / 2000.0) ** 2)
    s.th = np.where(Ld <= 1.0, -15.0 * (np.cos(np.pi * Ld) + 1) / 2, 0.0) / s.pi_c[:, None]
    s.run(900.0, dt)
    idx = np.where((xx > 0) & (s.th[0] < -1.0))[0]
    return s, xx[idx.max()] / 1000.0

s, front = straka(200.0, 2.0)
sym = np.abs(s.th - s.th[:, ::-1]).max()
report("density current at 200 m: front 15.0-16.0 km at 900 s (reference about 15.5), symmetric",
       15.0 <= front <= 16.0 and sym < 1e-9 and np.isfinite(s.u).all(),
       f"front {front:.2f} km, theta' min {s.th.min():.2f} K, max|u| {np.abs(s.u).max():.1f} m/s, symmetry {sym:.1e}")
s2, front2 = straka(200.0, 2.0, ns=12)
report("halving the acoustic step changes the front by < 0.1 km",
       abs(front2 - front) < 0.1, f"ns 6: {front:.2f} km; ns 12: {front2:.2f} km")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

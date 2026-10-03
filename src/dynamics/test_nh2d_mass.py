"""
Tests for the mass-coordinate 2-D non-hydrostatic core (nh2d_mass.py; CAM stage S1b).

Run:  python test_nh2d_mass.py      (about 1 minute)
"""
import numpy as np

import nh2d_mass
from nh2d_mass import NH2DMass

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


T0 = 250.0
th_iso = lambda z: T0 * np.exp(9.80665 * z / (1004.5 * T0))

# 1. An atmosphere at rest over a 1000 m mountain stays at rest ---------------
nx, dx = 120, 2000.0
x = (np.arange(nx) + 0.5) * dx - nx * dx / 2
m = NH2DMass(nx, 40, dx, 20000.0, th_iso, terrain=1000.0 * 1e8 / (x ** 2 + 1e8), ns=6)
mass0 = m.total_mass()
for _ in range(180):
    m.step(10.0)
umax, wmax = np.abs(m.u_centres()).max(), np.abs(m.w_field()).max()
report("rest over a 1000 m bell for 30 min: |u|, |w| < 1e-9 m/s; mass unchanged",
       umax < 1e-9 and wmax < 1e-9 and abs(m.total_mass() / mass0 - 1) < 1e-13,
       f"max|u| {umax:.1e}, max|w| {wmax:.1e}, mass change {m.total_mass() / mass0 - 1:.1e}")

# 2. Linear hydrostatic mountain wave, lower troposphere ----------------------
U0, h0, a = 10.0, 10.0, 10000.0
nx, dx = 300, 2000.0
x = (np.arange(nx) + 0.5) * dx - nx * dx / 2
mw = NH2DMass(nx, 60, dx, 30000.0, th_iso, terrain=h0 * a ** 2 / (x ** 2 + a ** 2), ns=6,
              sponge_depth=10000.0, sponge_rate=1 / 300.0, u0=U0)
mass0 = mw.total_mass()
for _ in range(1800):
    mw.step(10.0)
N = 9.80665 / np.sqrt(1004.5 * T0)
H = nh2d_mass.RD * T0 / 9.80665
l = np.sqrt(N ** 2 / U0 ** 2 - 1 / (4 * H ** 2))
z = mw.phi / 9.80665
X = np.broadcast_to(x, z.shape)
delta = lambda xx, zz: h0 * a * np.exp(zz / (2 * H)) * (a * np.cos(l * zz) - xx * np.sin(l * zz)) / (xx ** 2 + a ** 2)
w_an = U0 * (delta(X + 1.0, z) - delta(X - 1.0, z)) / 2.0
sel = (np.abs(X) <= 60000.0) & (z >= 1000.0) & (z <= 3000.0)
cc = np.corrcoef(mw.w_field()[sel], w_an[sel])[0, 1]
slope = np.polyfit(w_an[sel], mw.w_field()[sel], 1)[0]
report("mountain wave against the analytic solution, 1-3 km, 5 h: correlation >= 0.95, slope 0.85-1.10; mass conserved",
       cc >= 0.95 and 0.85 <= slope <= 1.10 and abs(mw.total_mass() / mass0 - 1) < 1e-12,
       f"correlation {cc:.3f}, slope {slope:.3f}, mass change {mw.total_mass() / mass0 - 1:.1e}")

# 3. Straka et al. (1993) density current at 200 m ------------------------------
L = 51200.0
ms = NH2DMass(int(L / 200.0), 32, 200.0, 6400.0, lambda zz: 300.0 + 0 * zz, K=75.0, ns=6)
xs = (np.arange(ms.nx) + 0.5) * 200.0 - L / 2
zc = 0.5 * (ms.phib[1:] + ms.phib[:-1]) / 9.80665
Ld = np.sqrt((np.broadcast_to(xs, zc.shape) / 4000.0) ** 2 + ((zc - 3000.0) / 2000.0) ** 2)
dT = np.where(Ld <= 1.0, -15.0 * (np.cos(np.pi * Ld) + 1) / 2, 0.0)
ms.Th = ms.mu[None, :] * (ms.thb + dT / (ms.pb / nh2d_mass.P0) ** nh2d_mass.KAPPA)
for _ in range(450):
    ms.step(2.0)
thp = ms.Th / ms.mu[None, :] - ms.thb
front = xs[np.where((xs > 0) & (thp[0] < -1.0))[0].max()] / 1000.0
sym = np.abs(thp - thp[:, ::-1]).max()
report("density current at 200 m: front 15.0-15.8 km, min theta' above -12 K (never below the initial -15), symmetric",
       15.0 <= front <= 15.8 and thp.min() > -12.0 and sym < 1e-9,
       f"front {front:.2f} km, min theta' {thp.min():.2f} K, symmetry {sym:.1e}")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

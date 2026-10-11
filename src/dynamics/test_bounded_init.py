"""
Tests for the bounded-domain initialisation (P-76, test E4): the mirrored
lowpass and the exact divergence removal.

Run:  python test_bounded_init.py      (a few seconds)
"""
import numpy as np

import initialization
import subgrid
from grid import CGrid

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


ny, nx = 97, 110
gr = CGrid(nx, ny, 12000.0, 12000.0, f0=1e-4, beta=0.0, edge_mode="replicate")
J, I = np.indices((ny, nx))
D = np.minimum.reduce([J, I, ny - 1 - J, nx - 1 - I])


def rms(a, m):
    return float(np.sqrt(np.mean(a[m] ** 2)))


# 1. A smooth non-periodic field passes the mirrored lowpass unchanged --------
th = 300 + 10 * (J / ny) + 2 * np.sin(2 * np.pi * I / 40) * np.cos(2 * np.pi * J / 55) + 3 * (I / nx) ** 2
ref = th.mean()
per = ref + initialization.spectral_lowpass(th - ref, gr)
mir = ref + initialization.spectral_lowpass(th - ref, gr, mirror=True)
e_per, e_mir = rms(per - th, D < 3), rms(mir - th, D < 3)
report("mirrored lowpass leaves a smooth 10 K gradient alone at the edges (outer 3 cells < 0.1 K); "
       "the periodic one does not",
       e_mir < 0.1 and e_per > 1.0 and rms(mir - th, D >= 18) < 1e-3,
       f"outer-3-cell rms change: periodic {e_per:.2f} K, mirrored {e_mir:.3f} K")

# 2. Grid-scale noise is still removed by the mirrored filter -------------------
rng = np.random.default_rng(1)
noise = rng.standard_normal((ny, nx))
left = initialization.spectral_lowpass(noise, gr, mirror=True)
report("mirrored lowpass still removes most 2-4 dx noise (variance kept < 30 %)",
       left.var() / noise.var() < 0.3, f"variance kept {left.var() / noise.var():.3f}")

# 3. Exact divergence removal: zero model divergence in EVERY cell -----------
psi = 2e6 * np.sin(np.pi * (J + 3) / 70.0) * np.cos(np.pi * I / 90.0) + 1e6 * (J / ny)
phi = 4e5 * np.cos(np.pi * I / 50.0) * np.sin(np.pi * J / 45.0)
u0 = -(np.roll(psi, -1, 0) - psi) / gr.dy + gr.dx_backward(phi)
v0 = (np.roll(psi, -1, 1) - psi) / gr.dx + gr.dy_backward(phi)
u0[-1] = u0[-2]
v0[:, -1] = v0[:, -2]
u3 = np.stack([u0, 0.5 * u0])
v3 = np.stack([v0, -0.3 * v0])


def div(u, v):
    return gr.dx_forward(u) + gr.dy_forward(v)


uh, vh, _ = subgrid.balance_initial_state(u3, v3, gr, verbose=False)
ue, ve, info = subgrid.balance_initial_state(u3, v3, gr, verbose=False, method="exact")
dh = np.abs(div(uh, vh))
de = np.abs(div(ue, ve))
report("exact removal: |div| < 1e-14 1/s in every cell including the edges; hybrid leaves > 1e-6 "
       "in the outer cells",
       de.max() < 1e-14 and dh[:, D < 3].max() > 1e-6,
       f"max|div| exact {de.max():.1e}, hybrid outer 3 cells {dh[:, D < 3].max():.1e}, "
       f"hybrid interior {dh[:, D >= 18].max():.1e}")

# 4. The western and southern boundary-face winds are not changed -------------
report("exact removal keeps the western and southern boundary-face winds",
       np.abs(ue[:, :, 0] - u3[:, :, 0]).max() == 0 and np.abs(ve[:, 0, :] - v3[:, 0, :]).max() == 0,
       f"max change u[:, :, 0] {np.abs(ue[:, :, 0] - u3[:, :, 0]).max():.1e}, "
       f"v[:, 0, :] {np.abs(ve[:, 0, :] - v3[:, 0, :]).max():.1e}")

# 5. A non-divergent field is left alone --------------------------------------
u5, v5 = ue.copy(), ve.copy()
u6, v6, _ = subgrid.balance_initial_state(u5, v5, gr, verbose=False, method="exact")
report("exact removal is a projection: applied twice it changes nothing (< 1e-10 m/s)",
       max(np.abs(u6 - u5).max(), np.abs(v6 - v5).max()) < 1e-10,
       f"max change {max(np.abs(u6 - u5).max(), np.abs(v6 - v5).max()):.1e} m/s")

# 6. The 3 km grid (389 x 439): exact there too, and fast ---------------------
# E5 (2026-10-10): the first, SuperLU-based solver was exact at 12 km on the
# desktop but returned residuals up to 10^9 times the divergence at 3 km on the
# server's SciPy 1.9.1; every 3 km run started with 10^7 m/s winds.
import time
g3 = CGrid(439, 389, 3000.0, 3000.0, f0=1e-4, beta=0.0, edge_mode="replicate")
J3, I3 = np.indices((389, 439))
psi3 = 2e6 * np.sin(np.pi * (J3 + 3) / 272.0) * np.cos(np.pi * I3 / 351.0)
phi3 = 4e5 * np.cos(np.pi * I3 / 197.0) * np.sin(np.pi * J3 / 179.0)
u3k = -(np.roll(psi3, -1, 0) - psi3) / 3000.0 + g3.dx_backward(phi3)
v3k = (np.roll(psi3, -1, 1) - psi3) / 3000.0 + g3.dy_backward(phi3)
d0 = np.abs(g3.dx_forward(u3k) + g3.dy_forward(v3k)).max()
t0 = time.time()
u3e, v3e = subgrid.remove_divergence_exact(np.stack([u3k] * 3), np.stack([v3k] * 3), g3)
t3 = time.time() - t0
d3 = np.abs(g3.dx_forward(u3e) + g3.dy_forward(v3e)).max()
report("exact removal on the 3 km grid (389 x 439, 3 levels): residual < 1e-9 of the divergence, max|u| bounded, < 10 s",
       d3 < 1e-9 * d0 and np.abs(u3e).max() < 10 * np.abs(u3k).max() and t3 < 10.0,
       f"max|div| {d0:.2e} -> {d3:.2e} 1/s; max|u| {np.abs(u3k).max():.0f} -> {np.abs(u3e).max():.0f} m/s; {t3:.2f} s")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

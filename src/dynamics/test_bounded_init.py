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

# J3 (2026-10-15/17, P-83): a level with NET convergence. The exact solve keeps
# the western and southern faces, so the whole net inflow leaves through the
# northern and eastern ones and the mean wind shifts; exact-sym takes it evenly
# from all four edges and keeps the mean wind.
D0 = -1.5e-5                                    # 1/s, the size seen at 250-300 hPa
xs = (I - 0.5 * (nx - 1)) * 12000.0; ys = (J - 0.5 * (ny - 1)) * 12000.0
psi = 2e6 * np.sin(np.pi * (J + 3) / 70.0) * np.cos(np.pi * I / 90.0)
u_rot = -(np.roll(psi, -1, 0) - psi) / 12000.0; v_rot = (np.roll(psi, -1, 1) - psi) / 12000.0
u_rot[-1] = u_rot[-2]; v_rot[:, -1] = v_rot[:, -2]
u_in = 20.0 + 0.5 * D0 * xs + u_rot; v_in = 10.0 + 0.5 * D0 * ys + v_rot
mean_div = float(np.mean(gr.dx_forward(u_in) + gr.dy_forward(v_in)))
ue_, ve_ = subgrid.remove_divergence_exact(u_in, v_in, gr)
us_, vs_ = subgrid.remove_divergence_exact(u_in, v_in, gr, symmetric=True)
ds = np.abs(gr.dx_forward(us_) + gr.dy_forward(vs_)).max()
# Reference: the same flow WITHOUT the net convergence, through the exact solve
# (the discrete rotational part has edge divergence of its own, handled the same way).
ur_, vr_ = subgrid.remove_divergence_exact(20.0 + u_rot, 10.0 + v_rot, gr)
shift_e = (float((ue_ - ur_).mean()), float((ve_ - vr_).mean()))
shift_s = (float((us_ - ur_).mean()), float((vs_ - vr_).mean()))
err_s = float(np.sqrt(np.mean((us_ - ur_) ** 2 + (vs_ - vr_) ** 2)))
err_e = float(np.sqrt(np.mean((ue_ - ur_) ** 2 + (ve_ - vr_) ** 2)))
report("exact-sym removes a net convergence evenly: divergence ~0, and the result equals the "
       "flow without the convergence (< 0.05 m/s); exact shifts the mean wind",
       ds < 1e-12 and err_s < 0.05 and max(map(abs, shift_e)) > 2.0,
       f"mean div {mean_div:.2e} 1/s; against the no-convergence flow: exact mean (du, dv) "
       f"({shift_e[0]:+.2f}, {shift_e[1]:+.2f}), rms {err_e:.2f} m/s; exact-sym "
       f"({shift_s[0]:+.3f}, {shift_s[1]:+.3f}), rms {err_s:.3f} m/s; exact-sym max|div| {ds:.1e}")

# J4 (2026-10-17, P-83): the minimum-norm solve with no face held fixed.
rng4 = np.random.default_rng(4)
u4 = 20.0 + 8.0 * np.sin(2 * np.pi * J / 61.0) + 5.0 * np.cos(2 * np.pi * I / 47.0) + rng4.standard_normal((ny, nx))
v4 = 5.0 + 6.0 * np.cos(2 * np.pi * (I + J) / 53.0) + rng4.standard_normal((ny, nx))
uf4, vf4 = subgrid.remove_divergence_free(u4, v4, gr)
ue4, ve4 = subgrid.remove_divergence_exact(u4, v4, gr)
df4 = np.abs(gr.dx_forward(uf4) + gr.dy_forward(vf4)).max()
uf5, vf5 = subgrid.remove_divergence_free(uf4, vf4, gr)
n_free = float(np.sqrt(np.sum((uf4 - u4) ** 2 + (vf4 - v4) ** 2)))
n_exact = float(np.sqrt(np.sum((ue4 - u4) ** 2 + (ve4 - v4) ** 2)))
report("exact-free: |div| < 1e-13 in every cell, a projection (twice = once), and its correction is no "
       "larger than exact's (minimum norm)",
       df4 < 1e-13 and max(np.abs(uf5 - uf4).max(), np.abs(vf5 - vf4).max()) < 1e-10 and n_free <= n_exact * (1 + 1e-12),
       f"max|div| {df4:.1e} 1/s; second pass changes {max(np.abs(uf5 - uf4).max(), np.abs(vf5 - vf4).max()):.1e} m/s; "
       f"correction norm free {n_free:.1f} vs exact {n_exact:.1f} m/s")
uf6, vf6 = subgrid.remove_divergence_free(u_in, v_in, gr)
shift_f = (float((uf6 - u_in).mean()), float((vf6 - v_in).mean()))
moved = max(np.abs(uf6[:, 0] - u_in[:, 0]).max(), np.abs(vf6[0, :] - v_in[0, :]).max())
report("exact-free on a net convergence: the mean wind is kept (< 0.5 m/s, exact shifts it by several), "
       "and the western and southern faces are allowed to change",
       max(map(abs, shift_f)) < 0.5 and moved > 0.1 and max(map(abs, shift_e)) > 2.0,
       f"mean (du, dv): exact-free ({shift_f[0]:+.2f}, {shift_f[1]:+.2f}), exact ({shift_e[0]:+.2f}, {shift_e[1]:+.2f}) m/s; "
       f"largest W/S face change {moved:.2f} m/s")
t0 = time.time()
u3f, v3f = subgrid.remove_divergence_free(np.stack([u3k] * 3), np.stack([v3k] * 3), g3)
t3f = time.time() - t0
d3f = np.abs(g3.dx_forward(u3f) + g3.dy_forward(v3f)).max()
report("exact-free on the 3 km grid (389 x 439, 3 levels): residual < 1e-9 of the divergence, < 10 s",
       d3f < 1e-9 * d0 and t3f < 10.0, f"max|div| {d0:.2e} -> {d3f:.2e} 1/s; {t3f:.2f} s")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

"""
Tests for nest.py (12 km forecast -> 3 km boundary frames; CAM stage S5).

Run:  python test_nest.py
"""
import numpy as np

import nest

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


ny, nx, NY, NX = 10, 12, 40, 48           # a 4x refinement, as 12 -> 3 km


def linear(ny_, nx_, sy, sx):
    """3 + 2 X + 5 Y at the points of a staggering, X, Y in domain fractions."""
    y = (np.arange(ny_) + (0.5 if sy == "c" else 0.0)) / ny_
    x = (np.arange(nx_) + (0.5 if sx == "c" else 0.0)) / nx_
    return 3.0 + 2.0 * x[None, :] + 5.0 * y[:, None]


# 1. Bilinear regridding is exact for a linear field, for every staggering,
#    wherever the destination point lies inside the source points.
errs = {}
for sy, sx in (("c", "c"), ("c", "f"), ("f", "c")):
    a = nest.regrid(linear(ny, nx, sy, sx), (NY, NX), sy, sx)
    want = linear(NY, NX, sy, sx)
    fy = nest.fractional_index(ny, NY, sy); fx = nest.fractional_index(nx, NX, sx)
    inside = ((fy >= 0) & (fy <= ny - 1))[:, None] & ((fx >= 0) & (fx <= nx - 1))[None, :]
    errs[sy + sx] = float(np.abs(a - want)[inside].max())
report("regridding is exact for linear fields (centres, u faces, v faces)",
       max(errs.values()) < 1e-12, ", ".join(f"{k} {v:.1e}" for k, v in errs.items()))

# 2. Same grid, same terrain, same levels: the frame is the input
sig = (np.linspace(0, 1, 11)[:-1] + np.linspace(0, 1, 11)[1:]) / 2
rng = np.random.default_rng(0)
u = rng.normal(0, 5, (10, ny, nx)); v = rng.normal(0, 5, (10, ny, nx))
th = 290 + 40 * (1 - sig)[:, None, None] + rng.normal(0, 1, (10, ny, nx))
pi = 80000 + rng.normal(0, 500, (ny, nx))
h = np.abs(rng.normal(300, 200, (ny, nx)))
uu, vv, tt, pp = nest.frame_from_coarse(u, v, th, pi, sig, 20000.0, h, h, sig)
d = max(np.abs(uu - u).max(), np.abs(vv - v).max(), np.abs(tt - th).max(), np.abs(pp - pi).max())
report("on the same grid, terrain and levels the frame equals the input", d < 1e-9, f"max difference {d:.1e}")

# 3. Hydrostatic: an isothermal atmosphere at rest over flat 12 km ground, a
#    1000 m hill on the 3 km grid
T = 260.0
ps0 = 101325.0
sig_s = sig
p_s = 20000.0 + sig_s[:, None, None] * (ps0 - 20000.0) * np.ones((10, ny, nx))
th_s = T * (nest.P0 / p_s) ** nest.KAPPA
Y, X = np.meshgrid((np.arange(NY) + 0.5) / NY, (np.arange(NX) + 0.5) / NX, indexing="ij")
hill = 1000.0 * np.exp(-((X - 0.5) ** 2 + (Y - 0.5) ** 2) / 0.02)
sig_d = (np.linspace(0, 1, 21)[:-1] + np.linspace(0, 1, 21)[1:]) / 2
uu, vv, tt, pp = nest.frame_from_coarse(0 * p_s, 0 * p_s, th_s, np.full((ny, nx), ps0 - 20000.0),
                                        sig_s, 20000.0, np.zeros((ny, nx)), hill, sig_d)
ps_true = ps0 * np.exp(-nest.G * hill / (nest.RD * T))
rel = float(np.abs((pp + 20000.0) / ps_true - 1).max())
p_d = 20000.0 + sig_d[:, None, None] * pp[None]
th_true = T * (nest.P0 / p_d) ** nest.KAPPA
inside = p_d <= p_s[-1, 0, 0]                    # above the lowest coarse level
dth = float(np.abs(tt - th_true)[inside].max())
report("over a 1000 m hill the surface pressure is hydrostatic (0.3 %) and theta is "
       "isothermal (0.05 K) above the lowest coarse level",
       rel < 3e-3 and dth < 0.05 and np.abs(uu).max() == 0,
       f"surface pressure max relative error {rel:.1e}; theta max error {dth:.2f} K")

# 4. Below the lowest coarse level the lowest value is kept
pdeep = np.array([[[90000.0]], [[110000.0]]])
q = nest.remap_column(np.array([[[1.0]], [[2.0]]]), np.array([[[50000.0]], [[100000.0]]]), pdeep)
report("values beyond the end levels are the end values (no extrapolation)",
       abs(q[1, 0, 0] - 2.0) < 1e-12 and 1.0 < q[0, 0, 0] < 2.0, f"{q.ravel()}")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

"""
The compiled kernels of the NH core against the NumPy reference (CAM stage S4a).

Run:  python test_nh3d_c.py

Needs a C compiler (gcc). Without one every test is reported as SKIPPED and
the file exits 0: the NumPy core is unaffected, and the server, where the CAM
runs, has gcc.
"""
import numpy as np

import cnh
import nh3d
from grid import CGrid
from sigma import SigmaLevels

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def case(backend, edge, nx=48, ny=40, nz=20, steps=6, ns=6, threads=0, wind=10.0):
    g = CGrid(nx, ny, 3000.0, 3000.0, f0=1e-4, beta=0.0, edge_mode=edge)
    X, Y = g.Xc - g.Lx / 2, g.Yc - g.Ly / 2
    m = nh3d.NH3D(g, SigmaLevels(nz), terrain=1200.0 * np.exp(-(X ** 2 + Y ** 2) / 30e3 ** 2),
                  ns=ns, backend=backend, threads=threads)
    th = m.Th / m.mu[None]
    th = th + 3.0 * np.exp(-((X + 25e3) ** 2 + (Y - 10e3) ** 2) / 15e3 ** 2)[None]
    m.Th = th * m.mu[None]
    m.U = m.U + wind * m.h_to_u(m.mu)[None]
    m.V = m.V + 0.5 * wind * m.h_to_v(m.mu)[None]
    for _ in range(steps):
        m.step(20.0)
    return m


def rel_diff(a, b):
    out = {}
    for k in ("mu", "U", "V", "W", "Th", "phi"):
        x, y = getattr(a, k), getattr(b, k)
        out[k] = float(np.abs(x - y).max() / max(np.abs(x).max(), 1e-30))
    return out


if not cnh.available():
    try:
        cnh.load("f64")
    except cnh.Unavailable as e:
        why = str(e).splitlines()[0]
    print(f"  [SKIPPED] no compiled kernels on this machine: {why}")
    raise SystemExit(0)

# 0. The tendencies alone, from the same disturbed state, both edge modes
for edge in ("periodic", "replicate"):
    a = case("numpy", edge, steps=2)
    b = case("c", edge, steps=0)
    for k in ("mu", "U", "V", "W", "Th", "phi"):
        setattr(b, k, getattr(a, k).copy())
    Fa = a.tendencies(a.mu, a.U, a.V, a.W, a.Th, a.phi)
    Fb = b.tendencies(b.mu, b.U, b.V, b.W, b.Th, b.phi)
    names = ("dmu", "FU", "FV", "FW", "FTh", "Fphi")
    d = {n: float(np.abs(x - y).max() / max(np.abs(x).max(), 1e-30)) for n, x, y in zip(names, Fa, Fb)}
    report(f"C tendencies match NumPy, {edge} edges (max relative difference < 1e-11)",
           max(d.values()) < 1e-11, ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 1. Periodic edges, bubble over terrain in a mean wind, 6 steps (36 substeps)
a, b = case("numpy", "periodic"), case("c", "periodic")
d = rel_diff(a, b)
report("C step (tendencies, set-up, substeps) matches NumPy over 6 steps, periodic edges (< 1e-11)",
       max(d.values()) < 1e-11 and np.abs(b.W).max() > 1e-3,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; max|w| {np.abs(b.W).max() / b.mu.max():.2e}")

# 2. Replicate (limited-area) edges, more substeps
a, b = case("numpy", "replicate", ns=12), case("c", "replicate", ns=12)
d = rel_diff(a, b)
report("C step matches NumPy over 6 steps, replicate edges, ns = 12",
       max(d.values()) < 1e-11,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 3. Thread count does not change the answer (columns are independent)
a, b = case("c", "replicate", threads=1), case("c", "replicate", threads=4)
d = rel_diff(a, b)
report("C results are bit-identical on 1 and 4 threads",
       max(d.values()) == 0.0, ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 4. A wrong dtype is refused, not silently converted
try:
    k = cnh.load("f64")
    z = np.zeros((2, 3, 4), dtype=np.float32)
    k.ac_pd(z, z, z, np.zeros((3, 3, 4), np.float32), z, 0.1, z)
    ok, why = False, "float32 arrays were accepted by the float64 kernel"
except TypeError as e:
    ok, why = True, f"refused: {str(e)[:70]}"
report("a float32 array passed to a float64 kernel is refused", ok, why)

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

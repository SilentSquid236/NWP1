"""
Reference integration for the backend port (numpy vs torch).

A realistic small case: 44 x 40 cells at 12 km, 20 sigma levels, an 800 m
ridge, a baroclinic jet, stable stratification and seeded noise, with
mixing, drag and convective adjustment on (the production defaults).
Integrated for N steps at a fixed dt; writes u, v, theta, pi at the end.

    python src/dynamics/backend_reference.py OUT.npz [--steps 300] [--backend numpy|torch]
"""
import argparse
import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

from grid import CGrid                       # noqa: E402
from sigma import SigmaLevels, P0, KAPPA     # noqa: E402
from primitive_sigma import PrimitiveSigma   # noqa: E402


def build(nx=44, ny=40):
    gr = CGrid(nx, ny, 12e3, 12e3, f0=1.0e-4, beta=1.6e-11, edge_mode="replicate")
    lev = SigmaLevels(20, p_top=20000.0)
    y, x = np.meshgrid(np.arange(ny) * 12e3, np.arange(nx) * 12e3, indexing="ij")
    terrain = 800.0 * np.exp(-(((x - 0.4 * nx * 12e3) / 60e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=terrain)
    ps = 101000.0 * np.exp(-9.81 * terrain / (287.0 * 285.0))
    m.pi = ps - lev.p_top
    p = lev.pressure(m.pi)
    T = 288.0 * (p / 101325.0) ** 0.19 - 6.0 * (y / y.max() - 0.5)[None]
    m.theta = T * (P0 / p) ** KAPPA
    sig = lev.sigma.reshape(-1, 1, 1)
    m.u = 30.0 * np.exp(-((sig - 0.3) / 0.25) ** 2) * np.ones((1, ny, nx))
    m.v = np.zeros_like(m.u)
    rng = np.random.default_rng(7)
    m.u = m.u + 0.3 * rng.standard_normal(m.u.shape)
    m.v = m.v + 0.3 * rng.standard_normal(m.v.shape)
    m.theta = m.theta + 0.05 * rng.standard_normal(m.theta.shape)
    return m


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("out")
    ap.add_argument("--steps", type=int, default=300)
    ap.add_argument("--dt", type=float, default=20.0)
    ap.add_argument("--backend", default="numpy")
    ap.add_argument("--threads", type=int, default=None)
    ap.add_argument("--nx", type=int, default=44)
    ap.add_argument("--ny", type=int, default=40)
    a = ap.parse_args()
    m = build(a.nx, a.ny)
    if a.backend != "numpy":
        m.to_backend(a.backend, threads=a.threads)
    import time
    t0 = time.time()
    for _ in range(a.steps):
        m.step(a.dt)
    wall = time.time() - t0
    out = {k: np.asarray(m.as_numpy(getattr(m, k))) for k in ("u", "v", "theta", "pi")}
    np.savez(a.out, wall=wall, steps=a.steps, **out)
    print(f"{a.backend}: {a.steps} steps in {wall:.2f} s, max|u| {np.abs(out['u']).max():.3f}")


if __name__ == "__main__":
    main()

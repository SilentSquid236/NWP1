
"""Time PrimitiveSigma.tendencies on 12 km and 3 km sized arrays at several torch thread counts (CAM stage S3 groundwork)."""
import sys, time, numpy as np
sys.path[:0] = ["src", "src/dynamics"]
import torch
from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from backend import set_threads
def build(nx, ny, nz, dx):
    gr = CGrid(nx, ny, dx, dx, f0=1e-4, beta=0.0, edge_mode="replicate")
    lev = SigmaLevels(nz)
    m = PrimitiveSigma(gr, lev)
    m.pi[:] = 101325.0 - lev.p_top
    p = lev.pressure(m.pi)
    y = np.linspace(0, 1, ny)[None, :, None]
    m.theta = (288.0 - 55.0 * (1 - p / p.max()) - 2.0 * np.cos(2 * np.pi * y)) / (p / P0) ** KAPPA
    m.u = 10.0 + np.random.default_rng(0).normal(0, 0.5, p.shape); m.v = np.zeros(p.shape)
    return m
for (name, nx, ny, nz, dx) in (("12 km, 20 levels", 110, 97, 20, 12e3), ("3 km, 40 levels", 440, 388, 40, 3e3)):
    for n in (8, 16, 32, 64):
        set_threads(n)
        m = build(nx, ny, nz, dx); m.to_backend("torch", threads=n); m.advection = "upwind3"
        args = (m.u, m.v, m.theta, m.pi)
        m.tendencies(*args)
        t0 = time.perf_counter(); k = 0
        while time.perf_counter() - t0 < 8.0:
            m.tendencies(*args); k += 1
        dt = (time.perf_counter() - t0) / k
        cells = nx * ny * nz
        print(f"{name:18s} threads {n:2d}: {dt*1e3:8.1f} ms per tendency call, {dt/cells*1e9:6.1f} ns per cell", flush=True)

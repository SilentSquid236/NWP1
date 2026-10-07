"""
Cost of the non-hydrostatic core (nh3d.NH3D) on sized grids -- CAM stage S3.

    python tools/bench_nh3d.py --nx 440 --ny 388 --dx 3000 --nz 40 --dt 20 --steps 5

Dynamics only (no physics hooks): a 1500 m terrain bump and a warm bubble so
every term is non-zero. Reports the set-up time, seconds per large step, the
cost per cell per step, the peak memory, and the functions that take the most
time (cProfile, own time). The 24 h projection is steps * (86400 / dt).

Run it on one core (OMP_NUM_THREADS=1): NumPy element-wise work is single
threaded, so this is the cost of the core as it is written today.
"""
import argparse
import cProfile
import io
import os
import pstats
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path[:0] = [str(ROOT / "src" / "dynamics"), str(ROOT / "src")]

import nh3d                       # noqa: E402
from grid import CGrid            # noqa: E402
from sigma import SigmaLevels     # noqa: E402


def peak_rss_gb():
    try:
        import resource
        return resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1024 ** 2
    except ImportError:            # Windows
        return float("nan")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nx", type=int, default=110)
    p.add_argument("--ny", type=int, default=97)
    p.add_argument("--dx", type=float, default=12000.0)
    p.add_argument("--nz", type=int, default=20)
    p.add_argument("--dt", type=float, default=60.0)
    p.add_argument("--ns", type=int, default=6)
    p.add_argument("--steps", type=int, default=5)
    p.add_argument("--profile", action="store_true")
    p.add_argument("--backend", choices=("numpy", "c"), default="numpy")
    p.add_argument("--threads", type=int, default=0, help="OpenMP threads for --backend c")
    a = p.parse_args()

    t0 = time.time()
    g = CGrid(a.nx, a.ny, a.dx, a.dx, f0=1e-4, beta=0.0, edge_mode="periodic")
    X, Y = g.Xc - g.Lx / 2, g.Yc - g.Ly / 2
    L = 20.0 * a.dx
    m = nh3d.NH3D(g, SigmaLevels(a.nz), terrain=1500.0 * np.exp(-(X ** 2 + Y ** 2) / L ** 2),
                  ns=a.ns, backend=a.backend, threads=a.threads)
    th = m.Th / m.mu[None]
    th = th + 2.0 * np.exp(-((X + 3 * L) ** 2 + Y ** 2) / L ** 2)[None]
    m.Th = th * m.mu[None]
    m.U = m.U + 10.0 * m.h_to_u(m.mu)[None]
    t_setup = time.time() - t0

    m.step(a.dt)                                   # warm-up
    prof = cProfile.Profile() if a.profile else None
    if prof:
        prof.enable()
    t0 = time.time()
    for _ in range(a.steps):
        m.step(a.dt)
    t_step = (time.time() - t0) / a.steps
    if prof:
        prof.disable()
    assert np.isfinite(m.U).all(), "non-finite state"

    cells = a.nx * a.ny * a.nz
    steps24 = 86400.0 / a.dt
    print(f"grid {a.nx}x{a.ny}x{a.nz} dx {a.dx/1000:g} km  dt {a.dt:g} s  ns {a.ns}  "
          f"backend {a.backend}  threads {getattr(m, 'threads', 1) if a.backend == 'c' else 1}")
    print(f"  set-up {t_setup:.1f} s  step {t_step:.3f} s  "
          f"{t_step / cells * 1e9:.0f} ns per cell-step  "
          f"24 h = {steps24:.0f} steps = {t_step * steps24 / 3600:.2f} h  "
          f"peak RSS {peak_rss_gb():.2f} GB")
    if prof:
        s = io.StringIO()
        pstats.Stats(prof, stream=s).sort_stats("tottime").print_stats(14)
        lines = [l for l in s.getvalue().splitlines() if l.strip()]
        start = next(i for i, l in enumerate(lines) if "ncalls" in l)
        print("\n".join(lines[start:start + 15]))


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Would more cores make a forecast faster? Measure before porting anything.

    python tools/bench_threads.py

The dynamical core is element-wise NumPy, which uses one core. torch is
already installed on the server and runs element-wise work on several
threads, so the question is how much that buys ON ARRAYS THIS SMALL: one
model field is 20 x 97 x 110 = 213k values (1.7 MB), and thread start-up is
paid on every operation, thousands of times per forecast hour.

Times three operation shapes the profile (research log 2026-09-22) says
dominate a step -- a 5-point Laplacian (hyperdiffusion, advection), a chain
of element-wise arithmetic (tendencies, Richardson number), and a vertical
cumulative sum (hydrostatic geopotential) -- in NumPy and in torch at 1..N
threads. N never exceeds the resources.py ceiling (50 % of cores, less when
the machine is busy). Takes well under a minute.
"""

import sys
import time
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
import resources
PLAN = resources.apply()

import numpy as np

SHAPE = (20, 97, 110)
REPS = 200


def lap_np(a, o):
    o[..., 1:-1, 1:-1] = (a[..., 1:-1, 2:] + a[..., 1:-1, :-2] + a[..., 2:, 1:-1]
                          + a[..., :-2, 1:-1] - 4.0 * a[..., 1:-1, 1:-1])
    return o


def chain_np(a, b, c):
    return (a * b + c) * np.exp(-0.01 * a) - 0.5 * b * b / (1.0 + c * c)


def csum_np(a):
    return np.cumsum(a[::-1], axis=0)[::-1]


def timed(fn, reps=REPS):
    fn()
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t0) / reps * 1e3


def main():
    ceiling = int(PLAN["torch_threads"])
    print(f"cores {PLAN['total_cores']}, thread ceiling {ceiling} "
          f"(resources.py), arrays {SHAPE} float64\n")
    rng = np.random.default_rng(0)
    a, b, c = (rng.random(SHAPE) for _ in range(3))
    o = np.empty_like(a)
    base = {"laplacian": timed(lambda: lap_np(a, o)),
            "chain": timed(lambda: chain_np(a, b, c)),
            "cumsum": timed(lambda: csum_np(a))}
    print(f"{'':12s}{'numpy':>9s}")
    for k, v in base.items():
        print(f"{k:12s}{v:8.3f} ms")

    try:
        import torch
    except ImportError:
        print("\ntorch not importable here; nothing more to measure")
        return 0
    ta, tb, tc = (torch.from_numpy(x) for x in (a, b, c))
    to = torch.empty_like(ta)

    def lap_t():
        to[..., 1:-1, 1:-1] = (ta[..., 1:-1, 2:] + ta[..., 1:-1, :-2] + ta[..., 2:, 1:-1]
                               + ta[..., :-2, 1:-1] - 4.0 * ta[..., 1:-1, 1:-1])

    def chain_t():
        return (ta * tb + tc) * torch.exp(-0.01 * ta) - 0.5 * tb * tb / (1.0 + tc * tc)

    def csum_t():
        return torch.flip(torch.cumsum(torch.flip(ta, [0]), 0), [0])

    threads = [t for t in (1, 2, 4, 8, 16, 26, 52) if t <= ceiling]
    print(f"\ntorch {torch.__version__}: time (ms) and speed-up over numpy")
    print(f"{'threads':>8s}" + "".join(f"{k:>18s}" for k in base))
    for n in threads:
        torch.set_num_threads(n)
        row = [timed(f) for f in (lap_t, chain_t, csum_t)]
        print(f"{n:8d}" + "".join(f"{t:9.3f} ({base[k]/t:4.1f}x)"
                                  for t, k in zip(row, base)))
    print("\nA whole-core port is worth it only if the best row is well above "
          "2x on all three: every model operation would pay the same overhead.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

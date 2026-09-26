#!/usr/bin/env python3
"""
Is the torch backend right and fast on THIS machine?  (Run it on the server.)

    python tools/check_backend.py            # full 110 x 97 grid, 60 steps
    python tools/check_backend.py --threads 4,8,12,16

Runs the same realistic case (src/dynamics/backend_reference.py) with NumPy
and with torch at several thread counts, and prints the time, the speed-up
and the largest relative difference from NumPy. Every torch row must show a
difference at round-off (below 1e-9); the fastest thread count is the one to
give forecast.py --threads. Uses at most the threads asked for.
"""
import argparse
import sys
import time
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "dynamics"))

import backend                               # noqa: E402
import backend_reference as br              # noqa: E402


def run(name, steps, threads=None):
    m = br.build(110, 97)
    if name != "numpy":
        m.to_backend(name, threads=threads)
    m.step(20.0)                             # warm-up (allocation, thread pool)
    t0 = time.time()
    for _ in range(steps):
        m.step(20.0)
    wall = time.time() - t0
    return {k: np.asarray(backend.to_numpy(getattr(m, k))) for k in ("u", "v", "theta", "pi")}, wall


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("--steps", type=int, default=60)
    ap.add_argument("--threads", default="1,4,8,12,16")
    a = ap.parse_args()
    if backend.TORCH is None:
        sys.exit("PyTorch is not importable here; only the numpy backend can run.")
    import torch
    print(f"torch {torch.__version__}, numpy {np.__version__}, "
          f"{torch.get_num_threads()} threads by default")
    ref, t_np = run("numpy", a.steps)
    print(f"{'backend':>14} {'s/step':>8} {'speed-up':>9} {'max rel diff':>13}")
    print(f"{'numpy':>14} {t_np / a.steps:8.3f} {1.0:9.1f} {0.0:13.1e}")
    for th in [int(x) for x in a.threads.split(",") if x.strip()]:
        out, t = run("torch", a.steps, threads=th)
        rel = max(float(np.abs(out[k] - ref[k]).max() / np.abs(ref[k]).max()) for k in ref)
        print(f"{'torch x' + str(th):>14} {t / a.steps:8.3f} {t_np / t:9.1f} {rel:13.1e}")
    print(f"\nA 24 h forecast is ~{24 * 3600 / 17.1:.0f} steps: numpy "
          f"~{24 * 3600 / 17.1 * t_np / a.steps / 60:.0f} min at these speeds.")


if __name__ == "__main__":
    main()

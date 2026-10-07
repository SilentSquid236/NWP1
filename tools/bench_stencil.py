"""
Backend options for the CAM core, on one representative stencil -- CAM stage S3.

The NH core's cost is dominated by short, memory-bound 3-D stencils (the
acoustic substep, the flux divergences). This times ONE such stencil

    out = a + dt * ( (b[i+1]-b[i-1]) * c  +  (b[j+1]-b[j-1]) * c
                     + 0.25 * (a[k+1] - 2 a[k] + a[k-1]) )

on a 3 km x 40 level array (440 x 388 x 40) with each candidate back end:

  numpy   float64 / float32        (today's code; one core)
  torch   float64 / float32, 1..N threads   (the hydrostatic core's backend)
  C       float64 / float32, OpenMP 1..N threads, compiled here with gcc -O3
          and called through ctypes (no package is installed)

All back ends are checked against the NumPy float64 answer before timing.
Usage:  python tools/bench_stencil.py [--threads 1,8,26,52] [--reps 20]
"""
import argparse
import ctypes
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

C_SRC = r"""
#include <stddef.h>
#define KERNEL(NAME, T)                                                        \
void NAME(const T *a, const T *b, const T *c, T *out,                          \
          int nz, int ny, int nx, T dt) {                                      \
    _Pragma("omp parallel for collapse(2) schedule(static)")                  \
    for (int k = 0; k < nz; k++)                                               \
    for (int j = 0; j < ny; j++) {                                             \
        int km = k > 0 ? k - 1 : k, kp = k < nz - 1 ? k + 1 : k;              \
        int jm = (j - 1 + ny) % ny, jp = (j + 1) % ny;                         \
        const T *ar = a + ((size_t)k * ny + j) * nx;                           \
        const T *br = b + ((size_t)k * ny + j) * nx;                           \
        const T *bn = b + ((size_t)k * ny + jp) * nx;                          \
        const T *bs = b + ((size_t)k * ny + jm) * nx;                          \
        const T *au = a + ((size_t)kp * ny + j) * nx;                          \
        const T *ad = a + ((size_t)km * ny + j) * nx;                          \
        const T *cr = c + ((size_t)k * ny + j) * nx;                           \
        T *o = out + ((size_t)k * ny + j) * nx;                                \
        for (int i = 0; i < nx; i++) {                                         \
            int im = i > 0 ? i - 1 : nx - 1, ip = i < nx - 1 ? i + 1 : 0;      \
            o[i] = ar[i] + dt * ((br[ip] - br[im]) * cr[i]                     \
                                 + (bn[i] - bs[i]) * cr[i]                     \
                                 + (T)0.25 * (au[i] - 2 * ar[i] + ad[i]));     \
        }                                                                      \
    }                                                                          \
}
KERNEL(stencil_f64, double)
KERNEL(stencil_f32, float)
"""


def numpy_kernel(a, b, c, dt):
    up = np.concatenate([a[1:], a[-1:]], axis=0)
    dn = np.concatenate([a[:1], a[:-1]], axis=0)
    dx = np.roll(b, -1, axis=2) - np.roll(b, 1, axis=2)
    dy = np.roll(b, -1, axis=1) - np.roll(b, 1, axis=1)
    return a + dt * (dx * c + dy * c + 0.25 * (up - 2 * a + dn))


def torch_kernel(torch):
    def k(a, b, c, dt):
        up = torch.cat([a[1:], a[-1:]], dim=0)
        dn = torch.cat([a[:1], a[:-1]], dim=0)
        dx = torch.roll(b, -1, 2) - torch.roll(b, 1, 2)
        dy = torch.roll(b, -1, 1) - torch.roll(b, 1, 1)
        return a + dt * (dx * c + dy * c + 0.25 * (up - 2 * a + dn))
    return k


def build_c():
    d = Path(tempfile.mkdtemp(prefix="nwp_stencil_"))
    (d / "k.c").write_text(C_SRC)
    so = d / "k.so"
    r = subprocess.run(["gcc", "-O3", "-march=native", "-fopenmp", "-shared", "-fPIC",
                        str(d / "k.c"), "-o", str(so)], capture_output=True, text=True)
    if r.returncode:
        raise RuntimeError(r.stderr)
    return ctypes.CDLL(str(so))


def timeit(fn, reps):
    fn()
    t0 = time.perf_counter()
    for _ in range(reps):
        fn()
    return (time.perf_counter() - t0) / reps


def main():
    p = argparse.ArgumentParser()
    p.add_argument("--nx", type=int, default=440)
    p.add_argument("--ny", type=int, default=388)
    p.add_argument("--nz", type=int, default=40)
    p.add_argument("--threads", default="1,8,26,52")
    p.add_argument("--reps", type=int, default=20)
    a_ = p.parse_args()
    threads = [int(t) for t in a_.threads.split(",")]
    shape = (a_.nz, a_.ny, a_.nx)
    n = np.prod(shape)
    rng = np.random.default_rng(0)
    A, B, C = (rng.standard_normal(shape) for _ in range(3))
    dt = 0.1
    ref = numpy_kernel(A, B, C, dt)
    rows = []

    def row(name, nthr, sec, err):
        rows.append((name, nthr, sec))
        print(f"  {name:14s} threads {nthr:3d}   {sec*1e3:8.2f} ms   {sec/n*1e9:6.2f} ns/cell   "
              f"max err {err:.1e}", flush=True)

    print(f"stencil on {a_.nx}x{a_.ny}x{a_.nz} = {n/1e6:.2f} M cells")
    for dtype in (np.float64, np.float32):
        a, b, c = A.astype(dtype), B.astype(dtype), C.astype(dtype)
        out = numpy_kernel(a, b, c, dtype(dt))
        row(f"numpy {np.dtype(dtype).name}", 1,
            timeit(lambda: numpy_kernel(a, b, c, dtype(dt)), a_.reps),
            float(np.abs(out - ref).max()))

    try:
        import torch
        k = torch_kernel(torch)
        for dtype in (torch.float64, torch.float32):
            a, b, c = (torch.from_numpy(x).to(dtype) for x in (A, B, C))
            for t in threads:
                torch.set_num_threads(t)
                out = k(a, b, c, dt)
                row(f"torch {str(dtype).split('.')[-1]}", t,
                    timeit(lambda: k(a, b, c, dt), a_.reps),
                    float(np.abs(out.double().numpy() - ref).max()))
    except ImportError:
        print("  torch: not installed")

    try:
        lib = build_c()
    except Exception as e:            # noqa: BLE001
        print(f"  C: could not build ({e})")
        return
    for dtype, fname, ct in ((np.float64, "stencil_f64", ctypes.c_double),
                             (np.float32, "stencil_f32", ctypes.c_float)):
        f = getattr(lib, fname)
        a, b, c = (np.ascontiguousarray(x.astype(dtype)) for x in (A, B, C))
        out = np.empty_like(a)
        ptr = lambda x: x.ctypes.data_as(ctypes.c_void_p)
        call = lambda: f(ptr(a), ptr(b), ptr(c), ptr(out),
                         ctypes.c_int(a_.nz), ctypes.c_int(a_.ny), ctypes.c_int(a_.nx), ct(dt))
        omp = ctypes.CDLL("libgomp.so.1")
        for t in threads:
            omp.omp_set_num_threads(ctypes.c_int(t))
            call()
            row(f"C {np.dtype(dtype).name}", t, timeit(call, a_.reps),
                float(np.abs(out.astype(np.float64) - ref).max()))


if __name__ == "__main__":
    main()

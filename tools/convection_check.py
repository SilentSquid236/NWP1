#!/usr/bin/env python3
"""
Where is a forecast state statically unstable, and how many sweeps does the
convective adjustment need to remove it?

    python tools/convection_check.py FORECAST.npz [--snapshot -1]

Takes one snapshot (default: the last) and reports, for theta increasing
upward (index 0 is the lid, so an interface k/k+1 is unstable when
theta[k] < theta[k+1]):
  * the unstable fraction of interfaces and the number of columns affected;
  * how the unstable interfaces split by level and by distance from the
    nearest lateral edge (the relaxation zone is the outer 15 cells by default);
  * the largest violation in K;
  * sweeps used and the fraction left unstable when
    convection.dry_convective_adjustment runs with caps of 20, 50, 200 and
    1000 sweeps (the model uses 20).
NumPy only; nothing is written.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src" / "dynamics"))

from convection import dry_convective_adjustment      # noqa: E402
from sigma import SigmaLevels                          # noqa: E402

TOL = 1e-10


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("forecast")
    ap.add_argument("--snapshot", type=int, default=-1)
    a = ap.parse_args()
    z = np.load(a.forecast, allow_pickle=True)
    i = a.snapshot
    th = np.asarray(z["theta"][i], dtype=float)
    u = np.asarray(z["u"][i], dtype=float)
    v = np.asarray(z["v"][i], dtype=float)
    pi = np.asarray(z["pi"][i], dtype=float)
    lev = SigmaLevels(th.shape[0], p_top=float(z["p_top"]))
    nz, ny, nx = th.shape
    t_h = float(z["times_s"][i]) / 3600.0

    bad = th[:-1] < th[1:] - TOL                           # (nz-1, ny, nx)
    viol = np.where(bad, th[1:] - th[:-1], 0.0)
    cols = bad.any(axis=0)
    jj, ii = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    edge = np.minimum(np.minimum(jj, ii), np.minimum(ny - 1 - jj, nx - 1 - ii))
    print(f"{a.forecast}  snapshot t+{t_h:.2f} h")
    print(f"unstable interfaces {bad.mean():.3e} ({int(bad.sum())} of {bad.size}); "
          f"columns with any: {int(cols.sum())} of {cols.size}; "
          f"largest violation {viol.max():.3g} K")
    if bad.any():
        per_level = bad.sum(axis=(1, 2))
        print("by interface: " + "  ".join(f"L{k:02d}/{k + 1:02d} {int(n)}"
                                           for k, n in enumerate(per_level) if n))
        e = np.broadcast_to(edge, bad.shape)[bad]
        bins = [(0, 4), (5, 9), (10, 14), (15, 19), (20, 10 ** 6)]
        print("by edge distance (cells): " + "  ".join(
            f"{lo}-{hi if hi < 10 ** 6 else '':}: {int(((e >= lo) & (e <= hi)).sum())}"
            for lo, hi in bins))
    for cap in (20, 50, 200, 1000):
        _, _, _, info = dry_convective_adjustment(th, u, v, pi, lev, max_sweeps=cap)
        print(f"adjustment cap {cap:4d}: sweeps {info['sweeps']:4d}, "
              f"unstable after {info['unstable_after']:.3e}")


if __name__ == "__main__":
    main()

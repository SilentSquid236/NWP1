#!/usr/bin/env python3
"""
Compare two forecast files hour by hour (e.g. the numpy and torch backends).

    python tools/compare_forecasts.py A.npz B.npz

For each output time both files share (within 0.01 h), prints the largest
relative difference over u, v, theta and pi, then each file's stop status
and wall time. Needs only NumPy.
"""
import sys

import numpy as np


def main(pa, pb):
    a, b = np.load(pa, allow_pickle=True), np.load(pb, allow_pickle=True)
    ta, tb = a["times_s"] / 3600.0, b["times_s"] / 3600.0
    print(f"{'hour':>6}  max relative difference (u, v, theta, pi)")
    for i, h in enumerate(tb):
        j = int(np.argmin(np.abs(ta - h)))
        if abs(ta[j] - h) > 0.01:
            continue
        d = max(float(np.abs(a[k][j] - b[k][i]).max() / max(np.abs(a[k][j]).max(), 1e-30))
                for k in ("u", "v", "theta", "pi"))
        print(f"{h:6.2f}  {d:.1e}")
    for name, z in ((pa, a), (pb, b)):
        print(f"{name}: {z['stopped']}, {float(z['wall_s']) / 60:.1f} min")


if __name__ == "__main__":
    if len(sys.argv) != 3:
        sys.exit(__doc__)
    main(sys.argv[1], sys.argv[2])

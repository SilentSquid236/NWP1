#!/usr/bin/env python3
"""Where does a forecast start to run away?  Reads a saved forecast.npz.

    python tools/locate_growth.py DATA/tensors_3d/obs_<stamp>/forecast.npz

For each pair of consecutive snapshots it prints the largest change in u and
v: its size, level, row, column, distance from the nearest edge (in cells)
and terrain height there, and how many points changed by more than 5 m/s,
with their minimum and median edge distance.  The relaxation zone is the
outer 10 cells, so an edge distance of about 10 is its inner boundary.
It also prints the first-guess line from availability.json, if present.
No model code is imported; only NumPy is needed.
"""
import json
import sys
from pathlib import Path

import numpy as np


def edge_dist(shape, r, c):
    ny, nx = shape
    return int(min(r, c, ny - 1 - r, nx - 1 - c))


def describe(d, terrain):
    """Largest |d| (d is [L, Y, X]) and the spread of big changes."""
    a = np.abs(np.nan_to_num(d, nan=0.0, posinf=1e9, neginf=1e9))
    k, r, c = np.unravel_index(int(np.argmax(a)), a.shape)
    ny, nx = a.shape[1:]
    tr = terrain[min(r, terrain.shape[0] - 1), min(c, terrain.shape[1] - 1)]
    big = np.argwhere(a > 5.0)
    if len(big):
        ed = np.minimum.reduce([big[:, 1], big[:, 2],
                                ny - 1 - big[:, 1], nx - 1 - big[:, 2]])
        spread = f"{len(big):5d} pts>5  edge min {ed.min():2d} med {int(np.median(ed)):2d}"
    else:
        spread = "    0 pts>5"
    return (f"{a[k, r, c]:6.1f} at L{k:02d} r{r:3d} c{c:3d} "
            f"edge {edge_dist((ny, nx), r, c):2d} z {tr:5.0f} m  {spread}"
            + (f"  NON-FINITE {int((~np.isfinite(d)).sum())}" if not np.isfinite(d).all() else ""))


def main(path):
    path = Path(path)
    f = np.load(path, allow_pickle=True)
    t = f["times_s"] / 3600.0
    u, v, th = f["u"], f["v"], f["theta"]
    terrain = np.asarray(f["terrain"], dtype=float)
    print(f"{path}  snapshots {len(t)}  u {u.shape[1:]}  v {v.shape[1:]}")
    av = path.parent / "availability.json"
    if av.exists():
        print("first guess:", json.loads(av.read_text()).get("first_guess"))
    print(f"{'hours':>11}  {'max|u|':>6} {'max|v|':>6} {'min th':>6}   "
          "largest change in u / v")
    for i in range(1, len(t)):
        du, dv = u[i] - u[i - 1], v[i] - v[i - 1]
        print(f"{t[i-1]:5.2f}-{t[i]:5.2f}  {np.nanmax(np.abs(u[i])):6.1f} "
              f"{np.nanmax(np.abs(v[i])):6.1f} {np.nanmin(th[i]):6.1f}   "
              f"u {describe(du, terrain)}")
        print(f"{'':35}v {describe(dv, terrain)}")


if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(__doc__)
    main(sys.argv[1])

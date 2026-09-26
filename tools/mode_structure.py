#!/usr/bin/env python3
"""
Where does a forecast's fastest-growing disturbance live, and is the flow
there unstable to begin with?

    python tools/mode_structure.py A.npz B.npz [--hours 1,2,3,4] [--box ROW COL]
                                                [--at 0,2,4]

A and B are two runs of the SAME case that differ only by round-off, e.g.
the numpy and torch backends (P-60). Their difference is then the leading
unstable mode seeded by round-off (a bred vector): it shows where the
instability lives hours before it is visible in the fields.

Part 1 prints, for each requested hour, where |B - A| is largest for u, v
and theta (level, row, column, edge distance, lat/lon), which three levels
carry most of the wind difference in a box around that point, how rough
the u difference is there ("rough": the share of its variance a 3x3 mean
removes, about 1 for a 2dx checkerboard and about 0 for a smooth mode),
and the e-folding time of the wind difference since the previous hour.

Part 2 prints, level by level, the stability of A's flow in a box around
the mode (the last requested hour's u/v maximum, or --box ROW COL), at
A's first snapshot or at each of the --at hours:
  eta/f  absolute vorticity over f; below 0 is inertial instability
  Ri     gradient Richardson number with the level below; below 0.25 is
         shear instability, below 0 static instability
  N2     buoyancy frequency squared with the level below (1/s^2)
and, over the whole domain, how many points per level have eta/f < 0 and
how many have Ri < 0.25. Rows are stored south to north, level 0 is the lid.
Needs only NumPy.
"""
import argparse
import sys

import numpy as np

RD, CP, G, P0 = 287.05, 1004.6, 9.80665, 1.0e5
KAPPA = RD / CP
OMEGA = 7.2921e-5


def edge(shape, r, c):
    ny, nx = shape
    return int(min(r, c, ny - 1 - r, nx - 1 - c))


def roughness(d, r, c, h):
    """Share of the variance of d (2D) in a box that a 3x3 mean removes.

    About 1 for a 2dx checkerboard, about 0 for a field smooth on the grid
    scale. It tells a grid-scale mode from a resolved one.
    """
    ny, nx = d.shape
    r0, r1 = max(r - h, 1), min(r + h + 1, ny - 1)
    c0, c1 = max(c - h, 1), min(c + h + 1, nx - 1)
    raw = d[r0:r1, c0:c1]
    sm = sum(d[r0 + i:r1 + i, c0 + j:c1 + j] for i in (-1, 0, 1) for j in (-1, 0, 1)) / 9.0
    v = float(raw.var())
    return 1.0 - float(sm.var()) / v if v > 0 else float("nan")


def centre_winds(u, v):
    """C-grid u (west faces) and v (south faces) to cell centres."""
    uc = u.copy()
    uc[..., :-1] = 0.5 * (u[..., :-1] + u[..., 1:])
    vc = v.copy()
    vc[..., :-1, :] = 0.5 * (v[..., :-1, :] + v[..., 1:, :])
    return uc, vc


def stability(u, v, theta, pi, sigma, p_top, lat, lon):
    """eta/f per level; N2 and Ri per interface (level k with k+1 below)."""
    uc, vc = centre_winds(u, v)
    dlat = np.gradient(lat, axis=0) * 111_132.0
    dlon = np.gradient(lon, axis=1) * 111_320.0 * np.cos(np.radians(lat))
    zeta = np.gradient(vc, axis=2) / dlon - np.gradient(uc, axis=1) / dlat
    f = 2 * OMEGA * np.sin(np.radians(lat))
    eta_f = 1.0 + zeta / f
    p = sigma[:, None, None] * pi[None] + p_top
    T = theta * (p / P0) ** KAPPA
    dz = RD * 0.5 * (T[:-1] + T[1:]) / G * np.log(p[1:] / p[:-1])
    n2 = G * (theta[:-1] - theta[1:]) / (0.5 * (theta[:-1] + theta[1:]) * dz)
    sh2 = ((uc[:-1] - uc[1:]) ** 2 + (vc[:-1] - vc[1:]) ** 2) / dz ** 2
    ri = n2 / np.maximum(sh2, 1e-12)
    speed = np.hypot(uc, vc)
    return eta_f, n2, ri, speed, p


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--hours", default="1,2,3,4")
    ap.add_argument("--box", nargs=2, type=int, metavar=("ROW", "COL"))
    ap.add_argument("--half", type=int, default=6, help="box half-width, cells")
    ap.add_argument("--at", default=None,
                    help="hours of A's snapshots for Part 2 (default: the first)")
    a = ap.parse_args()

    A, B = np.load(a.a, allow_pickle=True), np.load(a.b, allow_pickle=True)
    ta, tb = A["times_s"] / 3600.0, B["times_s"] / 3600.0
    lat, lon = A["lat"], A["lon"]
    shape = lat.shape
    h = a.half

    print("Part 1 -- the difference B - A")
    print(f"{'hour':>5} {'var':>5} {'max|d|':>9}  where (level, row, col, edge, lat, lon)"
          f"   top levels of wind difference in box   e-fold")
    prev = None
    centre = None
    for hr in [float(x) for x in a.hours.split(",")]:
        i = int(np.argmin(np.abs(ta - hr)))
        j = int(np.argmin(np.abs(tb - ta[i])))
        if abs(ta[i] - hr) > 0.13 or abs(tb[j] - ta[i]) > 0.01:
            print(f"{hr:5.2f}  no common snapshot")
            continue
        du = B["u"][j] - A["u"][i]
        dv = B["v"][j] - A["v"][i]
        dth = B["theta"][j] - A["theta"][i]
        wind = np.hypot(du, dv)
        for name, d in (("u", du), ("v", dv), ("theta", dth)):
            k, r, c = np.unravel_index(int(np.argmax(np.abs(d))), d.shape)
            extra = ""
            if name == "u":
                kw, rw, cw = np.unravel_index(int(np.argmax(wind)), wind.shape)
                centre = (rw, cw)
                box = wind[:, max(rw - h, 0):rw + h + 1, max(cw - h, 0):cw + h + 1]
                prof = np.sqrt((box ** 2).mean(axis=(1, 2)))
                top = np.argsort(prof)[::-1][:3]
                extra = "  " + " ".join(f"L{t:02d} {prof[t] / prof.max():4.0%}" for t in top)
                extra += f"  rough {roughness(du[kw], rw, cw, h):.2f}"
                wmax = float(wind.max())
                if prev is not None and wmax > 0 and prev[1] > 0 and wmax != prev[1]:
                    tau = (ta[i] - prev[0]) * 60.0 / np.log(wmax / prev[1])
                    extra += f"   {tau:5.0f} min"
                prev = (ta[i], wmax)
            print(f"{ta[i]:5.2f} {name:>5} {abs(d[k, r, c]):9.1e}  L{k:02d} r{r:3d} c{c:3d} "
                  f"edge {edge(shape, r, c):2d} {lat[r, c]:5.2f}N {lon[r, c]:7.2f}{extra}")

    if a.box:
        centre = tuple(a.box)
    if centre is None:
        sys.exit("no hour matched; nothing to centre Part 2 on")
    r0, c0 = centre
    sl = (slice(max(r0 - h, 0), r0 + h + 1), slice(max(c0 - h, 0), c0 + h + 1))
    idx = [0] if a.at is None else sorted({int(np.argmin(np.abs(ta - float(x))))
                                           for x in a.at.split(",")})
    for n in idx:
        eta_f, n2, ri, speed, p = stability(A["u"][n], A["v"][n], A["theta"][n], A["pi"][n],
                                            A["sigma"], float(A["p_top"]), lat, lon)
        print(f"\nPart 2 -- stability in A at t+{ta[n]:.2f} h, box of "
              f"+-{h} cells round r{r0} c{c0} ({lat[r0, c0]:.2f}N {lon[r0, c0]:.2f}, "
              f"edge {edge(shape, r0, c0)})")
        print("  lev    p hPa  max wind  min eta/f  min Ri(k,k+1)  min N2(k,k+1) | domain: eta/f<0  Ri<0.25")
        nz = eta_f.shape[0]
        for k in range(nz):
            pb = float(p[k][sl].mean()) / 100.0
            line = (f"  L{k:02d} {pb:8.0f} {float(speed[k][sl].max()):9.1f} "
                    f"{float(eta_f[k][sl].min()):10.2f}")
            if k < nz - 1:
                line += (f" {float(ri[k][sl].min()):14.2f} {float(n2[k][sl].min()):14.1e} |"
                         f" {int((eta_f[k] < 0).sum()):8d} {int((ri[k] < 0.25).sum()):8d}")
            else:
                line += f" {'':14} {'':14} | {int((eta_f[k] < 0).sum()):8d}"
            print(line)

if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""
Which term of the equations feeds the fastest-growing mode?

    python tools/mode_budget.py A.npz B.npz --hour 3 [--box ROW COL] [--half 8]

A and B are two runs of the same case that differ only by round-off (the
numpy and torch backends). At the chosen hour their difference d = B - A
is the leading unstable mode (tools/mode_structure.py; the bred-vector
idea of Toth and Kalnay 1993). This tool scales d
to about 1 mm/s, so it stays linear, and evaluates every tendency term of
the core at A and at A + d. It prints each term's contribution to the
growth of the mode's kinetic energy in a box round the mode:

    P_term = sum over the box of  du * dT_u + dv * dT_v

where dT is the term's tendency at A + d minus its tendency at A.
Horizontal and vertical advection are split into the base flow carrying
the mode (U.grad d) and the mode carrying the base flow (d.grad U), and
the parts are checked against the full difference.

Two checks come first:
- the terms must add to the core's own tendency (`PrimitiveSigma.tendencies`);
- the implied energy growth rate, sum(P) / sum(du^2 + dv^2), must be
  compared with 2 / (e-folding time) from `mode_structure.py`.

Rows run south to north; level 0 is the lid.
"""
import argparse
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))
sys.path.insert(0, str(ROOT / "src" / "dynamics"))

from grid import CGrid                                         # noqa: E402
from sigma import (SigmaLevels, hydrostatic_geopotential,       # noqa: E402
                   continuity, vertical_advection, pressure_gradient_force,
                   KAPPA, P0)
from subgrid import hyperdiffusion                              # noqa: E402
from turbulence import vertical_mixing                          # noqa: E402
from surface import surface_drag                                # noqa: E402
from primitive_sigma import PrimitiveSigma                      # noqa: E402


def build_grid(lat, lon):
    """The grid forecast.build_grid makes, from the file's own lat/lon."""
    ny, nx = lat.shape
    lat_min = float(lat[0, 0] - 0.5 * (lat[1, 0] - lat[0, 0]))
    lat_max = float(lat[-1, 0] + 0.5 * (lat[-1, 0] - lat[-2, 0]))
    lon_min = float(lon[0, 0] - 0.5 * (lon[0, 1] - lon[0, 0]))
    lon_max = float(lon[0, -1] + 0.5 * (lon[0, -1] - lon[0, -2]))
    lat0 = 0.5 * (lat_min + lat_max)
    dy = (lat_max - lat_min) * 111_132.0 / ny
    dx = (lon_max - lon_min) * 111_320.0 * np.cos(np.radians(lat0)) / nx
    omega = 7.2921e-5
    f0 = 2 * omega * np.sin(np.radians(lat0))
    beta = 2 * omega * np.cos(np.radians(lat0)) / 6_371_000.0
    return CGrid(nx, ny, dx, dy, f0=f0, beta=beta, edge_mode="replicate")


def terms(m, u, v, th, pi):
    """Every tendency term of PrimitiveSigma.tendencies, kept apart."""
    gr, lev = m.grid, m.lev
    phi = hydrostatic_geopotential(th, pi, lev, phi_surface=m.phi_s)
    _, sd = continuity(u, v, pi, lev, gr)
    T_ref = np.mean(th * (lev.pressure(pi) / P0) ** KAPPA, axis=(1, 2))
    fx, fy = pressure_gradient_force(phi, th, pi, lev, gr,
                                     reference=T_ref if m.ref_pgf else None)
    v_at_u, u_at_v = gr.v_to_u(v), gr.u_to_v(u)
    t = {
        "horizontal advection": (-m._horiz_adv(u, u, v_at_u), -m._horiz_adv(v, u_at_v, v)),
        "vertical advection": (-vertical_advection(u, sd, lev), -vertical_advection(v, sd, lev)),
        "Coriolis": (gr.f_u * v_at_u, -gr.f_v * u_at_v),
        "pressure gradient": (fx, fy),
    }
    z = np.zeros_like(u)
    t["hyperdiffusion"] = ((hyperdiffusion(u, gr, m.hyper), hyperdiffusion(v, gr, m.hyper))
                           if m.hyper > 0 else (z, z))
    if m.drag:
        ddu, ddv, _ = surface_drag(u, v, th, pi, lev, z0=m.z0, theta_s=m.theta_surface)
        t["surface drag"] = (ddu, ddv)
    if m.mixing:
        mu, mv, _, _ = vertical_mixing(u, v, th, pi, lev, ri_crit=m.ri_crit,
                                       k_max=m.k_max, mixing_length=m.mixing_length)
        t["vertical mixing"] = (mu, mv)
    if m.sponge_levels > 0:
        s = m._sponge
        t["sponge"] = (-s * (u - m._u_ref), -s * (v - m._v_ref))
    return t, sd


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("a")
    ap.add_argument("b")
    ap.add_argument("--hour", type=float, default=3.0)
    ap.add_argument("--box", nargs=2, type=int, metavar=("ROW", "COL"))
    ap.add_argument("--half", type=int, default=8)
    a = ap.parse_args()

    A, B = np.load(a.a, allow_pickle=True), np.load(a.b, allow_pickle=True)
    ta, tb = A["times_s"] / 3600.0, B["times_s"] / 3600.0
    i = int(np.argmin(np.abs(ta - a.hour)))
    j = int(np.argmin(np.abs(tb - ta[i])))
    if abs(tb[j] - ta[i]) > 0.01:
        sys.exit("no common snapshot at that hour")
    sA = {k: np.asarray(A[k][i], dtype=float) for k in ("u", "v", "theta", "pi")}
    d = {k: np.asarray(B[k][j], dtype=float) - sA[k] for k in sA}
    wind = np.hypot(d["u"], d["v"])
    scale = 1e-3 / float(wind.max())
    d = {k: scale * x for k, x in d.items()}

    lev = SigmaLevels(len(A["sigma"]), p_top=float(A["p_top"]))
    if not np.allclose(lev.sigma, A["sigma"]):
        sys.exit("the file's sigma levels do not match SigmaLevels")
    grid = build_grid(A["lat"], A["lon"])
    m = PrimitiveSigma(grid, lev, terrain=np.asarray(A["terrain"], dtype=float))
    m.set_reference(sA["u"], sA["v"])

    if a.box:
        r0, c0 = a.box
    else:
        _, r0, c0 = np.unravel_index(int(np.argmax(wind)), wind.shape)
    h = a.half
    box = (slice(None), slice(max(r0 - h, 0), r0 + h + 1), slice(max(c0 - h, 0), c0 + h + 1))

    s1 = {k: sA[k] + d[k] for k in sA}
    tA, sdA = terms(m, sA["u"], sA["v"], sA["theta"], sA["pi"])
    t1, sd1 = terms(m, s1["u"], s1["v"], s1["theta"], s1["pi"])

    # Check 1: the separated terms rebuild the core's own tendency difference.
    duA, dvA, _, _ = m.tendencies(sA["u"], sA["v"], sA["theta"], sA["pi"])
    du1, dv1, _, _ = m.tendencies(s1["u"], s1["v"], s1["theta"], s1["pi"])
    core = (du1 - duA, dv1 - dvA)
    mine = (sum(t1[k][0] - tA[k][0] for k in tA), sum(t1[k][1] - tA[k][1] for k in tA))
    err = max(float(np.abs(mine[q] - core[q]).max()) for q in (0, 1))
    ref = max(float(np.abs(core[q]).max()) for q in (0, 1))
    print(f"hour {ta[i]:.2f}, box +-{h} round r{r0} c{c0} "
          f"({float(A['lat'][r0, c0]):.2f}N {float(A['lon'][r0, c0]):.2f}), mode scaled to 1 mm/s")
    print(f"check 1: terms rebuild the core's tendency difference to {err / ref:.1e} (relative)")

    du, dv = d["u"][box], d["v"][box]
    E = float((du ** 2 + dv ** 2).sum())
    P = {k: float((du * (t1[k][0] - tA[k][0])[box] + dv * (t1[k][1] - tA[k][1])[box]).sum())
         for k in tA}

    # Advection split: base flow carrying the mode, and the mode carrying the base flow.
    gr = m.grid
    uA, vA = sA["u"], sA["v"]
    dU, dV = d["u"], d["v"]
    split = {
        "  U.grad d (base carries mode)": (-m._horiz_adv(dU, uA, gr.v_to_u(vA)),
                                           -m._horiz_adv(dV, gr.u_to_v(uA), vA)),
        "  d.grad U (mode carries base)": (-m._horiz_adv(uA, dU, gr.v_to_u(dV)),
                                           -m._horiz_adv(vA, gr.u_to_v(dU), dV)),
        "  sigma_dot(A) d/dsigma d": (-vertical_advection(dU, sdA, lev),
                                      -vertical_advection(dV, sdA, lev)),
        "  d(sigma_dot) d/dsigma U": (-vertical_advection(uA, sd1 - sdA, lev),
                                      -vertical_advection(vA, sd1 - sdA, lev)),
    }
    Ps = {k: float((du * x[0][box] + dv * x[1][box]).sum()) for k, x in split.items()}

    tot = sum(P.values())
    print(f"check 2: energy growth rate sum(P)/E = {tot / E:.2e} 1/s "
          f"(an amplitude e-folding of {2.0 / (tot / E) / 60.0:.0f} min if positive)\n")
    print(f"  {'term':34} {'P / E (1/s)':>12}  {'share of sum |P|':>16}")
    ab = sum(abs(x) for x in P.values())
    for k, x in sorted(P.items(), key=lambda kv: -kv[1]):
        print(f"  {k:34} {x / E:12.2e}  {x / ab:16.0%}")
        if k == "horizontal advection":
            for s in list(Ps)[:2]:
                print(f"  {s:34} {Ps[s] / E:12.2e}")
        if k == "vertical advection":
            for s in list(Ps)[2:]:
                print(f"  {s:34} {Ps[s] / E:12.2e}")
    prof = np.sqrt((du ** 2 + dv ** 2).sum(axis=(1, 2)))
    print("\n  levels carrying the mode's energy: " +
          " ".join(f"L{k:02d} {prof[k] ** 2 / (prof ** 2).sum():.0%}" for k in np.argsort(prof)[::-1][:5]))


if __name__ == "__main__":
    main()

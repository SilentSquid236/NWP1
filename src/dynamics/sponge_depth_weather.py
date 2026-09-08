"""
Does a DEEP sponge suppress developing weather?

WHAT THE LAST EXPERIMENT SETTLED, AND WHAT IT DID NOT

The divergent-only sponge failed outright -- reflection 55 m/s against 21 for
a plain 8-level layer -- because a mountain wave is not purely divergent, and
damping only its divergent part leaves the rest to reflect off the lid.

But the same run overturned the premise. A deep PLAIN sponge does not flatten
the jet at all:

    sponge levels |  5   |  8   |  12
    jet drift 24h | 0.24%| 0.22%| 0.22%
    reflection    | 36.4 | 21.3 | 15.1

The reason the sponge had to stay shallow was inherited from P-16, where it
relaxed toward the HORIZONTAL MEAN. It now relaxes toward a frozen reference
state, and for a steady jet the jet IS the reference, so there is nothing to
damp. That constraint was assumed, never re-measured after the fix.

WHICH LEAVES THE REAL QUESTION, AND THE ONE A STEADY JET CANNOT ANSWER

A frozen reference is exactly right for a jet that does not change and exactly
wrong for one that does. Relaxing 12 of 20 levels toward a state frozen at
hour zero should suppress DEVELOPMENT -- the weather -- even though it leaves
a steady jet alone.

So the measurement is baroclinic growth rate against sponge depth, the same
diagnostic that rejected divergence damping (1.21x/day at 0.00, 0.33x at
0.01). A deep sponge is only usable if growth survives it.

PREDICTION, BEFORE THE RUN: growth falls with sponge depth, and 12 levels is
unusable despite having the least reflection.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma


def growth(sponge_levels, days=2):
    gr = CGrid(48, 48, 60e3, 60e3, f0=1.0e-4, beta=1.6e-11)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev, sponge_levels=sponge_levels)

    k_y = 2 * np.pi / gr.Ly
    m.pi = np.full((gr.ny, gr.nx), 101325.0) - lev.p_top
    p = lev.pressure(m.pi)
    T = 258.0 + 6.0 * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA

    phi = m.geopotential()
    for k in range(lev.nz):
        m.u[k] = -0.5 * (gr.dy_forward(phi[k])
                         + gr.dy_backward(phi[k])) / gr.f0
    jet = float(np.abs(m.u).max())

    seed = 0.5 * np.sin(4 * np.pi * gr.Xc / gr.Lx) * np.sin(k_y * gr.Yc)
    m.theta += seed

    def eddy(mm):
        up = mm.u - mm.u.mean(axis=2, keepdims=True)
        vp = mm.v - mm.v.mean(axis=2, keepdims=True)
        return float((up ** 2 + vp ** 2).sum())

    dt = m.max_dt()
    m.run(24 * 3600, dt=dt)
    e1 = eddy(m)
    m.run(24 * 3600, dt=dt)
    e2 = eddy(m)
    if not np.isfinite(m.u).all() or e1 <= 0:
        return float("nan"), float("nan"), jet
    g = e2 / e1
    return g, float(np.abs(m.v).max()), jet


if __name__ == "__main__":
    print("baroclinic development against sponge depth "
          "(48x48x20, 60 km, 2 days)\n")
    print(f"{'sponge':>7} {'eddy energy x/day':>18} {'max|v|':>9} {'verdict':>12}")
    for nsp in (0, 5, 8, 12):
        g, vmax, jet = growth(nsp)
        verdict = ("healthy" if g > 1.2 else
                   "suppressed" if np.isfinite(g) else "diverged")
        print(f"{nsp:7d} {g:18.2f} {vmax:9.1f} {verdict:>12}", flush=True)

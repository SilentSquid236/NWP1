"""
Is development SUPPRESSED, or has it already SATURATED?

The ratio used so far -- eddy energy on day 2 over day 1 -- is the same
diagnostic that rejected divergence damping, and it has a blind spot. A wave
that grew fast and saturated during day 1 reports a ratio near 1 and looks
identical to a wave that never grew.

The numbers that raised the doubt:

    sponge |  0   |  2   |  3   |  5
    ratio  | 3.18 | 0.93 | 0.76 | 0.85
    max|v| | 12.7 | 12.8 |  7.7 |  2.1

At two levels the ratio collapses but max|v| is UNCHANGED. Those two facts
cannot both mean suppression. At five levels max|v| is 2.1 and it clearly is.

So take the whole curve instead of two points on it: eddy energy every six
hours, which separates "never grew" from "grew and levelled off".
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA
from primitive_sigma import PrimitiveSigma


def curve(nsp, hours=48, every=6):
    gr = CGrid(48, 48, 60e3, 60e3, f0=1.0e-4, beta=1.6e-11)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev, sponge_levels=nsp)
    k_y = 2 * np.pi / gr.Ly
    m.pi = np.full((gr.ny, gr.nx), 101325.0) - lev.p_top
    p = lev.pressure(m.pi)
    m.theta = (258.0 + 6.0 * np.cos(k_y * gr.Yc)) / (p / P0) ** KAPPA
    phi = m.geopotential()
    for k in range(lev.nz):
        m.u[k] = -0.5 * (gr.dy_forward(phi[k])
                         + gr.dy_backward(phi[k])) / gr.f0
    m.theta += 0.5 * np.sin(4 * np.pi * gr.Xc / gr.Lx) * np.sin(k_y * gr.Yc)

    def eddy(mm):
        up = mm.u - mm.u.mean(axis=2, keepdims=True)
        vp = mm.v - mm.v.mean(axis=2, keepdims=True)
        return float((up ** 2 + vp ** 2).sum())

    dt = m.max_dt()
    out = []
    for _ in range(hours // every):
        m.run(every * 3600, dt=dt)
        if not np.isfinite(m.u).all():
            out.append(float("nan"))
            break
        out.append(eddy(m))
    return out


if __name__ == "__main__":
    hrs = list(range(6, 49, 6))
    print("eddy kinetic energy against time (arbitrary units)\n")
    print("sponge | " + " ".join(f"{h:>9}h" for h in hrs))
    for nsp in (0, 2, 3, 5):
        c = curve(nsp)
        c += [float("nan")] * (len(hrs) - len(c))
        print(f"{nsp:6d} | " + " ".join(f"{v:10.2e}" for v in c), flush=True)

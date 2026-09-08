"""
The production question: can the radiative lid replace the sponge?

Partial absorption would still be a win if it costs nothing. The sponge's
price was development -- every setting turned a growing baroclinic wave into a
decaying one (P-49). The radiative boundary damps nothing, so in principle it
should not have that price at all.

Two measurements, the same two that have decided every sponge question:

  1. baroclinic development, 48 h -- eddy kinetic energy every 6 h
  2. survival over 2500 m terrain, 12 h

A configuration that develops weather AND survives terrain replaces the
sponge. One that develops but dies over terrain does not.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from lid_test import build_on


def curve(sponge, radiative, hours=48, every=6):
    gr = CGrid(48, 48, 60e3, 60e3, f0=1.0e-4, beta=1.6e-11)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev, sponge_levels=sponge,
                       radiative_top=radiative)
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


def terrain(sponge, radiative, hours=12):
    lev = SigmaLevels(20)
    m = build_on(lev, 2500.0, sponge)
    m.radiative_top = radiative
    dt = m.max_dt()
    done = 0
    for _ in range(hours):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    return done, (float(np.abs(m.u).max()) if np.isfinite(m.u).all()
                  else float("nan"))


if __name__ == "__main__":
    hrs = list(range(6, 49, 6))
    configs = [("sponge 5, rigid", 5, False),
               ("no sponge, rigid", 0, False),
               ("no sponge, RADIATIVE", 0, True),
               ("sponge 5, RADIATIVE", 5, True)]

    print("baroclinic development (eddy kinetic energy)\n")
    print("config                | " + " ".join(f"{h:>8}h" for h in hrs)
          + " | 48h/6h")
    ratios = {}
    for name, sp, rad in configs:
        c = curve(sp, rad)
        c += [float("nan")] * (len(hrs) - len(c))
        r = c[-1] / c[0] if np.isfinite(c[-1]) and c[0] else float("nan")
        ratios[name] = r
        print(f"{name:21} | " + " ".join(f"{v:9.2e}" for v in c)
              + f" | {r:6.2f}", flush=True)

    print("\n2500 m terrain, 12 h\n")
    print(f"{'config':21} {'survived':>9} {'max|u|':>8}")
    for name, sp, rad in configs:
        d, u = terrain(sp, rad)
        print(f"{name:21} {d:6d}/12 {u:8.1f}", flush=True)

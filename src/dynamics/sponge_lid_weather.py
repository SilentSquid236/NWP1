"""
The sponge is not too strong. It is in the wrong place.

Everything tried so far leaves development suppressed at five levels:

    rate 15 min / 1 h / 6 h      -> 0.85 / 0.78 / 0.82 x per day
    frozen vs running reference  -> 0.85 / 0.80
    against NO sponge            -> 3.18

Neither the strength nor the target moves it, which rules out both and points
somewhere else. Baroclinic instability is a COUPLED mode: an upper-level wave
and a lower-level wave holding each other up. Damp either end hard enough and
the mode stops growing, whatever you damp it toward.

And with p_top = 200 hPa the lid IS the tropopause. The top five of twenty
sigma levels sit around 200-300 hPa -- precisely where the upper half of a
baroclinic wave lives. Real models put the sponge in the stratosphere, above
the weather. Ours has nowhere to put it that is not in the weather.

PREDICTION, BEFORE THE RUN: raising the lid moves the sponge off the
baroclinic wave and development recovers toward 3.18, with the sponge still
absorbing the mountain wave.

THIS REOPENS A DECISION. The note in `SigmaLevels.__init__` says a 200 hPa lid
beats 50 hPa, re-measured on corrected initial states as recently as
2026-09-03. That measurement was SURVIVAL over terrain, and survival is not
the only thing a lid affects. If development says the opposite, the model has
been trading weather for stability without anyone noticing.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA
from primitive_sigma import PrimitiveSigma
import lid_test


def baroclinic(nsp, p_top):
    gr = CGrid(48, 48, 60e3, 60e3, f0=1.0e-4, beta=1.6e-11)
    lev = SigmaLevels(20, p_top=p_top)
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
    m.run(24 * 3600, dt=dt); e1 = eddy(m)
    m.run(24 * 3600, dt=dt); e2 = eddy(m)
    if not np.isfinite(m.u).all() or e1 <= 0:
        return float("nan"), float("nan"), p[:, 0, 0]
    return e2 / e1, float(np.abs(m.v).max()), p[:, 0, 0]


def terrain(nsp, p_top, hours=6):
    m = lid_test.build_on(SigmaLevels(20, p_top=p_top), 2500.0, nsp)
    u0 = m.u.copy()
    dt = m.max_dt()
    done = 0
    for _ in range(hours):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            break
        done += 1
    if not np.isfinite(m.u).all():
        return float("nan"), done
    return float(np.abs(m.u - u0).max()), done


if __name__ == "__main__":
    print("development and reflection against LID HEIGHT, 5-level sponge\n")
    print(f"{'p_top':>8} {'sponge base':>12} {'eddy x/day':>11} "
          f"{'refl max|du|':>13} {'survived':>9}")
    for p_top in (20000.0, 10000.0, 5000.0):
        for nsp in (0, 5):
            g, vmax, p = baroclinic(nsp, p_top)
            refl, done = terrain(nsp, p_top)
            base = p[4] / 100 if nsp else float("nan")
            print(f"{p_top/100:7.0f}h {base:11.0f}h {g:11.2f} "
                  f"{refl:13.2f} {done:6d}/6"
                  + ("   (no sponge)" if nsp == 0 else ""), flush=True)

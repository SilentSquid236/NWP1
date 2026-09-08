"""
The sponge's problem is WHAT it relaxes toward, not how hard.

Rate ladder so far, at five levels:

    timescale   | 15 min | 1 h
    eddy x/day  |  0.85  | 0.78
    reflection  | 36.39  | 37.84

Weakening the rate does not help, and the reason is arithmetic: even a
one-hour relaxation applied for forty-eight hours is forty-eight e-foldings.
Anything short of "off" is complete damping over a forecast, so the rate is
not the knob.

What the sponge relaxes TOWARD is. Right now it is a state frozen at hour
zero, so it damps every deviation from the initial condition -- and a
developing baroclinic wave is exactly a deviation from the initial condition.
The sponge cannot tell weather from a gravity wave because it is not looking
at anything that distinguishes them.

THE DISTINGUISHING PROPERTY IS TIMESCALE. A vertically propagating mountain
wave oscillates on minutes; baroclinic development takes a day. So relax
toward a RUNNING LOW-PASS of the state rather than a frozen snapshot: the
reference follows the slow evolution and lags the fast oscillation, so the
deviation the sponge sees is the wave and not the weather.

    ref <- ref + (dt/tau_ref) (state - ref)

PREDICTION, BEFORE THE RUN, with tau_ref = 6 h:

  1. development recovers toward the 3.18x/day of no sponge at all
  2. reflection stays close to the frozen sponge's 36 m/s, not worse
  3. the thermal-wind jet still persists

If 1 improves and 2 collapses, the sponge has simply been turned off and the
right conclusion is that this does not work.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA
from primitive_sigma import PrimitiveSigma
import lid_test


def running_reference(m, tau_ref):
    """Replace the frozen reference with a low-pass of the evolving state."""
    m.set_reference()
    orig_step = m.step

    def step(dt):
        orig_step(dt)
        a = min(dt / tau_ref, 1.0)
        m._u_ref += a * (m.u - m._u_ref)
        m._v_ref += a * (m.v - m._v_ref)

    m.step = step
    return m


def baroclinic(nsp, tau_ref=None):
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
    if tau_ref:
        running_reference(m, tau_ref)

    def eddy(mm):
        up = mm.u - mm.u.mean(axis=2, keepdims=True)
        vp = mm.v - mm.v.mean(axis=2, keepdims=True)
        return float((up ** 2 + vp ** 2).sum())

    dt = m.max_dt()
    m.run(24 * 3600, dt=dt); e1 = eddy(m)
    m.run(24 * 3600, dt=dt); e2 = eddy(m)
    if not np.isfinite(m.u).all() or e1 <= 0:
        return float("nan"), float("nan")
    return e2 / e1, float(np.abs(m.v).max())


def reflection(nsp, tau_ref=None, hours=6):
    m = lid_test.build_on(SigmaLevels(20), 2500.0, nsp)
    if tau_ref:
        running_reference(m, tau_ref)
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
    print("frozen vs running sponge reference\n")
    print(f"{'levels':>7} {'reference':>18} {'eddy x/day':>11} "
          f"{'refl max|du|':>13} {'survived':>9}")
    for nsp, tau, label in ((0, None, "no sponge"),
                            (5, None, "frozen"),
                            (5, 6 * 3600.0, "running, 6 h"),
                            (8, 6 * 3600.0, "running, 6 h"),
                            (5, 1 * 3600.0, "running, 1 h")):
        g, vmax = baroclinic(nsp, tau)
        refl, done = reflection(nsp, tau)
        print(f"{nsp:7d} {label:>18} {g:11.2f} {refl:13.2f} {done:6d}/6",
              flush=True)

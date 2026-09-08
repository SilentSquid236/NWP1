"""
Three levels survives 2500 m clean terrain but NOT the decisive test.

Reducing the default to 3 took `test_primitive_sigma.py` from 12/12 to 11/12
on the noisy 400 m case. Two different cases want two different minima, so the
production default has to be the minimum that satisfies BOTH, not whichever
one was measured last.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import (SigmaLevels, RD, G0, P0, KAPPA,
                   hydrostatic_geopotential, pressure_gradient_force)
from primitive_sigma import PrimitiveSigma
from subgrid import balance_initial_state
from initialization import filter_initial_state
from lid_test import build_on


def decisive(nsp):
    """The decisive test's own case: 400 m terrain, 1.2 m/s noise, filtered."""
    gr = CGrid(90, 88, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    rng = np.random.default_rng(0)
    h = 400.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                         ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h, sponge_levels=nsp)
    k_y = 2 * np.pi / gr.Ly
    m.pi = 101325.0 * np.exp(-G0 * h / (RD * 280.0)) - lev.p_top
    p = lev.pressure(m.pi)
    m.theta = (288.0 - 55.0 * (1 - p / p.max())
               - 1.5 * np.cos(k_y * gr.Yc)) / (p / P0) ** KAPPA
    phi = hydrostatic_geopotential(m.theta, m.pi, lev, phi_surface=m.phi_s)
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u[:] = -fy / gr.f0
    m.v[:] = fx / gr.f0
    m.u += rng.normal(0, 1.2, m.u.shape)
    m.v += rng.normal(0, 1.2, m.v.shape)
    m.u, m.v, m.theta = filter_initial_state(m.u, m.v, m.theta, gr)
    m.u, m.v, _ = balance_initial_state(m.u, m.v, gr, verbose=False)
    dt = m.max_dt()
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            break
        done += 1
    return done, (float(np.abs(m.u).max()) if np.isfinite(m.u).all()
                  else float("nan"))


def tall(nsp):
    m = build_on(SigmaLevels(20), 2500.0, nsp)
    dt = m.max_dt()
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    return done, (float(np.abs(m.u).max()) if np.isfinite(m.u).all()
                  else float("nan"))


if __name__ == "__main__":
    print("smallest sponge that survives BOTH production cases\n")
    print(f"{'sponge':>7} {'400 m + noise':>15} {'2500 m clean':>14}")
    for nsp in (3, 4, 5):
        a, ua = decisive(nsp)
        b, ub = tall(nsp)
        print(f"{nsp:7d} {a:11d}/12    {b:10d}/12", flush=True)

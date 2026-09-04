"""
Boundary-layer parameterization sweep.

QUESTION BEING ASKED
    Richardson-number vertical mixing raised 12-hour survival from 2-3/12 to
    6/12 on the hard cases (grid-scale noise, tall terrain) but no further.
    Surface drag is the other half of a boundary-layer scheme. Does adding it
    move the survival counts, and is the effect additive with mixing or
    redundant with it?

WHY A MATRIX AND NOT A PATCH
    The failure mode has already survived nine single-candidate patches (see
    docs/STABILITY.md). Turning both knobs independently over the same
    terrain/noise grid is the only way to attribute a change to a scheme
    rather than to a coincidence of settings.

Run:  python sweep_boundary_layer.py [--quick]
"""

import sys
import time

import numpy as np

from grid import CGrid
from sigma import (SigmaLevels, RD, G0, P0, KAPPA,
                   hydrostatic_geopotential, pressure_gradient_force)
from primitive_sigma import PrimitiveSigma
from subgrid import balance_initial_state

TERRAINS = [0.0, 1000.0, 2500.0]
NOISES = [0.0, 1.2]
HOURS = 12


def build(terrain_height, noise, mixing, drag, seed=0, dT=1.5):
    """
    The decisive test's initial state, corrected.

    TWO INITIALIZATION ERRORS WERE FIXED HERE (see docs/RESEARCH_LOG.md).

    1. The meridional temperature contrast was 6 K, which by thermal wind
       implies a 166 m/s jet. The old test then CLIPPED the wind at 60 m/s,
       destroying geostrophic balance over 34% of the domain and injecting a
       sharp, scale-selective imbalance. The measured grid-scale (2dx) growth
       in the clean case came from that clip, not from the dynamics. 1.5 K
       gives a 41 m/s jet -- a realistic Northeast winter jet -- with no clip.

    2. The balanced wind was taken as -d(phi)/dy / f. On sigma surfaces the
       horizontal force has TWO terms that largely cancel over sloping
       ground; keeping only the first implies an 845 m/s "balanced" wind over
       2500 m terrain. The wind now comes from the full PGF, giving 41 m/s
       with a physical 15 m/s cross-mountain ageostrophic component.

    Every terrain row of the earlier mixing baseline was measured against
    initial states carrying these artifacts and should be disregarded.
    """
    gr = CGrid(90, 88, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    rng = np.random.default_rng(seed)

    h = terrain_height * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                                  ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h, mixing=mixing, drag=drag)

    k_y = 2 * np.pi / gr.Ly
    p_s = 101325.0 * np.exp(-G0 * h / (RD * 280.0))
    m.pi = p_s - lev.p_top
    p = lev.pressure(m.pi)
    T = 288.0 - 55.0 * (1 - p / p.max()) - dT * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA

    phi = hydrostatic_geopotential(m.theta, m.pi, lev, phi_surface=m.phi_s)
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u[:] = -fy / gr.f0
    m.v[:] = fx / gr.f0

    if noise > 0:
        m.u += rng.normal(0, noise, m.u.shape)
        m.v += rng.normal(0, noise, m.v.shape)
        m.u, m.v, _ = balance_initial_state(m.u, m.v, gr, verbose=False)

    return m


def survive(m, hours=HOURS):
    """Hours completed before the run stops being finite or physical."""
    dt = m.max_dt()
    done = 0
    for _ in range(hours):
        try:
            m.run(3600, dt=dt)
        except Exception:
            break
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150.0:
            break
        done += 1
    umax = float(np.abs(m.u).max()) if np.isfinite(m.u).all() else float("nan")
    return done, umax, float(dt)


def main():
    quick = "--quick" in sys.argv
    hours = 3 if quick else HOURS
    rows = []
    print(f"\nBoundary-layer sweep  ({hours} h ceiling, 90x88x20)\n" + "=" * 78)
    print(f"{'terrain':>8} {'noise':>6} {'mixing':>7} {'drag':>5} "
          f"{'survived':>9} {'max|u|':>8} {'dt':>6} {'wall':>7}")
    for hgt in TERRAINS:
        for noise in NOISES:
            for mixing in (False, True):
                for drag in (False, True):
                    t0 = time.time()
                    m = build(hgt, noise, mixing, drag)
                    done, umax, dt = survive(m, hours)
                    wall = time.time() - t0
                    rows.append((hgt, noise, mixing, drag, done, umax))
                    print(f"{hgt:8.0f} {noise:6.1f} {str(mixing):>7} "
                          f"{str(drag):>5} {done:6d}/{hours} {umax:8.1f} "
                          f"{dt:6.1f} {wall:6.1f}s", flush=True)
    print("=" * 78)

    # Marginal effect of each scheme, averaged over the terrain/noise grid.
    def mean_where(mix, drg):
        v = [r[4] for r in rows if r[2] is mix and r[3] is drg]
        return sum(v) / len(v)

    print("\nmean hours survived (averaged over terrain x noise):")
    print(f"  neither        {mean_where(False, False):.2f}")
    print(f"  mixing only    {mean_where(True,  False):.2f}")
    print(f"  drag only      {mean_where(False, True):.2f}")
    print(f"  both           {mean_where(True,  True):.2f}")
    np.save("sweep_boundary_layer.npy", np.array(
        [(h, n, float(mx), float(dg), d, u) for h, n, mx, dg, d, u in rows]))


if __name__ == "__main__":
    main()

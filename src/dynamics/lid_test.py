"""
Is the model lid too low for a mountain wave?

MEASURED SO FAR (sponge_edge_test.py, 2500 m terrain, clean, filtered):

    sponge levels |  0   |  5   |  8   | 12
    peak growth k |  0   |  5   |  8   | 18
    max|du| (m/s) | 60.8 | 36.4 | 21.3 | 15.1

The growth peak tracks the SPONGE BASE exactly. That is partial reflection off
the absorbing layer's lower edge, not a physical instability. Deepening the
sponge suppresses it, but a 12-level sponge on a 20-level model is damping the
free troposphere -- the same mistake as the sponge that relaxed toward the
horizontal mean and flattened the jet.

The absorbing layer has to be deep compared with the wave it absorbs. This
script measures the mountain wave's vertical wavelength, lz = 2*pi*U/N, and
then tests whether raising the lid -- which buys sponge depth without spending
troposphere -- removes the reflection.

A PRIOR NEGATIVE RESULT IS BEING RE-OPENED. `SigmaLevels.__init__` records
that raising the lid from 200 to 50 hPa cost about an hour of stability. That
measurement was made on initial states now known to carry a clipped
166 m/s jet and a one-term geostrophic wind over terrain, so it does not
survive as evidence and is being re-run.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels, RD, G0, P0, KAPPA
from turbulence import richardson
from initialization import filter_initial_state
from subgrid import balance_initial_state
import probe_failure as P


def wave_scales(m):
    """N, U and the hydrostatic mountain wave's vertical wavelength."""
    _, N2, _, dz = richardson(m.u, m.v, m.theta, m.pi, m.lev)
    N = float(np.sqrt(np.maximum(np.nanmean(N2[:15]), 1e-8)))
    U = float(np.abs(m.u).mean())
    depth = float(np.nansum(np.nanmean(dz, axis=(1, 2))))
    return N, U, 2 * np.pi * U / N, depth


def prep(p_top, sponge_levels, nz=20, hgt=2500.0):
    lev = SigmaLevels(nz, p_top=p_top)
    m = P.build(hgt, 0.0, dT=1.5, clip=None, sponge_levels=sponge_levels)
    # rebuild on the requested vertical grid
    m = P.build(hgt, 0.0, dT=1.5, clip=None, sponge_levels=sponge_levels,
                levels=lev) if False else m
    return m, lev


def build_on(lev, hgt, sponge_levels):
    """Same initial state as probe_failure.build, on a chosen vertical grid."""
    from grid import CGrid
    from primitive_sigma import PrimitiveSigma
    from sigma import hydrostatic_geopotential, pressure_gradient_force
    gr = CGrid(90, 88, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                       ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h, sponge_levels=sponge_levels)
    k_y = 2 * np.pi / gr.Ly
    m.pi = 101325.0 * np.exp(-G0 * h / (RD * 280.0)) - lev.p_top
    p = lev.pressure(m.pi)
    T = 288.0 - 55.0 * (1 - p / p.max()) - 1.5 * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA
    phi = hydrostatic_geopotential(m.theta, m.pi, lev, phi_surface=m.phi_s)
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u[:] = -fy / gr.f0
    m.v[:] = fx / gr.f0
    m.u, m.v, m.theta = filter_initial_state(m.u, m.v, m.theta, gr)
    m.u, m.v, _ = balance_initial_state(m.u, m.v, gr, verbose=False)
    return m


if __name__ == "__main__":
    print("2500 m terrain, clean, filtered, mixing+drag on, 12 h ceiling\n")
    lev0 = SigmaLevels(20)
    m0 = build_on(lev0, 2500.0, 5)
    N, U, lz, depth = wave_scales(m0)
    print(f"wave scales at 200 hPa lid: N={N:.4f} 1/s, U={U:.1f} m/s, "
          f"lambda_z={lz/1000:.1f} km, model depth={depth/1000:.1f} km")
    print(f"  -> an absorbing layer should be a large fraction of "
          f"{lz/1000:.0f} km; five levels here is "
          f"{depth*5/20/1000:.1f} km\n")

    print(f"{'p_top':>8} {'nz':>4} {'sponge':>7} {'sponge km':>10} "
          f"{'survived':>9} {'max|u|':>8} {'peak k':>7}")
    for p_top in (20000.0, 10000.0, 5000.0):
        for nsp in (5, 8):
            lev = SigmaLevels(20, p_top=p_top)
            m = build_on(lev, 2500.0, nsp)
            _, _, _, dep = wave_scales(m)
            u0 = m.u.copy()
            dt = m.max_dt()
            done = 0
            for _ in range(12):
                m.run(3600, dt=dt)
                if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
                    break
                done += 1
            fin = np.isfinite(m.u).all()
            if fin:
                amp = np.sqrt(((m.u - u0) ** 2).mean(axis=(1, 2)))
                peak = int(np.argmax(amp))
                umax = float(np.abs(m.u).max())
            else:
                peak, umax = -1, float("nan")
            print(f"{p_top/100:7.0f}h {20:4d} {nsp:7d} {dep*nsp/20/1000:10.1f} "
                  f"{done:6d}/12 {umax:8.1f} {peak:7d}", flush=True)

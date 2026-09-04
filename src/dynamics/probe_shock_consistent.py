"""
The synthetic analysis itself was not hydrostatically self-consistent.

`probe_shock.build` perturbed temperature by -1.5*cos(k_y*y) K and, separately,
perturbed the geopotential heights by -45*cos(k_y*y) m. Those two numbers were
chosen independently, so they do not satisfy hydrostatic balance with each
other. A real analysis does. Everything measured from that state -- the 140 m
geopotential spread, the 2.7 m/s/h acceleration, the 100 m/s "balanced" wind --
was the model responding to an inconsistency in the test data.

This builds the heights BY INTEGRATING the temperature field, so the synthetic
analysis has the property a real one has, and re-runs the same measurements.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import (SigmaLevels, hydrostatic_geopotential,
                   pressure_gradient_force, RD, G0, P0, KAPPA)
from primitive_sigma import PrimitiveSigma
from interpolate import pressure_to_sigma, interp_log_p

LEVELS_HPA = [1000, 975, 950, 925, 900, 850, 800, 750, 700, 650,
              600, 550, 500, 450, 400, 350, 300, 275, 250, 200]
P_PA = np.array(LEVELS_HPA, dtype=float) * 100.0


def consistent_analysis(gr, dT=1.5, p_sea=101325.0):
    """
    A synthetic analysis whose heights are the hydrostatic integral of its own
    temperatures -- the property that makes it usable as a test of anything
    else.

        dz = -(R T / g) d(ln p),  integrated upward from 1000 hPa.
    """
    ny, nx = gr.ny, gr.nx
    T0, L = 288.15, 0.0065
    k_y = 2 * np.pi / gr.Ly

    T = np.empty((P_PA.size, ny, nx))
    for k, p in enumerate(P_PA):
        T_std = T0 * (p / 101325.0) ** (RD * L / G0)
        T[k] = T_std - dT * np.cos(k_y * gr.Yc)

    # Height of the 1000 hPa surface, hydrostatically below sea level pressure.
    z = np.empty_like(T)
    z[0] = (RD * T[0] / G0) * np.log(p_sea / P_PA[0])
    for k in range(1, P_PA.size):
        T_layer = 0.5 * (T[k - 1] + T[k])
        z[k] = z[k - 1] + (RD * T_layer / G0) * np.log(P_PA[k - 1] / P_PA[k])
    return T, z


def report(label, gr, lev, h, pi, u, v, th, hours=6):
    m = PrimitiveSigma(gr, lev, terrain=h)
    m.pi, m.u, m.v, m.theta = pi.copy(), u.copy(), v.copy(), th.copy()
    u0, v0 = m.u.copy(), m.v.copy()
    du, dv, _, _ = m.tendencies(m.u, m.v, m.theta, m.pi)
    m.run(hours * 3600, dt=m.max_dt())
    ok = np.isfinite(m.u).all()
    drift = float(max(np.abs(m.u - u0).max(), np.abs(m.v - v0).max())) \
        if ok else float("nan")
    print(f"  {label:<24} |du/dt| {np.abs(du).max()*3600:7.3f}  "
          f"|dv/dt| {np.abs(dv).max()*3600:7.3f} m/s/h   "
          f"wind change over {hours} h {drift:7.3f} m/s", flush=True)


if __name__ == "__main__":
    for hgt in (0.0, 1500.0, 2500.0):
        gr = CGrid(60, 60, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
                   edge_mode="replicate")
        lev = SigmaLevels(20)
        h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 150e3) ** 2 +
                           ((gr.Yc - gr.Ly / 2) / 150e3) ** 2))
        T, z = consistent_analysis(gr)
        zero = np.zeros_like(T)
        pi, u, v, th = pressure_to_sigma(zero, zero, T, z, P_PA, h, lev)

        phi = hydrostatic_geopotential(th, pi, lev, phi_surface=G0 * h)
        p_dst = lev.pressure(pi)
        phi_a = G0 * interp_log_p(z, P_PA, p_dst, extrapolate="linear")
        d = (phi - phi_a) / G0
        spread = float((d.max(axis=(1, 2)) - d.min(axis=(1, 2))).max())

        fx, fy = pressure_gradient_force(phi, th, pi, lev, gr)
        ub, vb = -fy / gr.f0, fx / gr.f0

        print(f"\nterrain {hgt:.0f} m   "
              f"geopotential error vs analysis: {spread:.2f} m spread   "
              f"implied balanced wind {np.abs(ub).max():.1f} m/s")
        report("at rest", gr, lev, h, pi, u, v, th)
        report("balanced", gr, lev, h, pi, ub, vb, th)

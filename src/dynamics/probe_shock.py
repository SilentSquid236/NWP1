"""
Where does the 9 m/s of initialization shock come from?

A converted analysis started at rest over a 1500 m mountain develops 9.1 m/s
of spurious wind in six hours. A state hydrostatically consistent with the
model's OWN discretisation develops 0.004 m/s over 4000 m terrain, so the
coordinate is not the cause (P-14's confirming measurement).

Hypothesis to test, stated first: the interpolated theta reproduces the
analysis's stratification but not the model's discrete hydrostatic integral,
so the model's geopotential differs from the analysis's geopotential by an
amount that VARIES HORIZONTALLY -- and a horizontally varying geopotential
error is a pressure-gradient force with nothing balancing it.

The measurement that decides it: compare the model's own Phi against the
analysis heights interpolated to the same pressures. If the difference is
horizontally uniform it is harmless; if it varies, its gradient is the shock.
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


def us_standard(p):
    T0, L = 288.15, 0.0065
    T = T0 * (p / 101325.0) ** (RD * L / G0)
    return T, (T0 - T) / L


def build(hgt=1500.0, ny=60, nx=60):
    gr = CGrid(nx, ny, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 150e3) ** 2 +
                       ((gr.Yc - gr.Ly / 2) / 150e3) ** 2))
    T1, z1 = us_standard(P_PA)
    T = np.repeat(np.repeat(T1[:, None, None], ny, 1), nx, 2)
    z = np.repeat(np.repeat(z1[:, None, None], ny, 1), nx, 2)
    k_y = 2 * np.pi / gr.Ly
    T = T - 1.5 * np.cos(k_y * gr.Yc)
    z = z - 1.5 * np.cos(k_y * gr.Yc) * 30.0
    u = np.zeros_like(T)
    return gr, lev, h, T, z, u


if __name__ == "__main__":
    gr, lev, h, T, z, u = build()
    pi, us_, vs_, th = pressure_to_sigma(u, u, T, z, P_PA, h, lev)

    m = PrimitiveSigma(gr, lev, terrain=h)
    m.pi, m.u, m.v, m.theta = pi, us_, vs_, th

    # The model's own geopotential, and the analysis's, at the same pressures.
    phi_model = hydrostatic_geopotential(th, pi, lev, phi_surface=G0 * h)
    p_dst = lev.pressure(pi)
    phi_analysis = G0 * interp_log_p(z, P_PA, p_dst, extrapolate="linear")

    d = phi_model - phi_analysis
    print("geopotential: model minus analysis, by level (0 = lid)")
    print(f"{'k':>3} {'p (hPa)':>9} {'mean (m)':>10} {'spread (m)':>11} "
          f"{'mean/g':>8}")
    for k in (0, 4, 8, 12, 16, 19):
        dk = d[k] / G0
        print(f"{k:3d} {p_dst[k].mean()/100:9.0f} {dk.mean():10.2f} "
              f"{dk.max()-dk.min():11.2f} {dk.mean():8.2f}")

    # The part that matters: the HORIZONTAL VARIATION, which is what a
    # pressure-gradient force sees. A uniform offset accelerates nothing.
    spread = (d / G0).max(axis=(1, 2)) - (d / G0).min(axis=(1, 2))
    print(f"\nmax horizontal spread of the error: {spread.max():.2f} m "
          f"at level {int(np.argmax(spread))}")

    # And the acceleration it implies, at rest.
    fx, fy = pressure_gradient_force(phi_model, th, pi, lev, gr)
    du, dv, dth, dpi = m.tendencies(m.u, m.v, m.theta, m.pi)
    print(f"initial acceleration at rest: max|du/dt| {np.abs(du).max():.2e} "
          f"m/s^2 = {np.abs(du).max()*3600:.1f} m/s per hour")

    # Control: a state built hydrostatically BY THE MODEL over the same
    # terrain, which is what P-14 measured at 0.004 m/s.
    m2 = PrimitiveSigma(gr, lev, terrain=h)
    T0 = 260.0
    ps2 = 101325.0 * np.exp(-G0 * h / (RD * T0))
    m2.pi = ps2 - lev.p_top
    m2.theta = T0 / (lev.pressure(m2.pi) / P0) ** KAPPA
    du2, _, _, _ = m2.tendencies(m2.u, m2.v, m2.theta, m2.pi)
    print(f"control (model-built isothermal state): "
          f"max|du/dt| {np.abs(du2).max():.2e} m/s^2 = "
          f"{np.abs(du2).max()*3600:.3f} m/s per hour")

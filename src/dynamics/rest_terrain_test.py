"""
DISCRIMINATOR: is the tall-terrain failure the mountain WAVE, or the
COORDINATE itself?

A motionless, hydrostatically balanced atmosphere over a mountain has no wave,
no shear, no jet -- nothing to go unstable. Anything that grows is the
sigma-coordinate pressure-gradient truncation error, which scales with terrain
slope and has no physical source to blame.

`test_sigma.py` already checks this at 1200 m for 6 hours. This pushes it to
the heights and durations where the model actually fails.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma

print("motionless isothermal atmosphere over a mountain, 12 h")
print(f"{'height':>8} {'max slope':>10} {'max|u| 6h':>11} {'max|u| 12h':>11} "
      f"{'dp_s 12h':>10}")
for hgt in (1200.0, 2500.0, 4000.0):
    gr = CGrid(90, 88, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                       ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h)
    T0 = 260.0
    p_s = 101325.0 * np.exp(-G0 * h / (RD * T0))
    m.pi = p_s - lev.p_top
    m.theta = T0 / (lev.pressure(m.pi) / P0) ** KAPPA
    pi0 = m.pi.copy()
    slope = float(np.abs(np.diff(h, axis=0)).max() / gr.dy)

    dt = m.max_dt()
    m.run(6 * 3600, dt=dt)
    u6 = float(np.abs(m.u).max()) if np.isfinite(m.u).all() else float("nan")
    m.run(6 * 3600, dt=dt)
    ok = np.isfinite(m.u).all()
    u12 = float(np.abs(m.u).max()) if ok else float("nan")
    dps = float(np.abs(m.pi - pi0).max()) if ok else float("nan")
    print(f"{hgt:8.0f} {slope:10.5f} {u6:11.3f} {u12:11.3f} {dps:10.1f}",
          flush=True)

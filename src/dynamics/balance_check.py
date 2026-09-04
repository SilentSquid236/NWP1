"""
Is the tall-terrain INITIAL STATE balanced, or is it another badly posed test?

The 4000 m case starts with max|v| = 36 m/s of cross-mountain flow, against
20 m/s at 2500 m and 9 m/s at 1000 m. That is derived from the sigma PGF, so
it is geostrophic by construction -- but geostrophic balance is not the same
as a state the model will hold, and the clipped-jet episode is a standing
reminder to check rather than assume.

The direct test: evaluate the tendencies at t = 0. A balanced state has small
residual tendencies; an unbalanced one is already accelerating before a single
step is taken.

Also reported: the nondimensional mountain height Nh/U, which says which
dynamical regime each case is in. Below ~1 the flow goes over the mountain and
the wave is linear; above ~1 the flow is partly blocked and the wave breaks.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from turbulence import richardson
from lid_test import build_on

print(f"{'terrain':>8} {'max|u0|':>8} {'max|v0|':>8} {'|du/dt|':>10} "
      f"{'|dv/dt|':>10} {'accel m/s/h':>12} {'N h / U':>9}")
for hgt in (0.0, 1000.0, 2500.0, 4000.0):
    m = build_on(SigmaLevels(20), hgt, 8)
    du, dv, dth, dpi = m.tendencies(m.u, m.v, m.theta, m.pi)
    _, N2, _, _ = richardson(m.u, m.v, m.theta, m.pi, m.lev)
    N = float(np.sqrt(np.maximum(np.nanmean(N2[:15]), 1e-8)))
    U = float(np.abs(m.u).max())
    a = float(max(np.abs(du).max(), np.abs(dv).max()))
    print(f"{hgt:8.0f} {np.abs(m.u).max():8.2f} {np.abs(m.v).max():8.2f} "
          f"{np.abs(du).max():10.2e} {np.abs(dv).max():10.2e} "
          f"{a*3600:12.1f} {N*hgt/max(U,1e-9):9.2f}", flush=True)

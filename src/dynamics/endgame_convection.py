"""
PREDICTION, WRITTEN BEFORE THE RUN.

Over 2500 m terrain the run dies at hour 12 with the minimum Richardson number
going negative at hour 10 (static instability: the mountain wave overturns).
With dry convective adjustment the prediction is:

  1. the minimum Ri floors near 0 instead of going negative
  2. the count of statically unstable interfaces stops growing
  3. the run completes 12 hours
  4. the wind stays near 41 m/s -- convection removes the overturning, it does
     not remove the wave, so this is NOT expected to look like extra damping

Failure of (4) while (3) succeeds would mean the scheme is buying stability by
destroying the flow, which is the failure mode the sponge already had once.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from turbulence import richardson
from lid_test import build_on

for conv in (False, True):
    m = build_on(SigmaLevels(20), 2500.0, 8)
    m.convection = conv
    u0 = m.u.copy()
    dt = m.max_dt()
    print(f"\n2500 m, clean, filtered, sponge=8, convection={conv}")
    print(f"{'hour':>5} {'max|u|':>8} {'min Ri':>8} {'Ri<0':>6} "
          f"{'mixed %':>8} {'sweeps':>7}")
    for hr in range(1, 19):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            print(f"{hr:5d}   non-finite")
            break
        Ri, _, _, _ = richardson(m.u, m.v, m.theta, m.pi, m.lev)
        info = m._conv_info or {}
        print(f"{hr:5d} {np.abs(m.u).max():8.2f} {np.nanmin(Ri[:14]):8.3f} "
              f"{int((Ri[:14] < 0).sum()):6d} "
              f"{info.get('unstable_before', 0)*100:7.2f}% "
              f"{info.get('sweeps', 0):7d}", flush=True)

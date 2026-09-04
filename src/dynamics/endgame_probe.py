"""
What finishes the 2500 m run at hour 11-12?

Already excluded by measurement:
  * the sigma coordinate      -- a motionless atmosphere over 4000 m terrain
                                 drifts 0.009 m/s in 12 h, error linear in slope
  * initial-state artifacts   -- state is filtered and rebalanced
  * sponge reflection alone   -- deepening the sponge moves the growth peak and
                                 halves its amplitude, but survival stays 11/12
  * lid height                -- raising it to 100 or 50 hPa is neutral to worse

So the last hours are watched directly, hour by hour.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from turbulence import richardson
from initialization import filter_initial_state
from subgrid import balance_initial_state
from lid_test import build_on

m = build_on(SigmaLevels(20), 2500.0, 8)
u0 = m.u.copy()
dt = m.max_dt()
print(f"2500 m, clean, filtered, sponge=8, dt={dt:.2f}s\n")
print(f"{'hour':>5} {'max|u|':>8} {'max|v|':>8} {'max|sd|':>9} {'min Ri':>8} "
      f"{'Ri<.25':>7} {'peak k':>7} {'dom wl':>8} {'min p_s':>9}")
for hr in range(1, 13):
    m.run(3600, dt=dt)
    if not np.isfinite(m.u).all():
        print(f"{hr:5d}   non-finite")
        break
    d = m.u - u0
    amp = np.sqrt((d ** 2).mean(axis=(1, 2)))
    Ri, _, _, _ = richardson(m.u, m.v, m.theta, m.pi, m.lev)
    k = np.abs(np.fft.rfft(d - d.mean(axis=2, keepdims=True), axis=2))
    kamp = np.sqrt((k ** 2).mean(axis=(0, 1)))
    dom = d.shape[2] / max(int(np.argmax(kamp[1:])) + 1, 1)
    print(f"{hr:5d} {np.abs(m.u).max():8.2f} {np.abs(m.v).max():8.2f} "
          f"{np.abs(m.sigma_dot()).max():9.2e} {np.nanmin(Ri[:14]):8.3f} "
          f"{int((Ri[:14] < 0.25).sum()):7d} {int(np.argmax(amp)):7d} "
          f"{dom:7.1f}d {np.nanmin(m.surface_pressure):9.0f}", flush=True)

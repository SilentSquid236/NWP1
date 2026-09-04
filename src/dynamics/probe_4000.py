"""
What fails at 4000 m that does not fail at 2500 m?

Convective adjustment took 2500 m from 11/12 to 12/12 and out to 16 hours.
At 4000 m it changes nothing: 6/12 either way. So the mode is different, and
it is watched the same way -- hour by hour, every quantity that could be the
one running away.

Ruled out already for tall terrain generally:
  * the sigma coordinate  -- a motionless atmosphere over 4000 m drifts
                             0.009 m/s in 12 h, error linear in slope
  * initial-state artifacts, lid height, sponge depth
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from turbulence import richardson
from lid_test import build_on

for hgt in (2500.0, 4000.0):
    m = build_on(SigmaLevels(20), hgt, 8)
    u0 = m.u.copy()
    dt = m.max_dt()
    print(f"\n{hgt:.0f} m terrain, clean, filtered, sponge=8, "
          f"convection on, dt={dt:.2f}s")
    print(f"{'hour':>5} {'max|u|':>8} {'max|v|':>8} {'max|sd|':>9} "
          f"{'min Ri':>8} {'peak k':>7} {'dom wl':>8} {'min p_s':>9} "
          f"{'conv %':>7}")
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
        info = m._conv_info or {}
        print(f"{hr:5d} {np.abs(m.u).max():8.2f} {np.abs(m.v).max():8.2f} "
              f"{np.abs(m.sigma_dot()).max():9.2e} "
              f"{np.nanmin(Ri[:14]):8.3f} {int(np.argmax(amp)):7d} "
              f"{dom:7.1f}d {np.nanmin(m.surface_pressure):9.0f} "
              f"{info.get('unstable_before',0)*100:6.2f}%", flush=True)

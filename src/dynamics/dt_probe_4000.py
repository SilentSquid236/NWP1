"""Same question, instrumented: dt/2 at 4000 m, printing every hour."""
import numpy as np, time
np.seterr(all="ignore")
from sigma import SigmaLevels
from turbulence import richardson
from lid_test import build_on

m = build_on(SigmaLevels(20), 4000.0, 8)
dt = m.max_dt() * 0.5
print(f"4000 m, dt={dt:.2f}s (half the CFL estimate)")
print(f"{'hour':>5} {'max|u|':>8} {'max|v|':>8} {'max|sd|':>9} {'min Ri':>8} "
      f"{'conv %':>7} {'sweeps':>7} {'wall':>7}")
for hr in range(1, 13):
    t0 = time.time()
    m.run(3600, dt=dt)
    if not np.isfinite(m.u).all():
        print(f"{hr:5d}   non-finite")
        break
    Ri, _, _, _ = richardson(m.u, m.v, m.theta, m.pi, m.lev)
    info = m._conv_info or {}
    print(f"{hr:5d} {np.abs(m.u).max():8.2f} {np.abs(m.v).max():8.2f} "
          f"{np.abs(m.sigma_dot()).max():9.2e} {np.nanmin(Ri[:14]):8.3f} "
          f"{info.get('unstable_before',0)*100:6.2f}% "
          f"{info.get('sweeps',0):7d} {time.time()-t0:6.0f}s", flush=True)

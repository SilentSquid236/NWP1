"""
Is the eddy-diffusivity CEILING the binding constraint at 4000 m?

WHAT THE MEASUREMENTS SAY SO FAR

  * Nh/U -- the nondimensional mountain height -- is 0.38 at 1000 m, 0.96 at
    2500 m and 1.19 at 4000 m. The model is clean below 1, needs convective
    adjustment at 1, and fails above it. That is the boundary between a
    mountain wave that propagates and one that breaks.
  * The 4000 m trajectory is IDENTICAL at dt and dt/2 through hour 6 (max|u|
    53.86 both, min Ri 8.190 both, max|sigma_dot| 3.85e-05 both). The solution
    is converged in time; the timestep is not the cause.
  * Convective adjustment holds min Ri at 0.007 instead of letting it go
    negative, and the fraction of overturning interfaces climbs steadily
    (0.09% -> 0.37% by hour 6) until the adjustment is firing domain-wide.

The adjustment removes static instability but conserves momentum within the
column. What a breaking wave also does is generate turbulence, and
`turbulence.K_MAX` caps the eddy diffusivity at 100 m^2/s. Observed values in
a breaking mountain wave are 10^2 to 10^3 m^2/s, so the cap -- not the physics
-- may be limiting the dissipation.

A ceiling is a parameter, not a scheme, so this is a sensitivity measurement.
If survival is flat in K_MAX the ceiling is innocent and the missing piece is
elsewhere.
"""
import numpy as np
np.seterr(all="ignore")

import turbulence
from sigma import SigmaLevels
from lid_test import build_on

print("4000 m terrain, clean, filtered, sponge=8, convection on, 12 h")
print(f"{'K_MAX':>8} {'survived':>9} {'max|u|':>8} {'min Ri':>8}")
for kmax in (100.0, 300.0, 1000.0):
    turbulence.K_MAX = kmax
    m = build_on(SigmaLevels(20), 4000.0, 8)
    dt = m.max_dt()
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    fin = np.isfinite(m.u).all()
    if fin:
        Ri, _, _, _ = turbulence.richardson(m.u, m.v, m.theta, m.pi, m.lev)
        mr = float(np.nanmin(Ri[:14]))
    else:
        mr = float("nan")
    print(f"{kmax:8.0f} {done:6d}/12 "
          f"{(np.abs(m.u).max() if fin else float('nan')):8.1f} {mr:8.3f}",
          flush=True)

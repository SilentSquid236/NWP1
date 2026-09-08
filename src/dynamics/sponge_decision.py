"""
THE DECISION-RELEVANT MEASUREMENT.

Everything measured so far says the sponge cannot separate a mountain wave
from a baroclinic wave, because in this model they occupy the same levels:

    absorbing well (reflection 21-36)  ->  development 0.85-0.96 x/day
    development fine (1.89-2.45)       ->  reflection 50-60, i.e. not absorbing

So the question stops being "how do we tune it" and becomes "is it worth its
cost at all". The sponge exists to stop the run dying over terrain. That was
true before convective adjustment existed. It may not be true now.

The measurement: 12-hour survival over 2500 m terrain against sponge depth,
with everything else at production settings.

If sponge=0 survives 12/12, the recommendation is to turn it off by default
and accept the reflection, because the reflection does not end runs and the
sponge costs a factor of four in development.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from lid_test import build_on

print("2500 m terrain, clean, filtered, convection on, 12 h ceiling\n")
print(f"{'sponge':>7} {'survived':>9} {'max|u|':>8} {'refl max|du|':>13}")
for nsp in (0, 3, 5, 8):
    m = build_on(SigmaLevels(20), 2500.0, nsp)
    u0 = m.u.copy()
    dt = m.max_dt()
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    ok = np.isfinite(m.u).all()
    print(f"{nsp:7d} {done:6d}/12 "
          f"{(np.abs(m.u).max() if ok else float('nan')):8.1f} "
          f"{(np.abs(m.u-u0).max() if ok else float('nan')):13.2f}",
          flush=True)

"""
Is the tall-terrain failure the boundary, or the timestep?

    2500 m terrain, 12 h    rigid  radiative
    no sponge                9/12       3/12

The boundary flux is applied EXPLICITLY to the thinnest layer in the column
(sigma stretch 1.4 puts the smallest dsigma at the top), and its own CFL limit
only enters `max_dt` once a flux has been computed -- which is never on the
first call, because `_top_flux` is None until the first tendency evaluation.
So the first steps run at a timestep chosen without knowing about the
boundary at all.

That is a specific, checkable claim: halve dt and the run should go further.
If it does not, the boundary is genuinely destabilising over terrain and the
sponge stays.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from lid_test import build_on

print("2500 m terrain, no sponge, radiative lid, 12 h\n")
print(f"{'dt factor':>10} {'dt (s)':>8} {'recheck':>8} {'survived':>9} "
      f"{'max|u|':>8}")
for factor, recheck in ((1.0, 50), (0.5, 50), (0.25, 10)):
    m = build_on(SigmaLevels(20), 2500.0, 0)
    m.radiative_top = True
    dt = m.max_dt() * factor
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt, recheck_steps=recheck)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    ok = np.isfinite(m.u).all()
    print(f"{factor:10.2f} {dt:8.2f} {recheck:8d} {done:6d}/12 "
          f"{(np.abs(m.u).max() if ok else float('nan')):8.1f}", flush=True)

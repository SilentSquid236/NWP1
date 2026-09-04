"""
Is the 4000 m failure a timestep problem or a physics problem?

max|sigma_dot| grows an order of magnitude before the run dies (3.9e-05 ->
2.9e-04), which is what a vertical-CFL violation would look like. The adaptive
timestep should be catching that, but 'should' is not a measurement.

Halving and quartering dt is decisive: if the run extends roughly in
proportion, it is CFL; if it dies at the same forecast hour, the timestep is
innocent and the cause is physical.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from lid_test import build_on

print("4000 m terrain, clean, filtered, sponge=8, convection on")
print(f"{'dt factor':>10} {'dt (s)':>8} {'recheck':>8} {'survived':>9} "
      f"{'max|u|':>8}")
for factor, recheck in ((1.0, 50), (0.5, 50), (1.0, 5)):
    m = build_on(SigmaLevels(20), 4000.0, 8)
    dt = m.max_dt() * factor
    done = 0
    for _ in range(12):
        m.run(3600, dt=dt, recheck_steps=recheck)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    fin = np.isfinite(m.u).all()
    print(f"{factor:10.2f} {dt:8.2f} {recheck:8d} {done:6d}/12 "
          f"{(np.abs(m.u).max() if fin else float('nan')):8.1f}", flush=True)

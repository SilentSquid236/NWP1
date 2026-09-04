"""Terrain rows re-measured with dry convective adjustment."""
import numpy as np
np.seterr(all="ignore")
from sigma import SigmaLevels
from lid_test import build_on

print("clean, filtered, mixing+drag on, sponge=8, 12 h ceiling")
print(f"{'terrain':>8} {'convection':>11} {'survived':>9} {'max|u|':>8}")
for hgt in (1000.0, 2500.0, 4000.0):
    for conv in (False, True):
        m = build_on(SigmaLevels(20), hgt, 8)
        m.convection = conv
        dt = m.max_dt(); done = 0
        for _ in range(12):
            m.run(3600, dt=dt)
            if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
                break
            done += 1
        fin = np.isfinite(m.u).all()
        print(f"{hgt:8.0f} {str(conv):>11} {done:6d}/12 "
              f"{(np.abs(m.u).max() if fin else float('nan')):8.1f}", flush=True)

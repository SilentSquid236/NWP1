"""
Where does the terrain case fail, now that the initialization artifacts are
gone? Baseline first (does filtering alone fix it?), then locate the mode.
"""
import numpy as np
np.seterr(all="ignore")
from initialization import filter_initial_state
from subgrid import balance_initial_state
import probe_failure as P


def prep(hgt, noise=0.0, filt=True, **kw):
    m = P.build(hgt, noise, dT=1.5, clip=None, **kw)
    if filt:
        m.u, m.v, m.theta = filter_initial_state(m.u, m.v, m.theta, m.grid)
        m.u, m.v, _ = balance_initial_state(m.u, m.v, m.grid, verbose=False)
    return m


def survive(m, hours=12):
    dt = m.max_dt(); done = 0
    for _ in range(hours):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
            break
        done += 1
    fin = np.isfinite(m.u).all()
    return done, (np.abs(m.u).max() if fin else float("nan"))


if __name__ == "__main__":
    import sys
    if sys.argv[1] == "baseline":
        print("clean states, mixing+drag on, with and without the filter")
        for hgt in (1000.0, 2500.0):
            for filt in (False, True):
                d, u = survive(prep(hgt, filt=filt))
                print("  terrain %5.0f m  filter=%-5s -> %2d/12 h, max|u| %.1f"
                      % (hgt, filt, d, u), flush=True)
    elif sys.argv[1] == "locate":
        P.locate(prep(2500.0), hours=6.0, snaps=6)

"""
P-40: does a higher eddy-diffusivity ceiling cost anything where it matters?

2500 m is the agreed target terrain (P-01) and already runs 12/12 at the
K_MAX = 100 default. The 4000 m ladder says 150 buys two hours there. A
default is only worth changing if it does not regress the case the project is
actually for.

L11 is the reason this exists: the sponge default was cut from 5 to 3 once on
a plausible reading, and a regression suite -- not self-review -- is what
caught that it took the decisive noisy case from 12/12 to 11/12.

Run:  python kmax_production.py
"""

import numpy as np
np.seterr(all="ignore")

import turbulence
from sigma import SigmaLevels
from lid_test import build_on

BAND = slice(8, 14)


def main():
    print("2500 m terrain, sponge 5, clean, filtered, 12 h ceiling")
    print(f"{'K_MAX':>7} {'survived':>9} {'max|u|':>8} {'jet':>7} "
          f"{'minRi':>9} {'maxK':>8} {'N2_mid':>11}")

    for k in (100.0, 150.0, 300.0):
        m = build_on(SigmaLevels(20), 2500.0, 5, k_max=k)
        dt = m.max_dt()
        done = 0
        for _ in range(12):
            m.run(3600, dt=dt)
            if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150:
                break
            done += 1

        if np.isfinite(m.u).all():
            Ri, N2, _, _ = turbulence.richardson(m.u, m.v, m.theta, m.pi, m.lev)
            print(f"{k:7.0f} {done:6d}/12 {np.abs(m.u).max():8.1f} "
                  f"{np.abs(m.u[BAND]).max():7.1f} "
                  f"{np.nanmin(Ri[BAND]):9.3f} {m._K_last.max():8.1f} "
                  f"{np.nanmean(N2[BAND]):11.3e}", flush=True)
        else:
            print(f"{k:7.0f} {done:6d}/12   (non-finite)", flush=True)


if __name__ == "__main__":
    main()

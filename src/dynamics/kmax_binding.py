"""
P-40: does the eddy-diffusivity ceiling buy forecast hours, and HOW?

WHAT IS ALREADY KNOWN

The first K_MAX ladder was a broken experiment (P-51): it assigned
`turbulence.K_MAX` at runtime, which nothing reads, so three settings returned
byte-identical survival and were recorded as an elimination. Re-run through
the constructor, 4000 m terrain with an 8-level sponge gives

    K_MAX  100 -> 6/12      300 -> 8/12      1000 -> 8/12

Two forecast hours, and a saturation somewhere between 300 and 1000.

WHY A SURVIVAL COUNT IS NOT ENOUGH TO ACCEPT THIS

L2. A scheme can buy stability by dissipating the thing that kills the run, or
by flattening the flow until nothing is left to break. Both look identical on
a survival count, and the sponge already failed exactly this way (P-16, P-49).
So the discriminators are measured alongside survival, and one of them is
written to be able to fail.

PREDICTIONS, WRITTEN BEFORE THE RUNS (2026-09-11)

If the extra hours come from DISSIPATING the breaking wave:

  P1  survival rises between 100 and 300
  P2  max|u| does NOT fall as the ceiling rises.  <-- the one that can only
      fail. A higher ceiling that flattens the jet is suppression, and the
      two extra hours would be worthless.
  P3  the overturning fraction at a given hour falls as the ceiling rises:
      more mixing is available exactly where Ri <= 0
  P4  mid-level stratification away from the mountain is unchanged. A scheme
      mixing out the whole column would show N^2 falling domain-wide.

If instead the hours come from SUPPRESSION, P2 fails: max|u| falls, the jet
weakens, and the run survives because it has been damped into something
that cannot break.

On the 300 -> 1000 saturation, the proposed mechanism is:

  P5  the REALIZED diffusivity never reaches 1000, so the ceiling stops
      binding somewhere above 300 and raising it further changes nothing.
      Falsifiable: if realized K does reach ~1000 and survival is still 8/12,
      this explanation is wrong and the saturation is something else.

A first probe at hour 1 with k_max = 1000 gave a realized max K of 51.3 --
the ceiling is not binding at all early on, which is consistent with P5 and
also says the ladder only starts to separate once the wave breaks.

Run:  python kmax_binding.py            # the full ladder
      python kmax_binding.py 100 300    # selected rungs
"""

import sys
import time

import numpy as np
np.seterr(all="ignore")

import turbulence
from sigma import SigmaLevels
from lid_test import build_on

HOURS = 12
TERRAIN = 4000.0
SPONGE = 8
LADDER = (100.0, 150.0, 200.0, 250.0, 300.0, 1000.0)

# Levels 0-13 are lid to mid-troposphere on this 20-level grid; the sponge
# occupies 0-7, so Ri there is not a statement about the atmosphere.
RI_BAND = slice(8, 14)


def column_stats(m):
    """Everything measured at the end of a forecast hour."""
    Ri, N2, S2, dz = turbulence.richardson(m.u, m.v, m.theta, m.pi, m.lev)
    K = m._K_last
    conv = m._conv_info or {}

    finite = np.isfinite(m.u).all()
    return {
        "max|u|": float(np.abs(m.u).max()) if finite else float("nan"),
        "mean|u|": float(np.abs(m.u).mean()) if finite else float("nan"),
        # The jet is the free-troposphere flow away from the sponge. If the
        # ceiling is buying survival by flattening it, this is where it shows.
        "jet": float(np.abs(m.u[RI_BAND]).max()) if finite else float("nan"),
        "minRi": float(np.nanmin(Ri[RI_BAND])) if finite else float("nan"),
        # Stratification well away from the breaking layer. P4 watches this.
        "N2_mid": float(np.nanmean(N2[RI_BAND])) if finite else float("nan"),
        "maxK": float(np.max(K)) if K is not None else float("nan"),
        "meanK": float(np.mean(K)) if K is not None else float("nan"),
        "overturn%": 100.0 * float(conv.get("unstable_before", np.nan)),
    }


def at_ceiling(K, k_max):
    """Fraction of interfaces the ceiling is actually clipping."""
    if K is None:
        return float("nan")
    return 100.0 * float(np.mean(K >= 0.999 * k_max))


def run_one(k_max, hours=HOURS, verbose=True):
    m = build_on(SigmaLevels(20), TERRAIN, SPONGE, k_max=k_max)
    dt = m.max_dt()
    trace, survived = [], 0
    t0 = time.time()

    for h in range(1, hours + 1):
        m.run(3600, dt=dt)
        dead = (not np.isfinite(m.u).all()) or np.abs(m.u).max() > 150
        s = column_stats(m)
        s["hour"] = h
        s["clip%"] = at_ceiling(m._K_last, k_max)
        trace.append(s)
        if verbose:
            print(f"    h{h:02d}  max|u| {s['max|u|']:6.1f}  jet {s['jet']:6.1f}"
                  f"  minRi {s['minRi']:8.3f}  maxK {s['maxK']:7.1f}"
                  f"  clip {s['clip%']:5.2f}%  overturn {s['overturn%']:5.3f}%"
                  f"  N2mid {s['N2_mid']:.3e}", flush=True)
        if dead:
            break
        survived = h

    return {
        "k_max": k_max,
        "survived": survived,
        "trace": trace,
        "wall_s": time.time() - t0,
    }


def main():
    ladder = ([float(a) for a in sys.argv[1:]] if len(sys.argv) > 1
              else list(LADDER))

    print(f"P-40: {TERRAIN:.0f} m terrain, sponge {SPONGE}, clean, filtered, "
          f"{HOURS} h ceiling")
    print("predictions are in this file's docstring, written before the runs\n")

    results = []
    for k_max in ladder:
        print(f"  K_MAX = {k_max:.0f}")
        r = run_one(k_max)
        results.append(r)
        print(f"    -> {r['survived']}/{HOURS} h  ({r['wall_s']:.0f} s wall)\n",
              flush=True)

    print("\nSUMMARY  (last surviving hour of each run)")
    print(f"{'K_MAX':>7} {'survived':>9} {'max|u|':>8} {'jet':>7} {'minRi':>9} "
          f"{'maxK':>8} {'clip%':>7} {'overturn%':>10} {'N2_mid':>11}")
    for r in results:
        last = r["trace"][r["survived"] - 1] if r["survived"] else r["trace"][-1]
        print(f"{r['k_max']:7.0f} {r['survived']:6d}/{HOURS} "
              f"{last['max|u|']:8.1f} {last['jet']:7.1f} {last['minRi']:9.3f} "
              f"{last['maxK']:8.1f} {last['clip%']:7.2f} "
              f"{last['overturn%']:10.3f} {last['N2_mid']:11.3e}")

    # The comparison that decides dissipation vs suppression is between runs
    # at the SAME hour, not at each run's own last hour -- a run that lives
    # longer is compared at a later, harder time otherwise.
    common = min(r["survived"] for r in results) if results else 0
    if common:
        print(f"\nAT A COMMON HOUR (h{common}), which is what P2 and P4 need")
        print(f"{'K_MAX':>7} {'max|u|':>8} {'jet':>7} {'minRi':>9} "
              f"{'maxK':>8} {'clip%':>7} {'overturn%':>10} {'N2_mid':>11}")
        for r in results:
            s = r["trace"][common - 1]
            print(f"{r['k_max']:7.0f} {s['max|u|']:8.1f} {s['jet']:7.1f} "
                  f"{s['minRi']:9.3f} {s['maxK']:8.1f} {s['clip%']:7.2f} "
                  f"{s['overturn%']:10.3f} {s['N2_mid']:11.3e}")

    return results


if __name__ == "__main__":
    main()

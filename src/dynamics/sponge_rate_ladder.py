"""
The sponge RATE was never measured. It should have been.

`sponge_rate` defaults to 1/900 s -- a fifteen-minute relaxation time, applied
to the top five of twenty levels, relaxing toward a state frozen at hour zero.
Nothing justifies that number; it was chosen, not derived.

What the depth experiment just measured:

    sponge levels |   0  |  5   |  8
    eddy energy   | 3.18 | 0.85 | 0.96   (x per day, >1 = weather develops)
    max|v|        | 12.7 |  2.1 |  1.6

The DEFAULT configuration turns 3.18x/day growth into 0.85x/day DECAY. A
frozen reference damps whatever deviates from the initial state, and a
developing baroclinic wave is exactly that. The sponge is not just absorbing
gravity waves near the lid; it is absorbing the weather.

This ladder measures both sides of the trade at once -- development, which the
sponge must not suppress, and mountain-wave reflection, which is why it exists
-- across rates from fifteen minutes to six hours.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
from sponge_depth_weather import growth
from sponge_divergent import terrain_case
from primitive_sigma import PrimitiveSigma
import lid_test


def terrain_with_rate(nsp, rate, hours=6):
    m = lid_test.build_on(SigmaLevels(20), 2500.0, nsp)
    prof = np.zeros_like(m._sponge)
    for k in range(nsp):
        frac = (nsp - k) / nsp
        prof[k, 0, 0] = rate * 0.5 * (1 - np.cos(np.pi * frac))
    m._sponge = prof
    u0 = m.u.copy()
    dt = m.max_dt()
    done = 0
    for _ in range(hours):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            break
        done += 1
    if not np.isfinite(m.u).all():
        return float("nan"), done
    return float(np.abs(m.u - u0).max()), done


def growth_with_rate(nsp, rate):
    import sponge_depth_weather as S
    orig = PrimitiveSigma.__init__

    def patched(self, *a, **kw):
        kw.setdefault("sponge_rate", rate)
        orig(self, *a, **kw)

    PrimitiveSigma.__init__ = patched
    try:
        return S.growth(nsp)
    finally:
        PrimitiveSigma.__init__ = orig


if __name__ == "__main__":
    print("sponge rate ladder: development vs reflection\n")
    print(f"{'levels':>7} {'timescale':>11} {'eddy x/day':>11} "
          f"{'refl max|du|':>13} {'survived':>9}")
    for nsp in (5, 8):
        for tau, label in ((900.0, "15 min"), (3600.0, "1 h"),
                           (6 * 3600.0, "6 h")):
            g, vmax, jet = growth_with_rate(nsp, 1.0 / tau)
            refl, done = terrain_with_rate(nsp, 1.0 / tau)
            print(f"{nsp:7d} {label:>11} {g:11.2f} {refl:13.2f} "
                  f"{done:6d}/6", flush=True)

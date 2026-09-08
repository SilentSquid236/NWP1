"""
Measure absorption where a reflected wave actually shows up.

The first assertion compared the disturbance AT THE LID and read the radiative
boundary as worse (32.2 m/s against 27.2). That is the wrong place to look: a
rigid lid holds the wave still and a radiative lid lets it move through, so
the top level is expected to be MORE active, not less. A reflected wave shows
up below, coming back down.

So: perturbation energy in the lower two thirds of the column, hour by hour,
rigid against radiative. If the boundary works, the two agree while the wave
is still on its way up and separate once it would have returned.
"""
import numpy as np
np.seterr(all="ignore")

from test_radiation import mountain_case


def run(label, sign, sponge, hours=10):
    m = mountain_case(sign=sign, sponge=sponge)
    u0 = m.u.copy()
    nz = m.lev.nz
    lower = slice(nz // 3, nz)
    dt = m.max_dt()
    print(f"\n{label}  (dt {dt:.1f} s)")
    print(f"{'hour':>5} {'E lower':>11} {'E top third':>12} "
          f"{'max|du|':>9} {'max|w~|':>10}")
    for hr in range(1, hours + 1):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            print(f"{hr:5d}   non-finite")
            return
        du = m.u - u0
        e_low = float((du[lower] ** 2 + m.v[lower] ** 2).sum())
        e_top = float((du[:nz // 3] ** 2 + m.v[:nz // 3] ** 2).sum())
        print(f"{hr:5d} {e_low:11.3e} {e_top:12.3e} "
              f"{np.abs(du).max():9.2f} "
              f"{np.abs(m.sigma_dot()).max():10.2e}", flush=True)


if __name__ == "__main__":
    run("rigid lid, no sponge", None, 0)
    run("radiative lid, no sponge", -1.0, 0)
    run("rigid lid + 5-level sponge", None, 5)

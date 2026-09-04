"""
Does the growth over tall terrain track the SPONGE BASE?

The 2500 m case grows at large scale (90dx, 45dx, 30dx -- not grid scale),
peaking at level k=5 of 20 counting from the top. The sponge covers k=0..4.
The growth peak sits exactly at the sponge's lower edge, which is what partial
reflection off an absorbing layer looks like.

DISCRIMINATOR: move the sponge base. If the peak follows it, the sponge is the
cause. If it stays at k=5 regardless, the sponge is innocent and the growth is
anchored to something physical (the mountain wave's own structure).

Also reports the minimum Richardson number aloft, which distinguishes wave
BREAKING (Ri < 0.25, a real physical process the model should represent) from
a numerical reflection.
"""
import numpy as np
np.seterr(all="ignore")

from turbulence import richardson
from initialization import filter_initial_state
from subgrid import balance_initial_state
import probe_failure as P


def prep(hgt=2500.0, sponge_levels=5, **kw):
    m = P.build(hgt, 0.0, dT=1.5, clip=None, sponge_levels=sponge_levels, **kw)
    m.u, m.v, m.theta = filter_initial_state(m.u, m.v, m.theta, m.grid)
    m.u, m.v, _ = balance_initial_state(m.u, m.v, m.grid, verbose=False)
    return m


def run(m, hours=6.0):
    u0 = m.u.copy()
    dt = m.max_dt()
    for _ in range(int(hours)):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            return None
    d = m.u - u0
    lev_amp = np.sqrt((d ** 2).mean(axis=(1, 2)))
    Ri, _N2, _S2, _dz = richardson(m.u, m.v, m.theta, m.pi, m.lev)
    return lev_amp, float(np.abs(d).max()), Ri


print("2500 m terrain, clean, filtered, mixing+drag on, 6 h")
print(f"{'sponge':>7} {'base k':>7} {'peak k':>7} {'max|du|':>9} "
      f"{'min Ri aloft':>13} {'Ri<0.25 pts':>12}")
for nsp in (0, 5, 8, 12):
    m = prep(sponge_levels=nsp)
    out = run(m)
    if out is None:
        print(f"{nsp:7d} {nsp:7d}   diverged before 6 h", flush=True)
        continue
    lev_amp, dmax, Ri = out
    aloft = Ri[:12]                      # above the boundary layer
    print(f"{nsp:7d} {nsp:7d} {int(np.argmax(lev_amp)):7d} {dmax:9.2f} "
          f"{np.nanmin(aloft):13.3f} {int((aloft < 0.25).sum()):12d}",
          flush=True)
    print("        level rms: "
          + " ".join(f"{a:5.2f}" for a in lev_amp), flush=True)

"""
P-50: WHERE does the radiative lid fail over terrain?

Eliminated already: the timestep (dt and dt/2 both 3/12) and radiating the
terrain's steady hydrostatic imprint (fixed by the running low-pass; survival
unchanged).

A specific suspect remains, and it is written down in another module's
docstring. `remove_divergence_spectral` says: "The FFT assumes periodicity,
which this domain does not have, so expect some error in the outermost cells.
The relaxation zone overwrites those anyway." The radiation condition uses the
same FFT on a `replicate` domain -- and its output is applied straight to the
prognostic surface pressure, where nothing overwrites the edges.

So the question is spatial, and the probe records WHERE the disturbance is:
in the outermost cells, or over the mountain.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels, hydrostatic_geopotential
from radiation import radiative_top_flux
from lid_test import build_on

EDGE = 5


def where(a):
    """Split a 2D field into edge frame and interior."""
    edge = np.concatenate([a[:EDGE].ravel(), a[-EDGE:].ravel(),
                           a[:, :EDGE].ravel(), a[:, -EDGE:].ravel()])
    inner = a[EDGE:-EDGE, EDGE:-EDGE]
    return float(np.abs(edge).max()), float(np.abs(inner).max())


m = build_on(SigmaLevels(20), 2500.0, 0)
m.radiative_top = True
u0 = m.u.copy()
dt = m.max_dt()
ny, nx = m.pi.shape
print(f"2500 m terrain, no sponge, radiative lid, dt={dt:.2f}s\n")
print(f"{'hour':>5} {'|F| edge':>10} {'|F| inner':>10} {'ratio':>7} "
      f"{'dpi edge':>10} {'dpi inner':>10} {'max|u|':>8} {'min p_s':>9}")
pi0 = m.pi.copy()
for hr in range(1, 13):
    m.run(3600, dt=dt)
    if not np.isfinite(m.u).all():
        print(f"{hr:5d}   non-finite")
        break
    phi = hydrostatic_geopotential(m.theta, m.pi, m.lev, phi_surface=m.phi_s)
    F = radiative_top_flux(phi, m.theta, m.pi, m.lev, m.grid, sign=-1.0,
                           reference=m._phi_top_ref)
    fe, fi = where(F)
    de, di = where(m.pi - pi0)
    print(f"{hr:5d} {fe:10.3e} {fi:10.3e} {fe/max(fi,1e-30):7.2f} "
          f"{de:10.1f} {di:10.1f} {np.abs(m.u).max():8.2f} "
          f"{np.nanmin(m.surface_pressure):9.0f}", flush=True)

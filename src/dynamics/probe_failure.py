"""
Probe the noisy-state failure instead of patching it.

Records, every step: min/max surface pressure, max|u|, max|sigma_dot|, and the
per-level maximum of every momentum tendency term separately, so the term that
runs away is identified rather than guessed.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import (SigmaLevels, RD, G0, P0, KAPPA, continuity,
                   hydrostatic_geopotential, pressure_gradient_force,
                   vertical_advection)
from primitive_sigma import PrimitiveSigma
from subgrid import balance_initial_state


def build(hgt=0.0, noise=1.2, seed=0, dT=6.0, clip=60.0, **kw):
    gr = CGrid(90, 88, 12e3, 12e3, f0=9.81e-5, beta=1.69e-11,
               edge_mode="replicate")
    lev = SigmaLevels(20)
    rng = np.random.default_rng(seed)
    h = hgt * np.exp(-(((gr.Xc - gr.Lx/2)/250e3)**2 + ((gr.Yc - gr.Ly/2)/250e3)**2))
    m = PrimitiveSigma(gr, lev, terrain=h, **kw)
    k_y = 2*np.pi/gr.Ly
    m.pi = 101325.0*np.exp(-G0*h/(RD*280.0)) - lev.p_top
    p = lev.pressure(m.pi)
    T = 288.0 - 55.0*(1 - p/p.max()) - dT*np.cos(k_y*gr.Yc)
    m.theta = T/(p/P0)**KAPPA
    # GEOSTROPHIC WIND FROM THE FULL SIGMA PGF, not from d(phi)/dy alone.
    #
    # On sigma surfaces the horizontal force has two terms and they largely
    # cancel over sloping ground. Using d(phi)/dy by itself keeps only the
    # first, so over a 2500 m mountain it implies an 845 m/s "balanced" wind:
    # an initialization artifact that has nothing to do with the model.
    phi = hydrostatic_geopotential(m.theta, m.pi, lev, phi_surface=m.phi_s)
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u[:] = -fy/gr.f0
    m.v[:] = fx/gr.f0
    m.u = gr.v_to_u(m.u) if m.u.shape != m.v.shape else m.u
    if clip is not None:
        m.u = np.clip(m.u, -clip, clip)
        m.v = np.clip(m.v, -clip, clip)
    if noise > 0:
        m.u += rng.normal(0, noise, m.u.shape)
        m.v += rng.normal(0, noise, m.v.shape)
        m.u, m.v, _ = balance_initial_state(m.u, m.v, gr, verbose=False)
    return m


def terms(m):
    """Momentum tendency, term by term, at the current state."""
    gr, lev = m.grid, m.lev
    u, v, theta, pi = m.u, m.v, m.theta, m.pi
    phi = hydrostatic_geopotential(theta, pi, lev, phi_surface=m.phi_s)
    dpi_dt, sd = continuity(u, v, pi, lev, gr)
    T_ref = (theta*(lev.pressure(pi)/P0)**KAPPA).mean(axis=(1, 2))
    fx, fy = pressure_gradient_force(phi, theta, pi, lev, gr,
                                     reference=T_ref if m.ref_pgf else None)
    v_at_u = gr.v_to_u(v)
    out = {
        "adv_h":   m._horiz_adv(u, u, v_at_u),
        "adv_v":   vertical_advection(u, sd, lev),
        "coriolis": gr.f_u*v_at_u,
        "pgf":     fx,
    }
    return out, dpi_dt, sd


def probe(m, max_steps=100000, label=""):
    dt = m.max_dt()
    print(f"\n{label}  dt={dt:.2f}s")
    print(f"{'step':>6} {'hour':>6} {'min p_s':>10} {'max|u|':>8} "
          f"{'max|sd|':>9} {'|adv_h|':>9} {'|adv_v|':>9} {'|cor|':>9} {'|pgf|':>9}")
    hist = []
    for n in range(max_steps):
        if n % 25 == 0 or n < 5:
            t, dpi, sd = terms(m)
            ps = m.surface_pressure
            row = (n, m.time/3600, np.nanmin(ps), np.nanmax(np.abs(m.u)),
                   np.nanmax(np.abs(sd)),
                   *[np.nanmax(np.abs(t[k])) for k in
                     ("adv_h", "adv_v", "coriolis", "pgf")])
            hist.append(row)
            print(f"{row[0]:6d} {row[1]:6.2f} {row[2]:10.1f} {row[3]:8.2f} "
                  f"{row[4]:9.2e} {row[5]:9.2e} {row[6]:9.2e} {row[7]:9.2e} "
                  f"{row[8]:9.2e}", flush=True)
        if not (np.isfinite(m.u).all() and np.isfinite(m.pi).all()):
            print(f"  ** non-finite at step {n}, t={m.time/3600:.3f} h")
            break
        if np.nanmin(m.surface_pressure) <= 0:
            j, i = np.unravel_index(np.nanargmin(m.surface_pressure),
                                    m.pi.shape)
            print(f"  ** surface pressure <= 0 at step {n}, "
                  f"t={m.time/3600:.3f} h, (j,i)=({j},{i}), "
                  f"p_s={np.nanmin(m.surface_pressure):.1f} Pa")
            break
        m.step(dt)
        if n % 50 == 0:
            dt = min(dt, m.max_dt())
    return hist


if __name__ == "__main__":
    import sys
    which = sys.argv[1] if len(sys.argv) > 1 else "noisy"
    if which == "noisy":
        probe(build(0.0, 1.2), 700, "flat + noise 1.2")
    elif which == "clean":
        probe(build(0.0, 0.0), 700, "flat, no noise")
    elif which == "terrain":
        probe(build(2500.0, 0.0), 700, "2500 m terrain, no noise")


# ---------------------------------------------------------------------------
# WHERE does the growth live? Level, latitude band, and horizontal scale.
# ---------------------------------------------------------------------------
def locate(m, hours=1.0, snaps=5):
    """Track the growing part of the wind in level / y / wavenumber."""
    dt = m.max_dt()
    u0 = m.u.copy()
    t_end = hours*3600
    marks = np.linspace(t_end/snaps, t_end, snaps)
    print(f"\nLOCATING GROWTH   dt={dt:.2f}s, {hours} h")
    for t in marks:
        while m.time < t - 1e-9:
            m.step(min(dt, t - m.time))
            if m.step_count % 50 == 0:
                dt = min(dt, m.max_dt())
            if not np.isfinite(m.u).all():
                break
        if not np.isfinite(m.u).all():
            print(f"  t={m.time/3600:.2f} h  non-finite")
            break
        d = m.u - u0
        lev_amp = np.sqrt((d**2).mean(axis=(1, 2)))
        y_amp = np.sqrt((d**2).mean(axis=(0, 2)))
        k = np.abs(np.fft.rfft(d - d.mean(axis=2, keepdims=True), axis=2))
        k_amp = np.sqrt((k**2).mean(axis=(0, 1)))
        nx = d.shape[2]
        # wavelength in grid cells for the three strongest nonzero wavenumbers
        top = np.argsort(k_amp[1:])[::-1][:3] + 1
        print(f"  t={m.time/3600:4.2f} h  max|du|={np.abs(d).max():7.2f}  "
              f"level rms (top->bottom): "
              f"{' '.join(f'{a:5.2f}' for a in lev_amp[::4])}")
        print(f"           peak level k={int(np.argmax(lev_amp))}/{len(lev_amp)-1}"
              f" (0=top), peak y-row={int(np.argmax(y_amp))}/{len(y_amp)-1}, "
              f"dominant wavelengths = "
              f"{', '.join(f'{nx/kk:.1f}dx' for kk in top)}")

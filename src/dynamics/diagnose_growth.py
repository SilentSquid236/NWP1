"""
Instability diagnosis: measure where the energy enters.

Rather than adding a sink and seeing whether the model survives longer, this
decomposes the kinetic-energy tendency term by term, level by level, and
wavenumber by wavenumber. The output identifies which term, which levels, and
which SCALE is growing, which narrows the cause to something specific instead
of leaving a list of suspects.

Four questions, each answered by a separate measurement:

  1. WHICH TERM?     dKE/dt = integral of u . (du/dt)_term, evaluated per term.
  2. WHICH LEVELS?   the same, resolved vertically.
  3. WHICH SCALE?    amplitude spectrum of u over time -- grid-scale growth
                     means a numerical mode; large-scale means a physical
                     instability or a real imbalance.
  4. ROTATIONAL OR DIVERGENT?  Helmholtz split of the growing part. Gravity
                     waves are divergent; balanced flow is rotational.

Run:  python diagnose_growth.py
"""

import numpy as np

from grid import CGrid
from sigma import (SigmaLevels, hydrostatic_geopotential, continuity,
                   vertical_advection, pressure_gradient_force,
                   RD, G0, P0, KAPPA)
from primitive_sigma import PrimitiveSigma
from subgrid import hyperdiffusion, remove_divergence_spectral


# ---------------------------------------------------------------------------
# 1 & 2. Energy budget by term and by level
# ---------------------------------------------------------------------------

def tendency_terms(m):
    """
    Evaluate each momentum tendency term separately, on the model's current
    state. Returns {name: (du, dv)}.
    """
    gr, lev = m.grid, m.lev
    u, v, theta, pi = m.u, m.v, m.theta, m.pi

    phi = hydrostatic_geopotential(theta, pi, lev, phi_surface=m.phi_s)
    _, sd = continuity(u, v, pi, lev, gr)
    fx, fy = pressure_gradient_force(phi, theta, pi, lev, gr)

    v_at_u = gr.v_to_u(v)
    u_at_v = gr.u_to_v(u)

    terms = {
        "horiz advection": (-m._horiz_adv(u, u, v_at_u),
                            -m._horiz_adv(v, u_at_v, v)),
        "vert advection": (-vertical_advection(u, sd, lev),
                           -vertical_advection(v, sd, lev)),
        "coriolis": (gr.f_u * v_at_u, -gr.f_v * u_at_v),
        "pressure gradient": (fx, fy),
    }
    if m.hyper > 0:
        terms["hyperdiffusion"] = (hyperdiffusion(u, gr, m.hyper),
                                   hyperdiffusion(v, gr, m.hyper))
    return terms


def energy_budget(m):
    """
    dKE/dt contribution from each term, total and per level.

    Coriolis should contribute ~zero: it rotates the wind without doing work.
    A non-negligible Coriolis contribution means the staggered interpolation
    is not energy-neutral, which is a specific, findable bug.
    """
    ds = m.lev.dsigma.reshape(-1, 1, 1)
    weight = m.pi * ds                       # mass per unit area in each layer

    out = {}
    for name, (du, dv) in tendency_terms(m).items():
        contrib = (m.u * du + m.v * dv) * weight
        out[name] = {
            "total": float(contrib.sum()),
            "by_level": contrib.sum(axis=(1, 2)),
        }
    return out


# ---------------------------------------------------------------------------
# 3. Which scale is growing
# ---------------------------------------------------------------------------

def spectrum(field, grid):
    """
    Isotropic amplitude spectrum, binned by total wavenumber index.
    Index 1 is the domain scale; the last bin is the 2dx grid scale.
    """
    f = np.fft.fft2(field - field.mean(), axes=(-2, -1))
    power = np.abs(f) ** 2

    ny, nx = field.shape[-2:]
    kx = np.fft.fftfreq(nx) * nx
    ky = np.fft.fftfreq(ny) * ny
    KX, KY = np.meshgrid(kx, ky)
    K = np.sqrt(KX ** 2 + KY ** 2)

    nbins = min(nx, ny) // 2
    bins = np.zeros(nbins)
    for i in range(nbins):
        mask = (K >= i) & (K < i + 1)
        if mask.any():
            bins[i] = power[..., mask].sum()
    return bins


def scale_growth(m, hours=2, samples=8, level=None):
    """
    Track the amplitude spectrum of u over time and report growth per band.

    Grid-scale growth (largest wavenumbers) indicates a numerical mode.
    Large-scale growth indicates a physical instability or a genuine imbalance
    in the initial state.
    """
    k = m.lev.nz - 4 if level is None else level
    dt_total = hours * 3600 / samples

    spectra = [spectrum(m.u[k], m.grid)]
    for _ in range(samples):
        m.run(dt_total)
        if not np.isfinite(m.u).all():
            break
        spectra.append(spectrum(m.u[k], m.grid))

    return np.array(spectra)


# ---------------------------------------------------------------------------
# 4. Rotational vs divergent
# ---------------------------------------------------------------------------

def helmholtz_energy(u, v, grid):
    """Kinetic energy split into rotational and divergent parts, per level."""
    rot_u = np.empty_like(u)
    rot_v = np.empty_like(v)
    for k in range(u.shape[0]):
        rot_u[k], rot_v[k] = remove_divergence_spectral(u[k], v[k], grid)

    div_u, div_v = u - rot_u, v - rot_v
    return (float((rot_u ** 2 + rot_v ** 2).sum()),
            float((div_u ** 2 + div_v ** 2).sum()))


# ---------------------------------------------------------------------------
# Driver
# ---------------------------------------------------------------------------

def build_case(jet_dT=3.0, noise=1.2, terrain_m=400.0, nx=90, ny=88, dx=12e3,
               seed=0):
    """The failing configuration, reproduced exactly."""
    gr = CGrid(nx, ny, dx, dx, f0=9.81e-5, beta=1.69e-11, edge_mode="replicate")
    lev = SigmaLevels(20)
    rng = np.random.default_rng(seed)

    h = terrain_m * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                             ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h)

    k_y = 2 * np.pi / gr.Ly
    p_s = 101325.0 * np.exp(-G0 * h / (RD * 280.0))
    m.pi = p_s - lev.p_top
    p = lev.pressure(m.pi)
    T = 288.0 - 55.0 * (1 - p / p.max()) - jet_dT * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA

    phi = m.geopotential()
    for k in range(lev.nz):
        dphidy = 0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k]))
        m.u[k] = -dphidy / gr.f0                 # NO clipping: keep balance

    if noise:
        m.u += rng.normal(0, noise, m.u.shape)
        m.v += rng.normal(0, noise, m.v.shape)
        for k in range(lev.nz):
            m.u[k], m.v[k] = remove_divergence_spectral(m.u[k], m.v[k], gr)
    return m


def main():
    print("=" * 70)
    print("INSTABILITY DIAGNOSIS")
    print("=" * 70)

    m = build_case()
    print(f"\nCase: {m}")
    print(f"  jet {np.abs(m.u).max():.1f} m/s, dt {m.max_dt():.1f} s, "
          f"Lamb {m.lamb_wave_speed:.0f} m/s")

    # --- Q1: which term supplies the energy? ------------------------------
    print("\n" + "-" * 70)
    print("Q1. dKE/dt by term, at t=0")
    print("-" * 70)
    b = energy_budget(m)
    scale = max(abs(v["total"]) for v in b.values()) or 1.0
    for name, d in sorted(b.items(), key=lambda kv: -abs(kv[1]["total"])):
        bar = "#" * int(40 * abs(d["total"]) / scale)
        print(f"  {name:>20}  {d['total']:+.3e}  {bar}")
    print("\n  Coriolis should be ~0 (it does no work). A large value means")
    print("  the staggered averaging is not energy-neutral.")

    # --- Q2: which levels? -------------------------------------------------
    print("\n" + "-" * 70)
    print("Q2. Vertical distribution of the dominant term")
    print("-" * 70)
    dom = max(b.items(), key=lambda kv: abs(kv[1]["total"]))
    lv = dom[1]["by_level"]
    top5 = np.argsort(-np.abs(lv))[:5]
    print(f"  dominant term: {dom[0]}")
    for k in sorted(top5):
        print(f"    level {k:2d} (sigma={m.lev.sigma[k]:.3f}, "
              f"p~{m.lev.pressure(m.pi.mean())[k]/100:6.1f} hPa): {lv[k]:+.3e}")

    # --- Q4 before/after: rotational vs divergent -------------------------
    print("\n" + "-" * 70)
    print("Q4. Rotational vs divergent kinetic energy over time")
    print("-" * 70)
    rot0, div0 = helmholtz_energy(m.u, m.v, m.grid)
    print(f"  {'time':>6} {'rotational':>14} {'divergent':>14} {'div/rot':>10}")
    print(f"  {'t=0':>6} {rot0:>14.4e} {div0:>14.4e} {div0/rot0:>10.2e}")

    for hr in (0.25, 0.5, 1.0, 2.0):
        m.run(hr * 3600 if hr == 0.25 else 0.25 * 3600)
        if not np.isfinite(m.u).all():
            print(f"  {'':>6} DIVERGED before t={hr} h")
            break
        r, d = helmholtz_energy(m.u, m.v, m.grid)
        print(f"  {hr:>5.2f}h {r:>14.4e} {d:>14.4e} {d/r:>10.2e}")

    # --- Q3: which scale? --------------------------------------------------
    print("\n" + "-" * 70)
    print("Q3. Which spatial scale grows")
    print("-" * 70)
    m2 = build_case()
    sp = scale_growth(m2, hours=1.0, samples=6)
    if len(sp) > 1:
        first, last = sp[0], sp[-1]
        nb = len(first)
        print(f"  {'band':>22} {'growth factor':>15}")
        bands = [(1, 3, "domain scale"), (3, 8, "synoptic"),
                 (8, 20, "mesoscale"), (nb - 8, nb - 3, "near grid scale"),
                 (nb - 3, nb, "2dx grid scale")]
        for lo, hi, label in bands:
            lo, hi = max(1, lo), min(nb, hi)
            if hi <= lo:
                continue
            a, bb = first[lo:hi].sum(), last[lo:hi].sum()
            g = bb / a if a > 0 else np.inf
            flag = "  <-- fastest" if g == max(
                (last[max(1, l):min(nb, h)].sum() /
                 max(first[max(1, l):min(nb, h)].sum(), 1e-300))
                for l, h, _ in bands) else ""
            print(f"  {label:>22} {g:>15.3f}{flag}")
        print("\n  Grid-scale growth => numerical mode.")
        print("  Large-scale growth => physical instability or imbalance.")
    else:
        print("  diverged before a spectrum could be measured")

    print("\n" + "=" * 70)


if __name__ == "__main__":
    main()

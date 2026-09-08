"""
Validation for the radiative upper boundary.

The first test is the one that decides the sign, and it is a measurement
rather than an argument: the same mountain wave is integrated with sign = -1
and sign = +1, and whichever LOSES wave energy is the radiating one. Getting
this backwards turns the boundary into a source that pumps energy in at
exactly the rate it should let it out.

Run:  python test_radiation.py
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import (SigmaLevels, continuity, hydrostatic_geopotential,
                   RD, G0, P0, KAPPA)
from primitive_sigma import PrimitiveSigma
from radiation import (radiative_top_flux, brunt_vaisala,
                       horizontal_wavenumber)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def mountain_case(sign=None, sponge=0, hgt=1500.0, nx=64, ny=64):
    """Uniform flow over a ridge -- the cleanest mountain-wave generator."""
    gr = CGrid(nx, ny, 12e3, 12e3, f0=9.81e-5, beta=0.0, edge_mode="replicate")
    lev = SigmaLevels(20)
    h = hgt * np.exp(-(((gr.Xc - gr.Lx / 2) / 60e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h, sponge_levels=sponge,
                       radiative_top=sign is not None,
                       radiation_sign=(sign or -1.0))
    T0 = 260.0
    m.pi = 101325.0 * np.exp(-G0 * h / (RD * T0)) - lev.p_top
    p = lev.pressure(m.pi)
    T = T0 - 40.0 * (1 - p / p.max())
    m.theta = T / (p / P0) ** KAPPA
    m.u[:] = 15.0
    return m


def wave_energy(m, u0):
    """Perturbation kinetic energy: what the boundary should be removing."""
    du = m.u - u0
    return float((du ** 2 + m.v ** 2).sum())


# ---------------------------------------------------------------------------
def test_rigid_lid_is_the_default():
    """Nothing changes unless the boundary is asked for."""
    gr = CGrid(32, 32, 20e3, 20e3)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev)
    rng = np.random.default_rng(0)
    m.u[:] = rng.normal(0, 5, m.u.shape)
    m.v[:] = rng.normal(0, 5, m.v.shape)
    m.pi[:] = 90000.0
    _, sd = continuity(m.u, m.v, m.pi, lev, gr)
    report("with no radiative top, sigma_dot vanishes at both ends",
           np.abs(sd[0]).max() < 1e-15 and np.abs(sd[-1]).max() < 1e-15,
           f"max|sigma_dot| top {np.abs(sd[0]).max():.2e}, "
           f"ground {np.abs(sd[-1]).max():.2e}")


def test_flux_opens_the_top_and_not_the_ground():
    """Mass may leave through the lid; it may never leak through the surface."""
    gr = CGrid(32, 32, 20e3, 20e3)
    lev = SigmaLevels(20)
    rng = np.random.default_rng(1)
    u = rng.normal(0, 5, (20, 32, 32))
    v = rng.normal(0, 5, (20, 32, 32))
    pi = np.full((32, 32), 90000.0)
    F = np.full((32, 32), 3.0)
    _, sd = continuity(u, v, pi, lev, gr, top_flux=F)
    ok = (abs(np.abs(sd[0]).max() - 3.0 / 90000.0) < 1e-12
          and np.abs(sd[-1]).max() < 1e-15)
    report("a top flux opens the lid and leaves the ground closed", ok,
           f"sigma_dot top {sd[0].max():.3e} (expected "
           f"{3.0/90000.0:.3e}), ground {np.abs(sd[-1]).max():.2e}")


def test_flux_is_zero_for_a_uniform_column():
    """
    A horizontally uniform state has nothing to radiate: |k| = 0 for the only
    mode present. A boundary that leaks mass out of a state at rest would be a
    slow, invisible drift in surface pressure.
    """
    gr = CGrid(32, 32, 20e3, 20e3)
    lev = SigmaLevels(20)
    pi = np.full((32, 32), 90000.0)
    p = lev.pressure(pi)
    theta = (260.0 - 40.0 * (1 - p / p.max())) / (p / P0) ** KAPPA
    phi = hydrostatic_geopotential(theta, pi, lev)
    F = radiative_top_flux(phi, theta, pi, lev, gr)
    report("a horizontally uniform column radiates nothing",
           np.abs(F).max() < 1e-9,
           f"max|F| {np.abs(F).max():.2e} Pa/s")


def test_brunt_vaisala_is_physical():
    """N for a normal troposphere is about 0.01 1/s."""
    gr = CGrid(16, 16, 20e3, 20e3)
    lev = SigmaLevels(20)
    pi = np.full((16, 16), 90000.0)
    p = lev.pressure(pi)
    theta = (260.0 - 40.0 * (1 - p / p.max())) / (p / P0) ** KAPPA
    N = brunt_vaisala(theta, pi, lev)
    report("Brunt-Vaisala frequency is physical", 0.005 < N < 0.03,
           f"N = {N:.4f} 1/s (troposphere is ~0.01, stratosphere ~0.02)")


def test_sign_is_measured_not_assumed():
    """
    THE DECIDING TEST. Integrate the same mountain wave with both signs and
    keep whichever loses perturbation energy. One radiates; the other pumps.
    """
    out = {}
    for sign in (-1.0, +1.0):
        m = mountain_case(sign=sign)
        u0 = m.u.copy()
        e = []
        dt = m.max_dt()
        for _ in range(4):
            m.run(3600, dt=dt)
            if not np.isfinite(m.u).all():
                e.append(float("nan"))
                break
            e.append(wave_energy(m, u0))
        out[sign] = e

    rigid = mountain_case(sign=None)
    u0 = rigid.u.copy()
    er = []
    dt = rigid.max_dt()
    for _ in range(4):
        rigid.run(3600, dt=dt)
        if not np.isfinite(rigid.u).all():
            er.append(float("nan"))
            break
        er.append(wave_energy(rigid, u0))

    a, b = out[-1.0], out[+1.0]
    fin = [x for x in (a[-1], b[-1], er[-1]) if np.isfinite(x)]
    ok = len(fin) >= 2 and min(a[-1], b[-1]) < er[-1]
    report("one sign radiates and the other does not", bool(ok),
           f"perturbation energy at 4 h -- rigid lid {er[-1]:.3e}; "
           f"sign -1 {a[-1]:.3e}; sign +1 {b[-1]:.3e}")


def test_it_does_not_leak_a_resting_atmosphere():
    """
    A balanced atmosphere at rest over terrain must not lose mass through the
    lid. This is what the raw form got wrong: over a mountain the top-level
    geopotential perturbation is dominated by the terrain's steady hydrostatic
    imprint, and radiating a steady anomaly pumps mass out of those columns for
    the whole run.
    """
    gr = CGrid(48, 48, 12e3, 12e3, f0=9.81e-5, edge_mode="replicate")
    lev = SigmaLevels(20)
    h = 1500.0 * np.exp(-(((gr.Xc - gr.Lx / 2) / 100e3) ** 2 +
                          ((gr.Yc - gr.Ly / 2) / 100e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h, sponge_levels=0,
                       radiative_top=True)
    T0 = 260.0
    m.pi = 101325.0 * np.exp(-G0 * h / (RD * T0)) - lev.p_top
    m.theta = T0 / (lev.pressure(m.pi) / P0) ** KAPPA
    m0 = m.total_mass()
    m.run(3 * 3600, dt=m.max_dt())
    ok = np.isfinite(m.pi).all()
    drift = abs(m.total_mass() - m0) / m0 if ok else float("nan")
    umax = float(np.abs(m.u).max()) if ok else float("nan")
    report("a resting atmosphere over terrain does not leak through the lid",
           bool(ok and drift < 1e-3 and umax < 5.0),
           f"mass drift {drift:.2e} relative, spurious max|u| {umax:.3f} m/s "
           f"over {h.max():.0f} m terrain in 3 h")


def test_development_is_not_suppressed():
    """
    THE POINT OF THE MODULE. The sponge turns a growing baroclinic wave into a
    decaying one at every setting (P-49). A boundary that damps nothing must
    not.

    Measured over 48 h, eddy kinetic energy at 48 h over 6 h:

        sponge 5, rigid lid   0.34   decays
        no sponge, rigid lid  1.82
        no sponge, radiative  2.50

    This runs a short version of the same thing and requires the radiative
    configuration to still be growing where the sponge is not.
    """
    # 36 h over 18 h, not 24 over 12: both configurations are still shedding
    # the initial transient at 12 h, so an early window does not separate them
    # (0.61 against 0.73, which is no test at all). The curves diverge once
    # growth takes over.
    from radiation_vs_sponge import curve
    c_sponge = curve(5, False, hours=36, every=6)
    c_rad = curve(0, True, hours=36, every=6)
    r_sponge = c_sponge[-1] / c_sponge[2]
    r_rad = c_rad[-1] / c_rad[2]
    report("development survives the radiative lid where a sponge kills it",
           r_rad > 1.0 and r_rad > r_sponge * 1.5,
           f"eddy energy 36 h / 18 h -- sponge {r_sponge:.2f} (decaying), "
           f"radiative {r_rad:.2f}")


if __name__ == "__main__":
    print("\nRadiative upper boundary\n" + "=" * 66)
    for fn in (test_rigid_lid_is_the_default,
               test_flux_opens_the_top_and_not_the_ground,
               test_flux_is_zero_for_a_uniform_column,
               test_brunt_vaisala_is_physical,
               test_sign_is_measured_not_assumed,
               test_it_does_not_leak_a_resting_atmosphere,
               test_development_is_not_suppressed):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 66)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

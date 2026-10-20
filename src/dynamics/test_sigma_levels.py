"""
Tests for the sigma level placement (P-84, test L1) and for the explicit
vertical-mixing limit in PrimitiveSigma.max_dt.

Run:  python test_sigma_levels.py
"""

import numpy as np

from grid import CGrid
from sigma import SigmaLevels, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from turbulence import richardson, eddy_diffusivity

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def lowest_level_height(lev, T=285.0, ps=101325.0):
    p1 = lev.p_top + lev.sigma[-1] * (ps - lev.p_top)
    return RD * T / G0 * np.log(ps / p1)


def test_default_grid_unchanged():
    lev = SigmaLevels(20)
    ref = np.linspace(0.0, 1.0, 21) ** 1.4
    ok = np.array_equal(lev.sigma_half, ref) and lev.form == "power" and lev.stretch == 1.4
    report("default grid is bit-identical to (h/N)^1.4", ok, repr(lev))


def test_power_grid_is_thickest_at_the_ground():
    lev = SigmaLevels(20)
    ok = lev.dsigma[-1] > lev.dsigma[0] and abs(lev.dsigma[-1] - (1 - 0.95 ** 1.4)) < 1e-12
    report("power 1.4: thickest layer at the ground (the P-84 finding)", ok,
           f"bottom dsigma {lev.dsigma[-1]:.4f}, top {lev.dsigma[0]:.4f}")


def test_surface_form():
    lev = SigmaLevels(20, stretch=1.4, form="surface")
    mono = bool(np.all(np.diff(lev.sigma_half) > 0))
    ends = lev.sigma_half[0] == 0.0 and lev.sigma_half[-1] == 1.0
    bottom = abs(lev.dsigma[-1] - (1.0 / 20) ** 1.4) < 1e-12
    mirror = np.allclose(lev.dsigma[::-1], SigmaLevels(20).dsigma)
    ok = mono and ends and bottom and mirror
    report("surface 1.4: monotonic, 0..1, thinnest at the ground, mirror of power", ok,
           f"bottom dsigma {lev.dsigma[-1]:.4f}, top {lev.dsigma[0]:.4f}")


def test_bad_form_rejected():
    try:
        SigmaLevels(20, form="tanh")
        ok = False
    except ValueError:
        ok = True
    report("unknown form is rejected", ok, "ValueError raised" if ok else "accepted")


def test_lowest_level_heights():
    cases = {(20, 1.4, "power"): 235, (20, 1.0, "power"): 169, (20, 1.4, "surface"): 51,
             (30, 1.32, "surface"): 38, (40, 1.19, "surface"): 42}
    got = {k: lowest_level_height(SigmaLevels(k[0], stretch=k[1], form=k[2])) for k in cases}
    ok = all(abs(got[k] - v) < 1.5 for k, v in cases.items())
    report("lowest full level heights (285 K, 1013 hPa)", ok,
           ", ".join(f"{k[0]}/{k[2]}/{k[1]}: {got[k]:.0f} m" for k in cases))


def sheared_isothermal(model, T0=260.0, u_top=30.0):
    """Isothermal column, u rising from 0 to u_top over the lowest 5 % of sigma."""
    model.pi = np.full(model.pi.shape, 101325.0 - model.lev.p_top)
    p = model.lev.pressure(model.pi)
    model.theta = T0 / (p / P0) ** KAPPA
    prof = u_top * np.minimum(1.0, (1.0 - model.lev.sigma) * 20.0)
    model.u = np.broadcast_to(prof[:, None, None], model.u.shape).copy()
    model.v = np.zeros_like(model.v)


def test_mixing_limit_does_not_bind_on_default_grid():
    gr = CGrid(16, 16, 25e3, 25e3)
    lev = SigmaLevels(20)
    a = PrimitiveSigma(gr, lev)
    b = PrimitiveSigma(gr, lev, mixing=False)
    for m in (a, b):
        sheared_isothermal(m)
    ok = a.max_dt() == b.max_dt()
    report("default grid: the mixing limit leaves dt exactly unchanged", ok,
           f"dt {a.max_dt():.3f} s with mixing, {b.max_dt():.3f} s without")


def test_mixing_limit_binds_on_thin_bottom_layers():
    gr = CGrid(16, 16, 25e3, 25e3)
    lev = SigmaLevels(20, stretch=1.4, form="surface")
    m = PrimitiveSigma(gr, lev)
    n = PrimitiveSigma(gr, lev, mixing=False)
    for x in (m, n):
        sheared_isothermal(x)
    Ri, _, S2, dz = richardson(m.u, m.v, m.theta, m.pi, lev)
    K = eddy_diffusivity(Ri, S2, m.ri_crit, m.k_max, m.mixing_length)
    expect = 0.4 * float(dz.min()) ** 2 / float(K.max())
    ok = m.max_dt() < n.max_dt() and abs(m.max_dt() - expect) < 1e-9
    report("surface grid with strong low-level shear: the mixing limit binds", ok,
           f"dt {m.max_dt():.2f} s (expected {expect:.2f}), {n.max_dt():.2f} s without mixing; "
           f"K max {float(K.max()):.0f} m^2/s, dz min {float(dz.min()):.0f} m")


if __name__ == "__main__":
    print("\nsigma levels and the vertical-mixing limit\n")
    for t in (test_default_grid_unchanged, test_power_grid_is_thickest_at_the_ground,
              test_surface_form, test_bad_form_rejected, test_lowest_level_heights,
              test_mixing_limit_does_not_bind_on_default_grid,
              test_mixing_limit_binds_on_thin_bottom_layers):
        t()
    print(f"\n{sum(results)}/{len(results)} passed")
    raise SystemExit(0 if all(results) else 1)

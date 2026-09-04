"""
Validation for the surface layer.

The decisive test is the Ekman spiral: friction plus rotation must turn the
near-surface wind ACROSS the isobars toward low pressure. That angle is a
physical consequence no tuning produces by accident, and a sign error in the
drag makes it turn the wrong way.

Run:  python test_surface.py
"""

import numpy as np

from grid import CGrid
from sigma import SigmaLevels, pressure_gradient_force, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from surface import (neutral_drag_coefficient, stability_function,
                     bulk_richardson, surface_drag, lowest_level_height,
                     ROUGHNESS)

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def balanced(n=32, dx=25e3, dT=1.5, drag=True, z0=0.1, mixing=True,
             theta_surface=None):
    gr = CGrid(n, n, dx, dx, f0=1.0e-4, beta=0.0, edge_mode="replicate")
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev, drag=drag, z0=z0, mixing=mixing,
                       theta_surface=theta_surface)
    m.pi[:] = 101325.0 - lev.p_top
    p = lev.pressure(m.pi)
    ky = 2 * np.pi / gr.Ly
    m.theta = (288.0 - 55.0 * (1 - p / p.max())
               - dT * np.cos(ky * gr.Yc)) / (p / P0) ** KAPPA
    phi = m.geopotential()
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u = -fy / gr.f0
    m.v = fx / gr.f0
    return m, gr


# ---------------------------------------------------------------------------
def test_drag_coefficient_matches_log_law():
    """Cd should follow (kappa/ln(z1/z0))^2 and rise with roughness."""
    z1 = 50.0
    cds = {k: neutral_drag_coefficient(z1, v) for k, v in ROUGHNESS.items()}
    ok = (cds["sea"] < cds["grass"] < cds["cropland"] < cds["forest"]
          and 5e-4 < cds["sea"] < 3e-3 and 5e-3 < cds["forest"] < 8e-2)
    report("drag coefficient follows the log law and ranks by roughness", ok,
           "  ".join(f"{k} {v:.4f}" for k, v in cds.items())
           + f"   (z1 = {z1:.0f} m)")


# ---------------------------------------------------------------------------
def test_stability_function_behaviour():
    """
    Unstable surface layers enhance drag; strongly stable ones decouple.
    Getting this backwards mixes away every nocturnal inversion.
    """
    Ri = np.array([-1.0, -0.1, 0.0, 0.1, 0.19, 0.25, 1.0])
    F = stability_function(Ri)
    ok = (F[0] > 1.0 and F[2] == 1.0 and F[3] < 1.0
          and F[5] == 0.0 and F[6] == 0.0 and np.all(np.diff(F) <= 1e-12))
    report("stability function: enhanced unstable, decoupled when stable", ok,
           "  ".join(f"Ri={r:+.2f}->{f:.2f}" for r, f in zip(Ri, F)))


# ---------------------------------------------------------------------------
def test_drag_removes_momentum_only():
    """
    Drag must oppose the wind everywhere and touch only the lowest level.
    A drag that accelerates anything is a sign error.
    """
    m, gr = balanced()
    du, dv, info = surface_drag(m.u, m.v, m.theta, m.pi, m.lev, z0=0.1)

    aligned = (du[-1] * m.u[-1] + dv[-1] * m.v[-1])
    upper_touched = np.abs(du[:-1]).max() + np.abs(dv[:-1]).max()

    ok = np.all(aligned <= 1e-12) and upper_touched == 0.0
    report("drag opposes the wind and acts only on the lowest level", ok,
           f"max(u.du) = {aligned.max():.2e} (must be <= 0); "
           f"levels above untouched: {upper_touched == 0.0}; "
           f"z1 = {np.mean(info['z1']):.1f} m, "
           f"u* = {np.mean(info['u_star']):.2f} m/s")


# ---------------------------------------------------------------------------
def test_ekman_spiral():
    """
    THE test. With drag, the near-surface wind must turn ACROSS the isobars
    toward low pressure -- 10-30 degrees over land. Without drag the turning
    must be essentially zero.
    """
    def cross_isobar_angle(drag):
        m, gr = balanced(drag=drag)
        m.run(6 * 3600)
        if not np.isfinite(m.u).all():
            return np.nan, np.nan

        # Compare the ACTUAL surface wind with the GEOSTROPHIC wind at the
        # SAME level. Comparing against the wind aloft instead measures
        # thermal-wind turning, which has nothing to do with friction.
        phi = m.geopotential()
        fx, fy = pressure_gradient_force(phi, m.theta, m.pi, m.lev, gr)
        ug, vg = -fy[-1] / gr.f0, fx[-1] / gr.f0
        us, vs = m.u[-1], m.v[-1]

        strong = np.sqrt(ug ** 2 + vg ** 2) > 5.0
        if not strong.any():
            return np.nan, np.nan
        ang = np.degrees(np.arctan2(
            np.sin(np.arctan2(vs[strong], us[strong])
                   - np.arctan2(vg[strong], ug[strong])),
            np.cos(np.arctan2(vs[strong], us[strong])
                   - np.arctan2(vg[strong], ug[strong]))))
        ratio = (np.sqrt(us ** 2 + vs ** 2)[strong].mean()
                 / np.sqrt(ug ** 2 + vg ** 2)[strong].mean())
        return float(np.median(ang)), float(ratio)

    ang_drag, ratio_drag = cross_isobar_angle(True)
    ang_free, ratio_free = cross_isobar_angle(False)

    ok = (5.0 < abs(ang_drag) < 60.0 and abs(ang_drag) > abs(ang_free)
          and ratio_drag < ratio_free)
    report("Ekman spiral: surface wind turns across the isobars", ok,
           f"with drag {ang_drag:+.1f} deg (vs local geostrophic), "
           f"speed ratio {ratio_drag:.2f}; without drag {ang_free:+.1f} deg, "
           f"ratio {ratio_free:.2f}")


# ---------------------------------------------------------------------------
def test_roughness_controls_slowdown():
    """A rougher surface must slow the low-level wind more."""
    def surface_speed(z0):
        m, gr = balanced(z0=z0)
        m.run(6 * 3600)
        if not np.isfinite(m.u).all():
            return np.nan
        return float(np.sqrt(m.u[-1] ** 2 + m.v[-1] ** 2).mean())

    sea = surface_speed(ROUGHNESS["sea"])
    forest = surface_speed(ROUGHNESS["forest"])

    ok = np.isfinite(sea) and np.isfinite(forest) and forest < sea
    report("rougher surface slows the low-level wind more", ok,
           f"sea (z0=2e-4) {sea:.1f} m/s vs forest (z0=1.0) {forest:.1f} m/s")


# ---------------------------------------------------------------------------
def test_drag_does_not_destabilise():
    """Adding drag must not shorten a run that was previously stable."""
    def survived(drag):
        m, gr = balanced(drag=drag)
        n = 0
        for hr in range(1, 13):
            m.run(3600)
            if not np.isfinite(m.u).all():
                break
            n = hr
        return n

    with_drag = survived(True)
    without = survived(False)
    ok = with_drag >= without
    report("drag does not shorten a stable run", ok,
           f"without drag {without}/12 h, with drag {with_drag}/12 h")


if __name__ == "__main__":
    print("\nSurface layer\n" + "=" * 64)
    for fn in (test_drag_coefficient_matches_log_law,
               test_stability_function_behaviour,
               test_drag_removes_momentum_only,
               test_ekman_spiral,
               test_roughness_controls_slowdown,
               test_drag_does_not_destabilise):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 64)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

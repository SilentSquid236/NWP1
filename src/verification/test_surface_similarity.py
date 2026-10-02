"""
Tests for the surface-layer similarity diagnostics (surface_similarity.py) and the
operator methods that use them.

Run:  python src/verification/test_surface_similarity.py
"""
import sys
from pathlib import Path

import numpy as np
np.seterr(all="ignore")

HERE = Path(__file__).resolve().parent
for p in (HERE, HERE.parent, HERE.parent.parent, HERE.parent / "dynamics"):
    sys.path.insert(0, str(p))

import config
from grid import CGrid
from sigma import SigmaLevels, P0, KAPPA, RD, G0
from sigma_operator import SigmaInterpolator
from surface_similarity import (surface_diagnostics, solve_zeta,
                                richardson_from_zeta, psi_m, psi_h)

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


# 1. Neutral: the log law the model's drag assumes -----------------------------
f, th2, ze, ri = surface_diagnostics(6.0, 290.0, 290.0, 230.0)
expect = np.log(10.0 / 0.1) / np.log(230.0 / 0.1)
report("neutral: the 10 m factor is ln(10/z0)/ln(z1/z0) and theta_2 = theta_g",
       abs(float(f) - expect) < 1e-6 and abs(float(th2) - 290.0) < 1e-9 and abs(float(ze)) < 1e-6,
       f"factor {float(f):.6f} against {expect:.6f}; zeta {float(ze):.1e}")

# 2. The stability solve inverts the Richardson relation -----------------------
ri_in = np.array([-5.0, -1.0, -0.2, -0.01, 0.01, 0.1, 0.5, 2.0, 8.0])
z = solve_zeta(ri_in, 230.0)
back = richardson_from_zeta(z, 230.0)
report("solve_zeta inverts Ri_b(zeta) on -5..8",
       np.max(np.abs(back - ri_in)) < 1e-6 and np.all(np.diff(z) > 0),
       "max |Ri(zeta(Ri)) - Ri| = %.1e; zeta " % np.max(np.abs(back - ri_in))
       + " ".join(f"{x:.3g}" for x in z))

# 3. Stability changes the 10 m wind the right way -----------------------------
dths = np.array([-3.0, -1.0, 0.0, 1.0, 3.0, 6.0])        # theta_1 - theta_g
F = np.array([float(surface_diagnostics(5.0, 290.0, 290.0 - d, 230.0)[0]) for d in dths])
report("10 m factor: above neutral by day, below at night, falling as the inversion grows",
       F[0] > F[1] > F[2] > F[3] > F[4] > F[5] > 0.0,
       "theta_1 - theta_g " + " ".join(f"{d:+.0f}:{x:.3f}" for d, x in zip(dths, F)))

# 4. theta_2 lies between the ground and the lowest level ----------------------
ok, worst = True, ""
for d in dths:
    for U in (1.0, 5.0, 12.0):
        _, t2, _, _ = surface_diagnostics(U, 290.0, 290.0 - d, 230.0)
        lo, hi = min(290.0, 290.0 - d), max(290.0, 290.0 - d)
        if not (lo - 1e-9 <= float(t2) <= hi + 1e-9):
            ok, worst = False, f"d {d} U {U}: {float(t2):.3f} outside [{lo}, {hi}]"
report("theta_2 lies between theta_g and theta_1 for every case", ok, worst or "18 cases")

# 5. Stability functions: zero at neutral, the right signs ---------------------
report("psi_m, psi_h: 0 at neutral, positive unstable, negative stable",
       abs(float(psi_m(0.0))) < 1e-12 and abs(float(psi_h(0.0))) < 1e-12
       and float(psi_m(-1.0)) > 0 and float(psi_h(-1.0)) > 0
       and float(psi_m(1.0)) < 0 and float(psi_h(1.0)) < 0,
       f"psi_m(-1) {float(psi_m(-1.0)):.3f}, psi_m(1) {float(psi_m(1.0)):.3f}, "
       f"psi_h(-1) {float(psi_h(-1.0)):.3f}, psi_h(1) {float(psi_h(1.0)):.3f}")

# 6. Through the operator: a neutral ground gives the neutral factor ------------
D = config.DOMAIN
ny, nx = 40, 44
lat0 = 0.5 * (D["lat_min"] + D["lat_max"])
dy = (D["lat_max"] - D["lat_min"]) * 111_132.0 / ny
dx = (D["lon_max"] - D["lon_min"]) * 111_320.0 * np.cos(np.radians(lat0)) / nx
gr = CGrid(nx, ny, dx, dy, f0=9.81e-5, beta=1.69e-11, edge_mode="replicate")
lev = SigmaLevels(20)
terrain = np.full((ny, nx), 300.0)
pi = 101325.0 * np.exp(-G0 * terrain / (RD * 280.0)) - lev.p_top
op = SigmaInterpolator(gr, D, lev.sigma, lev.p_top, pi, terrain)
p = lev.pressure(pi)
th = np.full(p.shape, 295.0)                              # neutral column
u = np.full(p.shape, 6.0); v = np.zeros(p.shape)
ex = ((lev.p_top + pi) / P0) ** KAPPA
lat, lon = lat0, 0.5 * (D["lon_min"] + D["lon_max"])
d = op.surface_similarity(th, u, v, 295.0 * ex, lat, lon)
fN, agl = op.wind_10m_factor(th, lat, lon)
cold = op.surface_similarity(th, u, v, (295.0 - 4.0) * ex, lat, lon)
T2, info = op.station_temperature_similarity(th, u, v, (295.0 - 4.0) * ex, lat, lon, 300.0)
# Tolerance 1e-4: tg and the surface pressure are interpolated to the point
# separately, so a "neutral" ground comes out at Ri_b ~ 1e-5, not 0, and the
# stable side changes the factor by ~3e-5 there.
report("operator: neutral ground gives the neutral factor; a cold ground lowers it and the 2 m T",
       abs(d["wind_factor"] - fN) < 1e-4 and cold["wind_factor"] < fN
       and T2 < 295.0 * float(ex[0, 0]) and info.get("temp_operator") == "similarity",
       f"neutral {d['wind_factor']:.4f} against log law {fN:.4f}; cold ground {cold['wind_factor']:.4f}; "
       f"2 m T {T2:.2f} K against {295.0 * float(ex[0, 0]):.2f} K at a neutral ground")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

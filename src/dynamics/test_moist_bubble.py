"""
The Bryan and Fritsch (2002) moist benchmark for the NH core (CAM stage S7).

A warm bubble rises through (a) a dry atmosphere of constant theta = 300 K and
(b) a saturated atmosphere with total water q_t = 0.02 that is neutral for a
reversible (no rain) moist ascent. With a consistent formulation the two
bubbles rise alike: Bryan and Fritsch's point is that moist models with
inconsistent thermodynamics (energy or mass not conserved, loading left out)
do not reproduce the dry result.

Here the moist base state is built with the model's own thermodynamics (the
Kessler saturation adjustment of moisture.py, with warm rain off), lifting a
saturated parcel level by level, so it is neutral for this model rather than
for Bryan and Fritsch's more complete one. The perturbation is the same
relative theta_rho perturbation in both, 2 cos^2(pi L / 2) / 300 with radius
2 km at 2 km height; the moist bubble is kept saturated at its q_t.

Gate (S7 plan, written before the run): after 1000 s the moist bubble's
maximum w and the height of its maximum theta_rho perturbation are within
10 % of the dry bubble's.

Domain 20 km x 10 km (lid at 250 hPa), dx = 200 m, 50 sigma levels, 2-D
(ny = 4, periodic), dt = 1.5 s. Run:  python test_moist_bubble.py  (~1-2 min)
"""
import sys
import time

import numpy as np

import nh3d
import moisture as mo
from grid import CGrid
from sigma import SigmaLevels

P0, KAPPA, RD, G, CP = nh3d.P0, nh3d.KAPPA, nh3d.RD, nh3d.G, nh3d.CP
NX, NY, DX, NZ = 100, 4, 200.0, 50
PS, PTOP, QT = 1.0e5, 2.5e4, 0.02
T_SFC_MOIST = 293.0          # saturated surface temperature of the moist base state
lev = SigmaLevels(NZ, p_top=PTOP, stretch=1.0)
backend = sys.argv[1] if len(sys.argv) > 1 else "numpy"


class _Hydro:
    drag = mixing = convection = False
    sponge_levels = 0
    # read by the C physics path even when drag and mixing are off
    z0, theta_surface, ri_crit, k_max, mixing_length = 0.1, None, 0.25, 200.0, 30.0

    def __init__(self, g):
        self.grid, self.lev, self.terrain = g, lev, np.zeros((g.ny, g.nx))
        self._sponge = np.zeros((lev.nz, 1, 1))

    def max_dt(self):
        return 1.5


def adjust(theta, qv, qc, p, n=8):
    """Saturation adjustment only (repeated Newton steps), no rain."""
    one = np.ones_like(np.asarray(theta, dtype=float))
    rho = p / (RD * theta * (p / P0) ** KAPPA)
    for _ in range(n):
        theta, qv, qc, _, _ = mo.kessler(theta, qv, qc, 0 * one, p * one, rho * one, 1.0 * one, 1.0,
                                         warm_rain=False)
    return theta, qv, qc


def moist_base(p_full):
    """theta, qv, qc on the full levels (index 0 = top) from a saturated parcel
    lifted from the surface with q_t fixed, 40 adjustment substeps per layer."""
    th = np.array([T_SFC_MOIST * (P0 / PS) ** KAPPA])
    qv = mo.qsat(np.array([T_SFC_MOIST]), np.array([PS]))
    qc = QT - qv
    out = np.empty((3, p_full.size))
    p_prev = PS
    for k in range(p_full.size - 1, -1, -1):
        for pp in np.linspace(p_prev, p_full[k], 41)[1:]:
            th, qv, qc = adjust(th, qv, qc, np.array([pp]), n=2)
        th, qv, qc = adjust(th, qv, qc, np.array([p_full[k]]), n=8)
        out[:, k] = th[0], qv[0], qc[0]
        p_prev = p_full[k]
    return out


def run(moist):
    g = CGrid(NX, NY, DX, DX, f0=0.0, beta=0.0, edge_mode="periodic")
    m = nh3d.NHModel(_Hydro(g), theta_ref=(lambda z: 300.0 + 0 * z), ns=6, dt_max=1.5,
                     backend=backend)
    m.warm_rain = False
    mu = np.full((NY, NX), PS - PTOP)
    p = PTOP + np.asarray(lev.sigma)[:, None, None] * mu[None]
    # height of the full levels in the base state (hydrostatic, from theta_rho)
    if moist:
        b = moist_base(p[:, 0, 0])
        th_b, qv_b, qc_b = (np.broadcast_to(a[:, None, None], p.shape).copy() for a in b)
    else:
        th_b = np.full(p.shape, 300.0)
        qv_b = qc_b = np.zeros(p.shape)
    thr_b = th_b * mo.rho_factor(qv_b, qc_b)
    m.set_state(np.zeros(p.shape), np.zeros(p.shape), th_b, mu,
                q=({"qv": qv_b, "qc": qc_b} if moist else None))
    phi = m.core.phi
    z = 0.5 * (phi[:-1] + phi[1:]) / G
    x = (np.arange(NX) + 0.5) * DX - NX * DX / 2
    L = np.sqrt((x[None, None, :] / 2000.0) ** 2 + ((z - 2000.0) / 2000.0) ** 2)
    dth = np.where(L <= 1.0, 2.0 * np.cos(np.pi * L / 2) ** 2, 0.0)
    thr_t = thr_b * (1.0 + dth / 300.0)
    if moist:
        # find theta (saturated at q_t) with theta * factor = the target theta_rho
        th = th_b.copy()
        for _ in range(30):
            T = th * (p / P0) ** KAPPA
            qv = np.minimum(mo.qsat(T, p), QT)
            f = mo.rho_factor(qv, QT - qv)
            th = th + (thr_t - th * f) / f * 0.8
        T = th * (p / P0) ** KAPPA
        qv = np.minimum(mo.qsat(T, p), QT)
        qc = QT - qv
        m.set_state(np.zeros(p.shape), np.zeros(p.shape), th, mu, q={"qv": qv, "qc": qc})
    else:
        m.set_state(np.zeros(p.shape), np.zeros(p.shape), thr_t, mu)
    t0 = time.time()
    while m.time < 1000.0 - 1e-6:
        m.step(min(1.5, 1000.0 - m.time))
    c = m.core
    w = c.winds()[2]
    thr = c.Th / c.mu[None]
    zc = 0.5 * (c.phi[:-1] + c.phi[1:]) / G
    dthr = thr - thr_b
    k, j, i = np.unravel_index(np.argmax(dthr), dthr.shape)
    return {"wmax": float(w.max()), "z_top": float(zc[k, j, i]), "dthr_max": float(dthr.max()),
            "secs": time.time() - t0, "qc_min": float(m.q("qc").min()) if moist else 0.0,
            "base_surface": (float(th_b[-1, 0, 0]), float(qv_b[-1, 0, 0]), float(qc_b[-1, 0, 0])),
            "base_top": (float(th_b[0, 0, 0]), float(qv_b[0, 0, 0]), float(qc_b[0, 0, 0]))}


if __name__ == "__main__":
    dry = run(False)
    wet = run(True)
    for name, r in (("dry", dry), ("moist", wet)):
        print(f"  {name:5s}: max w {r['wmax']:.2f} m/s, max theta_rho' {r['dthr_max']:.2f} K at "
              f"{r['z_top']:.0f} m, {r['secs']:.0f} s; base surface {r['base_surface']}, "
              f"top {r['base_top']}")
    rw = wet["wmax"] / dry["wmax"] - 1
    rz = wet["z_top"] / dry["z_top"] - 1
    ok = abs(rw) <= 0.10 and abs(rz) <= 0.10
    print(f"  [{'PASS' if ok else 'FAIL'}] moist vs dry after 1000 s: max w {rw:+.1%}, "
          f"height of max theta_rho' {rz:+.1%} (gate: both within 10 %)")
    raise SystemExit(0 if ok else 1)

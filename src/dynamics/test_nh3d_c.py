"""
The compiled kernels of the NH core against the NumPy reference (CAM stage S4a).

Run:  python test_nh3d_c.py

Needs a C compiler (gcc). Without one every test is reported as SKIPPED and
the file exits 0: the NumPy core is unaffected, and the server, where the CAM
runs, has gcc.
"""
import numpy as np

import cnh
import nh3d
from grid import CGrid
from sigma import SigmaLevels

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def case(backend, edge, nx=48, ny=40, nz=20, steps=6, ns=6, threads=0, wind=10.0):
    g = CGrid(nx, ny, 3000.0, 3000.0, f0=1e-4, beta=0.0, edge_mode=edge)
    X, Y = g.Xc - g.Lx / 2, g.Yc - g.Ly / 2
    m = nh3d.NH3D(g, SigmaLevels(nz), terrain=1200.0 * np.exp(-(X ** 2 + Y ** 2) / 30e3 ** 2),
                  ns=ns, backend=backend, threads=threads)
    th = m.Th / m.mu[None]
    th = th + 3.0 * np.exp(-((X + 25e3) ** 2 + (Y - 10e3) ** 2) / 15e3 ** 2)[None]
    m.Th = th * m.mu[None]
    m.U = m.U + wind * m.h_to_u(m.mu)[None]
    m.V = m.V + 0.5 * wind * m.h_to_v(m.mu)[None]
    for _ in range(steps):
        m.step(20.0)
    return m


def rel_diff(a, b):
    out = {}
    for k in ("mu", "U", "V", "W", "Th", "phi"):
        x, y = getattr(a, k), getattr(b, k)
        out[k] = float(np.abs(x - y).max() / max(np.abs(x).max(), 1e-30))
    return out


if not cnh.available():
    try:
        cnh.load("f64")
    except cnh.Unavailable as e:
        why = str(e).splitlines()[0]
    print(f"  [SKIPPED] no compiled kernels on this machine: {why}")
    raise SystemExit(0)

# 0. The tendencies alone, from the same disturbed state, both edge modes
for edge in ("periodic", "replicate"):
    a = case("numpy", edge, steps=2)
    b = case("c", edge, steps=0)
    for k in ("mu", "U", "V", "W", "Th", "phi"):
        setattr(b, k, getattr(a, k).copy())
    Fa = a.tendencies(a.mu, a.U, a.V, a.W, a.Th, a.phi)
    Fb = b.tendencies(b.mu, b.U, b.V, b.W, b.Th, b.phi)
    names = ("dmu", "FU", "FV", "FW", "FTh", "Fphi")
    d = {n: float(np.abs(x - y).max() / max(np.abs(x).max(), 1e-30)) for n, x, y in zip(names, Fa, Fb)}
    report(f"C tendencies match NumPy, {edge} edges (max relative difference < 1e-11)",
           max(d.values()) < 1e-11, ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 1. Periodic edges, bubble over terrain in a mean wind, 6 steps (36 substeps)
a, b = case("numpy", "periodic"), case("c", "periodic")
d = rel_diff(a, b)
report("C step (tendencies, set-up, substeps) matches NumPy over 6 steps, periodic edges (< 1e-11)",
       max(d.values()) < 1e-11 and np.abs(b.W).max() > 1e-3,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; max|w| {np.abs(b.W).max() / b.mu.max():.2e}")

# 2. Replicate (limited-area) edges, more substeps
a, b = case("numpy", "replicate", ns=12), case("c", "replicate", ns=12)
d = rel_diff(a, b)
report("C step matches NumPy over 6 steps, replicate edges, ns = 12",
       max(d.values()) < 1e-11,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 3. Thread count does not change the answer (columns are independent)
a, b = case("c", "replicate", threads=1), case("c", "replicate", threads=4)
d = rel_diff(a, b)
report("C results are bit-identical on 1 and 4 threads",
       max(d.values()) == 0.0, ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 4. A wrong dtype is refused, not silently converted
try:
    k = cnh.load("f64")
    z = np.zeros((2, 3, 4), dtype=np.float32)
    k.ac_pd(z, z, z, np.zeros((3, 3, 4), np.float32), z, 0.1, z)
    ok, why = False, "float32 arrays were accepted by the float64 kernel"
except TypeError as e:
    ok, why = True, f"refused: {str(e)[:70]}"
report("a float32 array passed to a float64 kernel is refused", ok, why)

# 5. Physics (drag + Richardson mixing + sponge) in C equals the NumPy physics
import turbulence


class _Hy:
    def __init__(self, g, nz):
        self.grid, self.lev = g, SigmaLevels(nz)
        X, Y = g.Xc - g.Lx / 2, g.Yc - g.Ly / 2
        self.terrain = 800.0 * np.exp(-(X ** 2 + Y ** 2) / 30e3 ** 2)
        self.drag, self.mixing, self.convection = True, True, False
        self.z0 = np.where(X > 0, 1.0, 0.0002)
        self.theta_surface = None
        self.ri_crit, self.k_max, self.mixing_length = turbulence.RI_CRIT, turbulence.K_MAX, turbulence.MIXING_LENGTH
        self.sponge_levels = 3
        self._sponge = np.zeros((nz, 1, 1)); self._sponge[:3] = 1e-3
    def max_dt(self):
        return 20.0


def nh_pair(edge):
    g = CGrid(40, 36, 3000.0, 3000.0, f0=1e-4, beta=0.0, edge_mode=edge)
    out = []
    for be in ("numpy", "c"):
        hy = _Hy(g, 20)
        m = nh3d.NHModel(hy, backend=be, threads=4)
        rng = np.random.default_rng(5)
        mu = 9.0e4 * np.ones((36, 40)) - 50 * hy.terrain
        th = 290 + 30 * (1 - hy.lev.sigma)[:, None, None] + rng.normal(0, 0.8, (20, 36, 40))
        u = 8 + 6 * (1 - hy.lev.sigma)[:, None, None] + rng.normal(0, 2, (20, 36, 40))
        v = rng.normal(0, 3, (20, 36, 40))
        m.set_state(u, v, th, mu)
        m._u_ref, m._v_ref = 0.5 * u, 0.5 * v
        out.append(m)
    return out


for edge in ("periodic", "replicate"):
    a, b = nh_pair(edge)
    Pa = a._physics(a.core.mu, a.core.U, a.core.V, a.core.W, a.core.Th, a.core.phi)
    Pb = b._physics(b.core.mu, b.core.U, b.core.V, b.core.W, b.core.Th, b.core.phi)
    d = {k: float(np.abs(Pa[k] - Pb[k]).max() / max(np.abs(Pa[k]).max(), 1e-30)) for k in ("U", "V", "Th", "W")}
    report(f"C physics (drag, mixing, sponge, w damping) matches NumPy, {edge} edges (< 1e-11)",
           max(d.values()) < 1e-11 and np.abs(Pa["Th"]).max() > 0,
           ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 6. Zone-only relaxation with the time interpolation in C equals relax_with(driver.at(t))
from boundaries import BoundaryDriver


class _Relax:
    def __init__(self, g, w=6):
        j, i = np.meshgrid(np.arange(g.ny), np.arange(g.nx), indexing="ij")
        dd = np.minimum(np.minimum(j, g.ny - 1 - j), np.minimum(i, g.nx - 1 - i))
        self.alpha2d = np.where(dd < w, 0.1 * (1 - dd / w), 0.0)


a, b = nh_pair("replicate")
rng = np.random.default_rng(9)
fr = []
for _ in range(2):
    fr.append({"u": a.u + rng.normal(0, 2, a.u.shape), "v": a.v + rng.normal(0, 2, a.v.shape),
               "theta": a.theta + rng.normal(0, 1, a.theta.shape), "pi": a.pi + rng.normal(0, 300, a.pi.shape)})
drv = BoundaryDriver([0.0, 3600.0], fr)
rl = _Relax(a.grid)
for m in (a, b):
    m._last_dt = 20.0
a.relax_with(rl, drv.at(1234.0))
b.relax_with_driver(rl, drv, 1234.0)
d = {k: float(np.abs(getattr(a.core, k) - getattr(b.core, k)).max() / max(np.abs(getattr(a.core, k)).max(), 1e-30))
     for k in ("mu", "U", "V", "W", "Th", "phi")}
report("zone-only C relaxation (time interpolation inside) matches relax_with(driver.at(t)) (< 1e-12)",
       max(d.values()) < 1e-12, ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 7. CAM stage S5e: blocked column kernels are bit-identical to the column ones
#    (nx = 45 is one block of 32 plus a partial block of 13)
import os

def _with_layout(layout, fn):
    old = os.environ.get("NWP_NH_LAYOUT")
    os.environ["NWP_NH_LAYOUT"] = layout
    try:
        return fn()
    finally:
        if old is None:
            os.environ.pop("NWP_NH_LAYOUT", None)
        else:
            os.environ["NWP_NH_LAYOUT"] = old

for edge in ("periodic", "replicate"):
    a = _with_layout("column", lambda: case("c", edge, nx=45, ny=31, steps=4, threads=4))
    b = _with_layout("block", lambda: case("c", edge, nx=45, ny=31, steps=4, threads=4))
    d = rel_diff(a, b)
    report(f"blocked column kernels are bit-identical to the S4 column kernels, {edge} edges",
           max(d.values()) == 0.0 and a._ck.layout == "column" and b._ck.layout == "block",
           ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

# 8. PAV convective adjustment in C equals the NumPy PAV (with the theta/u/v round trip)
from convection import dry_convective_adjustment_pav

a, b = nh_pair("replicate")
ca, cb_ = a.core, b.core
th, u, v, info = dry_convective_adjustment_pav(a.theta, a.u, a.v, ca.mu, a.lev, mix_momentum=True)
Th_ref, U_ref, V_ref = th * ca.mu[None], u * ca.h_to_u(ca.mu)[None], v * ca.h_to_v(ca.mu)[None]
ncol = cb_._ck.pav_col(cb_._xm1, cb_._ym1, 1e-10, True,
                       np.ascontiguousarray(b.lev.dsigma, dtype=float), cb_.mu, cb_.U, cb_.V, cb_.Th)
d = {k: float(np.abs(x - y).max() / max(np.abs(x).max(), 1e-30))
     for k, x, y in (("Th", Th_ref, cb_.Th), ("U", U_ref, cb_.U), ("V", V_ref, cb_.V))}
report("C PAV adjustment equals the NumPy PAV, momentum mixed (< 1e-14), same columns",
       max(d.values()) < 1e-14 and ncol == info["columns"] and ncol > 0,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; columns C {ncol}, NumPy {info['columns']}")

# 9. The largest wind for max_dt, in C, equals NumPy's
au, av = cb_._ck.uv_maxabs(cb_._xm1, cb_._ym1, cb_.mu, cb_.U, cb_.V)
nu, nv = float(np.abs(b.u).max()), float(np.abs(b.v).max())
report("C max |u|, |v| equal NumPy's exactly", au == nu and av == nv,
       f"u {au:.6f} vs {nu:.6f}, v {av:.6f} vs {nv:.6f}")

# 14. CAM stage S6: drag with the land surface's theta_g (Louis 1979) in C equals NumPy
for edge in ("periodic", "replicate"):
    a, b = nh_pair(edge)
    rng = np.random.default_rng(11)
    th1 = (a.core.Th / a.core.mu[None])[-1]
    ts = th1 + rng.normal(0, 3.0, th1.shape)          # stable and unstable columns
    for m in (a, b):
        m.land_surface = object()                      # marks "land surface on" (Louis drag)
        m.h.theta_surface = ts
    Pa = a._physics(a.core.mu, a.core.U, a.core.V, a.core.W, a.core.Th, a.core.phi)
    Pb = b._physics(b.core.mu, b.core.U, b.core.V, b.core.W, b.core.Th, b.core.phi)
    d = {k: float(np.abs(Pa[k] - Pb[k]).max() / max(np.abs(Pa[k]).max(), 1e-30)) for k in ("U", "V", "Th")}
    nst, nun = int((ts < th1).sum()), int((ts > th1).sum())
    report(f"C drag with Louis stability (theta_g) matches NumPy surface_drag(stability='louis'), {edge} (< 1e-11)",
           max(d.values()) < 1e-11 and nst > 0 and nun > 0,
           ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; stable cols {nst}, unstable {nun}")

# P-76: the edge relaxation with pi_mode "mean" and "none" in C equals relax_with
for mode in ("mean", "none"):
    a, b = nh_pair("replicate")
    rng = np.random.default_rng(9)
    fr = []
    for _ in range(2):
        fr.append({"u": a.u + rng.normal(0, 2, a.u.shape), "v": a.v + rng.normal(0, 2, a.v.shape),
                   "theta": a.theta + rng.normal(0, 1, a.theta.shape), "pi": a.pi + rng.normal(300, 300, a.pi.shape)})
    drv = BoundaryDriver([0.0, 3600.0], fr)
    rl = _Relax(a.grid); rl.pi_mode = mode
    for m in (a, b):
        m._last_dt = 20.0
    mu0 = a.core.mu.copy()
    a.relax_with(rl, drv.at(1234.0))
    b.relax_with_driver(rl, drv, 1234.0)
    d = {k: float(np.abs(getattr(a.core, k) - getattr(b.core, k)).max() / max(np.abs(getattr(a.core, k)).max(), 1e-30))
         for k in ("mu", "U", "V", "W", "Th", "phi")}
    dmu = a.core.mu - mu0
    zone = rl.alpha2d > 0
    shape_ok = (np.abs(dmu).max() == 0) if mode == "none" else \
        (np.ptp(dmu[zone] / rl.alpha2d[zone]) < 1e-6 * np.abs(dmu).max() / rl.alpha2d.max() + 1e-12 and np.abs(dmu).max() > 0)
    report(f"P-76 edge relaxation, surface pressure '{mode}': C matches relax_with (< 1e-12); mu change "
           + ("zero" if mode == "none" else "uniform per unit weight"),
           max(d.values()) < 1e-12 and shape_ok, ", ".join(f"{k} {v:.1e}" for k, v in d.items())
           + f"; max |dmu| {np.abs(dmu).max():.2f} Pa")

# P-79 / E6: a narrower zone for theta and q. C (target trick) equals relax_with (NumPy)
a, b = nh_pair("replicate")
rng = np.random.default_rng(19)
fr = []
for _ in range(2):
    fr.append({"u": a.u + rng.normal(0, 2, a.u.shape), "v": a.v + rng.normal(0, 2, a.v.shape),
               "theta": a.theta + rng.normal(0, 1, a.theta.shape), "pi": a.pi + rng.normal(300, 300, a.pi.shape)})
drv = BoundaryDriver([0.0, 3600.0], fr)
rl = _Relax(a.grid)
rl.alpha2d_thermo = _Relax(a.grid, w=3).alpha2d          # inside the dynamic zone (w = 6) and never larger
for m in (a, b):
    m._last_dt = 20.0
th0 = a.theta.copy()
a.relax_with(rl, drv.at(1234.0))
b.relax_with_driver(rl, drv, 1234.0)
d = {k: float(np.abs(getattr(a.core, k) - getattr(b.core, k)).max() / max(np.abs(getattr(a.core, k)).max(), 1e-30))
     for k in ("mu", "U", "V", "W", "Th", "phi")}
outer_dyn_only = (rl.alpha2d > 0) & (rl.alpha2d_thermo == 0)
dth = np.abs(a.theta - th0)[:, outer_dyn_only].max()
report("P-79 separate theta/q zone (2-D): C matches relax_with (< 1e-12); theta unchanged where only the dynamics are relaxed",
       max(d.values()) < 1e-12 and dth < 1e-12 and outer_dyn_only.any(),
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; max |dtheta| outside the theta zone {dth:.1e} K "
       f"({int(outer_dyn_only.sum())} columns)")

# P-79 / E6c: per-level theta/q weights (narrow near the ground, wide aloft): C equals NumPy
a, b = nh_pair("replicate")
rl = _Relax(a.grid)
nzl = a.lev.nz
rl.alpha2d_thermo = _Relax(a.grid, w=3).alpha2d
rl.thermo_levels = np.linspace(0.0, 1.0, nzl)                  # 0 at the top (wide), 1 at the ground (narrow)
for m in (a, b):
    m._last_dt = 20.0
th0 = a.theta.copy()
a.relax_with(rl, drv.at(2000.0))
b.relax_with_driver(rl, drv, 2000.0)
d = {k: float(np.abs(getattr(a.core, k) - getattr(b.core, k)).max() / max(np.abs(getattr(a.core, k)).max(), 1e-30))
     for k in ("mu", "U", "V", "W", "Th", "phi")}
dth_low = np.abs(a.theta - th0)[-1][outer_dyn_only].max()
dth_top = np.abs(a.theta - th0)[0][outer_dyn_only].max()
report("P-79 per-level theta zone (3-D): C matches relax_with (< 1e-12); lowest level untouched outside the narrow zone, top level relaxed",
       max(d.values()) < 1e-12 and dth_low < 1e-12 and dth_top > 1e-6,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; |dtheta| lowest {dth_low:.1e}, top {dth_top:.2e} K")

# 15-19. CAM stage S7: moisture in C against the NumPy reference
import moisture as mo


def moist_pair(edge, rain=False):
    a, b = nh_pair(edge)
    rng = np.random.default_rng(21)
    p = a.lev.p_top + np.asarray(a.lev.sigma)[:, None, None] * a.core.mu[None]
    T = a.theta * (p / nh3d.P0) ** nh3d.KAPPA
    qv = mo.qv_from_rh(np.clip(80 + rng.normal(0, 15, T.shape), 20, 100), T, p)
    qv = qv * np.where(rng.random(T.shape) < 0.2, 1.03, 1.0)          # some supersaturation
    qc = np.where(rng.random(T.shape) < 0.3, 1.5e-3 * rng.random(T.shape), 0.0)
    qr = np.where(rng.random(T.shape) < 0.3, 2e-3 * rng.random(T.shape), 0.0) if rain else 0 * qc
    qr[:, :2] = -1e-5                                                   # negatives to fill
    # one copy of the state for both: a.theta read after a's moisture is on
    # (with the negative qr below) would hand b a different theta
    u0, v0, th0, mu0 = a.u.copy(), a.v.copy(), a.theta.copy(), a.core.mu.copy()
    for m in (a, b):
        m.set_state(u0, v0, th0, mu0, q={"qv": qv, "qc": qc, "qr": qr})
        m.core.Q[2] = qr * m.core.mu[None]                             # keep the negatives
        m._u_ref, m._v_ref = 0.5 * u0, 0.5 * v0
    return a, b


for edge in ("periodic", "replicate"):
    a, b = moist_pair(edge)
    Pa = a._physics(a.core.mu, a.core.U, a.core.V, a.core.W, a.core.Th, a.core.phi)
    Pb = b._physics(b.core.mu, b.core.U, b.core.V, b.core.W, b.core.Th, b.core.phi)
    d = {k: float(np.abs(Pa[k] - Pb[k]).max() / max(np.abs(Pa[k]).max(), 1e-30)) for k in ("U", "V", "Th", "W", "Q")}
    report(f"S7 moist physics in C (dry theta in, q mixed with theta's K) matches NumPy, {edge} (< 1e-11)",
           max(d.values()) < 1e-11 and np.abs(Pa["Q"]).max() > 0,
           ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

for edge in ("periodic", "replicate"):
    a, b = moist_pair(edge, rain=True)
    q = np.ascontiguousarray(a.q())
    Fa = a.core.__class__._adv_q.__get__(a.core)(q, a.core.U, a.core.V)
    Fb = b.core._adv_q(q, np.ascontiguousarray(b.core.U), np.ascontiguousarray(b.core.V))
    d = float(np.abs(Fa - Fb).max() / np.abs(Fa).max())
    report(f"S7 moisture advection adv_q in C matches NH3D._adv_q, {edge} (< 1e-12)",
           d < 1e-12 and a.core._ck is None, f"max rel diff {d:.1e}")

for edge in ("periodic", "replicate"):
    a, b = moist_pair(edge, rain=True)
    rng = np.random.default_rng(31)
    ex = rng.normal(0, 1e-3, a.core.Q.shape)
    Qs = a.core.Q * (1 + 0.01 * rng.random(a.core.Q.shape))
    outs = []
    for m in (a, b):
        m.core._flux_avg = (np.ascontiguousarray(a.core.U), np.ascontiguousarray(a.core.V))
        m.core._extra_last = {"Q": ex}
        outs.append(m.core._q_stage(a.core.Q, Qs, a.core.mu * 1.001, 7.0))
    d = float(np.abs(outs[0] - outs[1]).max() / np.abs(outs[0]).max())
    report(f"S7 fused RK-stage Q update q_stage in C matches NH3D._q_stage, {edge} (< 1e-12)",
           d < 1e-12, f"max rel diff {d:.1e}")

a, b = moist_pair("replicate", rain=True)
ca, cb_ = a.core, b.core
qa = list(ca.Q / ca.mu[None, None])
th, u, v, info = dry_convective_adjustment_pav(ca.Th / ca.mu[None], a.u, a.v, ca.mu, a.lev,
                                               mix_momentum=True, scalars=qa)
Q_ref = np.stack(info["scalars"]) * ca.mu[None, None]
Qc_ = np.ascontiguousarray(cb_.Q)
ncol = cb_._ck.pav_col(cb_._xm1, cb_._ym1, 1e-10, True, np.ascontiguousarray(b.lev.dsigma, dtype=float),
                       cb_.mu, cb_.U, cb_.V, cb_.Th, Qm=Qc_)
d = {k: float(np.abs(x - y).max() / max(np.abs(x).max(), 1e-30))
     for k, x, y in (("Th", th * ca.mu[None], cb_.Th), ("Q", Q_ref, Qc_))}
report("S7 PAV in C mixes q over the same blocks as the NumPy PAV (< 1e-14), same columns",
       max(d.values()) < 1e-14 and ncol == info["columns"] and ncol > 0,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; columns {ncol}")

for warm in (True, False):
    a, b = moist_pair("replicate", rain=True)
    for m in (a, b):
        m.warm_rain = warm
        m._step_moisture(20.0)
    dTh = float(np.abs(a.core.Th - b.core.Th).max() / np.abs(a.core.Th).max())
    dQ = float(np.abs(a.core.Q - b.core.Q).max() / np.abs(a.core.Q).max())
    dr = float(np.abs(a.rain_step - b.rain_step).max() / max(np.abs(a.rain_step).max(), 1e-30))
    report(f"S7 Kessler step moist_col in C matches NHModel._step_moisture (warm_rain={warm}; existing rain falls either way) (< 1e-11)",
           max(dTh, dQ, dr) < 1e-11 and abs(a.q_lost - b.q_lost) <= 1e-9 * max(abs(a.q_lost), 1e-30)
           and a.rain_step.max() > 0,
           f"Th {dTh:.1e}, Q {dQ:.1e}, rain {dr:.1e} (max {a.rain_step.max():.3f} mm); "
           f"hole-filling loss {a.q_lost:.3e} vs {b.q_lost:.3e} kg/m^2")

# P-81 simple ice: the C kernel matches moisture.kessler(ice=True), and ice changes the step
for warm in (True, False):
    a, b = moist_pair("replicate", rain=True)
    a0, _ = moist_pair("replicate", rain=True)
    for m in (a, b):
        m.warm_rain = warm
        m.ice = True
        m._step_moisture(20.0)
    a0.warm_rain = warm
    a0._step_moisture(20.0)
    dTh = float(np.abs(a.core.Th - b.core.Th).max() / np.abs(a.core.Th).max())
    dQ = float(np.abs(a.core.Q - b.core.Q).max() / np.abs(a.core.Q).max())
    dr = float(np.abs(a.rain_step - b.rain_step).max() / max(np.abs(a.rain_step).max(), 1e-30))
    eff = float(np.abs(a.core.Q - a0.core.Q).max() / np.abs(a0.core.Q).max())
    report(f"P-81 simple ice: moist_col in C matches NHModel._step_moisture with ice (warm_rain={warm}) (< 1e-11), "
           f"and ice changes the step (> 1e-6)",
           max(dTh, dQ, dr) < 1e-11 and eff > 1e-6,
           f"Th {dTh:.1e}, Q {dQ:.1e}, rain {dr:.1e}; ice vs no ice {eff:.1e}")

a, b = moist_pair("replicate", rain=True)
rng = np.random.default_rng(9)
fr = []
for _ in range(2):
    fr.append({"u": a.u + rng.normal(0, 2, a.u.shape), "v": a.v + rng.normal(0, 2, a.v.shape),
               "theta": a.theta + rng.normal(0, 1, a.theta.shape), "pi": a.pi + rng.normal(0, 300, a.pi.shape),
               "qv": a.q("qv") * (1 + 0.2 * rng.random(a.theta.shape)), "qc": 0 * a.theta,
               "qr": a.q("qr").clip(0)})
drv = BoundaryDriver([0.0, 3600.0], fr)
rl = _Relax(a.grid)
for m in (a, b):
    m._last_dt = 20.0
a.relax_with(rl, drv.at(1234.0))
b.relax_with_driver(rl, drv, 1234.0)
d = {k: float(np.abs(getattr(a.core, k) - getattr(b.core, k)).max() / max(np.abs(getattr(a.core, k)).max(), 1e-30))
     for k in ("mu", "U", "V", "W", "Th", "phi", "Q")}
report("S7 moist zone relaxation in C matches relax_with: q exactly (< 1e-12); Th, phi to the "
       "theta_rho interpolation order (< 1e-5)",
       d["Q"] < 1e-12 and d["mu"] < 1e-12 and d["U"] < 1e-12 and d["Th"] < 1e-5 and d["phi"] < 1e-5,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()))

a, b = moist_pair("replicate", rain=True)
for m in (a, b):
    m.land_surface = None
    for _ in range(5):
        m.step(20.0)
d = {k: float(np.abs(getattr(a.core, k) - getattr(b.core, k)).max() / max(np.abs(getattr(a.core, k)).max(), 1e-30))
     for k in ("mu", "U", "V", "W", "Th", "phi", "Q")}
report("S7 five full moist NHModel steps (dynamics, q transport, physics, Kessler): C vs NumPy (< 1e-9)",
       max(d.values()) < 1e-9 and np.isfinite(b.core.Q).all() and b.rain.max() > 0,
       ", ".join(f"{k} {v:.1e}" for k, v in d.items()) + f"; rain max {b.rain.max():.3f} mm")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

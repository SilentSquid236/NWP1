"""
Tests for warm-rain moisture (moisture.py and the moist NH core; CAM stage S7).

Run:  python test_moisture.py      (about 1 min; the bubble pair is in
      test_moist_bubble.py)
"""
import numpy as np

import nh3d
import moisture as mo
from grid import CGrid
from sigma import SigmaLevels

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


G = nh3d.G
lev = SigmaLevels(20)

# 1. Saturation, RH and dewpoint ---------------------------------------------
es20 = float(mo.es_bolton(293.15))
qs = float(mo.qsat(293.15, 1.0e5))
td = float(mo.dewpoint(mo.qv_from_rh(50.0, 293.15, 1.0e5), 1.0e5))
td100 = float(mo.dewpoint(mo.qsat(280.0, 8.5e4), 8.5e4))
report("Bolton es(20 C) = 2337 Pa within 2 Pa; qs(20 C, 1000 hPa) 14.9 g/kg; dewpoint inverts es; "
       "Td at 50 % RH, 20 C is 9.3 C",
       abs(es20 - 2337.0) < 2.0 and abs(qs - 0.01488) < 1e-4 and abs(td100 - 280.0) < 1e-9
       and abs(td - 282.43) < 0.1,
       f"es {es20:.1f} Pa, qs {qs * 1e3:.3f} g/kg, Td(qs(280 K)) {td100:.9f} K, Td(50 %) {td - 273.15:.2f} C")

# 2. Saturation adjustment conserves water and moist enthalpy ------------------
nz = 5
p = np.linspace(5.0e4, 9.5e4, nz)[:, None]
pii = (p / nh3d.P0) ** nh3d.KAPPA
T = np.linspace(255.0, 290.0, nz)[:, None]
th0 = T / pii
qv0 = 1.05 * mo.qsat(T, p)                    # 5 % supersaturated
qc0 = np.zeros_like(qv0); qr0 = np.zeros_like(qv0)
rho = p / (nh3d.RD * T)
dz = np.full_like(p, 500.0)
th, qv, qc, qr, rain = mo.kessler(th0, qv0, qc0, qr0, p, rho, dz, 1.0)
Tn = th * pii
ss = np.abs(qv / mo.qsat(Tn, p) - 1).max()
dh = np.abs(nh3d.CP * (Tn - T) - mo.LV * (qv0 - qv)).max() / (mo.LV * (qv0 - qv)).max()
report("a 5 % supersaturated column adjusts to within 0.5 % of saturation in one step; water and "
       "cp T + Lv qv conserved",
       ss < 5e-3 and np.abs(qv + qc + qr - qv0).max() < 1e-15 and dh < 1e-12 and rain.max() == 0,
       f"max |qv/qs - 1| {ss:.1e}; water change {np.abs(qv + qc + qr - qv0).max():.1e}; "
       f"enthalpy residual {dh:.1e}; condensed {(qc.max()) * 1e3:.2f} g/kg")

# 3. Rain falls out: column water budget closes exactly ------------------------
nz = 30
p = np.linspace(2.5e4, 1.0e5, nz)[:, None] * np.ones((1, 3))
T = 300.0 * (p / 1e5) ** nh3d.KAPPA
pii = (p / nh3d.P0) ** nh3d.KAPPA
rho = p / (nh3d.RD * T)
dz = np.full_like(p, 250.0); dz[-5:] = 40.0                # thin layers near the ground: substeps
qv = mo.qsat(T, p)                                           # saturated: no evaporation
qc = np.zeros_like(p)
qr = np.zeros_like(p); qr[5:12] = 4e-3
W0 = (rho * dz * (qv + qc + qr)).sum(axis=0)
tot_rain = 0.0
for _ in range(60):
    th_, qv, qc, qr, rain = mo.kessler(T / pii, qv, qc, qr, p, rho, dz, 20.0)
    T = th_ * pii
    tot_rain = tot_rain + rain
W1 = (rho * dz * (qv + qc + qr)).sum(axis=0) + tot_rain
report("rain falling through 30 layers (40 m near the ground, fall substeps) for 20 min: column "
       "water + surface rain conserved, qr >= 0",
       np.abs(W1 / W0 - 1).max() < 1e-3 and qr.min() >= 0 and tot_rain.min() > 1.0,
       f"budget error {np.abs(W1 / W0 - 1).max():.1e} (saturation adjustment changes qv); "
       f"rain {tot_rain[0]:.2f} mm; min qr {qr.min():.1e}")

# 4. Moisture advection: a uniform q stays uniform; total conserved -----------
gr = CGrid(24, 20, 3000.0, 3000.0, f0=1e-4, beta=0.0, edge_mode="periodic")
X, Y = gr.Xc - gr.Lx / 2, gr.Yc - gr.Ly / 2
m = nh3d.NH3D(gr, lev, terrain=800.0 * np.exp(-(X ** 2 + Y ** 2) / 1.5e4 ** 2))
rng = np.random.default_rng(1)
m.U = m.h_to_u(m.mu)[None] * (10.0 + rng.normal(0, 2, m.U.shape))
m.V = m.h_to_v(m.mu)[None] * (-5.0 + rng.normal(0, 2, m.V.shape))
blob = np.exp(-(X ** 2 + (Y - 1.5e4) ** 2) / 1.0e4 ** 2)[None] * np.ones((lev.nz, 1, 1))
m.q_names = mo.NAMES
m.Q = np.stack([0.01 * np.ones_like(m.Th), 1e-3 * blob, np.zeros_like(m.Th)]) * m.mu[None, None]
tot0 = (m.Q[1] * m.ds[:, None, None]).sum()
for _ in range(30):
    m.step(10.0)
q = m.Q / m.mu[None, None]
dev = np.abs(q[0] / 0.01 - 1).max()
dtot = abs((m.Q[1] * m.ds[:, None, None]).sum() / tot0 - 1)
report("advection over 800 m terrain in a noisy 3-D flow (30 steps): uniform qv stays uniform to "
       "1e-12; a blob's total is conserved to 1e-12; qr stays 0",
       dev < 1e-12 and dtot < 1e-12 and np.abs(q[2]).max() == 0,
       f"max |qv/q0 - 1| {dev:.1e}; blob total change {dtot:.1e}; min blob q {q[1].min():.1e} "
       f"(negative overshoot, filled after the step in NHModel)")

# 5. Hole filling conserves the column ----------------------------------------
qq = np.array([[1e-3, 2e-3], [-2e-4, 1e-3], [5e-4, -3e-3]])     # (level, column)
w = np.array([0.2, 0.3, 0.5])[:, None]
qf, lost = mo.fill_negative(qq, w)
m_before, m_after = (qq[:, 0] * w[:, 0]).sum(), (qf[:, 0] * w[:, 0]).sum()
owed = -(qq[:, 1] * w[:, 0]).sum()
report("hole filling: negatives removed, column mass kept where it can be, owed mass reported",
       qf.min() >= 0 and abs(m_after - m_before) < 1e-18 and abs(lost[1] - owed) < 1e-18
       and np.all(qf[:, 1] == 0),
       f"column 1 mass {m_before:.3e} -> {m_after:.3e}; column 2 (net negative) lost {lost[1]:.2e}")


# 6. A moist atmosphere at rest stays at rest; dry theta round-trips -----------
class _Hydro:
    drag = mixing = convection = False
    sponge_levels = 0

    def __init__(self, g):
        self.grid, self.lev, self.terrain = g, lev, np.zeros((g.ny, g.nx))
        self._sponge = np.zeros((lev.nz, 1, 1))

    def max_dt(self):
        return 20.0


g6 = CGrid(16, 4, 3000.0, 3000.0, f0=0.0, beta=0.0, edge_mode="periodic")
h6 = _Hydro(g6)
m6 = nh3d.NHModel(h6)
mu6 = np.full((4, 16), 1.0e5 - lev.p_top)
p6 = lev.p_top + np.asarray(lev.sigma)[:, None, None] * mu6[None]
th6 = 300.0 * (1.0e5 / p6) ** 0 * (p6 / 1.0e5) ** -0.03
T6 = th6 * (p6 / nh3d.P0) ** nh3d.KAPPA
qv6 = mo.qv_from_rh(70.0, T6, p6)
z6 = np.zeros((lev.nz, 4, 16))
m6.set_state(z6, z6, th6, mu6, q={"qv": qv6})
rt = np.abs(m6.theta - th6).max()
for _ in range(90):
    m6.core.step(20.0)
w6 = np.abs(m6.core.winds()[2]).max()
report("a moist (70 % RH) atmosphere at rest: theta round-trips through theta_rho to 1e-12 K; "
       "|w| stays below 1 cm/s for 30 min",
       rt < 1e-12 and w6 < 0.01,
       f"theta round trip {rt:.1e} K; max|w| {w6:.1e} m/s; qv 0-level {qv6[-1].max() * 1e3:.1f} g/kg")

# 7. Moist NHModel step: Kessler wired in, no rain from a sub-saturated column --
m6.step(20.0)
qn = m6.q()
report("an NHModel step with moisture: sub-saturated air makes no cloud or rain; qv unchanged to 1e-12",
       np.abs(qn[1]).max() == 0 and np.abs(qn[2]).max() == 0 and np.abs(qn[0] / qv6 - 1).max() < 1e-12
       and m6.rain.max() == 0,
       f"max qc {np.abs(qn[1]).max():.1e}, max qr {np.abs(qn[2]).max():.1e}, "
       f"qv change {np.abs(qn[0] / qv6 - 1).max():.1e}")

# 8. The PAV adjustment mixes q over the blocks it mixes theta (P-74) -------
from convection import dry_convective_adjustment_pav
th8 = np.array([310.0, 305.0, 300.0, 302.0, 301.0])[:, None, None]      # lid .. ground; levels 2-3 unstable
qv8 = np.array([1e-3, 3e-3, 6e-3, 9e-3, 12e-3])[:, None, None]
lev5 = SigmaLevels(5)
pi8 = np.full((1, 1), 9.0e4)
z8 = np.zeros_like(th8)
th_, _, _, inf8 = dry_convective_adjustment_pav(th8, z8, z8, pi8, lev5, scalars=[qv8])
q_ = inf8["scalars"][0]
w8 = np.asarray(lev5.dsigma)[:, None, None]
report("the PAV adjustment mixes qv over the adjusted block (column water conserved) and leaves the rest",
       abs((q_ * w8).sum() / (qv8 * w8).sum() - 1) < 1e-14 and np.ptp(q_[2:4]) < 1e-15 and np.allclose(q_[[0, 1, 4]], qv8[[0, 1, 4]], rtol=1e-14, atol=0),
       f"q (g/kg) {np.round(qv8.ravel() * 1e3, 2)} -> {np.round(q_.ravel() * 1e3, 3)}; "
       f"theta {np.round(th_.ravel(), 2)}")

# 9. Simple ice (P-81): off is the warm scheme; on conserves, slows snow, melts, saturates over ice
from nh2d import P0 as _P0, KAPPA as _KAP
_rng = np.random.default_rng(1)
_nz, _ny, _nx = 30, 5, 6
_p = np.linspace(2.5e4, 1.0e5, _nz)[:, None, None] * np.ones((1, _ny, _nx))
_T = np.linspace(225.0, 295.0, _nz)[:, None, None] + _rng.normal(0, 1, (_nz, _ny, _nx))
_th = _T * (_P0 / _p) ** _KAP
_rho = _p / (287.05 * _T)
_dz = np.full_like(_p, 300.0)
_qs = mo.EPS * mo.es_bolton(_T) / (_p - mo.es_bolton(_T))
_qv = _qs * _rng.uniform(0.8, 1.1, _p.shape); _qc = _rng.uniform(0, 2e-3, _p.shape); _qr = _rng.uniform(0, 3e-3, _p.shape)
warm = _T > mo.T_ICE + 5.0
_a = mo.kessler(_th, _qv, _qc, _qr, _p, _rho, _dz, 20.0)
_b = mo.kessler(_th, _qv, _qc, _qr, _p, _rho, _dz, 20.0, ice=True)
W0 = (_rho * _dz * (_qv + _qc + _qr)).sum(0); W1 = (_rho * _dz * (_b[1] + _b[2] + _b[3])).sum(0) + _b[4]
_Tw = np.full_like(_p, 290.0); _thw = _Tw * (_P0 / _p) ** _KAP
_aw = mo.kessler(_thw, _qv, _qc, _qr, _p, _p / (287.05 * _Tw), _dz, 20.0)
_bw = mo.kessler(_thw, _qv, _qc, _qr, _p, _p / (287.05 * _Tw), _dz, 20.0, ice=True)
report("simple ice: in an all-warm column ice=True equals ice=False exactly; with cold layers water + surface "
       "precipitation is conserved to 1e-14",
       all(np.array_equal(x, y) for x, y in zip(_aw, _bw)) and np.abs(W1 / W0 - 1).max() < 1e-14,
       f"water change {np.abs(W1 / W0 - 1).max():.1e}")
_qr2 = np.where(_T < 263.0, 2e-3, 0.0); _z0 = np.zeros_like(_qv)
_w = mo.kessler(_th, _z0, _z0, _qr2, _p, _rho, _dz, 20.0, warm_rain=False)
_i = mo.kessler(_th, _z0, _z0, _qr2, _p, _rho, _dz, 20.0, warm_rain=False, ice=True)
_below = _qr2 == 0
report("simple ice: snow from the cold layers moves down more slowly than rain would (one 20 s step)",
       (_i[3] * _below).sum() < 0.7 * (_w[3] * _below).sum(),
       f"moved below the cold layers: rain {(_w[3] * _below).sum():.4f}, snow {(_i[3] * _below).sum():.4f} (kg/kg summed)")
# melting: one cold layer of snow above warm layers; the warm layer that receives it cools by LF/cp per kg
_pm = np.array([6.0e4, 6.5e4, 7.0e4])[:, None, None]; _Tm = np.array([265.0, 275.0, 280.0])[:, None, None]
_thm = _Tm * (_P0 / _pm) ** _KAP; _rhom = _pm / (287.05 * _Tm); _dzm = np.full_like(_pm, 500.0)
_qrm = np.array([2e-3, 0.0, 0.0])[:, None, None]; _zm = np.zeros_like(_qrm)
_m = mo.kessler(_thm, _zm, _zm, _qrm, _pm, _rhom, _dzm, 60.0, warm_rain=False, ice=True)
_dT = (_m[0][1] - _thm[1]) * (_pm[1] / _P0) ** _KAP
report("simple ice: snow falling into the layer above freezing cools it (melting), the layers below are untouched by melting",
       float(_dT[0, 0]) < 0 and float(_dT[0, 0]) > -1.0,
       f"temperature change of the first warm layer {float(_dT[0, 0]):+.4f} K")
_ad = mo.kessler(_th, _qv * 1.2, _qc, _z0, _p, _rho, _dz, 20.0, warm_rain=False, ice=True)
_Tn = (_p / _P0) ** _KAP * _ad[0]
_cold = (_T < mo.T_ICE - 5.0) & (_ad[2] > 1e-6)
_rhi = _ad[1] / (mo.EPS * mo.es_ice(_Tn) / (_p - mo.es_ice(_Tn)))
# one Newton step, as in the warm scheme: from 20 % supersaturation it lands within ~1 %
report("simple ice: 20 %-supersaturated cold cloud is adjusted to saturation over ice (one Newton step, within 2 %)",
       _cold.any() and np.abs(_rhi[_cold] - 1).max() < 2e-2,
       f"RH over ice after the step: {float(_rhi[_cold].min()):.4f}-{float(_rhi[_cold].max()):.4f}")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

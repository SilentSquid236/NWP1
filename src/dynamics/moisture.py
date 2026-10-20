"""
Warm-rain moisture for the non-hydrostatic core (CAM stage S7).

WHAT IS CARRIED
Three mixing ratios (kg/kg): water vapour qv, cloud water qc and rain qr.
The core carries them mass-weighted, Q = mu q, stacked in one array
(3, nz, ny, nx) in the order of NAMES, and advects them in flux form with the
same third-order scheme as Theta (nh3d.NH3D._adv_q).

HOW MOISTURE ENTERS THE DYNAMICS: DENSITY POTENTIAL TEMPERATURE
With moisture the core's Theta is mu * theta_rho, where

    theta_rho = theta (1 + qv / EPS) / (1 + qv + qc + qr)

is the density potential temperature (Emanuel 1994, sec. 4.3). The gas law
for the moist air is p alpha = R_d T_rho, with alpha the inverse density of
the whole air (dry air, vapour and condensate) and T_rho = theta_rho Pi. The
core's equation of state is p = P0 (R_d Theta / (mu P0 alpha))^gamma, so with
theta_rho in place of theta it is exact for moist air, and both the buoyancy
of vapour and the loading of cloud and rain enter through it. Nothing in the
dry solver (NumPy or C) changes. mu stays the mass of the whole air column:
the analysis pressure already includes the vapour.

theta itself is never stored: theta = theta_rho / rho_factor(qv, qc, qr).
theta_rho is a function of quantities that are each conserved following the
air, so advecting it is the same as advecting theta and q separately, up to
the scheme's own mixing (second order in the differences).

APPROXIMATIONS, KNOWN
- mu does not change when rain falls out or the surface evaporates (at most
  ~0.3 hPa of column mass in a day of heavy rain).
- Constant cp = CP of dry air and Lv = LV; no ice.
- The mixing ratios are relative to dry air in the microphysics formulas,
  but the mass weight mu is the whole air column (a ~1 % difference).

WARM RAIN: Kessler (1969) as in Klemp and Wilhelmson (1978)
Transcribed from WRF's module_mp_kessler.F (the reference warm-rain scheme):
- rain falls first, with the terminal velocity 36.34 (rho qr 1e-3)^0.1364
  sqrt(rho_0 / rho), in substeps that keep the fall Courant number <= 0.75;
- autoconversion 1e-3 s^-1 (qc - 1e-3) and accretion 2.2 qc qr^0.875
  (implicit in qr);
- saturation adjustment: one Newton step on qv - qs(T, p), with Bolton's
  (1980) saturation vapour pressure es = 611.2 exp(17.67 (T - 273.15) /
  (T - 29.65)) Pa;
- rain evaporation in sub-saturated air (Klemp and Wilhelmson 1978, eq. 2.14).
Here index 0 is the lid and nz-1 the ground (WRF: 1 is the ground).
"""
import numpy as np

NAMES = ("qv", "qc", "qr")
EPS = 0.622          # R_d / R_v
LV = 2.5e6           # J/kg, latent heat of vaporisation
SVP1, SVP2, SVP3, SVPT0 = 611.2, 17.67, 29.65, 273.15     # Bolton (1980), Pa and K
KC1, KC2, KC3, KC4 = 1.0e-3, 1.0e-3, 2.2, 0.875           # Kessler autoconversion / accretion
MAX_CR_SED = 0.75    # largest fall Courant number per sedimentation substep
# Simple ice (P-81, 2026-10-13; after Dudhia 1989 / WRF WSM3): the same three
# species, read by temperature. Below T_ICE, qc is cloud ice and qr is snow:
# saturation over ice with the latent heat of sublimation, snow falling at
# VS0 sqrt(rho_0 / rho), ice-to-snow aggregation and accretion with Lin et al.
# (1983)'s efficiency exp(ICE_EXP (T - T_ICE)), and snow melting (cooling by LF)
# as it falls into a layer above freezing. Off unless kessler(..., ice=True).
LF = 3.34e5                                   # J/kg, latent heat of fusion
SVPI2, SVPI3, SVPIT0 = 21.8745584, 7.66, 273.16   # saturation over ice (Murray 1967), with SVP1
T_ICE = 273.15
VS0 = 1.0            # m/s, snow fall speed at the ground density
KI1, KI2 = 1.0e-3, 1.0e-4                     # ice -> snow aggregation rate (1/s) and threshold
ICE_EXP = 0.025      # 1/K, Lin et al. (1983) aggregation efficiency


def es_bolton(T):
    """Saturation vapour pressure over water (Pa), Bolton (1980)."""
    return SVP1 * np.exp(SVP2 * (T - SVPT0) / (T - SVP3))


def qsat(T, p):
    """Saturation mixing ratio (kg/kg) at temperature T (K) and pressure p (Pa)."""
    es = es_bolton(T)
    return EPS * es / np.maximum(p - es, 1.0)


def rho_factor(qv, qc=0.0, qr=0.0):
    """theta_rho / theta."""
    return (1.0 + qv / EPS) / (1.0 + qv + qc + qr)


def qv_from_rh(rh_percent, T, p):
    """Vapour mixing ratio from relative humidity (%) as e / es."""
    rh = np.clip(np.asarray(rh_percent, dtype=float), 0.0, 100.0) / 100.0
    e = rh * es_bolton(T)
    return EPS * e / np.maximum(p - e, 1.0)


def dewpoint(qv, p):
    """Dewpoint (K) from the vapour mixing ratio and pressure (Bolton inverse)."""
    e = np.maximum(qv * p / (EPS + qv), 1e-3)
    lnr = np.log(e / SVP1)
    return (SVP3 * lnr - SVP2 * SVPT0) / (lnr - SVP2)


def fill_negative(q, w):
    """
    Remove negative mixing ratios column by column, conserving the column
    mass: negatives are set to zero and the positive values are scaled down
    by the mass they owed. q: (nz, ...) mixing ratio; w: layer weights
    broadcastable to q (mass per level, e.g. ds[:, None, None]). Returns the
    filled q and the mass that could not be repaid (a column whose total was
    negative), as a column sum of w q.
    """
    neg = np.minimum(q, 0.0)
    pos = np.maximum(q, 0.0)
    owe = -(neg * w).sum(axis=0)
    have = (pos * w).sum(axis=0)
    with np.errstate(divide="ignore", invalid="ignore"):
        scale = np.where(have > owe, 1.0 - owe / np.where(have > 0, have, 1.0), 0.0)
    lost = np.where(have > owe, 0.0, owe - have)
    return pos * scale[None], lost


def _vt(qr, rho, rho_bottom):
    """Rain terminal velocity (m/s), Klemp and Wilhelmson (1978)."""
    return 36.34 * np.maximum(0.001 * rho * qr, 0.0) ** 0.1364 * np.sqrt(rho_bottom / rho)


def es_ice(T):
    """Saturation vapour pressure over ice (Pa), Murray (1967)."""
    return SVP1 * np.exp(SVPI2 * (T - SVPIT0) / (T - SVPI3))


def kessler(theta, qv, qc, qr, p, rho, dz, dt, warm_rain=True, ice=False):
    """
    One Kessler step on columns (index 0 = lid). theta: potential temperature
    (K); p (Pa), rho (kg/m^3, of the whole air) and dz (m, layer depth) as
    from the model. Returns (theta, qv, qc, qr, rain), rain in mm (kg/m^2)
    reaching the ground in this step. Inputs are not modified.
    warm_rain=False turns autoconversion and accretion off (saturation
    adjustment only: the reversible case of Bryan and Fritsch 2002).
    ice=True: the simple-ice reading of qc and qr below T_ICE (module notes).
    """
    from nh2d import CP, P0, KAPPA
    theta, qv, qc, qr = (np.array(a, dtype=float, copy=True) for a in (theta, qv, qc, qr))
    nz = qr.shape[0]
    rho_b = rho[-1]
    pii = (p / P0) ** KAPPA
    # phase of each layer, from the temperature at the start of the step
    cold = (pii * theta < T_ICE) if ice else np.zeros(qr.shape, dtype=bool)
    vsnow = VS0 * np.sqrt(rho_b[None] / rho)
    def fallspeed(q):
        v = _vt(q, rho, rho_b[None])
        return np.where(cold, vsnow, v) if ice else v
    # --- sedimentation (first, as in WRF) -----------------------------------
    vt = fallspeed(qr)
    crmax = (vt * dt / dz).max(axis=0)
    nfall = np.maximum(1, np.floor(0.5 + crmax / MAX_CR_SED).astype(int))
    dtf = dt / nfall
    rain = np.zeros(qr.shape[1:])
    for n in range(int(nfall.max())):
        act = n < nfall
        flux = rho * qr * vt                                   # kg/m^2/s out of each layer's bottom
        div = flux.copy()
        div[1:] -= flux[:-1]
        qr = np.where(act[None], qr - dtf[None] * div / (rho * dz), qr)
        if ice:
            # snow from a frozen layer entering a layer above freezing melts there
            melt = np.zeros_like(qr)
            melt[1:] = np.where(cold[:-1] & ~cold[1:], flux[:-1], 0.0) * dtf[None] / (rho[1:] * dz[1:])
            theta = np.where(act[None], theta - LF / (CP * pii) * melt, theta)
        rain = rain + np.where(act, flux[-1] * dtf, 0.0)
        vt = fallspeed(qr)
    # --- autoconversion and accretion ----------------------------------------
    if warm_rain:
        factorn = 1.0 / (1.0 + KC3 * dt * np.maximum(qr, 0.0) ** KC4)
        qrprod = qc * (1.0 - factorn) + factorn * KC1 * dt * np.maximum(qc - KC2, 0.0)
        if ice:
            E = np.exp(ICE_EXP * (pii * theta - T_ICE))
            fi = 1.0 / (1.0 + KC3 * E * dt * np.maximum(qr, 0.0) ** KC4)
            qsprod = qc * (1.0 - fi) + fi * E * KI1 * dt * np.maximum(qc - KI2, 0.0)
            qrprod = np.where(cold, qsprod, qrprod)
        qc = np.maximum(qc - qrprod, 0.0)
        qr = np.maximum(qr + qrprod, 0.0)
    # --- saturation adjustment and rain evaporation ---------------------------
    temp = pii * theta
    gam = LV / (CP * pii)
    es = es_bolton(temp)
    f5 = SVP2 * (SVPT0 - SVP3) * LV / CP
    den = (temp - SVP3) ** 2
    if ice:
        LS = LV + LF
        gam = np.where(cold, LS / (CP * pii), gam)
        es = np.where(cold, es_ice(temp), es)
        f5 = np.where(cold, SVPI2 * (SVPIT0 - SVPI3) * LS / CP, f5)
        den = np.where(cold, (temp - SVPI3) ** 2, den)
    qvs = EPS * es / (p - es)
    prod = (qv - qvs) / (1.0 + p / (p - es) * qvs * f5 / den)
    rcgs = 0.001 * rho
    rq = rcgs * qr
    ern = np.minimum.reduce([
        dt * (((1.6 + 124.9 * rq ** 0.2046) * rq ** 0.525) / (2.55e8 / (p * qvs) + 5.4e5))
        * (np.maximum(qvs - qv, 0.0) / (rcgs * qvs)),
        np.maximum(-prod - qc, 0.0),
        qr])
    product = np.maximum(prod, -qc)
    theta = theta + gam * (product - ern)
    qv = np.maximum(qv - product + ern, 0.0)
    qc = qc + product
    qr = qr - ern
    return theta, qv, qc, qr, rain

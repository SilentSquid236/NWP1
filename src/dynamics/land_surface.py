"""
A ground temperature with its own energy budget, coupled to the lowest layer (P-59 step 2).

WHY

Three open problems come from the same missing physics: the ground has no
temperature of its own.
- P-59: there is no diurnal cycle in 2 m temperature.
- P-69: at night the model's 10 m wind is 1.6–1.8 times the observed. By day
  it is about right.
- The night warm bias left by the prescribed flux (`diurnal.py`).

Real nights are calm near the ground because the ground cools by longwave
radiation, a stable layer forms, and the surface layer decouples from the
flow above. NWP1's drag is always neutral (`theta_surface=None`) and its air
is never cooled from below, so neither can happen.

The prescribed flux was a stand-in. The user declined it as the default
(2026-10-02) because the aim is a convection-allowing model, and every CAM
gets its surface fluxes from a surface energy budget. This module is that
budget, in its simplest established form. Stage S5 of `docs/CAM_DESIGN.md`
will build on it.

WHAT THIS IS

The force-restore ground temperature of Deardorff (1978), over land:

    dTg/dt = C_G (Rn - H - LE)  -  (2 pi / tau_d) (Tg - T2)
    C_G    = 2 sqrt(pi) / (rho_c d1),   d1 = sqrt(k_s tau_d)

- Rn = (1 - albedo) SW_down + eps (LW_down - sigma Tg^4) is the net radiation.
  - SW_down = S0 tau max(sin e, 0) is clear-sky sunshine, as in `diurnal.py`.
  - LW_down = eps_a sigma Ta^4, from the lowest-level air temperature.
- H = rho cp C_h U (theta_g - theta_1) (p_s / P0)^kappa is the sensible heat
  flux.
  - C_h is the neutral transfer coefficient times the Louis (1979) heat
    stability function.
  - U = max(|V_1|, U_MIN) keeps free convection going in calm air.
- LE = EF max(Rn, 0) is latent heat. The model is dry, so evaporation is a
  fixed share of daytime net radiation. It is the one constant that stands in
  for moisture, and it goes when moisture arrives.
- T2 is the deep-soil temperature, held at its initial value for the 24 h run.

Over water the surface temperature is held at its initial value. Lakes and
the ocean barely change in a day.

The initial Tg, and both T2 and the water temperature, come from the
lowest-level air temperature at the cycle time, brought down to the ground
with the standard lapse rate. Nothing from after the cycle time is used.

The flux heats or cools the lowest layer exactly as in `diurnal.py`:
dtheta_1/dt = g H / (cp dp_1) (P0 / p_1)^kappa.

THE STABLE SIDE OF THE HEAT FLUX. Louis's heat function falls off much faster
than his momentum function on the stable side (F_h = 0.003 at Ri_b = 5). In
the night column test (`test_land_surface.py`), with the Louis heat function,
the heat flux collapsed to -1.1 W/m^2 at Ri_b 2.0 after 6 h, with the ground
5.5 K below the air: a runaway decoupling. Clear nights over land have fluxes
of tens of W/m^2. So by default (`long_tail=True`) the heat flux uses the
momentum function's long tail on the stable side as well, and Louis's own
form on the unstable side. This is a choice made from that measurement. It
is not taken from a reference.

The drag uses theta_g as its surface temperature, with the Louis (1979)
momentum function (`louis_momentum`). That function has a long stable tail,
so drag weakens at night but never switches off. The older cut-off form in
`surface.stability_function` drops drag to zero at Ri_b = 0.2, which with a
cooling ground decouples the surface completely. The ground would then cool
without limit.

CONSTANTS (values from the literature, not fitted to any case here)

| name | value | source / reasoning |
|---|---|---|
| albedo | 0.18 land | mixed forest and cropland |
| eps (ground emissivity) | 0.95 | vegetated land |
| eps_a (clear-sky air emissivity) | 0.79 | Brutsaert (1975) at e = 12 hPa and T = 285 K, a late-September Northeast value; there are no clouds |
| rho_c (soil heat capacity) | 2.0e6 J/m^3/K | moist loam |
| k_s (soil thermal diffusivity) | 5.0e-7 m^2/s | moist loam |
| tau_d | 86400 s | the forcing period |
| EF (evaporative fraction) | 0.5 | a Bowen ratio of 1, typical of an early-autumn mid-latitude land surface |
| U_MIN | 1.0 m/s | a floor on the wind speed for the heat flux |

LIMITS, STATED PLAINLY
- No clouds: every day and every night is clear. Overcast nights will cool
  too much and overcast days warm too much.
- One soil type, one albedo and one roughness everywhere on land.
- The Great Lakes stand above 0 m and are classed as land by the
  terrain-above-sea-level mask, as in `diurnal.py`.
- EF is a constant, not a moisture budget.
"""
import math
from datetime import timedelta

import numpy as np

from backend import xp_of
from sigma import CP, P0, G0, KAPPA, RD
from diurnal import S0, TAU, sin_solar_elevation

SIGMA_SB = 5.670374e-8          # Stefan-Boltzmann, W/m^2/K^4
ALBEDO = 0.18
EPS_G = 0.95
EPS_A = 0.79
RHO_C = 2.0e6
K_SOIL = 5.0e-7
TAU_D = 86400.0
EF = 0.5
U_MIN = 1.0
LAPSE = 0.0065

LOUIS_B, LOUIS_D = 5.0, 5.0
LOUIS_CM, LOUIS_CH = 7.4, 5.3


def brutsaert_emissivity(e_hpa, T):
    """Clear-sky air emissivity, Brutsaert (1975): 1.24 (e/T)^(1/7), e in hPa."""
    return 1.24 * (e_hpa / T) ** (1.0 / 7.0)


def louis_momentum(Ri, a2, z_over_z0):
    """Louis (1979) momentum stability factor F_m (multiplies the neutral Cd)."""
    xp = xp_of(Ri)
    pos = xp.maximum(Ri, 0.0)
    neg = xp.minimum(Ri, 0.0)
    st = 1.0 / (1.0 + 2.0 * LOUIS_B * pos / xp.sqrt(1.0 + LOUIS_D * pos))
    un = 1.0 - 2.0 * LOUIS_B * neg / (1.0 + 3.0 * LOUIS_B * LOUIS_CM * a2
                                      * xp.sqrt(-neg * z_over_z0))
    return xp.where(Ri >= 0.0, st, un)


def louis_heat(Ri, a2, z_over_z0):
    """Louis (1979) heat stability factor F_h (multiplies the neutral C_h)."""
    xp = xp_of(Ri)
    pos = xp.maximum(Ri, 0.0)
    neg = xp.minimum(Ri, 0.0)
    st = 1.0 / (1.0 + 3.0 * LOUIS_B * pos * xp.sqrt(1.0 + LOUIS_D * pos))
    un = 1.0 - 3.0 * LOUIS_B * neg / (1.0 + 3.0 * LOUIS_B * LOUIS_CH * a2
                                      * xp.sqrt(-neg * z_over_z0))
    return xp.where(Ri >= 0.0, st, un)


def force_restore_coefficient(rho_c=RHO_C, k_s=K_SOIL, tau_d=TAU_D):
    """C_G = 2 sqrt(pi) / (rho_c d1), d1 = sqrt(k_s tau_d); K per J/m^2."""
    return 2.0 * math.sqrt(math.pi) / (rho_c * math.sqrt(k_s * tau_d))


class ForceRestoreSurface:
    """
    Ground temperature over land, fixed surface temperature over water.

    Arrays are NumPy (ny, nx). The model passes NumPy views of its lowest
    level and converts the returned theta tendency for the torch backend.
    """

    def __init__(self, lat, lon, land, start, T_air1, z1, z0=0.1,
                 albedo=ALBEDO, ef=EF, eps_a=EPS_A, long_tail=True):
        self.lat = np.asarray(lat, dtype=float)
        self.lon = np.asarray(lon, dtype=float)
        self.land = np.asarray(land, dtype=bool)
        self.start = start
        self.z0 = np.asarray(z0, dtype=float)       # scalar or (ny, nx) map
        self.albedo, self.ef, self.eps_a = float(albedo), float(ef), float(eps_a)
        self.long_tail = bool(long_tail)
        T0 = np.asarray(T_air1, dtype=float) + LAPSE * np.asarray(z1, dtype=float)
        self.Tg = T0.copy()          # ground (land) or water surface temperature
        self.T2 = T0.copy()          # deep soil; also the fixed water temperature
        self.C_G = force_restore_coefficient()
        self.last = {}

    def sw_down(self, t_seconds):
        when = self.start + timedelta(seconds=float(t_seconds))
        return S0 * TAU * np.maximum(sin_solar_elevation(self.lat, self.lon, when), 0.0)

    def step(self, t_seconds, dt, u1, v1, theta1, ps, p1, dp1, z1):
        """
        Advance Tg by dt; return d(theta_1)/dt (K/s) and theta_g for the drag.

        u1, v1, theta1: lowest-level wind and potential temperature;
        ps: surface pressure (Pa); p1, dp1: lowest full-level pressure and
        layer thickness (Pa); z1: lowest full-level height above ground (m).
        """
        ex_s = (ps / P0) ** KAPPA
        Ta = theta1 * (p1 / P0) ** KAPPA
        theta_g = self.Tg / ex_s
        U = np.maximum(np.hypot(u1, v1), U_MIN)
        a2 = (0.4 / np.log(np.maximum(z1, 2 * self.z0) / self.z0)) ** 2
        Ri = G0 * z1 * (theta1 - theta_g) / (theta1 * U ** 2)
        fh = louis_heat(Ri, a2, z1 / self.z0)
        if self.long_tail:
            fh = np.where(Ri > 0.0, louis_momentum(Ri, a2, z1 / self.z0), fh)
        ch = a2 * fh
        rho = ps / (RD * self.Tg)
        H = rho * CP * ch * U * (theta_g - theta1) * ex_s        # W/m^2, up
        sw = self.sw_down(t_seconds + 0.5 * dt)
        lw_dn = self.eps_a * SIGMA_SB * Ta ** 4
        Rn = (1.0 - self.albedo) * sw + EPS_G * (lw_dn - SIGMA_SB * self.Tg ** 4)
        LE = self.ef * np.maximum(Rn, 0.0)
        dTg = self.C_G * (Rn - H - LE) - (2.0 * math.pi / TAU_D) * (self.Tg - self.T2)
        self.Tg = np.where(self.land, self.Tg + dt * dTg, self.T2)
        dtheta = G0 * H / (CP * dp1) * (P0 / p1) ** KAPPA
        self.last = {"H": H, "Rn": Rn, "LE": LE, "Ri": Ri}
        return dtheta, self.Tg / ex_s

    def __repr__(self):
        return (f"ForceRestoreSurface(start {self.start:%Y-%m-%d %HZ}, albedo {self.albedo:g}, "
                f"EF {self.ef:g}, eps_a {self.eps_a:g}, land {self.land.mean():.0%})")

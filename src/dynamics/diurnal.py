"""
Prescribed diurnal surface sensible heat flux from solar geometry (P-59).

WHY

The dry core had no diurnal cycle: no surface heat flux, no radiation, no
sun. Verification shows the result as the largest error in the project, a
temperature bias that tracks the sun (about -8 K in the afternoon, +6 K
before dawn), and near-surface wind that decays in the interior because
nothing mixes momentum down by day.

WHAT THIS IS

The smallest step that puts the sun in, using no information from after
the cycle time and no one else's model output:

    H = f * S0 * tau * max(sin e, 0)  +  H_night        over land, W/m^2
    H = 0                                                over water

e is the solar elevation from the date, time and position. S0 = 1361 W/m^2
(Kopp and Lean 2011). tau = 0.75 is a clear-sky transmission, because the
dry model has no clouds. f = 0.2 is the share of the incoming sunshine that
goes into sensible heat, and H_night = -30 W/m^2 is a steady longwave
cooling of the ground. The sum is continuous: it crosses zero about an hour
after sunrise and peaks near 120 W/m^2 at a late-September noon in the
Northeast.

The flux heats (or cools) the lowest model layer:

    dT/dt = g H / (cp dp1),    dtheta/dt = dT/dt * (P0 / p1)^kappa

with dp1 the layer's pressure thickness. Nothing else is done here. By day
the heated layer becomes unstable and the convective adjustment mixes it,
and with it the momentum, upward, which is a crude convective boundary
layer. By night the cooled layer is stable and stays where it is.

LIMITS, STATED PLAINLY

- No clouds, so every day is clear. Overcast days will be too warm.
- The equation of time is ignored (solar noon is up to 16 min off).
- Land is "terrain above 0 m". The Great Lakes stand above sea level, so
  they are treated as land and heated. The ocean, at 0 m, gets no flux.
- f, tau and H_night are constants, not a surface energy budget. They are
  the knobs a later scheme replaces.
- The relaxation zone pulls the edges back to the frozen initial state, so
  the diurnal cycle there is damped.

Declination: Cooper (1969).
"""
import math
from datetime import datetime, timedelta

import numpy as np

from sigma import CP, P0, G0, KAPPA

S0 = 1361.0          # W/m^2, total solar irradiance (Kopp and Lean 2011)
TAU = 0.75           # clear-sky transmission
F_SENSIBLE = 0.2     # share of incoming sunshine into sensible heat
H_NIGHT = -30.0      # W/m^2, steady longwave cooling of the ground


def solar_declination(day_of_year):
    """Declination in radians (Cooper 1969)."""
    return math.radians(23.45) * math.sin(2.0 * math.pi * (284 + day_of_year) / 365.0)


def sin_solar_elevation(lat_deg, lon_deg, when):
    """sin(solar elevation) at `when` (a naive UTC datetime); arrays broadcast."""
    lat = np.radians(np.asarray(lat_deg, dtype=float))
    dec = solar_declination(when.timetuple().tm_yday)
    hours = when.hour + when.minute / 60.0 + (when.second + when.microsecond / 1e6) / 3600.0
    omega = np.radians(15.0 * (hours + np.asarray(lon_deg, dtype=float) / 15.0 - 12.0))
    return np.sin(lat) * math.sin(dec) + np.cos(lat) * math.cos(dec) * np.cos(omega)


def surface_heat_flux(sin_e, land, f_sensible=F_SENSIBLE, tau=TAU, h_night=H_NIGHT):
    """Sensible heat flux into the atmosphere, W/m^2 (positive upward)."""
    return np.where(land, f_sensible * S0 * tau * np.maximum(sin_e, 0.0) + h_night, 0.0)


class DiurnalHeating:
    """
    The lowest-layer theta tendency from the prescribed flux.

    lat, lon, land are (ny, nx) NumPy arrays; start is the cycle time as a
    naive UTC datetime. Works on NumPy; the model converts the 2-D result
    for the torch backend.
    """

    def __init__(self, lat, lon, land, start, f_sensible=F_SENSIBLE, tau=TAU,
                 h_night=H_NIGHT):
        self.lat = np.asarray(lat, dtype=float)
        self.lon = np.asarray(lon, dtype=float)
        self.land = np.asarray(land, dtype=bool)
        self.start = start
        self.f_sensible, self.tau, self.h_night = float(f_sensible), float(tau), float(h_night)

    def flux(self, t_seconds):
        when = self.start + timedelta(seconds=float(t_seconds))
        s = sin_solar_elevation(self.lat, self.lon, when)
        return surface_heat_flux(s, self.land, self.f_sensible, self.tau, self.h_night)

    def theta_tendency(self, t_seconds, pi, lev):
        """d(theta)/dt of the lowest level, K/s, shape (ny, nx)."""
        pi = np.asarray(pi, dtype=float)
        dp1 = float(np.asarray(lev.dsigma)[-1]) * pi
        p1 = lev.p_top + float(np.asarray(lev.sigma)[-1]) * pi
        dT = G0 * self.flux(t_seconds) / (CP * dp1)
        return dT * (P0 / p1) ** KAPPA

    def __repr__(self):
        return (f"DiurnalHeating(start {self.start:%Y-%m-%d %HZ}, f {self.f_sensible:g}, "
                f"tau {self.tau:g}, night {self.h_night:g} W/m2, land {self.land.mean():.0%})")

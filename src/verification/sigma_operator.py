"""
Observation operator for a sigma-coordinate forecast.

WHY A SEPARATE OPERATOR

`GridInterpolator` was written against pressure levels, where every column
shares the same vertical coordinate. In sigma coordinates the pressure of a
level depends on the column's surface pressure, so a single 1D array of level
pressures does not exist and `vertical()` has nothing to interpolate against.

It also carries a trap for surface observations. `at_observation(..., pressure
=None)` falls back to `field3d[0]`, and in this project index 0 is the MODEL
LID. A 2 m thermometer would have been scored against the 200 hPa field --
a plausible-looking number, off by seventy kelvin, that no test of the
interpolation itself would catch. That fallback is why this class overrides
the method rather than adding to it.

WHAT A SURFACE OBSERVATION ACTUALLY COMPARES AGAINST

An ASOS thermometer sits 2 m above a station whose elevation is known. The
model's lowest level sits some hundreds of metres above a grid-cell mean
terrain that is not the station's elevation -- on a 12 km grid a valley
station can be several hundred metres below its own cell. Three corrections,
in order, each of which is a real physical step and not a fudge:

  1. interpolate horizontally to the station's lat/lon
  2. bring the value from the model's lowest-level HEIGHT down to the station
     elevation along a lapse rate
  3. report it as the 2 m value

Step 2 uses the standard 6.5 K/km lapse rate, which is wrong exactly when it
matters most -- a clear calm night with a surface inversion, when the valley
station is coldest and the correction is largest. That limitation is recorded
here rather than hidden, because it will show up in the archive as a
persistent cold-season bias at low-elevation stations and should be read as
the operator's error, not the model's.
"""

import numpy as np

from obs_operator import GridInterpolator, elevation_correct_temperature

RD = 287.05
G0 = 9.80665
P0 = 100_000.0
KAPPA = RD / 1004.6


class SigmaInterpolator(GridInterpolator):
    """
    Interpolate a sigma-coordinate model state to observation locations.

    Constructed from the forecast's own vertical description, so an archive
    written from one run cannot be silently scored with another run's levels.
    """

    def __init__(self, grid, domain, sigma, p_top, pi, terrain):
        super().__init__(grid, domain, levels=None)
        self.sigma = np.asarray(sigma, dtype=float)
        self.p_top = float(p_top)
        self.pi = np.asarray(pi, dtype=float)
        self.terrain = np.asarray(terrain, dtype=float)
        self.nz = self.sigma.size

    # --- geometry ----------------------------------------------------------

    def pressure_column(self, lat, lon):
        """Pressure of every model level at a point, bottom value included."""
        pi = self.horizontal(self.pi, lat, lon)
        if pi is None:
            return None
        return self.p_top + self.sigma * pi

    def surface_pressure(self, lat, lon):
        pi = self.horizontal(self.pi, lat, lon)
        return None if pi is None else self.p_top + pi

    def model_elevation(self, lat, lon):
        """Grid-cell terrain height at a point."""
        return self.horizontal(self.terrain, lat, lon)

    def lowest_level_height(self, theta3d, lat, lon):
        """
        Height above sea level of the lowest model level.

        Hydrostatic, from the model's own layer thickness rather than a
        constant, so it follows the vertical grid instead of assuming one.
        """
        p = self.pressure_column(lat, lon)
        ps = self.surface_pressure(lat, lon)
        h = self.model_elevation(lat, lon)
        if p is None or ps is None or h is None:
            return None
        th = self.horizontal(theta3d[-1], lat, lon)
        if th is None:
            return None
        T = th * (p[-1] / P0) ** KAPPA
        return h + RD * T / G0 * np.log(ps / p[-1])

    # --- interpolation -----------------------------------------------------

    def at_pressure(self, field3d, lat, lon, pressure):
        """
        Interpolate a 3D field to a pressure, in log(p), using THIS column's
        pressures. Returns None outside the column's own range rather than
        extrapolating.
        """
        p = self.pressure_column(lat, lon)
        if p is None:
            return None
        if pressure > p[-1] or pressure < p[0]:
            return None

        col = []
        for k in range(field3d.shape[0]):
            v = self.horizontal(field3d[k], lat, lon)
            if v is None:
                return None
            col.append(v)

        # Index 0 is the lid, so p increases with index and log(p) does too.
        return float(np.interp(np.log(pressure), np.log(p), col))

    def at_surface(self, field3d, lat, lon):
        """The lowest MODEL level -- index -1, not index 0."""
        return self.horizontal(field3d[-1], lat, lon)

    def at_observation(self, field3d, lat, lon, pressure):
        """
        Override. `pressure is None` means a surface observation, which is the
        BOTTOM of the column here; the base class would have returned the lid.
        """
        if pressure is None:
            return self.at_surface(field3d, lat, lon)
        return self.at_pressure(field3d, lat, lon, pressure)

    # --- derived quantities ------------------------------------------------

    def temperature(self, theta3d, lat, lon, pressure=None):
        """Temperature, converting from the model's potential temperature."""
        p = self.pressure_column(lat, lon)
        if p is None:
            return None
        if pressure is None:
            th = self.at_surface(theta3d, lat, lon)
            return None if th is None else th * (p[-1] / P0) ** KAPPA
        th = self.at_pressure(theta3d, lat, lon, pressure)
        return None if th is None else th * (pressure / P0) ** KAPPA

    def station_temperature(self, theta3d, lat, lon, station_elev):
        """
        Model temperature brought to a station's elevation.

        Returns (value, info) so the size of the correction is recorded
        alongside the number. A match that needed 400 m of lapse correction is
        not the same evidence as one that needed 5 m, and the archive should be
        able to tell them apart later.
        """
        T = self.temperature(theta3d, lat, lon)
        z1 = self.lowest_level_height(theta3d, lat, lon)
        if T is None or z1 is None:
            return None, {}
        if station_elev is None:
            return T, {"elev_correction_m": 0.0, "model_level_height_m": z1}
        Tc = elevation_correct_temperature(T, z1, station_elev)
        return Tc, {"elev_correction_m": float(station_elev - z1),
                    "model_level_height_m": float(z1),
                    "elev_correction_K": float(Tc - T)}

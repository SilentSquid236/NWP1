"""
The observation operator H: map a gridded model state to observation space.

This is the piece that makes model and observation comparable at all, and it
is used identically for verification (score the difference) and assimilation
(correct using the difference). Getting it right once serves both.

Two interpolations, each with a specific requirement:

HORIZONTAL -- bilinear. Exact for linear fields, which is the property the
tests exploit: interpolating a plane must reproduce the plane exactly.

VERTICAL -- linear in log(pressure), NOT in pressure or height. Atmospheric
variables are close to linear in log(p); interpolating linearly in p puts
meaningful errors into the mid-troposphere where the levels are far apart.
This is a small detail that quietly corrupts results if missed.

ELEVATION -- surface observations need a lapse-rate correction between the
station height and the model's grid-cell height. A valley station and its
grid cell can differ by hundreds of metres, which is several degrees. Left
uncorrected, that difference is read as forecast error (or worse, injected as
an analysis increment).
"""

import numpy as np

STANDARD_LAPSE = 6.5e-3        # K/m, temperature decrease with height


class GridInterpolator:
    """
    Interpolates model fields to arbitrary lat/lon (and pressure) points.

    The model grid is Cartesian in metres; observations arrive in lat/lon, so
    we need a mapping. An equirectangular projection about the domain centre
    is accurate to well under a grid cell over a domain this size, and keeps
    the transform invertible and cheap.
    """

    def __init__(self, grid, domain, levels=None):
        self.grid = grid
        self.domain = domain
        self.levels = levels

        self.lat0 = 0.5 * (domain["lat_min"] + domain["lat_max"])
        self.lon0 = 0.5 * (domain["lon_min"] + domain["lon_max"])
        self.m_per_deg_lat = 111_132.0
        self.m_per_deg_lon = 111_320.0 * np.cos(np.radians(self.lat0))

    # --- coordinate mapping -------------------------------------------------

    def lonlat_to_xy(self, lat, lon):
        x = (lon - self.domain["lon_min"]) * self.m_per_deg_lon
        y = (lat - self.domain["lat_min"]) * self.m_per_deg_lat
        return x, y

    def in_domain(self, lat, lon, margin_cells=1):
        x, y = self.lonlat_to_xy(lat, lon)
        gx, gy = self.grid.dx * margin_cells, self.grid.dy * margin_cells
        return (gx <= x <= self.grid.Lx - gx) and (gy <= y <= self.grid.Ly - gy)

    # --- interpolation ------------------------------------------------------

    def horizontal(self, field2d, lat, lon):
        """
        Bilinear interpolation of a 2D field. Returns None outside the domain
        rather than extrapolating -- an extrapolated "observation match" is a
        fabricated number, and silently scoring against it is worse than
        having no match.
        """
        if not self.in_domain(lat, lon):
            return None

        x, y = self.lonlat_to_xy(lat, lon)
        fx, fy = x / self.grid.dx - 0.5, y / self.grid.dy - 0.5   # cell centres
        i0, j0 = int(np.floor(fx)), int(np.floor(fy))
        tx, ty = fx - i0, fy - j0

        ny, nx = field2d.shape
        if not (0 <= i0 < nx - 1 and 0 <= j0 < ny - 1):
            return None

        f00 = field2d[j0, i0]
        f10 = field2d[j0, i0 + 1]
        f01 = field2d[j0 + 1, i0]
        f11 = field2d[j0 + 1, i0 + 1]

        return float((1 - tx) * (1 - ty) * f00 + tx * (1 - ty) * f10 +
                     (1 - tx) * ty * f01 + tx * ty * f11)

    def vertical(self, column, pressure):
        """
        Interpolate a column (bottom -> top) to a pressure, linear in log(p).

        Returns None above the model top or below its base -- extrapolating a
        sounding beyond the model's own range invents data.
        """
        p = self.levels.p
        if pressure > p[0] or pressure < p[-1]:
            return None

        lp = np.log(p)
        target = np.log(pressure)
        # p decreases with index, so log(p) does too -- flip for np.interp.
        return float(np.interp(target, lp[::-1], np.asarray(column)[::-1]))

    def at_observation(self, field3d, lat, lon, pressure):
        """Interpolate a 3D field (nz, ny, nx) to a point."""
        if pressure is None:
            return self.horizontal(field3d[0], lat, lon)

        col = []
        for k in range(field3d.shape[0]):
            v = self.horizontal(field3d[k], lat, lon)
            if v is None:
                return None
            col.append(v)
        return self.vertical(col, pressure)


def elevation_correct_temperature(model_T, model_elev, station_elev,
                                  lapse=STANDARD_LAPSE):
    """
    Adjust a model temperature from the grid-cell height to the station height.

    Going DOWN from model height to a lower station means warming, hence the
    sign: dz = station - model, and T_station = T_model - lapse * dz.

    The standard lapse rate is an approximation that fails in strong
    inversions -- precisely the calm clear nights when valley stations are
    coldest and the correction is largest. Worth remembering when a station
    shows persistent cold bias that this correction does not remove.
    """
    dz = station_elev - model_elev
    return model_T - lapse * dz

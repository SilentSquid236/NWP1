"""
The analysis grid, and terrain that does not come from a weather model.

THE GRID IS THE FORECAST'S GRID

`forecast.build_grid()` derives dx and dy from the domain extent and the array
shape, assuming cells uniform in latitude and longitude with dx evaluated at
the central latitude. The analysis must put its values on exactly those cell
centres, or the model starts from a field that is shifted by up to half a
cell -- a small error that is also a spurious pressure gradient everywhere.

TERRAIN (P-53)

HRRR terrain is gone with the rest of HRRR. Terrain is static geography, not
weather, so any elevation model will do; ETOPO (NOAA, 1 arc-minute) is served
by NOAA's ERDDAP as plain CSV, which needs nothing but urllib. It is fetched
ONCE, block-averaged onto the grid, and cached under <data>/static/. Ocean
depths are clipped to sea level: the atmosphere's lower boundary over water
is the water surface.

Block-averaging is smoother than the point-sampling the HRRR path used at
stride 4, so peaks are lower than in every terrain measurement made before
2026-09-22. A capability number measured on one is not a number for the other.
"""

import io
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]
sys.path.insert(0, str(ROOT))

ERDDAP = "https://coastwatch.pfeg.noaa.gov/erddap/griddap/etopo180.csv"
M_PER_DEG_LAT = 111_132.0
M_PER_DEG_LON = 111_320.0


def grid_shape(domain, spacing_m=12_000.0):
    """(ny, nx) for a target spacing -- 12 km matches the old stride 4."""
    lat0 = 0.5 * (domain["lat_min"] + domain["lat_max"])
    ly = (domain["lat_max"] - domain["lat_min"]) * M_PER_DEG_LAT
    lx = ((domain["lon_max"] - domain["lon_min"]) * M_PER_DEG_LON
          * np.cos(np.radians(lat0)))
    return int(round(ly / spacing_m)), int(round(lx / spacing_m))


def cell_centres(domain, ny, nx):
    """2-D lat, lon of cell centres, row 0 at the southern edge."""
    dlat = (domain["lat_max"] - domain["lat_min"]) / ny
    dlon = (domain["lon_max"] - domain["lon_min"]) / nx
    lat = domain["lat_min"] + (np.arange(ny) + 0.5) * dlat
    lon = domain["lon_min"] + (np.arange(nx) + 0.5) * dlon
    lon2, lat2 = np.meshgrid(lon, lat)
    return lat2, lon2


def spacing_m(domain, ny, nx):
    """(dy, dx) exactly as forecast.build_grid computes them."""
    lat0 = 0.5 * (domain["lat_min"] + domain["lat_max"])
    dy = (domain["lat_max"] - domain["lat_min"]) * M_PER_DEG_LAT / ny
    dx = ((domain["lon_max"] - domain["lon_min"]) * M_PER_DEG_LON
          * np.cos(np.radians(lat0)) / nx)
    return dy, dx


def bilinear(field, lat, lon, domain, clip=False):
    """
    Sample a (ny, nx) cell-centred field at points; NaN outside the domain.

    clip=True returns the nearest edge value instead of NaN. That is only for
    evaluating a FIRST GUESS at an observation beyond the edge (a sounding in
    the ring), so the observation can still correct the edge; it is never
    used to score anything.
    """
    ny, nx = field.shape[-2:]
    fy = (np.asarray(lat) - domain["lat_min"]) / (
        domain["lat_max"] - domain["lat_min"]) * ny - 0.5
    fx = (np.asarray(lon) - domain["lon_min"]) / (
        domain["lon_max"] - domain["lon_min"]) * nx - 0.5
    inside = (fy >= -0.5) & (fy <= ny - 0.5) & (fx >= -0.5) & (fx <= nx - 0.5)
    fy = np.clip(fy, 0, ny - 1 - 1e-9)
    fx = np.clip(fx, 0, nx - 1 - 1e-9)
    j0, i0 = np.floor(fy).astype(int), np.floor(fx).astype(int)
    j1, i1 = np.minimum(j0 + 1, ny - 1), np.minimum(i0 + 1, nx - 1)
    wy, wx = fy - j0, fx - i0
    v = ((1 - wy) * (1 - wx) * field[..., j0, i0] + (1 - wy) * wx * field[..., j0, i1]
         + wy * (1 - wx) * field[..., j1, i0] + wy * wx * field[..., j1, i1])
    return v if clip else np.where(inside, v, np.nan)


# ---------------------------------------------------------------------------
# Terrain
# ---------------------------------------------------------------------------

def etopo_url(domain, stride=2):
    """ERDDAP griddap CSV for the domain. stride 2 = 2 arc-minutes (~3.7 km)."""
    q = (f"altitude%5B({domain['lat_min']}):{stride}:({domain['lat_max']})%5D"
         f"%5B({domain['lon_min']}):{stride}:({domain['lon_max']})%5D")
    return f"{ERDDAP}?{q}"


def parse_erddap_csv(text):
    """ERDDAP CSV (two header rows: names, units) -> lat, lon, value arrays."""
    data = np.loadtxt(io.StringIO(text), delimiter=",", skiprows=2)
    return data[:, 0], data[:, 1], data[:, 2]


def block_average(lat, lon, val, domain, ny, nx):
    """Mean of all source points falling in each grid cell."""
    j = np.floor((lat - domain["lat_min"]) / (domain["lat_max"] - domain["lat_min"]) * ny)
    i = np.floor((lon - domain["lon_min"]) / (domain["lon_max"] - domain["lon_min"]) * nx)
    ok = (j >= 0) & (j < ny) & (i >= 0) & (i < nx)
    j, i, v = j[ok].astype(int), i[ok].astype(int), val[ok]
    s = np.zeros((ny, nx)); n = np.zeros((ny, nx))
    np.add.at(s, (j, i), v)
    np.add.at(n, (j, i), 1)
    if (n == 0).any():
        raise ValueError(f"{int((n == 0).sum())} grid cells received no terrain "
                         f"samples; the source is coarser than the grid")
    return s / n


def limit_slope(terrain, domain, max_slope=0.0086):
    """
    Smooth terrain until its slope is at most `max_slope` (the model's own
    measure, sigma.terrain_slope). Returns (terrain, passes, slope).

    WHY 0.0086 (P-56). It is the steepest idealised terrain this core has
    been measured to survive 12/12 h over (2500 m ridge). Raw ETOPO at 12 km
    reaches 0.0536, and every real-terrain run died within 4 h; the same
    state smoothed to 0.0083 lived to 13.7 h. Operational models filter their
    orography for the same reason.
    """
    sys.path.insert(0, str(ROOT / "src" / "dynamics"))
    from grid import CGrid
    from sigma import smooth_terrain, terrain_slope
    ny, nx = terrain.shape
    dy, dx = spacing_m(domain, ny, nx)
    g = CGrid(nx, ny, dx, dy, edge_mode="replicate")
    if terrain_slope(terrain, g) <= max_slope:
        return terrain, 0, terrain_slope(terrain, g)
    out, n, s = smooth_terrain(terrain, g, target_slope=max_slope, max_passes=200)
    return np.clip(out, 0.0, None), n, s


def load_terrain(domain, ny, nx, cache_dir, fetcher=None, verbose=True):
    """
    Terrain (m) on the grid: from the cache, or fetched once from ERDDAP.

    The cache file name carries the shape, so a change of resolution fetches
    afresh instead of silently reusing the wrong grid.
    """
    cache_dir = Path(cache_dir)
    path = cache_dir / f"terrain_etopo_{ny}x{nx}.npz"
    if path.exists():
        z = np.load(path)
        return z["terrain"].astype(float), str(path)
    if fetcher is None:
        from netpolicy import PoliteFetcher
        fetcher = PoliteFetcher()
    url = etopo_url(domain)
    if verbose:
        print(f"  terrain        : fetching ETOPO once from NOAA ERDDAP", flush=True)
    lat, lon, alt = parse_erddap_csv(fetcher.get_text(url, timeout=300))
    terrain = np.clip(block_average(lat, lon, alt, domain, ny, nx), 0.0, None)
    cache_dir.mkdir(parents=True, exist_ok=True)
    np.savez_compressed(path, terrain=terrain.astype(np.float32),
                        source=np.array("ETOPO via NOAA ERDDAP etopo180"),
                        url=np.array(url))
    return terrain, str(path)

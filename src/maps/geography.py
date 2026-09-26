"""
Map geography: a Lambert conformal projection and the state and coast lines.

Nothing here needs cartopy, shapely or pyproj. The shared server allows no
installs (CLAUDE.md constraint 1), and matplotlib plus NumPy are enough to
draw a projected map.

PROJECTION. Lambert conformal conic on a sphere (formulas as in Snyder 1987), standard parallels 39 and
45 N, centred at 42 N 74 W: the projection US regional products use, with
scale error under 1 % across 37-47.5 N.

BOUNDARIES. Natural Earth (2026) 1:50m coastline, lakes, country borders and
state/province lines. That is static geography, like the ETOPO terrain, not
weather. They are fetched ONCE as GeoJSON (parsed with the json module),
clipped to the domain plus a margin, and cached as a small npz of
NaN-separated polylines under <data>/static/. A copy placed at
src/maps/boundaries_ne50m.npz is used first, so a machine with no network
can still draw lines. If neither exists and the fetch fails, maps are drawn
without lines, with a warning, instead of failing the run.
"""

import json
import math
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[2]

NE_BASE = ("https://raw.githubusercontent.com/nvkelso/natural-earth-vector/"
           "master/geojson/")
LAYERS = {
    "coast": "ne_50m_coastline.geojson",
    "lakes": "ne_50m_lakes.geojson",
    "countries": "ne_50m_admin_0_boundary_lines_land.geojson",
    "states": "ne_50m_admin_1_states_provinces_lines.geojson",
}
CACHE_NAME = "boundaries_ne50m.npz"
MARGIN_DEG = 3.0


# ---------------------------------------------------------------------------
# Projection
# ---------------------------------------------------------------------------

class LambertConformal:
    """Spherical Lambert conformal conic. x, y in metres; forward only."""

    R_EARTH = 6_371_000.0

    def __init__(self, lat1=39.0, lat2=45.0, lat0=42.0, lon0=-74.0):
        p1, p2, p0 = (math.radians(v) for v in (lat1, lat2, lat0))
        t = lambda p: math.tan(math.pi / 4 + p / 2)
        if abs(lat1 - lat2) < 1e-9:
            self.n = math.sin(p1)
        else:
            self.n = (math.log(math.cos(p1) / math.cos(p2))
                      / math.log(t(p2) / t(p1)))
        self.F = math.cos(p1) * t(p1) ** self.n / self.n
        self.rho0 = self.R_EARTH * self.F / t(p0) ** self.n
        self.lon0 = lon0
        self.params = dict(lat1=lat1, lat2=lat2, lat0=lat0, lon0=lon0)

    def __call__(self, lon, lat):
        lon = np.asarray(lon, dtype=float)
        lat = np.asarray(lat, dtype=float)
        rho = (self.R_EARTH * self.F
               / np.tan(np.pi / 4 + np.radians(lat) / 2) ** self.n)
        th = self.n * np.radians(lon - self.lon0)
        return rho * np.sin(th), self.rho0 - rho * np.cos(th)

    def scale_factor(self, lat):
        """Map scale factor k at latitude lat (1 on the standard parallels)."""
        p = np.radians(np.asarray(lat, dtype=float))
        rho = self.R_EARTH * self.F / np.tan(np.pi / 4 + p / 2) ** self.n
        return self.n * rho / (self.R_EARTH * np.cos(p))

    def __repr__(self):
        return "LambertConformal(" + ", ".join(
            f"{k}={v:g}" for k, v in self.params.items()) + ")"


# ---------------------------------------------------------------------------
# Boundaries
# ---------------------------------------------------------------------------

def _lines_of(geometry):
    """Every coordinate sequence in a GeoJSON geometry, as lists of (lon, lat)."""
    if geometry is None:
        return []
    kind, c = geometry.get("type"), geometry.get("coordinates")
    if kind == "LineString":
        return [c]
    if kind in ("MultiLineString", "Polygon"):
        return list(c)
    if kind == "MultiPolygon":
        return [ring for poly in c for ring in poly]
    if kind == "GeometryCollection":
        return [l for g in geometry.get("geometries", []) for l in _lines_of(g)]
    return []


def clip_polylines(lines, box):
    """
    Keep the parts of each polyline inside box = (lon0, lon1, lat0, lat1).

    Returns lon, lat arrays with NaN between pieces, ready for one ax.plot.
    A segment leaving the box ends its piece, so no line is drawn across the
    gap to where it re-enters.
    """
    lon0, lon1, lat0, lat1 = box
    out_x, out_y = [], []
    for line in lines:
        a = np.asarray(line, dtype=float)
        if a.ndim != 2 or len(a) < 2:
            continue
        a = a[:, :2]
        inside = ((a[:, 0] >= lon0) & (a[:, 0] <= lon1)
                  & (a[:, 1] >= lat0) & (a[:, 1] <= lat1))
        if not inside.any():
            continue
        # runs of consecutive inside points
        edges = np.flatnonzero(np.diff(np.r_[0, inside.astype(int), 0]))
        for s, e in zip(edges[::2], edges[1::2]):
            if e - s >= 2:
                out_x.extend(a[s:e, 0]); out_x.append(np.nan)
                out_y.extend(a[s:e, 1]); out_y.append(np.nan)
    return np.asarray(out_x, dtype=np.float32), np.asarray(out_y, dtype=np.float32)


def domain_box(domain, margin=MARGIN_DEG):
    return (domain["lon_min"] - margin, domain["lon_max"] + margin,
            domain["lat_min"] - margin, domain["lat_max"] + margin)


def build_boundaries(domain, fetch_text):
    """Fetch every layer and clip it. fetch_text(url) -> str."""
    box = domain_box(domain)
    layers = {}
    for name, fname in LAYERS.items():
        gj = json.loads(fetch_text(NE_BASE + fname))
        lines = [l for f in gj.get("features", [])
                 for l in _lines_of(f.get("geometry"))]
        layers[name] = clip_polylines(lines, box)
    return layers


def save_boundaries(path, layers):
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    arrays = {}
    for k, (x, y) in layers.items():
        arrays[f"{k}_lon"], arrays[f"{k}_lat"] = x, y
    np.savez_compressed(path, source=np.array("Natural Earth 1:50m via "
                                              "github.com/nvkelso/natural-earth-vector"),
                        **arrays)


def read_boundaries(path):
    z = np.load(path)
    names = sorted({k.rsplit("_", 1)[0] for k in z.files if k.endswith("_lon")})
    return {n: (z[f"{n}_lon"], z[f"{n}_lat"]) for n in names}


def load_boundaries(domain, cache_dir=None, fetch=True, verbose=True):
    """
    Boundary polylines {layer: (lon, lat)}: the bundled copy, the cache, or a
    one-time fetch. Returns {} (and warns) if none of those works.
    """
    bundled = Path(__file__).with_name(CACHE_NAME)
    if bundled.exists():
        return read_boundaries(bundled)
    if cache_dir is None:
        sys.path.insert(0, str(ROOT))
        import config
        cache_dir = Path(config.DATA_ROOT) / "static"
    cached = Path(cache_dir) / CACHE_NAME
    if cached.exists():
        return read_boundaries(cached)
    if not fetch:
        return {}
    try:
        sys.path.insert(0, str(ROOT))
        from netpolicy import PoliteFetcher
        fetcher = PoliteFetcher()
        if verbose:
            print("  boundaries     : fetching Natural Earth 1:50m once", flush=True)
        layers = build_boundaries(domain,
                                  lambda u: fetcher.get_text(u, timeout=300))
        save_boundaries(cached, layers)
        return layers
    except Exception as e:                    # maps still render without lines
        print(f"  WARNING: no state/coast lines ({type(e).__name__}: {e})",
              flush=True)
        return {}

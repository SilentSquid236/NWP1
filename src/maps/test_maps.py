"""
Tests for the forecast maps: projection, boundary clipping, diagnostics,
and (once rendering exists) that every product renders.

The diagnostics are checked against a standard atmosphere, whose heights and
sea-level pressure are known in closed form:

    T(p) = T0 (p / p0)^(R L / g),   z(p) = T0 / L (1 - (p / p0)^(R L / g))

Run:  python src/maps/test_maps.py
"""

import sys
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
ROOT = HERE.parents[1]
for p in (ROOT, ROOT / "src", ROOT / "src" / "dynamics", ROOT / "src" / "analysis",
          ROOT / "src" / "verification"):
    sys.path.insert(0, str(p))

from maps import geography, derive          # noqa: E402
from sigma import SigmaLevels, RD, G0       # noqa: E402

T0, PS0, L = 288.15, 101325.0, derive.LAPSE
EXP = RD * L / G0
results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def std_z(p):
    return T0 / L * (1.0 - (p / PS0) ** EXP)


def std_state(terrain):
    """theta, pi on 20 sigma levels for a standard atmosphere over terrain."""
    lev = SigmaLevels(20, p_top=20000.0)
    ps = PS0 * (1.0 - L * terrain / T0) ** (1.0 / EXP)
    pi = ps - lev.p_top
    p = lev.pressure(pi)
    T = T0 * (p / PS0) ** EXP
    theta = T * (derive.P0 / p) ** derive.KAPPA
    return lev, theta, pi


def test_projection_is_conformal():
    P = geography.LambertConformal()
    lat, lon, d = 44.3, -70.1, 1e-4
    x0, y0 = P(lon, lat)
    xe, ye = P(lon + d, lat)
    xn, yn = P(lon, lat + d)
    R = P.R_EARTH * np.radians(d)
    se = np.hypot(xe - x0, ye - y0) / (R * np.cos(np.radians(lat)))
    sn = np.hypot(xn - x0, yn - y0) / R
    k = P.scale_factor(np.array([37.0, 39.0, 42.0, 45.0, 47.5]))
    report("projection: conformal, scale 1 on 39/45 N, < 0.5 % error in the domain",
           abs(se - sn) < 1e-6 and abs(k[1] - 1) < 1e-12 and abs(k[3] - 1) < 1e-12
           and np.all(np.abs(k - 1) < 0.005),
           f"east/north scale {se:.7f}/{sn:.7f}; k over 37-47.5 N {k.min():.4f}-{k.max():.4f}")


def test_clip_breaks_lines_at_the_box():
    line = [(-90, 40), (-75, 40), (-74, 41), (-60, 41), (-73, 42), (-72, 42)]
    x, y = geography.clip_polylines([line], (-80, -70, 35, 45))
    pieces = np.split(x, np.flatnonzero(np.isnan(x)) + 1)
    pieces = [q[~np.isnan(q)] for q in pieces if (~np.isnan(q)).any()]
    report("boundaries: clipping keeps inside runs, breaks where a line leaves",
           len(pieces) == 2 and list(pieces[0]) == [-75, -74] and list(pieces[1]) == [-73, -72],
           f"pieces {[list(map(float, q)) for q in pieces]}")


def test_geojson_geometry_kinds():
    gj = [{"type": "LineString", "coordinates": [[0, 0], [1, 1]]},
          {"type": "MultiLineString", "coordinates": [[[0, 0], [1, 1]], [[2, 2], [3, 3]]]},
          {"type": "Polygon", "coordinates": [[[0, 0], [1, 0], [1, 1], [0, 0]]]},
          {"type": "MultiPolygon", "coordinates": [[[[0, 0], [1, 0], [0, 0]]], [[[5, 5], [6, 5], [5, 5]]]]}]
    n = [len(geography._lines_of(g)) for g in gj]
    report("boundaries: every GeoJSON geometry kind yields its lines", n == [1, 2, 1, 2], f"lines {n}")


def test_heights_match_standard_atmosphere():
    terrain = np.zeros((6, 7)); terrain[:, 4:] = 1000.0
    lev, theta, pi = std_state(terrain)
    p, T, Z, ps = derive.column_state(theta, pi, lev.sigma, lev.p_top, terrain)
    err = {}
    for Lh in (850, 700, 500, 250):
        z, _ = derive.to_pressure(Z, p, Lh * 100.0, "Z", ps=ps, T=T, Z=Z, terrain=terrain)
        err[Lh] = float(np.abs(z - std_z(Lh * 100.0)).max())
    z1000, below = derive.to_pressure(Z, p, 1e5, "Z", ps=ps, T=T, Z=Z, terrain=terrain)
    err[1000] = float(np.abs(z1000 - std_z(1e5)).max())
    report("heights: 1000-250 hPa within 5 m of the standard atmosphere, flat and 1000 m ground",
           max(err.values()) < 5.0 and below[:, 4:].all() and not below[:, :4].any(),
           "max error (m) " + ", ".join(f"{k}: {v:.2f}" for k, v in err.items())
           + "; 1000 hPa flagged below ground only over the hill")


def test_mslp_recovers_sea_level_pressure():
    terrain = np.linspace(0, 1500, 12).reshape(3, 4)
    lev, theta, pi = std_state(terrain)
    p, T, Z, ps = derive.column_state(theta, pi, lev.sigma, lev.p_top, terrain)
    m = derive.mslp(ps, T, Z, terrain)
    report("MSLP: a standard atmosphere over 0-1500 m reduces to 1013.25 hPa",
           float(np.abs(m - PS0).max()) < 30.0,
           f"MSLP {m.min()/100:.2f}-{m.max()/100:.2f} hPa (surface {ps.min()/100:.0f}-{ps.max()/100:.0f})")


def test_interpolation_exact_at_model_levels():
    terrain = np.zeros((4, 5))
    lev, theta, pi = std_state(terrain)
    p, T, Z, ps = derive.column_state(theta, pi, lev.sigma, lev.p_top, terrain)
    errs = [float(np.abs(derive.to_pressure(T, p, float(p[k, 0, 0]), "T", ps=ps)[0] - T[k]).max())
            for k in range(lev.nz)]
    report("interpolation: reproduces the model-level value at its own pressure",
           max(errs) < 1e-9, f"max error {max(errs):.1e} K over {lev.nz} levels")


def test_vorticity():
    ny, nx, dx, dy = 20, 24, 12e3, 12e3
    y, x = np.meshgrid(np.arange(ny) * dy, np.arange(nx) * dx, indexing="ij")
    u0 = derive.relative_vorticity(np.full((ny, nx), 10.0), np.full((ny, nx), -3.0), dx, dy)
    W = 1e-4
    z = derive.relative_vorticity(-W * y, W * x, dx, dy)
    report("vorticity: 0 for uniform flow, 2 Omega for solid-body rotation",
           np.abs(u0).max() < 1e-15 and np.abs(z - 2 * W).max() < 1e-12,
           f"uniform max|zeta| {np.abs(u0).max():.1e}; solid body {z.min():.3e}-{z.max():.3e} (2W = {2*W:.1e})")


def test_destagger():
    ny, nx = 5, 6
    u = np.tile(np.arange(nx, dtype=float), (ny, 1))          # u at west faces x = i
    v = np.tile(np.arange(ny, dtype=float)[:, None], (1, nx))
    uc, vc = derive.destagger(u, v)
    report("destagger: face winds averaged to cell centres (edges replicate)",
           np.allclose(uc[:, :-1], np.arange(nx - 1) + 0.5) and np.allclose(uc[:, -1], nx - 1)
           and np.allclose(vc[:-1, 0], np.arange(ny - 1) + 0.5),
           f"u row {uc[0].tolist()}")


def test_wind_rotation_follows_the_meridians():
    from maps import render
    lat = np.array([[40.0, 40.0], [46.0, 46.0]]); lon = np.array([[-81.0, -67.0], [-81.0, -67.0]])
    m = render.Mapper(lat, lon, {})
    worst = 0.0
    for j in range(2):
        for i in range(2):
            la, lo, d = lat[j, i], lon[j, i], 1e-4
            x0, y0 = m.P(lo, la)
            for (u, v), (dlo, dla) in (((1.0, 0.0), (d, 0.0)), ((0.0, 1.0), (0.0, d))):
                x1, y1 = m.P(lo + dlo, la + dla)
                want = np.arctan2(y1 - y0, x1 - x0)
                ux, vy = m.rotate(np.array(u), np.array(v), m.theta[j, i])
                worst = max(worst, abs(np.angle(np.exp(1j * (np.arctan2(vy, ux) - want)))))
    report("wind barbs: east/north winds rotated onto the projected parallels and meridians",
           worst < 1e-4, f"worst direction error {np.degrees(worst):.2e} deg at the domain corners")


def test_every_product_renders():
    """Every forecast product, from a standard-atmosphere forecast, to PNG."""
    import tempfile
    from maps import render
    from geo import cell_centres, spacing_m
    import config
    ny, nx = 20, 24
    lat, lon = cell_centres(config.DOMAIN, ny, nx)
    dy, dx = spacing_m(config.DOMAIN, ny, nx)
    terrain = np.zeros((ny, nx)); terrain[5:9, 5:12] = 600.0
    lev, theta, pi = std_state(terrain)
    u = np.full(theta.shape, 12.0); v = np.full(theta.shape, -4.0)
    f = {"theta": theta[None], "pi": pi[None], "u": u[None], "v": v[None],
         "sigma": lev.sigma, "p_top": lev.p_top, "terrain": terrain,
         "hours": np.array([1.0]), "lat": lat, "lon": lon, "dx": dx, "dy": dy}
    d = derive.snapshot(f, 0)
    d["dT1"] = np.full(terrain.shape, -0.8)
    m = render.Mapper(lat, lon, {"coast": (np.array([-75.0, -72.0, np.nan], np.float32),
                                           np.array([40.0, 41.0, np.nan], np.float32))})
    sub = render.time_labels(__import__("datetime").datetime(2026, 9, 21, 12), 1)
    sizes = {}
    with tempfile.TemporaryDirectory() as tmp:
        for key, fn in render.FORECAST_PRODUCTS.items():
            fig = fn(m, d, sub, ("left", "right"))
            path = Path(tmp) / f"{key}.png"
            fig.savefig(path, dpi=render.DPI)
            render.plt.close(fig)
            sizes[key] = path.stat().st_size
        pairs = [{"lat": 42.0, "lon": -74.0, "err": 1.5}, {"lat": 43.0, "lon": -72.0, "err": -2.0}]
        fig = render.product_error(m, {"terrain": terrain}, sub, ("", ""), pairs, "t")
        fig.savefig(Path(tmp) / "err.png", dpi=render.DPI); render.plt.close(fig)
        sizes["err_t"] = (Path(tmp) / "err.png").stat().st_size
    report("rendering: every forecast product and the error map write a PNG",
           len(sizes) == len(render.FORECAST_PRODUCTS) + 1 and min(sizes.values()) > 10_000,
           f"{len(sizes)} images, {min(sizes.values()) // 1000}-{max(sizes.values()) // 1000} kB")


def test_viewer_embeds_every_product():
    import tempfile, json as _json
    from maps import viewer
    prods = {"mslp": {"group": "Surface", "label": "MSLP", "hours": [0, 1, 2]},
             "err_t": {"group": "Verification", "label": "T error", "hours": [1, 2]}}
    with tempfile.TemporaryDirectory() as tmp:
        p = Path(tmp) / "index.html"
        viewer.write_viewer(p, "run", prods, "forecast complete", "note")
        html = p.read_text(encoding="utf-8")
    j = html[html.index("const M = ") + 10: html.index(";\nconst groups")]
    man = _json.loads(j)
    report("viewer: self-contained page embeds every product and hour, no external URLs",
           man["products"] == prods and man["all_hours"] == [0, 1, 2]
           and "http" not in html.replace("http-equiv", ""),
           f"{len(man['products'])} products, hours {man['all_hours']}, {len(html) // 1000} kB")


def test_forecast_hours_match_real_snapshot_times():
    import make_maps
    t = np.ceil(np.arange(1, 17) * 3600 / 17.1 - 1e-9) * 17.1 / 3600   # as run_forecast saves them
    got = make_maps.match_hours(t)
    report("hours: every forecast hour is found although snapshots land seconds late",
           sorted(got) == list(range(1, 17)) and all(abs(t[j] - h) < 0.01 for h, j in got.items()),
           f"snapshot at {t[0]:.4f} h, {t[1]:.4f} h ... -> hours {min(got)}-{max(got)} ({len(got)})")


def test_valid_time_labels():
    from datetime import datetime
    from maps import render
    a = render.time_labels(datetime(2026, 9, 23, 6), 12)
    b = render.time_labels(datetime(2026, 12, 1, 0), 18)
    report("titles: init, forecast hour, and valid time in UTC and US Eastern (EDT/EST)",
           a[0].startswith("Init: 06z Sep 23 2026") and "Forecast hour: 12" in a[0]
           and a[1].startswith("Valid: 18z Wed Sep 23 2026") and "2 PM EDT Wed" in a[1]
           and b[1].startswith("Valid: 18z Tue Dec 1 2026") and "1 PM EST Tue" in b[1],
           f"{a[0]} | {a[1]} ; {b[1]}")


def test_pixel_geometry_for_hover():
    """A marker drawn at a lat/lon lands on the pixel the viewer computes for it."""
    import tempfile
    from PIL import Image  # Pillow (Pillow Contributors 2026)
    from maps import render
    from geo import cell_centres
    import config
    lat, lon = cell_centres(config.DOMAIN, 97, 110)
    m = render.Mapper(lat, lon, {})
    fig, ax = m.frame("t", ("a", "b"))
    pts = [(-74.0, 42.0), (-79.5, 38.5), (-68.5, 46.5)]
    for lo, la in pts:
        x, y = m.P(lo, la); ax.plot([x], [y], "s", color=(1, 0, 0), ms=5, zorder=20)
    with tempfile.TemporaryDirectory() as tmp:
        fig.savefig(Path(tmp) / "g.png", dpi=render.DPI); render.plt.close(fig)
        im = np.asarray(Image.open(Path(tmp) / "g.png").convert("RGB")).astype(int)
    F = m.pixel_frame(); worst = 0.0
    red = (im[..., 0] > 200) & (im[..., 1] < 60) & (im[..., 2] < 60)
    for lo, la in pts:
        x, y = m.P(lo, la)
        px = F["ax_left"] + (x - F["extent"][0]) / (F["extent"][1] - F["extent"][0]) * F["ax_w"]
        py = F["ax_top"] + (F["extent"][3] - y) / (F["extent"][3] - F["extent"][2]) * F["ax_h"]
        yy, xx = np.nonzero(red[int(py) - 8:int(py) + 9, int(px) - 8:int(px) + 9])
        worst = max(worst, np.hypot(xx.mean() - 8 - (px - int(px)), yy.mean() - 8 - (py - int(py))))
        # the JavaScript inverse projection, replicated
        P = F["proj"]; dy = P["rho0"] - y; rho = np.sign(P["n"]) * np.hypot(x, dy)
        la2 = np.degrees(2 * np.arctan((P["R"] * P["F"] / rho) ** (1 / P["n"])) - np.pi / 2)
        lo2 = P["lon0"] + np.degrees(np.arctan2(x, dy)) / P["n"]
        worst = max(worst, 100 * abs(la2 - la), 100 * abs(lo2 - lo))
    report("hover: rendered marker within 1.5 px of the viewer's pixel, inverse projection exact",
           worst < 1.5, f"worst offset {worst:.2f} px (or 0.01 deg units)")


def test_pack_roundtrip():
    import base64
    import make_maps
    a = np.linspace(-40, 40, 600).reshape(20, 30); b = a * 100; b[3, 4] = np.nan
    o = make_maps.pack([("a", a), ("b", b)])
    q = np.frombuffer(base64.b64decode(o["b64"]), dtype="<i2").reshape(2, 20, 30)
    ra = o["offset"][0] + q[0] * o["scale"][0]
    rb = np.where(q[1] == -32768, np.nan, o["offset"][1] + q[1] * o["scale"][1])
    err = max(np.abs(ra - a).max() / 80, np.nanmax(np.abs(rb - b)) / 8000)
    report("viewer data: int16 packing round-trips within 1/65000 of range, NaN kept",
           err <= 1 / 65000 + 1e-12 and np.isnan(rb[3, 4]), f"relative error {err:.1e}")


TESTS = [test_projection_is_conformal, test_clip_breaks_lines_at_the_box,
         test_geojson_geometry_kinds, test_heights_match_standard_atmosphere,
         test_mslp_recovers_sea_level_pressure, test_interpolation_exact_at_model_levels,
         test_vorticity, test_destagger, test_wind_rotation_follows_the_meridians,
         test_every_product_renders, test_viewer_embeds_every_product,
         test_forecast_hours_match_real_snapshot_times, test_valid_time_labels,
         test_pixel_geometry_for_hover, test_pack_roundtrip]

if __name__ == "__main__":
    print("=" * 62)
    print("Forecast maps")
    print("=" * 62)
    for fn in TESTS:
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

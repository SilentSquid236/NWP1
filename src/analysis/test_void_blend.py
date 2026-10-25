"""Tests for build.void_blend (test J5, P-82).   python src/analysis/test_void_blend.py"""
import numpy as np
import build, geo

D = {"lat_min": 37.0, "lat_max": 47.5, "lon_min": -82.0, "lon_max": -66.0}
ny, nx = 30, 40
lat, lon = geo.cell_centres(D, ny, nx)
res = []
def report(name, ok, detail):
    res.append(ok); print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")
def fields(val):
    return {v: np.full((3, ny, nx), float(val)) for v in ("TMP", "RH", "UGRD", "VGRD")}
def so_at(sites):
    return {"TMP": [(0, "S%d" % i, la, lo, 280.0, 1.0) for i, (la, lo) in enumerate(sites)]}

def test_near_unchanged_far_blended():
    bg, cold = fields(10.0), fields(0.0)
    out, info = build.void_blend(bg, cold, so_at([(42.0, -74.0)]), lat, lon, D, 300.0, 600.0, 0.5)
    d = np.hypot((lat - 42.0) * 111.132, (lon + 74.0) * 111.32 * np.cos(np.radians(0.5 * (lat + 42.0))))
    near, far = d <= 300.0, d >= 600.0
    ok = np.allclose(out["UGRD"][:, near], 10.0) and np.allclose(out["UGRD"][:, far], 5.0)
    report("near a sounding unchanged; beyond d1 the w_min mix", ok,
           f"near {out['UGRD'][0][near].min():.3f}-{out['UGRD'][0][near].max():.3f}, far {out['UGRD'][0][far].min():.3f}-{out['UGRD'][0][far].max():.3f}")

def test_linear_between():
    bg, cold = fields(1.0), fields(0.0)
    out, _ = build.void_blend(bg, cold, so_at([(42.0, -74.0)]), lat, lon, D, 300.0, 600.0, 0.0)
    d = np.hypot((lat - 42.0) * 111.132, (lon + 74.0) * 111.32 * np.cos(np.radians(0.5 * (lat + 42.0))))
    mid = (d > 300.0) & (d < 600.0)
    want = (600.0 - d[mid]) / 300.0
    err = float(np.max(np.abs(out["TMP"][0][mid] - want)))
    report("weight linear in distance between d0 and d1", err < 1e-9, f"max error {err:.1e} over {int(mid.sum())} cells")

def test_no_soundings_untouched():
    bg, cold = fields(3.0), fields(0.0)
    out, info = build.void_blend(bg, cold, {}, lat, lon, D, 300.0, 600.0, 0.0)
    report("no soundings: first guess untouched", (not info["applied"]) and np.all(out["VGRD"] == 3.0), str(info))

def test_domain_only_ignores_outside():
    bg, cold = fields(2.0), fields(0.0)
    out, info = build.void_blend(bg, cold, so_at([(50.0, -90.0)]), lat, lon, D, 300.0, 600.0, 0.0, domain_only=True)
    report("domain_only: a sounding outside the box counts as none", not info["applied"], str(info))

def test_large_scale_background_follows_gradient():
    """J5b: the Barnes cold background keeps a north-south gradient the soundings show."""
    levels_pa = np.array([85000.0, 70000.0, 50000.0])
    so = {"TMP": [(k, "S", la, -74.0, t0 - 0.8 * (la - 40.0), 1.0)
                  for k, t0 in enumerate((280.0, 270.0, 255.0)) for la in (38.0, 40.0, 42.0, 44.0, 46.0)]}
    cold, _ = build.first_guess_from_soundings(so, levels_pa, (ny, nx))
    flat = float(np.ptp(cold["TMP"][0]))
    cold, _ = build.analyse_upper(cold, so, lat, lon, D, L_km=1000.0)
    north = float(cold["TMP"][0][lat > 46.0].mean()); south = float(cold["TMP"][0][lat < 38.5].mean())
    report("large-scale background: uniform mean has no gradient, Barnes L=1000 km recovers it",
           flat == 0.0 and north < south - 3.0, f"uniform ptp {flat:.2f} K; Barnes north {north:.2f} vs south {south:.2f} K")

for t in (test_large_scale_background_follows_gradient, test_near_unchanged_far_blended, test_linear_between, test_no_soundings_untouched, test_domain_only_ignores_outside):
    t()
print(f"{sum(res)}/{len(res)} passed")

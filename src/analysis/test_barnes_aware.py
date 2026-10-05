"""
Tests for barnes.barnes_increments_aware (CAM stage S5, station gap-filling).

Run:  python test_barnes_aware.py
"""
import numpy as np

import barnes

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


rng = np.random.default_rng(1)
lat_g, lon_g = np.meshgrid(np.linspace(40, 44, 37), np.linspace(-78, -72, 41), indexing="ij")
lat_o = rng.uniform(40, 44, 60); lon_o = rng.uniform(-78, -72, 60)
d = rng.normal(0, 2, 60); d[3] = np.nan
s = rng.uniform(1, 2, 60)

# 1. No factors: identical to the reference Barnes, also when chunked
a = barnes.barnes_increments(lat_g, lon_g, lat_o, lon_o, d, s, 120.0, passes=3)
b = barnes.barnes_increments_aware(lat_g, lon_g, lat_o, lon_o, d, s, 120.0, passes=3, chunk=97)
report("with no factors it equals barnes_increments (chunked, 3 passes)",
       np.abs(a - b).max() < 1e-12, f"max difference {np.abs(a - b).max():.1e}")

# 2. Elevation: a valley station's increment reaches a nearby ridge much less
lg, lo = np.array([[42.0, 42.0]]), np.array([[-75.05, -74.95]])
zg = np.array([[150.0, 1150.0]])           # valley point, ridge point (same distance)
inc = barnes.barnes_increments_aware(lg, lo, [42.0], [-75.0], [5.0], [1.5], 40.0,
                                     z_g=zg, z_o=[150.0], H_m=300.0)
report("elevation weighting: a valley station sets the valley, not the ridge 1000 m above",
       inc[0, 0] > 3.5 and inc[0, 1] < 0.05,
       f"valley point {inc[0, 0]:.2f} K, ridge point {inc[0, 1]:.3f} K (station innovation 5 K)")

# 3. Land/sea: a land station's increment over the water is reduced
# Near a station (weight >> lam) the factor barely matters; it shortens the
# reach. So compare two points 50 km away, where the weight is comparable to lam.
lo50 = np.array([[-75.6, -74.4]])
inc_l = barnes.barnes_increments_aware(lg, lo50, [42.0], [-75.0], [5.0], [1.5], 40.0,
                                       land_g=np.array([[True, False]]), land_o=[True],
                                       coast_factor=0.5)
report("land/sea weighting: 50 km from a land station, the increment over water is smaller",
       inc_l[0, 1] < 0.75 * inc_l[0, 0],
       f"land point {inc_l[0, 0]:.2f}, water point {inc_l[0, 1]:.2f}")

# 4. Unknown station elevation is not penalised
inc_n = barnes.barnes_increments_aware(lg, lo, [42.0], [-75.0], [5.0], [1.5], 40.0,
                                       z_g=zg, z_o=[np.nan], H_m=300.0)
plain = barnes.barnes_increments(lg, lo, [42.0], [-75.0], [5.0], [1.5], 40.0)
report("a station with no elevation is weighted as plain Barnes",
       np.abs(inc_n - plain).max() < 1e-12, f"difference {np.abs(inc_n - plain).max():.1e}")

print(f"\n{sum(results)}/{len(results)} passed")
raise SystemExit(0 if all(results) else 1)

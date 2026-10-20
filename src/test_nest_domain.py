"""
Test N1 plumbing: a 3 km grid nested inside a LARGER 12 km grid.

    python src/test_nest_domain.py

1. Same box: spans_for gives None and the regrid is the old one, bit for bit.
2. A field linear in latitude and longitude comes through the offset regrid
   exactly (bilinear reproduces linear fields), on centres and on faces.
3. The 12 km analysis as the 3 km first guess (build.regrid_features with
   src_domain): the same linear check, and the old path when the boxes agree.
4. domain_of inverts cell_centres; spans_for refuses a 3 km box outside.
5. config.NWP_DOMAIN: unset gives the production box; set gives the new one.
"""
import os
import subprocess
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parents[1]
for p in ("src", "src/dynamics", "src/analysis", ""):
    sys.path.insert(0, str(ROOT / p))
import nest   # noqa: E402
import geo    # noqa: E402
import build  # noqa: E402

D3 = {"lat_min": 37.0, "lat_max": 47.5, "lon_min": -82.0, "lon_max": -66.0}
DB = {"lat_min": 34.0, "lat_max": 50.5, "lon_min": -86.0, "lon_max": -62.0}
results = []


def report(name, ok, detail=""):
    results.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}" + (f"  ({detail})" if detail else ""))


def points(dom, ny, nx, sy, sx):
    """lat, lon of a grid's points: 'c' centres, 'f' low faces."""
    dla = (dom["lat_max"] - dom["lat_min"]) / ny
    dlo = (dom["lon_max"] - dom["lon_min"]) / nx
    la = dom["lat_min"] + (np.arange(ny) + (0.5 if sy == "c" else 0.0)) * dla
    lo = dom["lon_min"] + (np.arange(nx) + (0.5 if sx == "c" else 0.0)) * dlo
    return np.meshgrid(la, lo, indexing="ij")


def lin(la, lo):
    return 3.0 * la - 2.0 * lo + 7.0


# 1. same box
rng = np.random.default_rng(1)
a = rng.normal(size=(5, 97, 110))
same = nest.spans_for(D3, dict(D3)) is None
b0 = nest.regrid(a, (389, 439), "c", "f")
b1 = nest.regrid(a, (389, 439), "c", "f", spans=nest.spans_for(D3, D3))
report("same box: spans_for is None and the regrid is bit-identical", same and np.array_equal(b0, b1))

# 2. offset regrid of a linear field
ny, nx, NY, NX = 153, 165, 389, 439
sp = nest.spans_for(DB, D3)
worst = 0.0
for sy, sx in (("c", "c"), ("c", "f"), ("f", "c")):
    la_s, lo_s = points(DB, ny, nx, sy, sx)
    la_d, lo_d = points(D3, NY, NX, sy, sx)
    got = nest.regrid(lin(la_s, lo_s), (NY, NX), sy, sx, spans=sp)
    worst = max(worst, float(np.abs(got - lin(la_d, lo_d)).max()))
report("offset regrid reproduces a linear field on centres and faces (< 1e-9)", worst < 1e-9, f"max error {worst:.1e}")

# 3. regrid_features with a source domain
la_s, lo_s = geo.cell_centres(DB, ny, nx)
f = np.broadcast_to(lin(la_s, lo_s), (2, 3, ny, nx)).copy()
out = build.regrid_features(f, D3, NY, NX, src_domain=DB)
la_d, lo_d = geo.cell_centres(D3, NY, NX)
e3 = float(np.abs(out - lin(la_d, lo_d)).max())
g = rng.normal(size=(2, 3, 97, 110))
old = build.regrid_features(g, D3, NY, NX)
new = build.regrid_features(g, D3, NY, NX, src_domain=geo.domain_of(*geo.cell_centres(D3, 97, 110)))
report("regrid_features: linear field from the larger box exact (< 1e-9); same box bit-identical",
       e3 < 1e-9 and np.array_equal(old, new), f"max error {e3:.1e}")

# 4. domain_of, and a 3 km box outside the 12 km one
dd = geo.domain_of(*geo.cell_centres(DB, ny, nx))
e4 = max(abs(dd[k] - DB[k]) for k in DB)
try:
    nest.spans_for(D3, DB)
    refused = False
except ValueError:
    refused = True
report("domain_of inverts cell_centres (< 1e-9 deg); a 3 km box outside the 12 km box is refused",
       e4 < 1e-9 and refused, f"{e4:.1e} deg")

# 5. config.NWP_DOMAIN
code = f"import sys; sys.path.insert(0, {str(ROOT)!r}); import config; d = config.DOMAIN; print(d['lat_min'], d['lat_max'], d['lon_min'], d['lon_max'])"
env0 = {k: v for k, v in os.environ.items() if k != "NWP_DOMAIN"}
r0 = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=env0, capture_output=True, text=True)
r1 = subprocess.run([sys.executable, "-c", code], cwd=ROOT, env=dict(env0, NWP_DOMAIN="34,50.5,-86,-62"),
                    capture_output=True, text=True)
report("NWP_DOMAIN unset gives 37-47.5 N, 82-66 W; set gives the new box",
       r0.stdout.split() == ["37.0", "47.5", "-82.0", "-66.0"] and r1.stdout.split() == ["34.0", "50.5", "-86.0", "-62.0"],
       f"{r0.stdout.strip()} | {r1.stdout.strip()} {r1.stderr.strip()[-120:]}")

print(f"\n{sum(results)}/{len(results)} passed")
sys.exit(0 if all(results) else 1)

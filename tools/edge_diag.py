"""
Edge diagnostic for the lateral boundaries (P-76): does the model's MSLP field
develop a ring of packed isobars along the relaxation zone?

For each forecast file, by distance from the nearest edge (in cells) and by
hour: the mean MSLP change since hour 0, and the mean MSLP gradient magnitude
(hPa per 100 km), with the MSLP derived exactly as the maps derive it
(maps/derive.py: column_state + mslp). The relaxation weight at each distance
(forecast.py's Relaxation3D, width and alpha as given) is written alongside.

    python tools/edge_diag.py --out edge.csv label:path/forecast.npz [...]
"""
import argparse
import csv
import sys
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
for p in (ROOT / "src", ROOT / "src" / "dynamics", ROOT / "src" / "maps", ROOT):
    sys.path.insert(0, str(p))

from derive import column_state, mslp          # noqa: E402


def edge_distance(ny, nx):
    j, i = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    return np.minimum(np.minimum(j, ny - 1 - j), np.minimum(i, nx - 1 - i))


def relax_weight(ny, nx, width=15, alpha=0.1):
    """Relaxation3D's weights (cosine profile), from forecast.py's grid."""
    try:
        from grid import CGrid
        from forecast import Relaxation3D
        g = CGrid(nx, ny, 3000.0, 3000.0, edge_mode="replicate")
        return np.asarray(Relaxation3D(g, width=width, alpha_max=alpha).alpha2d)
    except Exception:                              # noqa: BLE001
        return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("files", nargs="+", help="label:forecast.npz")
    ap.add_argument("--out", required=True)
    ap.add_argument("--hours", default="0,1,3,6,12,18,24")
    ap.add_argument("--maxd", type=int, default=45)
    a = ap.parse_args()
    want = [float(h) for h in a.hours.split(",")]
    rows = []
    for spec in a.files:
        lab, path = spec.split(":", 1)
        z = np.load(path, allow_pickle=True)
        th, pi = z["theta"], z["pi"]
        terrain = np.asarray(z["terrain"], float)
        ny, nx = terrain.shape
        lat = np.asarray(z["lat"], float) if "lat" in z.files and z["lat"].shape != () else None
        import config
        D = config.DOMAIN
        lat0 = 0.5 * (D["lat_min"] + D["lat_max"])
        dy_km = (D["lat_max"] - D["lat_min"]) * 111.132 / ny
        dx_km = (D["lon_max"] - D["lon_min"]) * 111.320 * np.cos(np.radians(lat0)) / nx
        dist = edge_distance(ny, nx)
        w = relax_weight(ny, nx)
        times = np.asarray(z["times_s"], float) / 3600.0
        # hour 0 is not stored by forecast.py; use the first snapshot as the reference
        ms = {}
        for h in want:
            k = int(np.argmin(np.abs(times - h)))
            if abs(times[k] - h) > 0.01 and h != 0:
                continue
            p, T, Z, ps = column_state(np.asarray(th[k], float), np.asarray(pi[k], float),
                                       z["sigma"], float(z["p_top"]), terrain)
            ms[float(times[k])] = mslp(ps, T, Z, terrain) / 100.0
        tk = sorted(ms)
        ref = ms[tk[0]]
        for t in tk:
            m = ms[t]
            gy, gx = np.gradient(m, dy_km, dx_km)
            grad = np.hypot(gx, gy) * 100.0          # hPa / 100 km
            dm = m - ref
            for d in range(0, a.maxd + 1):
                sel = dist == d
                rows.append([lab, f"{t:.2f}", d, f"{dm[sel].mean():.4f}", f"{np.abs(dm[sel]).mean():.4f}",
                             f"{grad[sel].mean():.4f}", f"{np.percentile(grad[sel], 90):.4f}",
                             "" if w is None else f"{w[sel].mean():.4f}"])
            inner = dist > a.maxd
            rows.append([lab, f"{t:.2f}", "interior", f"{dm[inner].mean():.4f}",
                         f"{np.abs(dm[inner]).mean():.4f}", f"{grad[inner].mean():.4f}",
                         f"{np.percentile(grad[inner], 90):.4f}", "0"])
        print(f"{lab}: {path} {ny}x{nx}, hours {tk}")
    with open(a.out, "w", newline="") as f:
        cw = csv.writer(f)
        cw.writerow(["run", "hour", "edge_dist_cells", "mean_dmslp_hpa", "mean_abs_dmslp_hpa",
                     "mean_grad_hpa_per_100km", "p90_grad_hpa_per_100km", "relax_weight"])
        cw.writerows(rows)


if __name__ == "__main__":
    main()

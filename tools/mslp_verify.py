"""
Sea-level pressure against the verifying analyses, and the edge-ring index
(P-76).

For each arm and cycle, the forecast MSLP at 6, 12, 18 and 24 h (derived as
the maps derive it: maps/derive.py) is compared with the p_msl of the
observation analysis valid at that time (tensors_3d or tensors_3km
obs_<valid>/obs_analysis_f00.npz). Scores are kept for the relaxation zone
(edge distance < --zone cells) and the interior separately. The ring index is
the mean MSLP gradient in the zone divided by the interior's. "persistence"
as an arm path uses the cycle's own analysis p_msl.

    python tools/mslp_verify.py --tensors data/tensors_3d --csv out.csv \
        full:runs/full/obs_{s}.npz none:runs/none/obs_{s}.npz persistence:persistence
"""
import argparse
import csv
import sys
from datetime import datetime, timedelta
from pathlib import Path

import numpy as np

ROOT = Path(__file__).resolve().parent.parent
for p in (ROOT / "src", ROOT / "src" / "dynamics", ROOT / "src" / "maps", ROOT):
    sys.path.insert(0, str(p))

import config                                   # noqa: E402
from derive import column_state, mslp           # noqa: E402

STAMPS = [f"202609{d:02d}_{h:02d}" for d in range(26, 31) for h in (0, 6, 12, 18)]


def spacing(ny, nx):
    D = config.DOMAIN
    lat0 = 0.5 * (D["lat_min"] + D["lat_max"])
    return ((D["lat_max"] - D["lat_min"]) * 111.132 / ny,
            (D["lon_max"] - D["lon_min"]) * 111.320 * np.cos(np.radians(lat0)) / nx)


def edge_distance(ny, nx):
    j, i = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    return np.minimum(np.minimum(j, ny - 1 - j), np.minimum(i, nx - 1 - i))


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arms", nargs="+", help="label:path-with-{s} | label:persistence")
    ap.add_argument("--tensors", required=True, help="the analyses' tensors directory")
    ap.add_argument("--csv", required=True)
    ap.add_argument("--zone", type=int, default=15)
    ap.add_argument("--stamps", default=",".join(STAMPS))
    ap.add_argument("--leads", default="6,12,18,24")
    a = ap.parse_args()
    T = Path(a.tensors)
    leads = [int(x) for x in a.leads.split(",")]
    rows, pooled = [], {}

    def analysis(stamp):
        f = T / f"obs_{stamp}" / "obs_analysis_f00.npz"
        if not f.exists():
            return None
        return np.asarray(np.load(f, allow_pickle=True)["p_msl"], float) / 100.0

    for s in a.stamps.split(","):
        rt = datetime.strptime(s, "%Y%m%d_%H")
        an0 = analysis(s)
        for spec in a.arms:
            lab, pat = spec.split(":", 1)
            z = None
            if pat != "persistence":
                f = Path(pat.format(s=s))
                if not f.exists():
                    continue
                z = np.load(f, allow_pickle=True)
                times = np.asarray(z["times_s"], float) / 3600.0
                ter = np.asarray(z["terrain"], float)
            for L in leads:
                ver = analysis((rt + timedelta(hours=L)).strftime("%Y%m%d_%H"))
                if ver is None:
                    continue
                if z is None:
                    if an0 is None:
                        continue
                    m = an0
                else:
                    k = int(np.argmin(np.abs(times - L)))
                    if abs(times[k] - L) > 0.01:
                        continue
                    p, Tt, Z, ps = column_state(np.asarray(z["theta"][k], float), np.asarray(z["pi"][k], float),
                                                z["sigma"], float(z["p_top"]), ter)
                    m = mslp(ps, Tt, Z, ter) / 100.0
                if m.shape != ver.shape:
                    raise SystemExit(f"{lab} {s}: forecast {m.shape} vs analysis {ver.shape}")
                ny, nx = m.shape
                dist = edge_distance(ny, nx)
                dy, dx = spacing(ny, nx)
                gy, gx = np.gradient(m, dy, dx)
                g = np.hypot(gx, gy)
                zone, inner = dist < a.zone, dist >= a.zone
                ring = float(g[zone].mean() / g[inner].mean())
                e = m - ver
                for reg, sel in (("zone", zone), ("interior", inner), ("all", np.ones_like(zone))):
                    r_ = float(np.sqrt(np.mean(e[sel] ** 2)))
                    b_ = float(np.mean(e[sel]))
                    rows.append([lab, s, L, reg, f"{r_:.4f}", f"{b_:.4f}", f"{ring:.4f}"])
                    pooled.setdefault((lab, L, reg), []).append((r_ ** 2, b_, ring))
    with open(a.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["arm", "cycle", "lead_h", "region", "rmse_hpa", "bias_hpa", "ring_index"])
        w.writerows(rows)
    labs = [s.split(":", 1)[0] for s in a.arms]
    print("MSLP vs verifying analysis: RMSE / bias (hPa), ring index (zone/interior gradient); n cycles")
    for reg in ("interior", "zone", "all"):
        print(f"  {reg}")
        for L in leads:
            cells = []
            for lab in labs:
                v = pooled.get((lab, L, reg))
                if v:
                    rm = np.sqrt(np.mean([x[0] for x in v])); bi = np.mean([x[1] for x in v]); ri = np.mean([x[2] for x in v])
                    cells.append(f"{lab} {rm:.2f}/{bi:+.2f} r{ri:.2f} n{len(v)}")
            print(f"    +{L:2d} h  " + "   ".join(cells))


if __name__ == "__main__":
    main()

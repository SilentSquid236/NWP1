"""
Edge-ring profile of sea-level pressure (P-76, test E4).

mslp_verify.py measures the ring as zone/interior gradient with a fixed
15-cell zone. That cannot compare runs whose relaxation zones differ in
width. This tool uses one fixed inner region for every arm (edge distance
>= --inner cells) and reports, at each lead:

  peak   the largest mean MSLP gradient of any --band-cell band outside the
         inner region, divided by the inner region's mean gradient (1 = no
         ring; the verifying analysis itself is reported as "analysis");
  inner  MSLP RMSE / bias (hPa) against the verifying analysis, inner region;
  outer  the same outside it.

The band profile (hPa/100 km by edge distance) is written to the CSV.

    python tools/ring_profile.py --tensors data/tensors_3km --inner 60 --csv out.csv \
        w15:runs/w15/obs_{s}.npz w60:runs/w60/obs_{s}.npz --stamps 20260926_12,...
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


def spacing(ny, nx):
    D = config.DOMAIN
    lat0 = 0.5 * (D["lat_min"] + D["lat_max"])
    return ((D["lat_max"] - D["lat_min"]) * 111.132 / ny,
            (D["lon_max"] - D["lon_min"]) * 111.320 * np.cos(np.radians(lat0)) / nx)


def edge_distance(ny, nx):
    j, i = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    return np.minimum(np.minimum(j, ny - 1 - j), np.minimum(i, nx - 1 - i))


def profile(m, dist, band, inner, dy, dx):
    gy, gx = np.gradient(m, dy, dx)
    g = np.hypot(gx, gy) * 100.0                 # hPa per 100 km
    gin = float(g[dist >= inner].mean())
    bands = [float(g[(dist >= b) & (dist < b + band)].mean()) for b in range(0, inner, band)]
    return gin, bands


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arms", nargs="+", help="label:path-with-{s}")
    ap.add_argument("--tensors", required=True)
    ap.add_argument("--csv", required=True)
    ap.add_argument("--inner", type=int, default=60)
    ap.add_argument("--band", type=int, default=5)
    ap.add_argument("--stamps", required=True)
    ap.add_argument("--leads", default="1,3,6,12")
    a = ap.parse_args()
    T = Path(a.tensors)
    leads = [int(x) for x in a.leads.split(",")]

    def analysis(stamp):
        f = T / f"obs_{stamp}" / "obs_analysis_f00.npz"
        if not f.exists():
            return None
        return np.asarray(np.load(f, allow_pickle=True)["p_msl"], float) / 100.0

    rows, pooled = [], {}
    for s in a.stamps.split(","):
        rt = datetime.strptime(s, "%Y%m%d_%H")
        for L in leads:
            ver = analysis((rt + timedelta(hours=L)).strftime("%Y%m%d_%H"))
            entries = []
            if ver is not None:
                entries.append(("analysis", ver))
            for spec in a.arms:
                lab, pat = spec.split(":", 1)
                f = Path(pat.format(s=s))
                if not f.exists():
                    continue
                z = np.load(f, allow_pickle=True)
                times = np.asarray(z["times_s"], float) / 3600.0
                k = int(np.argmin(np.abs(times - L)))
                if abs(times[k] - L) > 0.01:
                    continue
                ter = np.asarray(z["terrain"], float)
                p, Tt, Z, ps = column_state(np.asarray(z["theta"][k], float), np.asarray(z["pi"][k], float),
                                            z["sigma"], float(z["p_top"]), ter)
                entries.append((lab, mslp(ps, Tt, Z, ter) / 100.0))
            for lab, m in entries:
                ny, nx = m.shape
                dist = edge_distance(ny, nx)
                dy, dx = spacing(ny, nx)
                gin, bands = profile(m, dist, a.band, a.inner, dy, dx)
                peak = max(bands) / gin
                row = {"arm": lab, "cycle": s, "lead_h": L, "inner_grad": round(gin, 3),
                       "peak_ratio": round(peak, 3)}
                if ver is not None and lab != "analysis":
                    din, dout = (m - ver)[dist >= a.inner], (m - ver)[dist < a.inner]
                    row.update(inner_rmse=round(float(np.sqrt(np.mean(din ** 2))), 3),
                               inner_bias=round(float(din.mean()), 3),
                               outer_rmse=round(float(np.sqrt(np.mean(dout ** 2))), 3),
                               outer_bias=round(float(dout.mean()), 3))
                for b, v in zip(range(0, a.inner, a.band), bands):
                    row[f"g{b:03d}"] = round(v, 3)
                rows.append(row)
                pooled.setdefault((lab, L), []).append(row)
    keys = sorted({k for r in rows for k in r}, key=lambda k: (k[0] == "g", k))
    with open(a.csv, "w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=keys)
        w.writeheader()
        w.writerows(rows)
    print(f"ring profile: inner region edge distance >= {a.inner} cells, bands of {a.band}")
    print("arm        lead  n  peak  inner_grad  inner RMSE/bias   outer RMSE/bias   band profile (hPa/100 km)")
    for (lab, L), rs in sorted(pooled.items(), key=lambda x: (x[0][1], x[0][0])):
        mean = lambda k: np.mean([r[k] for r in rs if k in r]) if any(k in r for r in rs) else np.nan
        prof = " ".join(f"{mean(f'g{b:03d}'):4.1f}" for b in range(0, a.inner, a.band))
        print(f"{lab:10s} {L:4d} {len(rs):2d} {mean('peak_ratio'):5.2f} {mean('inner_grad'):6.2f}"
              f"      {mean('inner_rmse'):5.2f}/{mean('inner_bias'):+5.2f}     "
              f"{mean('outer_rmse'):5.2f}/{mean('outer_bias'):+5.2f}    {prof}")


if __name__ == "__main__":
    main()

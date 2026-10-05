"""
Compare verification archives of several forecast arms on PAIRED matches,
by lead, day/night and terrain class (CAM stage S5, campaign verification).

    python tools/cam_compare.py --grid <3 km forecast.npz> --csv out.csv \
        cam3:<archive> nh12:<archive> ad12:<archive> [am12:<archive>]

Each archive is a verification tree <root>/<YYYYMMDD_HH>/matches.jsonl as
written by verify.py. Only matches present in EVERY arm are scored (same
cycle, station, variable, valid time and lead), so the arms are compared on
the same observations.

Classes, fixed before any result was seen:
- day / night: the sun above or below the horizon at the station and valid
  time (NOAA low-precision solar position);
- terrain, from the 3 km grid around the station's nearest grid point
  (relief = max - min of the terrain within +-3 cells, a 21 km box):
  flat < 150 m, hilly 150-400 m, mountain > 400 m.

  CHANGED 2026-10-05, after the interim run. The relief was first taken
  from the forecast file's terrain. That terrain is slope-limited (0.0086),
  so no 21 km box can hold more than ~155 m: 365 of 367 stations came out
  flat and none mountain. --terrain-npz now takes the unlimited ETOPO
  terrain on the same grid (data/static/terrain_etopo_389x439.npz). The
  class edges are unchanged.

Writes one CSV row per (variable, arm, lead, period, terrain): n, bias, rmse;
lead -1 means all leads pooled. Prints the pooled tables.
"""
import argparse
import collections
import csv
import glob
import json
import math
import os
from datetime import datetime

import numpy as np

VARS = ("TMP", "UGRD", "VGRD")
RELIEF_EDGES = (150.0, 400.0)
RELIEF_HALF = 3


def key(cycle, m):
    return (cycle, m["station"], m["variable"], m["time"], int(round(m["lead_hours"])))


def solar_elevation(lat, lon, when):
    """Degrees; NOAA low-precision formula (good to ~0.5 deg)."""
    doy = when.timetuple().tm_yday
    hour = when.hour + when.minute / 60.0 + when.second / 3600.0
    g = 2 * math.pi / 365.0 * (doy - 1 + (hour - 12) / 24.0)
    decl = (0.006918 - 0.399912 * math.cos(g) + 0.070257 * math.sin(g) - 0.006758 * math.cos(2 * g)
            + 0.000907 * math.sin(2 * g) - 0.002697 * math.cos(3 * g) + 0.00148 * math.sin(3 * g))
    eqt = 229.18 * (0.000075 + 0.001868 * math.cos(g) - 0.032077 * math.sin(g)
                    - 0.014615 * math.cos(2 * g) - 0.040849 * math.sin(2 * g))
    tst = hour * 60 + eqt + 4 * lon
    ha = math.radians(tst / 4.0 - 180.0)
    la = math.radians(lat)
    cz = math.sin(la) * math.sin(decl) + math.cos(la) * math.cos(decl) * math.cos(ha)
    return math.degrees(math.asin(max(-1.0, min(1.0, cz))))


def terrain_classes(grid_npz, stations, terrain_npz=None):
    """station -> (elevation_m, relief_m, class) from the 3 km grid."""
    z = np.load(grid_npz, allow_pickle=True)
    lat, lon, ter = (np.asarray(z[k], float) for k in ("lat", "lon", "terrain"))
    if terrain_npz:
        zt = np.load(terrain_npz, allow_pickle=True)
        name = "terrain" if "terrain" in zt else [k for k in zt.keys() if np.ndim(zt[k]) == 2][0]
        raw = np.asarray(zt[name], float)
        if raw.shape != ter.shape:
            raise SystemExit(f"--terrain-npz shape {raw.shape} != grid {ter.shape}")
        print(f"terrain for the classes: {terrain_npz} [{name}], {raw.min():.0f}-{raw.max():.0f} m "
              f"(grid file {ter.min():.0f}-{ter.max():.0f} m)")
        ter = raw
    if lat.ndim == 1:
        lon, lat = np.meshgrid(lon, lat)
    ny, nx = ter.shape
    flat_lat, flat_lon = lat.ravel(), lon.ravel()
    out = {}
    for s, (sla, slo) in stations.items():
        d2 = (flat_lat - sla) ** 2 + ((flat_lon - slo) * math.cos(math.radians(sla))) ** 2
        j, i = divmod(int(np.argmin(d2)), nx)
        box = ter[max(0, j - RELIEF_HALF):j + RELIEF_HALF + 1, max(0, i - RELIEF_HALF):i + RELIEF_HALF + 1]
        rel = float(box.max() - box.min())
        cls = "flat" if rel < RELIEF_EDGES[0] else "hilly" if rel < RELIEF_EDGES[1] else "mountain"
        out[s] = (float(ter[j, i]), rel, cls)
    return out


def stats(e):
    e = np.asarray(e, float)
    return len(e), float(e.mean()) if len(e) else float("nan"), \
        float(np.sqrt(np.mean(e ** 2))) if len(e) else float("nan")


def main():
    p = argparse.ArgumentParser()
    p.add_argument("arms", nargs="+", help="label:archive_root")
    p.add_argument("--grid", required=True, help="a 3 km forecast.npz (lat, lon, terrain)")
    p.add_argument("--csv", required=True)
    p.add_argument("--stations-csv", default=None, help="write the station classes here")
    p.add_argument("--terrain-npz", default=None,
                   help="unlimited terrain on the same grid for the relief classes")
    a = p.parse_args()

    labels, roots = [], []
    for spec in a.arms:
        lab, root = spec.split(":", 1)
        labels.append(lab)
        roots.append(root)
    cyc_sets = [set(os.path.basename(os.path.dirname(f))
                    for f in glob.glob(os.path.join(r, "*", "matches.jsonl"))) for r in roots]
    cycles = sorted(set.intersection(*cyc_sets))
    print(f"arms {labels}; {len(cycles)} cycles in all arms: {cycles[0] if cycles else '-'} .. "
          f"{cycles[-1] if cycles else '-'}")

    # paired matches: keep only the error (and, for the first arm, the position)
    per_arm, where = [], {}
    for ai, r in enumerate(roots):
        d = {}
        for c in cycles:
            with open(os.path.join(r, c, "matches.jsonl")) as fh:
                for line in fh:
                    m = json.loads(line)
                    if m["variable"] in VARS:
                        k = key(c, m)
                        d[k] = float(m["error"])
                        if ai == 0:
                            where[k] = (m["station"], m["time"], m["lat"], m["lon"])
        per_arm.append(d)
    common = set.intersection(*(set(d) for d in per_arm))
    sizes = [len(d) for d in per_arm]
    print(f"matches per arm {sizes}; paired in all arms {len(common)}")

    stations = {}
    for k in common:
        s, _, la, lo = where[k]
        stations.setdefault(s, (la, lo))
    tc = terrain_classes(a.grid, stations, a.terrain_npz)
    rel = np.array([v[1] for v in tc.values()])
    print("station relief (m) quantiles 10/50/90/max: " +
          " / ".join(f"{q:.0f}" for q in np.percentile(rel, [10, 50, 90, 100])))
    if a.stations_csv:
        with open(a.stations_csv, "w", newline="") as f:
            w = csv.writer(f)
            w.writerow(["station", "lat", "lon", "elev_m", "relief_m", "class"])
            for s, (la, lo) in sorted(stations.items()):
                w.writerow([s, la, lo, *tc[s]])
    ncls = collections.Counter(v[2] for v in tc.values())
    print(f"stations {len(stations)}: " + ", ".join(f"{k} {ncls[k]}" for k in ("flat", "hilly", "mountain")))

    sun_cache = {}
    groups = collections.defaultdict(lambda: [[] for _ in labels])
    cyc_groups = collections.defaultdict(lambda: [[] for _ in labels])
    for k in common:
        st, tm, la, lo = where[k]
        var, lead = k[2], k[4]
        sk = (st, tm)
        if sk not in sun_cache:
            when = datetime.fromisoformat(tm.replace("Z", "+00:00")).replace(tzinfo=None)
            sun_cache[sk] = "day" if solar_elevation(la, lo, when) > 0 else "night"
        per = sun_cache[sk]
        cls = tc[st][2]
        for ai, d in enumerate(per_arm):
            e = d[k]
            for L in (lead, -1):
                for pr in (per, "all"):
                    for tcl in (cls, "all"):
                        groups[(var, L, pr, tcl)][ai].append(e)
            cyc_groups[(var, k[0])][ai].append(e)

    with open(a.csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["variable", "arm", "lead", "period", "terrain", "n", "bias", "rmse"])
        for (var, L, pr, tcl), lists in sorted(groups.items()):
            for lab, e in zip(labels, lists):
                n, b, r = stats(e)
                w.writerow([var, lab, L, pr, tcl, n, f"{b:.4f}", f"{r:.4f}"])
        for (var, c), lists in sorted(cyc_groups.items()):
            for lab, e in zip(labels, lists):
                n, b, r = stats(e)
                w.writerow([var, lab, f"cycle:{c}", "all", "all", n, f"{b:.4f}", f"{r:.4f}"])

    ref = labels[0]
    for var in VARS:
        print(f"\n{var}: RMSE by lead (paired), and {labels[1:]} minus {ref}")
        print("  lead      n  " + "  ".join(f"{l:>8s}" for l in labels))
        for L in sorted({k[1] for k in groups if k[0] == var and k[1] >= 0}):
            lists = groups[(var, L, "all", "all")]
            rm = [stats(e)[2] for e in lists]
            print(f"  {L:4d}  {len(lists[0]):6d}  " + "  ".join(f"{x:8.3f}" for x in rm))
        for pr in ("all", "day", "night"):
            for tcl in ("all", "flat", "hilly", "mountain"):
                lists = groups.get((var, -1, pr, tcl))
                if not lists or not lists[0]:
                    continue
                s = [stats(e) for e in lists]
                print(f"  pooled {pr:5s} {tcl:8s} n {s[0][0]:6d}  " +
                      "  ".join(f"{lab} {b:+.2f}/{r:.3f}" for lab, (_, b, r) in zip(labels, s)))
        diffs = []
        for c in cycles:
            lists = cyc_groups.get((var, c))
            if lists and lists[0]:
                diffs.append([stats(e)[2] for e in lists])
        if diffs:
            d = np.array(diffs)
            for ai in range(1, len(labels)):
                better = int(np.sum(d[:, 0] < d[:, ai]))
                print(f"  per cycle: {ref} better than {labels[ai]} in {better} of {len(d)} cycles; "
                      f"mean RMSE diff {ref} - {labels[ai]} {np.mean(d[:, 0] - d[:, ai]):+.3f}")


if __name__ == "__main__":
    main()

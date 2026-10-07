
"""Test AJ part 2: lowest-level and profile wind excess over the analysis, per variant.
usage: aj_cmp.py VAR1 VAR2 ...   (each VAR is a dir holding obs_YYYYMMDD_HH.npz forecasts)"""
import numpy as np, datetime as dt, os, sys
D = "/data5/pierce/AINWP/data/tensors_3d"
name = lambda t: f"obs_{t:%Y%m%d_%H}"
topo = np.load("/data5/pierce/AINWP/data/static/terrain_etopo_97x110.npz")["terrain"]; land = topo > 0
ny, nx = land.shape; jj, ii = np.mgrid[0:ny, 0:nx]
interior = np.minimum(np.minimum(jj, ny-1-jj), np.minimum(ii, nx-1-ii)) >= 12
regions = {"all": np.ones_like(land), "int_land": interior & land, "int_sea": interior & ~land}
cache = {}
def ana(t):
    if t not in cache:
        f = f"{D}/{name(t)}/forecast_persist.npz"
        cache[t] = None if not os.path.exists(f) else (lambda z: (z["u"][0], z["v"][0]))(np.load(f))
    return cache[t]
for var in sys.argv[1:]:
    files = sorted(x for x in os.listdir(var) if x.startswith("obs_") and x.endswith(".npz"))
    acc = {L: {"rat": {r: [] for r in regions}, "prof": []} for L in (1, 6, 12, 18, 24)}
    for fn in files:
        C = dt.datetime.strptime(fn[4:15], "%Y%m%d_%H")
        z = np.load(os.path.join(var, fn)); hrs = list(np.rint(z["times_s"]/3600).astype(int))
        for L in acc:
            a = ana(C + dt.timedelta(hours=L)) if L > 1 else ana(C)
            if a is None or L not in hrs: continue
            n = hrs.index(L); sm = np.hypot(z["u"][n], z["v"][n]); sa = np.hypot(*a)
            for r, m in regions.items(): acc[L]["rat"][r].append(sm[-1][m].mean() / sa[-1][m].mean())
            acc[L]["prof"].append([(sm[k] - sa[k])[interior].mean() for k in range(sm.shape[0])])
    print(f"== {var}  ({len(files)} forecasts)")
    for L, d in acc.items():
        if not d["prof"]: continue
        P = np.mean(d["prof"], 0)
        print(f"  L={L:2d} n={len(d['prof']):2d} ratio " + " ".join(f"{r} {np.mean(v):.3f}" for r, v in d["rat"].items())
              + "  | excess lowest {:+.2f}  2nd {:+.2f}  lv10 {:+.2f}  lv13-18 {:+.2f}".format(P[-1], P[-2], P[9], P[12:18].mean()))

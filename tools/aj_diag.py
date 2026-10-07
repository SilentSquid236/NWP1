
"""Test AJ diagnostic: where does the model's wind grow relative to the analysis?
Model (forecast.npz of cycle C, lead L) vs analysis valid at C+L (forecast_persist.npz
of cycle C+L, which holds that cycle's analysis converted to sigma)."""
import numpy as np, datetime as dt, os, json
D = "/data5/pierce/AINWP/data/tensors_3d"
cyc = [dt.datetime(2026,9,26)+dt.timedelta(hours=6*i) for i in range(20)]
name = lambda t: f"obs_{t:%Y%m%d_%H}"
st = np.load("/data5/pierce/AINWP/data/static/terrain_etopo_97x110.npz")
print("static keys", list(st.keys()))
k = [x for x in st.keys() if 'terr' in x or 'elev' in x or 'height' in x][0]
topo = st[k]; land = topo > 0
ny, nx = land.shape
jj, ii = np.mgrid[0:ny, 0:nx]
edge = np.minimum(np.minimum(jj, ny-1-jj), np.minimum(ii, nx-1-ii))
interior = edge >= 12
regions = {"all": np.ones_like(land), "edge<12": ~interior, "int_land": interior & land,
           "int_sea": interior & ~land, "int_land_hi": interior & land & (topo > 500)}
def spd(u, v): return np.hypot(u, v)
cache = {}
def ana(t):
    if t not in cache:
        f = f"{D}/{name(t)}/forecast_persist.npz"
        if not os.path.exists(f): cache[t] = None
        else:
            z = np.load(f); cache[t] = (z["u"][0], z["v"][0])
    return cache[t]
rows = []
prof = {L: [] for L in (6,12,18,24)}
for C in cyc:
    z = np.load(f"{D}/{name(C)}/forecast.npz")
    hrs = np.rint(z["times_s"]/3600).astype(int)
    a0 = ana(C)
    for L in (6,12,18,24):
        a = ana(C+dt.timedelta(hours=L))
        if a is None or L not in hrs: continue
        n = list(hrs).index(L)
        um, vm = z["u"][n], z["v"][n]
        sm, sa, sp = spd(um, vm), spd(*a), spd(*a0)
        r = {"cycle": name(C), "L": L}
        for rn, m in regions.items():
            r[rn+"_mod"] = float(sm[-1][m].mean()); r[rn+"_ana"] = float(sa[-1][m].mean()); r[rn+"_per"] = float(sp[-1][m].mean())
        rows.append(r)
        prof[L].append([float((sm[k_][interior]).mean() - (sa[k_][interior]).mean()) for k_ in range(sm.shape[0])])
    # hourly lowest-level interior-land speed for the model
print("rows", len(rows))
import csv
with open("aj_rows.csv", "w", newline="") as fh:
    w = csv.DictWriter(fh, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
out = {}
for L in (6,12,18,24):
    R_ = [r for r in rows if r["L"] == L]
    out[L] = {rn: (np.mean([r[rn+"_mod"] for r in R_]), np.mean([r[rn+"_ana"] for r in R_]), np.mean([r[rn+"_per"] for r in R_])) for rn in regions}
    print(f"L={L:2d} n={len(R_)}")
    for rn, (m, a, p) in out[L].items():
        print(f"   {rn:12s} model {m:6.2f}  analysis {a:6.2f}  persist {p:6.2f}  model/ana {m/a:5.3f}")
    P = np.array(prof[L]).mean(0)
    print("   interior model-analysis by level (top..bottom):", " ".join(f"{x:+.2f}" for x in P))

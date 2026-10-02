"""Test AE summary: compare verification archives (model runs and persistence)."""
import json, sys, glob, collections
import numpy as np

def load(arch):
    f = glob.glob(arch + "/*/matches.jsonl")
    if not f: return None
    rows = [json.loads(l) for l in open(f[0])]
    return rows

def per_lead_T(rows):
    d = collections.defaultdict(list)
    for m in rows:
        if m["variable"] == "TMP": d[int(round(m["lead_hours"]))].append(m["error"])
    return {k: (np.mean(v), np.sqrt(np.mean(np.square(v))), len(v)) for k, v in sorted(d.items())}

def interior_day_speed_bias(rows, lat2, lon2, edge):
    w = collections.defaultdict(dict)
    for m in rows:
        if m["variable"] in ("UGRD", "VGRD"):
            w[(m["station"], m["valid_time"][:13])][m["variable"]] = m
    out = []
    for (st, vt), r in w.items():
        if len(r) < 2: continue
        u, v = r["UGRD"], r["VGRD"]
        k = np.argmin((lat2 - u["lat"]) ** 2 + ((lon2 - u["lon"]) * np.cos(np.radians(u["lat"]))) ** 2)
        if edge.flat[k] < 15: continue
        hour = int(vt[11:13]); solar = (hour + u["lon"] / 15.0) % 24
        if 10 <= solar <= 17:
            out.append(np.hypot(u["forecast"], v["forecast"]) - np.hypot(u["observation"], v["observation"]))
    return (float(np.mean(out)), len(out)) if out else (float("nan"), 0)

if __name__ == "__main__":
    grid_npz = sys.argv[1]
    z = np.load(grid_npz); lat2, lon2 = z["lat"], z["lon"]
    ny, nx = lat2.shape
    jj, ii = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    edge = np.minimum(np.minimum(jj, ny - 1 - jj), np.minimum(ii, nx - 1 - ii))
    for spec in sys.argv[2:]:          # case:base_archive:heat_archive:persist_archive
        case, base, heat, pers = spec.split(":")
        R = {n: load(a) for n, a in (("base", base), ("heat", heat), ("persist", pers))}
        if any(v is None for v in R.values()):
            print(f"== {case}: missing archive(s): {[n for n, v in R.items() if v is None]}"); continue
        T = {n: per_lead_T(r) for n, r in R.items()}
        leads = sorted(set(T["base"]) & set(T["heat"]) & set(T["persist"]))
        print(f"== {case}  (T bias / RMSE, K)")
        print("lead   base            heat            persist         heat-base  heat-persist")
        for L in leads:
            b, h, p = T["base"][L], T["heat"][L], T["persist"][L]
            print(f"{L:4d}  {b[0]:+6.2f}/{b[1]:5.2f}   {h[0]:+6.2f}/{h[1]:5.2f}   {p[0]:+6.2f}/{p[1]:5.2f}   {h[1]-b[1]:+6.2f}    {h[1]-p[1]:+6.2f}")
        mean = {n: float(np.mean([T[n][L][1] for L in leads])) for n in T}
        swing = {n: float(max(T[n][L][0] for L in leads) - min(T[n][L][0] for L in leads)) for n in T}
        wins = sum(T["heat"][L][1] < T["persist"][L][1] for L in leads)
        sb = {n: interior_day_speed_bias(r, lat2, lon2, edge) for n, r in R.items()}
        print(f"  mean RMSE over leads: base {mean['base']:.2f}, heat {mean['heat']:.2f} ({mean['heat']-mean['base']:+.2f}), persist {mean['persist']:.2f}")
        print(f"  bias swing: base {swing['base']:.2f}, heat {swing['heat']:.2f} ({(swing['heat']/swing['base']-1)*100:+.0f} %), persist {swing['persist']:.2f}")
        print(f"  heat beats persistence (T RMSE) at {wins} of {len(leads)} leads")
        print(f"  interior daytime speed bias (local solar 10-17 h): base {sb['base'][0]:+.2f} (n {sb['base'][1]}), heat {sb['heat'][0]:+.2f}, persist {sb['persist'][0]:+.2f}")

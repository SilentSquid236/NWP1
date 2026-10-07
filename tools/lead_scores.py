
import json, glob, os, collections, numpy as np, csv, sys
V = "/data5/pierce/AINWP/data/verification_tests"
arch = {"persistence": "adp_w", "AD (dry, no sun)": "ad_w", "AK (land surface)": "ak/sim", "AM (candidate)": "am/sim"}
rows = []
for name, a in arch.items():
    acc = collections.defaultdict(list)
    for f in glob.glob(f"{V}/{a}/*/matches.jsonl"):
        for line in open(f):
            m = json.loads(line)
            if m.get("pressure") is not None or m["variable"] not in ("TMP", "UGRD", "VGRD"): continue
            acc[(m["variable"], int(round(m["lead_hours"])))].append(m["error"])
    for (var, L), e in sorted(acc.items()):
        e = np.asarray(e); rows.append([name, var, L, len(e), float(e.mean()), float(np.sqrt((e**2).mean()))])
with open("lead_scores.csv", "w", newline="") as fh:
    w = csv.writer(fh); w.writerow(["config", "variable", "lead_h", "n", "bias", "rmse"]); w.writerows(rows)
print(len(rows))

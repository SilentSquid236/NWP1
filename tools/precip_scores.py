"""
Hourly precipitation scores (S7): categorical skill of APCP1H pairs.

For each arm (label:archive), reads <archive>/<cycle>/matches.jsonl and keeps
the APCP1H pairs (ASOS p01i against the forecast's rain in the hour ending at
the snapshot; src/verify.py). Pairs are matched across arms by (cycle,
station, valid time), so every arm is scored on the same hours.

For each threshold (mm/h) it counts hits a, false alarms b, misses c and
correct negatives d, and reports
  frequency bias  (a+b)/(a+c)           1 = forecasts rain as often as it falls
  POD             a/(a+c)
  FAR             b/(a+b)
  ETS             (a-ar)/(a+b+c-ar), ar = (a+b)(a+c)/n   0 = no skill over chance
plus the mean forecast and observed hourly rain and the RMSE.

A model that never rains (persistence of the initial state) scores bias 0,
POD 0 and ETS 0: the reference any rain forecast has to beat.

    python tools/precip_scores.py cam3M:runs/v_cam3M nh12:runs/v_nh12 pers:runs/v_pers \
        --thresholds 0.25,1,2.5 --csv out.csv
"""
import argparse
import csv
import json
from pathlib import Path

import numpy as np


def load(archive):
    out = {}
    for f in sorted(Path(archive).glob("*/matches.jsonl")):
        cyc = f.parent.name
        with open(f) as fh:
            for line in fh:
                m = json.loads(line)
                if m.get("variable") != "APCP1H":
                    continue
                out[(cyc, m["station"], m["time"])] = (m["forecast"], m["observation"], m.get("lead_hours"))
    return out


def scores(f, o, thr):
    fy, oy = f >= thr, o >= thr
    a = int((fy & oy).sum()); b = int((fy & ~oy).sum()); c = int((~fy & oy).sum()); d = int((~fy & ~oy).sum())
    n = a + b + c + d
    ar = (a + b) * (a + c) / n if n else 0.0
    nan = float("nan")
    return {"n": n, "obs_events": a + c, "hits": a, "false_alarms": b, "misses": c,
            "freq_bias": (a + b) / (a + c) if a + c else nan,
            "pod": a / (a + c) if a + c else nan,
            "far": b / (a + b) if a + b else nan,
            "ets": (a - ar) / (a + b + c - ar) if (a + b + c - ar) else nan}


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("arms", nargs="+", help="label:archive_dir")
    ap.add_argument("--thresholds", default="0.25,1,2.5")
    ap.add_argument("--csv")
    a = ap.parse_args()
    arms = [s.split(":", 1) for s in a.arms]
    data = {lab: load(p) for lab, p in arms}
    keys = set.intersection(*[set(d) for d in data.values()]) if data else set()
    keys = sorted(keys)
    print(f"APCP1H pairs per arm: {[len(d) for d in data.values()]}; paired in all arms: {len(keys)}")
    if not keys:
        return
    obs = np.array([data[arms[0][0]][k][1] for k in keys])
    print(f"observed: mean {obs.mean():.3f} mm/h, wet hours (>= 0.25 mm) {int((obs >= 0.25).sum())}, "
          f"max {obs.max():.1f} mm/h")
    rows = []
    print(f"{'arm':10s} {'thr':>5s} {'events':>6s} {'bias':>6s} {'POD':>6s} {'FAR':>6s} {'ETS':>6s}   "
          f"mean fc  RMSE")
    for lab, _ in arms:
        f = np.array([data[lab][k][0] for k in keys])
        for thr in [float(x) for x in a.thresholds.split(",")]:
            s = scores(f, obs, thr)
            s.update(arm=lab, threshold_mm=thr, mean_fc=float(f.mean()), mean_obs=float(obs.mean()),
                     rmse=float(np.sqrt(np.mean((f - obs) ** 2))))
            rows.append(s)
            print(f"{lab:10s} {thr:5.2f} {s['obs_events']:6d} {s['freq_bias']:6.2f} {s['pod']:6.2f} "
                  f"{s['far']:6.2f} {s['ets']:6.3f}   {s['mean_fc']:.3f}  {s['rmse']:.3f}")
    if a.csv:
        with open(a.csv, "w", newline="") as fh:
            w = csv.DictWriter(fh, fieldnames=list(rows[0]))
            w.writeheader()
            w.writerows(rows)


if __name__ == "__main__":
    main()

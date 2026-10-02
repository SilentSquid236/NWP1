"""
Surface scores split by time of day (tests AJ and AK).

    python tools/phase_scores.py ARCHIVE [ARCHIVE ...]

Each ARCHIVE is a verification archive (<stamp>/matches.jsonl). The valid
hour decides the phase, in UTC for the Northeast (EDT = UTC-4):
day 14-22 UTC (10-18 EDT), night 02-10 UTC (22-06 EDT), and transition
otherwise. Printed per phase and lead group:
- 2 m temperature: bias and RMSE;
- 10 m wind: the ratio of the mean forecast speed to the mean observed speed
  over matched u-v pairs, and the vector RMSE.
"""
import collections
import glob
import json
import os
import sys

import numpy as np


def phase(hour):
    if 14 <= hour <= 22:
        return "day"
    if 2 <= hour <= 10:
        return "night"
    return "trans"


def lead_group(L):
    return "L01-06" if L <= 6 else ("L07-12" if L <= 12 else ("L13-18" if L <= 18 else "L19-24"))


def scores(root):
    T = collections.defaultdict(list)
    W = collections.defaultdict(lambda: [[], [], []])     # fspd, ospd, vec err^2
    for f in glob.glob(os.path.join(root, "*", "matches.jsonl")):
        wind = {}
        for line in open(f):
            m = json.loads(line)
            if m.get("pressure") is not None:
                continue
            L = int(round(m["lead_hours"]))
            ph = phase(int(m["valid_time"][11:13]))
            if m["variable"] == "TMP":
                for k in ((ph, lead_group(L)), (ph, "all"), ("all", lead_group(L)), ("all", "all")):
                    T[k].append(m["error"])
            elif m["variable"] in ("UGRD", "VGRD"):
                wind.setdefault((m["station"], m["valid_time"][:16], L), {})[m["variable"]] = \
                    (m["forecast"], m["observation"])
        for (st, vt, L), d in wind.items():
            if len(d) < 2:
                continue
            fu, ou = d["UGRD"]; fv, ov = d["VGRD"]
            ph = phase(int(vt[11:13]))
            for k in ((ph, lead_group(L)), (ph, "all"), ("all", lead_group(L)), ("all", "all")):
                W[k][0].append(np.hypot(fu, fv)); W[k][1].append(np.hypot(ou, ov))
                W[k][2].append((fu - ou) ** 2 + (fv - ov) ** 2)
    return T, W


def main(roots):
    groups = ("L01-06", "L07-12", "L13-18", "L19-24", "all")
    for root in roots:
        T, W = scores(root)
        n = len(glob.glob(os.path.join(root, "*", "matches.jsonl")))
        print(f"== {root} ({n} cycles)")
        print("  2 m T bias / RMSE (K)")
        for ph in ("day", "night", "trans", "all"):
            cells = []
            for g in groups:
                e = np.asarray(T.get((ph, g), []))
                cells.append(f"{g} {e.mean():+5.2f}/{np.sqrt((e ** 2).mean()):4.2f}" if e.size else f"{g} -")
            print(f"    {ph:5s} " + "  ".join(cells))
        print("  10 m wind speed ratio (forecast/observed) / vector RMSE (m/s)")
        for ph in ("day", "night", "trans", "all"):
            cells = []
            for g in groups:
                fs, os_, e2 = (np.asarray(x) for x in W.get((ph, g), [[], [], []]))
                cells.append(f"{g} {fs.mean() / os_.mean():4.2f}/{np.sqrt(e2.mean()):4.2f}" if fs.size else f"{g} -")
            print(f"    {ph:5s} " + "  ".join(cells))


if __name__ == "__main__":
    if len(sys.argv) < 2:
        raise SystemExit(__doc__)
    main(sys.argv[1:])

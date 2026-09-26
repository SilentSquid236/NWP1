#!/usr/bin/env python3
"""
Verification scores by variable and forecast hour, for one or two archives.

    python tools/score_by_lead.py ARCHIVE_A [ARCHIVE_B] [--max-lead 24]

Each ARCHIVE is a verification root written by `src/verify.py --archive`
(it holds <run>/matches.jsonl). With two archives the scores are printed
side by side, with B's RMSE minus A's. Use it to compare two forecasts of
the same case, e.g. with and without a setting. Only snapshots within 0.05 h
of a whole forecast hour are scored, so a run written every 15 minutes
and one written hourly compare like for like. Needs only NumPy.
"""
import argparse
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def load(root):
    recs = []
    for f in sorted(Path(root).glob("*/matches.jsonl")):
        with open(f, encoding="utf-8") as fh:
            recs.extend(json.loads(line) for line in fh if line.strip())
    if not recs:
        raise SystemExit(f"no matches.jsonl under {root}")
    table = defaultdict(list)
    for r in recs:
        lead = float(r["lead_hours"])
        if abs(lead - round(lead)) > 0.05:      # whole forecast hours only
            continue
        table[(r["variable"], round(lead))].append(
            float(r["forecast"]) - float(r["observation"]))
    return table


def stats(errs):
    e = np.asarray(errs)
    return len(e), float(e.mean()), float(np.sqrt((e ** 2).mean()))


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[1])
    ap.add_argument("a")
    ap.add_argument("b", nargs="?")
    ap.add_argument("--max-lead", type=float, default=48)
    a = ap.parse_args()
    A = load(a.a)
    B = load(a.b) if a.b else None
    for var in sorted({k[0] for k in A}):
        print(f"\n{var}")
        head = f"  {'lead':>4} {'n':>5} {'bias A':>8} {'rmse A':>8}"
        if B:
            head += f" {'n':>5} {'bias B':>8} {'rmse B':>8} {'B - A':>8}"
        print(head)
        for lead in sorted({k[1] for k in A if k[0] == var} |
                           ({k[1] for k in B if k[0] == var} if B else set())):
            if lead > a.max_lead:
                continue
            line = f"  {lead:4d}"
            sa = stats(A[(var, lead)]) if (var, lead) in A else None
            line += (f" {sa[0]:5d} {sa[1]:+8.2f} {sa[2]:8.2f}" if sa else f" {'':5} {'':8} {'':8}")
            if B:
                sb = stats(B[(var, lead)]) if (var, lead) in B else None
                line += (f" {sb[0]:5d} {sb[1]:+8.2f} {sb[2]:8.2f}" if sb else f" {'':5} {'':8} {'':8}")
                if sa and sb:
                    line += f" {sb[2] - sa[2]:+8.2f}"
            print(line)


if __name__ == "__main__":
    main()

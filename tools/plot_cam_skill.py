"""
Figure: 3 km against 12 km skill by lead (CAM stage S5, test S5g).

    python tools/plot_cam_skill.py cam_compare.csv out.png [--title-cycles N]

Reads the CSV written by tools/cam_compare.py (paired matches only). Three
panels: 2 m temperature, 10 m u, 10 m v RMSE by lead, all stations, day and
night pooled. The 3 km arm is the focal series.
"""
import argparse
import csv
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ARMS = [("cam3", "3 km NH (dry)"), ("nh12", "12 km NH (dry)"),
        ("ad12", "12 km driver (hydrostatic, test AD)"),
        ("am12", "12 km production (land surface)")]
COLORS = {"cam3": "#C0392B", "nh12": "#2E86C1", "ad12": "#7F8C8D", "am12": "#6C3483"}
PANELS = [("TMP", "2 m temperature RMSE (K)"), ("UGRD", "10 m u RMSE (m/s)"), ("VGRD", "10 m v RMSE (m/s)")]


def load(path):
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r["period"] != "all" or r["terrain"] != "all" or r["lead"].startswith("cycle"):
                continue
            out[(r["variable"], r["arm"], int(r["lead"]))] = (int(r["n"]), float(r["rmse"]))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("csv"); p.add_argument("png")
    p.add_argument("--cycles", type=int, required=True, help="number of cycles (for the caption line)")
    a = p.parse_args()
    d = load(a.csv)
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.9), layout="constrained")
    for ax, (var, ylab) in zip(axes, PANELS):
        for arm, label in ARMS[::-1]:
            leads = sorted(L for (v, ar, L) in d if v == var and ar == arm and L >= 0)
            if not leads:
                continue
            y = [d[(var, arm, L)][1] for L in leads]
            focal = arm == "cam3"
            ax.plot(leads, y, color=COLORS[arm], lw=1.8 if focal else 1.0,
                    alpha=1.0 if focal else 0.85, zorder=3 if focal else 2,
                    marker="o" if focal else None, ms=2.2, label=label)
        ax.set_xlabel("Lead (h)")
        ax.set_ylabel(ylab)
        ax.set_xticks([0, 6, 12, 18, 24])
        ax.margins(x=0.03, y=0.06)
        n = d.get((var, "cam3", -1), (0, 0))[0]
        ax.set_title(f"{ylab.split(' RMSE')[0]} (n = {n/1000:.0f}k)", loc="left")
        ax.set_ylabel("RMSE (" + ylab.split("(")[1])
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], frameon=False, loc="outside lower center", ncol=4, fontsize=6.5)
    fig.suptitle(f"3 km vs 12 km, {a.cycles} cycles (26-30 Sep 2026), lower is better",
                 x=0.01, ha="left")
    fig.savefig(a.png, dpi=300)
    print("wrote", a.png)


if __name__ == "__main__":
    main()

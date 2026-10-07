"""
Figure: the 3 km CAM with the land surface against 12 km (CAM stage S6, test S6c).

    python tools/plot_cam_land.py leads.csv daynight.csv out.png --cycles 20

leads.csv: rows of tools/cam_compare.py's CSV with period = terrain = "all"
(by lead). daynight.csv: the pooled 2 m T rows for period day and night.
Panels: 2 m T, 10 m u, 10 m v RMSE by lead; pooled 2 m T bias by day and night.
"""
import argparse
import csv

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

ARMS = [("cam3L", "3 km NH + land surface"), ("nh12L2", "12 km NH + land surface"),
        ("am12", "12 km production (hydrostatic + land surface)"), ("cam3", "3 km NH, dry (S5g)")]
COLORS = {"cam3L": "#C0392B", "nh12L2": "#2E86C1", "am12": "#6C3483", "cam3": "#AAB7B8"}
PANELS = [("TMP", "2 m temperature", "K"), ("UGRD", "10 m u", "m/s"), ("VGRD", "10 m v", "m/s")]


def load(path):
    out = {}
    with open(path, newline="") as f:
        for r in csv.DictReader(f):
            if r["lead"].startswith("cycle"):
                continue
            out[(r["variable"], r["arm"], int(r["lead"]), r["period"])] = (int(r["n"]), float(r["bias"]), float(r["rmse"]))
    return out


def main():
    p = argparse.ArgumentParser()
    p.add_argument("leads"); p.add_argument("daynight"); p.add_argument("png")
    p.add_argument("--cycles", type=int, required=True)
    a = p.parse_args()
    d = load(a.leads); dn = load(a.daynight)
    fig, axes = plt.subplots(1, 4, figsize=(9.6, 2.9), layout="constrained",
                             gridspec_kw={"width_ratios": [1, 1, 1, 0.8]})
    for ax, (var, name, unit) in zip(axes[:3], PANELS):
        for arm, label in ARMS[::-1]:
            leads = sorted(L for (v, ar, L, pr) in d if v == var and ar == arm and L >= 0 and pr == "all")
            y = [d[(var, arm, L, "all")][2] for L in leads]
            focal = arm == "cam3L"
            ax.plot(leads, y, color=COLORS[arm], lw=1.8 if focal else 1.0, zorder=3 if focal else 2,
                    marker="o" if focal else None, ms=2.2, label=label,
                    ls="--" if arm == "cam3" else "-")
        n = d[(var, "cam3L", -1, "all")][0]
        ax.set_title(f"{name} (n = {n/1000:.0f}k)", loc="left")
        ax.set_xlabel("Lead (h)"); ax.set_ylabel(f"RMSE ({unit})")
        ax.set_xticks([0, 6, 12, 18, 24]); ax.margins(x=0.03, y=0.06)
    ax = axes[3]
    x = np.arange(len(ARMS)); w = 0.38
    for k, (per, hatch) in enumerate((("day", None), ("night", "///"))):
        vals = [dn[("TMP", arm, -1, per)][1] for arm, _ in ARMS]
        ax.bar(x + (k - 0.5) * w, vals, w, color=[COLORS[a_] for a_, _ in ARMS],
               edgecolor="black", linewidth=0.4, hatch=hatch, label=per)
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xticks(x); ax.set_xticklabels(["3 km\n+land", "12 km\n+land", "12 km\nprod.", "3 km\ndry"])
    ax.set_ylabel("2 m T bias (K)"); ax.set_title("2 m T bias, day | night (hatched)", loc="left")
    ax.margins(y=0.08)
    h, l = axes[0].get_legend_handles_labels()
    fig.legend(h[::-1], l[::-1], frameon=False, loc="outside lower center", ncol=4, fontsize=6.5)
    fig.suptitle(f"3 km CAM with the land surface vs 12 km, {a.cycles} cycles (26-30 Sep 2026), paired matches; "
                 "RMSE lower is better", x=0.01, ha="left")
    fig.savefig(a.png, dpi=300)
    print("wrote", a.png)


if __name__ == "__main__":
    main()

"""
Visualise the terrain instability: where it starts and how it spreads.

Produces a vertical cross-section through the terrain at successive forecast
hours, for flat ground and for a mountain, so the difference is visible rather
than inferred from survival counts.

    python visualize_instability.py

Colour choices follow the data's job: sigma_dot is SIGNED (up/down), so it
gets a diverging map with a neutral midpoint; wind speed is a MAGNITUDE, so it
gets a single-hue sequential ramp. No rainbow anywhere -- a rainbow map invents
boundaries where the field is smooth, which is exactly the error we are trying
to see through.
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.colors import TwoSlopeNorm

from grid import CGrid
from sigma import SigmaLevels, pressure_gradient_force, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma

INK = "#1a1a1a"
MUTED = "#6b6b6b"
GRID = "#e4e4e4"
TERRAIN = "#8a8175"


def build(terrain_m, n=36, dx=25e3):
    gr = CGrid(n, n, dx, dx, f0=9.81e-5, beta=1.69e-11, edge_mode="replicate")
    lev = SigmaLevels(20)
    h = terrain_m * np.exp(-(((gr.Xc - gr.Lx / 2) / 250e3) ** 2 +
                             ((gr.Yc - gr.Ly / 2) / 250e3) ** 2))
    m = PrimitiveSigma(gr, lev, terrain=h)
    ky = 2 * np.pi / gr.Ly
    p_s = 101325.0 * np.exp(-G0 * h / (RD * 280.0))
    m.pi = p_s - lev.p_top
    p = lev.pressure(m.pi)
    m.theta = (288.0 - 55.0 * (1 - p / p.max())
               - 3.0 * np.cos(ky * gr.Yc)) / (p / P0) ** KAPPA
    phi = m.geopotential()
    fx, fy = pressure_gradient_force(phi, m.theta, m.pi, lev, gr)
    m.u = -fy / gr.f0
    m.v = fx / gr.f0
    return m, h


def collect(terrain_m, hours=4):
    m, h = build(terrain_m)
    snaps = [(0.0, m.sigma_dot().copy(), m.u.copy())]
    series = [(0.0, float(np.abs(m.sigma_dot()).max()))]
    for hr in range(1, hours + 1):
        m.run(3600)
        ok = np.isfinite(m.u).all()
        sd = m.sigma_dot() if ok else np.full_like(m.u, np.nan)
        snaps.append((float(hr), sd.copy(), m.u.copy()))
        series.append((float(hr), float(np.abs(sd).max()) if ok else np.nan))
        if not ok:
            break
    return m, h, snaps, series


def main():
    hours = 4
    flat_m, flat_h, flat_s, flat_t = collect(0.0, hours)
    mtn_m, mtn_h, mtn_s, mtn_t = collect(2500.0, hours)

    gr, lev = mtn_m.grid, mtn_m.lev
    j = gr.ny // 2                      # cross-section row through the peak
    x_km = gr.xc / 1000.0
    sigma = lev.sigma

    ncol = min(len(flat_s), len(mtn_s), hours + 1)
    fig, axes = plt.subplots(2, ncol, figsize=(3.1 * ncol, 6.4), sharey=True,
                             constrained_layout=True)
    if ncol == 1:
        axes = axes.reshape(2, 1)

    # One symmetric scale across every panel, so panels are comparable.
    vals = [np.abs(s[1][:, j, :]) for s in (flat_s + mtn_s)
            if np.isfinite(s[1]).any()]
    vmax = float(np.nanpercentile(np.concatenate([v.ravel() for v in vals]), 99))
    vmax = max(vmax, 1e-9)
    norm = TwoSlopeNorm(vmin=-vmax, vcenter=0.0, vmax=vmax)

    for row, (snaps, terr, label) in enumerate(
            [(flat_s, flat_h, "flat ground"),
             (mtn_s, mtn_h, "2500 m mountain")]):
        for col in range(ncol):
            ax = axes[row, col]
            t, sd, u = snaps[col] if col < len(snaps) else (np.nan, None, None)

            if sd is None or not np.isfinite(sd).any():
                ax.text(0.5, 0.5, "diverged", ha="center", va="center",
                        transform=ax.transAxes, color="#b3261e", fontsize=11)
                ax.set_facecolor("#faf8f8")
            else:
                # sigma_dot is on half levels; plot the layer-centre average.
                field = 0.5 * (sd[:-1, j, :] + sd[1:, j, :])
                ax.pcolormesh(x_km, sigma, field, cmap="RdBu_r", norm=norm,
                              shading="auto", rasterized=True)

            # terrain silhouette, drawn in sigma space (ground is sigma = 1)
            ax.fill_between(x_km, 1.0, 1.0 - terr[j, :] / max(terr.max(), 1) * 0.0,
                            color="none")
            ax2 = ax.twinx()
            ax2.fill_between(x_km, 0, terr[j, :], color=TERRAIN, alpha=0.55, lw=0)
            ax2.set_ylim(0, max(3000.0, terr.max() * 1.05))
            ax2.set_yticks([])

            ax.set_ylim(1.0, 0.0)            # ground at the bottom
            ax.set_xlim(x_km[0], x_km[-1])
            for s in ax.spines.values():
                s.set_color(GRID)
            ax.tick_params(colors=MUTED, labelsize=8)
            if row == 0:
                ax.set_title(f"+{int(t)} h" if np.isfinite(t) else "",
                             color=INK, fontsize=10, pad=6)
            if col == 0:
                ax.set_ylabel(f"{label}\nσ  (1 = ground)",
                              color=INK, fontsize=9)
            if row == 1:
                ax.set_xlabel("x  (km)", color=MUTED, fontsize=9)

    sm = plt.cm.ScalarMappable(cmap="RdBu_r", norm=norm)
    cb = fig.colorbar(sm, ax=axes, location="bottom", fraction=0.045, pad=0.02,
                      shrink=0.55)
    cb.set_label("vertical velocity  σ̇   (1/s)   —  "
                 "red = sinking, blue = rising",
                 color=MUTED, fontsize=9)
    cb.ax.tick_params(colors=MUTED, labelsize=8)
    cb.outline.set_edgecolor(GRID)

    fig.suptitle("Vertical velocity through the terrain, hour by hour",
                 color=INK, fontsize=13, y=1.04)
    fig.savefig("instability_cross_section.png", dpi=150,
                bbox_inches="tight", facecolor="white")
    print("wrote instability_cross_section.png")

    # --- growth curve ------------------------------------------------------
    fig2, ax = plt.subplots(figsize=(7.2, 3.8), constrained_layout=True)
    for (series, color, label) in [(flat_t, "#2563a8", "flat ground"),
                                   (mtn_t, "#b3261e", "2500 m mountain")]:
        t = [s[0] for s in series]
        y = [s[1] for s in series]
        ax.plot(t, y, lw=2, color=color, marker="o", ms=5, label=label)
        if not np.isfinite(y[-1]):
            ax.plot(t[-2], y[-2], marker="X", ms=11, color=color, mew=0)
            ax.annotate("diverged", (t[-2], y[-2]), textcoords="offset points",
                        xytext=(8, 6), color=color, fontsize=9)

    ax.set_yscale("log")
    ax.set_xlabel("forecast hour", color=MUTED, fontsize=10)
    ax.set_ylabel("max |σ̇|   (1/s)", color=MUTED, fontsize=10)
    ax.set_title("Spurious vertical motion grows only over terrain",
                 color=INK, fontsize=12, loc="left")
    ax.grid(True, color=GRID, lw=0.8)
    ax.set_axisbelow(True)
    for s in ax.spines.values():
        s.set_color(GRID)
    ax.tick_params(colors=MUTED, labelsize=9)
    leg = ax.legend(frameon=False, fontsize=9, labelcolor=INK)
    fig2.savefig("instability_growth.png", dpi=150, bbox_inches="tight",
                 facecolor="white")
    print("wrote instability_growth.png")


if __name__ == "__main__":
    main()

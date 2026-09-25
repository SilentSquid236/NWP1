"""
Forecast product maps, drawn with matplotlib alone (no cartopy).

One PNG per product per hour, 1000 x 760 px, on a Lambert conformal map of
the domain with state, province, coast and lake lines. Each product keeps a
FIXED colour scale at every hour and in every run, so two maps are
comparable by eye, as on Pivotal-style sites.

Products (key: title):

  Surface     mslp   MSLP and 1000-500 hPa thickness
              tlow   temperature, lowest model level
              wlow   wind, lowest model level
  Upper air   t850   850 hPa temperature, heights and wind
              t700   700 hPa temperature, heights and wind
              z500   500 hPa heights and absolute vorticity
              w250   250 hPa wind speed and heights
  Analysis    an_sfc   2 m analysis temperature with the station reports
              an_mslp  analysis MSLP with the station reports
              an_850, an_500  analysis with the soundings at that level
  Verification  err_t, err_w  forecast minus observed, by station

Hour 0 of the forecast products is the analysis itself (a forecast's
snapshots start at its first output hour). The model is dry, so there is
no precipitation, reflectivity, dewpoint or CAPE product.
"""

import math
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.colors import BoundaryNorm, ListedColormap  # noqa: E402

from maps.geography import LambertConformal         # noqa: E402

MS_TO_KT = 1.0 / 0.514444
FIGSIZE, DPI = (10.0, 7.6), 100
AX_W, AX_H = 0.98, 0.80          # map frame as a fraction of the figure
MODEL_LABEL = "NWP1 12 km dry sigma model | analysis from observations only"

PRODUCTS = {
    # key: (group, menu label)
    "mslp": ("Surface", "MSLP & 1000-500 thickness"),
    "tlow": ("Surface", "Temperature (lowest level)"),
    "wlow": ("Surface", "Wind (lowest level)"),
    "t850": ("Upper air", "850 hPa temp, height, wind"),
    "t700": ("Upper air", "700 hPa temp, height, wind"),
    "z500": ("Upper air", "500 hPa height & vorticity"),
    "w250": ("Upper air", "250 hPa jet & height"),
    "an_sfc": ("Analysis", "2 m temperature & reports"),
    "an_mslp": ("Analysis", "MSLP & reports"),
    "an_850": ("Analysis", "850 hPa & soundings"),
    "an_500": ("Analysis", "500 hPa & soundings"),
    "err_t": ("Verification", "Temperature error by station"),
    "err_w": ("Verification", "Wind speed error by station"),
}


def temperature_cmap():
    """Blue-to-red temperature palette in 2-unit steps (turbo, dimmed ends)."""
    try:
        base = plt.get_cmap("turbo")               # matplotlib >= 3.3
    except ValueError:
        base = plt.get_cmap("jet")
    return ListedColormap(base(np.linspace(0.04, 0.96, 256)))


# ---------------------------------------------------------------------------
# The map frame
# ---------------------------------------------------------------------------

class Mapper:
    """Projected grid, boundaries and the frame every product is drawn in."""

    def __init__(self, lat, lon, boundaries, projection=None):
        self.P = projection or LambertConformal()
        self.lat, self.lon = lat, lon
        self.X, self.Y = self.P(lon, lat)
        self.theta = self.P.n * np.radians(lon - self.P.lon0)   # grid-north rotation
        dlat = (lat[1, 0] - lat[0, 0]) / 2
        dlon = (lon[0, 1] - lon[0, 0]) / 2
        # the grid's outline, edge to edge
        s = np.linspace(0, 1, 60)
        la0, la1 = lat[0, 0] - dlat, lat[-1, 0] + dlat
        lo0, lo1 = lon[0, 0] - dlon, lon[0, -1] + dlon
        blat = np.r_[np.full(60, la0), la0 + s * (la1 - la0), np.full(60, la1), la1 - s * (la1 - la0)]
        blon = np.r_[lo0 + s * (lo1 - lo0), np.full(60, lo1), lo1 - s * (lo1 - lo0), np.full(60, lo0)]
        self.bx, self.by = self.P(blon, blat)
        pad = 0.01 * (self.bx.max() - self.bx.min())
        self.extent = (self.bx.min() - pad, self.bx.max() + pad,
                       self.by.min() - pad, self.by.max() + pad)
        self.lines = {k: self.P(x, y) for k, (x, y) in boundaries.items()}
        # figure height chosen so the map fills the frame's width
        aspect = (self.extent[1] - self.extent[0]) / (self.extent[3] - self.extent[2])
        self.figsize = (FIGSIZE[0], min(11.0, FIGSIZE[0] * AX_W / aspect / AX_H))

    def rotate(self, u, v, theta=None):
        """East/north wind components -> map x/y components."""
        th = self.theta if theta is None else theta
        c, s = np.cos(th), np.sin(th)
        return u * c - v * s, u * s + v * c

    def frame(self, title, subtitle, left_note="", right_note=""):
        fig = plt.figure(figsize=self.figsize, dpi=DPI)
        H = self.figsize[1]
        # fixed-size header (0.75 in) and footer (1.0 in) whatever the height
        self.fy = lambda inches_from_top: 1.0 - inches_from_top / H
        ax = fig.add_axes([0.01, 1.0 / H, AX_W, 1.0 - 1.75 / H])
        ax.set_xlim(self.extent[:2]); ax.set_ylim(self.extent[2:])
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.8)
        fig.text(0.01, self.fy(0.27), title, fontsize=13, fontweight="bold", va="center")
        fig.text(0.01, self.fy(0.55), MODEL_LABEL, fontsize=9, va="center", color="0.25")
        fig.text(0.99, self.fy(0.27), subtitle[0], fontsize=11, ha="right", va="center")
        fig.text(0.99, self.fy(0.55), subtitle[1], fontsize=9, ha="right", va="center", color="0.25")
        if left_note:
            fig.text(0.01, 0.012, left_note, fontsize=8, color="0.3", va="bottom")
        if right_note:
            fig.text(0.99, 0.012, right_note, fontsize=8, color="0.3", ha="right", va="bottom")
        return fig, ax

    def terrain_base(self, ax, terrain):
        ax.contourf(self.X, self.Y, terrain, levels=[0, 150, 300, 500, 750, 1000, 1500],
                    colors=[str(g) for g in (0.97, 0.93, 0.89, 0.85, 0.81, 0.77)],
                    extend="max", zorder=0)

    def boundaries(self, ax, dark=False):
        c = "k" if dark else "0.15"
        style = {"coast": (0.8, c), "lakes": (0.7, c), "countries": (1.0, c), "states": (0.55, "0.3")}
        for k, (x, y) in self.lines.items():
            lw, col = style.get(k, (0.5, "0.3"))
            ax.plot(x, y, color=col, lw=lw, zorder=5)
        ax.plot(self.bx, self.by, color="0.35", lw=0.8, ls=(0, (4, 3)), zorder=5)
        ax.set_xlim(self.extent[:2]); ax.set_ylim(self.extent[2:])

    def colorbar(self, fig, mappable, label, ticks=None):
        H = fig.get_size_inches()[1]
        cax = fig.add_axes([0.12, 0.52 / H, 0.76, 0.17 / H])
        cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", ticks=ticks)
        cb.set_label(label, fontsize=9)
        cb.ax.tick_params(labelsize=8)
        return cb

    def barbs(self, ax, u, v, every=7, color="k", length=4.9):
        sl = (slice(every // 2, None, every), slice(every // 2, None, every))
        ux, vy = self.rotate(u[sl], v[sl], self.theta[sl])
        ax.barbs(self.X[sl], self.Y[sl], ux * MS_TO_KT, vy * MS_TO_KT,
                 length=length, linewidth=0.55, color=color, zorder=6)

    def contour(self, ax, F, levels, color="k", lw=0.9, fmt="%d", ls="solid",
                label=True, fs=7):
        cs = ax.contour(self.X, self.Y, F, levels=levels, colors=color,
                        linewidths=lw, linestyles=ls, zorder=4)
        if label and len(cs.levels):
            ax.clabel(cs, fmt=fmt, fontsize=fs, inline=True, inline_spacing=2)
        return cs

    def highs_lows(self, ax, F, size=9, edge=3, min_sep=8):
        """H and L at local extrema over a (2*size+1)^2 window, away from the edge."""
        Fp = np.pad(F, size, mode="edge")
        mx = np.full(F.shape, -np.inf); mn = np.full(F.shape, np.inf)
        ny, nx = F.shape
        for dj in range(-size, size + 1, 3):
            for di in range(-size, size + 1, 3):
                w = Fp[size + dj: size + dj + ny, size + di: size + di + nx]
                mx = np.maximum(mx, w); mn = np.minimum(mn, w)
        inner = np.zeros(F.shape, bool); inner[edge:-edge, edge:-edge] = True
        done = []
        for mask, letter, col in (((F >= mx) & inner, "H", "#1f3fbf"),
                                  ((F <= mn) & inner, "L", "#c01818")):
            for j, i in zip(*np.nonzero(mask)):
                if any(abs(j - a) + abs(i - b) < min_sep for a, b in done):
                    continue
                done.append((j, i))
                ax.text(self.X[j, i], self.Y[j, i], letter, color=col, fontsize=18,
                        fontweight="bold", ha="center", va="center", zorder=7)
                ax.text(self.X[j, i], self.Y[j, i], f"\n\n{F[j, i]:.0f}", color=col,
                        fontsize=7.5, fontweight="bold", ha="center", va="center", zorder=7)

    def points(self, lon, lat):
        return self.P(np.asarray(lon, float), np.asarray(lat, float))


# ---------------------------------------------------------------------------
# Titles
# ---------------------------------------------------------------------------

def time_labels(run_time, hour):
    valid = run_time + timedelta(hours=float(hour))
    return (f"Init {run_time:%Y-%m-%d %H}Z   F{int(round(hour)):03d}",
            f"Valid {valid:%a %Y-%m-%d %H:%M}Z")


def mask_below(F, below):
    return np.where(below, np.nan, F)


# ---------------------------------------------------------------------------
# Forecast products
# ---------------------------------------------------------------------------

def product_mslp(m, d, sub, note):
    fig, ax = m.frame("MSLP (hPa), 1000-500 hPa thickness (dam)", sub, *note)
    m.terrain_base(ax, d["terrain"])
    th = d["thick_1000_500"] / 10.0
    lv = np.arange(480, 606, 6)
    m.contour(ax, th, [v for v in lv if v < 540], color="#2050d0", lw=0.9, ls="dashed")
    m.contour(ax, th, [540], color="#1030a0", lw=1.8, ls="dashed")
    m.contour(ax, th, [v for v in lv if v > 540], color="#d03020", lw=0.9, ls="dashed")
    p = d["mslp"] / 100.0
    m.contour(ax, p, np.arange(940, 1064, 2), color="k", lw=0.9, fmt="%d")
    m.highs_lows(ax, p)
    m.boundaries(ax)
    return fig


def product_tlow(m, d, sub, note):
    agl = float(np.nanmean(d["z_low_agl"]))
    fig, ax = m.frame(f"Temperature (\u00b0F), lowest model level (\u2248{agl:.0f} m above ground)", sub, *note)
    tf = (d["T_low"] - 273.15) * 9 / 5 + 32
    lv = np.arange(-20, 112, 2)
    cf = ax.contourf(m.X, m.Y, tf, levels=lv, cmap=temperature_cmap(), extend="both", zorder=1)
    m.contour(ax, tf, np.arange(-20, 111, 10), color="0.2", lw=0.5, fs=6.5)
    m.contour(ax, tf, [32], color="#0b2a8a", lw=1.6, label=False)
    m.boundaries(ax)
    m.colorbar(fig, cf, "\u00b0F", ticks=np.arange(-20, 111, 10))
    return fig


def product_wlow(m, d, sub, note):
    agl = float(np.nanmean(d["z_low_agl"]))
    fig, ax = m.frame(f"Wind speed (kt) and barbs, lowest model level (\u2248{agl:.0f} m above ground)", sub, *note)
    spd = np.hypot(d["u_low"], d["v_low"]) * MS_TO_KT
    lv = np.arange(0, 72, 4)
    cf = ax.contourf(m.X, m.Y, spd, levels=lv, cmap="YlGnBu", extend="max", zorder=1)
    m.boundaries(ax)
    m.barbs(ax, d["u_low"], d["v_low"])
    m.colorbar(fig, cf, "kt", ticks=lv[::2])
    return fig


def _level_temp(m, d, sub, note, L, tlv, hstep):
    fig, ax = m.frame(f"{L} hPa temperature (\u00b0C), height (dam), wind (kt)", sub, *note)
    below = d[f"below{L}"]
    tc = mask_below(d[f"T{L}"] - 273.15, below)
    cf = ax.contourf(m.X, m.Y, tc, levels=tlv, cmap=temperature_cmap(), extend="both", zorder=1)
    m.contour(ax, tc, [0], color="#0b2a8a", lw=1.6, label=False)
    z = mask_below(d[f"Z{L}"] / 10.0, below)
    zr = np.nanmin(z), np.nanmax(z)
    if np.isfinite(zr[0]):
        m.contour(ax, z, np.arange(math.floor(zr[0] / hstep) * hstep, zr[1] + hstep, hstep),
                  color="k", lw=1.0)
    m.boundaries(ax)
    m.barbs(ax, mask_below(d[f"u{L}"], below), mask_below(d[f"v{L}"], below))
    if below.any():
        ax.contourf(m.X, m.Y, below.astype(float), levels=[0.5, 1.5], colors="none",
                    hatches=["////"], zorder=3)
    m.colorbar(fig, cf, "\u00b0C" + ("   (hatched: below ground)" if below.any() else ""),
               ticks=tlv[::4])
    return fig


def product_t850(m, d, sub, note):
    return _level_temp(m, d, sub, note, 850, np.arange(-30, 32, 1), 3)


def product_t700(m, d, sub, note):
    return _level_temp(m, d, sub, note, 700, np.arange(-40, 22, 1), 3)


def product_z500(m, d, sub, note):
    fig, ax = m.frame("500 hPa height (dam) and absolute vorticity (10\u207b\u2075 s\u207b\u00b9)", sub, *note)
    av = d["avort500"] * 1e5
    lv = np.array([12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 34, 38, 42, 48])
    cmap = plt.get_cmap("YlOrRd")
    cf = ax.contourf(m.X, m.Y, av, levels=lv, cmap=cmap, extend="max", zorder=1)
    z = d["Z500"] / 10.0
    m.contour(ax, z, np.arange(math.floor(np.nanmin(z) / 6) * 6, np.nanmax(z) + 6, 6), color="k", lw=1.1)
    m.boundaries(ax)
    m.barbs(ax, d["u500"], d["v500"], color="0.25")
    m.colorbar(fig, cf, "10\u207b\u2075 s\u207b\u00b9", ticks=lv)
    return fig


def product_w250(m, d, sub, note):
    fig, ax = m.frame("250 hPa wind speed (kt) and height (dam)", sub, *note)
    spd = np.hypot(d["u250"], d["v250"]) * MS_TO_KT
    lv = np.arange(50, 190, 10)
    cf = ax.contourf(m.X, m.Y, spd, levels=lv, cmap="RdPu", extend="max", zorder=1)
    z = d["Z250"] / 10.0
    m.contour(ax, z, np.arange(math.floor(np.nanmin(z) / 12) * 12, np.nanmax(z) + 12, 12), color="k", lw=1.0)
    m.boundaries(ax)
    m.barbs(ax, d["u250"], d["v250"], color="0.2")
    m.colorbar(fig, cf, "kt", ticks=lv[::2])
    return fig


FORECAST_PRODUCTS = {"mslp": product_mslp, "tlow": product_tlow, "wlow": product_wlow,
                     "t850": product_t850, "t700": product_t700, "z500": product_z500,
                     "w250": product_w250}


# ---------------------------------------------------------------------------
# Analysis with observations
# ---------------------------------------------------------------------------

def _station_plot(m, ax, obs, value_key, fmt, withheld, color="k"):
    """Value text + barb at each station. obs: dicts with lat, lon, value(s), u, v."""
    for o in obs:
        x, y = m.points(o["lon"], o["lat"])
        held = o["station"] in withheld
        c = "#7a1fa2" if held else color
        if o.get(value_key) is not None:
            ax.text(x, y, fmt(o[value_key]), fontsize=6.5, color=c, ha="right", va="bottom",
                    zorder=8, fontweight="bold" if not held else "normal", clip_on=True)
        if o.get("u") is not None and o.get("v") is not None:
            th = m.P.n * np.radians(o["lon"] - m.P.lon0)
            ux, vy = m.rotate(np.array([o["u"]]), np.array([o["v"]]), th)
            ax.barbs([x], [y], ux * MS_TO_KT, vy * MS_TO_KT, length=4.6, linewidth=0.5,
                     color=c, zorder=8, clip_on=True)
        ax.plot([x], [y], marker="o", ms=2.2, color=c, mfc="none" if held else c, zorder=8,
                clip_on=True)


def product_an_sfc(m, d, sub, note, obs, withheld):
    fig, ax = m.frame("Analysis 2 m temperature (\u00b0F) with station reports (\u00b0F, kt)", sub, *note)
    tf = (d["T_2m"] - 273.15) * 9 / 5 + 32
    lv = np.arange(-20, 112, 2)
    cf = ax.contourf(m.X, m.Y, tf, levels=lv, cmap=temperature_cmap(), extend="both", zorder=1, alpha=0.85)
    m.boundaries(ax)
    _station_plot(m, ax, obs, "tf", lambda v: f"{v:.0f}", withheld)
    m.colorbar(fig, cf, "\u00b0F   (purple, open: withheld from the analysis, used only to score it)",
               ticks=np.arange(-20, 111, 10))
    return fig


def product_an_mslp(m, d, sub, note, obs, withheld):
    fig, ax = m.frame("Analysis MSLP (hPa) with station reports (hPa)", sub, *note)
    m.terrain_base(ax, d["terrain"])
    p = d["mslp"] / 100.0
    m.contour(ax, p, np.arange(940, 1064, 2), color="k", lw=0.9)
    m.highs_lows(ax, p)
    m.boundaries(ax)
    _station_plot(m, ax, obs, "pmsl", lambda v: f"{v:.1f}", withheld, color="#1a4a9a")
    return fig


def _an_level(m, d, sub, note, L, sondes, tlv, hstep):
    fig, ax = m.frame(f"Analysis {L} hPa temperature (\u00b0C) and height (dam) with soundings", sub, *note)
    below = d[f"below{L}"]
    tc = mask_below(d[f"T{L}"] - 273.15, below)
    cf = ax.contourf(m.X, m.Y, tc, levels=tlv, cmap=temperature_cmap(), extend="both", zorder=1, alpha=0.85)
    z = mask_below(d[f"Z{L}"] / 10.0, below)
    m.contour(ax, z, np.arange(math.floor(np.nanmin(z) / hstep) * hstep, np.nanmax(z) + hstep, hstep),
              color="k", lw=1.0)
    m.boundaries(ax)
    _station_plot(m, ax, sondes, "tc", lambda v: f"{v:.0f}", set(), color="#8a1010")
    m.colorbar(fig, cf, f"\u00b0C   ({len(sondes)} soundings reported at {L} hPa)", ticks=tlv[::4])
    return fig


# ---------------------------------------------------------------------------
# Verification
# ---------------------------------------------------------------------------

def product_error(m, d, sub, note, pairs, kind):
    """pairs: dicts with lat, lon, err (forecast minus observed)."""
    if kind == "t":
        title, unit, lim = "Forecast minus observed temperature (\u00b0C), by station", "\u00b0C", 8
    else:
        title, unit, lim = "Forecast minus observed wind speed (kt), by station", "kt", 20
    fig, ax = m.frame(title, sub, *note)
    m.terrain_base(ax, d["terrain"])
    m.boundaries(ax)
    e = np.array([p["err"] for p in pairs], dtype=float)
    if len(e):
        x, y = m.points([p["lon"] for p in pairs], [p["lat"] for p in pairs])
        sc = ax.scatter(x, y, c=np.clip(e, -lim, lim), cmap="RdBu_r", vmin=-lim, vmax=lim,
                        s=26, edgecolors="0.2", linewidths=0.4, zorder=8)
        m.colorbar(fig, sc, f"{unit}   (red: forecast too high / too strong)",
                   ticks=np.linspace(-lim, lim, 9))
        ax.text(0.01, 0.985, f"n = {len(e)}   bias {e.mean():+.2f} {unit}   "
                f"RMSE {np.sqrt(np.mean(e ** 2)):.2f} {unit}",
                transform=ax.transAxes, fontsize=10, va="top", zorder=9,
                bbox=dict(fc="white", ec="0.6", lw=0.5, pad=3))
    return fig

"""
Forecast product maps, drawn with matplotlib alone (no cartopy).

One PNG per product per hour on a Lambert conformal map. The frame is the
largest rectangle inside the model domain, so every pixel of the map holds
model data. The grid is drawn past the frame and clipped. Each product keeps
a FIXED colour scale at every hour and in every run, so maps compare by eye.
Names, units (mb, kt, degF at the surface, degC aloft) and layout follow
Pivotal Weather's model pages.

Hour 0 of the forecast products is the analysis (a forecast's snapshots
start at its first output hour). The model is dry, so there is no
precipitation, reflectivity, dewpoint or CAPE product, and "surface" fields
are the LOWEST MODEL LEVEL, as each title says.
"""

import math
from datetime import datetime, timedelta

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt                     # noqa: E402
from matplotlib.colors import ListedColormap        # noqa: E402

from maps.geography import LambertConformal         # noqa: E402

MS_TO_KT = 1.0 / 0.514444
WIDTH_IN, DPI = 10.0, 100
HEADER_IN, FOOTER_IN = 0.62, 1.05
MODEL_LABEL = "NWP1 12 km (dry sigma model, observation-only analysis)"

# key: (group, menu label, hover fields). Hover fields name entries of the
# per-hour hover data written by make_maps.hover_fields().
PRODUCTS = {
    "mslp": ("Surface", "MSLP & 1000-500 mb thickness", ["mslp", "thick"]),
    "tlow": ("Surface", "Temperature (lowest level)", ["t_low_f"]),
    "tchg": ("Surface", "1-hr temperature change", ["t_low_f", "dt1_f"]),
    "wlow": ("Surface", "Wind (lowest level)", ["wind_low"]),
    "t925": ("Upper air", "925 mb temp, height, wind", ["t925", "z925", "wind925"]),
    "t850": ("Upper air", "850 mb temp, height, wind", ["t850", "z850", "wind850"]),
    "t700": ("Upper air", "700 mb temp, height, wind", ["t700", "z700", "wind700"]),
    "w700": ("Upper air", "700 mb vertical velocity", ["omega700", "z700"]),
    "z500": ("Upper air", "500 mb height, wind, vorticity", ["z500", "avort500", "wind500"]),
    "w300": ("Upper air", "300 mb height & wind", ["z300", "wind300"]),
    "w250": ("Upper air", "250 mb height & wind", ["z250", "wind250"]),
    "an_sfc": ("Analysis", "2 m temperature & station reports", ["t2m_f"]),
    "an_mslp": ("Analysis", "MSLP & station reports", ["mslp"]),
    "an_850": ("Analysis", "850 mb & soundings", ["t850", "z850", "wind850"]),
    "an_500": ("Analysis", "500 mb & soundings", ["t500", "z500", "wind500"]),
    "err_t": ("Verification", "Temperature error by station", []),
    "err_w": ("Verification", "Wind speed error by station", []),
}


def temperature_cmap():
    try:
        base = plt.get_cmap("turbo")               # matplotlib >= 3.3
    except ValueError:
        base = plt.get_cmap("jet")
    return ListedColormap(base(np.linspace(0.04, 0.96, 256)))


# ---------------------------------------------------------------------------
# Time labels (UTC and US Eastern, as Pivotal prints them)
# ---------------------------------------------------------------------------

def _nth_sunday(year, month, n):
    d = datetime(year, month, 1)
    d += timedelta(days=(6 - d.weekday()) % 7)
    return d + timedelta(weeks=n - 1)


def eastern(utc):
    """(local datetime, 'EDT'|'EST') for a UTC datetime, by the US DST rule."""
    y = utc.year
    start = _nth_sunday(y, 3, 2) + timedelta(hours=7)     # 2 am EST = 07 UTC
    end = _nth_sunday(y, 11, 1) + timedelta(hours=6)      # 2 am EDT = 06 UTC
    if start <= utc < end:
        return utc - timedelta(hours=4), "EDT"
    return utc - timedelta(hours=5), "EST"


def _hm(t):
    h = t.hour % 12 or 12
    return f"{h} {'AM' if t.hour < 12 else 'PM'}"


def time_labels(run_time, hour):
    """(left line, right line) of the header for a forecast hour."""
    valid = run_time + timedelta(hours=float(hour))
    loc, zone = eastern(valid)
    left = (f"Init: {run_time:%H}z {run_time:%b} {run_time.day} {run_time:%Y}"
            f"   Forecast hour: {int(round(hour))}")
    right = (f"Valid: {valid:%H}z {valid:%a} {valid:%b} {valid.day} {valid:%Y}"
             f"   ({_hm(loc)} {zone} {loc:%a})")
    return left, right


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
        la0, la1 = lat[0, 0] - dlat, lat[-1, 0] + dlat
        lo0, lo1 = lon[0, 0] - dlon, lon[0, -1] + dlon
        s = np.linspace(0, 1, 200)
        sx, sy = self.P(lo0 + s * (lo1 - lo0), np.full_like(s, la0))
        nx_, ny_ = self.P(lo0 + s * (lo1 - lo0), np.full_like(s, la1))
        wx, wy = self.P(np.full_like(s, lo0), la0 + s * (la1 - la0))
        ex, ey = self.P(np.full_like(s, lo1), la0 + s * (la1 - la0))
        # The largest axis-aligned rectangle inside the projected domain.
        # Shrink it by half a cell so the contour fill reaches every edge.
        half = 0.5 * 12_000.0
        self.extent = (wx.max() + half, ex.min() - half, sy.max() + half, ny_.min() - half)
        self.lines = {k: self.P(x, y) for k, (x, y) in boundaries.items()}
        aspect = (self.extent[1] - self.extent[0]) / (self.extent[3] - self.extent[2])
        self.map_w_in = WIDTH_IN * 0.98
        self.map_h_in = self.map_w_in / aspect
        self.figsize = (WIDTH_IN, self.map_h_in + HEADER_IN + FOOTER_IN)

    # geometry the viewer needs to turn a mouse position into a grid point
    def pixel_frame(self):
        W, H = self.figsize[0] * DPI, self.figsize[1] * DPI
        left = 0.01 * W
        top = HEADER_IN * DPI
        return {"img_w": int(W), "img_h": int(H), "ax_left": left, "ax_top": top,
                "ax_w": self.map_w_in * DPI, "ax_h": self.map_h_in * DPI,
                "extent": list(map(float, self.extent)),
                "proj": {"n": self.P.n, "F": self.P.F, "rho0": self.P.rho0,
                         "R": self.P.R_EARTH, "lon0": self.P.lon0}}

    def rotate(self, u, v, theta=None):
        """East/north wind components -> map x/y components."""
        th = self.theta if theta is None else theta
        c, s = np.cos(th), np.sin(th)
        return u * c - v * s, u * s + v * c

    def frame(self, title, sub, left_note="", right_note=""):
        fig = plt.figure(figsize=self.figsize, dpi=DPI)
        H = self.figsize[1]
        ax = fig.add_axes([0.01, FOOTER_IN / H, 0.98, self.map_h_in / H])
        ax.set_xlim(self.extent[:2]); ax.set_ylim(self.extent[2:])
        ax.set_aspect("equal"); ax.set_xticks([]); ax.set_yticks([])
        for sp in ax.spines.values():
            sp.set_linewidth(0.8)
        fy = lambda d: 1.0 - d / H
        fig.text(0.01, fy(0.20), title, fontsize=13, fontweight="bold", va="center")
        fig.text(0.01, fy(0.47), sub[0], fontsize=10, va="center")
        fig.text(0.99, fy(0.47), sub[1], fontsize=10, ha="right", va="center")
        fig.text(0.01, 0.06 / H, MODEL_LABEL + (f"   |   {left_note}" if left_note else ""),
                 fontsize=8, color="0.3", va="bottom")
        if right_note:
            fig.text(0.99, 0.06 / H, right_note, fontsize=8, color="0.3", ha="right", va="bottom")
        return fig, ax

    def terrain_base(self, ax, terrain):
        ax.contourf(self.X, self.Y, terrain, levels=[-1, 150, 300, 500, 750, 1000, 1500],
                    colors=[str(g) for g in (0.97, 0.93, 0.89, 0.85, 0.81, 0.77)],
                    extend="max", zorder=0)

    def boundaries(self, ax):
        style = {"coast": (0.8, "k"), "lakes": (0.7, "k"), "countries": (1.0, "k"),
                 "states": (0.55, "0.25")}
        for k, (x, y) in self.lines.items():
            lw, col = style.get(k, (0.5, "0.3"))
            ax.plot(x, y, color=col, lw=lw, zorder=5)
        ax.set_xlim(self.extent[:2]); ax.set_ylim(self.extent[2:])

    def colorbar(self, fig, mappable, label, ticks=None):
        H = fig.get_size_inches()[1]
        cax = fig.add_axes([0.12, 0.72 / H, 0.76, 0.16 / H])
        cb = fig.colorbar(mappable, cax=cax, orientation="horizontal", ticks=ticks)
        cb.set_label(label, fontsize=9)
        cb.ax.tick_params(labelsize=8)
        return cb

    def barbs(self, ax, u, v, every=7, color="k", length=4.9):
        sl = (slice(every // 2, None, every), slice(every // 2, None, every))
        ux, vy = self.rotate(u[sl], v[sl], self.theta[sl])
        ax.barbs(self.X[sl], self.Y[sl], ux * MS_TO_KT, vy * MS_TO_KT,
                 length=length, linewidth=0.55, color=color, zorder=6, clip_on=True)

    def contour(self, ax, F, levels, color="k", lw=0.9, fmt="%d", ls="solid",
                label=True, fs=7):
        levels = [l for l in levels if np.nanmin(F) <= l <= np.nanmax(F)] \
            if np.isfinite(np.nanmin(F)) else []
        if not levels:
            return None
        cs = ax.contour(self.X, self.Y, F, levels=levels, colors=color,
                        linewidths=lw, linestyles=ls, zorder=4)
        if label:
            ax.clabel(cs, fmt=fmt, fontsize=fs, inline=True, inline_spacing=2)
        return cs

    def inside(self, x, y):
        e = self.extent
        return (x >= e[0]) & (x <= e[1]) & (y >= e[2]) & (y <= e[3])

    def highs_lows(self, ax, F, size=9, min_sep=8):
        """H and L at local extrema over a (2*size+1)^2 window, inside the frame."""
        Fp = np.pad(F, size, mode="edge")
        mx = np.full(F.shape, -np.inf); mn = np.full(F.shape, np.inf)
        ny, nx = F.shape
        for dj in range(-size, size + 1, 3):
            for di in range(-size, size + 1, 3):
                w = Fp[size + dj: size + dj + ny, size + di: size + di + nx]
                mx = np.maximum(mx, w); mn = np.minimum(mn, w)
        ok = self.inside(self.X, self.Y)
        ok[:3, :] = ok[-3:, :] = False; ok[:, :3] = ok[:, -3:] = False
        done = []
        for mask, letter, col in (((F >= mx) & ok, "H", "#1f3fbf"),
                                  ((F <= mn) & ok, "L", "#c01818")):
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

    def thin(self, obs, min_km=45.0, keep=None):
        """Greedy thinning of station plots to at least min_km apart, inside the frame."""
        out, xy = [], []
        order = sorted(obs, key=lambda o: (o["station"] in (keep or set()),
                                           -sum(v is not None for v in o.values())))
        for o in order:
            x, y = self.points(o["lon"], o["lat"])
            if not self.inside(x, y):
                continue
            if all((x - a) ** 2 + (y - b) ** 2 >= (min_km * 1e3) ** 2 for a, b in xy):
                out.append(o); xy.append((float(x), float(y)))
        return out


def mask_below(F, below):
    return np.where(below, np.nan, F)


def _height_levels(z, step):
    lo, hi = np.nanmin(z), np.nanmax(z)
    if not np.isfinite(lo):
        return []
    return list(np.arange(math.floor(lo / step) * step, hi + step, step))


# ---------------------------------------------------------------------------
# Surface products
# ---------------------------------------------------------------------------

def product_mslp(m, d, sub, note):
    fig, ax = m.frame("MSLP (mb) and 1000-500 mb thickness (dam)", sub, *note)
    m.terrain_base(ax, d["terrain"])
    th = d["thick_1000_500"] / 10.0
    lv = np.arange(480, 606, 6)
    m.contour(ax, th, [v for v in lv if v < 540], color="#2050d0", lw=0.9, ls="dashed")
    m.contour(ax, th, [540], color="#1030a0", lw=1.8, ls="dashed")
    m.contour(ax, th, [v for v in lv if v > 540], color="#d03020", lw=0.9, ls="dashed")
    p = d["mslp"] / 100.0
    m.contour(ax, p, np.arange(940, 1064, 2), color="k", lw=0.9)
    m.highs_lows(ax, p)
    m.boundaries(ax)
    return fig


def _t_f(K):
    return (K - 273.15) * 9 / 5 + 32


def product_tlow(m, d, sub, note):
    agl = float(np.nanmean(d["z_low_agl"]))
    fig, ax = m.frame(f"Temperature (\u00b0F), lowest model level (\u2248{agl:.0f} m AGL)", sub, *note)
    tf = _t_f(d["T_low"])
    cf = ax.contourf(m.X, m.Y, tf, levels=np.arange(-20, 112, 2), cmap=temperature_cmap(),
                     extend="both", zorder=1)
    m.contour(ax, tf, np.arange(-20, 111, 10), color="0.2", lw=0.5, fs=6.5)
    m.contour(ax, tf, [32], color="#0b2a8a", lw=1.6, label=False)
    m.boundaries(ax)
    m.colorbar(fig, cf, "\u00b0F", ticks=np.arange(-20, 111, 10))
    return fig


def product_tchg(m, d, sub, note):
    fig, ax = m.frame("1-hr temperature change (\u00b0F), lowest model level", sub, *note)
    dt = d["dT1"] * 9 / 5
    lv = np.arange(-8, 8.25, 0.5)
    cf = ax.contourf(m.X, m.Y, dt, levels=lv, cmap="RdBu_r", extend="both", zorder=1)
    m.contour(ax, dt, [-4, -2, 2, 4], color="0.25", lw=0.5, fs=6.5)
    m.boundaries(ax)
    m.colorbar(fig, cf, "\u00b0F change over the previous hour"
               + ("   (F001: from the analysis to the first model hour)" if round(d.get("hour", -1)) == 1 else ""),
               ticks=np.arange(-8, 9, 2))
    return fig


def product_wlow(m, d, sub, note):
    agl = float(np.nanmean(d["z_low_agl"]))
    fig, ax = m.frame(f"Wind speed (kt) and barbs, lowest model level (\u2248{agl:.0f} m AGL)", sub, *note)
    spd = np.hypot(d["u_low"], d["v_low"]) * MS_TO_KT
    lv = np.arange(0, 72, 4)
    cf = ax.contourf(m.X, m.Y, spd, levels=lv, cmap="YlGnBu", extend="max", zorder=1)
    m.boundaries(ax)
    m.barbs(ax, d["u_low"], d["v_low"])
    m.colorbar(fig, cf, "kt", ticks=lv[::2])
    return fig


# ---------------------------------------------------------------------------
# Upper air
# ---------------------------------------------------------------------------

def _level_temp(m, d, sub, note, L, tlv, hstep):
    fig, ax = m.frame(f"{L} mb temperature (\u00b0C), height (dam), wind (kt)", sub, *note)
    below = d[f"below{L}"]
    tc = mask_below(d[f"T{L}"] - 273.15, below)
    cf = ax.contourf(m.X, m.Y, tc, levels=tlv, cmap=temperature_cmap(), extend="both", zorder=1)
    m.contour(ax, tc, [0], color="#0b2a8a", lw=1.6, label=False)
    z = mask_below(d[f"Z{L}"] / 10.0, below)
    m.contour(ax, z, _height_levels(z, hstep), color="k", lw=1.0)
    m.boundaries(ax)
    m.barbs(ax, mask_below(d[f"u{L}"], below), mask_below(d[f"v{L}"], below))
    if below.any():
        ax.contourf(m.X, m.Y, below.astype(float), levels=[0.5, 1.5], colors="none",
                    hatches=["////"], zorder=3)
    m.colorbar(fig, cf, "\u00b0C" + ("   (hatched: below ground)" if below.any() else ""),
               ticks=tlv[::4])
    return fig


def product_t925(m, d, sub, note):
    return _level_temp(m, d, sub, note, 925, np.arange(-30, 36, 1), 3)


def product_t850(m, d, sub, note):
    return _level_temp(m, d, sub, note, 850, np.arange(-30, 32, 1), 3)


def product_t700(m, d, sub, note):
    return _level_temp(m, d, sub, note, 700, np.arange(-40, 22, 1), 3)


def product_w700(m, d, sub, note):
    fig, ax = m.frame("700 mb vertical velocity (-\u00b5b/s, positive = rising) and height (dam)", sub, *note)
    w = -d["omega700"] * 10.0                         # Pa/s -> -microbar/s
    lv = np.array([-30, -20, -15, -10, -6, -3, -1, 1, 3, 6, 10, 15, 20, 30])
    cf = ax.contourf(m.X, m.Y, w, levels=lv, cmap="PuOr", extend="both", zorder=1)
    z = d["Z700"] / 10.0
    m.contour(ax, z, _height_levels(z, 3), color="k", lw=1.0)
    m.boundaries(ax)
    m.colorbar(fig, cf, "-\u00b5b/s  (kinematic, from the model's divergence)", ticks=lv)
    return fig


def product_z500(m, d, sub, note):
    fig, ax = m.frame("500 mb height (dam), wind (kt), absolute vorticity (10\u207b\u2075 s\u207b\u00b9)", sub, *note)
    av = d["avort500"] * 1e5
    lv = np.array([12, 14, 16, 18, 20, 22, 24, 26, 28, 30, 34, 38, 42, 48])
    cf = ax.contourf(m.X, m.Y, av, levels=lv, cmap="YlOrRd", extend="max", zorder=1)
    z = d["Z500"] / 10.0
    m.contour(ax, z, _height_levels(z, 6), color="k", lw=1.1)
    m.boundaries(ax)
    m.barbs(ax, d["u500"], d["v500"], color="0.25")
    m.colorbar(fig, cf, "10\u207b\u2075 s\u207b\u00b9", ticks=lv)
    return fig


def _jet(m, d, sub, note, L):
    fig, ax = m.frame(f"{L} mb height (dam) and wind (kt)", sub, *note)
    spd = np.hypot(d[f"u{L}"], d[f"v{L}"]) * MS_TO_KT
    lv = np.arange(50, 190, 10)
    cf = ax.contourf(m.X, m.Y, spd, levels=lv, cmap="RdPu", extend="max", zorder=1)
    z = d[f"Z{L}"] / 10.0
    m.contour(ax, z, _height_levels(z, 12), color="k", lw=1.0)
    m.boundaries(ax)
    m.barbs(ax, d[f"u{L}"], d[f"v{L}"], color="0.2")
    m.colorbar(fig, cf, "kt", ticks=lv[::2])
    return fig


def product_w300(m, d, sub, note):
    return _jet(m, d, sub, note, 300)


def product_w250(m, d, sub, note):
    return _jet(m, d, sub, note, 250)


FORECAST_PRODUCTS = {"mslp": product_mslp, "tlow": product_tlow, "tchg": product_tchg,
                     "wlow": product_wlow, "t925": product_t925, "t850": product_t850,
                     "t700": product_t700, "w700": product_w700, "z500": product_z500,
                     "w300": product_w300, "w250": product_w250}


# ---------------------------------------------------------------------------
# Analysis with observations
# ---------------------------------------------------------------------------

def _station_plot(m, ax, obs, value_key, fmt, withheld, color="k", min_km=45.0):
    shown = m.thin(obs, min_km=min_km, keep=withheld)
    for o in shown:
        x, y = m.points(o["lon"], o["lat"])
        held = o["station"] in withheld
        c = "#7a1fa2" if held else color
        if o.get(value_key) is not None:
            ax.text(x, y, fmt(o[value_key]), fontsize=7, color=c, ha="right", va="bottom",
                    zorder=8, fontweight="bold" if not held else "normal", clip_on=True)
        if o.get("u") is not None and o.get("v") is not None:
            th = m.P.n * np.radians(o["lon"] - m.P.lon0)
            ux, vy = m.rotate(np.array([o["u"]]), np.array([o["v"]]), th)
            ax.barbs([x], [y], ux * MS_TO_KT, vy * MS_TO_KT, length=4.6, linewidth=0.5,
                     color=c, zorder=8, clip_on=True)
        ax.plot([x], [y], marker="o", ms=2.4, color=c, mfc="none" if held else c, zorder=8,
                clip_on=True)
    return len(shown)


def product_an_sfc(m, d, sub, note, obs, withheld):
    fig, ax = m.frame("Analysis 2 m temperature (\u00b0F) with station reports (\u00b0F, kt)", sub, *note)
    tf = _t_f(d["T_2m"])
    cf = ax.contourf(m.X, m.Y, tf, levels=np.arange(-20, 112, 2), cmap=temperature_cmap(),
                     extend="both", zorder=1, alpha=0.85)
    m.boundaries(ax)
    n = _station_plot(m, ax, obs, "tf", lambda v: f"{v:.0f}", withheld)
    m.colorbar(fig, cf, f"\u00b0F   ({n} of {len(obs)} stations shown; purple open = withheld)",
               ticks=np.arange(-20, 111, 10))
    return fig


def product_an_mslp(m, d, sub, note, obs, withheld):
    fig, ax = m.frame("Analysis MSLP (mb) with station reports (mb)", sub, *note)
    m.terrain_base(ax, d["terrain"])
    p = d["mslp"] / 100.0
    m.contour(ax, p, np.arange(940, 1064, 2), color="k", lw=0.9)
    m.highs_lows(ax, p)
    m.boundaries(ax)
    n = _station_plot(m, ax, [o for o in obs if o.get("pmsl") is not None], "pmsl",
                      lambda v: f"{v:.1f}", withheld, color="#1a4a9a", min_km=55.0)
    ax.text(0.01, 0.01, f"{n} stations plotted (thinned to 55 km)", transform=ax.transAxes,
            fontsize=8, color="0.3", zorder=9)
    return fig


def _an_level(m, d, sub, note, L, sondes, tlv, hstep):
    fig, ax = m.frame(f"Analysis {L} mb temperature (\u00b0C) and height (dam) with soundings", sub, *note)
    below = d[f"below{L}"]
    tc = mask_below(d[f"T{L}"] - 273.15, below)
    cf = ax.contourf(m.X, m.Y, tc, levels=tlv, cmap=temperature_cmap(), extend="both",
                     zorder=1, alpha=0.85)
    z = mask_below(d[f"Z{L}"] / 10.0, below)
    m.contour(ax, z, _height_levels(z, hstep), color="k", lw=1.0)
    m.boundaries(ax)
    n = _station_plot(m, ax, sondes, "tc", lambda v: f"{v:.0f}", set(), color="#8a1010", min_km=0.0)
    m.colorbar(fig, cf, f"\u00b0C   ({n} soundings in the frame at {L} mb)", ticks=tlv[::4])
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
                        s=28, edgecolors="0.2", linewidths=0.4, zorder=8, clip_on=True)
        m.colorbar(fig, sc, f"{unit}   (red: forecast too high / too strong)",
                   ticks=np.linspace(-lim, lim, 9))
        ax.text(0.01, 0.985, f"n = {len(e)}   bias {e.mean():+.2f} {unit}   "
                f"RMSE {np.sqrt(np.mean(e ** 2)):.2f} {unit}",
                transform=ax.transAxes, fontsize=10, va="top", zorder=9,
                bbox=dict(fc="white", ec="0.6", lw=0.5, pad=3))
    return fig

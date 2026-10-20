#!/usr/bin/env python3
"""Score NWP1's verified forecasts against persistence (stdlib only, Python >= 3.6).

    python3 tools/scoreboard.py [DATA_ROOT] [--days 7]

Reads the verification archives that `src/verify_pending.py` fills from live
ASOS observations:

    data/verification/<cycle>/matches.jsonl       12 km production runs
    data/verification_3km/<cycle>/matches.jsonl   3 km runs

and compares each forecast with two persistence forecasts built from the SAME
observations, so every comparison is paired (same station, same valid time):

    P0   station persistence: the station's own observation at the run's start
         time, held for the whole forecast. Very hard to beat in the first hours.
    P24  24 h persistence: the station's observation 24 h before the valid time
         ("tomorrow at this hour is like today"). Keeps the day-night cycle.

Writes
    data/status/scoreboard.md        the status page's skill section
    data/verification_scores.csv     one row per grid, cycle, variable, lead band

Skill against a reference is 1 - MSE_model / MSE_reference (Murphy 1988): above 0 the model
beats it, 0 ties, below 0 the model is worse. "best persistence" is the smaller
of the P0 and P24 errors in that lead band, a deliberately strict reference.
"""
import csv
import datetime as dt
import glob
import json
import math
import os
import sys
from collections import defaultdict

VARS = {  # archive variable -> (label, scale, unit)
    "TMP": ("2 m temperature", 1.0, "K"),
    "DPT": ("2 m dewpoint", 1.0, "K"),
    "WIND": ("10 m wind (vector)", 1.0, "m/s"),
    "PMSL": ("sea-level pressure", 0.01, "hPa"),
}
BANDS = [(1, 6), (7, 12), (13, 18), (19, 24)]
GRIDS = [("12 km", "verification"), ("3 km", "verification_3km")]
LAG_ALERT_H = 48.0   # newest verified run start; normal lag is ~26-35 h


def hour_key(t):
    """Round a datetime to the nearest whole hour."""
    t = t.replace(second=0, microsecond=0)
    return (t + dt.timedelta(minutes=30)).replace(minute=0)


def parse_time(s):
    s = s.replace("Z", "")
    if "+" in s:
        s = s.split("+")[0]
    s = s.split(".")[0]
    for f in ("%Y-%m-%dT%H:%M:%S", "%Y-%m-%dT%H:%M", "%Y-%m-%d %H:%M:%S"):
        try:
            return dt.datetime.strptime(s, f)
        except ValueError:
            pass
    raise ValueError(f"unparsed time {s!r}")


def band_of(lead):
    for lo, hi in BANDS:
        if lo <= lead <= hi:
            return (lo, hi)
    return None


def read_archive(path_glob):
    """Yield (cycle, var, station, valid_hour, lead_h, forecast, observation)."""
    for f in sorted(glob.glob(path_glob)):
        cyc = os.path.basename(os.path.dirname(f))
        try:
            dt.datetime.strptime(cyc, "%Y%m%d_%H")
        except ValueError:
            continue
        with open(f) as fh:
            for line in fh:
                try:
                    m = json.loads(line)
                except ValueError:
                    continue
                if m.get("source") != "asos":
                    continue
                v = m.get("variable")
                if v not in ("TMP", "DPT", "UGRD", "VGRD", "PMSL"):
                    continue
                e = m.get("error")
                if e is None or not math.isfinite(e):
                    continue
                yield (cyc, v, m["station"], hour_key(parse_time(m["valid_time"])),
                       float(m["lead_hours"]), float(m["forecast"]), float(m["observation"]))


def build(rows):
    """Return per-(cycle, var, band) sums of squared error for model, P0, P24.

    Wind is scored as a vector: a (station, valid hour) counts only when both
    components were matched, for the model and for each persistence reference.
    """
    rows = list(rows)
    obs = {}
    for cyc, v, st, vh, lead, fc, ob in rows:
        obs[(v, st, vh)] = ob
    # (cycle, var, band) -> [n, se_model, sum_err_model, se_p0, se_p24]
    acc = defaultdict(lambda: [0, 0.0, 0.0, 0.0, 0.0])
    wind = defaultdict(dict)    # (cyc, st, vh, lead) -> {"UGRD": (fc, ob), ...}
    for cyc, v, st, vh, lead, fc, ob in rows:
        b = band_of(int(round(lead)))
        if b is None:
            continue
        t0 = dt.datetime.strptime(cyc, "%Y%m%d_%H")
        if v in ("UGRD", "VGRD"):
            wind[(cyc, st, vh, lead)][v] = (fc, ob)
            continue
        p0 = obs.get((v, st, t0))
        p24 = obs.get((v, st, vh - dt.timedelta(hours=24)))
        if p0 is None or p24 is None:
            continue
        a = acc[(cyc, v, b)]
        a[0] += 1
        a[1] += (fc - ob) ** 2
        a[2] += fc - ob
        a[3] += (p0 - ob) ** 2
        a[4] += (p24 - ob) ** 2
    for (cyc, st, vh, lead), d in wind.items():
        if len(d) < 2:
            continue
        b = band_of(int(round(lead)))
        t0 = dt.datetime.strptime(cyc, "%Y%m%d_%H")
        p0u, p0v = obs.get(("UGRD", st, t0)), obs.get(("VGRD", st, t0))
        k24 = vh - dt.timedelta(hours=24)
        p24u, p24v = obs.get(("UGRD", st, k24)), obs.get(("VGRD", st, k24))
        if None in (p0u, p0v, p24u, p24v):
            continue
        (fu, ou), (fv, ov) = d["UGRD"], d["VGRD"]
        a = acc[(cyc, "WIND", b)]
        a[0] += 1
        a[1] += (fu - ou) ** 2 + (fv - ov) ** 2
        a[2] += math.hypot(fu, fv) - math.hypot(ou, ov)    # speed bias
        a[3] += (p0u - ou) ** 2 + (p0v - ov) ** 2
        a[4] += (p24u - ou) ** 2 + (p24v - ov) ** 2
    return acc


def pooled(acc, cycles, var, band=None):
    n = se = s = s0 = s24 = 0.0
    for (cyc, v, b), a in acc.items():
        if v != var or cyc not in cycles or (band is not None and b != band):
            continue
        n += a[0]; se += a[1]; s += a[2]; s0 += a[3]; s24 += a[4]
    if n == 0:
        return None
    return {"n": int(n), "rmse": math.sqrt(se / n), "bias": s / n,
            "p0": math.sqrt(s0 / n), "p24": math.sqrt(s24 / n),
            "skill_best": 1.0 - se / min(s0, s24) if min(s0, s24) > 0 else float("nan"),
            "skill_p24": 1.0 - se / s24 if s24 > 0 else float("nan")}


def fmt(x, d=2):
    return "—" if x is None or not math.isfinite(x) else f"{x:.{d}f}"


def main(argv):
    data = os.environ.get("NWP_DATA_ROOT", "/data5/pierce/AINWP/data")
    days = 7.0
    args = list(argv)
    if "--days" in args:
        i = args.index("--days"); days = float(args[i + 1]); del args[i:i + 2]
    if args:
        data = args[0]
    now = dt.datetime.utcnow()
    md = ["## Forecast skill against live observations\n",
          f"Scored from ASOS observations by `tools/scoreboard.py`, {now:%Y-%m-%d %H:%M}Z. "
          "Each forecast is compared with two persistence forecasts built from the same "
          "observations: **P0** holds each station's value at the start time; **P24** uses "
          "the value 24 h before the valid time. Skill = 1 − MSE(model) / MSE(best of P0, P24): "
          "above 0 the model beats both.\n"]
    csv_rows = []
    alerts = []
    for grid, sub in GRIDS:
        acc = build(read_archive(os.path.join(data, sub, "2*_??", "matches.jsonl")))
        cycles = sorted({k[0] for k in acc})
        if not cycles:
            md.append(f"### {grid}\n\n_No verified cycles yet._\n")
            continue
        last = dt.datetime.strptime(cycles[-1], "%Y%m%d_%H")
        lag = (now - last).total_seconds() / 3600.0
        recent = {c for c in cycles
                  if (last - dt.datetime.strptime(c, "%Y%m%d_%H")).total_seconds() <= days * 86400}
        if lag > LAG_ALERT_H:
            alerts.append(f"{grid}: newest verified run is {lag:.0f} h old")
        md.append(f"### {grid} — {len(recent)} cycles to {last:%m-%d %H}Z "
                  f"(newest verified run started {lag:.0f} h ago)\n")
        md.append("| variable | lead | n | model RMSE | bias | P0 | P24 | skill vs best persistence |")
        md.append("|---|---|---|---|---|---|---|---|")
        for var, (label, sc, unit) in VARS.items():
            for b in BANDS:
                r = pooled(acc, recent, var, b)
                if r is None:
                    continue
                md.append(f"| {label} ({unit}) | {b[0]}–{b[1]} h | {r['n']} | {fmt(r['rmse'] * sc)} | "
                          f"{fmt(r['bias'] * sc, 2)} | {fmt(r['p0'] * sc)} | {fmt(r['p24'] * sc)} | "
                          f"{fmt(r['skill_best'])} |")
        md.append("")
        # Per-cycle temperature line: does each run beat 24 h persistence?
        md.append(f"2 m temperature by cycle, 1–24 h ({grid}): model RMSE / P24 RMSE (K)\n")
        cells = []
        for c in sorted(recent)[-12:]:
            r = pooled(acc, {c}, "TMP")
            if r:
                mark = "better" if r["rmse"] < r["p24"] else "**worse**"
                cells.append(f"{c[4:6]}-{c[6:8]} {c[9:11]}Z {r['rmse']:.1f}/{r['p24']:.1f} {mark}")
        md.append("; ".join(cells) + "\n")
        for c in cycles:
            for var in VARS:
                for b in BANDS:
                    r = pooled(acc, {c}, var, b)
                    if r:
                        csv_rows.append([grid, c, var, f"{b[0]}-{b[1]}", r["n"], f"{r['rmse']:.4f}",
                                         f"{r['bias']:.4f}", f"{r['p0']:.4f}", f"{r['p24']:.4f}",
                                         f"{r['skill_best']:.4f}"])
    if alerts:
        md.insert(2, "**Verification is behind:** " + "; ".join(alerts) + "\n")
    os.makedirs(os.path.join(data, "status"), exist_ok=True)
    tmp = os.path.join(data, "status", "scoreboard.md.tmp")
    with open(tmp, "w") as fh:
        fh.write("\n".join(md) + "\n")
    os.replace(tmp, os.path.join(data, "status", "scoreboard.md"))
    tmp = os.path.join(data, "verification_scores.csv.tmp")
    with open(tmp, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["grid", "cycle", "variable", "lead_h", "n", "rmse", "bias",
                    "rmse_p0", "rmse_p24", "skill_vs_best_persistence"])
        w.writerows(csv_rows)
    os.replace(tmp, os.path.join(data, "verification_scores.csv"))
    print(f"scoreboard: {len(csv_rows)} rows; alerts: {len(alerts)}")
    return 0


if __name__ == "__main__":
    sys.exit(main(sys.argv[1:]))

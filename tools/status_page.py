#!/usr/bin/env python3
"""Build a one-page Markdown status report for NWP1 (stdlib only, Python >= 3.6).

Run by cron on the server (tools/status_push.sh, every 10 min); the page is
pushed to the 'status' branch of the GitHub repo so it can be read from a
phone. It reads only NWP1's own logs and progress files:

  data/status/notes.md                 notes appended by the AI collaborator
  data/verification_tests/*/progress.txt   experiments touched in the last 48 h
  data/logs/compute_hours.csv          every cam_cycle / daily run (wall, status)
  data/logs/cam_<stamp>.log            live step rate of running forecasts

Usage: python3 status_page.py [DATA_ROOT]  (writes Markdown to stdout)
"""
import csv
import datetime as dt
import glob
import os
import re
import subprocess
import sys

DATA = sys.argv[1] if len(sys.argv) > 1 else os.environ.get(
    "NWP_DATA_ROOT", "/data5/pierce/AINWP/data")
NOW = dt.datetime.utcnow()
CRON_HOURS = (1, 7, 13, 19)        # cam_cycle.sh at HH:45 UTC
ALERT_RE = re.compile(r"DIVERGED|STOP|TESTS FAILED|Traceback|rc=[1-9]|status=[1-9]|deadline at")


def acked():
    """Substrings of alerts already dealt with (data/status/ack.txt, one per line)."""
    try:
        return [l.strip() for l in open(os.path.join(DATA, "status", "ack.txt")) if l.strip()]
    except OSError:
        return []


def age(path):
    try:
        return (NOW - dt.datetime.utcfromtimestamp(os.path.getmtime(path))).total_seconds() / 3600
    except OSError:
        return None


def tail(path, n):
    try:
        with open(path, errors="replace") as f:
            lines = f.read().replace("\r", "\n").splitlines()
    except OSError:
        return []
    return [l.rstrip() for l in lines if l.strip()][-n:]


def next_cron():
    best = None
    for d in (0, 1):
        day = (NOW + dt.timedelta(days=d)).date()
        for h in CRON_HOURS:
            t = dt.datetime(day.year, day.month, day.day, h, 45)
            if t > NOW and (best is None or t < best):
                best = t
    return best


def fmt_age(h):
    if h is None:
        return "?"
    return f"{h * 60:.0f} min ago" if h < 1.5 else f"{h:.1f} h ago"


def section_scores(out, alerts):
    """Forecast skill against live observations, written by tools/scoreboard.py
    after each verification run (daily.sh verify)."""
    p = os.path.join(DATA, "status", "scoreboard.md")
    if not os.path.exists(p):
        alerts.append("no forecast scoreboard yet (tools/scoreboard.py has not run)")
        return
    a = age(p)
    if a is not None and a > 14:
        alerts.append(f"forecast scoreboard is {a:.0f} h old (verification runs every 6 h)")
    with open(p) as fh:
        txt = fh.read().rstrip()
    if "Verification is behind" in txt:
        alerts.append("verification is behind: see Forecast skill")
    out.append(txt + "\n")


def section_notes(out):
    p = os.path.join(DATA, "status", "notes.md")
    entries, cur = [], []
    try:
        for line in open(p, errors="replace").read().splitlines():
            if line.startswith("- ") and cur:
                entries.append("\n".join(cur)); cur = []
            if line.strip():
                cur.append(line)
        if cur:
            entries.append("\n".join(cur))
    except OSError:
        pass
    out.append("## Notes from the AI collaborator (newest first)\n")
    if not entries:
        out.append("_No notes yet._\n")
    for e in reversed(entries[-15:]):
        out.append(e)
    out.append("")


def running_forecasts():
    """[(stamp, last progress line)] for forecast.py processes now running."""
    try:
        ps = subprocess.run(["ps", "-eo", "etimes,args"], stdout=subprocess.PIPE,
                            universal_newlines=True, timeout=10).stdout
    except Exception:
        return []
    rows = []
    for line in ps.splitlines():
        if "forecast.py" not in line or "--run-dir" not in line:
            continue
        m = re.search(r"--run-dir\s+(\S+)", line)
        rd = m.group(1).rstrip("/") if m else "?"
        stamp = os.path.basename(rd).replace("obs_", "")
        res = "3 km" if "tensors_3km" in rd else "12 km"
        log = os.path.join(DATA, "logs", f"cam_{stamp}.log")
        steps = [l for l in tail(log, 400) if "steps/s" in l]
        prog = re.sub(r"\s+", " ", steps[-1].strip()) if steps else "(no step line yet)"
        rows.append((stamp, res, prog))
    return rows


def section_server(out, alerts):
    out.append("## Server\n")
    try:
        la = open("/proc/loadavg").read().split()[:3]
        out.append(f"- load average {' / '.join(la)} (104 logical CPUs)")
    except OSError:
        pass
    try:
        mi = dict(l.split(":", 1) for l in open("/proc/meminfo"))
        tot = int(mi["MemTotal"].split()[0]) / 1e6
        av = int(mi["MemAvailable"].split()[0]) / 1e6
        out.append(f"- memory available {av:.0f} of {tot:.0f} GB")
    except Exception:
        pass
    nc = next_cron()
    out.append(f"- next cron cycle {nc:%H:%M}Z (in {(nc - NOW).total_seconds() / 60:.0f} min)")
    rows = running_forecasts()
    if rows:
        out.append("\n| running forecast | grid | latest progress |\n|---|---|---|")
        for s, res, p in rows:
            out.append(f"| {s} | {res} | {p} |")
    else:
        out.append("- no forecast running")
    out.append("")


def section_experiments(out, alerts):
    out.append("## Experiments (active in the last 24 h)\n")
    found = False
    for p in sorted(glob.glob(os.path.join(DATA, "verification_tests", "*", "progress.txt")),
                    key=os.path.getmtime, reverse=True):
        a = age(p)
        if a is None or a > 24:
            continue
        found = True
        name = os.path.basename(os.path.dirname(p))
        lines = tail(p, 2000)
        done = sum(1 for l in lines if "(completed" in l)
        cut = sum(1 for l in lines if "deadline at" in l)
        fin = any(l.startswith("all done") for l in lines)
        res = os.path.join(os.path.dirname(p), f"{name}_result.txt")
        if fin or os.path.exists(res):
            state = "finished"
        elif any(l.startswith("stopped") for l in lines):
            state = "stopped"
        elif lines and lines[-1].startswith("yield") and a < 0.5:
            state = "waiting for a free slot (it steps aside for cron)"
        else:
            state = "running"
        out.append(f"### {name.upper()} — {state}, last update {fmt_age(a)}\n")
        out.append(f"- forecasts completed: **{done}**" + (f", cut at deadline: {cut}" if cut else ""))
        for l in lines:
            if ALERT_RE.search(l) and "yield" not in l:
                alerts.append(f"{name}: {l[:160]}")
        out.append("\n```\n" + "\n".join(l[:150] for l in lines[-8:] if not l.startswith("yield")) + "\n```\n")
    if not found:
        out.append("_None active._\n")


def section_cron(out, alerts):
    out.append("## Operational cron cycles (last 3 days)\n")
    p = os.path.join(DATA, "logs", "compute_hours.csv")
    rows = []
    try:
        for r in csv.DictReader(open(p)):
            try:
                cyc = dt.datetime.strptime(r["cycle"][:13], "%Y-%m-%dT%H")
            except ValueError:
                continue
            # campaign re-runs use old cycle dates; operational ones are recent
            if (NOW - cyc).total_seconds() < 3 * 86400:
                rows.append(r)
    except OSError:
        pass
    if not rows:
        out.append("_No rows._\n")
        return
    out.append("| cycle | started (UTC) | script | wall (min) | status |\n|---|---|---|---|---|")
    for r in rows[-12:]:
        st = r.get("status", "")
        flag = "" if st in ("0", "") else " ⚠"
        out.append(f"| {r['cycle']}Z | {r['start_utc'][5:16].replace('T', ' ')} | {r['script']} "
                   f"| {r['wall_min']} | {st}{flag} |")
        if flag and (NOW - dt.datetime.strptime(r["start_utc"][:16], "%Y-%m-%dT%H:%M")).total_seconds() < 86400:
            alerts.append(f"cron {r['cycle']} {r['script']} status {st}")
    out.append("")


def main():
    body, alerts = [], []
    scores = []; section_scores(scores, alerts)
    section_notes(body)
    exp = []; section_experiments(exp, alerts)
    cron = []; section_cron(cron, alerts)
    srv = []; section_server(srv, alerts)
    ack = acked()
    alerts = [a for a in alerts if not any(k in a for k in ack)]
    out = [f"# NWP1 status — {NOW:%Y-%m-%d %H:%M}Z\n",
           "Rebuilt every 10 min on the school server from NWP1's own logs. "
           "Times are UTC. If this time is more than ~20 min old, the server or its cron is not running.\n"]
    if alerts:
        out.append("## Attention\n")
        out += [f"- {a}" for a in alerts[-10:]]
        out.append("")
    out += scores + body + exp + cron + srv
    sys.stdout.write("\n".join(out) + "\n")


if __name__ == "__main__":
    main()

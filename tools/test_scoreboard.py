#!/usr/bin/env python3
"""Tests for tools/scoreboard.py on a synthetic archive with known answers.

    python3 tools/test_scoreboard.py
"""
import datetime as dt
import json
import math
import os
import sys
import tempfile

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import scoreboard as sb  # noqa: E402

PASS = []


def check(name, ok, detail=""):
    PASS.append(ok)
    print(f"{'PASS' if ok else 'FAIL'}  {name}  {detail}")


def write_archive(root, sub, cyc, matches):
    d = os.path.join(root, sub, cyc)
    os.makedirs(d, exist_ok=True)
    with open(os.path.join(d, "matches.jsonl"), "w") as fh:
        for m in matches:
            fh.write(json.dumps(m) + "\n")


def m(var, st, valid, lead, fc, ob, src="asos"):
    return {"source": src, "variable": var, "station": st, "valid_time": valid.isoformat(),
            "lead_hours": lead, "forecast": fc, "observation": ob, "error": fc - ob}


def main():
    t0 = dt.datetime(2026, 10, 11, 0)
    prev = t0 - dt.timedelta(hours=24)
    with tempfile.TemporaryDirectory() as root:
        # Previous cycle supplies observations from t0-24 h to t0 (lead 0..24).
        pm = [m("TMP", "AAA", prev + dt.timedelta(hours=h), h, 280.0, 280.0 + h * 0.1)
              for h in range(0, 24)]   # stops before t0: the t0 observation comes from the current cycle
        # Current cycle: obs = 283 + 0.5*h; model = obs + 1 (constant +1 K error).
        cm = [m("TMP", "AAA", t0 + dt.timedelta(hours=h), h, 284.0 + 0.5 * h, 283.0 + 0.5 * h)
              for h in range(0, 25)]
        # One record from another source must be ignored.
        cm.append(m("TMP", "AAA", t0 + dt.timedelta(hours=3), 3, 999.0, 283.0, src="raob"))
        # Wind: obs (3, 4) every hour; model (3, 4) + (0, 2); P0 = obs at t0 = (3, 4).
        for h in range(0, 25):
            v = t0 + dt.timedelta(hours=h)
            cm.append(m("UGRD", "AAA", v, h, 3.0, 3.0))
            cm.append(m("VGRD", "AAA", v, h, 6.0, 4.0))
        for h in range(0, 24):    # previous-cycle wind obs for P24
            v = prev + dt.timedelta(hours=h)
            pm.append(m("UGRD", "AAA", v, h, 0.0, 0.0))
            pm.append(m("VGRD", "AAA", v, h, 0.0, 0.0))
        write_archive(root, "verification", prev.strftime("%Y%m%d_%H"), pm)
        write_archive(root, "verification", t0.strftime("%Y%m%d_%H"), cm)

        acc = sb.build(sb.read_archive(os.path.join(root, "verification", "2*_??", "matches.jsonl")))
        cyc = t0.strftime("%Y%m%d_%H")
        r = sb.pooled(acc, {cyc}, "TMP", (1, 6))
        check("model RMSE equals the constant +1 K error", abs(r["rmse"] - 1.0) < 1e-9, f"{r['rmse']:.6f}")
        check("model bias is +1 K", abs(r["bias"] - 1.0) < 1e-9, f"{r['bias']:.6f}")
        check("n counts leads 1-6 only, other sources ignored", r["n"] == 6, f"n={r['n']}")
        # P0 = obs at t0 = 283; error at lead h = -0.5 h -> RMSE over h=1..6
        p0 = math.sqrt(sum((0.5 * h) ** 2 for h in range(1, 7)) / 6)
        check("P0 RMSE matches the analytic value", abs(r["p0"] - p0) < 1e-9, f"{r['p0']:.6f} vs {p0:.6f}")
        # P24 = obs 24 h earlier = 280 + 0.1 h; error = (280 + 0.1h) - (283 + 0.5h) = -3 - 0.4h
        p24 = math.sqrt(sum((3 + 0.4 * h) ** 2 for h in range(1, 7)) / 6)
        check("P24 RMSE matches the analytic value", abs(r["p24"] - p24) < 1e-9, f"{r['p24']:.6f} vs {p24:.6f}")
        sk = 1 - 1.0 / min(p0, p24) ** 2
        check("skill uses the better persistence", abs(r["skill_best"] - sk) < 1e-9, f"{r['skill_best']:.6f}")
        w = sb.pooled(acc, {cyc}, "WIND", (7, 12))
        check("wind is scored as a vector (error 2 m/s)", abs(w["rmse"] - 2.0) < 1e-9, f"{w['rmse']:.6f}")
        check("wind P0 is exact for steady wind", abs(w["p0"]) < 1e-12, f"{w['p0']:.6f}")
        check("wind P24 uses the earlier calm (5 m/s error)", abs(w["p24"] - 5.0) < 1e-9, f"{w['p24']:.6f}")
        rc = sb.main([root])
        md = open(os.path.join(root, "status", "scoreboard.md")).read()
        check("main writes the status section", rc == 0 and "Forecast skill against live observations" in md)
        check("main writes the CSV", os.path.exists(os.path.join(root, "verification_scores.csv")))
    print(f"\n{sum(PASS)}/{len(PASS)} passed")
    return 0 if all(PASS) else 1


if __name__ == "__main__":
    sys.exit(main())

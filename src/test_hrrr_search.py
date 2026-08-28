"""
Regression test for the HRRR variable search regex.

This is the bug that cost a live run: HRRR index entries begin with a colon
(":TMP:850 mb:anl"), so a regex anchored with "^" matched nothing. Herbie then
downloaded no file, and the failure surfaced far away as a FileNotFoundError
from cfgrib about a path that was never created.

The index lines below are copied verbatim from a real HRRR prs file, so this
test would have caught it without any network access.

Run:  python src/test_hrrr_search.py
"""

import re
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ingest_hrrr import HRRR_SEARCH, CF_ALIASES, resolve_variable

# Verbatim from hrrr.t13z.wrfprsf00.grib2.idx
REAL_INDEX = [
    ":HGT:50 mb:anl", ":TMP:50 mb:anl", ":RH:50 mb:anl", ":DPT:50 mb:anl",
    ":SPFH:50 mb:anl", ":VVEL:50 mb:anl", ":UGRD:50 mb:anl", ":VGRD:50 mb:anl",
    ":HGT:850 mb:anl", ":TMP:850 mb:anl", ":RH:850 mb:anl",
    ":UGRD:850 mb:anl", ":VGRD:850 mb:anl",
    ":HGT:1000 mb:anl", ":TMP:1000 mb:anl",
    # Non-isobaric entries that must NOT match
    ":TMP:2 m above ground:anl", ":UGRD:10 m above ground:anl",
    ":PRES:surface:anl", ":REFC:entire atmosphere:anl",
    ":TMP:surface:anl", ":APCP:surface:0-0 day acc fcst",
    # Forecast-hour form, which must still match
    ":TMP:500 mb:1 hour fcst", ":UGRD:500 mb:1 hour fcst",
]

results = []


def report(name, ok, detail):
    results.append(ok)
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def test_matches_isobaric_entries():
    rx = re.compile(HRRR_SEARCH)
    matched = [e for e in REAL_INDEX if rx.search(e)]
    wanted = [e for e in REAL_INDEX
              if re.search(r":(TMP|RH|UGRD|VGRD|HGT):\d+ mb:", e)]
    ok = len(matched) == len(wanted) and len(matched) > 0
    report("matches the five variables on pressure levels", ok,
           f"{len(matched)} of {len(REAL_INDEX)} index entries matched "
           f"(expected {len(wanted)})")


def test_rejects_non_isobaric():
    rx = re.compile(HRRR_SEARCH)
    bad = [e for e in REAL_INDEX
           if rx.search(e) and " mb:" not in e]
    ok = not bad
    report("does not match surface or 2 m entries", ok,
           f"unwanted matches: {bad if bad else 'none'}")


def test_not_anchored_at_start():
    """
    THE regression. HRRR entries start with ':', so '^' before the variable
    name matches nothing. Assert the pattern has no such anchor.
    """
    anchored = HRRR_SEARCH.startswith("^")
    rx = re.compile(HRRR_SEARCH)
    hits = sum(1 for e in REAL_INDEX if rx.search(e))

    broken = re.compile(r"^(?:TMP|RH|UGRD|VGRD|HGT):\d+ mb:")
    broken_hits = sum(1 for e in REAL_INDEX if broken.search(e))

    ok = not anchored and hits > 0 and broken_hits == 0
    report("pattern is not anchored with ^ (the live-run bug)", ok,
           f"current pattern {HRRR_SEARCH!r} -> {hits} hits; "
           f"the old ^-anchored pattern -> {broken_hits} hits")


def test_matches_forecast_hours_too():
    """Analysis rows end ':anl', forecast rows ':N hour fcst'. Both must match."""
    rx = re.compile(HRRR_SEARCH)
    anl = [e for e in REAL_INDEX if rx.search(e) and e.endswith(":anl")]
    fcst = [e for e in REAL_INDEX if rx.search(e) and "fcst" in e]
    ok = len(anl) > 0 and len(fcst) == 2
    report("matches both analysis and forecast-hour entries", ok,
           f"{len(anl)} analysis, {len(fcst)} forecast rows matched")


def test_cf_short_names_resolve():
    """
    cfgrib renames GRIB variables to CF short names, so what we ASK for is not
    what comes back: TMP->t, RH->r, UGRD->u, VGRD->v, HGT->gh. These are the
    exact names a real HRRR prs download produced.
    """
    class FakeDS:
        def __init__(self, names): self.data_vars = list(names)
        def __contains__(self, k): return k in self.data_vars

    ds = FakeDS(["t", "u", "v", "gh", "r"])       # verbatim from a live run
    got = {c: resolve_variable(ds, c)
           for c in ("TMP", "RH", "UGRD", "VGRD", "HGT")}
    want = {"TMP": "t", "RH": "r", "UGRD": "u", "VGRD": "v", "HGT": "gh"}

    ok = got == want
    report("cfgrib CF short names resolve to our channels", ok,
           f"{got}")


def test_missing_variable_reports_usefully():
    """The error must name what was tried and what was available."""
    class FakeDS:
        def __init__(self, names): self.data_vars = list(names)
        def __contains__(self, k): return k in self.data_vars

    try:
        resolve_variable(FakeDS(["t", "u"]), "HGT")
        ok, msg = False, "no error raised"
    except KeyError as e:
        msg = str(e)
        ok = "gh" in msg and "['t', 'u']" in msg
    report("missing variable error names candidates and contents", ok,
           msg[:110])


if __name__ == "__main__":
    print("\nHRRR search regex and variable naming\n" + "=" * 62)
    for fn in (test_matches_isobaric_entries, test_rejects_non_isobaric,
               test_not_anchored_at_start, test_matches_forecast_hours_too,
               test_cf_short_names_resolve,
               test_missing_variable_reports_usefully):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

#!/usr/bin/env python3
"""
Token ledger: what this project has cost, and what it would cost via the API.

    python tools/tokens.py --report
    python tools/tokens.py --add 2026-09-16 --in 120000 --cache-read 8400000 \
                           --cache-write 950000 --out 145000 --note "P-50 probes"
    python tools/tokens.py --prices
    python tools/tokens.py --estimate

WHY THIS IS TWO NUMBERS AND NOT ONE

The project is being built through a Claude subscription, so the MARGINAL cost
of a token is zero -- the cost is the monthly plan, whatever is used. That is
the number that leaves the bank account.

It is not the number a paper wants. "What does it cost to build a regional NWP
model with an AI collaborator" is asked in units someone else can reproduce,
and the reproducible unit is API pricing. So the ledger reports:

    BILLED      what was actually paid (the subscription, prorated)
    API-EQUIV   what the same tokens would have cost at list API prices

API-EQUIV is a shadow price. It is the honest figure for the write-up, and it
must never be described as money spent.

WHAT IS MEASURED AND WHAT IS NOT

Each row is tagged `measured` or `estimated`, and the report keeps the two
totals apart. A session's real numbers come from the usage breakdown at the end
of that session; there is no way to recover them afterwards, which makes this
the same shape of problem as the verification archive (P-07): record it on the
day or lose it.

Rows already in the ledger for sessions before 2026-09-16 are ESTIMATES with a
stated method, not recovered measurements. They should be read as an
order of magnitude.
"""

import argparse
import csv
import sys
from datetime import date
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
LEDGER = ROOT / "docs" / "token_ledger.csv"

# ---------------------------------------------------------------------------
# Prices, dated -- because an undated number cannot expire (L7 in the learning
# log). Checked 2026-09-16 against a pricing aggregator; the authority is
# https://docs.claude.com/en/docs/about-claude/pricing and these should be
# re-verified there before any figure from this tool is published.
# ---------------------------------------------------------------------------
PRICES_ASOF = "2026-09-16"
PRICES = {                      # US dollars per million tokens
    "claude-opus-5":   {"in": 5.00, "out": 25.00, "cache_read": 0.50,
                        "cache_write_1h": 10.00, "cache_write_5m": 6.25},
    "claude-sonnet-5": {"in": 2.00, "out": 10.00, "cache_read": 0.20,
                        "cache_write_1h": 4.00, "cache_write_5m": 2.50},
    "claude-haiku-4.5": {"in": 1.00, "out": 5.00, "cache_read": 0.10,
                         "cache_write_1h": 2.00, "cache_write_5m": 1.25},
}
DEFAULT_MODEL = "claude-opus-5"

# The subscription actually being paid for. Set this to what is true; the
# billed total is prorated across the months the project has run.
PLAN_USD_PER_MONTH = 200.0
PLAN_NAME = "Max"

# Factor-of-two band on estimated rows; see --estimate for why it is this wide.
EST_UNCERTAINTY = 2.0

FIELDS = ["date", "model", "input", "cache_read", "cache_write", "output",
          "basis", "note"]


def cost(row):
    """API-equivalent cost of one ledger row, in dollars."""
    p = PRICES.get(row["model"], PRICES[DEFAULT_MODEL])
    n = lambda k: float(row.get(k) or 0) / 1e6
    return (n("input") * p["in"]
            + n("cache_read") * p["cache_read"]
            + n("cache_write") * p["cache_write_1h"]
            + n("output") * p["out"])


def load():
    if not LEDGER.exists():
        return []
    with open(LEDGER, newline="", encoding="utf-8") as f:
        return [r for r in csv.DictReader(f) if r.get("date")]


def save(rows):
    LEDGER.parent.mkdir(parents=True, exist_ok=True)
    with open(LEDGER, "w", newline="", encoding="utf-8") as f:
        # lineterminator: csv defaults to \r\n, which git rewrites on every
        # commit and shows the whole file as changed.
        w = csv.DictWriter(f, fieldnames=FIELDS, lineterminator="\n")
        w.writeheader()
        for r in sorted(rows, key=lambda r: r["date"]):
            w.writerow({k: r.get(k, "") for k in FIELDS})


def add(args):
    rows = load()
    rows.append({
        "date": args.add,
        "model": args.model,
        "input": args.input,
        "cache_read": args.cache_read,
        "cache_write": args.cache_write,
        "output": args.output,
        "basis": args.basis,
        "note": args.note or "",
    })
    save(rows)
    r = rows[-1]
    print(f"added {r['date']}  {args.basis}  "
          f"API-equivalent ${cost(r):,.2f}")
    return 0


def report(args):
    rows = load()
    if not rows:
        print(f"No ledger yet at {LEDGER.relative_to(ROOT)}.")
        print("Add a session:  python tools/tokens.py --add <date> --in N "
              "--cache-read N --cache-write N --out N")
        return 1

    print(f"Token ledger  ({len(rows)} sessions, prices as of {PRICES_ASOF})\n")
    print(f"{'date':>11} {'basis':>9} {'input':>11} {'cache rd':>12} "
          f"{'cache wr':>11} {'output':>10} {'API-equiv':>11}")

    tot = {k: 0.0 for k in ("input", "cache_read", "cache_write", "output")}
    by_basis = {}
    for r in rows:
        c = cost(r)
        for k in tot:
            tot[k] += float(r.get(k) or 0)
        by_basis[r["basis"]] = by_basis.get(r["basis"], 0.0) + c
        print(f"{r['date']:>11} {r['basis']:>9} {float(r['input'] or 0):11,.0f} "
              f"{float(r['cache_read'] or 0):12,.0f} "
              f"{float(r['cache_write'] or 0):11,.0f} "
              f"{float(r['output'] or 0):10,.0f} {c:11,.2f}")
        if r.get("note"):
            print(f"{'':>11} {r['note'][:86]}")

    total_tokens = sum(tot.values())
    api = sum(by_basis.values())
    print("\n" + "-" * 80)
    print(f"{'TOTAL':>11} {'':>9} {tot['input']:11,.0f} "
          f"{tot['cache_read']:12,.0f} {tot['cache_write']:11,.0f} "
          f"{tot['output']:10,.0f} {api:11,.2f}")
    print(f"\n{total_tokens/1e6:,.1f} M tokens total.")

    for basis in sorted(by_basis):
        print(f"  {basis:<10} API-equivalent  ${by_basis[basis]:,.2f}")

    # Uncertainty band on the estimated part. The weak parameter is the average
    # conversation size re-sent per step, uncertain by roughly a factor of two,
    # and cache reads dominate the token count -- so the band is +/- that
    # factor rather than anything tighter. A measured row has no band.
    est = by_basis.get("estimated", 0.0)
    if est > 0:
        print(f"  {'':<10} estimated part is uncertain by about 2x: "
              f"${est/EST_UNCERTAINTY:,.0f} to ${est*EST_UNCERTAINTY:,.0f}")

    months = months_spanned(rows)
    billed = PLAN_USD_PER_MONTH * months
    print(f"\nBILLED      ${billed:,.2f}   "
          f"({PLAN_NAME} plan, {months} month(s) at "
          f"${PLAN_USD_PER_MONTH:,.0f})")
    print(f"API-EQUIV   ${api:,.2f}   "
          f"what the same tokens would cost at list prices")
    if api > 0 and billed > 0:
        ratio = api / billed
        print(f"\nList-price tokens per dollar billed: {ratio:.2f}x.")
        if ratio < 1:
            print("Below 1x, so at these prices the plan is not being used to "
                  "its list value --\nthe cost here is wall-clock and human "
                  "attention, not tokens.")
    print("\nAPI-EQUIV is a shadow price for the write-up. It is not money "
          "spent.")
    if any(r["basis"] == "estimated" for r in rows):
        print("Rows marked `estimated` are not recovered measurements -- see "
              "--estimate for the method.")
    return 0


def months_spanned(rows):
    ds = sorted(r["date"] for r in rows)
    a = date.fromisoformat(ds[0])
    b = date.fromisoformat(ds[-1])
    return max(1, (b.year - a.year) * 12 + (b.month - a.month) + 1)


def prices(args):
    print(f"Prices as of {PRICES_ASOF}, US dollars per million tokens.\n")
    print(f"{'model':>18} {'input':>8} {'output':>8} {'cache rd':>9} "
          f"{'cache wr 1h':>12}")
    for m, p in PRICES.items():
        print(f"{m:>18} {p['in']:8.2f} {p['out']:8.2f} "
              f"{p['cache_read']:9.2f} {p['cache_write_1h']:12.2f}")
    print("\nRe-verify at https://docs.claude.com/en/docs/about-claude/pricing "
          "before publishing any figure derived from these.")
    return 0


def estimate(args):
    print("""How the pre-2026-09-16 rows were estimated, and why they are weak

There is no way to recover a past session's token counts from inside a later
one, so everything before the ledger existed is an estimate. The method:

  1. `docs/PROMPT_LOG.md` records every human prompt -- a hard count of turns.
  2. An agentic turn in this project is dominated by CACHE READS, because the
     whole conversation is re-sent each turn under a one-hour cache TTL. Input
     per turn therefore grows roughly linearly with turn number within a
     session, not with the size of the prompt.
  3. Output per turn is bounded by what was produced: the repository is about
     90,000 words of code and documentation, and drafts and rewrites mean the
     generated total is several times the surviving total.
  4. Long-running measurement sessions add very little -- a 12-hour forecast
     sweep costs wall-clock, not tokens.

So the estimate is turn-count times an assumed average context, and the
assumed average is the weak link: it is uncertain by a factor of about two,
which propagates straight to the total. Treat the estimated rows as an order
of magnitude and nothing better.

WHAT TO DO INSTEAD, from now on

At the end of each session, read the usage breakdown and record it:

    python tools/tokens.py --add <today> --in N --cache-read N \\
        --cache-write N --out N --basis measured --note "<what the session did>"

This is the same shape of problem as the verification archive (P-07): the
number exists on the day and not afterwards.""")
    return 0


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--report", action="store_true")
    ap.add_argument("--prices", action="store_true")
    ap.add_argument("--estimate", action="store_true")
    ap.add_argument("--add", metavar="YYYY-MM-DD")
    ap.add_argument("--model", default=DEFAULT_MODEL)
    ap.add_argument("--input", "--in", dest="input", type=float, default=0)
    ap.add_argument("--cache-read", dest="cache_read", type=float, default=0)
    ap.add_argument("--cache-write", dest="cache_write", type=float, default=0)
    ap.add_argument("--output", "--out", dest="output", type=float, default=0)
    ap.add_argument("--basis", choices=["measured", "estimated"],
                    default="measured")
    ap.add_argument("--note", default="")
    args = ap.parse_args()

    if args.add:
        return add(args)
    if args.prices:
        return prices(args)
    if args.estimate:
        return estimate(args)
    return report(args)


if __name__ == "__main__":
    raise SystemExit(main())

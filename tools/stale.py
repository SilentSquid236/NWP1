#!/usr/bin/env python3
"""
Find measurements that the code has moved past.

    python tools/stale.py                 # dated measurements, oldest first
    python tools/stale.py --undated       # numbers with no date anywhere near
    python tools/stale.py --grace 30      # how far behind is "well behind"
    python tools/stale.py --all           # every dated measurement, stale or not

WHY THIS EXISTS

A constraint can expire silently (L7 in docs/LEARNING_LOG.md). "A deep sponge
flattens the jet" was true when it was measured. The sponge was then changed to
relax toward a frozen reference, which for a steady jet *is* the jet, so the
constraint was false from that moment -- and it was carried for four more days,
ruling out an experiment that turned out to be free.

Nothing detected that, because nothing could: the number was still sitting in
the comment where it had always been, and the code underneath it had changed.

WHAT THIS CAN AND CANNOT DO

It CANNOT know whether a claim is still true. Deciding that needs the run.
It also cannot always tell a date that timestamps a claim from a date that is
merely data -- an example in a --help string, a CSV fixture row. The obvious
ones are filtered; the rest show up as a handful of hits with the offending
line printed, which costs a reader seconds.

It CAN know that nothing has re-checked it: the file carrying the measurement
was modified well after the measurement was taken. That is a weak signal and
it is deliberately weak. A false positive costs one re-run. The false negative
this is built against cost four days.

THE PART THAT MATTERS MORE THAN THE CHECK

`--undated`. A number with no date cannot expire, because nothing can tell
when the code moved past it. Dating the measurements is what makes the
staleness check worth anything at all; until then this tool is looking at a
small corner of the evidence and reporting that corner as if it were the
whole.

HOW THE DATE OF A FILE IS DECIDED

git, when it is available: the commit date of the last commit that touched the
file. That is the honest answer and it is what to use.

Otherwise the filesystem mtime, which is all the server has -- git is not
installed there and cannot be. MTIMES ARE MUCH WEAKER, and the weakness runs
one way: a fresh copy stamps every file at once, so everything looks as if it
changed today and every measurement older than the grace period is reported.
Measured on this repository 2026-09-10, immediately after a clone: 0 hits from
git dates, 12 from mtimes, and all twelve were files nobody had touched. So
on the server, read the output as "these are the numbers whose files are in
the batch that moved", not as twelve separate findings.
"""

import argparse
import re
import subprocess
import tokenize
import sys
from datetime import date, datetime, timedelta
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

SKIP_DIRS = {".git", "data", "__pycache__", ".ipynb_checkpoints"}
SUFFIXES = {".py", ".md", ".sh"}
# Generated files. Their numbers are checksums and byte counts, which are
# recomputed by the generator and cannot go stale in the sense meant here.
#
# APPEND-ONLY RECORDS are excluded for a different reason, and it is the
# tool's sharpest limitation. A research-log entry is frozen on the day it
# is written; that is the point of the log. But appending TODAY's entry
# updates the whole FILE's date, so every historical measurement in it then
# looks like a number sitting in a file that moved. On 2026-09-12 that was
# twelve of twelve hits -- all of them entries from August that nobody has
# any business re-running. A check with a 100% false-positive rate gets
# switched off, so these are out of scope until staleness can be scoped to
# a SECTION rather than a file.
#
# PROBLEMS.md and LEARNING_LOG.md stay in scope deliberately: both are
# edited in place when a status changes, so a number in them really can be
# left behind by the code.
SKIP_FILES = {"docs/MANIFEST.txt", "docs/STRUCTURE.md",
              "docs/RESEARCH_LOG.md", "docs/PROMPT_LOG.md"}

# An ISO date written anywhere: a comment, a log heading, a docstring.
DATE = re.compile(r"(?<!\d)(20\d\d)-(\d\d)-(\d\d)(?!\d)")

# A date that is DATA, not a timestamp on a claim: the example date in a
# command line, a CSV fixture row, a clock time. Reading one of those as
# "when this was measured" is how the first version of this tool dated every
# number in ingest_hrrr.py to the example date in its --start help string,
# and then called them all stale.
DATE_IS_DATA = re.compile(
    r"--start|e\.g\.|\d{4}-\d\d-\d\d[T ]\d\d:?\d\d|^[^,]*,[^,]*,[^,]*,")

# A line that REPORTS something measured, rather than merely containing a
# digit: units, survival counts, ratios, orders of magnitude. Deliberately
# narrow. A checker that flags loop bounds and array indices gets ignored,
# which is worse than not having one.
MEASUREMENT = re.compile(
    r"""
      \b\d+\s*/\s*\d+\s*(?:h\b|hours?\b)?     # survival counts: 11/12, 8/12 h
    | \d\s*(?:m/s|m\^2/s|m2/s|hPa|kPa|\bPa\b|\bkm\b|\bK\b|\bm\b|\bs\b)
    | \d\s*%                                  # percentages
    | \d[eE][-+]?\d                            # 3.90e-05
    | \bRi\s*[=~<>]+\s*-?\d                   # Ri = 0.23
    | \be-folding\b
    """,
    re.X,
)

# Lines that look like a measurement but are template, schema, or prose about
# how to write one -- not a number anybody took.
NOISE = re.compile(r"<[a-z ]+>|YYYY-MM-DD|\bregenerate\b|\bpython tools/", re.I)


def prose_lines(path):
    """
    Line numbers carrying PROSE, where a measurement gets written down.

    For Python that means comments and string literals only. Code is excluded
    on purpose: `if hours == 12` and `0.514444` in a unit-conversion assertion
    are not measurements of anything, and an earlier version of this tool that
    read every line reported 991 "numbers", which is the same as reporting
    none. Everything else -- markdown, shell, text -- is prose throughout.
    """
    if path.suffix != ".py":
        return None                      # None = every line counts
    keep = {}
    try:
        with path.open("rb") as fh:
            for tok in tokenize.tokenize(fh.readline):
                if tok.type not in (tokenize.COMMENT, tokenize.STRING):
                    continue
                if tok.start[0] == tok.end[0]:
                    # A one-line comment or literal: keep the token text, not
                    # the code it shares the line with. `m.h += 1.0 * ...  #
                    # small => linear` is a line of code with a comment on the
                    # end, and only the comment is prose.
                    keep.setdefault(tok.start[0], "")
                    keep[tok.start[0]] += " " + tok.string
                else:
                    for ln in range(tok.start[0], tok.end[0] + 1):
                        keep.setdefault(ln, None)      # None = whole line
    except (tokenize.TokenError, SyntaxError, OSError, UnicodeDecodeError):
        return None
    return keep


def files():
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file() or p.suffix not in SUFFIXES:
            continue
        rel = p.relative_to(ROOT)
        if SKIP_DIRS & set(rel.parts) or rel.as_posix() in SKIP_FILES:
            continue
        yield p


def file_date(path, use_git=True):
    """(date, source) for the file's last modification."""
    if use_git:
        try:
            out = subprocess.run(
                ["git", "-C", str(ROOT), "log", "-1", "--format=%cs", "--",
                 str(path.relative_to(ROOT))],
                capture_output=True, text=True, timeout=10)
            if out.returncode == 0 and out.stdout.strip():
                return date.fromisoformat(out.stdout.strip()), "git"
        except (OSError, ValueError, subprocess.SubprocessError):
            pass
    return date.fromtimestamp(path.stat().st_mtime), "mtime"


def scan(path, window):
    """
    (dated, undated) measurement lines in one file.

    A measurement is DATED if a date appears within `window` lines of it --
    that is the same "anywhere near it" a reader applies, and it is why the
    window is generous. Undated ones are reported separately rather than
    guessed at.
    """
    try:
        lines = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return [], []

    dates_at = {}
    for i, line in enumerate(lines):
        if DATE_IS_DATA.search(line):
            continue
        m = DATE.search(line)
        if m:
            try:
                dates_at[i] = date(int(m.group(1)), int(m.group(2)),
                                   int(m.group(3)))
            except ValueError:
                pass

    # A DATED SECTION dates everything inside it. In the research log a
    # heading is "## 2026-09-05 -- The sponge does not absorb weather", and
    # every number for the next two hundred lines belongs to that day. A
    # fixed line window would call almost all of them undated, which is both
    # wrong and the kind of noise that gets a checker switched off.
    section = {}
    if path.suffix == ".md":
        current = None
        for i, line in enumerate(lines):
            if line.startswith("#"):
                m = DATE.search(line)
                current = None if m is None else dates_at.get(i)
            if current is not None:
                section[i] = current

    prose = prose_lines(path)

    dated, undated = [], []
    for i, line in enumerate(lines):
        text = line
        if prose is not None:
            if (i + 1) not in prose:
                continue
            text = prose[i + 1] or line
        if not MEASUREMENT.search(text) or NOISE.search(text):
            continue
        near = [(abs(i - j), d) for j, d in dates_at.items()
                if abs(i - j) <= window]
        when = min(near)[1] if near else section.get(i)
        if when is not None:
            dated.append((i + 1, when, line.strip()))
        else:
            undated.append((i + 1, line.strip()))
    return dated, undated


def trim(s, n=88):
    s = " ".join(s.split())
    return s if len(s) <= n else s[:n - 1] + "…"


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--grace", type=int, default=14,
                    help="days a file may move after a measurement before it "
                         "is called stale (default 14)")
    ap.add_argument("--window", type=int, default=12,
                    help="how many lines away a date still counts as being "
                         "near a number (default 12)")
    ap.add_argument("--undated", action="store_true",
                    help="list numbers with no date near them instead")
    ap.add_argument("--all", action="store_true",
                    help="list every dated measurement, stale or not")
    ap.add_argument("--no-git", action="store_true",
                    help="use filesystem mtimes, as the server must")
    args = ap.parse_args()

    n_dated = n_undated = 0
    stale, undated_rows, source_used = [], [], set()

    for path in files():
        dated, undated = scan(path, args.window)
        n_dated += len(dated)
        n_undated += len(undated)
        rel = path.relative_to(ROOT)
        if undated:
            undated_rows.append((rel, undated))
        if not dated:
            continue
        mod, source = file_date(path, use_git=not args.no_git)
        source_used.add(source)
        for lineno, taken, text in dated:
            behind = (mod - taken).days
            if args.all or behind > args.grace:
                stale.append((behind, rel, lineno, taken, mod, text))

    if args.undated:
        for rel, rows in undated_rows:
            print(f"\n{rel}")
            for lineno, text in rows:
                print(f"  {lineno:>5}  {trim(text)}")
        print(f"\n{n_undated} number(s) with no date within "
              f"{args.window} lines, in {len(undated_rows)} file(s).")
        print("A number with no date cannot expire. Dating them is what makes "
              "the staleness check\nworth anything.")
        return 0

    stale.sort(reverse=True)
    for behind, rel, lineno, taken, mod, text in stale:
        print(f"\n{rel}:{lineno}")
        print(f"  measured {taken}, file last changed {mod} "
              f"({behind} days later)")
        print(f"  {trim(text)}")

    print(f"\n{len(stale)} of {n_dated} dated measurement(s) sit in files "
          f"changed more than {args.grace} days later.")
    print(f"{n_undated} more report numbers with no date near them "
          f"(--undated).")
    if source_used:
        print(f"file dates from: {', '.join(sorted(source_used))}")
    print("\nThis says nothing has re-checked them. It does not say they are "
          "wrong.")
    return 1 if stale else 0


if __name__ == "__main__":
    sys.exit(main())

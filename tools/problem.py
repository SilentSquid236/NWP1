#!/usr/bin/env python3
"""
Maintain docs/PROBLEMS.md.

    python tools/problem.py new "Sponge base reflects"     # append an OPEN entry
    python tools/problem.py list                           # open problems only
    python tools/problem.py list --all                     # every entry
    python tools/problem.py check                          # audit the register

`new` gives the next free P-nn and writes the skeleton into the OPEN section.
Closing a problem is a manual edit, deliberately: the entry has to gain the
measurement that confirms the fix, and a script cannot supply that.

`check` is the part worth running before a commit. It flags:
  * FIXED entries with no "Confirmed by" line -- a fix without a number
  * OPEN entries with no "Symptom"
  * duplicate or missing problem numbers
"""

import re
import sys
from datetime import date
from pathlib import Path

DOC = Path(__file__).resolve().parent.parent / "docs" / "PROBLEMS.md"

TEMPLATE = """## P-{num:02d} — {title}
**Category** ? · **First seen** {date} · **Status** OPEN

**Symptom.** <what is observed, with numbers>

**What is known.** <measurements so far>

**Ruled out.** <candidates already eliminated, and by what measurement>

---

"""

HEAD = re.compile(r"^## (P-(\d+)) — (.+)$", re.M)

# The ELIMINATED section is a table rather than a run of headed entries: a
# ruled-out candidate needs one line and one measurement, not five paragraphs.
# The checker still has to see those rows, or every one of them looks like a
# gap in the numbering.
ROW = re.compile(r"^\| P-(\d+) \| (.+?) \| (.+?) \|\s*$", re.M)


def entries(text):
    """(number, title, body) for every entry, headed or tabulated."""
    out = []
    marks = list(HEAD.finditer(text))
    for i, m in enumerate(marks):
        end = marks[i + 1].start() if i + 1 < len(marks) else len(text)
        out.append((int(m.group(2)), m.group(3), text[m.start():end]))
    for m in ROW.finditer(text):
        title = m.group(2).replace("**", "").strip()
        body = f"**Status** ELIMINATED\n**Confirmed by.** {m.group(3).strip()}"
        out.append((int(m.group(1)), title, body))
    return sorted(out)


def status_of(body):
    m = re.search(r"\*\*Status\*\*\s+([A-Z]+)", body)
    if m:
        return m.group(1)
    for s in ("FIXED", "ELIMINATED", "REVERTED", "ACCEPTED"):
        if f"**{s}**" in body or f"· **{s}" in body:
            return s
    return "FIXED" if "**Fixed**" in body else "OPEN"


def cmd_new(title):
    text = DOC.read_text(encoding="utf-8")
    nums = [n for n, _, _ in entries(text)]
    num = max(nums) + 1 if nums else 1
    block = TEMPLATE.format(num=num, title=title, date=date.today().isoformat())
    marker = "\n# FIXED\n"
    if marker not in text:
        raise SystemExit("could not find the FIXED section header")
    text = text.replace(marker, "\n" + block + marker, 1)
    DOC.write_text(text, encoding="utf-8")
    print(f"added P-{num:02d} — {title}  (OPEN)")
    return 0


def cmd_list(show_all):
    text = DOC.read_text(encoding="utf-8")
    for num, title, body in entries(text):
        st = status_of(body)
        if show_all or st == "OPEN":
            print(f"  P-{num:02d}  {st:<11} {title}")
    return 0


def cmd_check():
    text = DOC.read_text(encoding="utf-8")
    problems, seen = [], set()
    for num, title, body in entries(text):
        st = status_of(body)
        if num in seen:
            problems.append(f"P-{num:02d} duplicate number")
        seen.add(num)
        if st == "FIXED" and "Confirmed by" not in body:
            problems.append(f"P-{num:02d} '{title}' is FIXED with no "
                            f"'Confirmed by' measurement")
        if st == "OPEN" and "Symptom" not in body:
            problems.append(f"P-{num:02d} '{title}' is OPEN with no Symptom")
    missing = sorted(set(range(1, max(seen) + 1)) - seen) if seen else []
    if missing:
        problems.append("gaps in numbering: "
                        + ", ".join(f"P-{n:02d}" for n in missing))
    if problems:
        for p in problems:
            print("  " + p)
        print(f"\n{len(problems)} issue(s) in the register")
        return 1
    print(f"  register clean: {len(seen)} entries, "
          f"{sum(1 for _, _, b in entries(text) if status_of(b) == 'OPEN')} open")
    return 0


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1
    cmd = sys.argv[1]
    if cmd == "new" and len(sys.argv) > 2:
        return cmd_new(" ".join(sys.argv[2:]))
    if cmd == "list":
        return cmd_list("--all" in sys.argv)
    if cmd == "check":
        return cmd_check()
    print(__doc__)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

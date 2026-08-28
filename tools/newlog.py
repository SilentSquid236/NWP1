#!/usr/bin/env python3
"""
Append a dated, templated entry to docs/RESEARCH_LOG.md.

    python tools/newlog.py "Sigma coordinate prototype"

Writes the skeleton with today's date and opens nothing -- edit the file.
The point is to make logging cheap enough that it actually happens, and to
keep every entry in the same shape so the log stays scannable.
"""

import sys
from datetime import date
from pathlib import Path

TEMPLATE = """
## {date} — {title}

**Context.** <what prompted this>

**Hypothesis.** <what you expect, written BEFORE the result>

**Method.** <what was run; enough to reproduce>

**Result.** <numbers; a table if more than one>

**Interpretation.** <what it means, and what it does not>

**Status.** <kept / reverted / open — and why>

---
"""


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        return 1

    title = " ".join(sys.argv[1:])
    log = Path(__file__).resolve().parent.parent / "docs" / "RESEARCH_LOG.md"
    if not log.exists():
        print(f"{log} not found")
        return 1

    text = log.read_text()
    marker = "## Template for new entries"
    entry = TEMPLATE.format(date=date.today().isoformat(), title=title)

    if marker in text:
        # Insert before the template section so entries stay in one block.
        head, tail = text.split(marker, 1)
        log.write_text(head.rstrip() + "\n\n" + entry.strip() + "\n\n" + marker + tail)
    else:
        log.write_text(text.rstrip() + "\n\n" + entry.strip() + "\n")

    print(f"Added entry '{title}' ({date.today().isoformat()}) to {log}")
    print("Now fill it in -- hypothesis BEFORE result.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Check that the code and docs cite outside work the AMS way.

    python tools/check_refs.py

Scans every .py and .md file (except the generated MANIFEST and STRUCTURE
and the verbatim PROMPT_LOG) for author-year citations such as
"Davies (1976)", "(Wicker and Skamarock 2002)", "Koch et al. (1983)" or
"NOAA (2026b)". It reports:
  * a citation whose first author (or, for an organisation, the last word
    of its name) and year have no entry in docs/REFERENCES.md;
  * an "&" between author names, which AMS style writes as "and";
  * an entry in docs/REFERENCES.md that nothing cites.
Exit status 1 if anything is reported.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
REFS = ROOT / "docs" / "REFERENCES.md"
SKIP = {"docs/MANIFEST.txt", "docs/STRUCTURE.md", "docs/PROMPT_LOG.md", "docs/REFERENCES.md"}
NOT_NAMES = {"January", "February", "March", "April", "May", "June", "July", "August",
             "September", "October", "November", "December", "Jan", "Feb", "Mar", "Apr",
             "Jun", "Jul", "Aug", "Sep", "Sept", "Oct", "Nov", "Dec", "Python", "Since",
             "In", "From", "Of", "The", "On", "Until", "Before", "After", "By"}

CITE = re.compile(r"\b([A-Z][A-Za-z'\-]+)(?:,? and [A-Z][A-Za-z'\-]+| et al\.)?\s\(?"
                  r"((?:19|20)\d\d)([a-z])?(?:, ([a-z]))?(?=[);,.:\s]|$)")
AMP = re.compile(r"\b[A-Z][a-z]+ & [A-Z][a-z]+,? \(?(?:19|20)\d\d")
ENTRY = re.compile(r"^(?P<who>[^\n]+?), (?P<year>(?:19|20)\d\d)(?P<suf>[a-z])?: ", re.M)


def entries():
    """{(key, year+suffix)} for every reference; key = first surname or last word."""
    out = {}
    for m in ENTRY.finditer(REFS.read_text(encoding="utf-8")):
        who = m.group("who")
        person = re.match(r"([A-Z][A-Za-z'\-]+), [A-Z]\.", who)
        key = person.group(1) if person else who.split()[-1]
        out[(key, m.group("year") + (m.group("suf") or ""))] = who
    return out


def main():
    known = entries()
    used, problems = set(), []
    files = [p for p in ROOT.rglob("*") if p.suffix in (".py", ".md")
             and "data" not in p.relative_to(ROOT).parts[:1]
             and "__pycache__" not in p.parts]
    for p in sorted(files):
        rel = p.relative_to(ROOT).as_posix()
        if rel in SKIP:
            continue
        text = p.read_text(encoding="utf-8", errors="replace")
        for n, line in enumerate(text.splitlines(), 1):
            for m in AMP.finditer(line):
                problems.append(f"{rel}:{n}: '&' in a citation ({m.group(0)}); write 'and'")
            for m in CITE.finditer(line):
                name, year, s1, s2 = m.groups()
                if name in NOT_NAMES:
                    continue
                keys = [(name, year + s) for s in {s1 or "", s2 or ""}]
                if not any(k in known for k in keys) and (name, year) not in known:
                    problems.append(f"{rel}:{n}: '{m.group(0).rstrip(');,. ')}' "
                                    f"is not in docs/REFERENCES.md")
                    continue
                used.update(k for k in keys if k in known)
                if (name, year) in known:
                    used.add((name, year))
                if s1 and not s2:
                    used.add((name, year + s1))
    for k, who in sorted(known.items()):
        if k not in used:
            problems.append(f"docs/REFERENCES.md: '{who}, {k[1]}' is never cited")
    for msg in problems:
        print(msg)
    print(f"references: {len(known)} entries, {len(used)} cited, {len(problems)} problems")
    return 1 if problems else 0


if __name__ == "__main__":
    sys.exit(main())

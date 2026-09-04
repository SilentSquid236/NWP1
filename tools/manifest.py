#!/usr/bin/env python3
"""
Write or check a manifest of every file in the project.

    python tools/manifest.py            # write docs/MANIFEST.txt
    python tools/manifest.py --check    # compare this copy against it

WHY

`tools/pull.sh` fetches over curl, and curl only ever adds and overwrites. It
cannot tell you that a file is MISSING, because a file that never arrived
looks exactly like a file that was never in the repo. Neither can git, if the
file was never committed in the first place -- which is the more common cause
here: work delivered into a chat and saved locally is not in the repository
until someone pushes it.

`--check` answers the question pull.sh cannot: does this copy contain what the
project is supposed to contain, byte for byte.

Excludes data/ (the verification archive, which is machine-specific and
irreplaceable) and caches.
"""

import hashlib
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
OUT = ROOT / "docs" / "MANIFEST.txt"

SKIP_DIRS = {"__pycache__", ".git", "data", ".idea", ".vscode"}
SKIP_EXT = {".pyc", ".log", ".npz", ".npy", ".gz", ".zip"}


def files():
    for p in sorted(ROOT.rglob("*")):
        if not p.is_file():
            continue
        if any(part in SKIP_DIRS for part in p.parts):
            continue
        if p.suffix in SKIP_EXT:
            continue
        if p == OUT:
            continue
        yield p


def digest(p):
    return hashlib.sha256(p.read_bytes()).hexdigest()[:16]


def write():
    rows = [f"{digest(p)}  {p.stat().st_size:>8}  "
            f"{p.relative_to(ROOT).as_posix()}" for p in files()]
    OUT.parent.mkdir(parents=True, exist_ok=True)
    OUT.write_text(
        "# Manifest of every source, doc and tool file in the project.\n"
        "# Regenerate with: python tools/manifest.py\n"
        "# Check a copy with: python tools/manifest.py --check\n"
        "#\n"
        "# sha256(16)          bytes  path\n"
        + "\n".join(rows) + "\n")
    print(f"wrote {OUT.relative_to(ROOT)} — {len(rows)} files")
    return 0


def check():
    if not OUT.exists():
        print(f"No {OUT.relative_to(ROOT)}. Run without --check on a known-good "
              f"copy first.")
        return 1

    want = {}
    for line in OUT.read_text().splitlines():
        if line.startswith("#") or not line.strip():
            continue
        h, size, path = line.split(None, 2)
        want[path] = (h, int(size))

    have = {p.relative_to(ROOT).as_posix(): p for p in files()}

    missing = sorted(set(want) - set(have))
    extra = sorted(set(have) - set(want))
    differ = sorted(k for k in set(want) & set(have)
                    if digest(have[k]) != want[k][0])

    for k in missing:
        print(f"  MISSING  {k}  ({want[k][1]} bytes)")
    for k in differ:
        print(f"  DIFFERS  {k}")
    for k in extra:
        print(f"  EXTRA    {k}")

    if not (missing or differ or extra):
        print(f"  complete — {len(want)} files, all matching")
        return 0
    print(f"\n{len(missing)} missing, {len(differ)} differing, {len(extra)} extra")
    return 1


if __name__ == "__main__":
    raise SystemExit(check() if "--check" in sys.argv else write())

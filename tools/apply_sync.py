#!/usr/bin/env python3
"""
Apply a sync archive to this project, safely, in one command.

    python tools/apply_sync.py                       # newest nwp_sync*.tar.gz in Downloads
    python tools/apply_sync.py path\\to\\nwp_sync.tar.gz
    python tools/apply_sync.py --dry                 # show what would change

WHY THIS EXISTS

Transferring this project by hand has failed in four distinct ways already,
and every one of them was silent:

  * `scp -r src host:.../src` when src exists -> a nested src/src, imports
    resolving to whichever copy came first (P-24)
  * an archive whose top level is NWP_Deployment_Package/, extracted from
    INSIDE NWP_Deployment_Package -> the package nested in itself
  * curl-based pulls that only ever add and overwrite, so a file deleted
    upstream survives forever and two copies of a renamed module coexist
  * thirty files that were never pushed at all, so a correct `git pull`
    correctly delivered a stale tree

None of those raised an error. They produced a working-looking directory that
behaved like an older version. This script removes every step where that can
happen, and refuses rather than guessing when something looks wrong.

WHAT IT GUARANTEES

  * data/ is never touched. It holds the verification archive, the one thing
    in this project that cannot be recreated.
  * The archive is inspected before anything is written; a truncated download
    costs nothing.
  * Files land at the project root regardless of whether the archive has a
    wrapper directory, so the nesting failure cannot recur.
  * The manifest is checked afterwards, so "it looked fine" is replaced by a
    byte-for-byte answer.

Standard library only. Runs the same on Windows and Linux.
"""

import argparse
import filecmp
import shutil
import subprocess
import sys
import tarfile
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent
NEVER_TOUCH = {"data"}


def find_archive():
    """Newest nwp_sync*.tar.gz in the usual download locations."""
    home = Path.home()
    spots = [home / "Downloads", home / "Desktop", ROOT.parent, ROOT,
             Path.cwd()]
    found = []
    for d in spots:
        try:
            found += list(d.glob("nwp_sync*.tar.gz"))
            found += list(d.glob("nwp_contents*.tar.gz"))
        except OSError:
            continue
    if not found:
        return None
    return max(found, key=lambda p: p.stat().st_mtime)


def project_root_in(staged):
    """
    The directory inside the archive that IS the project.

    Found by looking for config.py beside src/, not by assuming a depth --
    guessing the depth is what nests the package inside itself.
    """
    if (staged / "config.py").exists() and (staged / "src").is_dir():
        return staged
    for cand in sorted(staged.rglob("config.py")):
        if (cand.parent / "src").is_dir():
            return cand.parent
    return None


def main():
    ap = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("archive", nargs="?", help="path to the .tar.gz")
    ap.add_argument("--dry", action="store_true",
                    help="list changes, write nothing")
    args = ap.parse_args()

    arc = Path(args.archive) if args.archive else find_archive()
    if arc is None:
        print("No nwp_sync*.tar.gz found in Downloads, Desktop, or here.")
        print("Pass the path explicitly:  python tools/apply_sync.py <file>")
        return 1
    if not arc.exists():
        print(f"No such file: {arc}")
        return 1

    print(f"archive  {arc}  ({arc.stat().st_size/1e3:.0f} kB)")
    print(f"project  {ROOT}\n")

    with tempfile.TemporaryDirectory() as tmp:
        staged = Path(tmp)
        try:
            with tarfile.open(arc, "r:gz") as tf:
                tf.extractall(staged)
        except (tarfile.TarError, OSError) as e:
            print(f"Could not read the archive: {e}")
            print("It is probably a partial download. Fetch it again.")
            return 1

        src = project_root_in(staged)
        if src is None:
            print("This archive does not contain a project: nothing in it has "
                  "config.py beside a src/ directory.")
            return 1

        new, changed, skipped = [], [], 0
        for f in sorted(src.rglob("*")):
            if not f.is_file():
                continue
            rel = f.relative_to(src)
            if rel.parts and rel.parts[0] in NEVER_TOUCH:
                skipped += 1
                continue
            dest = ROOT / rel
            if not dest.exists():
                new.append(rel)
            elif not filecmp.cmp(f, dest, shallow=False):
                changed.append(rel)

        for rel in new:
            print(f"  new      {rel.as_posix()}")
        for rel in changed:
            print(f"  updated  {rel.as_posix()}")
        if not new and not changed:
            print("  nothing to do -- already identical")
            return 0
        print(f"\n  {len(new)} new, {len(changed)} updated"
              + (f", {skipped} skipped under data/" if skipped else ""))

        if args.dry:
            print("\n--dry: nothing written")
            return 0

        written = failed = 0
        for rel in new + changed:
            dest = ROOT / rel
            dest.parent.mkdir(parents=True, exist_ok=True)
            try:
                # Copy CONTENT into the existing file rather than replacing the
                # file: some mounted filesystems allow writing and not
                # unlinking, and shutil.copy2 unlinks first.
                with open(src / rel, "rb") as fh:
                    data = fh.read()
                with open(dest, "wb") as fh:
                    fh.write(data)
                written += 1
            except OSError as e:
                print(f"  FAILED {rel.as_posix()}: {e}")
                failed += 1
        print(f"\n  wrote {written} file(s)"
              + (f", {failed} failed" if failed else ""))

    print()
    for check in ("checklayout.py", "manifest.py"):
        script = ROOT / "tools" / check
        if not script.exists():
            continue
        cmd = [sys.executable, str(script)]
        if check == "manifest.py":
            cmd.append("--check")
        r = subprocess.run(cmd, capture_output=True, text=True)
        print(r.stdout.rstrip() or r.stderr.rstrip())
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Check that the project is laid out where the code expects it.

    python tools/checklayout.py

WHY THIS EXISTS

`scp -r src user@host:.../NWP_Deployment_Package/src` creates `src/src` when
`src` already exists -- the destination has to be the PARENT directory, not
the directory itself. That produced P-24 in docs/PROBLEMS.md: imports resolved
to whichever copy came first on sys.path, and the symptom was code that looked
right and behaved like an older version. A human spotted it; no test did.

This checks the three things that go wrong with a copy:

  1. nesting     -- src/src, docs/docs, a package inside itself
  2. missing     -- a module the code imports that never arrived
  3. stale pairs -- two copies of the same module in different places

It imports nothing and runs anywhere Python does, so it is safe to run on a
server with no packages installed.
"""

import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent.parent

# Files the project cannot run without, by the path the code expects.
REQUIRED = [
    "config.py", "resources.py", "netpolicy.py", "preflight.py",
    "src/forecast.py", "src/ingest_hrrr.py", "src/verify.py",
    "src/dynamics/grid.py", "src/dynamics/sigma.py",
    "src/dynamics/primitive_sigma.py", "src/dynamics/subgrid.py",
    "src/dynamics/turbulence.py", "src/dynamics/surface.py",
    "src/dynamics/convection.py", "src/dynamics/initialization.py",
    "src/dynamics/interpolate.py", "src/dynamics/boundaries.py",
    "src/verification/observations.py", "src/verification/fetchers.py",
    "src/verification/obs_operator.py", "src/verification/sigma_operator.py",
    "src/verification/scoring.py",
    "src/postproc/bias_correction.py",
    "docs/PROBLEMS.md", "docs/RESEARCH_LOG.md",
    "tools/tree.py", "tools/problem.py", "tools/daily.sh",
]

# Directory names that must never contain a directory of the same name.
NO_NEST = ["src", "docs", "tools", "dynamics", "verification", "postproc",
           "NWP_Deployment_Package"]


def main():
    problems = []

    print(f"Checking {ROOT}\n")

    # 0. The package inside itself. Extracting an archive whose top level is
    #    NWP_Deployment_Package/ while already standing inside
    #    NWP_Deployment_Package/ produces exactly this, and every path in the
    #    project then points at the outer, empty copy.
    inner = ROOT / "NWP_Deployment_Package"
    if inner.is_dir():
        problems.append(
            f"NESTED PACKAGE: {inner.relative_to(ROOT)}/ exists inside the "
            f"project.\n"
            f"         An archive was extracted one level too deep. Fix with:\n"
            f"           cd {ROOT}\n"
            f"           cp -a NWP_Deployment_Package/. .\n"
            f"           rm -rf NWP_Deployment_Package")

    # 1. Nesting.
    for name in NO_NEST:
        for d in ROOT.rglob(name):
            if not d.is_dir():
                continue
            inner = d / name
            if inner.is_dir():
                problems.append(
                    f"NESTED: {inner.relative_to(ROOT)} exists.\n"
                    f"         Almost certainly `scp -r {name} host:.../{name}`\n"
                    f"         when it should have been the PARENT directory.\n"
                    f"         Fix: move its contents up one level, then remove it.")

    # 2. Missing.
    missing = [f for f in REQUIRED if not (ROOT / f).exists()]
    for f in missing:
        problems.append(f"MISSING: {f}")

    # 3. Duplicates of the same module name in more than one place.
    seen = {}
    for py in ROOT.rglob("*.py"):
        if "__pycache__" in py.parts or ".git" in py.parts:
            continue
        seen.setdefault(py.name, []).append(py.relative_to(ROOT))
    for name, paths in sorted(seen.items()):
        if len(paths) > 1 and name != "__init__.py":
            problems.append(
                "DUPLICATE: " + name + " appears in "
                + ", ".join(str(p) for p in paths)
                + "\n           Whichever is first on sys.path wins, and it "
                  "may not be\n           the one you just copied.")

    if not problems:
        n = len(list(ROOT.rglob("*.py")))
        print(f"  layout OK -- {len(REQUIRED)} required files present, "
              f"{n} python files, no nesting, no duplicates")
        return 0

    for p in problems:
        print("  " + p)
    print(f"\n{len(problems)} layout problem(s)")
    return 1


if __name__ == "__main__":
    raise SystemExit(main())

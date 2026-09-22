---
name: nwp-sync
description: Use when moving NWP1 work between places — the GitHub repository, the Windows desktop clone, the shared compute server, or the claude.ai project docs — or when files are missing, nested wrongly, or out of date after a transfer.
---

# Syncing NWP1

## The places, and which is the source of truth

| where | holds | reached by |
|---|---|---|
| GitHub `SilentSquid236/NWP1` | code, `docs/`, `tools/`, `skills/` — **source of truth** | git on Windows; `git pull` (or `tools/pull.sh`) on the server |
| Windows desktop `Desktop\NWP\NWP_Deployment_Package` | the clone with push credentials | git |
| shared server (the Xeon), `/data5/pierce/NWP` | the working copy that runs; `data/` is at `NWP_Deployment_Package/data` inside it | `git pull` — git is on the server since 2026-09-22; `tools/pull.sh` (curl) is the fallback |

`data/` exists only on the server and is **never overwritten** by any sync.

## Getting work into the repository

A session that cannot push produces a **patch**, not a tarball of files:

    git format-patch <base>..HEAD --stdout > nwp1-<topic>-<date>.patch

and it is applied from Windows, where credentials and file deletion both work:

    git checkout main && git pull
    git checkout -b <topic>
    git am nwp1-<topic>-<date>.patch
    git push -u origin <topic>

Do **not** write git history through a device bridge that cannot delete files:
git replaces `.git/index` by unlinking it, so `reset`, `am` and `checkout`
fail part-way and leave `.git/*.lock` files. Read-only git over such a bridge
is fine.

## Getting the repository onto the server

git has been on the server since 2026-09-22 (reported by the human; before
that it was absent and `tools/pull.sh` was the only route).

**Server layout, checked 2026-09-22.** Run git in **`/data5/pierce/NWP`** — not
`/data5/pierce/Data5/NWP`, which is an old partial copy. `/data5/pierce/NWP`
is already a git checkout: git 2.52.0, `origin` =
`https://github.com/SilentSquid236/NWP1.git`, branch `main` tracking
`origin/main`, clean apart from the untracked `NWP1-main/`,
`NWP_Deployment_Package/` and `nwp.tar.gz` (all checked 2026-09-22). The
conversion below is therefore **not needed there**; `git pull --ff-only` is.
`NWP_DATA_ROOT` (in `~/.bashrc`) is
`/data5/pierce/NWP/NWP_Deployment_Package/data`: the archive lives inside a
nested copy of the package, which has its own `.git`. **Never delete, move or
`git clean` that nested folder.** `manifest.py --check` in the root reports its
files as EXTRA; that is expected and is not a reason to remove them.

Kept for any copy that has no usable `.git` — convert it **in place, once**,
never by deleting or re-cloning the directory:

    git --version                  # record it; the recipe below needs >= 1.8
    git init
    git remote add origin https://github.com/SilentSquid236/NWP1.git
    git fetch origin main
    git reset origin/main          # moves the index only; no file is touched
    git branch -M main
    git branch --set-upstream-to=origin/main main
    git status                     # read this list before the next line
    git checkout -- .              # overwrites tracked files that differ
    python tools/manifest.py --check

After that, an update is

    git pull --ff-only
    python tools/manifest.py --check

The server **only pulls**. GitHub stays the source of truth and pushes happen
from Windows. `git fetch` is not rate-limited by `netpolicy.py`; the repository
is small (the 2026-09-22 zip is 0.67 MB), so that is within the bandwidth rule,
but do not fetch large refs or LFS content this way.

`data/` is in `.gitignore`, so git neither tracks nor overwrites it — and for
exactly that reason **`git clean -x` and `git clean -X` delete it, and
`git stash --all` sweeps it out of the working tree. Never run them on the
server.** Plain `git clean -f` (no `-x`)
leaves ignored files alone, but there is no reason to run it there either.

The conversion recipe has never been run: the server root turned out to be a
checkout already. Treat it as written, not exercised. The untracked copies are
unprotected against `git clean -f -d`, which would remove `NWP1-main/` and
`nwp.tar.gz` (a nested repository like `NWP_Deployment_Package/` survives a
single `-f` but not `-ff`). There is no reason to run `git clean` on the
server at all.

`tools/pull.sh` still works and is the fallback if git is unavailable or the
conversion has not been done:

    bash tools/pull.sh main --dry
    bash tools/pull.sh main
    python tools/manifest.py --check

## Applying an archive by hand

    python tools/apply_sync.py --dry
    python tools/apply_sync.py

It finds the project root by looking for `config.py` beside `src/`, so an
archive with a wrapper directory cannot nest the package inside itself; it
writes file contents in place (for mounts that forbid unlinking); and it never
touches `data/`.

## The silent failures this procedure prevents

- `scp -r src host:.../src` when `src` exists → `src/src`. Copy to the parent.
- An archive with `NWP_Deployment_Package/` at its top, extracted from inside
  that directory → the package nested in itself.
- curl-based pulls only add and overwrite — a file deleted upstream survives
  forever. `manifest.py --check` reports it as EXTRA.
- Work never pushed at all, so a correct pull delivers a stale tree.
- `git clean -x` or `-X` on the server deletes `data/` (and `git stash --all`
  removes it from the tree): it is ignored, and ignored is not protected.
- On the server the archive is inside `NWP_Deployment_Package/`, which looks
  like a nesting accident. Tidying it away deletes the archive.

None of these raises an error. Always finish with `tools/manifest.py --check`.

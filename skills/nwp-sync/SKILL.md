---
name: nwp-sync
description: Use when moving NWP1 work between places — the GitHub repository, the Windows desktop clone, the shared compute server, or the claude.ai project docs — or when files are missing, nested wrongly, or out of date after a transfer.
---

# Syncing NWP1

## The places, and which is the source of truth

| where | holds | reached by |
|---|---|---|
| GitHub `SilentSquid236/NWP1` | code, `docs/`, `tools/`, `skills/` — **source of truth** | git on Windows; `tools/pull.sh` on the server |
| Windows desktop `Desktop\NWP\NWP_Deployment_Package` | the clone with push credentials | git |
| shared server (data5 / the Xeon) | the working copy that runs, plus `data/` | `tools/pull.sh` (curl; **no git there**) |

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

None of these raises an error. Always finish with `tools/manifest.py --check`.

#!/usr/bin/env bash
#
# Pull the repository onto a machine with no git, using curl.
#
#   tools/pull.sh              # main branch, into this project directory
#   tools/pull.sh dev          # another branch
#   tools/pull.sh main --dry   # show what would change, touch nothing
#
# WHY CURL AND NOT GIT
#
# git is not installed on the server and cannot be (no sudo, and the admin
# policy forbids installs). GitHub serves any branch as a tarball over plain
# HTTPS, which needs nothing but curl and tar.
#
# WHAT THIS GETS RIGHT THAT A BARE curl | tar DOES NOT
#
#   * STAGING. The tarball is extracted to a temporary directory and inspected
#     BEFORE anything in the project is touched. A half-downloaded archive
#     then costs nothing.
#   * THE EXTRA TOP LEVEL. GitHub wraps everything in `<repo>-<branch>/`, and
#     the project may sit inside that again. Both are detected rather than
#     assumed -- guessing wrong is what produces the nested `src/src` that
#     took a human to spot (P-24 in docs/PROBLEMS.md).
#   * data/ IS NEVER TOUCHED. It holds the verification archive, which is the
#     one thing in this project that cannot be recreated (P-07).
#   * BANDWIDTH. Rate-limited, because this is a shared link.
#   * A LAYOUT CHECK afterwards, so a bad copy is caught now rather than as
#     confusing behaviour later.

set -uo pipefail

REPO="${NWP_REPO:-SilentSquid236/NWP1}"
BRANCH="${1:-main}"
DRY=0
[ "${2:-}" = "--dry" ] && DRY=1

DEST="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RATE="${NWP_MAX_MBPS:-8}M"
URL="https://codeload.github.com/$REPO/tar.gz/refs/heads/$BRANCH"

STAGE="$(mktemp -d "${TMPDIR:-/tmp}/nwppull.XXXXXX")"
trap 'rm -rf "$STAGE"' EXIT INT TERM

echo "repo    $REPO ($BRANCH)"
echo "dest    $DEST"
echo "rate    $RATE max"
echo

echo "downloading..."
if ! curl -fsSL --limit-rate "$RATE" --retry 3 --retry-delay 5 \
        -o "$STAGE/repo.tar.gz" "$URL"; then
    echo "FAILED to download $URL"
    echo "  * check the branch name exists"
    echo "  * a private repo needs a token:"
    echo "      curl -H \"Authorization: Bearer \$GITHUB_TOKEN\" -fsSL ... "
    exit 1
fi
echo "  $(du -h "$STAGE/repo.tar.gz" | cut -f1) downloaded"

tar -xzf "$STAGE/repo.tar.gz" -C "$STAGE" || { echo "FAILED to extract"; exit 1; }

# GitHub wraps everything in <repo>-<branch>/. Find the real project root by
# looking for config.py rather than assuming a depth.
SRC=""
for cand in "$STAGE"/*/ "$STAGE"/*/*/; do
    [ -f "${cand}config.py" ] && [ -d "${cand}src" ] && { SRC="${cand%/}"; break; }
done
if [ -z "$SRC" ]; then
    echo "FAILED: no directory in the archive contains both config.py and src/"
    echo "        archive top level:"
    ls -1 "$STAGE" | sed 's/^/          /'
    exit 1
fi
echo "  project root in archive: ${SRC#$STAGE/}"

echo
echo "changes:"
CHANGED=0
while IFS= read -r rel; do
    case "$rel" in
        data/*|.git/*|*/__pycache__/*|*.pyc) continue ;;
    esac
    if [ ! -f "$DEST/$rel" ]; then
        echo "  new      $rel"; CHANGED=$((CHANGED+1))
    elif ! cmp -s "$SRC/$rel" "$DEST/$rel"; then
        echo "  updated  $rel"; CHANGED=$((CHANGED+1))
    fi
done < <(cd "$SRC" && find . -type f | sed 's|^\./||' | sort)

if [ "$CHANGED" -eq 0 ]; then
    echo "  none -- already up to date"
    exit 0
fi
echo "  ($CHANGED file(s))"

if [ "$DRY" -eq 1 ]; then
    echo
    echo "--dry: nothing written"
    exit 0
fi

echo
echo "copying..."
# Copy the CONTENTS of the staged root into the destination. The trailing /.
# is what stops src/ landing inside src/.
cp -a "$SRC/." "$DEST/" || { echo "FAILED to copy"; exit 1; }
echo "  done ($CHANGED file(s)); data/ untouched"

echo
if [ -f "$DEST/tools/checklayout.py" ]; then
    python "$DEST/tools/checklayout.py"
else
    echo "note: tools/checklayout.py not present, skipping layout check"
fi

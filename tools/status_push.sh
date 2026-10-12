#!/bin/bash
# Rebuild the NWP1 status page and push it to the 'status' branch on GitHub.
# Cron, every 10 min (tagged # NWP1-cron). Installed copy: data/status/bin/.
#
# The branch holds one commit (amended each time), so it never grows. The
# deploy key (~/.ssh/nwp1_status_deploy) works for this one repository only;
# delete the key on GitHub (Settings -> Deploy keys) to stop the pushes.
set -u
DATA=${NWP_DATA_ROOT:-/data5/pierce/AINWP/data}
B=$DATA/status
KEY=$HOME/.ssh/nwp1_status_deploy
REPO=$B/repo
REMOTE=git@github.com:SilentSquid236/NWP1.git
export GIT_SSH_COMMAND="ssh -i $KEY -o IdentitiesOnly=yes -o BatchMode=yes -o ConnectTimeout=20"

mkdir -p $B
python3 $B/bin/status_page.py "$DATA" > $B/STATUS.md.tmp && mv -f $B/STATUS.md.tmp $B/STATUS.md || exit 1
[ -f $KEY ] || exit 0

if [ ! -d $REPO/.git ]; then
  mkdir -p $REPO && cd $REPO && git init -q && git checkout -q --orphan status
  git remote add origin $REMOTE
  git config user.name "NWP1 status page"
  git config user.email "nwp1-status@localhost"
fi
cd $REPO || exit 1
cp -f $B/STATUS.md README.md
git add README.md
if git rev-parse -q --verify HEAD > /dev/null; then
  git commit -q --amend -m "status $(date -u +%FT%TZ)"
else
  git commit -q -m "status $(date -u +%FT%TZ)"
fi
timeout 90 git push -q -f origin status 2>&1 | grep -v '^$' | tail -2

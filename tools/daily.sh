#!/usr/bin/env bash
#
# One day of the archive: ingest, forecast, verify.
#
#   tools/daily.sh                 # today's 00Z run
#   tools/daily.sh 2026-09-04T00   # a specific run
#
# Written to be run from cron on the shared server. Cron gives you a bare
# environment, no terminal, and no memory of the last run, so:
#
#   * every path is absolute, derived from this script's own location
#   * a LOCK FILE stops two runs overlapping. A 12-hour forecast takes longer
#     than the gap between some cron schedules, and two copies competing for
#     the same cores on a shared machine is exactly what the 50% ceiling in
#     resources.py exists to avoid.
#   * output goes to a dated log, kept, because a failure at 4 a.m. is only
#     diagnosable from what it wrote at the time
#   * python runs with -u. Without it Python BLOCK-BUFFERS stdout whenever it
#     is not a terminal, so a redirected log stays empty for many minutes and
#     a running job is indistinguishable from a frozen one. That is not a
#     cosmetic detail: it is how the first real run got reported as a freeze.
#   * the exit code is the FIRST failure, not the last command's
#
# WHY THIS RUNS EVERY DAY AND THE OTHER SCRIPTS DO NOT (P-07)
#
# Verification needs forecasts and observations paired in time. A day not
# archived cannot be recovered later: the observations remain downloadable,
# but the forecast that was valid for them was never made. Missing a day
# costs a day of evidence permanently.
#
# Suggested crontab entry (03:30 local, after the 00Z HRRR is complete):
#
#   30 3 * * *  /path/to/NWP_Deployment_Package/tools/daily.sh >/dev/null 2>&1
#
# Nothing here needs root, and nothing installs anything.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN="${1:-$(date -u +%Y-%m-%dT00)}"
STAMP="$(echo "$RUN" | tr -d ':-' | tr 'T' '_')"

DATA="${NWP_DATA_ROOT:-$ROOT/data}"
LOGDIR="$DATA/logs"
LOCK="$DATA/daily.lock"
LOG="$LOGDIR/daily_$STAMP.log"
RUNDIR="$DATA/tensors/analysis_$STAMP"

mkdir -p "$LOGDIR" "$DATA"

# Lock. `set -o noclobber` makes this atomic without needing flock, which is
# not present everywhere.
if ! (set -o noclobber; echo "$$ $(date -u +%FT%TZ)" > "$LOCK") 2>/dev/null; then
    echo "another run is active (lock: $LOCK, holder: $(cat "$LOCK" 2>/dev/null))"
    exit 75          # EX_TEMPFAIL: try again later, do not alarm
fi
trap 'rm -f "$LOCK"' EXIT INT TERM

# How to see where a step is, while it runs:
#     tail -f "$LOG"                     # what it has printed
#     kill -USR1 $(pgrep -f forecast.py) # traceback, run continues
#
STATUS=0
step() {
    local name="$1"; shift
    echo "=== $name  $(date -u +%FT%TZ)" >> "$LOG"
    if "$@" >> "$LOG" 2>&1; then
        echo "    ok" >> "$LOG"
    else
        local rc=$?
        echo "    FAILED rc=$rc" >> "$LOG"
        [ "$STATUS" -eq 0 ] && STATUS=$rc
        return $rc
    fi
}

{
    echo "daily archive run for $RUN"
    echo "root   $ROOT"
    echo "data   $DATA"
} >> "$LOG"

# 1. Ingest. --hours 13 gives a 12-hour forecast one boundary frame per hour
#    plus the initial state.
step ingest python -u "$ROOT/src/ingest_hrrr.py" \
    --start "$RUN" --hours 13 --stride 4 || true

# 2. Forecast. Runs even if some hours are missing -- a shorter forecast is
#    still worth archiving; no forecast at all is not.
step forecast python -u "$ROOT/src/forecast.py" \
    --run-dir "$RUNDIR" --hours 12 --output-every 1 || true

# 3. Verify and archive. THIS IS THE STEP THAT CANNOT BE DEFERRED, so it is
#    attempted even if the forecast step reported a failure: a forecast that
#    diverged at hour 8 still produced eight hours worth archiving.
if [ -f "$RUNDIR/forecast.npz" ]; then
    step verify python -u "$ROOT/src/verify.py" \
        --forecast "$RUNDIR/forecast.npz" --run-time "$RUN" || true
else
    echo "=== verify  SKIPPED: no forecast.npz produced" >> "$LOG"
    [ "$STATUS" -eq 0 ] && STATUS=1
fi

echo "=== done  $(date -u +%FT%TZ)  status=$STATUS" >> "$LOG"
exit "$STATUS"

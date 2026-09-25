#!/usr/bin/env bash
#
# One forecast cycle: observations -> analysis -> forecast -> archive.
# And, separately, verification of the cycles whose windows have closed.
#
#   bash tools/daily.sh                         # latest cycle, 24 h, observations
#   bash tools/daily.sh 2026-09-21T12           # a specific cycle (00/06/12/18Z)
#   bash tools/daily.sh 2026-09-21T12 12        # ... 12 hours long
#   bash tools/daily.sh 2026-09-21T12 24 hrrr   # the old HRRR-seeded baseline
#   bash tools/daily.sh verify                  # verify every closed window
#
# It is a BASH script: `python tools/daily.sh` fails at the first line of
# shell (prompt 92).
#
# WHAT CHANGED ON 2026-09-22 (prompts 94-101, P-53, P-55)
#
#   * The run starts from OBSERVATIONS at its own cycle time (src/ingest_obs.py)
#     and holds its edges to that analysis. No HRRR unless asked for.
#   * Four cycles a day, 12-24 h each, and each run inside 1.5 h of wall
#     clock: the forecast gets a deadline computed from what is left of the
#     budget after ingest, and writes the hours it reached.
#   * Verification is NOT part of a run -- the observations it needs do not
#     exist yet. `daily.sh verify` scores every run whose window has closed.
#   * Maps (2026-09-25): after the forecast, src/make_maps.py draws every
#     product for the hours reached into <rundir>/maps/, with index.html as
#     the viewer; `verify` adds the error maps. A map failure never changes
#     the cycle's status.
#   * One variable, RUNDIR, gives every step its directory. The old script
#     handed the forecast $DATA/tensors/... while ingest wrote to
#     $DATA/tensors_3d/... (P-55).
#
# Suggested crontab (UTC; soundings reach IEM a little after nominal time):
#
#   45 1,7,13,19 * * *  bash /data5/pierce/AINWP/tools/daily.sh >/dev/null 2>&1
#   30 3 * * *          bash /data5/pierce/AINWP/tools/daily.sh verify >/dev/null 2>&1
#
# Cron does not read ~/.bashrc, so NWP_DATA_ROOT is unset there and DATA falls
# back to $ROOT/data -- which is /data5/pierce/AINWP/data, the same place.
#
# Every path is absolute and derived from this script's location, a lock
# stops two jobs competing for cores, output goes to a dated log that is kept,
# python runs unbuffered (-u), and the exit code is the FIRST failure.
# Nothing here needs root, and nothing installs anything.

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="${NWP_DATA_ROOT:-$ROOT/data}"
TENSORS="$DATA/tensors_3d"          # == config.TENSOR_DIR; the ONE place runs live
LOGDIR="$DATA/logs"
LOCK="$DATA/daily.lock"
BUDGET_MIN="${NWP_CYCLE_BUDGET_MIN:-90}"
RESERVE_MIN=8                       # after the forecast: writing, archiving, maps

mkdir -p "$LOGDIR" "$DATA"

latest_cycle() {
    # The newest 00/06/12/18Z cycle at least 100 min old, so its soundings
    # and hourly reports have had time to reach the archive.
    local t h
    t="$(date -u -d '-100 min' +%Y-%m-%dT%H)"
    h="${t:11:2}"
    printf '%sT%02d' "${t:0:10}" $(( (10#$h / 6) * 6 ))
}

MODE="cycle"
if [ "${1:-}" = "verify" ]; then
    MODE="verify"
    RUN="$(date -u +%Y-%m-%dT%H)"
else
    RUN="${1:-$(latest_cycle)}"
fi
HOURS="${2:-24}"
SOURCE="${3:-obs}"
STAMP="$(echo "$RUN" | tr -d ':-' | tr 'T' '_')"
LOG="$LOGDIR/${MODE}_$STAMP.log"

if ! (set -o noclobber; echo "$$ $MODE $RUN $(date -u +%FT%TZ)" > "$LOCK") 2>/dev/null; then
    echo "another job is active (lock: $LOCK, holder: $(cat "$LOCK" 2>/dev/null))"
    exit 75          # EX_TEMPFAIL: try again later, do not alarm
fi
trap 'rm -f "$LOCK"' EXIT INT TERM

STATUS=0
T0=$(date +%s)
step() {
    local name="$1"; shift
    echo "=== $name  $(date -u +%FT%TZ)  (+$(( ($(date +%s) - T0) / 60 )) min)" >> "$LOG"
    "$@" >> "$LOG" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then
        echo "    ok" >> "$LOG"
    else
        echo "    FAILED rc=$rc" >> "$LOG"
        [ "$STATUS" -eq 0 ] && STATUS=$rc
    fi
    return $rc
}

{
    echo "NWP1 $MODE for $RUN"
    echo "root   $ROOT"
    echo "data   $DATA"
} >> "$LOG"

maps() {
    # Maps never change a cycle's status: a failed map is not a failed forecast.
    local keep=$STATUS
    step maps python -u "$ROOT/src/make_maps.py" "$@" || true
    STATUS=$keep
}

if [ "$MODE" = "verify" ]; then
    step verify python -u "$ROOT/src/verify_pending.py"
    # Error maps for every run verified since its maps were last drawn.
    for d in "$TENSORS"/obs_*; do
        [ -f "$d/verified.json" ] || continue
        [ -f "$d/maps/index.html" ] && [ "$d/maps/index.html" -nt "$d/verified.json" ] && continue
        maps --run-dir "$d" --products err_t,err_w
    done
    echo "=== done  $(date -u +%FT%TZ)  status=$STATUS" >> "$LOG"
    exit "$STATUS"
fi

echo "cycle  $RUN, $HOURS h, source $SOURCE, budget $BUDGET_MIN min" >> "$LOG"

# 1. Initial state.
if [ "$SOURCE" = "hrrr" ]; then
    RUNDIR="$TENSORS/analysis_$STAMP"
    step ingest python -u "$ROOT/src/ingest_hrrr.py" \
        --start "$RUN" --hours "$(( HOURS + 1 ))" --stride 4 || true
else
    RUNDIR="$TENSORS/obs_$STAMP"
    step ingest python -u "$ROOT/src/ingest_obs.py" --cycle "$RUN" || true
fi
echo "rundir $RUNDIR" >> "$LOG"

# 2. Forecast, inside what is left of the budget.
USED_MIN=$(( ($(date +%s) - T0) / 60 ))
LEFT_MIN=$(( BUDGET_MIN - USED_MIN - RESERVE_MIN ))
if [ "$LEFT_MIN" -lt 5 ]; then
    echo "=== forecast  SKIPPED: ingest used $USED_MIN of $BUDGET_MIN min" >> "$LOG"
    [ "$STATUS" -eq 0 ] && STATUS=1
else
    # rc 2 = deadline cut the run short, 3 = diverged; both still write the
    # hours reached, and both count as the run's status.
    step forecast python -u "$ROOT/src/forecast.py" \
        --run-dir "$RUNDIR" --hours "$HOURS" --output-every 1 \
        --deadline-min "$LEFT_MIN" || true
fi

if [ ! -f "$RUNDIR/forecast.npz" ]; then
    echo "=== no forecast.npz in $RUNDIR -- nothing to verify later" >> "$LOG"
    [ "$STATUS" -eq 0 ] && STATUS=1
fi

# 3. Maps and the viewer (<rundir>/maps/index.html), for the hours reached.
if [ -f "$RUNDIR/forecast.npz" ] || [ -f "$RUNDIR/obs_analysis_f00.npz" ]; then
    maps --run-dir "$RUNDIR"
fi

echo "=== done  $(date -u +%FT%TZ)  +$(( ($(date +%s) - T0) / 60 )) min  status=$STATUS" >> "$LOG"
exit "$STATUS"

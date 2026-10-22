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
#   10 5,11,17,23 * * * bash /data5/pierce/AINWP/tools/daily.sh verify >/dev/null 2>&1
#     (four times a day, in the gaps between the :45 cam_cycle runs, so each
#     run is scored ~5 h after its window closes; was 30 3 * * * until 2026-10-19)
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
# Array backend for the core (prompt 128). "torch" runs the same float64
# physics on THREADS cores, ~8x faster; switch the default only after
# tools/check_backend.py shows round-off agreement on this machine.
BACKEND="${NWP_BACKEND:-numpy}"
THREADS="${NWP_THREADS:-8}"
# P-82, test J1 (2026-10-15, user decision "adopt A3, then J2"): upper-air
# analysis options. NWP_UPPER_DOMAIN_ONLY=1 leaves soundings outside the domain
# out of the upper air (their innovations were taken against edge values);
# NWP_UPPER_L_KM sets its length scale (build.py default 450 km). Unset = as before.
UPPER_ARGS=""
[ "${NWP_UPPER_DOMAIN_ONLY:-0}" = 1 ] && UPPER_ARGS="--upper-domain-only"
[ -n "${NWP_UPPER_L_KM:-}" ] && UPPER_ARGS="$UPPER_ARGS --upper-L-km $NWP_UPPER_L_KM"
# NWP_INIT_BALANCE (P-83, tests J2-J4, 2026-10-17): the start-up divergence removal
# of the 12 km forecast (forecast.py --init-balance; unset = its default, exact).
# exact-free (J4) holds no boundary face fixed. tools/cam_cycle.sh reads the same
# switch for its 12 km moist run, and NWP_CAM_INIT_BALANCE for the 3 km run.
INIT_ARGS=""
[ -n "${NWP_INIT_BALANCE:-}" ] && INIT_ARGS="--init-balance $NWP_INIT_BALANCE"
# NWP_12KM_MOIST=1 (2026-10-20, prompt 232, user decision "promote the moist NH
# run"): the 12 km PRODUCT is the moist non-hydrostatic forecast forecast_nhm.npz
# (40 levels, warm rain and, with NWP_CAM_ICE=1, ice; the file cam_cycle.sh used
# to make only for the 3 km edges). It runs first, then the dry hydrostatic
# forecast.npz, which stays as the cycling first guess (ingest_obs.previous_run)
# until a test switches that too. Maps, the viewer and verification use the
# product file; verify_pending.py reads NWP_12KM_PRODUCT. Unset = as before.
MOIST12="${NWP_12KM_MOIST:-0}"
NH_THREADS="${NWP_CAM_THREADS:-26}"
ICE_ARG=""; [ "${NWP_CAM_ICE:-0}" = 1 ] && ICE_ARG="--ice"
if [ "$MOIST12" = 1 ]; then export NWP_12KM_PRODUCT="${NWP_12KM_PRODUCT:-forecast_nhm.npz}"; fi
product_of() {   # the product file of a run dir: the moist run when promoted and present
    if [ "$MOIST12" = 1 ] && [ -f "$1/forecast_nhm.npz" ]; then echo forecast_nhm.npz; else echo forecast.npz; fi
}

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
        maps --run-dir "$d" --forecast "$(product_of "$d")" --products err_t,err_w
    done
    # Skill against persistence from the live-observation archives, for the
    # status page (data/status/scoreboard.md) and data/verification_scores.csv.
    # Like the maps, a scoreboard failure never changes the verify status.
    keep=$STATUS
    step scoreboard python -u "$ROOT/tools/scoreboard.py" "$DATA" || true
    STATUS=$keep
    echo "=== done  $(date -u +%FT%TZ)  status=$STATUS" >> "$LOG"
    exit "$STATUS"
fi

echo "cycle  $RUN, $HOURS h, source $SOURCE, budget $BUDGET_MIN min, backend $BACKEND x$THREADS" >> "$LOG"

# 1. Initial state.
if [ "$SOURCE" = "hrrr" ]; then
    RUNDIR="$TENSORS/analysis_$STAMP"
    step ingest python -u "$ROOT/src/ingest_hrrr.py" \
        --start "$RUN" --hours "$(( HOURS + 1 ))" --stride 4 || true
else
    RUNDIR="$TENSORS/obs_$STAMP"
    step ingest python -u "$ROOT/src/ingest_obs.py" --cycle "$RUN" $UPPER_ARGS || true
fi
echo "rundir $RUNDIR" >> "$LOG"

# 2a. The moist NH 12 km forecast: the 12 km product when NWP_12KM_MOIST=1
# (same command as cam_cycle.sh step 1b, which then reuses this file).
if [ "$MOIST12" = 1 ] && [ -f "$RUNDIR/obs_analysis_f00.npz" ]; then
    rm -f "$RUNDIR/forecast_nhm.npz"     # never keep a stale file from an earlier attempt
    step "moist 12 km forecast (product)" python -u "$ROOT/src/forecast.py" --run-dir "$RUNDIR" \
        --hours "$HOURS" --output-every 1 --core nh --nh-backend c --advection upwind3 \
        --nh-physics step --levels 40 --threads "$NH_THREADS" --moist $ICE_ARG \
        --out "$RUNDIR/forecast_nhm.npz" $INIT_ARGS --deadline-min 30 || true
fi

# 2b. The dry hydrostatic forecast (the product when NWP_12KM_MOIST is unset;
# otherwise the next cycle's first guess), inside what is left of the budget.
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
        --backend "$BACKEND" --threads "$THREADS" $INIT_ARGS \
        --deadline-min "$LEFT_MIN" || true
fi

if [ ! -f "$RUNDIR/forecast.npz" ]; then
    echo "=== no forecast.npz in $RUNDIR -- nothing to verify later" >> "$LOG"
    [ "$STATUS" -eq 0 ] && STATUS=1
fi

# 3. Maps and the viewer (<rundir>/maps/index.html), for the hours reached.
# NWP_DEFER_MAPS=1 (set by cam_cycle.sh with NWP_CAM_MAPS=1, 2026-10-20): the
# caller draws them on spare CPUs so the 3 km run does not wait for them.
if [ "${NWP_DEFER_MAPS:-0}" = 1 ]; then
    echo "=== maps deferred to the caller (NWP_DEFER_MAPS=1)" >> "$LOG"
elif [ -f "$RUNDIR/forecast.npz" ] || [ -f "$RUNDIR/obs_analysis_f00.npz" ]; then
    maps --run-dir "$RUNDIR" --forecast "$(product_of "$RUNDIR")"
fi

# Compute ledger (docs/COMPUTE_LEDGER.md): one line per run, kept with the logs.
log_compute() {   # script, threads
    local f="$LOGDIR/compute_hours.csv" wall
    wall=$(( $(date +%s) - T0 ))
    [ -f "$f" ] || echo "start_utc,cycle,script,wall_min,threads,core_h,status" > "$f"
    echo "$(date -u -d @"$T0" +%FT%TZ),$RUN,$1,$(awk "BEGIN{printf \"%.1f\", $wall/60}"),$2,$(awk "BEGIN{printf \"%.2f\", $wall*$2/3600}"),$STATUS" >> "$f"
}
log_compute daily.sh "$THREADS"
echo "=== done  $(date -u +%FT%TZ)  +$(( ($(date +%s) - T0) / 60 )) min  status=$STATUS" >> "$LOG"
exit "$STATUS"

#!/usr/bin/env bash
#
# One CAM cycle (CAM stage S5): the 12 km cycle, then the 3 km nest inside it.
#
#   bash tools/cam_cycle.sh                     # latest cycle, 24 h
#   bash tools/cam_cycle.sh 2026-09-28T06       # a specific cycle
#   bash tools/cam_cycle.sh 2026-09-28T06 12    # ... 12 hours long
#   NWP_SKIP_12KM=1 bash tools/cam_cycle.sh 2026-09-28T06
#                                               # the 12 km run already exists
#
# What it does, in order, all inside ONE budget (NWP_CYCLE_BUDGET_MIN, 90):
#   1. tools/daily.sh for the cycle: the 12 km analysis and forecast, as in
#      production (its own lock and log).
#   2. The 3 km analysis, from the observations the 12 km cycle archived,
#      starting from the 12 km analysis regridded, with the terrain- and
#      coast-aware gap-filling of the surface stations (ingest_obs.py
#      --spacing-km 3 --gapfill; test S5b). It goes to $DATA/tensors_3km.
#   3. The 3 km forecast: the non-hydrostatic core in C (--core nh
#      --nh-backend c), 40 levels, physics once per step, edges from the 12 km
#      forecast every hour (--boundary-forecast; nest.py). Its deadline is
#      what is left of the budget.
#
# NWP_CAM_THREADS (default 26 = one socket) sets the OpenMP threads of the
# 3 km run. The S3b/S4b measurements: 26 -> 52 threads gains only 1.2x.
# Nothing here installs anything; the C kernels are compiled by gcc at first
# use and cached in $NWP_CBUILD (default $DATA/cbuild).

set -uo pipefail

ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
DATA="${NWP_DATA_ROOT:-$ROOT/data}"
T12="$DATA/tensors_3d"
T3="$DATA/tensors_3km"
LOGDIR="$DATA/logs"
LOCK="${NWP_CAM_LOCK:-$DATA/cam.lock}"   # a campaign running two cycles at once gives each its own
BUDGET_MIN="${NWP_CYCLE_BUDGET_MIN:-90}"
RESERVE_MIN=6   # forecast set-up (~1 min) + output write (~1.3 min since S5e writes
                # float32 at zlib 1; it was ~6 min at float64/zlib 6 in S5c2), with margin.
                # The deadline counts only the steps.
CAM_THREADS="${NWP_CAM_THREADS:-26}"
export NWP_CBUILD="${NWP_CBUILD:-$DATA/cbuild}"
# resources.py sets the BLAS/OpenMP caps of every Python step to a fraction of
# the machine (default 50 % = 52 threads). A quarter (26, one socket) for the
# CAM steps: test S5b's parallel analyses at the default took the load to 94.
export NWP_RESOURCE_FRACTION="${NWP_RESOURCE_FRACTION:-0.25}"
# OpenMP threads of the C kernels: one per physical core, consecutive core
# numbers (CAM S5d). On met-23030149 the NUMA nodes are the even and the odd
# CPU numbers, so 26 "close" threads alternate between the two sockets and use
# both memory controllers. Test S5f: 0.727 s/step this way against 0.808 s
# pinned to one socket. Two cycles at once must each be pinned to a socket
# (taskset -c 0,2,..,50 and 1,3,..,51): 0.81-0.82 s/step each, against 1.25
# when both use the default placement and land on the same cores.
export OMP_PLACES="${OMP_PLACES:-cores}" OMP_PROC_BIND="${OMP_PROC_BIND:-close}"
mkdir -p "$LOGDIR" "$T3" "$NWP_CBUILD"

latest_cycle() {          # as in daily.sh
    local t h
    t="$(date -u -d '-100 min' +%Y-%m-%dT%H)"
    h="${t:11:2}"
    printf '%sT%02d' "${t:0:10}" $(( (10#$h / 6) * 6 ))
}

RUN="${1:-$(latest_cycle)}"
HOURS="${2:-24}"
STAMP="$(echo "$RUN" | tr -d ':-' | tr 'T' '_')"
D12="$T12/obs_$STAMP"
D3="$T3/obs_$STAMP"
LOG="$LOGDIR/cam_$STAMP.log"

if ! (set -o noclobber; echo "$$ cam $RUN $(date -u +%FT%TZ)" > "$LOCK") 2>/dev/null; then
    echo "another CAM job is active (lock: $LOCK, holder: $(cat "$LOCK" 2>/dev/null))"
    exit 75
fi
trap 'rm -f "$LOCK"' EXIT INT TERM

STATUS=0
T0=$(date +%s)
step() {
    local name="$1"; shift
    echo "=== $name  $(date -u +%FT%TZ)  (+$(( ($(date +%s) - T0) / 60 )) min)" >> "$LOG"
    "$@" >> "$LOG" 2>&1
    local rc=$?
    if [ $rc -eq 0 ]; then echo "    ok" >> "$LOG"; else
        echo "    FAILED rc=$rc" >> "$LOG"; [ "$STATUS" -eq 0 ] && STATUS=$rc; fi
    return $rc
}

echo "NWP1 CAM cycle $RUN, $HOURS h, budget $BUDGET_MIN min, $CAM_THREADS threads" >> "$LOG"

# 1. The 12 km cycle
if [ -z "${NWP_SKIP_12KM:-}" ]; then
    step "12 km cycle (daily.sh)" bash "$ROOT/tools/daily.sh" "$RUN" "$HOURS" || true
fi
for f in "$D12/obs_analysis_f00.npz" "$D12/forecast.npz"; do
    if [ ! -f "$f" ]; then
        echo "=== no $f: the 3 km run needs the 12 km analysis and forecast" >> "$LOG"
        exit 1
    fi
done

# 2. The 3 km analysis
step "3 km analysis" python -u "$ROOT/src/ingest_obs.py" --cycle "$RUN" --from-raw \
    --raw-from "$D12" --out-root "$T3" --spacing-km 3 \
    --background-analysis "$D12/obs_analysis_f00.npz" --gapfill || exit "$STATUS"

# 3. The 3 km forecast, inside what is left of the budget
USED_MIN=$(( ($(date +%s) - T0) / 60 ))
LEFT_MIN=$(( BUDGET_MIN - USED_MIN - RESERVE_MIN ))
if [ "$LEFT_MIN" -lt 5 ]; then
    echo "=== 3 km forecast SKIPPED: $USED_MIN of $BUDGET_MIN min already used" >> "$LOG"
    exit 1
fi
step "3 km forecast" python -u "$ROOT/src/forecast.py" --run-dir "$D3" --hours "$HOURS" \
    --output-every 1 --core nh --nh-backend c --threads "$CAM_THREADS" \
    --nh-physics step --levels 40 --boundary-forecast "$D12/forecast.npz" \
    --deadline-min "$LEFT_MIN" || true

# Compute ledger (docs/COMPUTE_LEDGER.md): one line per run, kept with the logs.
log_compute() {   # script, threads
    local f="$LOGDIR/compute_hours.csv" wall
    wall=$(( $(date +%s) - T0 ))
    [ -f "$f" ] || echo "start_utc,cycle,script,wall_min,threads,core_h,status" > "$f"
    echo "$(date -u -d @"$T0" +%FT%TZ),$RUN,$1,$(awk "BEGIN{printf \"%.1f\", $wall/60}"),$2,$(awk "BEGIN{printf \"%.2f\", $wall*$2/3600}"),$STATUS" >> "$f"
}
log_compute cam_cycle.sh "$CAM_THREADS"
echo "=== done  $(date -u +%FT%TZ)  (+$(( ($(date +%s) - T0) / 60 )) min)  status=$STATUS" >> "$LOG"
exit "$STATUS"

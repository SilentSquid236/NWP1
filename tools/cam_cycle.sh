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
# NWP_CAM_MOIST=1 (E7, 2026-10-12): the 3 km run carries water vapour, cloud
# and rain (--moist), with its edges from a moist 12 km NH forecast run here
# (E4's nhMb set-up), because the production 12 km forecast has no qv/qc/qr. The
# moist 3 km takes ~96 min per 24 h, so the default budget rises to 150 min.
CAM_MOIST="${NWP_CAM_MOIST:-0}"
# NWP_CAM_ICE=1 (P-81, test R2, 2026-10-14): simple ice (--ice) in both moist runs.
CAM_ICE_ARG=""; [ "${NWP_CAM_ICE:-0}" = 1 ] && CAM_ICE_ARG="--ice"
# NWP_UPPER_DOMAIN_ONLY=1, NWP_UPPER_L_KM (P-82, test J1, 2026-10-15): the same
# upper-air analysis options as tools/daily.sh, for the 3 km analysis.
UPPER_ARGS=""
[ "${NWP_UPPER_DOMAIN_ONLY:-0}" = 1 ] && UPPER_ARGS="--upper-domain-only"
[ -n "${NWP_UPPER_L_KM:-}" ] && UPPER_ARGS="$UPPER_ARGS --upper-L-km $NWP_UPPER_L_KM"
# NWP_INIT_BALANCE / NWP_CAM_INIT_BALANCE (P-83, tests J2-J4, 2026-10-17): start-up
# divergence removal (forecast.py --init-balance) of the 12 km moist run and of the
# 3 km run. Unset = forecast.py's default (exact).
INIT12_ARGS=""; [ -n "${NWP_INIT_BALANCE:-}" ] && INIT12_ARGS="--init-balance $NWP_INIT_BALANCE"
INIT3_ARGS=""; [ -n "${NWP_CAM_INIT_BALANCE:-}" ] && INIT3_ARGS="--init-balance $NWP_CAM_INIT_BALANCE"
if [ "$CAM_MOIST" = 1 ]; then BUDGET_MIN="${NWP_CYCLE_BUDGET_MIN:-150}"
else BUDGET_MIN="${NWP_CYCLE_BUDGET_MIN:-90}"; fi
RESERVE_MIN=6   # forecast set-up (~1 min) + output write (~1.3 min since S5e writes
                # float32 at zlib 1; it was ~6 min at float64/zlib 6 in S5c2), with margin.
                # The deadline counts only the steps.
CAM_THREADS="${NWP_CAM_THREADS:-26}"
# NWP_CAM_ARGS: extra forecast.py options for the 3 km run, e.g. "--moist"
# (CAM stage S7); empty by default.
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
# Maps (2026-10-20, prompt 232): NWP_CAM_MAPS=1 draws the 12 km maps in the
# background as soon as the 12 km product exists, and the 3 km maps (with
# click-to-sounding data) after the 3 km forecast, both niced and pinned to
# NWP_MAP_CPUS, outside the "close" cores of the 26-thread 3 km run (0-25).
# A map failure never changes the cycle status. Unset = daily.sh draws the
# 12 km maps itself before the 3 km run, and there are no 3 km maps (as before).
CAM_MAPS="${NWP_CAM_MAPS:-0}"
MAP_CPUS="${NWP_MAP_CPUS:-50,51}"
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

echo "NWP1 CAM cycle $RUN, $HOURS h, budget $BUDGET_MIN min, $CAM_THREADS threads$([ "$CAM_MOIST" = 1 ] && echo ', moist')" >> "$LOG"

# 1. The 12 km cycle
if [ -z "${NWP_SKIP_12KM:-}" ]; then
    export NWP_CAM_THREADS="$CAM_THREADS"   # daily.sh: threads of its moist 12 km run
    [ "$CAM_MAPS" = 1 ] && export NWP_DEFER_MAPS=1
    step "12 km cycle (daily.sh)" bash "$ROOT/tools/daily.sh" "$RUN" "$HOURS" || true
fi
# Edges: the 12 km forecast of the same cycle; NWP_CAM_BOUNDARY overrides the
# file (CAM stage S6 campaign: the AM forecasts kept by test AM).
BFC="${NWP_CAM_BOUNDARY:-$D12/forecast.npz}"
for f in "$D12/obs_analysis_f00.npz" "$BFC"; do
    if [ ! -f "$f" ]; then
        echo "=== no $f: the 3 km run needs the 12 km analysis and forecast" >> "$LOG"
        exit 1
    fi
done

# 1a. The 12 km maps, in the background (NWP_CAM_MAPS=1).
draw_maps() {   # run dir, forecast file name
    nice -n 19 taskset -c "$MAP_CPUS" env OMP_NUM_THREADS=2 MKL_NUM_THREADS=2 OPENBLAS_NUM_THREADS=2 \
        NWP_RESOURCE_FRACTION=0.02 timeout 3600 python -u "$ROOT/src/make_maps.py" --run-dir "$1" --forecast "$2"
}
MAPS_PID=""
if [ "$CAM_MAPS" = 1 ] && [ -z "${NWP_SKIP_12KM:-}" ]; then
    P12=forecast.npz
    if [ "${NWP_12KM_MOIST:-0}" = 1 ] && [ -f "$D12/forecast_nhm.npz" ] && \
       [ "$(stat -c %Y "$D12/forecast_nhm.npz")" -ge "$T0" ]; then P12=forecast_nhm.npz; fi
    echo "=== 12 km maps ($P12) started in the background on CPUs $MAP_CPUS  $(date -u +%FT%TZ)" >> "$LOG"
    ( M0=$(date +%s); draw_maps "$D12" "$P12" > "$LOGDIR/maps12_$STAMP.log" 2>&1
      echo "=== 12 km maps finished rc=$?  $(date -u +%FT%TZ)  ($(( ($(date +%s) - M0) / 60 )) min)" >> "$LOG" ) &
    MAPS_PID=$!
fi

# 1b. Moist runs: a moist 12 km NH forecast for the 3 km edges (E4 nhMb set-up:
# bounded start-up is the forecast.py default). Same threads as the 3 km run.
CAM_MOIST_ARGS=""
if [ "$CAM_MOIST" = 1 ]; then
    CAM_MOIST_ARGS="--moist $CAM_ICE_ARG"
    if [ -z "${NWP_CAM_BOUNDARY:-}" ]; then
        BFM="$D12/forecast_nhm.npz"
        # NWP_12KM_MOIST=1 (2026-10-20): daily.sh (step 1) has just made this file as
        # the 12 km product; reuse it when it was written after this cycle started (T0).
        if [ "${NWP_12KM_MOIST:-0}" = 1 ] && [ -f "$BFM" ] && [ "$(stat -c %Y "$BFM")" -ge "$T0" ]; then
            echo "=== 12 km moist NH forecast: reusing the 12 km product $BFM" >> "$LOG"
        else
        rm -f "$BFM"          # never nest in a stale file from an earlier attempt
        step "12 km moist NH forecast (3 km edges)" python -u "$ROOT/src/forecast.py" --run-dir "$D12" \
            --hours "$HOURS" --output-every 1 --core nh --nh-backend c --advection upwind3 \
            --nh-physics step --levels 40 --threads "$CAM_THREADS" --moist $CAM_ICE_ARG --out "$BFM" \
            $INIT12_ARGS --deadline-min 30 || true
        fi
        if [ -f "$BFM" ]; then BFC="$BFM"
        else echo "=== moist 12 km forecast missing: 3 km edges from $BFC (no moisture at the edges)" >> "$LOG"; fi
    fi
fi

# 2. The 3 km analysis
step "3 km analysis" python -u "$ROOT/src/ingest_obs.py" --cycle "$RUN" --from-raw \
    --raw-from "$D12" --out-root "$T3" --spacing-km 3 \
    --background-analysis "$D12/obs_analysis_f00.npz" --gapfill $UPPER_ARGS || exit "$STATUS"

# 3. The 3 km forecast, inside what is left of the budget
USED_MIN=$(( ($(date +%s) - T0) / 60 ))
LEFT_MIN=$(( BUDGET_MIN - USED_MIN - RESERVE_MIN ))
if [ "$LEFT_MIN" -lt 5 ]; then
    echo "=== 3 km forecast SKIPPED: $USED_MIN of $BUDGET_MIN min already used" >> "$LOG"
    exit 1
fi
# Edge relaxation in physical units (P-79, tests E6-E6c, 2026-10-10): u, v and
# pi over 180 km (the 12 km zone) at the 12 km time scale; theta and q over
# 180 km above sigma 0.80 but only 45 km near the ground, where the 12 km
# forecast's boundary layer would cost the 3 km 2 m temperature. Override with
# NWP_CAM_RELAX (empty string = the old 15-cell, CFL-step relaxation).
CAM_RELAX="${NWP_CAM_RELAX---relax-width-km 180 --relax-thermo-width-km 45 --relax-thermo-sigma 0.80,0.90 --relax-dt-ref 17.1}"
step "3 km forecast" python -u "$ROOT/src/forecast.py" --run-dir "$D3" --hours "$HOURS" \
    --output-every 1 --core nh --nh-backend c --threads "$CAM_THREADS" \
    --nh-physics step --levels 40 --boundary-forecast "$BFC" \
    --deadline-min "$LEFT_MIN" $CAM_RELAX $CAM_MOIST_ARGS $INIT3_ARGS ${NWP_CAM_ARGS:-} || true

# Compute ledger (docs/COMPUTE_LEDGER.md): one line per run, kept with the logs.
log_compute() {   # script, threads
    local f="$LOGDIR/compute_hours.csv" wall
    wall=$(( $(date +%s) - T0 ))
    [ -f "$f" ] || echo "start_utc,cycle,script,wall_min,threads,core_h,status" > "$f"
    echo "$(date -u -d @"$T0" +%FT%TZ),$RUN,$1,$(awk "BEGIN{printf \"%.1f\", $wall/60}"),$2,$(awk "BEGIN{printf \"%.2f\", $wall*$2/3600}"),$STATUS" >> "$f"
}
log_compute cam_cycle.sh "$CAM_THREADS"

# 4. The 3 km maps and viewer (NWP_CAM_MAPS=1), after the 12 km maps are done.
if [ "$CAM_MAPS" = 1 ]; then
    [ -n "$MAPS_PID" ] && wait "$MAPS_PID"
    if [ -f "$D3/forecast.npz" ]; then
        M0=$(date +%s)
        echo "=== 3 km maps on CPUs $MAP_CPUS  $(date -u +%FT%TZ)" >> "$LOG"
        draw_maps "$D3" forecast.npz >> "$LOG" 2>&1
        echo "    3 km maps rc=$?  ($(( ($(date +%s) - M0) / 60 )) min)" >> "$LOG"
        echo "$(date -u -d @"$M0" +%FT%TZ),$RUN,maps_3km,$(awk "BEGIN{printf \"%.1f\", ($(date +%s)-$M0)/60}"),2,$(awk "BEGIN{printf \"%.2f\", ($(date +%s)-$M0)*2/3600}"),0" >> "$LOGDIR/compute_hours.csv"
    fi
fi
echo "=== done  $(date -u +%FT%TZ)  (+$(( ($(date +%s) - T0) / 60 )) min)  status=$STATUS" >> "$LOG"
exit "$STATUS"

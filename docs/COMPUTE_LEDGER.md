# Compute ledger

Machine time used by this project, kept from 2026-10-03 at the user's request
(prompt 174). One row per job in [`COMPUTE_LEDGER.csv`](COMPUTE_LEDGER.csv).

## What is counted, and how

- **Server jobs** (met-23030149, 2 × 26-core Xeon Gold 6230R; 104 threads with
  hyper-threading). Each row is a job submitted through the agent's job system.
  - `run_h` runs from the job script's start (its `cmd.sh` written on the server)
    to its last file write.
  - `wait_before_start_min` is the time between submission and start: mostly
    the user approving the job.
- **`threads_max`** is the most threads the job used at once, read from its
  script. Jobs that ran sequential thread sweeps, waited for another job, or
  ended with single-thread verification have `core_h_upper` = run_h ×
  threads_max as an **upper bound**, and say so in `note`.
- **Not counted:**
  - work run on the server before 2026-09-30 by hand or under `nohup`
    (tests up to AA, the September stability ladders), which left no job
    record;
  - the user's own runs;
  - desktop runs shorter than 5 min.
- **Production cycles.** From 2026-10-03, `tools/daily.sh` and
  `tools/cam_cycle.sh` append one line per cycle to
  `$DATA/logs/compute_hours.csv` (start, cycle, script, wall minutes, threads).
  Those runs will be in that file, not here.

## Totals to 2026-10-05 (51 server jobs)

| Work | Jobs | Run hours | Core-hours (upper bound) |
|---|---|---|---|
| Dry 12 km model (tests AB–AN) | 29 | 23.6 | 505 |
| CAM stages S2–S5 (AO, S3a–S5h, benchmarks, profiles) | 22 | 16.9 | 729 |
| **All server jobs** | **51** | **40.5** | **1235** |

Computed from the CSV rows (run_h 40.51, core-h 1234.60). The CAM rows are
the two CAM benchmarks of 2026-10-02 (857ef10a, 1975c024) and every job
from test AO (7efb640f) on; all other server rows count as the dry 12 km
model, including the wait job for test AA.

Approval waits added 3.4 h of wall clock in total (2.1 h of it one test job
submitted while the user was away, 2026-10-04).

The largest single jobs:
- the 3 km campaign S5g, ≤525 core-hours (10.1 h on both sockets);
- the holdout AN, ≤157 core-hours;
- the campaigns AL, AK and AM, 38–74 each.

One 3 km 24 h cycle at 26 threads costs **23 core-hours** alone on the
machine (S5c3: 53.6 min) and 26 core-hours when two run at once, one per
socket (S5g: 59–67 min each). Before the S5e kernels a cycle cost ~40
core-hours and did not finish. The S5g campaign is most of the CAM total.

# NWP1 — 3D Neural Weather Emulator

A Conv3d autoregressive emulator that steps a gridded atmospheric state
forward one hour at a time, trained on HRRR analysis fields.

State vector: **5 variables x 20 pressure levels** — `TMP, RH, UGRD, VGRD, HGT`.
Temperature is treated as the strict control variable and carries a 5x weight
in the loss, so thermodynamic error dominates optimisation.

## Layout

```
config.py                     paths + hyperparameters (env-var driven)
resources.py                  50% resource cap — import before torch
test_env.py                   environment smoke test — run this first
requirements.txt              ingestion machine only
src/
  ingest_hrrr.py              HRRR -> domain-subset tensors
  dynamics/                   the physics core (see its own README)
  verification/               observation operator, QC, scoring, archive
  postproc/                   adaptive bias correction
  nwp_emulator_3d.py          Conv3d encoder/decoder
  autoregressive_dataset.py   pairs state(T) -> state(T+1)
  train_autoregressive.py     training loop, checkpointing, resume
data/tensors_3d/              .pt tensors (gitignored)
data/checkpoints/             model weights (gitignored)
```

## Two-machine setup

Development and data ingestion happen on a local Linux/WSL box; training runs
on the shared Xeon. Code moves by git, data does not.

| | ingestion machine (WSL) | shared server |
|---|---|---|
| installs packages | yes | no — uses preinstalled stack |
| runs Herbie / cfgrib | yes | only if already available |
| runs training | small tests | yes |

Set the data root per machine:

```bash
# WSL
export NWP_DATA_ROOT=/mnt/c/Users/Epier/Desktop/NWP/NWP_Deployment_Package/data
# server
export NWP_DATA_ROOT=/data5/pierce/Data5/NWP/data
```

Add it to `~/.bashrc` on each so you don't have to think about it again.

## Workflow

```bash
python test_env.py                          # what's available here?
python src/ingest_hrrr.py --start 2026-08-01 --hours 6 --dry-run
python src/ingest_hrrr.py --start 2026-08-01 --hours 24
python src/train_autoregressive.py --dry-run  # verify wiring, one batch
python src/train_autoregressive.py --epochs 5
python src/train_autoregressive.py --resume   # continue from latest.pth
```

## Shared-server etiquette

The Xeon (104 cores, 376 GB) is shared with ~30 researchers, so **this
project never uses more than 50% of the machine's cores** unless explicitly
told otherwise — and below that ceiling it scales down further to whatever
other users are leaving free, re-checked between every epoch. The cap is
enforced in `resources.py`, which must be imported before torch — OpenMP and
MKL read their thread limits at import time, and torch will otherwise seize
every core on the box.

Check what the policy resolves to on any machine:

```bash
python resources.py
```

How the adaptive cap behaves at 104 cores (50% ceiling = 50 threads):

| 1-min load | others using | threads we take |
|---|---|---|
| 0-52 | ~0 | 50 (ceiling) |
| 78 | ~28 | 38 |
| 100 | ~50 | 27 |
| 150 | ~100 | 4 (floor) |

Adjustments smaller than 4 threads are ignored so the run doesn't oscillate.
Disable adaptation and pin to the ceiling with `NWP_ADAPTIVE=0`; set the floor
with `NWP_MIN_CORES`.

Raise the ceiling only when you know the box is quiet:

```bash
python src/train_autoregressive.py --resource-fraction 0.75
export NWP_RESOURCE_FRACTION=0.75      # or set it for the session
```

Core count is read via `sched_getaffinity`, so a cpuset or scheduler
allocation is respected rather than the host's full core count.

Memory is reported, not enforced — a hard RSS limit would kill a run mid-epoch
instead of slowing it. Watch the reported budget and reduce batch size or
domain if you approach it.

Long runs should go under `nohup` or `tmux` so they survive a dropped SSH
session.

## Status

Working: model, dataset pairing, training loop, checkpoint/resume, config.

Ingestion writes `data/tensors_3d/<run>/live_hrrr_f<HH>.pt`, each a dict with
`hrrr_features` holding a `[5, 20, Y, X]` tensor plus lat/lon and metadata.

**Analysis vs forecast mode matters.** `--mode analysis` (default) fetches
successive hourly F00 analyses, so T -> T+1 pairs teach real atmospheric
evolution. `--mode forecast` walks one run's F00..FNN, where every step after
F00 is HRRR's own forecast — pairs then teach the model to imitate HRRR rather
than the atmosphere. Train on analysis.

## Test suites

```bash
cd src/dynamics && python test_shallow_water.py    # 8
cd src/dynamics && python test_boundaries.py       # 6
cd src/dynamics && python test_primitive3d.py      # 8
cd src/dynamics && python test_subgrid.py          # 7
cd src/verification && python test_verification.py # 9
cd src/postproc && python test_bias_correction.py  # 7
```

45 tests. Physics tests assert analytic answers or convergence order, not
tolerances chosen to pass.

Still open:
- **normalisation** — geopotential height (~5000) and RH (0-100) differ by
  orders of magnitude, so the 5x temperature loss weight does not yet mean what
  it appears to. Needs per-channel standardisation before results are
  trustworthy.
- data assimilation layer — see `docs/DATA_ASSIMILATION.md`

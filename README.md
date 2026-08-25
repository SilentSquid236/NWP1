# NWP1 — 3D Neural Weather Emulator

A Conv3d autoregressive emulator that steps a gridded atmospheric state
forward one hour at a time, trained on HRRR analysis fields.

State vector: **5 variables x 15 pressure levels** — `TMP, RH, UGRD, VGRD, HGT`.
Temperature is treated as the strict control variable and carries a 5x weight
in the loss, so thermodynamic error dominates optimisation.

## Layout

```
config.py                     paths + hyperparameters (env-var driven)
resources.py                  50% resource cap — import before torch
test_env.py                   environment smoke test — run this first
requirements.txt              ingestion machine only
src/
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
python src/train_autoregressive.py --dry-run  # verify wiring, one batch
python src/train_autoregressive.py --epochs 5
python src/train_autoregressive.py --resume   # continue from latest.pth
```

## Shared-server etiquette

The Xeon is shared with ~30 researchers, so **this project never uses more
than 50% of the machine's cores** unless explicitly told otherwise. The cap is
enforced in `resources.py`, which must be imported before torch — OpenMP and
MKL read their thread limits at import time, and torch will otherwise seize
every core on the box.

Check what the policy resolves to on any machine:

```bash
python resources.py
```

Raise it only when you know the box is quiet:

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

**Not built yet: ingestion.** Nothing currently writes the
`data/tensors_3d/<run>/live_hrrr_f<HH>.pt` files the dataset expects. Each
must contain a dict with key `hrrr_features` holding a `[5, 15, Y, X]` tensor,
channels in the `config.CHANNELS` order.

Also open:
- domain subsetting — full CONUS at 1799x1059 is ~500 MB/hour in float32
- normalisation — raw geopotential and temperature differ by orders of
  magnitude, which will distort the weighted loss until channels are scaled
- data assimilation layer for balloon telemetry into the F00 state

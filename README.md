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
  ingest_hrrr.py              HRRR -> domain-subset .npz (init + boundaries)
  forecast.py                 end-to-end: initialise, integrate, output
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
python src/ingest_hrrr.py --start 2026-08-01 --hours 24
python src/forecast.py --run-dir data/tensors_3d/analysis_20260801_00 --hours 12
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

## Documentation

| document | what it is |
|---|---|
| `docs/AI_COLLABORATION.md` | the research study: AI-assisted model building, error taxonomy |
| `docs/RESEARCH_LOG.md` | dated record of every experiment, including failures |
| `docs/METHODOLOGY.md` | standards of evidence — how we decide something works |
| `docs/CAPABILITIES.md` | honest assessment of what the model can and cannot do |
| `docs/STABILITY.md` | the open stability problem, with all measurements |
| `docs/DATA_ASSIMILATION.md` | observation operator, QC, increments — design |
| `docs/POSTPROCESSING.md` | bias correction and learned post-processing — design |

New log entry: `python tools/newlog.py "Short title"`. Hypothesis before
result.

## What to expect

See **`docs/CAPABILITIES.md`** for an honest assessment. Short version: the dry
core is correct and verified, mid-tropospheric flow is reasonable for 12-24 h,
and surface weather is not — there is no moisture, radiation, boundary layer,
or terrain. Do not use it for forecasts anyone relies on.

Measured runtime for a 12 h forecast (single core, 20 levels): **2.1 min at
12 km, 22 min at 6 km, 4.2 h at 3 km**. The code is single-threaded numpy, so
the 104-core server does not speed it up. Run at 6-12 km — the hydrostatic
core resolves nothing extra at 3 km.

## Before a live run

```bash
python preflight.py --hours 12 --stride 4
```

Checks packages, config, disk, Herbie's cache location, compute budget,
estimated runtime, network reachability, and the fast test suites. Read-only
and quick — it exists to fail in seconds on a missing package rather than
twenty minutes into a download. Exits 0 for GO, 1 for NO-GO with the blocking
items named.

### Resolution: use `--stride`

The ingest subsamples the grid AFTER the domain cut, so the geographic extent
is unchanged and only the spacing coarsens. `build_grid()` derives dx from the
domain and array shape, so the dynamics follows automatically.

| stride | spacing | grid | stored | 12 h forecast |
|---|---|---|---|---|
| 1 | 3 km | 388x438 | 68 MB/h | 4.2 hours |
| 2 | 6 km | 194x219 | 17 MB/h | 22 min |
| **4** | **12 km** | **97x109** | **4.2 MB/h** | **2.1 min** |

**Start at stride 4.** The core is hydrostatic with no convection scheme, so
3 km resolves nothing the equations can represent — it is false precision at
120x the cost. Cost scales as stride^-3.

## Shared-resource policy

Two caps, same philosophy: take a modest share by default, adapt, and never
be the reason a colleague's session is slow.

**CPU** (`resources.py`) — ceiling of 50% of cores, scaling down to what other
users leave free, rechecked between epochs. Note the dynamics is single-
threaded numpy, so this constrains intent more than throughput today.

**Network** (`netpolicy.py`) — token-bucket rate limit (default 8 MB/s, ~6% of
a 1 Gb/s link), sequential requests only, an enforced gap between them,
exponential backoff, and a content-addressed cache so nothing is downloaded
twice. A 13-hour HRRR ingest is ~220 MB, reported up front before it runs.

```bash
export NWP_MAX_MBPS=25          # raise the cap when the link is quiet
export NWP_CACHE_DIR=/data5/pierce/NWP/cache
python src/ingest_hrrr.py --start 2026-08-01 --hours 13 --pause 2
```

Saturating a shared link is more disruptive than saturating a core: everyone's
ssh, file transfer, and data fetch degrades at once, and it is not obvious to
them why.

## Test suites

```bash
cd src/dynamics && python test_shallow_water.py    # 8
cd src/dynamics && python test_boundaries.py       # 6
cd src/dynamics && python test_primitive3d.py      # 8
cd src/dynamics && python test_subgrid.py          # 7
cd src/verification && python test_verification.py # 9
cd src/postproc && python test_bias_correction.py  # 7
cd src/verification && python test_fetchers.py     # 9
cd src && python test_forecast.py                  # 7
python test_netpolicy.py                           # 9
```

70 tests. Physics tests assert analytic answers or convergence order, not
tolerances chosen to pass.

Still open:
- **terrain** — pressure coordinates assume a flat lower boundary, and the
  Appalachians run through this domain. Sigma coordinates are the fix.
- **moisture** — the core is dry: no condensation, no latent heat.
- **data assimilation** — designed, not built. See `docs/DATA_ASSIMILATION.md`.
- **live fetcher verification** — parsers are tested offline against captured
  samples; the network calls themselves are untested.

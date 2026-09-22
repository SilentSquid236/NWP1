# NWP1 — a regional weather model built with an AI collaborator

Read this first. It is written for a fresh session with no memory of how the
project got here, and everything below was paid for at least once.

## What this is

A dry, hydrostatic, terrain-following (sigma) primitive-equation model for the
US Northeast (37.0–47.5 N, 82.0–66.0 W), built from scratch — not a wrapper
around an existing model. It is initialised from HRRR analyses and verified
against **observations only**.

It is also the subject of a research project: *AI to build* — how an AI
collaborator performs when building a numerical weather model, what it gets
wrong, how the errors are caught, and whether lessons transfer. The research
record in `docs/` is half the deliverable, not a side effect.

Human collaborator: Ethan Pierce (epierce2296@gmail.com).

## Constraints that are not negotiable

These bound every design decision and a proposal that breaks one is wrong
however good it is otherwise.

- **The compute server allows no installs.** Shared Xeon, ~30 users, no sudo,
  no conda/venv, admin policy forbids adding packages. Design around what is
  already installed (numpy, torch, herbie, cfgrib are there). Do not suggest
  workarounds. **This includes installing Claude Science itself on that
  server** — run it on the Windows desktop instead, or ask the admin first.
- **git is on the server since 2026-09-22** (it was absent until then). Code
  moves by `git pull` on the server; `tools/pull.sh` (curl + tar) is the
  fallback. The server only pulls — GitHub is the source of truth. `data/` is
  git-ignored, so **`git clean -x`/`-X` and `git stash --all` would take it
  away: never run them there.** git arriving does not relax the no-installs
  rule. Procedure: `skills/nwp-sync`.
- **Never exceed 50% of the server's cores** unless told otherwise, adapting to
  other users' load (`resources.py`). Never saturate shared bandwidth
  (`netpolicy.py`, 8 MB/s default).
- **No model output enters a forecast or scores one** (since 2026-09-22; before
  that, HRRR could seed). Each run starts from observations valid at its own
  cycle time; where soundings are missing, the upper air comes from this
  model's own previous forecast. Nothing observed after the cycle time enters
  the run. Verification is against observations only, after the forecast
  window has closed. `--source hrrr` survives only as a labelled baseline.
- **`data/` is never overwritten** by any sync or tool. It holds the
  verification archive, the one thing here that cannot be recreated.

## How the work is done — the method

**Probe the error, don't guess-check.** The single most consequential
instruction in the project (2026-09-02), given after nine failed patches
against a model that turned out not to be broken. When something fails:
record everything, locate where it fails in space, level, scale and time, and
only then change code. After **two failed adjustments in a row**, stop tuning
and write down the mechanism — on 2026-09-06 that meant writing a feedback
loop and estimating its gain (154·|k| per second), which explained a failure
no setting could fix.

**A fix must predict the outcome it was proposed to explain.** Write the
predictions before the run, including at least one that can only fail. Deepening
the sponge halved the reflection and changed survival by nothing; accepting it
would have been wrong.

**Suspect the test before the model.** Four of the worst defects in this
project were in test setups: a 166 m/s jet clipped at 60, a geostrophic wind
built from one of two PGF terms (845 m/s over terrain), an atmosphere started
at rest and then blamed for accelerating, and a K_MAX ladder that varied
nothing because the constant was bound at import.

**Identical results at different settings are a bug report.** Three settings
giving byte-identical survival counts was a disconnected knob, not a finding.

Full list with origins and callbacks: `docs/LEARNING_LOG.md`.

## Current state (2026-09-22)

The sigma core is complete through boundary layer (Richardson mixing, surface
drag), dry convective adjustment, initialization filtering, a pressure-to-sigma
converter, and a radiative upper boundary.

**Since 2026-09-22 a run starts from observations, not HRRR.** Four cycles a
day (00/06/12/18Z), 12–24 h each, each inside 1.5 h of wall clock:
`tools/daily.sh` → `src/ingest_obs.py` (every reliable source at or before the
cycle time; missing ones skipped) → `src/forecast.py` (edges held to the
initial analysis; deadline; in-hour blow-up stop) → archive. Verification is a
separate job once the window has closed: `bash tools/daily.sh verify` →
`src/verify_pending.py`. The observation fetchers have now met the live
services from the desktop (P-54 found and fixed); the server has not run any
of it yet.

Since 2026-09-22 the work runs in **Claude Science** (desktop app, Windows) with
the model `claude-opus-5-5`; the three skills are imported there and the key
facts in this file are in its project memory. That move is a seam in the study
(`docs/AI_COLLABORATION.md`, "Instrument changes").

Measured capability: 12/12 forecast hours on flat ground, 1000 m and 2500 m
terrain; 8/12 at 4000 m with the eddy-diffusivity ceiling at its new default of
200 m²/s. The agreed terrain target was 2 km and is met.

**Open** (`python tools/problem.py list`):

| | |
|---|---|
| P-07 | **the verification archive on the server has no data** — the only time-sensitive item; a day not archived is gone |
| P-56 | **the first observation-built forecast diverges at 3.75 h**: v at level 17, over 812 m terrain in northern Maine; located, not yet explained |
| P-53 | observation-only initial state with frozen edges — a design with a known cost (error spreads in from the edges) |
| P-55 | the old `daily.sh` gave the forecast the wrong run directory (predicted; fixed; unconfirmed on the server) |
| P-06 | fetchers met the live services from the desktop on 2026-09-22 (P-54 found there); the server has not run them |
| P-49 | every Rayleigh sponge setting turns growing baroclinic weather into decaying weather |
| P-50 | the radiative lid (fixes P-49's development: 0.34 → 3.12) still dies over tall terrain; fourth cause is top-level shear the sponge was masking |
| P-52 | a blowing-up run is only detected at forecast-hour boundaries, so it slows down instead of stopping |
| P-02, P-03 | sponge reflection; 2500 m wind growth past hour 15 |

## Where things are

| | |
|---|---|
| `src/dynamics/` | the model; `primitive_sigma.py` is the core |
| `src/analysis/` | observation analysis: `sources.py` (one adapter per source), `build.py` (first guess, Barnes, heights), `geo.py` (grid, ETOPO terrain) |
| `src/ingest_obs.py` | one cycle's initial state from observations |
| `src/forecast.py` | driver: analysis → sigma → forecast (obs or `--source hrrr` frames) |
| `src/verify.py`, `src/verify_pending.py` | observations → matched pairs → archive, once a window has closed |
| `docs/PROBLEMS.md` | what is wrong, what fixed it, what was ruled out and how |
| `docs/RESEARCH_LOG.md` | dated entries: hypothesis before the result, then the numbers |
| `docs/LEARNING_LOG.md` | lessons, and where each came back |
| `docs/PROMPT_LOG.md` | every human prompt, classified — the study's input record |
| `docs/TOKEN_COST.md` | what the project costs, billed vs API-equivalent |
| `docs/STRUCTURE.md` | generated tree; reading order at the top |
| `skills/` | the recurring procedures, as skills |

**On the server** (checked 2026-09-22): the project root is
**`/data5/pierce/NWP`**, and it already has a `.git`. The data root is **not**
where older docs said: `~/.bashrc` line 44 sets
`NWP_DATA_ROOT=/data5/pierce/NWP/NWP_Deployment_Package/data`, so the
verification archive sits inside a nested copy of the package that has its
own `.git`. **That nested folder looks like leftover junk and is not: removing
it removes `data/`.** `NWP1-main/`, `nwp.tar.gz` and the copies under
`/data5/pierce/Data5/` are older transfers. No crontab is installed.

Every module has a `test_*.py` beside it. The suites are the specification.

## Before claiming anything is done

    python tools/problem.py check      # no FIXED entry without a measurement
    python tools/checklayout.py        # no src/src, no duplicates
    python tools/manifest.py --check   # byte-for-byte against the manifest
    python tools/tree.py               # regenerate STRUCTURE.md
    # and every test suite the change could touch

Then: a research-log entry (hypothesis stated before the result), the problem
register updated in place, the prompt logged, and the session's token counts
recorded with `tools/tokens.py --add`.

## Things a new session gets wrong

- Index 0 in every vertical array is the **model lid**, not the ground. On the
  analysis's PRESSURE levels (`config.PRESSURE_LEVELS`) index 0 is 1000 hPa —
  the opposite way round.
- The analysis grid helpers are `src/analysis/geo.py`, not `grid.py`:
  `src/dynamics/grid.py` (CGrid) is on the same import path.
- The eddy-diffusivity ceiling and Ri_crit are **instance state**
  (`PrimitiveSigma(k_max=...)`). Setting `turbulence.K_MAX` does nothing.
- The sponge relaxes toward a **frozen reference**, not the horizontal mean.
- A test that starts an atmosphere with a temperature gradient and no wind will
  accelerate. That is geostrophic adjustment, not a bug.
- Measurements with no date cannot expire. Date them.

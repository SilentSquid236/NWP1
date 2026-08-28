# Research Log

A dated record of what was tried, what happened, and what it means. Entries
are append-only: superseded conclusions are struck through rather than
deleted, because the path to a result is part of the result.

**Negative results are recorded with the same weight as positive ones.** Most
of the useful findings below are things that did not work.

---

## 2026-08-24 — Environment: Windows is the wrong platform

**Context.** The project began on Windows with a PyTorch/Herbie/MetPy stack.

**What happened.** Three consecutive failures: PyTorch's `fbgemm.dll` missing
a C++ dependency; `conda` invisible to PowerShell; and `ecCodes` failing to
load its C library through three separate fixes (`ecmwflibs`,
`os.add_dll_directory`, `ECCODES_DIR`).

**Finding.** ecCodes is built for Linux HPC. Its Windows port is unreliable
enough that the standard advice is to stop using Windows.

**Resolution.** Moved to a shared Linux server. Later discovered the server
already had the complete stack — eccodes 2.43, cfgrib, herbie, metpy — so the
entire Windows effort was avoidable. **Lesson: inventory the target
environment before building for the development one.**

---

## 2026-08-25 — Neural emulator: a bug, then a dead end

**Design.** A Conv3d network mapping atmospheric state T → T+1, trained on
HRRR analyses. 5 variables × 15 pressure levels.

**Bug found.** `AutoregressiveDataset3D` sorted tensor filenames
lexicographically, so with more than ten forecast hours `f10` sorted between
`f1` and `f2`. Combined with the consecutive-hour check, this silently dropped
most training pairs. Fixed with numeric sorting; regression test added.

**The dead end.** ~~The emulator is the project's core.~~ A model trained on
HRRR output is bounded by HRRR: it learns that model's biases as if they were
physics, and cannot exceed its teacher. Scoring it against HRRR would measure
only how well it copied.

**Decision.** Abandon the emulator. Build a physics core; keep the neural
approach for stage 4, as learned parameterizations *inside* a physical model,
where it adds something the equations cannot express.

**This is the project's central methodological claim** and everything after
follows from it.

---

## 2026-08-25 — Shallow water core: build the small thing first

**Design.** 2D shallow-water equations, Arakawa C-grid, beta-plane,
Wicker–Skamarock RK3.

**Why not go straight to 3D.** Shallow water contains advection, Coriolis, the
pressure gradient and gravity waves — every hard part except vertical
structure — and has *analytic solutions to test against*. A bug found here
takes an afternoon; the same bug in a 3D moist model is nearly invisible,
because a numerical instability is indistinguishable from real convection.

**Validation.** 8 tests: rest stays at rest (exact), mass conserved (exact),
gravity wave speed 0.7% of √(gH), geostrophic balance 0.00% drift over 24 h,
CFL limit real, energy drift bounded.

**Immediate payoff.** 4 of 6 initial tests failed. Three were test configs
violating their own CFL limit. The fourth was subtler: a `tanh` jet is not
periodic in y, so the wrap-around seam created an artificial gradient that
destroyed geostrophic balance within hours. A sinusoid — periodic by
construction — took drift from 131% to 0.00%. **That class of bug would have
been invisible in a 3D model.**

---

## 2026-08-25 — Vector-invariant momentum: right change, wrong reason

**Prediction.** Rewriting momentum in vector-invariant (Sadourny) form would
fix the −4.3% energy drift.

**Result. The prediction was wrong.** Both forms lose identical energy:

| dt | advective | vector-invariant |
|---|---|---|
| dt_max | −4.325% | −4.320% |
| /8 | −0.017% | −0.010% |

Drift shrinks ~7.2× per halving of dt — it is RK3 **time** truncation, not the
spatial scheme. I had assumed a spatial cause without checking.

**What it actually bought.** Potential enstrophy conservation improved **14×**
in a vorticity-rich shear flow (0.0015% vs 0.0218% over 48 h). The first test
showed 0.0000% for both because a smooth blob barely perturbs enstrophy; it
took a rolling-up shear layer to make the diagnostic sensitive.

**Lesson.** When a change does not produce the predicted effect, look for the
effect it *does* produce before keeping or discarding it. Both outcomes are
now encoded as tests.

---

## 2026-08-26 — 3D primitive equations, and an invisible axis bug

**Design.** Dry hydrostatic primitive equations on 20 pressure levels.
Prognostic u, v, θ; geopotential from hydrostatic integration; omega from
continuity.

**Bug found (serious).** `grid.shift` used `axis=1` for x — correct for 2D
`(ny, nx)` fields, but in a 3D `(nz, ny, nx)` field axis 1 is *y*. The core
was differencing north–south when it meant east–west. It produced NaN rather
than plausible output, which was the lucky outcome. Operators now map onto the
last two dimensions, serving 2D and 3D unchanged.

**Negative result: flux-form theta transport is unstable here.** Flux form
conserves the domain integral exactly — but only if discrete continuity holds
exactly, and `diagnose_omega` applies a linear correction to force omega to
zero at the lid and ground. Multiplying that residual by θ (~300 K) is a large
spurious heating. Splitting about the mean profile did not rescue it. Reverted
to advective form (1.3e-05 drift per 12 h) and documented.

**Validation approach that mattered.** The thermal-wind state is balanced only
to discretisation accuracy, so a tolerance on the residual would be arbitrary.
Instead the test refines the grid: imbalance falls 6.58e-04 → 1.65e-04 →
4.14e-05, ratios 3.97 and 3.99 against a theoretical 4.0. **Truncation error
converges; bugs do not.**

---

## 2026-08-26 — Dissipation: continuous vs discrete eigenvalues

**Bug found.** The hyperdiffusion coefficient was derived from the continuous
`k⁴` with `k = π/dx`. The *discrete* Laplacian's response at 2Δx is `4/dx²`,
not `(π/dx)² = 9.87/dx²` — so the damping was **6× weaker than intended**,
e-folding in 18 hours instead of 3.

Caught because the test asserted the damping *time*, not merely that damping
existed. Now derived from the discrete operator eigenvalue: 3.00 h at the grid
scale, 8284 h at 16Δx.

**Second bug.** The stochastic perturbation's spectral filter could remove
every resolvable mode when `length_scale` approached the domain size,
returning a constant field with zero variance — perturbations silently doing
nothing. Now capped, with a hard error rather than a dead field.

**Reported and corrected.** An early baroclinic test showed "eddy energy
×1.6e29", which passed its threshold but was division by near-zero: the seed
was in θ, so eddy *wind* energy started at exactly zero. Re-measured as a
growth rate between day 1 and day 2.

---

## 2026-08-27 — Live data: four interface bugs in a row

Each surfaced only against real HRRR, and each is now pinned by an
offline-runnable test.

1. **`operator.py` shadowed the stdlib**, breaking `collections` and therefore
   `numpy`, with a circular-import traceback that never mentioned the file.
   Renamed `obs_operator.py`.
2. **Herbie's cache defaulted to `~/data`.** A failed write there does not
   raise; the file simply never appears, surfacing later as a
   `FileNotFoundError` from cfgrib. Redirected to the data root.
3. **The GRIB search regex was anchored with `^`.** HRRR index entries begin
   with a colon (`:TMP:850 mb:anl`), so it matched **zero of 708 messages**.
   Herbie downloaded nothing. This is the one bug that could not be caught
   offline — and it now is, via captured index lines.
4. **cfgrib renames variables to CF short names** (`TMP`→`t`, `HGT`→`gh`).
   Resolved through an alias table.

**Meta-observation.** Everything testable offline was tested and worked. Every
failure was at an interface with an external system whose conventions I had
assumed. That is a reusable prior for this kind of work.

---

## 2026-08-28 — Initialisation: analysis data is not balanced for our grid

**Symptom.** First forecast from real HRRR diverged in 1 hour with
`max|omega| = 131 Pa/s`, where the real atmosphere is order 1 Pa/s.

**Diagnosis.** Working back from the timestep showed initial omega was already
~50 Pa/s before a single step. HRRR winds are balanced for *HRRR's*
discretisation. Coarsened and differenced with our operators they carry ~100×
too much grid-scale divergence, and the column integral converts that to tens
of Pa/s of vertical motion.

**Hypotheses tested and rejected:**
- *Aliasing from strided coarsening* — block-averaging did not reduce it.
- *A-grid vs C-grid wind placement* — interpolating to faces changed nothing.
- *An operator bug* — the divergence operator returns **1.08e-18** on a
  discretely consistent rotational flow. Machine precision. Not the cause.

**Solution (kept).** Helmholtz split: solve `∇²χ = div` in Fourier space using
the eigenvalues of *our discrete* Laplacian, so cancellation is exact rather
than approximate. Subtracting `∇χ` leaves the rotational flow.

| | before | after |
|---|---|---|
| max\|div\| | 7.5e-04 1/s | 9.3e-05 1/s |
| implied omega | 60 Pa/s | 1.08 Pa/s |
| correlation with rotational flow | — | 0.997 |

Boundary frames are balanced too — otherwise relaxation re-injects at the
edges what was removed from the interior, every step.

---

## 2026-08-28 — Stability: a structural limit, not a tuning problem

**Symptom.** With a healthy initial state, omega grows ~4× per hour and the
run dies at 2–3 hours.

**Everything tried, measured.** Hyperdiffusion damping times 3 h → 0.5 h,
crossed with divergence damping 0.0 → 0.2. **No configuration reaches 6
hours**, and the strongest setting is worse than a moderate one.

**Negative result with a cost.** Divergence damping extends survival 1 h → 3 h
but suppresses the physics:

| div_damp | baroclinic growth/day |
|---|---|
| 0.00 | 1.21× |
| 0.01 | 0.33× |
| 0.10 | 0.49× |

Any level that helps stability also damps baroclinic development. Default off.

**Two of my own bugs surfaced here.** Divergence damping written as a tendency
with a coefficient in m²/s violated its own diffusion stability limit
(`nu·dt/dx² ≤ 0.25`) — rewritten as a dimensionless post-step filter, stable
for any dt. And the sponge layer first relaxed wind toward the horizontal
mean, which the thermal-wind test caught: that flattens a jet, which is
legitimate structure, not wave noise.

**Diagnosis.** The numerics are verified and all 29 idealised tests pass; the
instability appears only with realistic sheared, noisy states. Most likely
cause is the **vertical coordinate**: in pure pressure coordinates with a
rigid flat lower boundary, omega must vanish at both ends, and diagnosing it
from divergence while enforcing both conditions over-constrains the column.
The correction feeds vertical advection → divergence → omega, a tight loop
with no physical damping.

**Status.** Correct dry solver for smooth balanced states; not usable for
forecasts from real analyses. Fix is sigma coordinates — already needed for
terrain, now also for stability. Full detail in `docs/STABILITY.md`.

---

## Recording for the AI-collaboration study

Each entry should also note, where applicable:

- **defects introduced**, tagged with a category from `docs/AI_COLLABORATION.md`
  (A external-interface, B discrete-vs-continuous, C stability/dimensional,
  D array/language semantics, E test design, F wrong causal hypothesis)
- **how each was detected** — offline test, real data, targeted measurement,
  self-review, human observation
- **whether a stated hypothesis survived measurement**
- **human interventions**, separating direction-setting from correction

---

## Template for new entries

```markdown
## YYYY-MM-DD — Short title

**Context.** What prompted this.

**Hypothesis.** What was expected, stated before the result.

**Method.** What was run. Enough to reproduce.

**Result.** Numbers. Tables where there is more than one.

**Interpretation.** What it means, and what it does not mean.

**Status.** Kept / reverted / open. If reverted, why.
```

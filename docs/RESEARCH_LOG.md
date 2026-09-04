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

## 2026-08-28 — Sigma coordinate: structural fix for the stability failure

**Context.** The pressure-coordinate core is stable on smooth balanced states
and diverges within 2-3 hours on real HRRR analyses, at every damping setting
tried. Diagnosis pointed at the vertical coordinate, not the numerics.

**Hypothesis.** In pure pressure coordinates with a rigid flat lower boundary,
omega must vanish at both ends of the column. Diagnosing it from divergence
and then enforcing both conditions over-constrains the system; the correction
that pins omega = 0 at the ground redistributes error through the column every
step. In sigma coordinates the ground is sigma = 1 by definition and surface
pressure is PROGNOSTIC, so sigma_dot = 0 at both ends should fall out of the
formulation with no correction at all.

**Method.** Implemented `sigma.py`: stretched sigma levels, hydrostatic
integration, prognostic-pi continuity with diagnosed sigma_dot, flux-form
vertical advection, and the sigma pressure-gradient force. Seven validation
tests, two of which the pressure-coordinate version could not pass.

**Result.** 7/7.

| test | result |
|---|---|
| hydrostatic exact, isothermal | 5.9e-12 m error |
| **sigma_dot = 0 at lid and ground, NO correction** | **0.00e+00 at both** |
| column mass tendency sums to zero | 3.9e-14 vs 9.7e+07 Pa |
| vertical advection of a constant | 8.9e-16 K/s (flux form) |
| PGF cancellation over a 1500 m ridge | 0.505% residual |
| ground is sigma=1 at 0-3000 m terrain | exact |

**Interpretation.** The hypothesis holds for the boundary condition:
sigma_dot vanishes exactly at both ends with no correction, removing the
feedback loop identified as the likely instability mechanism. Whether that
yields a stable 12-hour forecast on real data is NOT yet established -- the
3D prognostic core has not been ported to this coordinate.

**Defects introduced (category C).** The pressure-gradient force was derived
wrong TWICE: first as `-grad(Phi) - R T grad(ln p_s)` (valid only at sigma=1
with p_top=0), then with the second term's sign flipped. Both caught by the
same test -- an isothermal atmosphere in exact balance over a ridge, where the
two large terms must cancel. Measuring which combination cancelled (A-B,
residual 2.9e-04, against A+B at 1.2e-01) settled what derivation had not.

**Detection.** Targeted measurement against an analytic balance. Neither
error would have been visible in a forecast; both would have produced a
plausible but wrong flow over terrain.

**Status.** Coordinate layer complete and validated. Next: port the 3D
prognostic core onto it, then re-run the real-data case that fails today.

---

## 2026-08-28 — Sigma 3D core: partial success, stability still open

**Context.** Port the 3D prognostic core onto the validated sigma coordinate
layer and re-run the case that dies at hour 3 in pressure coordinates.

**Hypothesis.** Prognostic surface pressure removes the over-constrained
lower boundary, so the divergence/vertical-velocity feedback disappears and
realistic states integrate stably.

**Result. 4/6 — the hypothesis is NOT confirmed.**

Verified working:

| test | result |
|---|---|
| rest over flat ground, 12 h | 0.00e+00 (exact) |
| surface pressure evolves under divergent flow | 255 Pa over 3 h |
| total mass conserved, 12 h | 0.00e+00 (exact) |
| thermal-wind jet persists 24 h | 0.12% drift, 0.07% spurious v |

Still failing:

| test | result |
|---|---|
| rest over a 1200 m mountain, 6 h | 9.1 m/s spurious wind |
| 12 h from a realistic noisy sheared state | diverges at 1-4 h |

**Two real findings.**

*The external mode.* With prognostic surface pressure the fastest signal is
the Lamb wave at sqrt(R·T) ~ 290 m/s, not the ~100 m/s internal wave. A rigid
lid suppresses that mode, so the pressure-coordinate CFL carried over made dt
3x too large: surface pressure went NEGATIVE within 20 steps. Fixed by
computing the wave speed from the temperature field. This is why operational
models sub-step or semi-implicitly treat the external mode rather than
resolving it explicitly — a cost that arrives with the free surface.

*Defect introduced (category C).* I added hyperdiffusion to the prognostic
`pi` tendency to damp grid-scale surface-pressure noise. Hyperdiffusion is
only conservative on a periodic domain; on a bounded domain it is a MASS
SOURCE. Measured: p_s inflating 1088 → 1243 hPa over three hours. Removed.
The existing mass-conservation test did not catch it because that test uses a
periodic domain with `hyper=0` — a gap in coverage, not a gap in the code.

**Test-design errors of my own (category E).** The realistic-state test built
its balanced wind and then CLIPPED it at ±60 m/s, which destroys geostrophic
balance exactly where it clips and imposes a large artificial imbalance. And
because the sigma column extends to 50 hPa rather than the pressure version's
200 hPa, the same temperature gradient produces roughly double the jet — my
"gentle" configurations were generating 110-210 m/s jets, far outside
anything realistic.

**Interpretation.** Sigma fixed what it was predicted to fix — the boundary
condition is now exact with no correction, mass conserves exactly, terrain is
representable, and balanced flow is better preserved than before (0.12% vs
3.93%). It did NOT deliver a stable integration from a noisy realistic state.

Remaining suspects, in order: the sigma pressure-gradient cancellation over
terrain (9 m/s spurious over 6 h is too large and is a known weakness with a
known fix — computing the PGF as departures from a reference state); the
explicitly resolved external mode; and my synthetic initial states being
unrepresentative of real analyses.

**Status.** Coordinate layer validated (7/7). 3D core 4/6, and the two
failures are the ones that matter. **Stopping the patch-and-retest loop here**
-- five consecutive fixes each moved the failure without removing it, which
is the signature of an unidentified root cause rather than a list of bugs.
Next step is diagnosis, not another damping term: instrument where the energy
enters, rather than guessing which sink to add.

---

## 2026-08-28 — Instability diagnosis: instrument first, then narrow

**Context.** Five consecutive fixes had each moved the sigma core's failure
without removing it. Stopped patching and built a diagnostic instead.

**Method.** `src/dynamics/diagnose_growth.py` answers four questions by
measurement rather than by hypothesis:

1. **Which term?** dKE/dt = integral of u . (du/dt) evaluated per tendency term.
2. **Which levels?** the same, resolved vertically.
3. **Which scale?** amplitude spectrum of u over time, binned by wavenumber.
4. **Rotational or divergent?** Helmholtz split of the growing part.

**Result.** The failure localised immediately.

| question | answer |
|---|---|
| which term | pressure gradient, +5.9e+04 (largest source) |
| which levels | the TOP of the model (sigma=0.027, ~76 hPa) |
| which scale | meso/synoptic grow 5.2-5.8x; **grid scale only 1.4x** |
| which component | divergent energy grows 3x; rotational flat |

Grid-scale growth would mean a numerical mode. It is not grid scale, so it is
not that.

**The decisive control.** Varying terrain and noise independently:

| noise | terrain | jet | PGF work | survived |
|---|---|---|---|---|
| 0.0 | 0 m | 164 m/s | **0.000e+00** | **6/6** |
| 0.0 | 400 m | 184 m/s | -2.7e-11 | 2/6 |
| 1.2 | 400 m | 173 m/s | -7.7e+03 | 1/6 |

**On flat ground the model is stable with a 164 m/s jet and exactly zero PGF
work** -- discrete geostrophic balance is perfect. Terrain breaks it with no
noise at all. Noise is a modest aggravator, not the cause.

**Two candidate fixes tested and REJECTED.**

*Reference-state PGF.* Rewrote the force as
`-grad(Phi + R T0 ln p) - R (T - T0) grad(ln p)`, so the large terms cancel
analytically rather than numerically. Standard remedy for sigma-coordinate
pressure-gradient error. **No effect on survival time.**

*Full-PGF geostrophic initialisation.* My initialiser balanced only against
-dPhi/dy, which is correct on flat ground (where grad(pi) = 0) but omits half
the force over terrain. Fixed to balance against the complete sigma PGF:
initial PGF work dropped to **-3.2e-14**, machine zero, confirming a genuinely
balanced state. Survival improved 2/6 -> 3/6 and **the run still dies.**

**Interpretation.** A perfectly balanced state over 400 m of terrain diverges
in ~3 hours. The instability is in the model's treatment of terrain, not in
the initial state, not in the pressure-gradient formulation, and not in
grid-scale noise. Energy enters through the pressure-gradient term at upper
levels in the divergent component -- consistent with error accumulating
upward through the hydrostatic integral, whose absolute magnitude is largest
at the top.

Untested candidates, now narrow: the hydrostatic integration over a
horizontally varying pi (layer thicknesses differ column to column, and the
integral starts from a terrain-following surface); the stretched sigma grid
interacting with terrain slope; and vertical resolution at the model top.

**Status.** Open. But the question has gone from "why does it blow up" to
"why does balanced flow over terrain leak energy into divergent modes at the
model top" -- which is answerable.

**Method note for the AI-collaboration study.** Five patch-and-retest cycles
produced no progress; one instrument produced a decisive localisation in a
single run. The instrument also **falsified two plausible fixes** that would
otherwise have been adopted on the strength of sounding right. The human
called for this change of approach.

---

## 2026-08-28 — Terrain sweep: the instability scales with slope

**Context.** The diagnostic localised energy entry to the pressure-gradient
term at upper levels, in the divergent component. Next question: does the
failure scale with terrain, and if so with height or with slope?

**Method.** Sweep terrain from flat to extreme with everything else fixed,
initialising against the FULL sigma pressure-gradient force so the state is
genuinely balanced at every height. Then vary the vertical discretisation
independently at fixed terrain.

**Result 1 — monotonic in slope.**

| terrain | max slope | jet | PGF work at t=0 | survived |
|---|---|---|---|---|
| 0 m | 0 | 176 m/s | **0.000e+00** | **4/4** |
| 250 m | 8.5e-04 | 176 m/s | -6.2e-14 | 3/4 |
| 1000 m | 3.4e-03 | 175 m/s | +3.3e-13 | 3/4 |
| 2500 m | 8.5e-03 | 174 m/s | -2.5e-13 | 2/4 |
| 5000 m | 1.7e-02 | 188 m/s | +4.5e-13 | 1/4 |

Initial PGF work is machine zero at EVERY terrain height, so the state is
perfectly balanced in all cases and the error is generated during
integration. Flat ground is stable with a 176 m/s jet.

**Result 2 — the vertical grid barely matters, but the lid does.**

At 2500 m terrain:

| stretch | p_top | nz | survived |
|---|---|---|---|
| 1.4 | 50 hPa | 20 | 2/4 |
| 1.0 (uniform) | 50 hPa | 20 | 2/4 |
| 1.4 | **200 hPa** | 20 | **3/4** |
| 1.0 | 200 hPa | 20 | 3/4 |
| 1.4 | 50 hPa | 30 | 2/4 |

Level stretching and level count: no effect. Model top: raising it from
200 to 50 hPa costs an hour. Consistent with energy entering aloft, where
R*T/p is largest and the hydrostatic integral has accumulated furthest.
**Default p_top changed to 200 hPa** -- a dry model with no stratospheric
physics gains nothing from a 50 hPa lid.

**Result 3 — Coriolis is not implicated.** The budget showed non-zero
Coriolis work, which should be identically zero. Checked directly: exactly
neutral on an f-plane (1.7e-16, machine precision) and the averaging
operators are exactly adjoint. On a beta-plane the error is 1.7e-08 relative
-- an e-folding time near 700 days. Negligible. A plausible-looking suspect,
eliminated in two minutes by measurement.

**Candidates tested and rejected across this and the previous entry:**
reference-state PGF (no effect), full-PGF balanced initialisation (+1 hour),
level stretching (none), level count (none), grid-scale damping (none),
noise (aggravator only, not cause), Coriolis energy error (negligible).

**Interpretation.** The instability is specific: **balanced flow over sloping
terrain leaks energy into divergent modes at the model top, at a rate that
scales with terrain slope.** Everything else has been eliminated by
measurement rather than by argument.

The remaining untested candidate is the hydrostatic integration itself. Phi
is built by integrating upward from a terrain-following surface, so each
column accumulates a different path, and horizontal differences of Phi at
upper levels are differences between separately accumulated integrals. That
is exactly where a slope-dependent, top-heavy, divergent error would come
from. Testing it means comparing against an integration formulated to keep
horizontal consistency -- not a damping term.

**Status.** Open, but narrow and specific. Six candidates eliminated, one
identified and untested.

---

## 2026-08-28 — Hydrostatic consistency: a 10^11 improvement that did not fix it

**Context.** The terrain sweep showed survival falling monotonically with
slope, with a perfectly balanced initial state. Remaining suspect: the
hydrostatic/pressure-gradient discretisation.

**Criterion.** If temperature is a function of pressure alone, the true
pressure-gradient force is EXACTLY zero over any terrain. A discretisation
that does not reproduce that is hydrostatically inconsistent, and its residual
is a spurious force proportional to slope. This is a sharp, machine-precision
test -- much stronger than "the terms nearly cancel".

**Where the error actually lived.** Resolving the residual vertically
overturned my earlier reading:

| level | pressure | residual |
|---|---|---|
| 0 (top) | 205 hPa | 3.0e-05 |
| 10 | 478 hPa | 1.2e-03 |
| 19 (surface) | 861 hPa | **2.1e-03** |

Largest at the BOTTOM, in the layers sitting on the terrain where sigma
surfaces are most steeply tilted -- the boundary layer. The energy budget had
pointed at the top because it weights by wind speed and the jet is aloft; the
FORCE error is near-surface. Two different questions with two different
answers, and I had conflated them.

**Four formulations measured** (residual at 3000 m terrain, must be zero):

| formulation | residual |
|---|---|
| sigma·grad(pi)/p at cell centres (what we had) | 2.1e-03 |
| same, coefficient averaged to velocity points | 2.3e-05 |
| flux form, grad(R T ln p) | 1.0e-14 |
| **grad(ln p) differenced directly, T on velocity points** | **8.2e-15** |

**The bug.** I had expanded `grad(ln p)` analytically as `sigma·grad(pi)/p`.
Algebraically identical; discretely not -- the expansion is evaluated at cell
centres while `grad(pi)` lives on the faces. Differencing `ln p` directly and
averaging T to the velocity point is exact. Improvement: **2.1e-03 to 8e-15**,
eleven orders of magnitude.

Also tried and rejected: the Simmons-Burridge half-level + alpha geopotential
construction. Both it and the naive integration give IDENTICAL residuals
(5.83e-05 vs 5.98e-05 at 500 m), so the geopotential was never the problem.
My first SB implementation was also wrong -- I differenced the half-level Phi
where the force needs the full-level gradient -- and made things 300x worse
before I caught it.

**Result: the fix did NOT resolve the instability.** With the PGF exact to
machine precision, terrain runs still fail:

| terrain | survived (8 h target) |
|---|---|
| 0 m | 7/8 |
| 1000 m | 3/8 |
| 2500 m | 3/8 |
| 5000 m | 1/8 |

Still monotonic in slope. So hydrostatic inconsistency was real, was worth
eleven orders of magnitude, and was **not the cause**.

**Interpretation.** The remaining error is not in the pressure-gradient force
for a horizontally uniform temperature. It must involve the terms that vanish
in that test: horizontal temperature gradients on tilted sigma surfaces, or
vertical advection through them. The next measurement should hold terrain
fixed and vary the horizontal temperature gradient, which the consistency test
by construction cannot see.

**Method note.** This is the second time a plausible, standard, textbook fix
produced a large measured improvement in the thing it targets and no
improvement in the thing that matters. Worth keeping both facts: the PGF is
now correct and should stay correct; and correctness there was not sufficient.

**Status.** PGF fix kept (test now asserts machine-zero consistency rather
than a 2% tolerance). Instability open. Seven candidates eliminated.

---

## 2026-08-28 — Visualisation: what the instability actually looks like

**Context.** Survival counts say a run failed; they do not say what failed.
Built `src/dynamics/visualize_instability.py` -- vertical cross-sections of
sigma_dot through the terrain at successive forecast hours, flat versus
mountain, on one shared colour scale.

**What the picture shows.**

*Flat ground:* vertical velocity is **horizontally uniform** -- clean
horizontal bands spanning the whole domain, alternating sign between hours.
The domain breathing as one column. Amplitude steady near 4e-05, harmless.

*2500 m mountain:* by +1 h **cellular structure appears on the mountain
flanks** -- alternating up/down columns of roughly 200 km wavelength -- and
extends upward through the depth of the model. By +3 h the field is
disordered, with narrow near-grid-scale vertical stripes, and the run dies.

Growth curve: flat stays flat; the mountain climbs steadily and jumps ~4x
between hours 2 and 3.

**Hypothesis it suggested.** Flow over terrain excites vertically propagating
gravity waves; the rigid lid reflects them; upgoing and reflected waves
interfere and grow. This fits every measured fact -- terrain-only,
slope-scaling, divergent, concentrated aloft. And the sigma core had **no
absorbing layer at all**: the sponge written for the pressure-coordinate
version was never ported.

**Test.** Added a Rayleigh sponge over the top 5 levels relaxing wind toward
the initial reference state (not the horizontal mean, which flattens jets --
the mistake the thermal-wind test caught previously).

**Result: no effect whatsoever.**

| terrain | no sponge | 5-level sponge |
|---|---|---|
| 0 m | 7/12 | 7/12 |
| 1000 m | 3/12 | 3/12 |
| 2500 m | 3/12 | 3/12 |
| 5000 m | 1/12 | 1/12 |

Identical to the hour. Lid reflection is not the mechanism. **Eighth
candidate eliminated.** The sponge is kept (it is correct and costs nothing)
but it is not a fix.

**What the visualisation earned regardless.** It converted "diverges at hour
3" into a specific picture: a ~200 km cellular disturbance forming on the
terrain flanks, not at the peak, and filling the column. That is a
length scale and a location, both of which are constraints on any future
hypothesis. The flank-not-peak detail matters -- a peak-centred error would
suggest the terrain representation itself; flanks suggest the SLOPE, which
matches the slope-scaling result exactly.

**Status.** Open. Eight candidates eliminated by measurement. The remaining
question is sharper than before: what generates a 200 km cellular divergent
disturbance on sloping coordinate surfaces when the pressure-gradient force
is exact to machine precision for horizontally uniform temperature?

---

## 2026-08-28 — Orographic response: slope, not height, and not wave breaking

**Context.** The cross-sections showed a ~200 km cellular disturbance forming
on the mountain FLANKS. Question: is the model producing a physical mountain
wave that then breaks, or a numerical artifact?

**Test 1 — wave breaking, rejected.** Orographic theory predicts breaking
above a nondimensional mountain height Nh/U ~ 0.85. Varying U at fixed
terrain:

| h | U | Nh/U | survived |
|---|---|---|---|
| 2500 m | 20 | **1.89** | 6/8 |
| 2500 m | 40 | 0.95 | 4/8 |
| 2500 m | 80 | **0.47** | 3/8 |

Breaking would make Nh/U = 1.89 the WORST case. It is the best. Survival
tracks U, not Nh/U. **Not wave breaking.**

**Test 2 — the scaling is forced ascent.** Varying height and width
independently, so slope and height decouple:

| w = U*slope | terrain | survived |
|---|---|---|
| 0.170 m/s | 1500 m, 300 km | 7/8 |
| 0.171 m/s | **3000 m**, 600 km | 6/8 |
| 0.334 m/s | 1500 m, 150 km | 6/8 |
| 0.341 m/s | 3000 m, 300 km | 4/8 |
| 0.343 m/s | **6000 m**, 600 km | 5/8 |
| 0.667 m/s | 3000 m, 150 km | 4/8 |

**A 6000 m mountain with a gentle slope outlives a 3000 m mountain with a
steep one.** Height is not the variable; slope is. Survival tracks the
terrain-forced vertical velocity w = U * slope.

**Test 3 — a wrong premise, corrected by measurement.** I tested "advective
consistency": with theta = theta(p) and uniform flow over terrain, is the
advective tendency zero? Measured 13 K/hour at 3000 m and U = 40, with the
vertical term at 1e-18. I first read that as a cancellation failure.

It is not. With flow along sigma surfaces, sigma_dot = 0 is CORRECT -- the air
follows the terrain up and over -- and the horizontal term is the physical
adiabatic cooling of rising air. Check: w = U * slope = 0.8 m/s, dry adiabatic
lapse 9.8 K/km, gives ~28 K/hour against 13 K/hour measured. Same order. The
tendency is real, and my test premise was wrong.

**Test 4 — smoothing does not help here.** Eight 1-2-1 passes move the slope
only 1.8e-02 -> 1.4e-02 and survival not at all (4/10 throughout). A 1-2-1
filter removes grid-scale roughness; this mountain is smooth at 150 km and its
slope is RESOLVED. Operational orography filtering works because raw terrain
carries grid-scale structure -- it cannot rescue a genuinely steep resolved
mountain.

**Interpretation.** The orography is forcing the model hard and physically:
0.7 m/s of ascent, driving ~13 K/hour of adiabatic cooling. That forcing is
real and correctly computed. What the model lacks is any means of handling the
response -- no gravity-wave drag, no turbulence, no boundary layer. The
disturbance on the flanks is where the forcing is strongest, which is exactly
where slope peaks.

**IMPORTANT CAVEAT.** The original failing 12-hour real-data forecast used
**flat terrain** -- `forecast.py` does not yet ingest orography. So everything
in this entry is a SEPARATE problem from the one that started the
investigation. The flat-ground sigma case still fails at 7/12 hours with no
terrain at all. Two distinct issues; do not conflate them.

**Tools added.** `terrain_slope`, `forced_ascent`, `smooth_terrain` in
`sigma.py` -- so orography can be characterised before a run rather than
diagnosed after one.

**Status.** Orographic behaviour now understood and quantified. The flat-ground
failure remains open and is the one that matters for the real forecast.

---

## 2026-08-28 — Flat-ground failure was my test case; missing physics is the rest

**Context.** With orography understood, the flat-ground failure (7/12 hours)
remained -- and that is the one blocking the real forecast, since
`forecast.py` does not yet ingest terrain.

**Finding 1: the flat-ground instability was an artefact of my test jet.**

| jet | Rossby number | survived |
|---|---|---|
| 97 m/s | 6.94 | 7/12 |
| **49 m/s** | 3.47 | **12/12** |
| 24 m/s | 1.73 | 12/12 |
| 3 m/s | 0.23 | 12/12 |

My synthetic jet was 97 m/s across a 250 km feature -- Rossby number ~7, so
the flow is nowhere near geostrophic and a "geostrophically balanced" initial
state is badly unbalanced in fact. The adjustment radiates as gravity waves,
surface pressure rings at +-16 hPa from the first hour, and the run eventually
goes nonlinear. **At realistic jet strengths the model integrates 12 hours
cleanly on flat ground.** Fourth test-design error of the session (category E).

Two diagnostic dead ends along the way, both my own premises: uniform flow
with no pressure gradient "growing" at 60 m/s turned out to be a textbook
INERTIAL OSCILLATION (predicted 60.6 m/s at t=6 h, measured 60.84); and
balancing surface pressure against the surface wind changed nothing, because
initial dp_s/dt was already 0.00 hPa/h.

**Finding 2: what remains is missing physics, not numerics.** With a realistic
jet, survival tracks vertical shear between adjacent levels:

| max shear per level | survived |
|---|---|
| 3.2 m/s | 12/12 |
| 4.2 m/s | 3/12 |
| 8.4 m/s | 2/12 |

Shear of ~8 m/s across a ~400 m layer with N ~ 0.015 gives Ri ~ 0.5,
approaching the Ri = 0.25 threshold below which shear instability is
physically expected. **The instability is real.** The model simply has no
turbulence to mix it away -- exactly the same conclusion as the orographic
case, where real forced ascent had no drag to absorb it.

**Action: Richardson-number vertical mixing** (`turbulence.py`). Standard
Louis-type formulation: `K = l^2 |S| f(Ri)`, with `f` falling to zero at
Ri_c = 0.25 and full strength where the column is statically unstable. Fluxes
vanish at lid and ground, so it redistributes momentum and heat within a
column without creating either.

| terrain | noise | no mixing | with mixing |
|---|---|---|---|
| 0 m | none | 12/12 | 12/12 |
| 0 m | 1.2 | 2/12 | **6/12** |
| 1000 m | none | 8/12 | **12/12** |
| 1000 m | 1.2 | 2/12 | **6/12** |
| 2500 m | 1.2 | 3/12 | **6/12** |

**The first change all session that improved anything.** Nine previous
candidates -- reference-state PGF, divergence damping, sponge layers,
hyperdiffusion strength, level stretching, level count, Coriolis, terrain
smoothing, balanced surface pressure -- changed nothing measurable. This one
helps everywhere, and it is physics the model was simply missing rather than a
numerical patch.

Not a complete fix: noisy cases reach 6/12, not 12/12.

**Interpretation.** Two of the three failures now have the same explanation:
the model produces physically correct responses (orographic ascent, shear
instability) and lacks the parameterized physics that would dissipate them.
That is a different and much more tractable problem than a numerical bug --
and it is what `docs/CAPABILITIES.md` predicted at the outset, listing the
boundary layer as a first-order omission.

**Status.** Flat ground + realistic jet: 12/12, clean. Terrain and noisy
states: improved but incomplete. Next candidate is surface drag, the other
half of a boundary-layer scheme, which would damp the near-surface shear that
the mixing scheme currently has to handle alone.

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

## 2026-09-01 — Sigma coordinate: structural fix for the stability failure

**Context.** <what prompted this>

**Hypothesis.** <what you expect, written BEFORE the result>

**Method.** <what was run; enough to reproduce>

**Result.** <numbers; a table if more than one>

**Interpretation.** <what it means, and what it does not>

**Status.** <kept / reverted / open — and why>

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

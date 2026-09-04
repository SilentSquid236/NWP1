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

## 2026-09-02 — Surface drag measured, and the instability traced to the test

**Context.** Richardson mixing had lifted 12-hour survival on the hard cases
from 2-3/12 to 6/12 and stopped there. Surface drag was the stated next
candidate: mixing redistributes momentum inside a column, only drag removes
it. Rather than adding drag and re-running the single decisive test, both
schemes were turned independently over the same terrain x noise grid, so a
change could be attributed to a scheme rather than to a coincidence.

**Hypothesis (stated before the runs).** Drag would raise survival on the
noisy and terrain cases by damping near-surface shear that mixing alone has to
handle.

**Method.** `src/dynamics/sweep_boundary_layer.py`: terrain 0 / 1000 / 2500 m
x noise 0.0 / 1.2 m/s x mixing off/on x drag off/on, 24 runs, 90x88x20,
12-hour ceiling, failure = non-finite or max|u| > 150 m/s.

**Result (first sweep).** Drag changed *nothing*. Twelve matched pairs, twelve
identical survival counts, max|u| agreeing to about 0.5 m/s. Mean hours
survived: neither 1.33, mixing 2.67, drag 1.33, both 2.67. The absolute
numbers were also worse than the 6/12 previously recorded, so the sweep was
not reproducing the earlier baseline either.

**Probe instead of patch.** `src/dynamics/probe_failure.py` recorded surface
pressure, wind, and every momentum term each step, then located the growing
mode in level, latitude and wavenumber.

| measurement | result |
|---|---|
| min surface pressure at failure | 101.0 kPa — never approached zero |
| max\|u\| trajectory | 63.8 -> 161 m/s over 1.2 h, then non-finite |
| growing scale, clean case | **2dx**, peaking two rows from the boundary |
| 2dx damping, interior | e-folding 10 800 s (3 h, as designed) |
| 2dx damping, replicate edge | e-folding 18 400 s (1.7x weaker) |
| observed 2dx growth | doubling in ~20 min |

Damping was losing the race by roughly a factor of nine. The question then
became what *generates* a 2dx mode in a clean, supposedly balanced state.

**Three defects, all in the initial state, none in the dynamics.**

1. The 6 K meridional temperature contrast implies a **166 m/s** jet by
   thermal wind (Ro = 3.2). The test then clipped the wind at +/-60 m/s,
   destroying geostrophic balance over **33.6%** of the domain. That clip is
   the 2dx source. Category E (test design) — the third time an unrealistic
   test jet has been mistaken for a model instability.
2. The balanced wind was taken as `-d(phi)/dy / f`. On sigma surfaces the
   horizontal force has two terms that largely cancel over sloping ground;
   keeping only the first implies an **845 m/s** "balanced" wind over 2500 m
   terrain. Every terrain row of the earlier mixing baseline was measured
   against that. Category B (discrete-vs-continuous / formulation).
3. 1.2 m/s of **white** noise puts 89% of its variance at wavelengths the grid
   cannot carry. Real analyses are filtered before they are integrated; this
   one was not. Category E.

**Result (second sweep, corrected initial states).** 1.5 K contrast, no clip
(41 m/s jet), wind from the full sigma PGF:

| terrain | noise | neither | mixing | drag | both |
|---|---|---|---|---|---|
| 0 m | 0.0 | 12/12 | 12/12 | 12/12 | 12/12 |
| 0 m | 1.2 | 1/12 | 1/12 | 1/12 | 1/12 |
| 1000 m | 0.0 | 7/12 | 12/12 | **11/12** | 12/12 |
| 1000 m | 1.2 | 1/12 | 1/12 | 1/12 | 1/12 |
| 2500 m | 0.0 | 5/12 | 6/12 | **7/12** | 8/12 |
| 2500 m | 1.2 | 1/12 | 1/12 | 1/12 | 1/12 |

Means: neither 4.50, mixing 5.50, drag 5.50, both 6.33. **Drag is worth as
much as mixing over terrain and the two are partly additive** — invisible in
the first sweep because the initialization artifact dominated everything.

**Noise threshold.** Flat ground, balanced 41 m/s jet, mixing and drag on:

| white noise | survived |
|---|---|
| 0.00 m/s | 12/12 |
| 0.15 m/s | 12/12 |
| 0.30 m/s | 12/12 |
| 0.60 m/s | 7/12 |
| 1.20 m/s | 1/12 |

A sharp threshold between 0.3 and 0.6 m/s, unmoved by any boundary-layer
setting. This is a *resolution* limit, not a physics gap.

**Initialization filter.** `src/dynamics/initialization.py`: raised-cosine
spectral lowpass, full response above 8dx, zero at or below 4dx, applied to
u, v and the theta deviation from the level mean. Order was measured, not
assumed:

| treatment | initial max\|div\| | survived |
|---|---|---|
| none | 3.90e-05 1/s | 1/12 |
| filter only | 9.93e-05 1/s | 11/12 |
| filter, then rebalance | 1.23e-05 1/s | **12/12**, max\|u\| 42.4 |

Note the middle row: filtering *raises* divergence yet survives ten hours
longer. Divergence is not the controlling variable — wavenumber content is.
Filtering changes u, v and theta separately, so it puts divergence back into a
balanced state; balancing afterwards removes it. Sub-4dx wind rms goes
0.808 -> 0.049 m/s.

**Interpretation.** The stability failure that survived nine single-candidate
patches was never in the dynamics. It was an unbalanced, over-strong,
unfiltered initial state, and the schemes added along the way were being
scored on their ability to survive an artifact. What it does *not* mean: the
boundary-layer schemes were wasted. With honest initial states both measure
as real, and drag turns out to be the stronger of the two over terrain.

What still stands open: 2500 m terrain reaches 8/12, not 12/12, even clean and
filtered. That is the next genuine question, and it is now uncontaminated.

**Status.** Kept. `surface.py`, `turbulence.py`, `initialization.py` and the
corrected `test_primitive_sigma.py` all retained. All suites green:
shallow water 8/8, boundaries 6/6, sigma 7/7, subgrid 7/7, surface 6/6,
initialization 5/5, sigma 3D core **6/6** (previously 5/6).

**For the collaboration study.** Defects introduced this session: two category
E (test design: clipped super-geostrophic jet, unfiltered white noise) and one
category B (geostrophic wind from one PGF term over sigma terrain). All three
were detected by *targeted measurement*, none by the test suite — the suite
reported them as a model failure. The stated hypothesis (drag improves the
noisy cases) did **not** survive: drag has no effect on noise at any
amplitude, and its real benefit is over terrain, which the hypothesis did not
predict. Human intervention was direction-setting and decisive: "taking a step
back and probing the error is a better idea than guess checking" is what
produced this entry rather than a tenth patch.

---

## 2026-09-02 — Prompt log and generated project structure

**Context.** The research log records what was done; nothing recorded what was
*asked*. For a study whose subject is AI-assisted building, the input side is
half the data and is the half that disappears fastest.

**Method.** `docs/PROMPT_LOG.md` — all 60 human prompts to date, verbatim, in
order, each tagged (direction / constraint / correction / methodological /
observation / administrative) with what it caused. `tools/tree.py` generates
`docs/STRUCTURE.md`: the tree shape is read from disk so it cannot drift, and
any file lacking an annotation is printed as `(unannotated)` rather than
passing silently.

**Result.** Median prompt length **11 words**. Distribution: direction 37%,
observation 22%, administrative 17%, methodological 15%, constraint 13%,
correction 5%. The four highest-leverage prompts average 19 words and none
names a technique.

**Interpretation.** The AI proposed nearly every equation, discretisation and
test in the repository, and none of the project's turns. It never proposed
starting from 2D, deferring moisture, banning HRRR from verification, capping
server usage, or probing rather than patching — the five decisions that most
determined how the work went. What this does *not* show is that the direction
was hard to produce: each of those is a short sentence. The scarce input was
knowing which sentence to say, and when.

**Status.** Kept. Append prompts as they arrive rather than in batches;
reconstructing intent afterwards is the self-report problem the methodology
section warns about.

---

## 2026-09-03 — Tall terrain: the wave was right, the aftermath was missing

**Context.** With the initialization artifacts gone, one case was still open:
2500 m terrain reached 11/12 hours clean and filtered, and no boundary-layer
setting moved it. Following the same method as last time, the failure was
located before anything was added.

**Hypothesis (stated before the runs).** The growth peaked at level k=5 of 20,
which is exactly the base of the 5-level sponge — so partial reflection off
the absorbing layer, fixable by deepening or raising the lid.

**Method and result — the hypothesis half survived.** Moving the sponge base
and watching where the growth peak went:

| sponge levels | peak growth level k | max\|du\| at 6 h | min Ri aloft |
|---|---|---|---|
| 0 | 0 (the lid) | 60.8 m/s | 0.83 |
| 5 | 5 | 36.4 m/s | 0.23 |
| 8 | 8 | 21.3 m/s | 0.94 |
| 12 | 18 (the surface) | 15.1 m/s | 1.96 |

The peak tracks the sponge base exactly, and the amplitude falls as the sponge
deepens. Reflection is real. But **survival barely moved** — 11/12 at sponge 5
and 11/12 at sponge 8 — so reflection was not what killed the run.

Raising the lid was tested too, which re-opened a prior negative result
recorded in `SigmaLevels.__init__` (that measurement had been taken on the
clipped-jet state and no longer counted as evidence). It survives re-testing:
200 hPa 11/12, 100 hPa 10-11/12, 50 hPa 9-10/12. Raising the lid is neutral to
worse. The note in the code now says so on valid data.

**Ruling out the coordinate.** A motionless isothermal atmosphere over a
mountain has no wave and no shear; anything that grows is sigma-coordinate
truncation error.

| terrain | max slope | spurious max\|u\| at 6 h | at 12 h |
|---|---|---|---|
| 1200 m | 0.0041 | 0.002 m/s | 0.003 m/s |
| 2500 m | 0.0086 | 0.004 m/s | 0.006 m/s |
| 4000 m | 0.0137 | 0.006 m/s | 0.009 m/s |

Linear in slope, and nine millimetres per second over twelve hours at 4000 m.
The coordinate is not the problem.

**The actual cause, watched hour by hour.** 2500 m, clean, filtered, sponge 8:

| hour | 1 | 3 | 6 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|---|
| max\|u\| | 41.3 | 41.3 | 41.3 | 41.3 | 41.6 | 42.1 | 42.2 | dead |
| min Ri | 11.5 | 2.07 | 0.94 | 0.33 | 0.23 | **-0.05** | **-1.15** | — |
| interfaces with Ri<0 | 0 | 0 | 0 | 0 | 0 | 1 | 18 | — |

**The wind never runs away.** It sits at 41 m/s from the first hour to the
last. What runs away is the stratification: Ri falls monotonically and goes
negative at hour 10. Ri < 0 is N² < 0 — potential temperature decreasing with
height. The mountain wave steepens as it propagates upward and **overturns**,
and the model had nothing that removes a statically unstable layer.

This is correct physics with a missing consequence. Mountain waves do break.
`turbulence.eddy_diffusivity` does treat Ri <= 0 as full-strength mixing, but
it is a diffusion capped at K = 100 m²/s, which relaxes a 600 m layer in
dz²/K = 3600 s. The wave steepens faster than an hour. Diffusion lost the race
exactly as hyperdiffusion lost the race against grid-scale noise last week —
the same failure shape, in a different scheme.

**What was added.** `src/dynamics/convection.py`: dry convective adjustment.
Wherever theta decreases with height, contiguous unstable segments are mixed
to their mass-weighted mean — neutral stratification, enthalpy conserved —
with the wind mixed over the same layers so momentum is conserved and
convective momentum transport is carried. Applied as a **post-step
adjustment**, not a tendency: an adjustment enforcing an inequality has no
meaningful time derivative, and inside the Runge-Kutta stages an intermediate
state would re-create the instability the final state must be free of.

Segment mixing replaced a first attempt at pairwise mixing, which is
conservative but converges like a diffusion: a fully inverted 20-level column
still had 0.26 K of spread after 200 sweeps. Segments settle it in one.
Conservation measured at 2.8e-16 relative for heat and 2.6e-16 for momentum.

**Prediction, written before the run, and the outcome.**

| prediction | outcome |
|---|---|
| min Ri floors near 0 instead of going negative | held: 0.09, 0.011, 0.020, 0.028, 0.013 at hours 10-14 |
| the count of Ri<0 interfaces stops growing | held: 0 for the whole run |
| the run completes 12 hours | held: reached **16** hours |
| the wind is NOT damped — convection removes overturning, not the wave | held: 41-42 m/s through hour 14 |

The fourth was the one worth stating. A scheme that bought stability by
flattening the flow would have looked identical in the first three, and that
is exactly how the first sponge implementation failed.

**Terrain rows, re-measured:**

| terrain | without convection | with convection |
|---|---|---|
| 1000 m | 12/12 | 12/12 |
| 2500 m | 11/12 | **12/12**, max\|u\| 42.4 |
| 4000 m | 6/12 | 6/12 |

**Interpretation.** The tall-terrain failure was the absence of a physical
process, not a numerical defect — the opposite of last week's finding, and
worth noting that the same debugging method produced both answers. What it
does not mean: terrain is finished. 4000 m fails at 6/12 with or without
convection, and the 2500 m run past hour 15 starts growing the wind rather
than the instability, so the mode there is different again.

Also unresolved, and now separated from the failure: sponge reflection is
measurably real, halving in amplitude between 5 and 8 levels, and it is
sitting there contaminating the upper levels whether or not it ends the run.

**Status.** Kept, on by default. All suites green: shallow water 8/8,
boundaries 6/6, sigma 7/7, subgrid 7/7, surface 6/6, initialization 5/5,
convection 5/5, sigma 3D core 6/6.

**For the collaboration study.** The stated hypothesis (sponge reflection)
was *partly* right and would have been accepted as the answer by a
patch-and-check loop — deepening the sponge does reduce the growth, visibly
and by a factor of two. Requiring it to move the survival count is what
exposed that it was the wrong cause. Category F (wrong causal hypothesis),
caught by insisting the fix predict the outcome it was proposed to explain.

---

## 2026-09-04 — 4000 m: five more candidates eliminated, and a regime boundary

**Context.** Convective adjustment took 2500 m terrain to 12/12. At 4000 m it
changed nothing — 6/12 with it and without. A scheme that fixes one case and
not the other is evidence the two cases fail differently, so the 4000 m
failure was probed rather than treated as more of the same.

**Hypotheses tested, in order, each with the measurement that settled it.**

**1. Is it the timestep?** `max|sigma_dot|` grows an order of magnitude before
the run dies (3.9e-05 → 2.9e-04), which is what a vertical-CFL violation looks
like. Halving dt:

| hour | 1 | 2 | 3 | 4 | 5 | 6 |
|---|---|---|---|---|---|---|
| max\|u\|, dt = 14.89 s | 53.86 | 54.39 | 53.04 | 54.30 | 54.98 | 54.75 |
| max\|u\|, dt = 7.45 s | 53.86 | 54.39 | 53.04 | 54.30 | 54.98 | 54.75 |
| min Ri, dt = 14.89 s | 8.190 | 1.966 | 1.101 | 0.428 | 0.032 | 0.007 |
| min Ri, dt = 7.45 s | 8.190 | 1.966 | 1.101 | 0.428 | 0.032 | 0.007 |

Identical to four significant figures for six hours. **The solution is
converged in time.** The timestep is not the cause, and `max|sigma_dot|`
growing is a symptom of the breaking, not of a CFL violation.

**2. Is the eddy-diffusivity ceiling binding?** The adjustment holds min Ri at
0.007 rather than letting it go negative, but the overturning fraction climbs
steadily (0.09% → 0.37% by hour 6) until it fires domain-wide. A breaking wave
generates turbulence, and `K_MAX` caps the diffusivity at 100 m²/s where
observed values in a breaking mountain wave are 10²–10³.

| K_MAX (m²/s) | 100 | 300 | 1000 |
|---|---|---|---|
| survived | 6/12 | 6/12 | 6/12 |

Flat. **The ceiling is innocent.** Tenth candidate eliminated by measurement.

**3. Is the initial state balanced?** Evaluating the tendencies at t = 0,
which is the check the clipped-jet episode should have had:

| terrain | max\|u₀\| | max\|v₀\| | peak acceleration | Nh/U |
|---|---|---|---|---|
| 0 m | 41.6 | 0.0 | 30.0 m/s/h | 0.00 |
| 1000 m | 41.5 | 9.3 | 31.5 m/s/h | 0.38 |
| 2500 m | 41.3 | 14.7 | 34.1 m/s/h | 0.96 |
| 4000 m | 53.8 | 33.5 | **95.2 m/s/h** | **1.19** |

The 4000 m state is measurably less balanced — three times the peak
acceleration, and 34 m/s of cross-mountain flow before a step is taken.

**The regime boundary.** Nh/U, the nondimensional mountain height, is the
parameter that orders every result in this and the previous entry:

| terrain | Nh/U | outcome |
|---|---|---|
| 1000 m | 0.38 | 12/12 with or without convection — linear wave |
| 2500 m | 0.96 | 11/12 without convection, **12/12 with** — wave at the overturning threshold |
| 4000 m | 1.19 | 6/12 regardless — blocked / breaking regime |

Nh/U ≈ 1 is the classical boundary between a mountain wave that propagates
over the obstacle and one where the low-level flow is partly blocked and the
wave breaks. The model reproduces that boundary without having been told about
it, which is a point in its favour, and it fails on the far side of it, which
is where a scheme it does not have would be needed.

**What is NOT missing.** Orographic gravity-wave drag, the standard
parameterization for this regime, would be wrong here: it parameterizes
*subgrid* orography, and this mountain is 250 km wide on a 12 km grid —
resolved by a factor of twenty. The model is explicitly simulating the wave.
Adding GWD would double-count it. Noting this because it is the scheme a
literature search suggests first, and it would have been a plausible-looking
mistake.

**Perspective on the failing case.** The highest terrain in the Northeast
domain is Mount Washington at 1917 m, and on a 12 km grid the cell-mean
elevation is well under 1500 m — Nh/U ≈ 0.5, comfortably inside the validated
envelope. **4000 m is a stress test of a mountain the domain does not
contain.** It stays on the list because it marks where the model's physics
runs out, not because the forecast needs it.

**An engineering note.** The convective adjustment initially swept the whole
domain every step and became the dominant cost of a 4000 m run — a 12-hour
integration that should take 20 minutes had not finished in 100. Compacting to
the columns that actually contain an inversion made the cost proportional to
the convection rather than to the domain: 1.7 ms on a stable state, 270 ms
when 0.3% of the domain is overturning. Worth recording because the symptom
looked like a hang, not like a performance bug.

**Status.** 4000 m open, and better bounded: not the coordinate, not the
timestep, not the sponge, not the lid, not the initialization, not the
convective adjustment, not the mixing ceiling. Nh/U > 1 with a resolved
mountain is the regime, and it is outside anything the operational domain
requires. All suites green.

**For the collaboration study.** Two of the three hypotheses this session were
mine and both were wrong (timestep, diffusivity ceiling — categories C and F).
The measurement that resolved each took under twenty minutes to design and
answered definitively. Ten candidates have now been eliminated by measurement
across this failure; the running cost of the probe-first method is roughly one
afternoon per eliminated candidate, against five patch cycles that eliminated
nothing.

---

## 2026-09-04 — A problem register, and an audit that found twenty gaps

**Context.** The research log is chronological, which is right for the study
but wrong for answering "what is broken now". A problem diagnosed across four
sessions was scattered across four entries, and the ruled-out candidates —
the expensive part of every diagnosis here — existed only as prose inside
whichever entry happened to mention them.

**Method.** `docs/PROBLEMS.md`: one entry per problem, updated in place,
across the whole project history. Statuses are OPEN, FIXED, ELIMINATED,
REVERTED and ACCEPTED — ACCEPTED existing so that "known, understood,
deliberately not fixed" is a stateable position rather than something that
looks like an oversight. `tools/problem.py` appends entries and audits the
register.

**The audit is the part that earned its keep.** The rule this project runs on
is that a fix must predict the outcome it was proposed to explain, so
`problem.py check` flags any FIXED entry without a "Confirmed by" measurement.
Run against the first draft it returned **20 issues**: thirteen fixes asserted
with no number attached, four open problems with no stated symptom, and a
numbering gap where the eliminated-candidates table was invisible to the
parser. Every one of the thirteen was a real fix — but "fixed" with no
measurement beside it is exactly the habit that produced nine failed patch
cycles, and writing them out forced the numbers to be found again.

**Result.** 45 entries: 7 open, 22 fixed, 12 eliminated, 2 reverted,
2 accepted. Register clean.

The distribution is worth noting on its own. **Twelve of the forty-five are
candidates that were investigated and were not the cause** — more than a
quarter of the register is negative results. That is the honest cost of the
probe-first method and the part that normally disappears from a repository
entirely.

**Interpretation.** The register makes one thing legible that the log did not:
of the seven open problems, only three are model defects (P-01, P-02, P-03)
and four are unbuilt work (P-04 to P-07). And of those four, P-07 — the
verification archive — is the only item in the project that gets permanently
more expensive every day it stays open, because a day not archived cannot be
obtained retroactively.

**Status.** Kept. Run `python tools/problem.py check` before a commit.

---

## 2026-09-04 — The sigma core is now reachable from real data

**Context.** The terrain target was agreed at 2 km, and the model does 2500 m
at 12/12, so P-01 moved from OPEN to ACCEPTED and the register's largest
remaining items were the two that had nothing to do with physics: the driver
still built a `Primitive3D` on pressure levels. Every result measured since
the coordinate change was unreachable from a real forecast — the only core
real data could run was the one that diverges in two to three hours.

**Method.** `src/dynamics/interpolate.py`, and a rewrite of the driver.

The conversion is three steps and the order matters. Terrain height gives
surface pressure by finding the pressure at which the analysis geopotential
height equals the terrain — an interpolation rather than a hydrostatic guess,
so it inherits the analysis's own stratification. Surface pressure gives the
target pressure of every sigma level. The analysis columns are then
interpolated to those pressures in **log(p)**, which matters more than it
sounds: the level spacing runs from 25 hPa near the ground to 50 hPa aloft,
and a field is far more nearly linear in log(p) than in p.

**Extrapolation was the part with a trap in it.** Sigma levels near the ground
over low terrain sit at pressures *below* the analysis's lowest level — 1000
hPa is about 100 m above sea level, not the surface — so something has to be
said about the layer beneath. Theta follows the lapse rate of the lowest two
levels rather than being held constant, because holding it constant makes the
near-surface layer exactly neutral, which the convective adjustment then reads
as marginal everywhere on step one. Wind is held constant, because
extrapolating a shear downward produces surface winds the drag scheme fights.

**Result.**

| test | result |
|---|---|
| field linear in log(p) reproduced | 0.00e+00 error |
| source levels recovered | 3.6e-15 |
| surface pressure vs standard atmosphere, 0–2500 m terrain | **7.6 Pa** |
| sea-level terrain extrapolated below 1000 hPa | 1013.3 hPa vs 1013.25 |
| converted analysis integrated 6 h | held |
| `test_interpolate.py` | 7/7 |
| `test_forecast.py`, rewritten for sigma | 11/11 |

**Two defects found on the way, both by tests written to have a known
answer.**

The bracket search in `surface_pressure_from_heights` had the height ordering
backwards. After sorting by descending pressure, index 0 is the highest
pressure and therefore the *lowest* height — heights increase with index — and
I had written the comparison the other way. Every column above the lowest
analysis level stayed pinned at that level's pressure: a **253 hPa** error
over 2500 m of terrain. Caught only because the test compares against a
standard atmosphere, where the answer is known in closed form; a
self-consistency check would have passed. Category D, and the fourth time an
ordering convention has been the defect.

The second was in a test rather than the code, and is worth recording because
it hid a real hazard. The relaxation test reported "edge moved +0 Pa" — it had
handed the model the same array it later compared against, and `apply` updates
in place. The test was measuring nothing. The driver had the same aliasing
hazard: assigning the analysis arrays to the model directly would let the
first relaxation step quietly rewrite the driving data. Both now assign
copies.

**A new problem, opened rather than absorbed.** A converted analysis started
at rest over a 1500 m mountain develops **9.1 m/s** of spurious wind in six
hours, where a state hydrostatically consistent with the model's own
discretisation develops 0.004 m/s. The interpolated theta reproduces the
analysis's stratification but not the model's hydrostatic integral over its
sigma layers, so the geopotential carries a gradient no wind balances. Against
a 20–40 m/s analysis wind that is a 25–30% error injected before the first
step, and it is now P-46 rather than a footnote.

**Interpretation.** The physics measured over the last week is now reachable
from real data, which it was not this morning. What this does **not** mean is
that a real forecast has been run: the surface-field GRIB search has never
touched HRRR, and on this project's record (P-20 through P-22) that is where
the next defect will be.

**Status.** Kept. P-01 accepted at the 2 km target, P-04 and P-05 closed, P-46
opened. Register: 46 entries, 5 open. All suites green, including
`test_forecast.py` 11/11 and `test_interpolate.py` 7/7.

---

## 2026-09-04 (later) — P-46 dies: three hypotheses, and the test was the defect

**Context.** P-46 was opened this morning on a real measurement: a converted
analysis started at rest over a 1500 m mountain developed 9.1 m/s of spurious
wind in six hours, where a state the model built itself develops 0.004 m/s.
Against a 20-40 m/s analysis wind that is a 25-30% error injected before the
first step, so it was the highest-value open item.

**Hypothesis 1, stated first: geopotential mismatch.** The interpolated theta
reproduces the analysis's stratification but not the model's discrete
hydrostatic integral, so the model's geopotential differs from the analysis's
by an amount that varies horizontally — and a horizontally varying
geopotential error is a pressure-gradient force with nothing balancing it.

The measurement supported it, at first. The error's horizontal spread grew
upward from 3 m at the ground to **140 m at the lid**, exactly the shape a
column-by-column integration error would have.

`hydrostatic_geopotential` is triangular, so it inverts exactly:

    T[-1] = (phi[-1] - phi_s) / (R ln(p_s/p[-1]))
    T[k]  = 2 (phi[k] - phi[k+1]) / (R ln(p[k+1]/p[k])) - T[k+1]

| | geopotential spread | acceleration at rest |
|---|---|---|
| interpolated theta | 140.08 m | 2.70 m/s per hour |
| exact hydrostatic inversion | **0.00 m** | **2.67 m/s per hour** |

**The error went to zero and the acceleration did not move.** The hypothesis
was wrong. The inversion also produced a statically unstable profile with an
8.5 K sawtooth and went NaN in six hours — that `- T[k+1]` makes the inverse
an alternating recursion, so an error at one level flips sign and persists
upward. Kept in `interpolate.py` as `hydrostatic_theta`, unused, with the
warning in its docstring, because the next person to have this idea should be
able to read why it does not work.

**Hypothesis 2: small-scale structure from the interpolation.** The sigma PGF
is a difference of two large terms that cancel only when theta is smooth.

| treatment | spurious wind after 6 h |
|---|---|
| interpolated, raw | 9.06 m/s |
| + horizontal spectral filter | 9.06 m/s |
| + one pass vertical smoothing | 9.06 m/s |
| + three passes | 9.03 m/s |

Nothing. Two hypotheses, both wrong, both eliminated by measurement rather
than argument.

**Hypothesis 3: the test.** The tendency breakdown had been sitting in the
output the whole time — the acceleration was 2.70 m/s/h in du/dt and
**35.45 in dv/dt**. The initial state has a meridional temperature gradient
and no wind. That is not a balanced state the model is corrupting; it is an
unbalanced state the model is correctly adjusting toward balance. 9 m/s over
six hours is geostrophic adjustment.

The decisive control: **on flat ground the same setup drifts 8.98 m/s.** There
is no terrain, no conversion over terrain, and almost all of the drift is
still there.

**And the synthetic analysis was itself inconsistent.** It perturbed
temperature by −1.5 K and geopotential height by −45 m of cos(k_y·y),
independently. Those two are not in hydrostatic balance with each other; a
real analysis is. Rebuilding the heights as the hydrostatic integral of the
temperatures:

| | inconsistent analysis | self-consistent |
|---|---|---|
| geopotential spread, 1500 m terrain | 140 m | **3.24 m** |
| implied balanced wind | 100.4 m/s | 61.9 m/s |

**What the conversion actually costs**, measured on a consistent analysis:

| terrain | geopotential error (horizontal spread) | drift in 6 h at rest |
|---|---|---|
| flat | **0.01 m** | 8.98 m/s (all adjustment) |
| 1500 m | 3.24 m | 10.19 m/s |
| 2500 m | 4.19 m | 11.54 m/s (**+2.56** over flat) |

**Interpretation.** P-46 was not a defect. The conversion is accurate to a
centimetre on flat ground and four metres of geopotential height over 2500 m
terrain, and terrain adds 2.6 m/s to a 9 m/s adjustment that a rest start
demands on its own. Category E, the fourth test-design error of the project —
and the first one found by the AI rather than by a human noticing an anomaly,
which is worth recording given that the score on that was previously 3-1
against.

What it does *not* mean: the two rejected hypotheses were wasted. Hypothesis 1
produced the exact-inversion operator and the measurement showing why exact is
not the same as usable, and hypothesis 2 established that the conversion
introduces no small-scale structure worth filtering. Both are now in the
eliminated table rather than available to be re-proposed.

**The tests were rewritten to measure the right thing.** The old decisive test
asserted "spurious wind under 25 m/s" from a rest start, which passes for the
wrong reason. It is replaced by two: one on the geopotential error, which is
what the conversion is responsible for, and one measuring terrain's
contribution against a flat control rather than against zero.
`test_interpolate.py` 8/8.

**Status.** P-46 eliminated. Register: 46 entries, 4 open — and every one of
the four now needs either a network or an idea, not a measurement.

---

## 2026-09-04 (evening) — The archive machinery, and a 74 K trap in the old operator

**Context.** P-07 is the only item in the register that gets permanently more
expensive every day it stays open. Observations stay downloadable for years;
the forecast that was valid for them is only makeable on the day. The driver
port closed the blocker this morning, so this is the piece that turns a
running model into evidence.

**The design decision, stated before the code.** The archive stores raw
observations verbatim and compressed, written **before** any parsing, QC or
matching is attempted, with the forecast copied beside them. Matched pairs are
derived. If the observation operator changes — and it will, because the
elevation correction is a standard 6.5 K/km lapse rate that is wrong on
exactly the calm clear nights when it is largest — every match can be
recomputed from the raw payload and the forecast. A failure anywhere
downstream must never cost the irreplaceable part.

**A 74 K trap found on the way.** `GridInterpolator.at_observation` falls back
to `field3d[0]` when an observation has no pressure — which is every surface
observation, and index 0 in this project is the **model lid**. An ASOS
thermometer would have been scored against the 200 hPa field.

| | value returned for a 2 m thermometer |
|---|---|
| sigma operator (lowest level) | 286.6 K |
| base class (`field3d[0]`, the lid) | 212.8 K |

That is not a crash and not an obviously wrong number in isolation — it is a
temperature. It would have appeared as a catastrophic, uniform cold bias and
been read as a model failure. `SigmaInterpolator` overrides the method rather
than extending it, and the test that catches it compares the two operators
directly so the trap cannot come back.

The operator also had to become column-aware: in sigma coordinates there is no
shared 1D array of level pressures, so `vertical()` had nothing to interpolate
against. Measured: the lowest level sits at 985 hPa over flat ground and
821 hPa over 1500 m terrain, in the same run.

**A defect of my own, and the failure mode is worth naming.** `verify.py`
first used variable names of its own — `"temperature"`, `"u_wind"`,
`"v_wind"` — while the fetchers, `config.CHANNELS` and `RANGE_LIMITS` all use
the GRIB-style `TMP`, `UGRD`, `VGRD`. Every observation fell through to
"variable not verified". The archive came out **empty while every other check
passed**: raw observations stored, forecast copied, metadata written, no
exception raised anywhere. Inventing a second vocabulary for something that
already had one produces a pipeline that reports success and archives nothing.
Category A in spirit — an interface assumption — though the interface was
internal.

**Result.** `test_sigma_operator.py` 7/7, `test_verify.py` 7/7. The archiver
round-trips raw observations byte-for-byte, records lead time on every pair
(skill decay is not recoverable from a pair alone), records the size of every
elevation correction (so a bias caused by the operator can be told apart from
a bias in the model), refuses to reach the network under `--report-only`,
refuses a pre-sigma forecast file, and — tested explicitly — **adds nothing on
a second run of the same day**. A daily job will be run twice; without that,
every score would be silently weighted by how often someone re-ran the script.

`tools/daily.sh` runs ingest → forecast → verify from cron: lock file (two
12-hour forecasts competing for the same cores is what the 50% ceiling
exists to prevent), dated logs kept because a 4 a.m. failure is only
diagnosable from what it wrote at the time, first-failure exit code, and
verification attempted even when the forecast step failed — a forecast that
diverged at hour 8 still produced eight hours worth archiving.

**Interpretation.** P-07 stays OPEN, and should. The machinery exists and is
tested, but the archive has no data in it, and nothing here has met the live
service. On this project's record that is where the next defect is: every
interface defect so far (P-20 to P-24) passed a full offline suite and
appeared on first contact.

**Status.** Machinery kept. Register: 46 entries, 4 open.

---

## 2026-09-04 — Configuration change: model reasoning effort raised to high

**Recorded because the study's subject is the AI, not only the model.** A
change to how the collaborator is configured is a change to the instrument,
and an instrument change part-way through an uncontrolled n=1 study has to be
in the record or every before/after comparison in this log is quietly
confounded.

**What changed.** The reasoning effort setting was raised to **high** by the
human collaborator on 2026-09-04. Session model identifier: `claude-opus-5`.
The serving model can differ from the configured one and can change
mid-session, so that identifier is what was configured, not a guarantee of
what answered any particular turn.

**What was in progress at the switch.** The sigma core complete through
convective adjustment; the driver ported to sigma; the verification archiver
built; P-46 eliminated. Register at 46 entries, 4 open.

**What this does NOT license.** Any claim that work after this point is better
than work before it. There is no control, no repeated trial, and no blind
comparison — the tasks on either side of the switch are different tasks. The
error taxonomy in `docs/AI_COLLABORATION.md` is the only quantitative record,
and it is small enough that a difference of two or three defects is noise.

**What it might reasonably be compared on, later, with all the above caveats.**
Defects per session by category; how often a stated hypothesis survives
measurement (currently roughly half); and the ratio of AI-proposed technique
to human-proposed direction, which `docs/PROMPT_LOG.md` tracks.

Anyone reading this log as evidence should treat 2026-09-04 as a seam and not
pool across it without saying so.

**Status.** Recorded. No code change.

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

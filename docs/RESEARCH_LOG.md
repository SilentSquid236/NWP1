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
| survived | ~~6/12~~ 6/12 | ~~6/12~~ **8/12** | ~~6/12~~ **8/12** |

~~Flat. **The ceiling is innocent.** Tenth candidate eliminated by
measurement.~~

**WITHDRAWN 2026-09-10.** The ladder was not varying the ceiling. It assigned
`turbulence.K_MAX` at runtime, which the mixing scheme never reads, because
the value is bound into `vertical_mixing`'s signature at import (P-51). Three
identical survival counts were the *symptom* of that, and were read as a
result. Re-run through the constructor the second and third rungs are 8/12,
so the ceiling is worth two forecast hours and this candidate was never
eliminated. The struck-through numbers are left in place because a wrong
conclusion drawn from a broken measurement is part of the record; see P-40 and
the 2026-09-10 entry.

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

## 2026-09-05 — The sponge does not absorb weather selectively, and never did

**Context.** P-02: growth over terrain peaks exactly at the sponge's lower
edge and moves when the edge moves. Real, measured, and not what ends a run.
The stated reason it could not simply be deepened was P-16 — an early sponge
that relaxed toward the horizontal mean and flattened a jet.

**Hypothesis, stated before the runs.** Gravity waves are divergent; balanced
flow is rotational; `remove_divergence_spectral` splits them exactly. Damp
only the divergent component and the layer can be deep without flattening
anything.

**Result: the hypothesis failed, and took the premise with it.**

| sponge | kind | jet drift 24 h | reflection max\|du\| 6 h | peak level |
|---|---|---|---|---|
| 5 | plain | 0.24% | 36.39 | 5 |
| 8 | plain | 0.22% | 21.26 | 8 |
| 12 | plain | 0.22% | 15.08 | 18 |
| 8 | divergent only | 0.22% | **55.47** | 0 (the lid) |
| 12 | divergent only | 0.22% | 53.74 | 0 |

Divergent-only absorbs *less* than plain: a mountain wave is not purely
divergent, and its rotational part reflected off the lid untouched.

But the same table shows a **12-level plain sponge does not flatten the jet at
all** — 0.22% drift, identical to 5. The constraint inherited from P-16 no
longer applied: that sponge relaxed toward the horizontal mean, this one
relaxes toward a frozen reference, and for a steady jet the jet *is* the
reference. The reason the sponge had to stay shallow had been assumed for
four days and was false.

**Which exposed the real question.** A frozen reference is exactly right for a
jet that does not change and exactly wrong for one that does. A steady-jet
test cannot see that. So: baroclinic development against sponge depth.

| sponge | 0 | 2 | 3 | 5 | 8 | 12 |
|---|---|---|---|---|---|---|
| eddy energy x/day | **3.18** | 0.93 | 0.76 | 0.85 | 0.96 | 1.02 |

**The default configuration turns growth into decay.** Not at twelve levels —
at *two*.

**A metric trap, caught before it was believed.** At sponge=2 the ratio
collapses to 0.93 while max|v| is 12.8 against 12.7 with no sponge. Those two
cannot both mean suppression: a wave that grew fast and saturated during day 1
reports a ratio near 1 and looks identical to one that never grew. So the
whole curve, not two points on it:

| sponge | 6 h | 12 h | 24 h | 36 h | 48 h | |
|---|---|---|---|---|---|---|
| 0 | 6.5e+03 | 4.2e+03 | 3.7e+03 | 4.2e+03 | **1.2e+04** | grows |
| 2 | 3.6e+03 | 2.9e+03 | 1.7e+03 | 1.3e+03 | 1.6e+03 | decays |
| 3 | 2.8e+03 | 2.1e+03 | 1.0e+03 | 8.6e+02 | 7.9e+02 | decays |
| 5 | 2.0e+03 | 1.3e+03 | 8.1e+02 | 7.1e+02 | 6.9e+02 | decays |

It is suppression. Monotonic decay with any sponge, growth without. The curve
also shows the no-sponge case oscillating by a factor of two between samples,
which means the day-2-over-day-1 ratio was partly luck of sampling — the same
diagnostic that rejected divergence damping, used for four days without anyone
checking its blind spot.

**Diagnosis.** The lid is at 200 hPa, so the sponge sits at **301 hPa** — in
the upper troposphere. Baroclinic instability is a coupled mode between an
upper and a lower wave; damping either end stops it, whatever it is damped
toward. Rayleigh damping cannot distinguish a mountain wave from a baroclinic
wave here because they occupy the same levels.

**Four more attempts, all measured, none works.** Rate from 15 min to 6 h:
0.85, 0.78, 0.82. Running low-pass reference instead of frozen: 0.80 at 6 h,
0.61 at 1 h. Lid at 100 and 50 hPa: 0.97 and 1.89 — development recovers as
predicted, but reflection climbs to 53 and 60, which is the no-sponge value.
More levels to buy both: 26 levels diverged, 30 levels gave 2.45 development
and 56.4 reflection — still not absorbing.

The pattern is the same in every row: **whenever the sponge absorbs
(reflection 21-36) development collapses; whenever development survives
(1.89-2.45) the sponge is not absorbing.**

**Interpretation.** This is not a tuning problem and there is no setting that
resolves it. The fix is a radiative upper boundary condition, which acts on
wave flux rather than amplitude and therefore has no such trade-off. Not
attempted.

What it does *not* mean: that current 12-hour forecasts are invalid. The
curves diverge most after 24 h; at 12 h sponge=3 holds 2.1e+03 against
4.2e+03. Affected, not invalidated — and the first thing to fix before
extending past a day.

**Change made, then reverted the same day.** The default was cut from 5 to 3
on the reasoning that a shallower layer must do less damage. It does not — by
the 48 h / 6 h ratio, three levels decays *faster* than five (0.28 against
0.34), which was visible in the table I had just written and did not read
carefully enough. The change regressed the decisive noisy case from 12/12 to
11/12 and broke the Ekman spiral test, so it went back to 5. Both suites
recovered on revert (6/6 each).

Depth is not the mechanism, so trading depth buys nothing. The setting stays
at 5 with the measurement in the code, so the next person to reach for it can
see that it has already been tried.

**Status.** P-49 opened, P-02 marked as the symptom of it. Register: 49
entries, 5 open.

**For the collaboration study.** Three failures worth counting. The stated
hypothesis (divergent-only damping) was wrong — category F. And a constraint
inherited from an earlier fix (P-16, "a deep sponge flattens the jet") was
carried for four days without being re-measured after the thing it described
had changed. That second one is not in the taxonomy yet: not a wrong
hypothesis, but a *stale* one — a fact that was true when recorded and
silently expired. Worth its own category if it recurs.

The third is the cleanest: I changed a default on a reading contradicted by
the table in the same commit. The measurement was correct, present, and
misread — and only the regression suite caught it. Category F, but a
sub-species worth naming: not a hypothesis that survived measurement and
failed later, one that the measurement already refuted at the moment it was
made.

---

## 2026-09-05 (later) — A radiative upper boundary: the mechanism works, the terrain case does not

**Context.** P-49 established that no Rayleigh sponge can separate a mountain
wave from a baroclinic wave, because amplitude is not what distinguishes them.
The standard answer is a boundary that passes vertical wave flux instead of
damping amplitude.

**What was built.** `src/dynamics/radiation.py` — the hydrostatic
Klemp-Durran / Bougeault condition, for each horizontal wavenumber:

    w(k) = |k| phi'(k) / N

converted to the mass flux this model's continuity equation wants, since at
sigma = 0 the pressure is constant and omega_top = pi * sigmadot_top:

    F(k) = -rho_top g (|k| / N) phi'(k)      [Pa/s]

`continuity` gained a `top_flux` argument. The algebra needed care: the
constant added to the partial integral must be F itself, not F(1-s), and the
two end values are what check it — sigma_dot comes out as F at the lid and
exactly 0 at the ground (measured 3.95e-20). Mass may leave through the top;
it may never leak through the surface.

**The sign was measured, not argued.** Get it backwards and the boundary is a
wave SOURCE pumping energy in at exactly the rate it should let it out. With
sign +1 the run goes non-finite within four hours; with sign -1 it does not.
That is the test, and it is in the suite.

**Result on the thing that mattered.** Eddy kinetic energy of a growing
baroclinic wave, 48 h, ratio of 48 h to 6 h:

| configuration | 6 h | 24 h | 48 h | ratio |
|---|---|---|---|---|
| sponge 5, rigid lid | 2.04e+03 | 8.06e+02 | 6.85e+02 | **0.34** decays |
| no sponge, rigid lid | 6.51e+03 | 3.72e+03 | 1.19e+04 | 1.82 |
| no sponge, **radiative** | 1.58e+04 | 6.68e+03 | 3.95e+04 | **2.50** |
| sponge 5, radiative | 1.48e+03 | 5.74e+02 | 5.48e+02 | 0.37 |

Better than a rigid lid, not merely unharmed — an amplitude e-folding of about
**1.8 days** against 6.7 days for the rigid lid, and 1-3 days is what
baroclinic waves actually do. The last row is worth noting too: with the
sponge still on, the sponge dominates and the boundary buys nothing. They are
not additive.

**And it fails the case the sponge was there for.** 2500 m terrain, 12 h:

| configuration | survived |
|---|---|
| sponge 5, rigid | 12/12 |
| no sponge, rigid | 9/12 |
| no sponge, radiative | **3/12** |

**Two candidate causes, both measured, neither is it.**

*The timestep.* The flux is explicit and acts on the thinnest layer in the
column, so CFL was the obvious suspect. dt and dt/2 give 3/12 and 3/12 —
identical. Not the timestep.

*Radiating a steady anomaly.* Over a mountain the top-level geopotential
perturbation is dominated by the terrain's own hydrostatic imprint, which is
balanced and does not propagate. Radiating it pumps mass out of those columns
continuously. So the condition now acts on the deviation from a 3-hour running
low-pass — only the transient part radiates. That **fixed the resting case**
(a balanced atmosphere over 1500 m terrain now drifts 1.4e-05 in mass and
1.07 m/s in wind over 3 h) and improved development from 1.82 to 2.50. Terrain
survival stayed at 3/12.

**Interpretation.** The mechanism is right and the measurement supporting it is
the strongest evidence in this log for any single change: the sponge decays
weather at 0.34 and this grows it at 2.50, with a physically correct e-folding
time. But it is not usable yet, and the honest statement is that the terrain
failure has an unknown cause, not a suspected one — the two obvious
explanations were tested and eliminated.

**Status.** Kept, **off by default**. The sponge remains production with its
cost recorded. P-50 opened for the terrain failure. `test_radiation.py` 7/7,
including the sign measurement, the resting-atmosphere leak test, and a
development test that requires growth where the sponge decays (36 h / 18 h:
sponge 0.75, radiative 1.91).

**A test-design note.** The first version of the development test compared
24 h to 12 h and reported 0.61 against 0.73 — no discrimination at all, because
both configurations are still shedding the initial transient at 12 h. The
curves only separate once growth takes over. An assertion window chosen for
speed rather than from the measured curve is not a test.

The same mistake appeared once more in this session, in the other direction:
the first radiative test asserted that the disturbance AT THE LID should be
smaller, and it was larger (32.2 against 27.2 m/s). That is the expected
behaviour — a rigid lid holds the wave still and a transparent one lets it
through, so the top level is *more* active, not less. A reflected wave shows up
below, which is where the measurement belongs.

---

## 2026-09-06 — Three causes behind one failure, and the analysis that found the third

**Context.** P-50: the radiative lid dies over 2500 m terrain at 3/12 hours.
The timestep and the terrain's steady imprint had both been eliminated the day
before, so the cause was genuinely unknown.

**The probe was spatial, because the suspect was.** Another module's docstring
already said it: `remove_divergence_spectral` notes that the FFT assumes
periodicity, the domain is not periodic, and the outermost cells carry the
error — but it gets away with it because the lateral relaxation overwrites
them. The radiation flux has no such protection; it goes straight into
prognostic surface pressure. So the probe recorded flux at the edges against
flux in the interior.

| hour | \|F\| edge / interior |
|---|---|
| 1 | 1.49 |
| 2 | 2.28 |
| 3 | **42.32** |

Interior flux was flat at 0.67-0.79 Pa/s the whole time. The boundary was not
radiating a wave; it was amplifying its own transform error. Windowing before
the transform and tapering after it fixed that — edge ratio down to 0.13.

**Which exposed a second problem underneath.** With the edges quiet, the
INTERIOR flux started growing: 1.25 → 4.59 Pa/s in an hour, death at hour 3 —
sooner than before. The transfer is proportional to \|k\|, so the shortest waves
get the most flux. A raised-cosine cutoff below 8 grid cells slowed it to 2.17
but did not stop it. That cutoff is worth keeping on its own terms: the
hydrostatic radiation condition is only valid where \|k\| << N/U, which fails
long before the grid scale, so it removes exactly the wavenumbers the
condition was never derived for.

**Then I stopped tightening knobs and did the loop analysis.** Two failed
adjustments in a row is the signal the project's own methodology says to stop
on.

    F  ->  dpi/dt  ->  pi  ->  phi_top  ->  F

phi_top is the hydrostatic integral from the ground up, so it moves when pi
moves, by about R T / p_s per pascal — 0.94 m²/s² per Pa here. With rho g / N
about 164 and \|k\| about 5e-5 for a 120 km wave, the loop gain is 154 \|k\| per
second: **an e-folding of about 130 seconds.** And the sign that radiates
waves correctly is the sign that makes this loop grow, so it was never
resolvable by choosing a sign. That number explains the growth rates observed
and nothing else did.

**The first fix for it was wrong, in an instructive way.** Subtracting the
column's hydrostatic response to pi made things worse — flux 25 Pa/s, death an
hour earlier. Over 2500 m terrain pi varies by 27000 Pa because of the
*mountain*, so pi' is ±13000 Pa against a phi' of about 150 m²/s². That
subtraction does not open the loop, it injects a terrain-shaped signal two
orders of magnitude larger than the wave. The loop runs through *changes* in
pi, so it is the change that has to be removed — pi against its own running
low-pass, matching the treatment phi_top already had.

**Result.** Flux bounded (1.52, 2.22, 1.11 Pa/s), wind steady at 39 m/s for
three hours, and the run still ends at hour 4. Three causes found and fixed, a
fourth remains, and nothing in the recorded diagnostics is running away before
the end.

**What improved anyway.** Development, the thing the boundary exists for, got
better with every fix: the suite's 36 h / 18 h eddy-energy ratio went 1.91 →
**3.12**, against the sponge's 0.75. All suites green — radiation 7/7, sigma
7/7, sigma 3D core 6/6.

**Interpretation.** A single symptom with three independent causes stacked
behind it, each of which had to be removed before the next was visible. Two
were found by measurement and one by arithmetic — and the arithmetic one was
the only one that could not have been found by trying settings, because no
setting fixes a loop whose two requirements have opposite signs.

**Status.** Kept, off by default. P-50 stays open with a much sharper
description than it had.

**For the collaboration study.** The useful entry here is not a defect but a
decision: after two failed adjustments in a row, switching from "try the next
setting" to "write down the feedback loop and estimate its gain" is what
produced the third cause. That is the same intervention the human made on
2026-09-02 ("probing the error is a better idea then guess checking"), applied
without being prompted this time. Worth noting for whether the habit
transfers.

---

## 2026-09-22 — Packaged for Claude Science, and two seams in the study

**Context.** The project is moving to Claude Science, Anthropic's research
workbench, where each project keeps its own memory and skills. A session there
starts with none of this conversation. So the question was not "which files"
but "what does a fresh session need in order not to repeat the mistakes this
record documents."

**What was built.**

- `CLAUDE.md` — the brief a new session reads first: what the model is, the
  five non-negotiable constraints, the method, current state, open problems,
  where things live, what must pass before anything is called done, and the
  specific things a new session gets wrong (index 0 is the lid; K_MAX is
  instance state; the sponge relaxes toward a frozen reference).
- Three skills under `skills/`: `nwp-debug` (probe-first diagnosis),
  `nwp-record-session` (the bookkeeping), `nwp-sync` (moving work between
  GitHub, the desktop and the server).
- `CLAUDE_SCIENCE.md` — the import steps, and what does not come across.

**Why skills, specifically.** The learning log's main finding was that most
lessons in this project survive only as habits — nothing enforces them, and a
habit is exactly what a fresh session lacks. A skill is a habit written down
where the next session will load it. The three skills carry the lessons with
the most callbacks: probe-don't-guess (L1), a fix must predict (L2), suspect
the test (L3), identical results are a bug (L12), and the transfer failures
that kept recurring silently.

**Every claim in `CLAUDE.md` was checked against the code before it went in**:
K_MAX 200 and Ri_crit 0.25 in `turbulence.py`; sponge 5 levels and radiative
lid off by default in `primitive_sigma.py`; the domain bounds in `config.py`.
A brief with a wrong fact in it is worse than no brief, because the new session
has no way to know.

**What could not be confirmed.** How Claude Science picks up an existing
repository or a `CLAUDE.md` file is not documented in the announcement or the
review checked. The import guide therefore makes its first step a test: ask the
new session to list the five constraints, and paste `CLAUDE.md` in by hand if
it cannot.

**A constraint the move could have broken.** Claude Science can run on a Linux
box, and the natural idea is to run it on the Xeon where the data is. That
server forbids installing packages, which covers this app. The guide says so
explicitly: install it on the Windows desktop.

**Two seams in the study, both recorded.** The session model was changed to
`claude-opus-5-5` this turn, and the environment is about to change. Both are
in the "Instrument changes" table in `docs/AI_COLLABORATION.md`. Counts before
and after should not be pooled without saying so.

**The package is built on the unmerged P-40 branch** (`6d26e1e`) plus the
token-ledger commit, so it carries the K_MAX default of 200 and P-52. Merging
`p40/ceiling-ladder` into `main` first keeps GitHub and the package agreeing.

**Status.** Kept.

---

## 2026-09-22 — First session in Claude Science; git reaches the server

**Context.** The first session in the new environment (prompts 83–84). The
import guide made its step 2 a test: if the new session cannot list the five
constraints in `CLAUDE.md`, it has not picked up the context.

**Hypothesis.** Stated in `CLAUDE_SCIENCE.md` before the move: the app might or
might not load `CLAUDE.md` by itself, and the five-constraint test would show
which.

**Method.** Prompt 83 asked the session to read `claude.md`. It checked project
memory, then looked for the file; it listed the constraints back from the file
it found. The model id was read from the session runtime rather than assumed.

**Result.**

| check | outcome |
|---|---|
| `CLAUDE.md` loaded without being asked | **no** — project memory was empty and nothing had been read |
| file at the path named | no `claude.md` at `Desktop\NWP\`; found at `NWP1\CLAUDE.md` (Windows paths are case-insensitive, so this was location, not case) |
| five constraints listed back | all five |
| model | `claude-opus-5-5`, matching the human's report |
| skills | three imported; `nwp-sync` revised first (below) |
| research record | read in place from the granted folder, not copied |
| GitHub `main` | still `b97bc06` in the desktop clone; no `CLAUDE.md` or `skills/` there; `package/claude-science` unmerged |

**git on the server (prompt 84).** The human reported git is now installed on
the server. Constraint 2 in `CLAUDE.md` was rewritten rather than deleted, and
the text that said git "cannot be" installed was corrected in `tools/pull.sh`
and `tools/stale.py`. `skills/nwp-sync` gained an in-place conversion of the
`pull.sh` copy to a git checkout. Writing it surfaced one new hazard: `data/` is
in `.gitignore`, which keeps git from overwriting it but also means
`git clean -x`/`-X` delete it and `git stash --all` removes it. "Ignored" had
been doing the work of "protected", and it is not the same thing. The
conversion recipe has **not been run** on the server; its status is that of
P-06, written and not exercised. The server's git version was not reported.

A side effect worth measuring later: `tools/stale.py` prefers commit dates, and
on 2026-09-10 mtimes produced 12 false hits against 0 from git dates. Once the
server copy is a checkout, the stale check there should stop reporting a
whole transfer batch as changed.

**Interpretation.** The step-2 test was needed; without it the session would
have started with none of the constraints. It does not show the brief was
inadequate: once read, it was sufficient. The constraints are now in project
memory, which gives a prediction that can fail: **the next Claude Science
session can state the five constraints before reading any file.**

**Found in the record.** The Aggregate table in `docs/PROMPT_LOG.md` is stale.
It lists DIR 22, OBS 13, ADM 10, MET 9, CON 8, COR 3. A regex count of the tag
column over the 76 numbered rows present before this session (prompts to 82)
gives DIR 28, ADM 15, MET 14, OBS 10, CON 7, COR 3. The difference in CON
suggests either a tag the regex missed or a miscount in the table. It was left
unedited for the human to resolve, because it feeds the study's main finding.

**A defect the move exposed.** Running the pre-done checks on Windows for the
first time, `tools/tree.py` and `tools/manifest.py` wrote `STRUCTURE.md` and
`MANIFEST.txt` with CRLF line endings. `.gitattributes` requires LF
everywhere; on Linux text mode writes LF by construction, so it had not
shown. Both now write LF explicitly (`newline="\n"`; `write_bytes` in
`manifest.py`, because `write_text(newline=)` needs Python 3.10 and the
server's version was not checked). After the fix no `.py`, `.md`, `.sh`,
`.txt` or `.csv` file in the tree contains CRLF. It is category A in the
taxonomy, an unobserved platform convention, found by a check rather than by
failing.

**For the collaboration study.** The environment seam and the model seam fall
on the same day but are not the same event: prompt 82's session already ran on
`claude-opus-5-5` in the old environment, so exactly one session separates
them. Prompt 84 is a CON that *removes* a constraint. It undoes prompt 15
("sadly git is not on the server"), which killed the git-based transfer plan
and led, via prompt 68, to `pull.sh`.

**Addendum, same day: where the server copy actually is.** The first `git
pull` was run in `/data5/pierce/Data5/NWP`, a path the AI inferred from the
data root in `config.py`; it failed (`not a git repository`, and no
`tools/manifest.py`). A read-only survey (`find`, `ls`, `~/.bashrc`, `crontab`,
`manifest.py --check`) found:

| item | finding |
|---|---|
| project root | `/data5/pierce/NWP` — top level, files dated 2026-09-22, manifest matches; **already has `.git`** (14:53). The nested `NWP_Deployment_Package/` is also complete (all four probe files) but dated 2026-09-04 |
| data root | `NWP_DATA_ROOT=/data5/pierce/NWP/NWP_Deployment_Package/data` (`~/.bashrc` line 44) |
| nested copies | `NWP_Deployment_Package/` (own `.git`, holds `data/`), `NWP1-main/`, `nwp.tar.gz` |
| other strays | `/data5/pierce/Data5/{NWP,NWP1-main,data,src/data}`, `/data5/pierce/config.py` |
| `manifest.py --check` in the root | 0 missing, 0 differing, 189 extra (the nested copies) |
| crontab | none — `tools/daily.sh` has never been scheduled |

The documented server data root (`/data5/pierce/Data5/NWP/data`, in
`README.md` and `config.py`) did not exist; both now give the real value. The
finding that matters most: the verification archive sits inside a directory
that every earlier lesson about nesting (P-24) would mark for deletion. The
AI's inference of the root was wrong, and it was a guess stated as "most
likely" rather than checked; the survey should have come first. That is L1
(probe, don't guess) applied to a path instead of a model.

**Server git state** (read-only: `git --version`, `remote -v`, `branch -vv`,
`log`, `status`): git 2.52.0; `origin` is the GitHub repository; `main` at
`44068f2` tracking `origin/main`, with PRs #2 (`package/claude-science`) and #3
(`p40/ceiling-ladder`) merged; working tree clean except the untracked
`NWP1-main/`, `NWP_Deployment_Package/` and `nwp.tar.gz`. So no conversion was
needed and `git pull --ff-only` works from the root.

**Correction to the table above:** "`package/claude-science` unmerged" was
wrong. It was read from the desktop clone, which had not been fetched since
`b97bc06`; GitHub had merged it. A local clone is a snapshot of the remote as
of its last fetch, not the remote.

---

## 2026-09-22 — No HRRR: every run starts from observations at its own cycle time

**Context.** Prompts 94–102. The human redirected the model away from HRRR
entirely, in seven short steps that are worth keeping in order, because the
AI's first reading of them was wrong and was corrected by the fourth:

1. "lets not use the hrrr but real data like radar surface obs soundings etc..."
2. Asked what should feed the lateral boundaries: "the model will have to
   interpolate values from surface observation radar data and soundings near
   the area".
3. "take in all data sources it can that are reliable ... If a source isnt
   availbe that hour then the model will skip it and use the data it has."
4. Four cycles, 00/06/12/18Z, each 12–24 h long.
5. Each run finishes in under 1.5 hours.
6. "the model should be based off inital conditions so the 00z run uses 00z
   conditions and models from there" — **this corrected the AI**, which had
   read (2) as boundaries built from observations *during* the forecast window,
   i.e. a hindcast, and had planned around a run that could only start a day
   late.
7. Where soundings can't be found, use the previous run's forecast. After the
   forecast window, check the forecast against surface observations.

**Decision, as built.**

| | before | now |
|---|---|---|
| initial state | HRRR analysis | every reliable observation valid at the cycle time; missing sources skipped and logged |
| upper air with no soundings (06Z, 18Z, missing sites) | — | the previous run's forecast valid at that time; cold start from a standard atmosphere if there is none |
| lateral boundaries | hourly HRRR analyses | held to the initial analysis for the whole run (nothing observed later may enter) |
| terrain | HRRR surface height | a static DEM (ETOPO via NOAA ERDDAP), fetched once |
| verification | ASOS after the run | surface observations, once the forecast window has closed; the cycle-time analysis scored only on withheld stations |
| budget | none | 1.5 h wall clock per run at ≤ 50 % of cores, verification excluded |

`ingest_hrrr.py` is kept, not deleted, and reachable with `--source hrrr`, so an
HRRR-seeded run remains available as a baseline for the same day.

**What this costs, stated before building it.** Frozen edges are the price of
a true forecast from observations. Air crosses about 860 km in 12 h at 20 m/s,
and the domain is about 1300 km across, so by hour 12 much of the interior is
downstream of an edge that stopped changing at hour 0. Scores at long lead
times will partly measure the boundaries, not the dynamics. The analysis
extends beyond the model domain so observed upstream air at least starts at
the edges. This is written down now so that a poor 24 h score is not later
read as a dynamics failure.

**First contact with the live services** (desktop, 2026-09-22; P-06 is about
exactly this):

| service | finding |
|---|---|
| IEM RAOB | **the fetcher's request is rejected.** `raob_url()` sends `ts1`/`ts2`; the service now requires `sts`/`ets` (HTTP 422), accepts **one station per request** (4-character limit), and wants the `K` prefix (`KOKX`). |
| IEM RAOB availability, 2026-09-21 12Z | of 21 active IDs in the domain plus a ~3° ring, 10 returned data: APX, BUF, DTX, GSO, GYX, IAD, MHX, OKX, PIT, RNK. CAR and ILN had 00Z but not 12Z; Albany (`KALB`/`KALY`/`_ALY`), Wallops and both Canadian sites (`CWMJ`, `CWQI`) returned nothing at either time. **Inside the domain: 6 soundings, not the 8 the docs assume.** RNK (Blacksburg) is inside the domain and missing from `NORTHEAST_RAOB`. |
| IEM ASOS | `sts`/`ets` accepted; one hour for 20 US and Canadian networks = 665 rows, 56 kB, 0.4 s; `mslp` and `alti` available. The network list contains look-alikes that are other countries (`DE__ASOS`, `MA__ASOS`, `MD__ASOS`, `PA__ASOS`) — selected by exact name, not prefix. |
| NDBC | 191 met-reporting buoys and fixed stations within 1° of the domain; 5-day files ~67 kB each. |
| NOAA ERDDAP `etopo180` | plain CSV of terrain height; no library needed. |

The RAOB result is one more interface defect found only on first contact with
a live service, in code that every offline test passed (category A).

**A prediction about the run in progress.** `tools/daily.sh` passes the
forecast `--run-dir $DATA/tensors/analysis_<stamp>`, but `ingest_hrrr.py`
writes to `config.TENSOR_DIR` = `$DATA/tensors_3d/analysis_<stamp>`. Predicted,
before the log comes back: ingest succeeds, **forecast fails with "No
live_hrrr_f*.npz in .../tensors/analysis_20260922_00"**, verify is SKIPPED, and
the run ends `status=1`. If the forecast step instead finds its files, this
reading of the two paths is wrong.

**Predictions for the observation-built runs** (written 2026-09-22, before any
analysis code exists; each is checked on 2026-09-21 12Z and 2026-09-21 18Z
unless stated):

| # | prediction | fails if |
|---|---|---|
| P1 | leave-one-out sounding temperature error, averaged over the soundings available: ≤ 2.0 K RMS at 500 hPa, ≤ 3.0 K at 850 hPa | either is exceeded |
| P2 | 500 hPa heights integrated hydrostatically from the surface-pressure anchor and the analysed virtual temperature match the soundings' reported heights to ≤ 30 m RMS | > 30 m (the integration or the anchor is wrong) |
| P3 | surface temperature at withheld ASOS stations (every fifth) after elevation correction: ≤ 2.5 K RMS | > 2.5 K |
| P4 | a forecast from the obs-built state, prepared by the same filter and rebalance, survives 12/12 h at stride-4 spacing; 24 h is a guess of ≥ 18 h | fails before hour 12 |
| P5 | at 12 h lead, surface-temperature error within 200 km of the edges exceeds the interior's | the interior is as bad or worse, which would mean the edges are not what limits skill |
| P6 | removing one source (buoys) changes the surface analysis by < 0.5 K everywhere more than 300 km from every buoy | a larger change far away, which would mean the influence radius is too large |
| P7 | the 18Z run, with no soundings, takes its upper air from the 12Z run's 6 h forecast and says so in its metadata; its 500 hPa temperature differs from the 12Z analysis by < 3 K RMS | it silently cold-starts, or differs by more |

P6 and P7 are the ones that can fail for reasons other than skill: they test
that the plumbing does what it claims.

**Results, 2026-09-22** (desktop, 12 threads of 24; live observations for
2026-09-21 12Z and 18Z; each prediction checked as written above):

| # | result | verdict |
|---|---|---|
| P1 | leave-one-out over the 6 in-domain soundings: 500 hPa **1.71 K**; 850 hPa **3.81 K** (BUF +5.9, GYX +6.5, the other four within 2.2) | 500 holds; **850 fails** |
| P2 | 500 hPa height from the pressure anchor vs reported: **8.4 m** RMS (6 soundings) | holds |
| P3 | withheld ASOS temperature: **1.02 K** RMS, bias +0.25 (71 stations, 12Z); 1.74 K, bias −0.70 at 18Z | holds |
| P4 | 24 h requested; **diverged at 3.75 h** (max\|u\| 372 m/s), stopped inside the hour by the new guard | **fails** |
| P5 | no run reached 12 h | not assessed |
| P6 | with buoys removed, surface T changes by ≤ **0.21 K** more than 300 km from every buoy (max 5.25 K near them) — but only 9 % of the grid is that far from a buoy, so the test is weak | holds, weakly |
| P7 | the 12Z run died before +6 h, so the 18Z run had no usable forecast; it logged that and **fell to a standard atmosphere** | not assessed; see below |

**P4, located before anything was changed** (`src/analysis/probe_obs_blowup.py`,
5-minute snapshots): max\|u\| sat at 46.1 m/s for 3.5 h (at the lid, 2 cells
from the edge, pinned by the frozen boundary). The runaway was the
**meridional wind**, first growing by > 20 % between 2.50 and 2.59 h at
**level 17 of 20 (near the ground)**, 14 cells from the edge, over **812 m of
terrain** in northern Maine; only 2 points grew by > 5 m/s, with a ~5.5 Δx
pattern along the row. So it is interior, low-level, over terrain and near
grid scale, and **not** edge-driven. It is the same shape as P-50, whose
runaway was also v. The initialisation reported that the divergence did not
reach its target (4.7e-4 → 9.2e-5 s⁻¹). Opened as P-56; nothing has been
tuned.

**The P7 finding changed the code.** A standard atmosphere has no jet and no
gradients. When the previous forecast cannot supply +6 h and there are no
soundings at the cycle, the first guess is now the previous run's
**analysis** (6 h old, built from real soundings), and the label says so.
The human's rule ("if upper air soundings cant be found use the previous runs
forecast") says nothing about a previous run that died. This is the AI's
extension of it, flagged for confirmation.

**Defects found on the way, all in code that had passed its tests:**

1. The buddy check compared a sounding's 250 hPa temperature with a
   neighbour's 1000 hPa one (pressure-blind), and on 2026-09-21 12Z rejected
   every upper-air value it had buddies for.
2. After that fix, it still accepted a "consensus" of one neighbouring
   sounding's many significant levels (KRNK vs KGSO) — it now counts stations,
   not values.
3. The raob request rejected by IEM (P-54).
4. `daily.sh` run directory (P-55; predicted from the code and fixed before
   the log came back, so still unconfirmed on the server).

Items 1 and 2 would have hit verification against soundings as well.

**Cost.** Ingest 131–148 s, of which NDBC's ~200 sequential 5-day files take
~120 s. Forecast 1.36 min per forecast hour on 12 km (10 threads), so 24 h is
~33 min. A 24 h cycle is ~36 min on the desktop; the server is not measured.
Deferred verification of 3 h: 3077 pairs, **TMP RMSE 2.23 K, bias −1.47 K;
u 2.47, v 2.19 m/s**. Its QC buddy check is O(n²) and took ~5 min on 50 000
observations. That is outside the run budget, but it grows with window length.

**Not run here.** `src/verification/test_verification.py` crashes inside
`numpy.corrcoef` in this desktop sandbox; plain `numpy.corrcoef` crashes the
same way (a delay-loaded DLL), so that suite is not assessed on the desktop.
`tools/daily.sh` could not be syntax-checked: Git's bash cannot start in the
sandbox.

**Status.** Built and kept. P-56 open (the first observation-built forecast
dies at 3.75 h). The pipeline has run end to end on the desktop only.

**Addendum, same evening: a new server root.** Git on the server "got messed
up" (prompt 106). The rebuild recipe was run in a new, empty folder,
`/data5/pierce/AINWP`: a clean clone with every tracked file listed as
"deleted", which was the index describing files not yet written, not a loss.
The human made `AINWP` the project root (prompt 108) with its own data root,
`AINWP/data` (prompt 109). `/data5/pierce/NWP` and its nested archive are kept
untouched. A side effect worth recording: the old layout had a latent cron
defect. Cron does not read `~/.bashrc`, so `NWP_DATA_ROOT` would have been
unset under cron, and `daily.sh` would have fallen back to `$ROOT/data`, a
different directory from the one every hand run used. In the new layout the
default and the variable are the same path.

**The model uses one core, and always has** (prompt 111). The dynamics is
element-wise NumPy (stencils, `np.roll`-style differences), and NumPy runs
those on a single thread. The thread caps `resources.py` sets (`OMP_…`,
`MKL_…`, reported as "torch threads 10 max") bound BLAS and torch, neither of
which the core calls. So the "50 % of cores" rule has never been the binding
limit, and the log line suggests a parallelism that does not exist. The
desktop's 1.36 min per forecast hour was therefore a single-core number, and
the 1.5 h budget should be judged against single-core speed on the Xeon.
Nothing changed; measure first (the log's steps/s line), then decide.

**First server run, 2026-09-22 18Z** (`AINWP`, 13 min end to end, status 3):

| | |
|---|---|
| machine | 104 cores, 376 GB (so the "50 %" cap is 52 cores; one is used) |
| ingest | ~2 min; ETOPO terrain fetched from the server (0–1161 m) |
| first guess | no soundings at 18Z and no earlier run in `AINWP`: standard atmosphere. **Initial max\|u\| 6.2 m/s** — the jet is simply absent |
| speed | dt 17.1 s, 2.1 steps/s, so **1.7 min per forecast hour**; 24 h ≈ 40 min, inside the 1.5 h budget on one core |
| outcome | max\|u\| 6.3 → 7.4 → 10.1 → 18.3 → 14.6 → 37.7 m/s over hours 1–6, then **435 m/s at 6.31 h**, stopped mid-hour |

So P-56 has two cases: a realistic 12Z state (46 m/s jet, 10 soundings)
dying at 3.75 h, and a near-calm 18Z state dying at 6.31 h. What the two share
is the observation-built lower atmosphere, ETOPO terrain and single-frame
(frozen) edges. What they do not share is a strong upper flow. **A calm start
blowing up rules out "the jet is too strong" as the cause.** The HRRR-seeded
runs that reached 12/12 h differed in all three shared respects, so a
discriminating run is needed before any change.

**Where a forecast hour goes** (prompt 113; desktop, cProfile, 1 h from the
2026-09-21 12Z analysis, 87.6 s in `step`). There is no single hotspot. The
largest own-time entries are `np.roll` 11.0 s, `hyperdiffusion` 10.3 s,
`tendencies` 10.1 s, `richardson` 7.1 s, `pressure_gradient_force` 5.1 s,
`vertical_advection` 5.0 s, and a long tail of small stencils. A 5-point
Laplacian on one 20 × 97 × 110 field takes 2.24 ms with `np.roll` and 1.29 ms
with slicing (1.7×), a single-core gain available without threads. Whether
threads help depends on how torch scales on arrays this small (213k values)
when every one of thousands of operations per hour pays thread start-up.
`tools/bench_threads.py` measures exactly that on the server. **Prediction,
written before it runs:** torch will be ≤ 1.5× NumPy at 1 thread and will not
exceed 3× at any thread count for the Laplacian and the vertical cumsum. It
fails if some thread count gives > 3× on all three operations.

**P-56 test A, predictions written before the run** (prompt 114: "lets do
that"). The 2026-09-21 12Z analysis, rebuilt from the same archived
observations with the surface blend switched off, so the lowest kilometre
comes from the soundings and the first guess only. Heights keep the same
sea-level-pressure anchor. Same terrain, same frozen edges, same 24 h
request, desktop.

- A1: if the blend is the cause, the run passes hour 6 (the blended run died
  at 3.75 h), and max|v| at level 17 near row 82, col 89 stays under 20 m/s
  through hour 4.
- A2: if the blend is not the cause, it dies before hour 5, with v running
  away first at a low level over terrain.

Either way, only this one thing changes.

**Test A result.** Diverged at **4.10 h** (blended: 3.75 h). The 5-minute probe
puts the first > 20 % growth at 2.84–2.92 h, **v at level 17**, row 47,
col 53: central New York near the Catskills, **712 m** of terrain, 47 cells
from every edge. **A2 holds; A1 fails.** The surface blend is not the cause.
Removing it moved the first growth from northern Maine to central New York
and delayed it by about 20 minutes, but the signature is the same: the
meridional wind, third level above the ground, over 700–800 m of terrain,
nowhere near an edge.

Two side observations:

- The blended run's lowest-level temperatures differed from the unblended
  run's by up to 10.4 K, yet at the failure column the blend changed T by
  only 0.5 K. Superadiabatic layers at 1000–975 hPa occur in both states
  (676 columns blended, 1243 unblended, of 10 670). The blend removes about
  half of them rather than adding any.
- The probe (4.2 h requested) and the forecast (24 h requested) of the same
  state diverged at 3.94 and 4.10 h. `run_forecast` sets dt = duration /
  n_steps, so a different requested length is a slightly different dt, and
  in this regime that is not a neutral change (compare the chunking caution
  under P-40).

**Next discriminator** (not yet run): the same observation state over flat
ground. If the flat run passes hour 6, the failure is the state's
interaction with terrain, and HRRR terrain (test B) and slope limits become
the targets. If it still dies, terrain is ruled out and the frozen edges
(test C) are next.

**Flat-ground result** (same observations, same blend, terrain set to 0,
24 h requested, desktop, 20.9 min). max|u| held at 46.0 m/s through **12 h**;
the minimum theta fell from 279.8 to 276.1 K between hours 8 and 10. At
**13.0 h** the state went non-finite in a single hour (NumPy warned of an
invalid value in a power: a negative base, i.e. a pressure or column depth
going negative), with no gradual wind growth first. **Terrain is implicated:**
the same state dies at 3.75 h over terrain and lives 12 h without it. That
leaves the frozen edges as a secondary suspect at most, and only for the
separate hour-13 failure, which has a different signature (sudden non-finite,
not a v runaway). The 12/12 h flat-ground capability measured before
2026-09-22 was never tested past hour 12, so hour 13 is new ground for any
initial state.

**Mechanism to test next, written down before changing anything** (two
suspects have now been checked, per the method: A the blend, F the
terrain). An HRRR state has sensible values at pressure levels BELOW the
ground (NCEP extrapolates them smoothly). The observation analysis does not:
under 700–800 m of terrain, the 1000–925 hPa levels hold a first guess plus
Barnes increments from distant stations, never observed there. If
`pressure_to_sigma` draws on those levels when it builds the lowest sigma
levels over terrain, the state there is unconstrained, and level 17 over
700–800 m is exactly where it would show. Test: replace below-ground
pressure-level values with a smooth downward extrapolation from the lowest
level above ground (T along a standard lapse, winds held), and rerun the
blended 12Z state over ETOPO. Prediction: it passes hour 6. Fails if it dies
before hour 5 at level 17 again.

**Result: refuted.** The below-ground extrapolation changed 7507 values, by
up to 7.55 K, all under the terrain (nothing above ground moved). The run
diverged at **3.75 h**, to the step the same as without it (553 m/s).
Survival did not move at all, so the switch is kept in the code but **off by
default**.

That is two mechanisms tested and refuted (the surface blend, below-ground
values), and one discriminator positive (flat ground survives 12 h). Per the
method, no third guess is made from the desktop. What still differs between
the observation runs that die over terrain and the HRRR runs that survived
over terrain is the terrain field itself (ETOPO block-averaged vs HRRR's),
and the balance of the flow near the ground over it. Test B (HRRR terrain
under the observation state) separates the two, and needs the server.

**Tests B and C, predictions written before they run** (server, the
2026-09-22 18Z case, which died at 6.31 h over ETOPO; HRRR used as a
diagnostic only, in directories named `*_test*`, never archived or verified):

- **B** (`src/analysis/probe_terrain_b.py`): the same archived observations,
  rebuilt over HRRR terrain block-averaged onto the same cells. If the
  terrain field is the cause, it passes hour 6.31. If it still dies before
  hour 7, the terrain field is ruled out and the flow's balance over terrain
  is what is left.
- **C**: an HRRR initial state with ONE frame, so its edges are frozen like
  the observation runs'. The HRRR runs that reached 12 h had hourly edges. If
  frozen edges are harmless over terrain, C survives 12 h; if C dies early,
  the frozen edges are part of P-56 after all.

**Tests on the server, 2026-09-23 00Z–01Z** (prompt 115).

*Benchmark.* torch 2.8.0 against NumPy on 20 × 97 × 110 float64:

| threads | Laplacian | element-wise chain | vertical cumsum |
|---|---|---|---|
| 1 | 1.0× | 1.3× | 1.4× |
| 8 | **13.8×** | 6.6× | **6.4×** |
| 26 | 12.4× | **8.9×** | 3.2× |

**The prediction (≤ 3× for the Laplacian and the cumsum) failed**, and
decisively. At 8 threads all three operations are 6–14× faster, so a torch
port of the core is worth doing. 8 threads is the sweet spot for the stencils,
well under the 52-core ceiling.

*Test B* did not run: `probe_terrain_b.py` had not been pushed. It is now in
26eb43a.

*Test C* (HRRR start, one frame, so frozen edges; HRRR terrain at stride 4 =
10 km, 127 × 122): **diverged at 3.16 h.** max|u| went 53.6 → 65 → 68.6 →
177 m/s, then non-finite. The initial divergence was 1.85e-3 → 5.06e-4 s⁻¹
after filtering: about 20× the same cycle's observation state (18Z,
9.74e-5 → 2.53e-5), and about 5.5× the 2026-09-21 12Z one (→ 9.21e-5).

**Correction to this entry and to P-56, found by checking the record rather
than memory.** I wrote above that "the HRRR-seeded runs that reached 12/12 h"
differed from the observation runs. **There were no such runs.** The
2026-09-04 entry says so plainly: the sigma core became reachable from real
data that day, but "a real forecast has [not] been run". Every 12/12 result
in this project is on idealised states over idealised terrain. Test C is the
first forecast this core has ever made from a real analysis, and it dies
sooner than the observation runs. So P-56 is not an observation problem: **no
real state over real terrain has survived.** The flat-ground run (12 h) is
the only real-state run that has. That is lesson L3 (suspect the test before
the model) applied to a baseline: the comparison run P-56 was measured
against had never happened.

**Mechanism, written down before the next run.** The project measured, on
2026-09-0x, that survival over idealised terrain tracks slope (forced ascent
w = U × slope), not height: 2500 m with a maximum slope of 0.0086 survived
12/12. Real terrain at 12 km is far steeper. ETOPO block-averaged has a **maximum
slope of 0.0536** by the model's own measure (`sigma.terrain_slope`,
one-sided differences, the measure the 0.0086 was taken with; 0.0316 with
centred differences), **6.2× beyond anything shown to survive**. The two
failure points sit on slopes of 0.011 and 0.013 (centred).

Test S: the same 2026-09-21 12Z observations (blend on), over ETOPO smoothed
with the model's own `smooth_terrain` until the maximum slope is ≤ 0.0086,
analysis rebuilt on that terrain. **Prediction: it passes hour 6.** It fails
if it dies before hour 5, with v at a low level over terrain again. If it
passes, the operational answer is to smooth real orography the way
operational centres do, and the slope limit becomes a measured parameter.

**Test S result: the prediction holds.** `smooth_terrain` needed 11 passes to
bring the maximum slope from 0.0536 to 0.0083; the highest cell fell from
1161 to 881 m, and the domain mean did not change. The run held max|u| at
46.1 m/s and **survived to 13.70 h** (unsmoothed: 3.75 h), 21.4 min on the
desktop. **Slope is the mechanism of the early failure.**

It still died, and the way it died matches the flat-ground run: minimum theta
drifting down from hour 8 (275.2 → 271.2 K by hour 13), max|sigma_dot|
growing roughly fourfold from hour 9, then a runaway at 13.7 h. Two runs, one
flat and one smoothed, both fail between hours 13 and 14, well after the
edges have had time to act on the interior (air crosses ~860 km in 12 h).
That makes the frozen edges the first suspect for this SECOND failure mode.
Test C does not separate the two, because its terrain was unsmoothed.

**Kept:** ingest now limits the slope of the terrain it writes, to 0.0086
by default (`--max-slope`), with the number of passes logged. The HRRR path
is left unsmoothed, as a baseline.

**Next two runs, predictions written first** (prompt 116):

- **Live server cycle**, the latest 00/06/12/18Z with slope limiting on,
  24 h requested. It passes hour 6 (every unsmoothed real run died by 6.31 h)
  and stops between hours 12 and 16 with the second failure mode (theta
  minimum falling first). It fails if it dies before hour 6.
- **Locating the hour-13 failure** (desktop, test S state, 10-minute
  snapshots). If the frozen edges drive it, the first runaway is within ~15
  cells of an edge, or on the inflow (western) side. It fails if the first
  growth is deep in the interior, far from every edge.

**Located: the prediction holds.** Test S state, 5-minute snapshots,
diverged at 13.72 h (22 min, desktop). The first > 20 % growth came between
**11.59 and 11.67 h**: u from 50 to 64.5 m/s at level 15, row 84, **col 10,
10 cells from the western edge**. That is exactly the inner boundary of the
10-cell relaxation zone, on the inflow side of a westerly flow. In those
5 minutes 667 points grew by > 5 m/s (u) and 473 (v), at edge distances of
9–10 cells minimum and 19 median, in a ~3.7 Δx pattern. theta changed by up
to 6.7 K and pi by 6.2 hPa: a violent, near-grid-scale event, not a slow
drift. max|v| then grew steadily, 23 → 28 → 39 → 48 → 62 m/s over
11.75–13.42 h.

**Mechanism, stated before anything is changed.** For 11 h the interior
evolves (the cold air visible in the falling theta minimum is advected
eastward), while the relaxation zone keeps pulling the western edge back to
the hour-0 state. On an inflow boundary, that is a growing mismatch
concentrated where the relaxation weight goes from full to zero, which is
where it broke. A frozen edge is only harmless while the interior still
resembles hour 0. This is the cost P-53 predicted, arriving as an
instability rather than a slow error.

Candidate responses, none tried yet, all consistent with "nothing observed
later may enter": a gentler and wider relaxation for a frozen driver, or
relaxation only where flow enters, with the outflow edge left free. Each
needs its own prediction before its run.

**Server live cycle (2026-09-23 06Z, slope-limited terrain, 0–881 m).**
It ran unattended through `daily.sh` and diverged at **16.31 h** (max|u|
438 m/s) after 27.9 min, at 2.0 steps/s. The predictions, scored:

| Prediction | Result |
|---|---|
| Passes hour 6 | **holds** (it passed hour 15) |
| Stops at 12–16 h | **marginal miss**: 16.31 h, just past the window |
| The theta minimum falls first | **refuted**: it held at 278.6 K throughout |

The failure was abrupt, not a drift. max|u| went 15.8 m/s at 15 h, 36.5 at
16 h, then 438 at 16.31 h, and max|sigma_dot| rose fourfold in that hour.
So a falling theta minimum is not a precursor in general. In the 12Z case it
was cold advection, which this weak-flow night did not have.

**A separate finding in the same log: the 06Z state has almost no wind.**
max|u| was 8.5 m/s before initialisation and 6.4 after, anywhere in the
column. A real September atmosphere over the Northeast has 20–40 m/s near
200–250 hPa. There are no soundings at 06Z, and on the fresh AINWP root no
00Z forecast existed, so the first guess must have come from lower in the
fallback chain. Its winds aloft are then whatever the surface observations
and the fallback supply, which is close to nothing. `availability.json`
records which one was used. Once cycles run back to back, each 06/18Z run
starts from the previous run's +6 h. This run could not, and its forecast
is skilful only near the ground at best.

Next discriminator (prediction first): `tools/locate_growth.py` on the
saved snapshots. If this is the same frozen-edge failure, the largest
change between 15 and 16 h lies 8–12 cells from an edge, at the inner
boundary of the relaxation zone. It fails if it lies more than 20 cells in.

**Located (server, prompt 118): the prediction holds in every interval.**
The first guess was `standard_atmosphere`, which has no wind at all, so the
near-windless state is explained. In all 15 hourly intervals, the largest
change in u and in v lies **10–13 cells from the edge** (usually exactly
10 or 11). In 28 of the 30 u/v lines it is at rows 10–16 and columns 10–14:
the south-west corner of the relaxation zone's inner boundary, over
540–880 m of the West Virginia Appalachians. The two exceptions are early
and small, still on column 10 but further north over low ground: v at row
75 (146 m, 2.0 m/s, hours 2–3) and u at row 77 (163 m, 2.9 m/s, hours 5–6).
The corner cell itself, row 10 col 10, has the largest change in u in
hours 1–4.

| Hours | Points changing > 5 m/s | Largest change u / v (m/s) |
|---|---|---|
| 1–5 | 0 | 2–4 / 2–5 |
| 7–10 | 3–6 | 8–9 / 7–8 |
| 12–13 | 33 | 12 / 10 |
| 14–15 | 79 / 71 | 14 / 17 |
| 15–16 | 158 / 148 | 36 / 68 |

The minimum edge distance of the fast-changing points is 10 or 11 in
every interval where there are any (never less than 10), and the median is
10–12. The disturbance never moves in; it
widens along the zone boundary until it runs away. This is not the 12Z
picture of an inflow mismatch: this flow is under 10 m/s and almost nothing
evolves to mismatch. Both failures do sit near a western corner, though.
In the 12Z case, row 84 is 12 cells from the northern edge.

**Two explanations remain, and they make different predictions:**

- **H1, the zone boundary itself.** The failure is made where the
  relaxation weight falls to zero. With the cosine profile, alpha per step
  is 0.0245 one cell inside the boundary and 0 at it: a pinned state next to
  a free one. The location should then move with the width.
- **H2, the place.** Steep terrain near the western corners (with the
  sigma-coordinate pressure-gradient error in an unbalanced, standard
  atmosphere state) makes it. The coincidence with column 10 is then
  chance, twice.

**Test W (server, before any result): the same 06Z state, 18 h,
15-minute snapshots, relaxation width 6 and width 15, run side by side.**

| | H1 predicts | H2 predicts |
|---|---|---|
| Width 6 | first growth 6–7 cells from the edge; dies **before** 16.3 h (steeper weight: 0.067 per step one cell in) | growth stays at rows 10–16, cols 10–14 (edge distance ≈ 10) |
| Width 15 | first growth 15–16 cells in; survives **past** 16.3 h (0.011 per step one cell in) | the corner is inside the zone and pinned; growth appears elsewhere over terrain, or not at all |

H1 is refuted if width 6 keeps its growth at edge distance ≥ 9. H2 is
refuted if both runs move their growth with the boundary. If both hold in
part (moves with the boundary, but only in the west), the answer is the
combination: the zone boundary where it crosses steep terrain.

**Test W result (server, prompt 119). H1 holds and H2 is refuted as the
sole cause: the growth moves with the zone boundary.**

| | Width 6 | Width 10 (06Z run) | Width 15 |
|---|---|---|---|
| alpha per step, one cell inside | 0.067 | 0.0245 | 0.011 |
| Outcome | diverged 16.20 h | diverged 16.31 h | reached 18.00 h, running away (max\|v\| 31.0 → 91.5 m/s in the last 15 min) |
| Edge distance of the largest change | 5–9 in 119 of 126 lines; 10–11 in 5; 32 and 43 in 2 (at 0.5–0.75 h, ≤ 1.1 m/s) | 10–13 in 30 of 30 | 15–19 in 105 of 142; 25–38 in 37 (all by 7 h, ≤ 1.0 m/s) |
| Minimum edge distance of points changing > 5 m/s | 6 in 24 of 29 intervals, 7–9 in 5 | 10–11 | 15 in 25 of 39, 16–17 in 14 |
| Where it ran away | southern edge, rows 8–9, cols 46–49, terrain 1–2 m (near the Delmarva coast) | south-west corner, 540–880 m | south-west corner, rows 15–21, cols 15–18, 730–850 m |

Scored against the predictions:

- Width 6 dies before 16.3 h: holds, but by 0.11 h, which is too small to
  mean anything.
- Width 6 grows 6–7 cells in: holds (mostly 6–8).
- Width 15 survives past 16.3 h: holds, by about 2 h. It was already
  running away when the run ended.
- Width 15 grows about 15 cells in: holds (mostly 15–17).
- H2 (it stays at the Appalachian corner): refuted. The width-6 run ran
  away over flat coastal ground at 1–2 m.

The south-west corner is still where the width-10 and width-15 runs fail,
so the place has some influence, but it is neither necessary nor the cause.

**What this establishes.** The failure is made at the inner boundary of
the relaxation zone, wherever that boundary is. A gentler weight ramp
delays it (by about 2 h from 0.0245 to 0.011 per step), and a steeper one
does not bring it much earlier. So the ramp's steepness is not the main
control. What all three runs share is the thing itself: a band pinned to
the hour-0 state next to an interior that is free to evolve.

**A defect found reading these logs: P-57.** At 18.00 h the width-15 log
printed "max|u| 17.9 m/s" while max|v| was 91.5 m/s. The in-loop guard and
the progress log both look only at u. The first P-56 runaway (3.75 h) started in v, and so did the last
width-15 hour, so a v runaway could pass the 150 m/s ceiling unreported
until u followed: the P-52 class again, "failure not detected". Fixed in `forecast.py`: the guard
now uses max(|u|, |v|) and finiteness of both, and the hourly line prints
max|v|. `test_forecast.py` passes 11/11. Earlier divergence times were detected on u, so they
may be late by however long u took to follow v. The locations and the
order of events come from snapshots and are unaffected.

**Next: test O, predictions first.** Same 06Z state, 24 h, 15-minute
snapshots, run side by side:

- **O1, no relaxation** (`--relax-width 0`). The edges then follow the
  model's own replicate condition, and no observed-later data is involved
  either way. If the pinned/free interface is what makes the failure,
  there is no growth band at a fixed edge distance and the run survives
  past 18 h. It fails if it dies by 18 h with its growth 5 or more cells
  in. If it dies with its growth at edge distance 0–2, that is a different
  failure (the free edge), not a refutation. It is recorded as a new
  problem, and it would say the zone is needed but in a different form.
- **O2, weak pinning** (width 10, `--relax-alpha 0.1`: 0.1 per step at the
  edge, 0.00245 one cell in). If the pull's strength matters, rather than
  only its presence, it survives past 18 h. If it fails, its growth is
  still about 10 cells in.

**Test O result (server, prompt 124).**

| | O1: no relaxation | 06Z run (width 10, alpha 1) | O2: width 10, alpha 0.1 |
|---|---|---|---|
| Outcome | diverged **6.67 h** (163 m/s) | diverged 16.31 h | diverged **21.67 h** (154 m/s) |
| Edge distance of the largest change | 0–2 in 43 of 50 lines, 3–4 in 7 | 10–13 in 30 of 30 | 8–13 in 149 of 170; the other 21 at 5–7, 14–18 or 26–43 (all by 7.75 h, ≤ 1.6 m/s) |
| Minimum edge distance of points changing > 5 m/s | 0 in 22 of 22 intervals | 10–11 | 9–10 in 45 of 54, 8 in 2, 11–12 in 7 |
| Where it ran away | southern edge, rows 0–4, cols 53–109, over the sea (0 m) | south-west corner, 540–880 m | south-west corner, rows 9–14, cols 9–14, 660–880 m |

Scored against the predictions:

- O1 survives past 18 h if the interface is the cause: **refuted.** It
  died much earlier, at the physical edge (edge distance 0). Growth started
  at 4 h along the southern boundary over the ocean. By the rule written
  beforehand, this is a different failure (the replicate edge), not a
  vindication of the zone. **Some relaxation is needed**: a free edge fails
  in under 7 h.
- O2 survives past 18 h if the pull's strength matters: **holds**, 21.67 h
  against 16.31 h.
- O2 fails about 10 cells in: **holds**, mostly 10–11. It runs away at the
  same south-west corner over the Appalachians.

**What the three tests say together.**

| Change from the 06Z default | Hours gained |
|---|---|
| width 10 → 6 | −0.1 |
| width 10 → 15 (gentler ramp, 0.011 one cell in) | ≈ +2 |
| alpha 1 → 0.1 (0.00245 one cell in) | **+5.4** |
| no zone at all | −9.6 |

The failure is made at the zone's inner edge, and it weakens as the pull
toward the frozen hour-0 state weakens. Without a zone, the edge itself
fails first. The south-west corner (the Appalachians under a
standard-atmosphere first guess) is where it breaks whenever the zone
boundary lies near it. The simplest reading: the frozen edge state is not
in the model's own balance over that terrain. The stronger the pull back
to it, the faster the mismatch at the zone's inner edge feeds the growth.

This is not yet a fix, and weakening the pull is tuning, not a cure: every
case still fails, only later. But the cycle needs 24 h, and 21.67 h is
close.

**Next: test P (predictions first).** Same 06Z state, 24 h, run side by
side:

- **P1**, width 15 and alpha 0.1, the two helpful changes together. If the
  effects add, it survives 24 h. It fails if it dies before 21.67 h (worse
  than alpha 0.1 alone).
- **P2**, width 10 and alpha 0.03, the dose-response check. If strength is
  the control, it outlives O2 (> 21.67 h). If it dies earlier than O2,
  weaker is not simply better. The edge would then be too free, as in O1,
  and its growth would move toward edge distance 0–4.

A 12Z case with real soundings is needed before any default changes. The
06Z state has a standard-atmosphere first guess, and a setting tuned on
one windless night is not a result.

**Test P result (server, prompt 130). Both runs completed 24 h, the first
real-data sigma forecasts in this project to do so.**

| | P1: width 15, alpha 0.1 | P2: width 10, alpha 0.03 |
|---|---|---|
| Outcome | **completed 24 h** (45.7 min, NumPy) | completed 24 h (45.6 min) |
| max\|u\|, max\|v\| over the run | 10.0, 11.1 m/s | 22.5, 18.5 m/s |
| 15-min intervals with a point changing > 5 m/s | **0 of 94** | 14 of 94, from 21 h (at most 9 points) |
| Largest 15-min change | 2.3 m/s | 15.2 m/s, rising at the end |
| Where the change is largest | spread over 10–43 cells in (112 of 188 lines at 14–20; no narrow band), all under 2.3 m/s | 9–12 cells in (149 of 188 lines), south-west corner, 730–870 m |

Scored against the predictions:

- P1 survives 24 h: **holds**, and with no sign of the zone-boundary growth.
- P2 outlives O2 (21.67 h): **holds**, it reached 24 h.
- P2's growth moves toward edge distance 0–4 if the pull is too weak:
  **refuted**. It stayed at the zone boundary (9–12 cells in) and was
  growing when the run ended, so it would probably have failed within a
  few hours.

Weaker pull keeps delaying the same zone-boundary failure. The wider, weak
zone is the only setting found that removes it for 24 h. **One case is
not a result**, though: this is the 06Z standard-atmosphere night, with
winds under 12 m/s and nothing much for an edge to fight.

**Next: test Q, a 12Z case with real soundings (predictions first).**
Build 2026-09-25 12Z from observations, then run three 24 h forecasts
side by side:

- **Q0n**, default zone, NumPy.
- **Q0t**, default zone, torch × 8. This is the real-case backend check:
  it should match Q0n to round-off early (≤ 1e-9 in the first hours).
  Differences may grow only where the run is already failing. It should
  be ≥ 2.5× faster in wall time.
- **Q1**, width 15 and alpha 0.1, torch × 8.

Predictions:

- Q0 fails before 24 h, somewhere in 10–17 h, near the zone's inner
  boundary. That is the 12Z 2026-09-21 case again, which died at 13.7 h.
- Q1 completes 24 h with no 15-minute interval where more than 20 points
  change by > 5 m/s. It fails if it dies before 24 h or shows the band at
  the zone boundary.

If Q1 holds, width 15 / alpha 0.1 becomes the default: two cases, one calm
and one with a jet. If it fails, the zone's form, not its strength, is
next.



---

## 2026-09-25 — Forecast maps and a Pivotal-style viewer

**Context.** Prompt 120: the forecasts need maps "like how a site like
pivotal weather has their maps". The user chose an HTML viewer and all four
product groups: surface, upper air, the analysis with its observations, and
forecast minus observed. The server has matplotlib but no cartopy, and
nothing may be installed (constraint 1).

**What was built.**
- `src/maps/geography.py`: a spherical Lambert conformal projection (39/45 N,
  centred 42 N 74 W) in NumPy. Natural Earth 1:50m coast, lakes, borders and
  state lines are fetched once as GeoJSON, clipped with the json module and
  NumPy, and cached like the ETOPO terrain.
- `src/maps/derive.py`: every product field from the model's own variables
  and constants. That covers heights by the core's own hydrostatic
  integration, ln-p interpolation to 1000/850/700/500/250 hPa, MSLP,
  thickness, de-staggered winds and absolute vorticity. Hour 0 is drawn from
  the analysis.
- `src/maps/render.py`: 7 forecast products, 4 analysis-with-reports products
  and 2 error products, each on a fixed colour scale.
- `src/maps/viewer.py`: one self-contained `index.html` per run.
- `src/make_maps.py`: renders it all. `daily.sh` runs it after each forecast,
  and `verify` adds the error maps. A map failure never changes a cycle's
  status.

**Checks.** `src/maps/test_maps.py` 11/11:
- Heights are within 1.7 m of the standard atmosphere at 1000–250 hPa, on
  flat ground and over 1000 m.
- MSLP reduces a standard atmosphere over 0–1500 m to 1013.25 hPa.
- Interpolation is exact at the model levels.
- Vorticity is 0 for uniform flow and 2Ω for solid-body rotation.
- The projection is conformal, with scale 0.9986–1.0033 over the domain.
- Barbs are rotated onto the projected meridians to within 3×10⁻⁵ degrees.
- Every product renders, and the viewer embeds every product and hour with
  no external URL.

A synthetic 24 h run (a moving low over analytic terrain, with the captured
2026-09-21 12Z station payloads as reports) rendered 179 images in 45 s on
the desktop, well inside the cycle's 8-minute reserve.

**What is not verified yet.** The state and coast lines have not been
drawn anywhere. The desktop sandbox had no network, so the Natural Earth
fetch happens on the server's first run, and the alignment of lines with
the grid is checked by eye there. The JavaScript viewer was not run in a
browser here (no browser or node in the sandbox). Its manifest is tested;
its behaviour is not.

**Honesty in the labels.** The model is dry and has no 2 m or 10 m
diagnosis. So "near-surface" maps say "lowest model level (≈ 236 m above
ground)" in their titles, and there are no precipitation, dewpoint, radar
or CAPE products. The analysis map shows 2 m temperature because that is
what the analysis fitted. Withheld stations are drawn open and in purple,
because they score the analysis and did not build it.

**Revision after the first server render (prompt 123).** The user listed:
too many station plots; the map should be filled edge to edge; only F000
shows; valid times off on some products; hover for values; click for a
SHARPpy-style sounding; follow Pivotal's parameters. Changes:

- **Only F000.** A real defect, P-58. The snapshots land up to one 17 s
  model step after each hour, and the map script demanded exact hours. It is fixed and tested
  against the model's own step rule.
- **Frame.** The map is now the largest rectangle inside the projected
  domain. The grid is drawn beyond it and clipped, so no blank corners
  remain.
- **Station plots.** They are thinned to 45 km (55 km for pressure), and
  the colour-bar label says how many are shown.
- **Times.** Titles follow Pivotal's pattern: "Init: 06z Sep 23 2026
  Forecast hour: 12 / Valid: 18z Wed Sep 23 2026 (2 PM EDT Wed)". Eastern
  time comes from the US DST rule, with no time-zone database needed. What
  the user saw as "off" has not been identified. The one known time error
  (P-58) would have hidden every hour but 0, so the user is asked which
  product it was.
- **Pivotal parameters.** Added 925 mb temperature/height/wind, 700 mb
  vertical velocity (kinematic, -µb/s), 300 mb jet and a 12-hr temperature
  change. The 12-hr change became a 1-hr change at prompt 126, drawn from
  F001 on (F001 is the analysis-to-first-hour step), on a ±8 °F scale in
  0.5 °F steps. Units are mb and kt, °F at the surface and °C aloft.
- **Hover and sounding.** Per-hour JavaScript data files (int16, a scale
  and offset per field, loaded with script tags so file:// works): 32–33
  map fields at every grid point, and model columns at every second point.
  That is about 35 MB of data per 24 h run and about 70 MB with the images.
  The page inverts the projection in JavaScript. The frame geometry is
  tested by drawing markers and finding them within 1 px of where the page
  computes them. The sounding has no dewpoint or CAPE except at hour 0,
  because the model is dry, and it says so.

`src/maps/test_maps.py` 15/15. The synthetic 24 h run drew 315 images in
58 s. The page's JavaScript is still untested in a browser.

**First real maps and first server verification (prompt 125: `maps.zip`
from the 06Z 2026-09-23 run).**

- **Maps.** All 17 hours (F000–F016, the run died at 16.31 h) and all 17
  products rendered, with 16 hours of error maps after `verify`.
- **Lines.** The state, coast and lake lines line up with the geography
  (Finger Lakes, Long Island, Cape Cod, Chesapeake), and the stations sit
  on land and coast where they should.
- **Viewer data.** The data files decode sensibly. At Albany at F006:
  terrain 242 m, MSLP 1030.7 mb, lowest-level 42.6 °F, 850 mb 6.1 °C. The
  sounding column has p_s 999.8 mb and 200 hPa −55.9 °C.
- **What the maps expose.** The 850 mb and 700 mb maps are nearly
  uniform, with calm winds. That is the standard-atmosphere first guess of
  a 06Z cycle with no soundings and no previous run, made visible at a
  glance.

The first server verification (≈350 surface stations an hour, read from
the error maps; `verification_20260923_06Z.csv`):

| Lead (valid) | T bias / RMSE (°C) | Wind bias / RMSE (kt) |
|---|---|---|
| F001 (07Z) | +1.68 / 2.63 | +0.86 / 3.84 |
| F005 (11Z) | +2.16 / 3.31 | +0.86 / 4.10 |
| F007 (13Z) | −0.67 / 2.14 | −0.78 / 3.82 |
| F010 (16Z) | −4.91 / 5.84 | −2.92 / 4.82 |
| F014 (20Z) | −6.88 / 7.72 | −2.38 / 5.02 |
| F016 (22Z) | −5.82 / 6.53 | −1.18 / 4.41 |

The error is the missing diurnal cycle (P-59): the core has no surface
heating or radiation. The model is too warm and too windy by night and
too cold and too calm by day, with the sign change at sunrise. No
persistence reference was scored, so how much of the error the model
adds, or removes, relative to holding the analysis fixed is not assessed.


---

## 2026-09-25 — The PyTorch backend

**Context.** Prompt 127 asked whether more CPU would speed things up;
prompt 128: "do the pytorch port". The core is element-wise NumPy on one
core, and the server benchmark (prompt 115) put torch at 6–14× on
model-sized arrays at 8 threads.

**Design: one physics source, two backends.** A second copy of the
physics in torch would drift from the first. Instead,
`src/dynamics/backend.py` gives the hot code a NumPy-named namespace:
NumPy itself for arrays, or a torch mapping (float64, CPU) for tensors.
Nine functions across `grid`, `sigma`, `turbulence`, `surface`,
`convection`, `primitive_sigma` and the forecast relaxation now ask
`xp_of(array)` which to use. NumPy constants (levels, Coriolis, terrain,
sponge, relaxation weights) pass through `xp.asarray`. That is free for
NumPy and a bounded cache for torch, because `ndarray * tensor` silently
turns the tensor back into NumPy. `grid.shift` now slices instead of
`np.roll` plus an edge fill: the same values, and expressible in both
backends. Not ported: stochastic physics and the radiative top (both off
in production). `to_backend("torch")` refuses them rather than running
them wrongly.

**Checks.**

| Check | Result |
|---|---|
| NumPy path after the refactor vs before, 300 steps, realistic case | **bit-identical** (all four fields) |
| All 11 dynamics suites + forecast, maps and analysis suites (NumPy) | 74 + 47 tests pass |
| torch vs NumPy, 300 steps (44×40) | 1.6e-12 relative (u), round-off |
| torch vs NumPy, 60 steps, full 110×97 grid | 1.7e-13 relative |
| `forecast.py` 1 h on a synthetic analysis, torch vs NumPy | 3.7e-13 relative (v); both "completed" |
| `test_backend.py` (new) | 3/3 |

**Speed on the desktop (full grid, 60 steps).**

| Backend | Time | Speed-up |
|---|---|---|
| NumPy | 20.1 s | 1.0× |
| torch × 1 | 5.9 s | 3.4× |
| torch × 4 | 2.4 s | **8.4×** |
| torch × 8 | 2.4 s | 8.4× |
| torch × 12 | 3.0 s | 6.7× |

A forecast hour in `forecast.py` took 10.5 s against 90 s. The server
benchmark's peak near 8 threads holds here: past 8 threads it gets slower.
Initialisation (filter and balance, about a minute) still runs in NumPy.

**Not yet done.** The server has torch 2.8 (the desktop has 2.13) and
Python 3.9. `tools/check_backend.py` runs the same comparison there. The
backend becomes the cycle default (`NWP_BACKEND=torch` in `daily.sh`)
only after that shows round-off agreement.

**Server check (prompt 129).** `test_backend.py` passes 3/3 on the Xeon
(torch 2.8.0+cpu, NumPy 1.24.4, Python 3.9). `check_backend.py`, full
grid, 60 steps:

| Backend | s/step | Speed-up | Max relative difference |
|---|---|---|---|
| NumPy | 0.237 | 1.0× | — |
| torch × 1 | 0.308 | 0.8× | 2.2e-12 |
| torch × 4 | 0.114 | 2.1× | 2.2e-12 |
| torch × 8 | 0.085 | **2.8×** | 2.1e-12 |
| torch × 12 | 0.087 | 2.7× | 2.1e-12 |
| torch × 16 | 0.089 | 2.7× | 2.1e-12 |

Agreement holds everywhere. The speed-up is a third of the desktop's, and
the estimate given to the user ("about 5 minutes for 24 h") was wrong.
That estimate was carried over from the desktop without asking why the
desktop's ratio was so large.

The per-component profile on the desktop shows the reason. There NumPy is
the slow side: a step takes 358 ms against the server's 237. Torch at 8
threads takes 41 ms on the desktop and 85 ms on the server, a slower,
shared core. The ratio therefore depends on both machines' NumPy as much as
on torch.

| Desktop, ms | NumPy | torch × 1 | torch × 8 |
|---|---|---|---|
| step | 358 | 97 | 41 |
| tendencies (×3 per step) | 112 | 30 | 15 |
| vertical mixing | 18.8 | 6.4 | 2.5 |
| surface drag | 10.1 | 3.2 | 0.9 |
| hyperdiffusion (×3 per tendency) | 8.4 | 1.5 | 0.5 |
| continuity | 6.2 | 1.4 | 0.6 |
| geopotential | 5.1 | 2.2 | 1.5 |
| convective adjustment | 1.9 | 1.1 | 0.4 |

At torch × 1 no single term dominates. What remains is the cost of many
small operations, which is what fusing them (`torch.compile`) would
attack. That needs a C++ compiler at run time: to be checked on the
server, and not testable on this desktop.

Expected for a real 24 h cycle on the server: the production forecast ran
at about 0.48 s/step in NumPy, so torch × 8 should give about 0.17 s/step,
or about 15 min instead of about 40 plus a minute of set-up. Measured on
the real case next.

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

## 2026-09-10 — The 8 September session lands in the repository, and two of its numbers do not survive

**Context.** The work of 2026-09-08 existed only as notes outside the
repository: the K_MAX ladder that returned identical numbers at every setting,
the lesson drawn from it, and a staleness checker described but never written.
The last commit on `main` was 2026-09-08 17:37 and contained none of it. This
entry records bringing it in, which turned out not to be a transcription job.

**Hypothesis (stated before anything was run).**

1. The identical ladder results come from Python's default-argument binding:
   `vertical_mixing(..., k_max=K_MAX)` freezes `K_MAX` at import, so
   `turbulence.K_MAX = x` at runtime changes nothing the model reads. If that
   is the mechanism, setting the global and calling with the default must
   return the import-time ceiling exactly, while passing the value explicitly
   must return the new one.
2. The 8 September note that the codebase contains "5 dated measurements
   against 19 undated" should reproduce, at least in direction and roughly in
   magnitude, once a checker is written.

**Method.**

- A neutral column (constant θ, so N² ≈ 0 and the Richardson function is at
  full strength) with 12 m/s of shear per level, run through
  `turbulence.vertical_mixing` twice: once with `turbulence.K_MAX` assigned at
  runtime and the default call, once with `k_max=` passed explicitly.
- `tools/stale.py` written and run against the whole repository, with file
  dates from git and again from filesystem mtimes.
- The three mixing parameters moved to instance state on `PrimitiveSigma`;
  `lid_test.build_on` given a `**model_kw` passthrough; `kmax_ladder.py`
  changed to vary the ceiling through the constructor.
- A guard-rail test added to `test_primitive_sigma.py`, and the suite re-run.

**Result.**

Hypothesis 1 — confirmed exactly:

| how the ceiling was set | max K |
|---|---|
| `turbulence.K_MAX = 400`, default call | 100.0 |
| `k_max=400` passed explicitly | 400.0 |

So every rung of P-40's ladder ran the same experiment, and "6/12, 6/12, 6/12"
was not an elimination. Re-run through the constructor on 2026-09-08 the same
ladder gives 6/12, **8/12**, **8/12**. P-40 is back in OPEN; the code defect
is P-51 and is FIXED.

Hypothesis 2 — **did not survive.** `tools/stale.py` as committed reports, on
this repository:

| | count |
|---|---|
| measurements carrying a date | 432 |
| measurements with no date within twelve lines | 476 |
| dated measurements in files changed >14 days later, git dates | 0 |
| the same, from filesystem mtimes on a fresh clone | 12 |

Not 5 and 19. The direction holds — there are more undated numbers than dated
ones — but the magnitudes are off by two orders of magnitude, because the
earlier figures came from a version of the checker that no longer exists and
whose scope is unknown. Nothing can be reconstructed from a number reported
without the script that produced it.

**Interpretation.** The staleness checker's first finding was about a claim in
the note that proposed it, which is the tightest possible demonstration of the
lesson it was built for (L7): a measurement about the project's measurements
went stale in two days.

What it does *not* mean: that the 8 September session was wrong. Its central
finding — the knob was not connected, and a recorded negative result was
therefore worthless — reproduces exactly and is the more important of the two.
The count was a supporting detail stated with more precision than it had.

There is a second correction. L5 ("a diffusion loses a race against a growing
mode") cited P-40 as evidence that raising K_MAX does nothing. That evidence
is withdrawn. The lesson's own argument is untouched — a rearrangement beats a
faster diffusion because the wave steepens in less than the diffusive
timescale — but a good pattern attracts confirming evidence and does not
check it, and this one collected a number that was never real.

**What this changes about the tall-terrain problem.** The eddy-diffusivity
ceiling is worth two forecast hours at 4000 m and saturates between 300 and
1000 m²/s, which is the shape of a constraint that binds and then stops
binding. Before that is called an improvement it has to pass the L2 test: the
prediction written down first is that max|u| must NOT fall, because a scheme
that buys stability by flattening the flow looks identical on a survival
count. That measurement has not been made.

**Status.** Kept. P-51 FIXED with a guard rail; P-40 reopened with the re-run
in it; `docs/LEARNING_LOG.md` and `tools/stale.py` now in the repository.

**For the collaboration study.** One defect, category D (language semantics
defeating a correct experimental design), detected by *a result being too
clean* — two numbers agreeing to one decimal place. That detection route
appears nowhere else in this project and nothing in the toolchain looks for
it; it was noticed, not caught. Both stated hypotheses were tested before
anything was rewritten, and one of them failed, which is the reason to write
them down. The human intervention was one sentence long ("can you sync now")
and administrative in tag, but it is what forced the reconciliation that found
both corrections — a transfer that had been recorded as done, and was not.

---

## 2026-09-12 — The ceiling was binding on one interface in a thousand, and the prediction that failed was the useful one

**Context.** P-40 was reopened on 2026-09-08 after the ladder that eliminated
it turned out not to have varied anything (P-51). The re-run said the
eddy-diffusivity ceiling was worth two forecast hours at 4000 m. Two forecast
hours is the largest single movement this problem has ever had, which is
exactly the kind of result that should be distrusted until it has been asked
to fail.

**Hypothesis, committed before any run finished** (in the docstring of
`src/dynamics/kmax_binding.py`, commit 907c5bb, made while the ladder was
still running):

| | prediction |
|---|---|
| P1 | survival rises between 100 and 300 |
| **P2** | **max\|u\| does NOT fall as the ceiling rises** — the one written to be able to fail |
| P3 | the overturning fraction at a given hour falls as the ceiling rises |
| P4 | mid-level stratification away from the mountain is unchanged |
| P5 | the realized diffusivity never reaches 1000, so the ceiling stops binding and raising it further changes nothing |

P2 is the discriminator. A ceiling that buys stability by flattening the jet
is suppression, and the two hours would be worthless — which is precisely how
the first sponge failed (P-16) and how P-49 still fails.

**Method.** 4000 m terrain, 8-level sponge, clean and filtered, 12-hour
ceiling, one run per rung at K_MAX = 100, 110, 125, 150, 200, 250, 300, 1000.
Recorded hour by hour: max|u|, the free-troposphere jet, min Ri, N² in the
mid-troposphere, realized max K, the fraction of interfaces the ceiling is
actually clipping, the overturning fraction, and — added partway through, for
reasons below — the mean timestep and the wall clock. Compared at a COMMON
hour rather than at each run's own last hour, since a run that lives longer is
otherwise being compared at a later and harder time.

**Result — survival.**

| K_MAX | 100 | 110 | 125 | 150 | 200 | 250 | 300 | 1000 |
|---|---|---|---|---|---|---|---|---|
| survived | 6/12 | 6/12 | 7/12 | 8/12 | 8/12 | 8/12 | 8/12 | 8/12 |

A monotone ramp between 100 and 150 and flat above it. The earlier reading
that it "saturates between 300 and 1000" was an artifact of three rungs.

**Result — the discriminators, at hour 6.**

| K_MAX | max\|u\| | jet | overturning | N² mid | clip% |
|---|---|---|---|---|---|
| 100 | 54.8 | 48.1 | 0.369% | 2.214e-04 | 0.11 |
| 150 | 54.7 | 48.4 | 0.375% | 2.214e-04 | 0.02 |
| 200 | 54.7 | 48.5 | 0.375% | 2.214e-04 | 0.01 |
| 300 | 54.7 | 48.5 | 0.379% | 2.214e-04 | 0.00 |
| 1000 | 54.7 | 48.5 | 0.378% | 2.214e-04 | 0.00 |

P2 held: the wind is 54.7 m/s at every setting and the jet is marginally
*stronger* with more mixing, not weaker. P4 held to four significant figures.
P5 held — at K_MAX 1000 the realized diffusivity peaks at 605 and the ceiling
never clips at all, so above ~600 the parameter is inert.

**P3 failed, and it is the most useful line in the table.** More available
mixing does not reduce the overturning; it goes very slightly the other way,
and the convective adjustment fires at the same rate at every setting. The
mechanism I proposed — that a higher ceiling mixes away the static instability
faster — is wrong. Something else is buying the hours.

**Interpretation.** The ceiling was truncating the diffusivity the scheme
itself asked for, on about **one interface in a thousand**, in the breaking
region, and that truncation cost two forecast hours. Above ~150 the formula
rarely asks for more, and above ~600 it never does.

What it does *not* mean: that we know why. Since the extra diffusivity changes
neither the overturning fraction nor the stratification nor the wind, the
remaining candidate is that it acts on MOMENTUM in the breaking layer —
removing shear that would otherwise concentrate into the runaway that ends the
run — rather than on buoyancy. That is a hypothesis, it is untested, and the
way to test it is the momentum budget in that layer, not another ladder.

**The production case is untouched.** 2500 m with a 5-level sponge, which is
the terrain this project is for: 12/12 at 100, 150 and 300, max|u| 44.0 / 43.9
/ 43.9, jet 23.5 at all three, N² identical. So the default can move without
paying anything where it matters. Raised to 200 — clear of the ramp, identical
on every discriminator to 150 and 300, and still at the bottom of the 10²–10³
m²/s observed in breaking mountain waves.

**A second defect, found by accident** (P-52). K_MAX = 110 spent 911 seconds
of wall clock in forecast hour 7 without finishing. Inside that hour the wind
reached **7176 m/s** with the timestep collapsed to a mean of 0.87 s. Every
sweep in this project tests for failure at the hour boundary, and `run()`
itself only stops on a non-finite surface pressure, so a run that is blowing
up but still finite keeps integrating — and gets slower as it does, because
the adaptive dt ratchets down and never back up. The ladder spent more time on
its two failing rungs than on the six that survived.

**And a caution about the instrument.** The first attempt at a stall guard
integrated each hour in ten-minute chunks so it could check between them. That
changed the answer: `run()` truncates its final step to land exactly on the
requested duration, so chunking changes the step sequence, and at K_MAX = 110
the chunked and whole-hour runs diverge (min Ri 0.012 against 0.022 at hour 6)
and then fail differently. The regime is marginal enough that subdividing the
integration differently is not a neutral act. The working guard is a callback
that only reads the clock, and it was checked to reproduce the un-chunked run
exactly before being trusted.

**Status.** Kept. P-40 closed with the default at 200; P-52 opened; the
mechanism behind the two hours left explicitly open rather than narrated.

**For the collaboration study.** Four predictions were written down and
committed before any result existed, and the one that failed is the one that
produced knowledge: P3's failure is what turned "the ceiling dissipates the
overturning" from a conclusion into an open question. The two that held (P2,
P4) only licensed the change; they did not teach anything. Note also that
the predictions were committed to git *while the runs were still going*,
which is the cheapest possible way to stop a hypothesis from drifting toward
the data — a habit worth keeping, and one that no tool enforces.

Defects this session: one category E/H (P-52, failure detected only at hour
boundaries), found not by a test but by noticing that a run was taking too
long. That is the second defect in this project detected by something being
anomalous rather than wrong — the first was P-51, found because two numbers
agreed too well.

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

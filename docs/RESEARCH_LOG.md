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

**Viewer confirmed in a browser (prompt 131).** The user opened the 06Z
2026-09-23 maps. The hover readout and the click-for-sounding panel both
work, and the valid times are correct on every product. The "valid times
are off" complaint (prompt 124) was made on the first map version. The
revision after it changed the hour matching (P-58) and the title times,
and the current build no longer shows the problem. Nothing on the viewer
is still unconfirmed.

**Torch on a real case (prompt 132, 06Z 2026-09-23, default settings).**
The torch × 8 run diverged at **16.31 h**, the same as NumPy. The largest
relative difference was 1.4e-12 at 1 h, 1.0e-11 at 7 h, 2.9e-10 at 10 h
and 1.4e-8 at 16 h. It grew only once the state became unstable, as
predicted. Wall time was 19.2 min against 27.9 min, **1.45×**. That is well
short of the 2.8× from `check_backend.py`, and not yet explained.

**Test Q result (prompt 133, 12Z 2026-09-25, first guess `sounding_mean`,
max|v| 33 m/s at the start).**

| Run | Setting | Backend | Outcome | Wall |
|---|---|---|---|---|
| Q0n | default zone | numpy | diverged 7.45 h | 9.1 min |
| Q0t | default zone | torch × 8 | diverged 7.45 h | 3.8 min (2.39×) |
| Q1 | width 15, alpha 0.1 | torch × 8 | diverged 12.65 h | 6.6 min |

The three ran side by side. Scored against the predictions:

- **Q0 fails in 10–17 h near the zone boundary: refuted.** It failed
  sooner, at 7.45 h, and not at the zone boundary.
- **Q1 completes 24 h with no bad interval: refuted.** It diverged at
  12.65 h. Its first interval with points changing by > 5 m/s was
  4.25–4.50 h (214 points).
- **Q0t matches Q0n to ≤ 1e-9 in the first hours: holds** through 3.75 h
  (8.0e-10), with 2.0e-9 at 4 h. Both diverged at 7.45 h.
- **Torch ≥ 2.5× faster: narrowly missed**, at 2.39×.

**What the failure is (P-60).** Both runs are quiet until 4.00–4.25 h,
then blow up at the same place:

- levels L03–L05 (about 270–330 hPa, jet level);
- rows 76–78, columns 86–90 (central Maine);
- 17–21 cells from the nearest edge.

The zone settings made no difference to the onset, only to how long the
wreck took to reach 150 m/s. From 5 h on, 5 100–11 500 points change by
> 5 m/s each 15 minutes, so Q1's extra 5 h are not forecast hours.

The numpy–torch difference is the useful measurement. It is round-off
amplified by whatever grows fastest, and it grows with an **e-folding time
of 24 min from the first snapshot**. So the mode is present in the
initial state, not created by the edge later. That is too fast for
inertial instability (at most about f, an e-folding near 3 h). It is fast
enough for shear or static instability, or a numerical mode. The onset
band also straddles the base of the wind sponge (L00–L04), the vertical
counterpart of P-56's lateral zone edge.

**Decisions.**

- Width 15 with alpha 0.1 does **not** become the default. It fixed one
  mechanism on one case and was no help against the other.
- **Torch is cleared for production.** Two real cases diverge at the same
  step as NumPy, with differences at round-off that grow only through the
  instability. Set `NWP_BACKEND=torch NWP_THREADS=8`.

**Next: test R (predictions first).**

1. `tools/mode_structure.py` on the Q0n/Q0t pair, and as a check of the
   tool on the 06Z pair.
   - The Q mode sits in the onset box (rows 70–84, columns 80–96,
     L03–L05) from hour 1.
   - The 06Z mode sits at the south-west zone boundary (edge distance
     9–13) by hour 10, where `locate_growth` put that failure. If it does
     not, the tool is wrong and Part 1 is not evidence.
2. Stability in the Q initial state at the onset box, which decides
   between the candidates:
   - Ri < 0.25 or N2 < 0 at L03–L05 there means the analysed flow is
     unstable: candidate (a).
   - Ri > 1 and eta/f > 0 at every level there points to (b).
3. Sponge depth, Q case, default zone, torch, 8 h:
   - sponge 8 levels (**R8**) and sponge 3 levels (**R3**), against 5.
   - If (b), the onset level moves with the sponge base: to about
     L06–L08 in R8 and L01–L03 in R3, and the onset time changes.
   - If (a), the onset stays at L03–L05 in central Maine at 4.0–4.5 h in
     both.
   - Either result eliminates one candidate. If both runs move, or
     neither does, both candidates are open.

**Test R result (prompt 134).**

*Tool check, 06Z pair: holds.* `mode_structure.py` puts the 06Z mode at
the south-west zone boundary (edge distance 10–11, rows 11–15,
columns 10–12, 38.2–38.7 N 80.2–80.5 W) at 6, 10, 14 and 16 h, where
`locate_growth` found that failure. It starts at L17–L18 and has risen
to L11–L14 by 14–16 h. At t+1 h the lowest interface there is statically
unstable (N2 −2.8e-5, Ri −1082 at L18/L19), and 394 points in the domain
have Ri < 0.25 at L18. That is noted, not pursued.

*The Q mode.* Its position by hour:

| Hour | Rows | Columns | Levels | Edge distance |
|---|---|---|---|---|
| 1 | 69 | 84–85 | L02–L04 | 24–25 |
| 3 | 72 | 89–90 | L03–L04 | 19–20 |
| 4 | 77–78 | 88–89 | L04 | 18–19 |

- The prediction (rows 70–84, columns 80–96, L03–L05 from hour 1) holds
  approximately. At hour 1 the mode is one row south and one level higher
  than predicted; from hour 3 it is inside the box.
- It is top-heavy. In the box, each level from the lid (L00) to L04
  carries 85–100 % of the peak wind difference at 1–3 h. The theta
  difference peaks lower, at L10 (525 hPa).
- It e-folds in 46, 24, 22 and 17 min over hours 0.5–1, 1–2, 2–3 and 3–4.

*Stability at the start (t+0.25 h, ± 6 cells round r77 c89).* eta/f is
at least 0.66 at every level, Ri at least 1.59 at every interface, and
N2 at least 6.3e-5 s⁻². No point in the domain has Ri < 0.25 at L00–L18.
By these grid-scale measures the initial flow in the box is stable.
**Candidate (a) is refuted as stated.**

*Sponge depth (Q case, default zone, torch, 8 h).*

| Run | Sponge | First 15 min with > 5 m/s changes | Where the largest changes are then | 8 h |
|---|---|---|---|---|
| Q0 | 5 levels | 4.00–4.25 h (32 points; 629 in the next) | L03–L05, r76–78 c86–90 | diverged 7.45 h |
| R8 | 8 levels | 4.25–4.50 h (351 points) | L00–L01, r74–75 c87–88 | completed, 8 216 points > 5 m/s at 5.25–5.50 h |
| R3 | 3 levels | 3.50–3.75 h (1 point, edge 10); 4.25–4.50 h (48–50) | L03–L04, r84 c83–87, then L02 r77 c93 | completed, 5 359 points > 5 m/s at 6.00–6.25 h |

Under (b) the onset level should have moved with the sponge base, to
L06–L08 in R8 and L01–L03 in R3, and the onset time should have changed.
Instead the onset time is 4.0–4.5 h in all three runs, in the same part
of Maine. R8's largest changes moved **up** to the lid, not down.
**Candidate (b) is refuted.** A deeper sponge slows the growth after
onset (R8 and R3 both reached 8 h, Q0 did not), but it does not move the
mode or delay it. The design paragraph's last sentence ("If both runs
move, or neither does, both candidates are open") contradicted the two
explicit predictions above it. The runs are scored on the explicit
predictions.

*Speed.* R8 and R3 ran side by side at 0.105 and 0.135 s/step (1828
steps in 3.2 and 4.1 min), which is close to `check_backend.py`. The 06Z
torch run's 0.335 s/step now looks like load on the machine, not the
code. A 24 h torch forecast is about 10–12 min.

**Two candidates refuted, so stop and re-examine the assumptions before
testing a third** (the rule this project runs on).

1. **Part 2 measured the wrong time.** Growth sped up from 46 to 17 min
   while the difference was still 1e-11–1e-8, well inside the linear
   range. A linear mode on a steady flow grows at a constant rate, so the
   flow under it changed during the first 4 h and became more unstable.
   Stability has to be measured at 2–4 h, not 0.25 h. Measurement next:
   Part 2 at 1–4 h round the mode's track.
   - Prediction: if a resolved physical instability is responsible, the
     minimum Ri at L00–L04 falls below 0.25 (or N2 below 0) in the box by
     3–4 h, before the 4.25 h onset.
   - If Ri stays above 1 and eta/f above 0 through 4 h, the mode is
     numerical.
2. **The mode reaches the lid.** It spans L00–L04, and the box's
   strongest wind is at the lid itself (32.6 m/s at L00, 206 hPa). The
   200 hPa lid cuts through the jet. On the 06Z night the wind at the
   lid was under 1 m/s, and that case had no such mode. This is
   candidate (c), the rigid lid under a jet. It is **not tested yet**.
   It gets written predictions after the measurement in item 1.

**Measurement R2 (prompt 135): stability round the mode at 1–4 h.** The
box is ± 8 cells round r73 c87, Q0n. Minimum Ri at each interface:

| Interface | 1 h | 2 h | 3 h | 4 h |
|---|---|---|---|---|
| L02/L03 | 1.12 | 0.92 | 0.90 | 0.90 |
| L03/L04 | 0.91 | 0.44 | 0.30 | **0.28** |
| L04/L05 | 1.62 | 0.95 | 0.58 | 0.39 |
| L05/L06 | 1.71 | 1.54 | 1.24 | 1.47 |

- N2 in the box stays at or above 5.8e-5 s⁻² at every level and every
  hour. eta/f stays above 0 at L00–L06 (its box minimum is 0.32–0.87).
- Over the whole domain, points with Ri < 0.25 at L03/L04 number 0, 0,
  19 and **579** at 1, 2, 3 and 4 h. The runaway begins at 4.0–4.5 h.

**Scored: neither branch.** The prediction assigned Ri < 0.25 to a
physical instability and Ri > 1 to a numerical mode. The result, a fall
from 0.91 to 0.28, lands in the gap between them, which the prediction
left unassigned. It was not sharp enough.

What it does show: the jet's lower flank (L03–L05) sharpens through the
same 4 h in which the growth speeds up. **The model has no vertical
dissipation there.** `eddy_diffusivity` is exactly zero for Ri ≥ 0.25,
and a desktop check confirms it: with the default settings a 1 h run is
bit-identical with mixing on and off. Horizontal hyperdiffusion is the
only thing opposing the sharpening, until Ri crosses 0.25 and mixing
switches on at up to l²|S|. That happens across hundreds of points at
3–4 h, just before the runaway.

Candidate (c), the lid, cannot be tested cleanly: the analysis stops at
200 hPa, the lid itself, so raising the lid would need extrapolated data.

**Candidate (e): no vertical dissipation in a sharpening shear layer.
Test S (predictions first).** Q case, default zone and sponge, torch,
8 h. Two new options, `--ri-crit` and `--no-mixing`; the default path
is bit-identical.

- **S1, `--ri-crit 1.0`.** Mixing acts wherever Ri < 1, weighted
  (1 − Ri)², so it is weak at 0.5–1.
- **S0, `--no-mixing`.** This is a control.

Predictions. Onset means the first 15-minute interval with more than 20
points changing by more than 5 m/s: 4.00–4.25 h in Q0, 4.25–4.50 h in R8
and R3.

- If (e): S1's box minimum Ri at L03/L04 stays at or above 0.5 through
  4 h, and S1's onset is after 5.0 h or absent within 8 h. (e) is
  refuted if S1's onset falls in 4.0–4.75 h.
- S0 is identical to Q0 until mixing first acts (about 3 h), and its
  onset stays in 4.0–4.75 h. If S0's onset moves out of that range, the
  switch-on of mixing below 0.25 is part of the trigger.

If S1 holds, an Ri_c of 1 is still not adopted from one case. It needs
the 06Z case and a second jet case, and a stable-regime mixing function
chosen for physical reasons rather than a threshold tuned to this
failure.

**Test S result (prompt 136).**

| Run | Setting | Onset (first 15 min with > 20 points changing > 5 m/s) | Where | 8 h |
|---|---|---|---|---|
| Q0 | Ri_c 0.25 | 4.00–4.25 h (32) | L03–L05, r76–78 c86–90 | diverged 7.45 h |
| S1 | Ri_c 1.0 | 4.00–4.25 h (43) | L05, r77 c88 | completed (5.4 min) |
| S0 | no mixing | 4.00–4.25 h (48) | L04, r77 c88 | diverged 7.37 h |

Minimum Ri in the box at L03/L04 for S1 was 0.91, 0.53, 0.41 and 0.37 at
1, 2, 3 and 4 h (Q0: 0.91, 0.44, 0.30, 0.28). Extending mixing to Ri < 1
barely slows the sharpening.

- **(e), S1 onset after 5 h with box Ri ≥ 0.5: refuted.** The onset is
  unchanged, and Ri is 0.41 by 3 h.
- **S0 identical to Q0 until mixing acts, onset in 4.0–4.75 h: holds.**
  Every `locate_growth` line matches Q0 through 2.50–2.75 h. The first
  difference is at 2.75–3.00 h (3.5 against 3.3 m/s) at the south-east
  zone-boundary point r10 c99, so mixing first acts about 2.75 h.

Vertical mixing does not control the P-60 onset: off, as it is, or
extended to Ri < 1, the onset is 4.00–4.25 h at r77 c88.

*Side finding (P-56, not P-60).* Without mixing, one point on the
south-east zone boundary runs away first (r10–12 c97–99, edge 10, L03–L05,
over the sea). Its v changes by 4.3, 10.6 and then 36.5 m/s per 15 min
over 3.25–4.00 h. With the default mixing the same point peaks at
3.3 m/s. Mixing had been holding a zone-boundary point in check.

**Three P-60 candidates refuted: (a), (b), (e). (c) cannot be tested
cleanly. Change of strategy.** Guessing a mechanism and testing it has
now failed three times. The next test instead splits the possibilities,
each arm separating one broad class from the rest.

**Test T (predictions first).** Q case, default settings otherwise,
torch, 8 h.

- **T1, `--dt-factor 0.5`** (timestep about 15.8 → 7.9 s; the Q runs took 1828 steps for 8 h).
  - A time-stepping instability would change its growth rate. T1's
    onset would be at least 1 h later (5.0–5.25 h or after), or absent
    within 8 h.
  - An onset within 4.0–4.75 h rules out time discretization.
- **T2, `--hyper-factor 4`** (2Δx damped in 45 min instead of 3 h; 4Δx
  in about 12 h).
- **Roughness** of the Q0n–Q0t u difference at 3–4 h, now printed by
  `mode_structure.py`. It is the share of variance a 3×3 mean removes.
  Calibration: 2–4Δx 0.9–1.0, 6Δx 0.56–0.80, 8Δx 0.35–0.58, 10Δx
  0.24–0.42, 20Δx 0.06–0.12. At 1 h the difference is still raw
  round-off (0.96 on the synthetic pair), so only 3–4 h counts.
  - A grid-scale mode (≤ about 6Δx): roughness ≥ 0.8 and T2's onset at
    least 1 h later.
  - A resolved mode (≥ about 10Δx): roughness ≤ 0.4 and T2's onset
    within 4.0–4.75 h.
  - **Declared in advance:** roughness 0.4–0.8 is inconclusive, and T2
    alone decides.

If neither T1 nor T2 moves the onset and the mode is smooth, it is a
resolved structure: something in the analysed state, or resolved
dynamics. The next step is then to look at the analysis over Maine, not
another switch.

**Test T result (prompt 138).**

| Run | Setting | Onset (first 15 min with > 20 points changing > 5 m/s) | Where | 8 h |
|---|---|---|---|---|
| Q0 | default | 4.00–4.25 h (32) | L03–L05, r76–78 c86–90 | diverged 7.45 h |
| T1 | timestep × 0.5 | 4.00–4.25 h (272) | L03–L04, r77 c89–90 | diverged 7.65 h (8.2 min) |
| T2 | hyperdiffusion × 4 | **5.75–6.00 h** (48) | L04, r76 c87 | completed (3.0 min) |

The roughness of the Q0n–Q0t u difference is 0.98, 0.99, 0.99 and 0.99
at 1, 2, 3 and 4 h.

- **Time discretization: ruled out.** Halving the timestep leaves the
  onset at 4.00–4.25 h, and the growth after onset is, if anything,
  faster.
- **Grid-scale mode: both parts of the prediction hold.** Roughness is
  0.99 at 3–4 h, above the 0.8 threshold and in the 2–4Δx band, and T2's
  onset is 1.75 h later. At 1 h a roughness of 0.98 is expected from
  round-off alone. By 4 h the difference has grown 5 000-fold into one
  box, so the 0.99 there describes the mode.
- **The delay fits a 2–3Δx scale (a rough check).** Hyperdiffusion × 4
  adds three times the default grid-scale damping rate (1/3 h) at 2Δx,
  0.56 of that at 3Δx and 0.25 at 4Δx. With a mean e-folding of 24 min
  and an onset at 4.1 h without it, the estimated onsets are 6.8 h
  (2Δx), 5.3 h (3Δx) and 4.6 h (4Δx). The observed 5.75–6.00 h lies
  between 2Δx and 3Δx.

So P-60 is a **grid-scale (2–3Δx) spatial-discretization mode** at jet
level. It is present from the start, does not depend on the timestep, is
damped but not removed by hyperdiffusion, and grows faster as the shear
under the jet sharpens. Stronger hyperdiffusion alone is not a fix: T2
still blows up after 6 h, and damping 2Δx in 45 min is tuning against a
growth rate, not a correction to its source.

**Next, measurement U: which term feeds the mode.** `tools/mode_budget.py`
takes the Q0n state and the Q0n state plus the scaled mode (1 mm/s, so
linear). It evaluates every tendency term of the core at both, and
prints each term's share of the mode's kinetic-energy growth in the box.
Advection is split into the base flow carrying the mode and the mode
carrying the base flow. It was checked on the synthetic pair, where the
terms rebuild the core's own tendency difference to 3.6e-12 and the
advection splits add up exactly.

- Check: the implied energy growth rate at 3 h must be within a factor 2
  of 2/(22 min) = 1.5e-3 s⁻¹. That is 0.75–3.0e-3 s⁻¹, from the measured
  2–3 h e-folding. If it is outside that range, the tendencies do not
  describe the mode. The relaxation zone and convective adjustment act
  outside them, but neither is active in the box: it is 18+ cells in and
  N2 > 0.
- Expectation, not a test: the largest positive term is the mode's
  vertical motion acting on the base flow's vertical shear,
  d(sigma_dot) dU/dsigma. That is where the sharpening is (Ri falls at
  L03–L05) and the mode is deep but grid-scale in the horizontal. The
  budget will show if the horizontal term d.grad U or the pressure
  gradient term dominates instead.

**References in AMS format (prompt 137).** The user asked that anything
taken from other people be referenced in American Meteorological Society
format. `docs/REFERENCES.md` now lists 38 works:
- methods and formulas: Davies relaxation, the Wicker–Skamarock RK3,
  Simmons–Burridge, Barnes and Koch et al., Louis, Manabe et al., the
  Magnus constants, Bolton, and others;
- data: IEM, NDBC, MRMS, ETOPO1 through ERDDAP, HRRR, Natural Earth;
- software: NumPy, PyTorch, Matplotlib, xarray, cfgrib, Herbie, Pillow;
- designs the maps imitate: Pivotal Weather, SHARPpy.

Each DOI was resolved against CrossRef or DataCite, and author initials,
titles, volumes and pages were taken from the registry record. Two
records needed correcting. CrossRef misspells Buizza et al.'s second
author ("Milleer"). The PyTorch paper has no DOI, so its entry was
checked against the NeurIPS proceedings and dblp.

In-text author–year citations were added where each work is used, and
the three "&" citations were changed to "and". The new
`tools/check_refs.py` enforces the rule. It flags a citation with no
entry, an "&", or an entry nothing cites, and it joins the pre-commit
checks.

The check also found Sadourny (1975), cited in the shallow-water core,
which the first inventory had missed.

**Measurement U result (prompt 139): the P-60 kinetic-energy budget.** Q0n
against Q0t; box ± 8 round the mode; P/E in 1/s.

| Term | 2 h | 3 h | 4 h |
|---|---|---|---|
| vertical advection: d(sigma_dot) dU/dsigma | **+2.13e-3** | **+1.46e-3** | **+2.52e-3** |
| vertical advection: sigma_dot d/dsigma d | −2.2e-6 | −2.4e-6 | +1.9e-6 |
| horizontal advection (both parts) | +4.3e-5 | −2.5e-5 | −4.4e-5 |
| pressure gradient | −1.02e-3 | +2.10e-3 | −5.25e-3 |
| sponge | −3.4e-4 | −4.1e-4 | −2.8e-4 |
| hyperdiffusion | −8.0e-5 | −8.5e-5 | −8.5e-5 |
| sum | 7.2e-4 | 3.04e-3 | −3.14e-3 |

- **Check 1 holds.** The terms rebuild the core's tendency difference to
  within 1.5e-12 at every hour.
- **Check 2, as specified, fails.** Energy growth at 3 h was to be 0.75–
  3.0e-3 s⁻¹; it is 3.04e-3. At 2 h it is 7.2e-4, just below the band,
  and at 4 h it is negative. The cause is visible in the table: the
  pressure-gradient term changes sign from hour to hour (−1.0, +2.1,
  −5.3e-3). An oscillating mode trades kinetic and potential energy
  through that term, so one snapshot of its kinetic-energy budget swings
  with the phase. The check was specified for the wrong quantity. It
  should be the time mean.
- **The expectation holds.** At every hour the one steady positive term
  is the mode's vertical motion acting on the base flow's vertical
  shear, d(sigma_dot) dU/dsigma, at 1.5–2.5e-3 s⁻¹. Without the
  oscillating pressure term, that source minus the sponge and
  hyperdiffusion sinks (about 4e-4) leaves about 1.6e-3 s⁻¹. The
  measured energy growth is 1.5–1.9e-3 s⁻¹, from the amplitudes
  1.3e-10, 2.0e-9 and 5.6e-8 at 2, 3 and 4 h.
- **The 06Z corner mode is different** (12 h, south-west corner, L15–L19).
  Its sources are horizontal advection with the mode carrying the base
  flow (d.grad U, 1.45e-4) and vertical mixing (+8.2e-5). That fits a
  low-level, mixing-active boundary-zone problem (P-56). It is not P-60.

**What this points to.** Horizontally the mode is 2–3Δx (test T). Its
energy comes from the vertical shear, through its own vertical motion,
while the grid-scale Ri stays at or above 0.28. Continuous theory does
not allow shear instability there (Miles 1961; Howard 1961). So the
restoring force the discrete mode feels is weaker than the one Ri
measures.

A known way that happens: the core uses the Lorenz grid (u, v and theta
on the same full levels, sigma_dot on half levels). The Lorenz grid
carries a vertical computational mode, a level-to-level zigzag in theta
that the hydrostatic pressure gradient barely feels (Arakawa and Konor
1996). **Candidate (f): the P-60 mode rides on the Lorenz computational
mode.**

**Measurement U2 (predictions first).** `mode_budget.py --hours 2,4`, the
mean over the nine 15-minute snapshots, plus the mode's vertical
structure at 4 h.

1. The time-mean energy growth is within a factor 2 of 1.68e-3 s⁻¹
   (0.84–3.4e-3). That is the check re-specified for the time mean.
2. The pressure-gradient term's time mean is small: its magnitude under
   30 % of the mean d(sigma_dot) dU/dsigma term. d(sigma_dot) dU/dsigma
   stays the largest positive mean term.
3. If (f) holds, theta in the mode zigzags. Its adjacent-level
   correlation is below −0.5 across the mode's levels (L00–L06) at 4 h.
   If that correlation is above 0 there, (f) is refuted: the mode is
   smooth in the vertical and grid-scale only in the horizontal. The
   next question is then how the C-grid's divergence, and so
   sigma_dot, responds at 2–3Δx. Between −0.5 and 0 is **declared
   inconclusive in advance**.

**Measurement U2 result (prompt 140).** Q0n against Q0t, eight snapshots
from 2.00 to 3.75 h.

1. **Time-mean energy growth 1.08e-3 s⁻¹, inside 0.84–3.4e-3: holds.**
2. **The mean pressure-gradient term is small: holds.** It is −2.8e-4,
   15 % of the mean d(sigma_dot) dU/dsigma term (1.84e-3, 71 % of the
   total). The pressure term swings from −4.7e-3 to +3.0e-3 between
   snapshots.
3. **(f), a zigzag in theta: refuted.** Theta's adjacent-level
   correlation is 1.00 at every pair from L00/L01 to L18/L19. For u it
   is 0.83–1.00 over L00–L09, and for sigma_dot 0.94–1.00.

The mode is smooth and deep in the vertical, and the reason is in
`continuity`. sigma_dot on each half level is built from the divergence
integrated from the lid down. A divergent disturbance confined to
L00–L09 therefore puts the same horizontal pattern into sigma_dot all the
way to the ground. Theta then follows sigma_dot acting on the stable
background theta.

So the P-60 mode is:
- **2–3Δx in the horizontal** (test T);
- **divergent** (its energy arrives through its own sigma_dot);
- **oscillating** (the pressure-gradient exchange);
- **deep** (lid to about 500 hPa);
- **independent of the timestep**;
- **fed by the jet's vertical shear**.

That is the description of a grid-scale gravity-wave mode drawing on the
shear.

**Candidate (g): a divergent grid-scale mode, which divergence damping
should remove. Test V (predictions first).** Divergence damping acts only
on the divergent part of the wind (Skamarock and Klemp 1992). It is new
in `subgrid.divergence_damping`, off by default. The new tests check that
a non-divergent flow gets zero tendency (5.5e-16), that a 2Δx wave decays
at exactly 4ν/dx², that torch equals numpy, and that the default core is
bit-identical. `--div-damp C` sets ν = C dx dy / dt.

At the Q case's dt of 15.8 s and dx of 12 km:

| C | ν (m²/s) | 2Δx e-fold | 3Δx | 10Δx | 20Δx |
|---|---|---|---|---|---|
| 0.01 (V1) | 9.1e4 | 6.6 min | 8.8 min | 69 min | 4.5 h |
| 0.003 (V2) | 2.7e4 | 22 min | 29 min | 3.8 h | 15 h |

The mode's amplitude grows at about 8.4e-4 s⁻¹ (20 min e-fold).

- **If (g):** V1 has no onset within 8 h (no 15-minute interval with more
  than 20 points changing by more than 5 m/s). V2's onset is at least
  2 h later than Q0's, or absent.
- **(g) refuted:** V1's onset falls in 4.0–4.75 h, as in Q0.
- **Cost check:** V1's max|u| and max|v| stay within 1 m/s of Q0's
  through 3.75 h, so the resolved flow is barely touched before the
  onset. On the desktop synthetic case, C = 0.01 changed max|u| by
  0.9 m/s in 1 h.
- **Measurement:** `mode_budget.py` now prints rms(divergence) /
  rms(vorticity) of the mode per level. For a divergent mode it should
  be above 1 at L00–L06.

If V1 holds, divergence damping is a treatment, not yet a default. It
needs the 06Z case, a second jet case, and verification scores no worse
than without it.

**Test V result (prompt 141): (g) holds.** Q case, default zone, torch,
8 h.

| Run | ν actually used | 2Δx / 10Δx e-fold | Largest 15-min change, 0–6.25 h | Intervals with > 5 m/s | 8 h |
|---|---|---|---|---|---|
| Q0 | 0 | — | 33.3 m/s | from 4.00 h | diverged 7.45 h |
| V1 | 5.88e4 m²/s | 10.2 / 107 min | 4.6 m/s | **none** | completed (2.6 min) |
| V2 | 1.77e4 m²/s | 34 / 355 min | 4.3 m/s | **none** | completed (2.5 min) |

- **V1 no onset within 8 h: holds.** Not one point changes by more than
  5 m/s in any interval of the printed 0.25–6.25 h, and the run completed
  8 h.
- **V2's onset at least 2 h later, or absent: holds** (absent, same
  evidence).
- **Cost check: holds.** max|u| and max|v| stay 27.6 and 33.2 m/s in
  both runs, the same as Q0 through 3.75 h (and to 6.25 h).
- **The mode is divergent: holds.** rms(divergence) / rms(vorticity) of
  the Q0n–Q0t mode at 3 h is 51.7, 15.8, 5.45, 2.40, 4.48 and 1.95 at
  L00–L05, and above 1 at 17 of the 19 level pairs.
- **What is left** (both runs, 5.25–6.00 h): v changes of 3.8–4.6 m/s at
  L03, rows 85–86, 10–11 cells from the northern edge. That is the
  default zone's inner boundary, the P-56 kind of growth, not P-60.

**A bug in `--div-damp`, found in these logs.** The option computed ν
from `model.max_dt()` before the initial state was loaded, when u = v = 0
and theta was unset. dt came out 24.5 s instead of 15.8 s, so ν was 0.64
of the intended C dx dy / dt: 5.88e4 rather than 9.1e4 m²/s for V1. The
design table above (9.1e4 and 2.7e4) is therefore not what ran. The
table here quotes the ν the logs print, which is what was used.

Fixed: ν is now computed from the loaded state's CFL step. On the desktop
synthetic case C = 0.01 now gives 9.45e4 m²/s, against 9.42e4 expected
from C dx dy / dt. Weaker damping than designed was enough: even V2's
1.77e4 removes the onset.

**Test W (predictions first): 24 h, both cases, with and without the
wide weak zone.** C = 0.0064, which with the fix reproduces V1's ν on the
Q case (about 5.8e4 m²/s; about 5.4e4 at the 06Z case's 17.1 s step).

| Run | Case | Zone |
|---|---|---|
| W1 | Q (12Z 2026-09-25, jet) | default |
| W2 | Q | width 15, alpha 0.1 |
| W3 | 06Z 2026-09-23 (calm) | default |
| W4 | 06Z | width 15, alpha 0.1 |

Predictions:

- **W1:** no interior P-60 onset through 24 h. Either it completes 24 h,
  or it fails from zone-boundary growth, with the median edge distance of
  the > 5 m/s points at 12 or less.
- **W2:** completes 24 h with no interval in which more than 20 points
  change by more than 5 m/s.
- **W3:** divergence damping does not touch the 06Z corner mode, which
  the budget showed is fed by d.grad U and mixing near the ground. W3
  diverges in 15–18 h at the south-west corner, 9–13 cells in.
- **W4:** completes 24 h, as P1 did.
- **Skill, via the new `tools/score_by_lead.py`:**
  - W2 against Q0n over hours 1–4, before Q0n's onset: temperature and
    wind RMSE within ± 0.2 (K, m/s) of Q0n.
  - W4 against the 06Z default run over hours 1–12: within ± 0.2 as well.
  - That is the damping, and the zone change, costing no skill before
    the default fails.

If W2 and W4 hold, the proposal is a new default of width 15, alpha 0.1
and C = 0.0064. It would stand on one calm case and one jet case. That
meets the two-case rule set for the relaxation change, and it would still
be checked on the first live cycles.

**Test W result (prompt 142).** C = 0.0064 gave ν = 5.85e4 m²/s on the Q
case and 5.38e4 on the 06Z case. Four runs side by side, torch × 8 each.

| Run | Case | Zone | 24 h | Points changing > 5 m/s in 15 min | Wall |
|---|---|---|---|---|---|
| W1 | Q (jet) | default | **completed** | at most 6 per interval, all 10–12 cells in (east zone boundary), from about 17 h; largest 13.7 m/s at 22.50–22.75 h | 16.8 min |
| W2 | Q | width 15, alpha 0.1 | **completed** | **none** in the whole run; max\|u\|, max\|v\| 28.3 and 33.1 at 23.5 h | 9.7 min |
| W3 | 06Z (calm) | default | **completed** | **none**; max 8.7 and 10.1 | 28.4 min |
| W4 | 06Z | width 15, alpha 0.1 | **completed** | **none**; largest change 0.5 m/s at 23.5 h | 28.1 min |

(The 06Z runs ran at 3.0 steps/s against 5.4–9.4 for the Q runs, with
four 8-thread runs sharing the machine.)

- **W1, no interior onset (completes, or fails only at the zone
  boundary): holds.** It completed. Every fast point is at edge distance
  10–12, the default zone's inner boundary, and none is inside.
- **W2, completes 24 h with no interval over 20 points: holds**, with no
  point at all.
- **W3, still fails at the 06Z south-west corner in 15–18 h: refuted.**
  It completed 24 h with no fast point. Divergence damping also removes
  the 06Z failure, which the budget had attributed to d.grad U and
  mixing near the ground. The 12 h budget described where the energy
  entered the mode, not what would stop it. A mode damped as divergent
  must have been divergent enough, and the prediction reasoned from the
  source term alone.
- **W4, completes 24 h: holds.**
- **Skill** (`tools/score_by_lead.py`; RMSE differences, damped minus
  undamped):
  - Q case, W2 against Q0n over hours 1–4, before Q0n's onset:
    temperature +0.00 to +0.04 K, u −0.03 to −0.06 m/s, v −0.06 to
    +0.11 m/s. **Within ± 0.2: holds.** At 6–7 h, as Q0n breaks up, W2
    is better by 0.47–1.65 (u) and 0.22–0.50 (v).
  - 06Z, W4 against the default run over hours 1–16: temperature −0.06
    to +0.02, u −0.02 to −0.17, v −0.02 to +0.04. **Within ± 0.2:
    holds.** The wind is slightly better throughout.
- **P-59 dominates the error.** In the Q case, which starts at 12Z and
  runs through the day, the temperature bias is −8.4 to −8.7 °C at hours
  6–8 (early afternoon) in both runs. The damping neither causes nor
  cures it.

**Decision (pre-registered in test W): new forecast defaults**, relaxation
width 15, alpha 0.1, divergence damping C = 0.0064. They are in
`src/forecast.py`, and `daily.sh` passes no flags, so the next cycle uses
them. `--relax-width 10 --relax-alpha 1 --div-damp 0` reproduces the old
model bit-identically (checked on the desktop). `src/test_forecast.py`
passes 11/11.

**P-60 is FIXED** on its case. **P-56 stays OPEN**: two of its cases pass,
but the two that first defined it have not been rerun. One is 12Z
2026-09-21, which was desktop-only and whose raw files went with the
wiped scratch data. The other is 18Z 2026-09-22, the first server cycle,
which died at 6.31 h before the slope limit existed.

**Test X (predictions first): the 18Z 2026-09-22 cycle rebuilt from its
archived raw observations** with today's ingest (slope-limited terrain),
24 h, torch. The original directory is copied aside first.

- **X1, new defaults:** completes 24 h with no 15-minute interval in which
  more than 20 points change by more than 5 m/s.
- **X0, old settings** (`--relax-width 10 --relax-alpha 1 --div-damp 0`):
  fails before 24 h, somewhere in 6–18 h. The earlier real cases failed
  at 6.31, 7.45, 13.7 and 16.3 h.
- **Skill:** X1 within ± 0.2 K and ± 0.2 m/s of X0 over X0's hours
  before its onset.
- If X1 fails, the new defaults are reverted to the old ones until the
  failure is understood, and P-56 gets the new case.

**Test X result (prompt 148).** First guess `standard_atmosphere`
(a calm case: initial max|u| 6.2 m/s); terrain slope-limited, 0–881 m.

- **X0, old settings, fails in 6–18 h: holds.** It diverged at 14.05 h,
  at the south-west zone corner as in the 06Z case. That is rows 10–13,
  columns 10–16, 9–12 cells in, levels L12–L19 over 700–840 m of
  terrain. The first points over 5 m/s appear about 9.75 h. It runs away
  at 11.75 h (102 points), with 290–700 points per interval after that
  (u and v).
- **X1, new defaults, completes 24 h with no bad interval: holds.** No
  point changed by more than 5 m/s in any 15-minute interval. max|u| and
  max|v| were 6.0 and 8.4 m/s at the end; 27.1 min at 3.1 steps/s.
- **Skill within ± 0.2 before X0's onset: holds.** X1 minus X0 RMSE is
  +0.03 to +0.10 K (temperature), −0.05 to +0.01 m/s (u) and −0.10 to
  +0.03 m/s (v), over leads 1–6.
- **But verification stopped at lead 6.25 h** (00:15Z), in both runs.
  Cause: the verification ASOS request sent dates only, and IEM ended it
  at 00Z (P-61). Fixed. Re-verifying X into fresh archives will give the
  full 24 h, and it confirms the fix.

**P-56 is FIXED.** All three real cases that can still be run now complete
24 h with the new defaults, and fail with the old ones (16.31, 7.45 and
14.05 h). The case that opened P-56, 12Z 2026-09-21, cannot be rerun. It
is to be reopened if a live cycle diverges.

**Test Y result (prompt 147): the calm cases are slow for a reason in the
state.** Each case ran alone for 1 h at torch × 8. The server had
3 users, load average 2–7 on 104 cores.

| Case | Normal | Denormals flushed | Result difference |
|---|---|---|---|
| calm (18Z 2026-09-22) | 2.9 steps/s | 3.0 | 0.0 |
| jet (12Z 2026-09-25) | 11.1 | 11.2 | 1.1e-13 |

- **Denormals are refuted:** flushing them changes nothing. Neither
  saved state contains a single subnormal value.
- **Machine load is refuted:** each run was alone on a quiet machine,
  and in test W the two cases differed while running at the same time.
- **Something in the calm state makes each step 3.8× more expensive.**
  The candidate is `dry_convective_adjustment`. It runs after every step
  and repeats its sweep, two Python loops over 20 levels, until no
  interface is unstable, up to 20 sweeps. The measured extra cost,
  about 255 ms a step, matches roughly 13–25 sweeps. Test Z counts the
  sweeps and times the adjustment (predictions: in the calm case a mean
  of 10 or more sweeps and at least 60 % of the step time; in the jet
  case 0 sweeps on most steps and under 10 %).

**Re-verification of test X (prompt 149): P-61 confirmed fixed, and P-62
found.** X was verified again into fresh archives, so the truncated
observation caches were not reused.

- X1 is now scored at every lead from 1 to 23 h, with 342–362
  temperature pairs per lead. X0 is scored to 14 h, its last snapshot.
  The prediction said leads 1–24: lead 24 is missing because the
  forecast file itself stops at 23.75 h.
- **P-62:** `run_forecast` compared `model.time`, a sum of 5 040 float
  steps that ends a few nanoseconds short of 24 h, with a 1e-9 s
  tolerance. So every run lost its final snapshot. The last step now
  always writes it, and a stub-model regression test was added
  (`test_forecast.py`, 12/12).
- X1 against X0 over leads 7–14: temperature RMSE +0.01 to +0.03, u
  −0.08 to −0.21 (better), v −0.01 to −0.10. Skill is unchanged or
  better throughout.
- **P-59, the night side:** bias rises from +1.3 °C (lead 5) to +6.2 °C
  (lead 17, 11Z), then falls to −1.4 °C (lead 23, 17Z).

**Test Z result (prompt 150): the convective adjustment is the cost, and it
never converges.**

| Case | Calls | Sweeps per call | Adjustment time | Unstable interfaces at each call |
|---|---|---|---|---|
| calm (18Z 2026-09-22) | 210 | **20 (the cap) every time** | 51.6 s of ~72 s stepping | mean 8.19e-3, max 9.16e-3 |
| jet (12Z 2026-09-25) | 229 | 0 every time | 0.2 s | 0 |

The prediction (calm: a mean of 10 or more sweeps and at least 60 % of the
time; jet: mostly 0 and under 10 %) holds. It also shows something the
prediction did not ask about: the calm case hits the cap on **every** call,
with a near-constant 0.8 % unstable at the start of each. So the
adjustment ends every step with instability left.

On the desktop the scheme converges on a column with 3 K noise, but only
after 69 sweeps (2.4 % still unstable at 20). Registered as P-63.

**Measurement Z2 (predictions first).** `tools/convection_check.py` on the
calm and jet 1 h states:

- **(a) Relaxation re-imposes it:** at least 80 % of the unstable
  interfaces lie in the relaxation zone (edge distance 14 or less).
- **(b) Interior needing more sweeps:** most lie deeper than that, and
  the adjustment converges (none left) with a higher cap, at more than
  20 sweeps.
- Either way, the jet state has none.
- Both can be partly true. The split by edge distance decides which fix
  comes first: adjusting the driving state once, or a converging and
  cheaper adjustment.

**Measurement Z2 result (prompt 159; run by the AI directly on the server):
all of the instability is in the relaxation zone.** The 1 h states from
test Z, `tools/convection_check.py`, one thread:

| Case | Unstable interfaces | Columns | Largest | By interface | By edge distance 0–4 / 5–9 / 10–14 / 15+ | Sweeps to converge |
|---|---|---|---|---|---|---|
| calm (18Z 2026-09-22) | 7.08e-3 (1 436 of 202 730) | 1 012 of 10 670 | 0.471 K | L18/19: 1 012, L17/18: 424 | 747 / 470 / 219 / **0** | 30 (the cap of 20 leaves 9.96e-4) |
| jet (12Z 2026-09-25) | 0 | 0 | — | — | — | 0 |

- **(a) holds, at 100 %** against the 80 % threshold: not one unstable
  interface lies deeper than 14 cells. The count falls with distance from
  the edge, as the relaxation weight does.
- The second half of (b) is also true: with a higher cap the adjustment
  converges, in 30 sweeps. So the cap of 20 is too low for this state, but
  the state is re-created every step, so raising the cap would only make each
  step more expensive.
- **Mechanism.** The observation-only runs have one driving frame, the
  analysis, held for the whole run. At 18Z that analysis carries an
  afternoon superadiabatic surface layer. Every step the adjustment mixes
  the edge columns to neutral, then the relaxation (weight up to 0.1 a step)
  pulls them back toward the unstable frame. The instability is re-created
  every step, for 24 hours, night included.
- **A side effect this implies, not yet measured.** The same cycle is a heat
  source. The relaxation restores the warm surface layer, the adjustment
  mixes that heat upward, and it repeats every step. The relaxation zone
  covers half the domain. So the cycle may contribute to the night warm
  bias of P-59 (+6.2 °C at lead 17 in this case). That is a hypothesis for
  test AA, not a finding.

**Fix.** Each boundary frame gets the same adjustment the model applies,
once, when the frame is built, run to convergence (`stabilise_frame` in
`src/forecast.py`, up to 1 000 sweeps; it raises an error if the frame is
still unstable). A frame that is already stable comes back unchanged,
array for array. After the fix, a weighted mean of two stable columns
(model and frame) is stable, so the relaxation can no longer create
instability. `--raw-boundaries` restores the old behaviour for comparison.
Two regression tests in `src/test_forecast.py` (14/14):

- a stable frame is returned bit-identical;
- relaxing an adjusted edge column toward the raw frame leaves 5.17e-2 of
  interfaces unstable, and toward the stabilised frame 0. The frame
  adjustment conserves each column's mass-weighted theta (drift 3e-16).

**Test AA (predictions first).** Torch × 8, new defaults, on the server.

- **AA1, calm (18Z 2026-09-22), 24 h, the fix on**, compared with X1 (same
  settings, fix off; archive `x1b`):
  1. The log reports 1 of 1 frames unstable, 1e-3 to 2e-2 of interfaces,
     converged in 100 sweeps or fewer.
  2. At least 8 steps/s over the run (X1: about 3). Without the adjustment
     cost, test Z's calm hour ran at about 10.
  3. `convection_check` on the 1 h snapshot finds exactly 0 unstable
     interfaces.
  4. Temperature and wind RMSE are not worse than X1 by more than 0.2 at any
     lead from 1 to 23.
  5. The P-59 side effect: mean temperature bias at leads 12–20, AA1 minus X1.
     A decrease of 0.2 K or more supports the heat-source hypothesis. A
     change within ±0.2 K is inconclusive, and an increase of 0.2 K or more
     refutes it.
- **AA2, jet (12Z 2026-09-25), 1 h, the fix on and `--raw-boundaries`:** the
  log reports 0 of 1 frames unstable, and the two outputs are bit-identical.
- **AA3, calm (06Z 2026-09-23), 1 h, the fix on:** at least 8 steps/s. If
  slower, `convection_check` shows whether it has an interior instability,
  which (b) would then have to fix.

**Test AA result (prompt 160; run and collected by the AI on the server,
commit 444fe9f): the frame fix is correct, and the slowdown has a second
cause.**

| Prediction | Result | |
|---|---|---|
| AA1-1: frame unstable, 1e-3 to 2e-2, ≤ 100 sweeps | 1 of 1, 9.14e-3, 34 sweeps | held |
| AA1-2: ≥ 8 steps/s | 3.2 (X1: 3.1) | **refuted** |
| AA1-3: 0 unstable interfaces at 1 h | 175 (8.63e-4), largest **7.73e-10 K**, all ≤ 14 cells from the edge, 174 of them at L17/18; 3 sweeps to clear | **refuted** — but the violations are now round-off size (Z2: 1 436, largest 0.471 K) |
| AA1-4: RMSE not worse than X1 by > 0.2 | T −0.01 to +0.02, u 0.00 to +0.01, v 0.00 to +0.01 | held |
| AA1-5: night bias, leads 12–20 | −0.00 to −0.01 K | inconclusive by the declared band; any heat-source effect is < 0.02 K |
| AA2: 0 frames unstable; bit-identical | 0 of 1; max diff 0.0 in u, v, theta, pi | held |
| AA3: ≥ 8 steps/s | 2.9; frame 2.04e-3 unstable, 34 sweeps; 1 h state 19 violations, largest 1.13e-9 K, all at the edge | **refuted** |

- The fix removes the real instability from the edges, costs nothing in
  skill, and leaves the stable jet case bit-identical. It stays in.
- The heat-source hypothesis of Z2 is not supported: the night bias moved
  by 0.01 K. The re-imposed layer was not what P-59 is made of.
- The speed is unchanged. That is the third explanation of P-63's cost to
  be refuted, after denormals and machine load (test Y). Under the stop
  rule, the next test splits the unstable interfaces into classes rather
  than guessing at a fourth mechanism.
- What the saved states cannot show: they are taken after the relaxation,
  not at the moment the adjustment runs. The violations left in them are
  near 1e-9 K, which points to one candidate. The adjustment leaves layers
  exactly neutral, each step moves them by round-off, and a tolerance of
  1e-10 K counts that as instability. That is a hypothesis; AB measures it.

**Test AB (predictions first; class-splitting).** A wrapper around
`dry_convective_adjustment` (as in test Z) records, at the start of every
call in the calm 1 h runs (18Z 2026-09-22 and 06Z 2026-09-23, torch × 8,
fix on): the unstable interfaces split by violation size (< 1e-8 K, 1e-8 to
1e-3 K, ≥ 1e-3 K), by edge distance (≤ 14 cells or deeper), and by level;
the sweeps used and the time; and, every 30th call, the sweeps needed to
converge with a cap of 1 000.

- **Round-off class:** at least 80 % of the unstable interfaces at the start
  of the calls are below 1e-8 K. The fix would be a physical tolerance in the
  adjustment (1e-10 K is 1 000 times below anything that matters).
- **Real class:** at least 80 % are 1e-3 K or larger. Then the step itself
  makes real instability every step, and the next measurement finds which
  tendency does it.
- Anything in between is mixed, and inconclusive for choosing a fix.
- Either way, sweeps and time per call should track the number of unstable
  interfaces. If the cap is hit on calls with only round-off violations,
  the sweep itself is not converging on them, which is a defect in the
  algorithm, not in the state.
- The prediction: **the round-off class.**

**Test AB result (prompt 160, same session; the predictions file on the
server is timestamped 21:15:43 EDT, before the run): not round-off, and the
defect is in the algorithm.**

| | calm 18Z 2026-09-22 | calm 06Z 2026-09-23 |
|---|---|---|
| sweeps per call | **20 (the cap) on all 210** | 20 on all 211 |
| adjustment time | 52.1 s of 74.8 s | 53.3 s of 76.3 s |
| unstable at call start | about 775 a call, steady | about 224 a call |
| size < 1e-8 / 1e-8 to 1e-3 / ≥ 1e-3 K | 0.9 % / **96.2 %** / 2.9 % | 0.0 % / **79.4 %** / 20.6 % |
| largest violation, first call / after | 4.92 K / about 1.6e-3 K | 2.57 K / about 1.7e-3 K |
| in the relaxation zone | 72.6 % | 100 % |
| by interface | 16/17: 14 469, 17/18: 75 538, 18/19: 72 798 | 17/18: 19 999, 18/19: 27 181 |
| left unstable after the call | every call (mean 8.84e-4) | every call (mean 9.34e-5) |
| sweeps to converge, cap 1 000 | 34 on the first call, then 23 | 34, then 24 |

- **The round-off prediction is refuted:** under 1 % of the violations are
  below 1e-8 K. By the declared bands the result is mixed (96 % in the
  middle band), so it does not choose between a tolerance and a tendency
  hunt.
- **What it shows instead.** The first call meets the real instability
  (4.92 K; in the 18Z case, 400 of the 1 858 interfaces are in the interior,
  where the initial state itself is superadiabatic) and mixes it to exactly
  neutral. From then on, every step moves those neutral layers by about
  1e-4 K. A steady few hundred interfaces become slightly unstable, and the
  sweep scheme needs 23–24 sweeps to clear them. At the cap of 20 it leaves
  some behind, and the next step starts the same way. The jet case never had
  an unstable layer, so it has no neutral layers, and it costs nothing.
- **So the cost is the algorithm, not the state.** A scheme that mixes only
  the currently unstable segments needs many sweeps on a column of chained
  near-neutral layers, because each mix can unsettle the interface next to
  it. That is candidate (b) of Z2.
- **Fix: pool-adjacent-violators (PAV).** One pass per column from the
  ground up merges each new layer into the block below while the block mean
  is higher (Ayer et al. 1955). It gives the stable state closest to the
  input in the mass-weighted least-squares sense: the least mixing that
  removes the instability.
- **An expectation that was checked and was wrong.** I expected PAV to reach
  the same state as the sweep scheme run to convergence. On the desktop it
  does not: 16 K apart on a column with 3 K noise, 2e-4 K on a near-neutral
  one. The sweep scheme merges every contiguous run of layers next to an
  unstable interface, across the stable interfaces in between, so it mixes
  more than it needs to. PAV therefore changes the model's answer as well as
  its cost, and test AC has to judge it on skill, not identity.
- The same desktop check reproduces the P-63 defect in a synthetic column: a
  neutral bottom block nudged by 1e-4 K. The 20-sweep cap leaves 2–9 % of
  interfaces unstable, and convergence takes 42–142 sweeps. PAV clears it in
  one pass, in 4–5 ms against 21–30 ms for 20 sweeps (1 600 columns).

## 2026-10-01 — P-63 fixed, a 20-cycle stability campaign, and a first sun

**Summary.**
- **P-63:** fixed (the frame fix of test AA and the PAV adjustment; the cause was found in tests AA and AB on 29 September). Calm cases now run
  3.3 times faster, with identical skill.
- **Stability:** the dry model completed 20 of 20 consecutive real cycles
  (26–30 Sep), 23 real 24 h forecasts in all since the new defaults.
- **Skill against persistence (P-59, P-68):** the dry model is worse than
  persistence at every lead, so it does not beat doing nothing. The prescribed
  surface heating brings temperature level with persistence (−0.24 K against
  no heating, 17 of 24 leads better than persistence). It is stable only with
  momentum mixing off (P-67).
- **Wind:** the verification compared a ~300 m wind with 10 m observations
  (P-68). With a 10 m operator, the model's real wind error remains (+0.34
  m/s u), and it grows with lead time.
- **Pipeline defects found and fixed:** cycling never used the previous
  forecast (P-64), the land mask (P-65), and a 110 times slower QC (P-66).

**For the collaboration study.** This was the first long block run without
the human: about 11 hours, after "keep going". It became possible because
the AI could reach the server directly (prompts 152–158). Every test still
had its predictions written before the run. Because the AI cannot commit,
the predictions were also left as timestamped files on the server
(`*_predictions.txt`, written by the same job before its first run). Two
of the AI's own claims were wrong and were caught by checks, not by the
human:
- that PAV reaches the converged sweep state (a unit test showed 16 K);
- a hand-copied AF summary (recomputed from the table: 1.578/1.599, not
  1.553/1.585).

Two failures were found by watching for anomalies rather than by tests:
- P-64 surfaced as a fallback message repeated in 18 logs;
- P-66 surfaced as a verification that was 600 times slower than its fetch.

**Test AC (predictions first; prompt 161, written 2026-10-01 before any
run).** `--conv-scheme pav` against the sweep scheme, both with the frame
fix, torch × 8. These runs use a development copy on the server (commit
444fe9f plus the uncommitted PAV and persistence changes, under
`scratch/dev`, data from the normal data root). The git checkout is not
touched, because the user commits.

- **AC1, calm (18Z 2026-09-22), 24 h, PAV** against AA1 (sweep):
  1. Completes 24 h.
  2. At least 8 steps/s (AA1: 3.2).
  3. Every hourly snapshot has exactly 0 unstable interfaces. The model's
     state is fully adjusted after each step, and the relaxation then mixes
     two stable columns.
  4. Temperature and wind RMSE are within ±0.2 of AA1 at every lead. PAV
     mixes less than the sweep scheme, so near-surface layers change, but
     only where the column was unstable.
- **AC2, jet (12Z 2026-09-25), 1 h, PAV:** bit-identical to AA2. That case
  never has an unstable interface, and PAV returns a stable state
  unchanged.
- **AC3, calm (06Z 2026-09-23), 24 h, PAV against AC3s, the sweep scheme,
  same settings:** both complete 24 h; PAV at 8 steps/s or more and the
  sweep scheme at 3.5 or less; RMSE within ±0.2 at every lead.
- **Persistence reference, for both calm cases:** the initial state held
  for 24 h (`--persistence`) and verified the same way. No prediction is
  needed: it is the baseline the model is scored against from now on.
- **Decision rule:** if 1–4 and AC2 and AC3 hold, PAV becomes the default.
  A skill change outside ±0.2 at any lead stops the adoption until it is
  understood.

**Test AC result (2026-10-01, run and collected by the AI): PAV holds on
everything scored.**

| Prediction | Result | |
|---|---|---|
| AC1 completes 24 h | yes | held |
| AC1 ≥ 8 steps/s | 24 h in 8.1 min, about 10.4 steps/s (AA1: 26.6 min) | held |
| AC1: 0 unstable interfaces in every hourly snapshot | 0 in 24 of 24 (the sweep run AC3s: unstable in 24 of 24, up to 9.4e-5) | held |
| AC1 RMSE within ±0.2 of AA1 | identical to the 0.01 printed, at every lead, for T, u and v | held |
| AC2 bit-identical to AA2 | max diff 0.0 in u, v, theta, pi | held |
| AC3 and AC3s complete; PAV ≥ 8, sweep ≤ 3.5 steps/s | 7.8 min (about 10.8 steps/s) against 27.2 min (about 3.1) | held |
| AC3 RMSE within ±0.2 of AC3s | re-verified after IEM 503 errors: identical to the 0.01 printed, at every lead, for T, u and v | held |

**The persistence reference (18Z 2026-09-22), the first time the model has
been scored against doing nothing:**
- **Temperature:** the model is better by 0.2–0.44 K at leads 5–20 and worse
  by 0.16–0.45 K at leads 1–4 and 21–24. Both errors are dominated by the
  missing diurnal cycle (P-59), RMSE up to 7.3 K for the model and 7.7 K for
  persistence.
- **Wind: persistence is better at almost every lead.** u RMSE model minus
  persistence rises from about 0 at lead 1 to +0.5 at 8–10 h and +1.07 at
  24 h; v is +0.3 to +0.6 from lead 3 on. The model's u bias grows to
  +2.4 m/s, against +1.4 for persistence. So the near-surface wind gets
  worse than the analysis as the dry model runs. This is one case; the
  campaign below tests whether it holds generally.
- **06Z 2026-09-23 against persistence is worse.** Temperature is 0.1–0.16 K
  better for 5 h, then worse at leads 7–23, by up to +1.30 K at lead 18. The
  model's next-afternoon cold bias is −6.81 K against −5.81 for
  persistence: without the sun, the model drifts colder than simply holding
  the night-time state. u is worse at leads 2–7 and 17–24 (up to +0.35), and
  v at leads 1–5 and 16–24 (up to +0.58).

**Decision: PAV is the default** (`--conv-scheme pav`; `sweep` reproduces the
old model). Every AC prediction held. P-63 is FIXED, and calm cases now
run at the speed of the jet case.

**Test AD (predictions first): a 20-cycle stability campaign.** Every 00,
06, 12 and 18Z cycle from 26 to 30 September 2026, ingested from
observations, in two chains (26–28 and 29–30 September; the first cycle of
each chain is a cold start, and the rest use the previous run's 6 h forecast
as in production). Each forecast is 24 h, PAV, new defaults, torch × 8, two
at a time. Each is verified with its persistence reference.

1. **Stability: all 20 forecasts complete 24 h,** with no divergence and no
   deadline stop. A failure is registered as a new problem, with its place
   and time from `tools/locate_growth.py`.
2. Every forecast takes 15 min or less of wall clock.
3. Every hourly snapshot of every run has 0 unstable interfaces.
4. All 20 ingests produce an analysis.
5. **Skill against persistence, averaged over the 20 cycles** (expectations
   from one case, not decisions): wind RMSE, model minus persistence, is
   above +0.2 m/s for u and v at leads 6–24 h. For temperature the
   difference is within ±0.5 K, which is inconclusive by design, because
   neither has a diurnal cycle.

**Where the wind error comes from (18Z 2026-09-22, model AC1 against
persistence, split by station position).** At the ~100 stations inside the
relaxation zone, the model is within +0.1 to +0.3 m/s (vector RMSE) of
persistence at every lead, which is expected, because the zone is relaxed
toward the frozen initial state. At the ~225 interior stations it is worse
by +0.6 to +0.9 at 4–10 h and by +1.1 to +1.5 at 20–24 h. The interior speed
bias is −0.6 to −0.9 m/s in the evening (persistence: +0.3) and −3.4 m/s the
next afternoon (persistence: −2.5). So the interior near-surface wind weakens
and does not recover by day. Vertical mixing is exactly zero above
Ri = 0.25, and nothing heats the ground. So by day nothing brings momentum
down to the lowest layer, while drag keeps taking it out. This is the wind
side of P-59.

**P-59 step 1: a prescribed diurnal surface heat flux (prompt 161).**
`src/dynamics/diurnal.py`. H = f·S0·τ·max(sin e, 0) + H_night over land,
0 over water, with:
- e the solar elevation (declination from Cooper 1969; the equation of time
  ignored);
- S0 = 1361 W/m² (Kopp and Lean 2011);
- τ = 0.75 clear sky, f = 0.2, H_night = −30 W/m².

It heats the lowest layer, and PAV then mixes heat and momentum upward,
which is a crude convective boundary layer. It uses only the date, the time
and the position, so the observation-only rule holds. The model is
unchanged unless `--surface-heating` is given. Tests (`test_diurnal.py`,
5/5):
- the Albany noon elevation is 46.64° at 16:55Z, as expected;
- the lowest-layer energy closes to 2e-16;
- the late-September daily mean is +16 W/m², with a peak of 117;
- in the model, land columns gain the prescribed heat to 2e-4 after mixing,
  and water changes by 1e-5 K;
- torch agrees with NumPy to 2e-13.

Limits: no clouds; land is "terrain above 0 m", so the Great Lakes are
heated; the relaxation zone damps the cycle near the edges.

**Test AE (predictions first).** PAV, new defaults, torch × 8, 24 h, with
and without `--surface-heating`, run from a second development copy
(`scratch/dev2`):

- **AE1** 18Z 2026-09-22 (baseline AC1), **AE3** 06Z 2026-09-23 (baseline
  AC3), and **AEj** the jet case, 12Z 2026-09-25 (baseline AEj0, PAV without
  heating). Persistence for each.
1. **Stability:** all complete 24 h.
2. **Temperature:** RMSE averaged over leads 1–24 lower than the baseline by
   at least 0.5 K in each case, and the swing of the bias (largest minus
   smallest lead-mean bias) smaller by at least 30 %. In the jet case the
   afternoon cold bias at leads 6–8 (Q case: −8.4 to −8.7 °C) improves by at
   least 3 K.
3. **Wind:** at the interior stations the speed bias by day (leads where
   local solar time is 10–17 h) is less negative by at least 0.3 m/s.
4. At least 6 steps/s.
5. Temperature RMSE better than persistence at 16 or more of the 24 leads,
   in each case.

Any improvement short of a threshold is called partial; a worsening refutes
it.

**Test AE result (2026-10-01): run with the faulty land mask (P-65, about
half the ocean heated as land), so it is a provisional reading, not the
test.** All three runs completed in 7.4–8.7 min.

| | 18Z 2026-09-22 | 06Z 2026-09-23 | jet, 12Z 2026-09-25 | prediction |
|---|---|---|---|---|
| T RMSE, mean over leads, base → heat | 4.18 → 4.01 (−0.17) | 4.22 → 3.55 (−0.67) | 4.84 → 4.05 (−0.79) | ≥ 0.5 lower: 2 of 3 |
| bias swing | 8.19 → 6.11 (−25 %) | 9.02 → 7.15 (−21 %) | 10.17 → 7.45 (−27 %) | ≥ 30 %: partial in all three |
| jet leads 6–8 cold bias | | | −8.3/−8.7/−8.7 → −6.1/−6.1/−5.8 | ≥ 3 K better: partial (1.9–2.3) |
| heat beats persistence | 15 of 24 leads | 14 of 24 | 24 of 24 | ≥ 16: jet only |
| interior daytime speed bias | −2.06 → −1.95 | −1.34 → −1.47 | +3.76 → +3.77 | ≥ +0.3: **refuted** |

- **Temperature:** the sun clearly helps. The 06Z case's next-afternoon cold
  bias goes from −6.8 to −4.4 K, and the jet case's afternoon error falls by
  more than 2 K. Where it hurts: the 18Z case is too warm at leads 4–9
  (late afternoon and evening, +0.3 to +0.7 K RMSE), and the 06Z case is too
  warm the next evening and night (leads 19–24, +1.2 to +2.5 K bias). The
  heat put into the mixed layer by day is not removed at night: the night
  flux cools only the lowest layer.
- **Wind: no change at all,** which refutes item 3. The verification takes
  the wind of the lowest model level, about 300 m above the ground, and
  compares it with 10 m observations. Momentum mixed down into a mixed layer
  changes the 300 m wind little. In the jet case the model is +3.8 m/s too
  fast by day, while persistence is +0.5. That is consistent with a 300 m
  wind scored against 10 m observations. The wind question needs an
  observation operator first; tuning the heating will not answer it.
- Parameters were **not** tuned on these three cases. AE2, with the correct
  mask, is the measurement.

**Two defects found while AE and AD ran (P-64, P-65).**
- **P-65:** the AE log reported "land 85%". The land mask was
  `terrain > 0` on the slope-limited run terrain, which spreads land
  heights over the sea (15.5 % exactly 0 m, against 31.9 % in the raw ETOPO
  grid). So AE heated and cooled about half the ocean. **The first AE results
  are contaminated and are reported as such;** AE2 repeats the test with the
  mask from the raw grid.
- **P-64:** every chained campaign cycle fell back from the previous
  forecast to the previous analysis or the sounding mean. Snapshots were
  stamped one step past the hour (6 h + 7.8 s), and the ingest asks for 6 h
  to within 3.6 s. So the documented cycling has never happened. Fixed in
  the working copy: whole steps per output interval, so snapshots land on
  the hour (the step becomes 17.06 s instead of 17.13), and a 72 s tolerance
  in the ingest.

**AE2 (predictions as for AE, items 1–5),** from `scratch/dev3` (444fe9f plus
every change of this session), with the raw-terrain land mask.

**AF (predictions first): does the forecast first guess help the analysis?**
Each chained campaign cycle (18 of them) is rebuilt from its archived
payloads (`--from-raw`) in a separate data root, now with the previous
run's 6 h forecast usable. The withheld-station score, surface temperature
at about 70 ASOS stations never assimilated, is compared with the
campaign's own analysis, which used the previous analysis or the sounding
mean.
- The forecast first guess is **better** if the mean withheld RMSE is lower
  by at least 0.1 K and lower in at least 12 of 18 cycles. It is **worse** if
  the mean is higher by at least 0.1 K. Anything else is no detectable
  difference.
- The expectation is **no detectable difference**: surface observations are
  assimilated on top of the first guess, and the withheld score is a
  surface score. At lead 6 the model and persistence scored within 0.3 K of
  each other in AC.
- This is one step of cycling. The previous forecasts come from non-cycled
  analyses.



---

**Test AD result, stability part (2026-10-01, 13:34 EDT): the dry model is
stable on 20 of 20 consecutive real cycles.**
- 1 held: 20 of 20 forecasts completed 24 h, with no divergence and no
  deadline stop.
- 2 held: 7.3–9.3 min of wall clock each, two at a time.
- 3 held: 0 unstable interfaces in all 480 hourly states.
- 4 held: all 20 ingests produced an analysis.

With the three cases of tests W, X and AC, that is 23 real 24 h forecasts
without a failure since the new defaults (d23e5fa) and PAV. The skill part
(item 5) waits for the verification windows to close.

One cycle, 12Z on 28 September, did use its previous forecast. Its log says
"previous_forecast …:+6h", because that run's 6 h snapshot happened to fall
within the 3.6 s window. That is P-64's mechanism seen directly.

**AE20 and a step control (predictions first).** Three cases are too few to
judge a physics change, so the heating is also run on all 20 campaign
cycles (`scratch/dev3`, correct land mask). Each is scored against the same
cycle without heating (the campaign forecast) and against persistence, on
the same observations (copied, `--report-only`). The predictions are those
of AE, items 1–5, applied to the RMSE pooled over the 20 cycles at each
lead, plus:
- in at least 15 of the 20 cycles, the heated run's T RMSE averaged over
  leads is lower than the unheated run's.

**Step control:** the 18Z 2026-09-22 case without heating from `dev3`, whose
step is 17.06 s instead of 17.13 (P-64), against AC1. RMSE within ±0.05 at
every lead.

**Verification was the bottleneck (P-66).** Each `verify.py` run took
10–14 min, almost all of it in an all-pairs Python buddy search in QC. The
replacement is a sorted-by-time NumPy search. On a real 24 h ASOS window
(50 170 observations) it takes 7.6 s against 849.1 s, with identical
verdicts (0 differ, 645 rejections each). At 14:36 EDT the fast version was
copied into `scratch/dev3` while the verification pool was running.
Verifications that start after that use it. Because the verdicts are
identical, the archives do not depend on which version matched them.

**AF result (prompt 161): no detectable difference, as expected.** With the
previous forecast as first guess, the 18 chained analyses have a mean
withheld-station RMSE of 1.578 K, against 1.599 K for the campaign's
(previous analysis or sounding mean): −0.021 K. The rebuild is lower in 10
of 18 cycles and higher in 7; one is identical, because it already used its
forecast. That is short of both bands (0.1 K and 12 of 18). By cycle the
difference runs from −0.253 K (27 Sep 06Z) to +0.250 K (28 Sep 00Z). So fixing P-64
neither helps nor hurts the surface analysis measurably. It restores the
documented design, and its effect aloft is not measured by this score.

**A stability failure with the heating (P-67).** In AE20, the heated run
of 00Z 2026-09-28 diverged at 22.46 h. The same cycle without heating
completed in AD. The growth is in the lowest levels over northern New York
and New England, and the unheated run has the same noisy patch, but
bounded. So the heating is not adopted, and it stays off by default.

**Test AG (predictions first; class-splitting, the 28 Sep 00Z cycle, 24 h,
heating on unless stated; each change made alone):**
- **AG1** `--no-conv-momentum` (the adjustment mixes theta only). If (a), the
  instant momentum mixing, drives it: completes 24 h with max|u| ≤ 20 m/s.
- **AG2** `--sh-night 0` (no night cooling). If (b), the night cooling over
  terrain, seeds it: completes 24 h.
- **AG3** `--sh-fraction 0` (night cooling only, no daytime heating). If (b)
  is enough by itself: max|u| above 25 m/s or divergence. If it completes
  below 20 m/s, the night cooling alone is harmless.
- If AG1 and AG2 both still fail, it is neither alone, and the next step is
  a closer look at the state over the patch, not a fourth guess.

**Test AG result (2026-10-01, 15:25 EDT): the night cooling and the
momentum mixing together, neither alone.**

| run (28 Sep 00Z, 24 h) | night cooling | day heating | momentum mixed | result |
|---|---|---|---|---|
| AE20 | yes | yes | yes | diverged 22.46 h |
| AG1 | yes | yes | **no** | completed, max\|u\| ≤ 15.6 m/s |
| AG2 | **no** | yes | yes | completed, max\|u\| ≤ 15.4 m/s |
| AG3 | yes | **no** | yes | **diverged 17.48 h** |

- AG1 holds (a) as necessary, and AG2 holds (b) as necessary. AG3 shows that
  the night cooling with momentum mixing fails without any daytime heating,
  and sooner, because the cooling then runs all day.
- So the failure is an interaction. Steady cooling of the lowest layer over
  land builds strong low-level cold pools. Where the adjustment fires in or
  around them, it mixes u and v instantly through the block. Without the
  cooling the mixing is harmless (AG2, and the 20 unheated campaign runs).
  Without the mixing the cooling is harmless (AG1).
- Momentum mixing also had no measurable effect on the verified wind in AE.
  The candidate configuration is therefore heating with momentum mixing off.
  That is a choice of configuration, not yet an explanation of the growth;
  the cold-pool dynamics remain unexplained.

**Test AH (predictions first): the candidate configuration on all 20
campaign cycles.** `--surface-heating --no-conv-momentum`, PAV, torch × 8,
verified on the campaign's observations.
1. All 20 complete 24 h.
2. Against the unheated runs (AD), with T RMSE pooled over the cycles at each
   lead: mean over leads lower by at least 0.3 K; better in at least 15 of 20
   cycles; bias swing smaller by at least 20 %.
3. Against persistence: T RMSE better at 16 or more of 24 leads.
4. The wind, against the unheated runs: within ±0.2 m/s at every lead.
   Momentum is not mixed, and the verified level is about 300 m.

**Test AD result, skill part (20 cycles, pooled at each lead).**

| | T RMSE, model − persistence | u | v |
|---|---|---|---|
| leads 1–24 | +0.04 to +0.54 K, mean **+0.21**; better at 0 of 24 leads | mean +0.62 m/s; leads 6–24 +0.68 | mean +0.60; leads 6–24 +0.60 |
| by cycle (mean over leads) | better in 3 of 20 | better in 2 of 20 | better in 5 of 20 |

- Item 5 for the wind holds as predicted: the model is worse than
  persistence by more than 0.2 m/s at leads 6–24. For temperature the mean
  (+0.21 K) is inside the ±0.5 K band declared inconclusive. But the model is
  never better at any lead: **the dry model, without a sun, does not beat
  doing nothing at the surface.**
- **The wind numbers are dominated by the observation operator (P-68).**
  Persistence at lead 1 is the analysis, which assimilated these
  observations. Yet its wind speed at the stations is 1.7–2.3 times the
  observed 10 m speed in nearly every cycle: 11.1 against 4.8 m/s at
  27 Sep 12Z, and 2.0 against 0.9 at 30 Sep 06Z. That gives a pooled lead-1
  bias of −2 m/s in u and v (the flow was mostly northeasterly). The
  verification takes the lowest model level, about 300 m above the ground,
  and scores it against 10 m anemometers. Until the operator brings the
  model wind down to 10 m, a wind comparison measures the operator first.

**AE20 result (heating on all 20 cycles; correct land mask; momentum mixed;
19 scored, one diverged, P-67):**
- **T against no heating:** better at 21 of 24 leads, by −0.247 K on
  average (leads 6–24: −0.263), and in 16 of 19 cycles. The largest gains are
  −0.40 to −0.45 K at leads 6–13.
- **T against persistence:** better at 17 of 24 leads, by −0.027 K on
  average. It wins by up to 0.33 K at leads 4–15 and loses by up to 0.60 K at
  leads 19–24.
- **Wind:** u is +0.15 m/s worse, v is −0.08 better.

**AE2 result (three cases, correct mask, momentum mixed).**
- T RMSE averaged over leads: −0.21 K (18Z 09-22), −0.63 K (06Z 09-23) and
  −0.81 K (the jet case).
- Bias swing: −23 %, −22 % and −24 %.
- Leads better than persistence: 16, 13 and 24 of 24.
- The jet case's cold bias at leads 6–8 improves by 2.0–2.3 K.
- Interior daytime speed bias: −2.06 → −1.87, −1.34 → −1.48, +3.76 → +3.70.

Against the predictions: item 2 holds in two of three cases for the RMSE
and is partial for the swing (all short of 30 %) and for the jet cold bias
(short of 3 K). Item 3, the wind, is refuted. Item 5 holds in two of three.
The correct land mask changed the result very little.

**Step control:** `dev3` (17.06 s step) against AC1 (17.13 s): identical to
the 0.01 printed at every lead. P-64's step change does not move the scores.

**Test AH result (2026-10-01, 16:58 EDT): heating with momentum mixing off,
20 cycles.**

| prediction | result | |
|---|---|---|
| 1. all 20 complete 24 h | 20 of 20, 7.4–9.3 min, including 28 Sep 00Z, which diverged with momentum mixed | held |
| 2. T RMSE against no heating, mean over leads, ≥ 0.3 K lower | −0.243 K (leads 6–24 −0.259); better at 21 of 24 leads | partial |
| 2. better in ≥ 15 of 20 cycles | 17 of 20 | held |
| 2. pooled bias swing ≥ 20 % smaller | 0.86 → 0.40 K (−53 %); pooling all four start times blurs the diurnal phase, so this is the weakest of the checks | held |
| 3. T better than persistence at ≥ 16 of 24 leads | 17 of 24 (mean −0.032 K): better by up to 0.33 K at leads 4–16, worse by up to 0.6 K at 19–24 | held |
| 4. wind within ±0.2 m/s of no heating at every lead | v: −0.01 to −0.15, better at all 24 leads; u: −0.02 to +0.25, outside the band at leads 21–23 | refuted, narrowly |

- With momentum mixed (AE20, 19 cycles), the scores are the same to about
  0.01 K: the momentum mixing was not doing anything measurable for skill,
  and it was what made 28 Sep 00Z diverge.
- `docs/diurnal_heating_skill.png` shows pooled T RMSE and bias by lead for
  persistence, the dry model, and the dry model with heating (19 cycles
  common to all runs).
- **The bias by lead tells the remaining story.** The unheated model drifts
  cold (to −0.5 K pooled). The heated model drifts warm (+0.1 to +0.5 K), with
  peaks at leads 5–6, 11–12 and 17–18. Those are the leads where one of the
  four start times reaches mid-afternoon. The heating is now a little too
  strong by day and too weak at removing heat by night (AE). The
  late-lead loss to persistence (leads 19–24) is where both lack what the
  real atmosphere does overnight.
- **Status:** the heating improves temperature, but not by the 0.3 K
  threshold. It is stable only with momentum mixing off. It stays off by
  default, and `--surface-heating --no-conv-momentum` is the candidate.
  Adoption is the user's call.

**P-68 step 1: a 10 m wind operator.** `verify.py --wind-operator log10m`
scales the lowest-level wind by ln(10/z0)/ln(z1/z0), where z1 is the
level's height above the model ground and z0 = 0.1 m. That is the model's
own neutral drag profile. The default is unchanged (`lowest`). The original
value is kept in each match (`forecast_lowest_level`, `wind_factor`). Test
(`test_sigma_operator.py`, 8/8): z1 is 237 m above flat ground (factor
0.593) and 226 m above 800 m terrain (0.596). Above about 50 m a neutral log
law over-reduces, so this is a first operator, not a surface-layer scheme.

**Test AI (predictions first): AD, its persistence and AH re-scored with
the 10 m operator** (`--report-only` on the archived observations, 20
cycles each).
1. Persistence at lead 1, pooled: the ratio of mean forecast to observed
   wind speed falls from its current value (to be printed alongside) to
   0.9–1.4.
2. Persistence at lead 1: u and v RMSE each lower by at least 0.5 m/s.
3. Temperature scores identical to the old archives, because the operator
   touches only wind.
4. No prediction on whether the model still loses to persistence in wind.
   That is the question the operator lets us ask.

**Test AI result (2026-10-01, 17:15 EDT): the 10 m operator removes most
of the wind error that was not the forecast's.**

| prediction | result | |
|---|---|---|
| 1. persistence lead-1 speed ratio (mean forecast / observed) falls to 0.9–1.4 | 1.86 → 1.10 (lead 12: 2.04 → 1.21; lead 24: 2.22 → 1.32) | held |
| 2. persistence lead-1 u and v RMSE each lower by ≥ 0.5 m/s | u 3.49 → 2.28 (−1.21); v 3.11 → 1.94 (−1.16); lower at all 24 leads in all 20 cycles | held |
| 3. temperature identical to the old archives | +0.000 at every lead in every cycle | held |

**What the operator then shows (item 4, no prediction).**
- **The dry model still loses to persistence in wind, by much less:** u
  +0.336 m/s (leads 6–24 +0.387), v +0.145 (+0.157). It is better in 2 and
  5 of 20 cycles.
- **Heating hardly changes the wind:** u +0.054, v −0.022 against the
  unheated model.
- **The model's near-surface wind strengthens with lead time.** Its speed
  ratio is 1.17 at lead 1, 1.39 at 12 h and 1.51 at 24 h, against 1.10,
  1.21 and 1.32 for persistence. Persistence's ratio also grows, so part of
  it is in the observations or the analysis. The extra growth in the model,
  +0.07 at lead 1 and +0.19 at 24 h, is a model error. This contradicts
  the 18Z 09-22 interior analysis, where the model wind weakened, so it
  varies with the case and needs its own look.
- **Recommendation:** make `--wind-operator log10m` the verification
  default. That is a change to the scorer, and it is the user's decision.
  With the old operator, wind scores measure the 300 m-to-10 m speed ratio
  more than the forecast.

## 2026-10-02 — The user's decisions, the night wind, and a ground temperature

**Summary.**
- **Decisions:** the goal is now a convection-allowing model
  (`docs/CAM_DESIGN.md`, a draft for review). Verification uses `log10m` by
  default. The P-64 fix stays. Dry fixes come first.
- **P-69 (night wind):** the model's 10 m wind excess was a night-time
  problem: a ground with no temperature and a neutral diagnosis. A
  force-restore ground temperature plus a similarity operator fix most of it
  (test AK).
- **P-67:** traced by a closed budget to the dispersion of the second-order
  centred advection at a northern-Vermont low-level jet. Third-order
  upwind-biased advection removes it with no change in skill (test AL).
- **The candidate configuration (test AM):** upwind3, the land surface and a
  land/sea roughness. On the 20 campaign cycles it beats persistence in 2 m
  T (−0.29 K), u (−0.045 m/s) and v (−0.069 m/s), pooled, and AD by 0.50 K,
  0.38 and 0.22 m/s. Late-lead temperature is still worse than persistence
  (no clouds). Nothing is switched on by default; the user decides.

**Context (prompts 162–163).** The user pulled c7380aa, said "keep going",
and asked for every approval needed before leaving for class. Six questions
were put to them, and they decided:

| question | decision |
|---|---|
| make `--surface-heating --no-conv-momentum` the default? | **no**: "we want a convective allowing model" |
| the follow-up: plan toward a 3 km non-hydrostatic CAM? | **yes**; a design and roadmap first, and no large rewrite until they have reviewed it |
| `--wind-operator log10m` as the verification default? | **yes** (done: `verify.py`) |
| keep the P-64 cycling fix? | **yes** |
| next: dry fixes, moisture, or both? | **dry fixes first** |
| run server jobs autonomously; set `NWP_BACKEND`/`NWP_THREADS` in `~/.bashrc`? | **yes** to both (done) |

On the server, all six test suites passed at c7380aa (convection 12/12 with torch).

**The CAM roadmap** is in `docs/CAM_DESIGN.md`, a draft for the user's
review. It covers what a convection-allowing model needs, what NWP1 keeps,
replaces and drops, a cost table against the 1.5 h budget, and stages S0–S7
with a gate test at each.
- A 3 km model over the whole domain is estimated at about 10–13 h per 24 h
  on today's settings, which does not fit.
- A 3 km nest of 450–600 km inside the 12 km model, or 4 km over the whole
  domain with float32 and more threads, might fit.
- These numbers are scaled from the measured 8.3 min run with an assumed cost
  factor, so stage S3 measures them before the grid is chosen.

**Test AJ (diagnostic; predictions first, `scratch/aj/aj_predictions.txt`): where does the near-surface wind grow?**
The model's wind at leads 6/12/18/24 h was compared with the analysis valid
at the same time: the converted analysis in the later cycle's
`forecast_persist.npz`. The AD campaign gave 70 pairs.

| prediction | result | |
|---|---|---|
| 1. domain lowest-level model/analysis ratio ≥ 1.05 at 24 h | 1.15 (6 h), 1.15, 1.22, 1.33 (24 h) | held |
| 2. drag-roughness hypothesis: interior land > 1.05 at 24 h, interior sea ≤ 1.00 | land 1.35; **sea 1.55**, the largest anywhere | **refuted** |
| 3. the excess is in the lowest levels (lowest > 2 × level 10) | lowest +1.41 (6 h) to +2.31 m/s (24 h); level 10 −0.11 to +0.31; mid-troposphere −0.7 to −1.3 | held |

One z0 everywhere cannot explain a growth that is largest over the sea, where
z0 = 0.1 m is about 500 times too rough. The edge zone grows least (1.03 to
1.19): it is relaxed toward the frozen initial state.

**The time-of-day split (no prediction; `tools/phase_scores.py` on the AD,
persistence and AH archives).** This changes what the problem is.

| 10 m wind, speed ratio forecast/observed | day (10–18 EDT) | night (22–06 EDT) |
|---|---|---|
| model (AD), leads 1–6 / 7–12 / 13–18 / 19–24 | 1.13 / 1.10 / 1.02 / 1.17 | 1.37 / 1.64 / 1.81 / 1.74 |
| persistence | 1.04 / 0.90 / 0.83 / 1.05 | 1.12 / 1.44 / 1.69 / 1.47 |

| 2 m T bias (K), all leads | day | night |
|---|---|---|
| model (AD) | −2.30 | +1.54 |
| persistence | −2.17 | +1.69 |
| model + prescribed heating (AH) | −1.35 | +1.49 |

- The model's mean 10 m speed is nearly the same day and night (3.9 m/s).
  The observed speed falls from 3.6 by day to 2.2–2.7 at night.
- By day the ratio is about 1.1. The "growth with lead time" found in test AI
  is mostly the share of night-time verification hours changing with lead,
  plus a night ratio that rises as the night goes on.
- The same holds for temperature. The prescribed heating fixed part of the
  day and none of the night.
- Both the wind and the temperature errors point at the same missing physics:
  the ground has no temperature. It does not cool at night, the surface layer
  never becomes stable, the drag stays neutral, and both 10 m operators assume
  a neutral profile.
- **Registered as P-69.**

**Test AJ2 (class-splitting; predictions first, `verification_tests/aj2/aj2_predictions.txt`).**
Four variants were run on 8 cycles (27 Sep 00Z to 28 Sep 18Z): control,
`--no-mixing`, `--no-conv-momentum`, and `--z0 1.0` (neutral drag × 2.0).
The aim is to see which removes the lowest-level excess over the analysis.
`forecast.py` gained `--z0` and `--no-drag` for this.

**Test AJ2 result (2026-10-02, 11:37 EDT).** 32 runs, in 9.6 min or less
each. The excess is the model's lowest-level wind minus the analysis valid
at the same time, in the interior, over the 8 cycles.

| variant | completed | lowest-level excess 6 h / 24 h (m/s) | u, v RMSE vs control (pooled, log10m) |
|---|---|---|---|
| control | 8/8 | +2.35 / +4.54 | — |
| `--no-mixing` | **3/8** (diverged at 6.5–22.1 h) | +2.40 / +4.03 (3 runs) | u −0.13, v +0.12 (3 cycles) |
| `--no-conv-momentum` | 8/8 | +2.33 / +4.52 | −0.01, −0.01 |
| `--z0 1.0` (neutral drag × 2.0) | **7/8** (28 Sep 06Z diverged at 18.99 h) | **+0.79 / +3.05** | **u −0.56, v −0.26** (24 of 24 leads, 7 of 7 cycles); T −0.08 K |

| prediction | result | |
|---|---|---|
| 1. control reproduces AD | identical to ±0.000 at every lead | held |
| 2. no mixing cuts the 6 h excess to ≤ 50 % | +2.40 against +2.35 | **refuted**; Richardson mixing is not the source |
| 3. no convective momentum within ±25 % | +2.33 against +2.35 | held |
| 4. z0 = 1.0 lowers the 24 h excess by 25–60 % | −33 % (−66 % at 6 h) | held |
| 5. z0 = 1.0 lowers pooled u and v RMSE | u −0.56, v −0.26 m/s | held |

**Stop rule, applied as written.** At 24 h no variant removes ≥ 50 % of the
excess. At 6 h the stronger drag removes 66 %. The 24 h picture is not a
surface one: on these 8 cycles the lower troposphere (levels 13–18) is
+2.5 m/s faster than the analysis by 24 h. So the late growth is in the
resolved dynamics or the frozen boundaries, and the next step there is a
momentum budget, not more parameters.
- These 8 windy cycles grow faster than the 20-cycle average: domain ratio
  1.86 at 24 h, against 1.33.

**What the drag result says.** The wind bias in u and v is negative at every
lead (u −1.6 to −2.8 m/s), while the speed is too high. The model wind is
too fast and turned, which is what too little surface friction does: too
little flow across the isobars toward low pressure. Doubling the drag
removes about a third of the u bias.
- That is the largest wind improvement of any change so far.
- No roughness is adopted from it. One value on 8 cycles, with one more
  divergence, is not a basis for a default. A land-use roughness (forest
  about 1 m, open sea about 0.0002 m, as tabulated in `surface.py`) is the
  physical version, and it is a candidate after test AK.

**The 28 Sep divergences share a place (`tools/locate28.py`).** The z0 = 1.0 run
of 28 Sep 06Z and the P-67 heated run of 28 Sep 00Z start the same way. The
fastest wind is at the second or third level from the ground, at
45.2–45.4 N, 72.2–72.8 W: northern Vermont, 19–21 cells from the northern
edge, over 170–310 m terrain between the Green Mountains and the border.
- The z0 = 1.0 run grew from 17 m/s to 41 m/s between 15 and 18 h, at
  21Z–00Z 28/29 Sep.
- The P-67 run grew from 22 to 31 m/s at 17 h (17Z) and from 33 to 51 m/s
  at 21–22 h.
- The control has its fastest low-level wind in the same box from 13 to 19 h
  (16–21 m/s) and survives.
- So P-67 is not specific to the heating. It is a weak spot in this place on
  this day, and three different changes push it over. Recorded under P-67.

**The P-67 tendency budget (`tools/budget28.py`, from the saved hourly states; the budget closes to 0.00).**
At each hour the fastest wind in the hot-spot box was taken, and every term
of the momentum tendency there was projected on the wind direction. The
point is a strict local speed maximum in its 3 × 3 × 3 neighbourhood at
every hour.

| run, hour | speed | horizontal advection | vertical advection | pressure gradient | mixing | hyperdiffusion |
|---|---|---|---|---|---|---|
| z0 = 1.0, 13 h | 15.1 | +7.5 | +5.0 | −4.9 | −0.5 | −2.0 |
| z0 = 1.0, 16 h | 24.0 | **+24.3** | −0.5 | −3.3 | 0.0 | −3.6 |
| z0 = 1.0, 18 h | 41.5 | **+66.0** | +44.7 | −23.5 | −51.7 | −7.1 |
| control, 17 h | 21.3 | +24.1 | +2.5 | −9.8 | −12.8 | −2.6 |
| control, 24 h | 19.2 | +27.3 | −1.0 | −5.8 | −0.3 | −2.0 |

(m/s per hour along the wind.)
- At a local maximum of speed, horizontal advection cannot raise the speed
  in the continuous equations: the along-wind part of (V·∇)V is V·∇|V|,
  which is zero there.
- The model's term is +7 to +66 m/s per hour. That is the dispersion error
  of the second-order centred advection (`_horiz_adv`), which is
  non-dissipative and leaves 2-dx structure to the hyperdiffusion. The
  hyperdiffusion removes 2–7 m/s per hour.
- The pressure gradient and the mixing oppose the growth. In the control
  they win; with stronger drag (more shear under the maximum) or with heating,
  they do not.
- This turns P-67 from "heating plus momentum mixing" into a numerical
  property of the core's advection, made visible by a strong low-level jet
  over northern Vermont.

**Test AL (predictions first, `verification_tests/al/al_predictions.txt`).**
This tests `--advection upwind3`: third-order upwind-biased horizontal
advection (Wicker and Skamarock 2002), written as the fourth-order centred
difference plus a |u|-scaled fourth-derivative damping. It is dissipative
only at the shortest scales, and it is the family the S1 non-hydrostatic
core would use.
- `test_advection.py` 4/4:
  - order 2.00 against 2.99;
  - one revolution of a 3-cell bump leaves a minimum of −0.286 for centred2
    and −0.024 for upwind3;
  - the 2-dx wave is invisible to centred2 and damped by upwind3.
- The test runs the 20 campaign cycles plus the two P-67 failures with the
  new scheme. It starts after AK.

**Test AL result (2026-10-02, 14:50 EDT).**

| prediction | result | |
|---|---|---|
| 1. 20 of 20 complete | 20 of 20 | held |
| 2. both P-67 failures complete, with the hot-spot wind under 30 m/s every hour | both complete 24 h. The fastest wind anywhere in the domain is 15.6 m/s (z0 = 1.0, 06Z) and 14.9 m/s (heated, 00Z), against 41.5 and 51.1 before the divergence with centred2 | held |
| 3. skill within ±0.05 K and no worse than +0.10 m/s | T +0.002 K, u +0.005, v −0.000 m/s | held |
| 4. run time ≤ 25 % longer | mean 9.8 min against about 8.3 (+18 %; max 11.3 against 9.3). AL ran 3 streams to AD's 2, so part of this is contention | held |
| 5. 24 h lowest-level excess over the analysis down ≥ 10 % | +2.32 against +2.31 m/s | **refuted** |

- **P-67 is fixed by the advection scheme.** The two configurations that
  diverged now run 24 h with no hot spot, and scores do not change. The
  20-cycle low-level wind excess is not an advection effect (prediction 5),
  which agrees with AJ2: that excess is surface physics and diagnosis (P-69).
- **Recommendation:** make `--advection upwind3` the default. It is a change
  to the core's numerics, so it is put to the user rather than switched
  silently. It is also the CAM roadmap's horizontal scheme.

**Test AM (predictions first, `verification_tests/am/am_predictions.txt`).**
The combined candidate on the 20 cycles: upwind3, the land surface, and a
land/sea roughness map (z0 = 1.0 m over land and 0.0002 m over water, the
forest and open-sea values in `surface.py`'s table, not fitted). With the
hot spot gone, AJ2's drag result can be tried without its divergence.
`--z0-land`/`--z0-sea` were added (`test_land_surface.py` 9/9: a map gives
each region its scalar drag).

**Test AM result (2026-10-02, 16:30 EDT).** All 20 runs completed, in
9.7–11.2 min each.

| prediction | result | |
|---|---|---|
| 1. 20 of 20 complete | 20 of 20 | held |
| 2. against AK (similarity): u RMSE down ≥ 0.20, v down ≥ 0.10 m/s | u **−0.274** (24 of 24 leads, 20 of 20 cycles); v **−0.051** (23 of 24 leads) | u held; v **refuted** (the right sign, half the size) |
| 3. T within ±0.10 K of AK | −0.023 K | held |
| 4. against persistence: u no worse than +0.05; v better | u −0.045 (15 of 24 leads); v −0.069 (22 of 24) | held |
| 5. speed ratio day ≤ 1.15, night 0.80–1.00 | 1.06 / 0.89 | held |

| 20 cycles, pooled RMSE, mean over leads | 2 m T (K) | u (m/s) | v (m/s) |
|---|---|---|---|
| AM against AD (dry, no sun) | −0.499 (19 of 24 leads) | −0.381 (24 of 24, 20 of 20 cycles) | −0.215 (23 of 24, 19 of 20) |
| AM against persistence | **−0.288** (18 of 24) | **−0.045** (15 of 24) | **−0.069** (22 of 24) |

- **The first configuration that beats persistence in all three surface
  variables, pooled over the 20 campaign cycles.** The margins in wind are
  small. Temperature loses at the last leads: at 23–24 h it is about 0.8–1 K
  worse than persistence, because the diurnal swing is too large with no
  clouds (day bias +1.32 K, night −2.13 K at leads 19–24).
- The combination is upwind3 + land surface + land/sea roughness. Every
  constant in it comes from a table or the literature. None was fitted to
  these cycles.
- Figure: `docs/am_skill_by_lead.png` (pooled RMSE by lead for persistence,
  AD, AK and AM; `tools/lead_scores.py`).
- These are 20 cycles in one week of late September. The candidate needs
  more weeks before anyone calls it an improvement in general. It has not
  been run on a stormy or a cloudy week on purpose.


**P-59 step 2: a ground temperature with its own energy budget (`src/dynamics/land_surface.py`).**
This is the force-restore scheme of Deardorff (1978): clear-sky sunshine,
longwave in and out, sensible heat by a bulk formula, latent heat as a fixed
share (EF = 0.5) of daytime net radiation (the model is dry), and a deep-soil
restoring term. The water temperature is held fixed. The drag now sees the
ground through the Louis (1979) momentum function, which weakens at night
without switching off.
- Every constant is a literature or climatological value, listed in the
  module. Nothing was fitted to these cases.
- One choice came from a measurement. Louis's heat function on the stable
  side collapsed the night flux to −1.1 W/m² in the column test, a runaway
  decoupling. So the heat flux uses the momentum function's long tail there
  (−13.8 W/m² in the same test). It is recorded as a choice, not a reference.
- `test_land_surface.py` 8/8:
  - the constants;
  - the Louis forms;
  - a night column (ground 4.7 K below the air);
  - a day column (peak H 195 W/m², Rn 462);
  - water held fixed;
  - drag stronger over a warm ground and weaker but not zero over a cold one;
  - coupled to the model, a night cools the lowest layer.
- `--land-surface` writes the ground temperature (`tg`) with each snapshot.
  A persistence run with `--land-surface` holds the ground temperature too.
- It is off by default. It replaces the prescribed heating rather than adding
  to it.

**The surface-layer operator (`src/verification/surface_similarity.py`, `verify.py --surface-operator similarity`).**
This is Monin–Obukhov similarity between the ground and the lowest level:
Paulson (1970) unstable, Beljaars and Holtslag (1991) stable. It gives the
10 m wind and the 2 m temperature. It is how most NWP models diagnose them,
here from a lowest level at about 230 m rather than 10–30 m, which is an
extrapolation and is stated as such.
- With a neutral ground it reproduces the `log10m` factor exactly.
- At night it takes 0.2–0.33 of the lowest-level wind instead of 0.59.
- `test_surface_similarity.py` 6/6.

On the server, a 3 h smoke run (27 Sep 12Z, torch × 2) completed. The ground
ran 2.5 K above the lowest-level air temperature at 13Z and 3.9 K at 15Z.
That includes about 1.5 K of lapse over 230 m.

**Test AK (predictions first, `verification_tests/ak/ak_predictions.txt`).**
The land surface on all 20 campaign cycles, scored both ways (standard, and
similarity), with a persistence that holds its ground temperature.

**Test AK result (2026-10-02, 13:25 EDT).** All 20 runs completed, in
7.8–10.2 min each.

| 20 cycles, pooled | 2 m T bias day / night (K) | 2 m T RMSE (K) | 10 m wind ratio day / night | wind vector RMSE (m/s) |
|---|---|---|---|---|
| persistence (standard) | −2.17 / +1.69 | 3.56 | 0.95 / 1.41 | 3.32 |
| AD, no sun (standard) | −2.30 / +1.54 | 3.76 | 1.10 / 1.63 | 3.67 |
| AK, land surface, standard operators | −0.98 / +1.77 | 3.49 | 1.07 / 1.82 | 3.96 |
| **AK, land surface, similarity operators** | **+0.51 / −1.33** | **3.25** | **1.22 / 0.92** | **3.48** |
| persistence holding its ground temperature (similarity) | −1.98 / +1.87 | 3.59 | 0.61 / 0.89 | 3.00 |

| prediction | result | |
|---|---|---|
| 1. 20 of 20 complete | 20 of 20 | held |
| 2. standard: day bias up ≥ 0.8 K; night changes < 0.5 K; RMSE down ≥ 0.15 K | +1.32 K; +0.23 K; −0.256 K (23 of 24 leads, 17 of 20 cycles) | held |
| 3. standard wind: night ratio above 1.63 | 1.82 | held (the drag weakens at night, as designed) |
| 4. similarity T: night bias within ±0.75 K; day within ±1.0 K; RMSE ≤ 3.40 K | night **−1.33**; day +0.51; RMSE 3.25 | night **refuted**; the others held |
| 5. similarity wind: night ratio ≤ 1.30; day 0.90–1.25; vector RMSE ≤ 3.50 | 0.92; 1.22; 3.48 | held |
| 6. similarity T RMSE below persistence's 3.56 K | 3.25 (better than persistence at 18 of 24 leads; mean −0.265 K) | held |

**What it means.**
- **The night wind excess (P-69) was mostly the diagnosis.** With a ground
  temperature and similarity, the night ratio goes from 1.63 to 0.92 and the
  day ratio stays at 1.22.
  - Against AD the wind RMSE falls in u by 0.107 m/s (18 of 24 leads) and in
    v by 0.163 m/s (22 of 24 leads, 17 of 20 cycles). That is the first
    model change that lowers the wind error.
  - Against persistence, u is still worse (+0.23 m/s) and v is level (−0.02).
- **Temperature beats persistence for the first time over a campaign**
  (3.25 against 3.56 K pooled), but not at every lead.
  - At leads 19–24 the diurnal swing is too large: night bias −2.61 K, day
    bias +1.59 K. AK is worse than AD at leads 23–24 (+0.58, +0.34 K) and
    worse than persistence there by about 1 K.
  - The stop rule applies to the night bias (colder than −1.0 K). The first
    suspect is the clear-sky assumption: every night in the model is clear,
    so it radiates too much on cloudy nights. No constant is retuned on this
    test.
- **Persistence scored by similarity is not a fair reference.** Its ground
  temperature is held fixed (day 0.61, night 0.89), so its daytime 10 m wind
  is reduced as if the cycle-time stability lasted all day. The standard
  persistence (`adp_w`) stays the reference.
- `--land-surface` stays off by default until the user decides, as with the
  heating. The `similarity` operator only works with it.


## 2026-10-02 (evening) — Threads measured, a holdout week, and the first non-hydrostatic core

**Context (prompts 164–165).**
- The thread cap may go above 26, as long as other users can still get cores.
- The user chose roadmap order B: the non-hydrostatic core first, then
  moisture in it.
- "You can run the check": a holdout test of the AM candidate.

**Threads measured (`tools/bench_cam.py`; server idle, nice 19).**
- At 12 km a 2 h forecast takes 49–60 s at 8 to 64 threads, slowest at 64.
- On 3 km, 40-level arrays the tendency call takes 1670 ms at 8 threads,
  1147 at 16, 1087 at 32 and 1364 at 64.
- The code is limited by memory bandwidth (about 160 ns per cell).
- So more threads buy 1.5 times, not the 2.5–4 times assumed in
  `docs/CAM_DESIGN.md`.
- A 3 km whole-domain hydrostatic core alone would take about 4.4 h per
  24 h. The 1.5 h budget is the binding constraint. The design document's
  section 3a records this, with float32 and domain decomposition as the
  untested remedies.

**Test AN (predictions first, `verification_tests/an/an_predictions.txt`): a holdout for the AM candidate.**
- 32 cycles the candidate never saw: 19–25 Sep, chained from a spin-up at
  18 Sep 18Z, and 1 Oct.
- Built from archived observations in a separate data root, then
  persistence, the default dry model and AM, then verification. It runs on
  the server; result below when it lands.
- The analyses build in 5–15 s each.

**Stage S1a: a 2-D non-hydrostatic core (`src/dynamics/nh2d.py`).**
- A compressible x–z core: height coordinate, flat ground, periodic in x.
- Exner-pressure form; Arakawa C grid; RK3 with acoustic sub-steps
  (forward u; vertically implicit w–pressure with off-centring 0.1); a
  pressure-based divergence damping; upwind3 horizontal advection.
- Splitting S1 into S1a (height coordinate, flat) and S1b (the Laprise mass
  coordinate, with terrain) is a deviation from the roadmap. It was made so
  that a failed benchmark points at one cause, the solver or the coordinate,
  not both.

**A lapse:** the first density-current run (100 m, dt 1 s) was made before
its predictions were written. The predictions for the other runs were
written after it and before them (`s1a_predictions.txt`, quoted here).

| run (Straka et al. 1993 density current, K = 75 m²/s, 900 s) | front (km) | min θ′ (K) | max\|u\| | max\|w\| | symmetry error (K) |
|---|---|---|---|---|---|
| reference, 25 m (Straka et al. 1993) | about 15.5 | about −9.8 | | | |
| 200 m, dt 2 s | 15.50 | −11.61 | 31.1 | 18.2 | 1.6e-13 |
| 100 m, dt 1 s (run before predictions) | 15.35 | −10.22 | 34.1 | 16.1 | 3.7e-13 |
| 100 m, dt 0.5 s | 15.35 | −10.20 | 34.1 | 16.1 | 1.5e-13 |
| 100 m, ns 12 | 15.35 | −10.20 | 34.1 | 16.1 | 3.9e-13 |
| 50 m, dt 0.5 s | 15.38 | −9.65 | 35.9 | 16.0 | 3.9e-13 |

| prediction | result | |
|---|---|---|
| 1. 200 m: front 14.0–16.0 km, min θ′ −8.0 to −10.8 K; 50 m: front 15.0–15.9 km, min θ′ −9.2 to −10.4 K | 200 m: front held, **min θ′ −11.61 refuted**; 50 m: both held | partly refuted |
| 2. dt 0.5 s: front within 0.10 km, min θ′ within 0.15 K | 0.00 km, 0.02 K | held |
| 3. ns 12: front within 0.10 km | 0.00 km | held |
| 4. symmetric to 1e-9 K | 3.9e-13 at most | held |

- The minimum θ′ converges toward the reference: −11.6 at 200 m, −10.2 at
  100 m, −9.65 at 50 m, against −9.8 at 25 m. The front is within 0.2 km of
  the reference at every resolution.
- The flow shows the three rotors of the published solution
  (`docs/nh2d_density_current.png`).
- The 200 m overshoot is the coarse resolution plus a non-monotone advection
  scheme.
- A sound pulse travels at 345.3 m/s against 345.8 m/s in theory.
- `test_nh2d.py` 4/4: the tridiagonal solve, the sound speed, the 200 m
  current, and the acoustic-step independence.

**Not yet done in S1.**
- The pressure equation is the linearised Exner form, so mass is not
  conserved exactly. The mass-conservation gate belongs to S1b, whose
  coordinate is in flux form.
- The mountain-wave gate is also S1b, which needs terrain.

**Stage S1b: the same solver in the mass coordinate, with terrain (`src/dynamics/nh2d_mass.py`).**
- The coordinate is Laprise (1992), as in WRF-ARW (Skamarock and Klemp 2008):
  NWP1's sigma made non-hydrostatic.
- Prognostic variables: column mass μ (flux form), U = μu, W = μw,
  Θ = μθ, and the geopotential φ. Pressure comes from the equation of state.
- The horizontal pressure gradient is in perturbation form against a
  hydrostatic reference, with the reference balance removed analytically.
- RK3 with acoustic sub-steps linearised about each stage state (Klemp et al.
  2007): U forward, then μ, Ω and Θ, then a vertically implicit
  W–φ tridiagonal solve.

**Predictions (`s1b_predictions.txt`)** were written after the two rest tests
and before the benchmark runs. The rest tests: an atmosphere at rest over a
1000 m bell stays at rest to 3e-12 m/s for 1 h, and mass changes by 3e-16.

| prediction | result | |
|---|---|---|
| 1. mountain wave (h 10 m, a 10 km, U 10 m/s, isothermal 250 K; 60 levels; 5 h): correlation with the analytic w ≥ 0.90 and slope 0.85–1.15 at 1–10 km | correlation 0.63, slope 0.55 | **refuted** |
| 2. mass change < 1e-12 | −1.3e-16 | held |
| 3. density current at 200 m: front within 0.2 km of S1a (15.50), min θ′ within 0.6 K of S1a (−11.61) | front 15.50; **min θ′ −18.2 K**, colder than the initial −15 K | **refuted** |

**What the refutations showed, and what was changed (recorded as a change
after a failed test, not a pass):**
- **Mountain wave.** The fit is good near the ground (1–3 km: correlation
  0.96, slope 0.94) and degrades with height. It did not improve from 5 to
  10 h, so the cause is not spin-up.
  - 60 levels give about 6 points per vertical wavelength (2π/l = 3.2 km).
    With 120 levels the fit is correlation 0.91 and slope 0.78 at 5 h, and
    0.96 and 0.90 at 10 h. So it was vertical resolution, plus upper levels
    still adjusting at 5 h.
  - With the advection fix below, 120 levels and 10 h: correlation 0.92,
    slope 0.86 (`docs/nh2d_mass_mountain_wave.png`).
- **Density current.** The θ flux was second-order centred in both
  directions, which undershoots at a sharp cold front. S1a used upwind3
  horizontally. Θ now uses a third-order upwind-biased flux in x and in the
  interior of the column.
  - At 200 m: front 15.30 km, min θ′ −10.70 K.
  - At 100 m: front 15.25 km, min θ′ −9.99 K, against about 15.5 km and
    −9.8 K in the reference and 15.35 km and −10.2 K for S1a.
  - Mass is conserved to round-off, and the field is symmetric to 1e-11.
- `test_nh2d_mass.py` 3/3, about 45 s: rest over the bell; the lower-
  troposphere wave against the analytic solution (60 levels, 5 h: correlation
  0.95, slope 0.95); the 200 m density current.

**Where S1 stands.**
- The roadmap's S1 gates (density current within the published spread,
  mountain wave against its analytic solution, mass to round-off) are met by
  S1b at 100 m and 120 levels.
- Two of the three predictions written beforehand failed. Each failure
  identified a real limitation, and both are recorded above.
- Next is S2: the 3-D version on the 12 km domain, compared with the
  hydrostatic model over the 20-cycle campaign.

---

## 2026-10-02 (night) — S2: the 3-D non-hydrostatic core in the forecast pipeline

**Context.** The user chose roadmap option B (non-hydrostatic core first,
moisture second). S1 ended with a 2-D mass-coordinate core that met its
gates. S2 puts the 3-D version behind the interface the forecast already
uses, so the same analysis, boundaries, physics and verification apply.

**What was built.**
- `src/dynamics/nh3d.py`, class `NH3D`: the S1b mass-coordinate
  (Laprise 1992) split-explicit core in 3-D on NWP1's sigma levels (index 0
  = lid) and C-grid. Prognostic μ, U, V, W, Θ, φ; acoustic substeps ns = 6;
  upwind3 horizontal θ flux (the S1b fix).
- Class `NHModel`, an adapter that looks like `PrimitiveSigma` to
  `forecast.py`: properties u, v, theta, pi, surface_pressure; `set_state`
  (W = 0 and φ integrated hydrostatically from α = RT/p); `max_dt` (60 s cap,
  0.7 advective CFL); `sigma_dot`; the PAV convective adjustment after each
  step; and the existing drag, Richardson mixing and sponge as tendencies,
  with W damped in the sponge.
- `forecast.py --core {hydrostatic,nh}`. The hydrostatic model is still built
  first. It supplies the grid, terrain and physics switches, and the
  prepared initial state is loaded into the NH core. `--core nh` is NumPy
  only and refuses `--land-surface` and `--surface-heating` for now.

**Checks.** `test_nh3d.py` 5/5:
1. Rest over a 1500 m mountain when the state is the core's own reference
   profile: |u|, |v|, |w| < 1e-9 m/s after 1 h. This holds by construction
   and only shows the code is self-consistent.
2. Inertial oscillation to 0.01 m/s.
3. The x and y directions give transposed fields to 1e-10 K.
4. The adapter's edge relaxation depends on elapsed time, not on the step
   count (see the defect below).
5. The September terrain test that the first hydrostatic core failed (it
   diverged at +3 h; `docs/instability_growth.png`). This time the state
   differs from the reference profile: isothermal 250 K at rest over a
   2500 m mountain. The spurious wind settles at 0.13–0.27 m/s and |w| at
   2–4 mm/s, with no growth over 24 h and mass conserved to 4e-16
   (`docs/nh3d_terrain_rest.png`).

`test_forecast.py` 15/15 is unchanged.

**First real case (desktop, NumPy, one core; 28 Sep 06Z, the P-67 case).**
The 24 h run completed in 18.0 min, with max|u| 13.4–15.6 m/s. The
reference is AL (hydrostatic core, upwind3, same physics).

*Defect found by the comparison (P-70, category C: a per-step constant
used where a per-time rate was meant).* `forecast.py`'s Davies weights are applied once per step. Their
values (width 15, α 0.1) were set with the hydrostatic core's ~17 s step,
so at the NH core's 60 s step the same weights relaxed ~3.5 times more
weakly per hour. The adapter now rescales them to the same e-folding time,
a_eff = 1 − (1 − a)^(dt/dt_ref), with dt_ref the hydrostatic core's stable
step for the run (16.7 s here).

| NH − AL, rms over the grid | 1 h | 24 h |
|---|---|---|
| u, per-step weights (m/s) | 0.38 | 0.73 |
| u, weights per unit time (m/s) | 0.025 | 0.21 |
| v, weights per unit time (m/s) | 0.024 | 0.19 |
| θ, weights per unit time (K) | 0.012 | 0.12 |
| surface pressure, per-step weights (hPa) | 0.84 | 0.59 |
| surface pressure, weights per unit time (hPa) | 0.078 | 0.085 |

Before the fix the u difference was 38 % of AL's own 24 h change
(1.92 m/s rms). After it, the difference is 11 %. The hydrostatic core's
24 h change in θ is 1.0 K rms. The largest local θ difference is 2.1 K,
where the two runs' convective adjustments differ.

**For the collaboration study.** This defect could only be found by
running the new core next to the old one on the same case. Every unit test
of each module passed. The weights were correct for the core they were
tuned on and silently wrong for any other step length. A comparison against
observations alone would have shown it only as slightly worse boundary
behaviour.

**Status.** S2 integration done; the 20-cycle gate (test AO) is the next
entry.

---

## 2026-10-02 (night) — Test AN: the AM candidate on a week it never saw

**Context.** AM (`--advection upwind3 --land-surface --z0-land 1.0
--z0-sea 0.0002`, scored with the similarity operator) beat persistence on
the 26–30 Sep development week by T −0.288 K, u −0.045 and v −0.069 m/s.
Every choice in it was made while looking at that week. The user approved a
holdout check before any of it becomes a default.

**Method.** New analyses were built from archived observations in a separate
data root (`scratch/holdroot`, code c7380aa, same ingest). They cover
19 Sep 00Z – 25 Sep 18Z (28 cycles, chained from an 18 Sep 18Z spin-up) and
1 Oct 00–18Z (4 cycles). Three arms: persistence, AD (the present default)
and AM. Predictions were written before any cycle was built
(`an/an_predictions.txt`).

**Result.**

| Prediction | Outcome |
|---|---|
| 1. All AD and AM runs complete 24 h | **Refuted.** AM 32/32. AD diverged twice: 22 Sep 18Z at 22.50 h and 1 Oct 12Z at 23.59 h |
| 2. AM (sim) beats persistence in T, at ≥ 12 of 24 leads | **Held.** −1.055 K, 20 of 24 leads, 27 of 30 cycles |
| 3. AM u and v no worse than persistence by > 0.05 m/s | **Held**, u narrowly: u +0.043 (13 of 24 leads), v −0.029 (18 of 24) |
| 4. AM beats AD in T, u and v | **Held.** T −1.089 K, u −0.308, v −0.081 m/s |

The comparisons use the 30 cycles every arm verified. The job verified
persistence and AM only after AD's verification succeeded, so the two cycles
where AD diverged were not scored for any arm. AM completed both of them.

Pooled RMSE, all leads (standard operators unless stated):

| Arm | 2 m T (K) | 10 m wind vector (m/s) | Speed ratio day / night |
|---|---|---|---|
| Persistence | 4.94 | 3.09 | 0.63 / 0.99 |
| AD | 4.59 | 3.36 | 0.89 / 1.37 |
| AM, standard operator | 3.99 | 3.47 | 0.75 / 1.46 |
| AM, similarity operator | 3.44 | 3.10 | 0.86 / 0.64 |

**Interpretation.**
- **AM's temperature margin is larger on the holdout than on the week it
  was developed on** (−1.06 against −0.29 K). This is not evidence that AM
  improved. Persistence did much worse on the holdout week (T RMSE 4.94
  against 3.56): its diurnal range was larger, and AM's advantage is the
  land surface's diurnal cycle. AM also beats persistence in T under the
  standard operator (3.99 against 4.94), so the result does not depend on
  the operator.
- **The wind margin is small, and the same leads lose as before.** u is
  worse than persistence at leads 1–4 and 19–24, by up to +0.28 m/s. T is
  worse at leads 21–24, by up to +0.47 K. This is the same late-lead
  weakness the development week showed. The model has no clouds and no
  moisture, so its second afternoon and night are where persistence wins
  (P-59).
- **AM's night 10 m wind is now too weak under the similarity operator
  (ratio 0.64).** In the development week it was 0.89. The stable-layer
  reduction may be too strong on clear, calm nights. This is recorded and
  not tuned.
- **The centred2 advection diverged twice more on unseen cycles, near the
  end of the run (P-67).** upwind3 completed both. These are the first
  failures of the unmodified default configuration; the two earlier P-67
  cases used z0 = 1.0 and surface heating. centred2 has now failed on four
  real cases, and upwind3 on none of the 54 it has run (20 campaign + 2
  P-67 cases + 32 holdout).

**Status.** The holdout supports the candidate. Predictions 2–4 held;
prediction 1 failed only for the old default.

**Decision (user, 2026-10-03): adopt the full AM configuration.** Defaults
now:
- `forecast.py`: `--advection upwind3`; the force-restore land surface on
  for the hydrostatic core (`--no-land-surface` turns it off; it is off
  automatically with `--surface-heating` or `--core nh`); and roughness
  1.0 m over land and 0.0002 m over water unless `--z0` is given.
- `verify.py`: `--surface-operator auto`. Forecasts that carry a ground
  temperature are scored with the similarity operator. Persistence, and
  forecasts with no ground temperature, keep the standard operator, so the
  persistence reference is unchanged (`resolve_surface_operator`; test in
  `test_verify.py`, 8/8).

The old configuration is `--advection centred2 --no-land-surface --z0 0.1`.
`daily.sh` passes none of these flags, so the server cycle picks up the new
defaults on its next `git pull`. P-67 is closed. P-59 stays open for the
late leads and the weak night wind.

---

## 2026-10-03 — Test AO: the S2 gate passes

**Method.** The non-hydrostatic core (`--core nh --advection upwind3 --z0 0.1
--no-land-surface`, i.e. AL's physics) was run on the 20 campaign cycles
(26 Sep 00Z – 30 Sep 18Z). Each run used NumPy on one server core at nice 19,
with all 20 in parallel. Runs were scored with the same operators as AL and
compared with AL (the hydrostatic core with the same physics). Predictions
were written on 2026-10-02 before any server run (`ao/ao_predictions.txt`).
Submission was delayed a night, and the job line was updated for the new
defaults; the predictions were not changed.

| Prediction | Outcome |
|---|---|
| 1. 20/20 complete, max wind < 60 m/s | **Held.** 20/20; highest max wind in a run 47.6 m/s (the 26 Sep jet) |
| 2. Mean over leads within ±0.10 K and ±0.10 m/s of AL | **Held.** T −0.015 K, u −0.006, v +0.012 m/s |
| 3. No lead differs by more than 0.25 | **Held.** Largest: T −0.03 K, u −0.03, v +0.02 m/s |
| 4. Every cycle under 40 min on one core | **Held.** 25.8–28.6 min |

**Interpretation.**
- At 12 km the two cores give the same forecast to within a few hundredths
  of a kelvin and a metre per second. That is the expected result, since
  non-hydrostatic effects are negligible at this grid spacing. It shows the
  NH core, its adapter and the physics hooks reproduce the model already
  verified. It says nothing yet about storm-scale skill.
- The NH core is slightly better in T at every lead (by up to 0.03 K). The
  difference is too small to interpret.
- **Cost is now the binding question.** On one core the NH core takes
  ~27 min for 24 h at 12 km. The hydrostatic core takes ~10 min at
  8 torch threads. At 3 km there are 16 times as many columns, and the
  step is about 4 times shorter. S3 (benchmarks: float32, threads, domain
  decomposition, a torch backend for `nh3d.py`) decides whether 3 km fits
  the 1.5 h budget, and on what domain.

**Status.** S2 passed. Next is S3.

---

## 2026-10-03 — S3a/S3b: what the 3 km core costs, and what could pay for it

**Context.** S2 passed with the NH core in NumPy on one core. Before the CAM
is designed further, S3 measures its cost at 3 km and the back ends that
could bring it inside the 1.5 h budget. Constraint: nothing may be installed
on the server. Probes found:
- torch 2.8 (CPU) works, but `torch.compile` does not: there is no C++
  compiler.
- numba is present but broken (fails to import against this NumPy).
- `gcc` 11.5 with OpenMP works.

The server has two 26-core Xeon Gold 6230R sockets (two NUMA nodes).
Predictions were written before each run (`s3/s3a_predictions.txt`,
`s3/s3b_predictions.txt`).

**S3a: the core as written** (`tools/bench_nh3d.py`, dynamics only, NumPy
float64, one core).

| Grid | Step | ns per cell-step | 24 h | Peak memory |
|---|---|---|---|---|
| 12 km × 20 (110×97), dt 60 s | 0.58 s | 2727 | 0.23 h | 0.30 GB |
| 12 km × 40, dt 60 s | 1.19 s | 2782 | 0.47 h | 0.40 GB |
| 4 km × 40 (330×291), dt 25 s | 14.2 s | 3689 | 13.6 h | 2.1 GB |
| 3 km × 20 (440×388), dt 20 s | 12.5 s | 3647 | 14.9 h | 1.9 GB |
| 3 km × 40, dt 20 s | 35.0 s | 5122 | **42.0 h** | 3.2 GB |

1. Cost per cell-step size-independent within ±30 %: **refuted.** It grows
   with the array: +35 % at 4 km and +88 % at 3 km × 40, as the working set
   leaves the cache.
2. 3 km × 40 over 20 h on one core: **held** (42 h, 28 times the budget).
3. Memory under 10 GB: **held** (3.2 GB).
4. The acoustic substep is the largest own-time entry at every size:
   **held** (about half the step). `np.roll` is next (about 15–20 %).

**S3b: one stencil, three back ends** (`tools/bench_stencil.py`). One
memory-bound stencil of the core's kind on the 3 km × 40 array. The C
version is plain C with OpenMP, compiled by the server's gcc and called
through ctypes.

| Back end | Threads | ns/cell | vs NumPy float64 |
|---|---|---|---|
| NumPy float64 (today) | 1 | 34.6 | 1× |
| NumPy float32 | 1 | 14.2 | 2.4× |
| torch float64 | 26 | 9.8 | 3.5× |
| torch float32 | 26 | 4.6 | 7.6× |
| C float64 | 1 | 3.2 | 11× |
| C float64 | 26 | 0.38 | 91× |
| C float32 | 26 | 0.16 | **216×** |
| C float32 | 52 | 0.17 | 204× |

1. C float32 at 26 threads at least 20× faster: **held** (216×).
2. torch float32 at most 4× faster: **refuted** (7.6×). The September
   thread benchmark understated torch on this stencil.
3. 26 → 52 threads gains less than 1.3×: **held**. 52 threads is slightly
   slower, as expected when one thread first touches the arrays, which then
   sit on one NUMA node.
4. All back ends agree with NumPy float64 (float32 to 4.6e-7): **held**.

**Interpretation.**
- NumPy evaluates each operation as a separate pass over memory. On a
  memory-bound stencil, most of the cost is moving arrays, not arithmetic.
  A fused C loop reads each input once. That alone is 11× on one core.
  OpenMP adds another 8–14× on one socket, until the socket's memory
  bandwidth is saturated: 0.16 ns/cell with four float32 arrays is about
  100 GB/s, close to the six DDR4 channels' peak.
- **Projection (figure `s3_backend_cost.png`).** A whole-core speed-up is
  smaller than one stencil's, because of the tridiagonal solves, the
  physics and the Python between kernels. If the core reaches 30×, the
  3 km × 40 dynamics take about 1.4 h. At 100× they take 0.4 h, which
  leaves room for moisture and physics. torch float32 at 7.6× gives
  5.5 h, which does not fit on the full domain.
- So the whole 3 km domain inside 1.5 h needs compiled kernels (C/OpenMP on
  one socket, about 26 threads). Without them the options are 4 km, a
  smaller 3 km domain, or a larger budget. That is a design decision for
  the user.

**Status.** S3a/S3b done. The back end, and with it the 3 km domain, waits
for the user's choice.

---

## 2026-10-03 — S4: the NH dynamics in C, 44 times faster, same answer

**Context.** The user chose compiled C kernels (prompt 170). Design:
- One C file, `src/dynamics/nh3d_kernels.c`, compiled by the machine's own
  gcc at first use (`src/dynamics/cnh.py`). It is cached by a hash of the
  source and flags and called through ctypes, so nothing is installed.
- `NH3D(backend="c")` and `forecast.py --nh-backend c --threads N`. NumPy
  stays the default and the reference.
- Every kernel transcribes a NumPy expression in the same order of
  operations, compiled with `-ffp-contract=off` so that rounding can match.
- Horizontal neighbours come from index tables built with the grid's own
  edge rule, so periodic and replicate edges are both exact.

There are eleven kernels in all:
- the acoustic substep in three passes (pressure and divergence damping;
  the U, V update; a column pass for μ, Ω, Θ, φ and the implicit w solve);
- the tendencies in three passes (column diagnostics; momentum; Θ, W, φ);
- the acoustic set-up in two passes;
- the final sum.

The physics hooks (drag, mixing, PAV) are still NumPy.

**Checks** (`test_nh3d_c.py` 6/6, run on the server; the desktop has no
working gcc and reports SKIPPED):
- tendencies equal NumPy to 1.8e-13 relative (dμ, FΘ and Fφ exactly);
- a full step (6 steps) equals NumPy to 4e-13 relative, with both edge
  modes and ns = 6 and 12;
- the result is bit-identical on 1 and 4 threads;
- a wrong dtype is refused.

**Speed** (predictions in `s4/s4a…s4c_predictions.txt`, written before each
run):

| Test | Configuration | Step | 24 h | Prediction |
|---|---|---|---|---|
| S3a | NumPy, 3 km × 40, 1 core | 35.0 s | 42.0 h | — |
| S4a | acoustic substeps only in C, 26 threads | 17.6 s | 21.1 h | 2–3× faster: **held** (2.0×) |
| S4b | all dynamics in C, 13 threads | 1.36 s | 1.63 h | — |
| S4b | all dynamics in C, 26 threads | 0.79 s | **0.95 h** | under 1.0 s: **held** (44× NumPy) |
| S4b | all dynamics in C, 52 threads | 0.64 s | 0.77 h | gain from 26 → 52 under 1.3×: **held** (1.23×) |
| S4b | 12 km × 20, 8 threads | 0.076 s | 1.8 min | under 0.05 s: **refuted** |

The column pass of the acoustic substep is now the largest single cost
(33 % of the step).

**Real case (S4c).** 28 Sep 06Z, 24 h, test AO's configuration, 8 threads.
The run completed in 4.7 min, against 27.8 min for AO's NumPy run of the
same cycle. Field differences from the NumPy run at 24 h are rms
3.2e-13 m/s in u, 3.4e-13 in v and 1.9e-13 K in θ, i.e. round-off. All
three predictions held. At 12 km most of the 4.7 min is now the NumPy
physics and Python. The physics is called in every RK stage.

**Interpretation.**
- The full 3 km × 40 level dynamics now fit inside the budget on one
  socket: 0.95 h for 24 h at 26 threads. That leaves about 0.5 h for
  physics, moisture, I/O and the analysis.
- Two costs come next:
  1. The NumPy physics is called three times per step and will dominate
     at 3 km. The standard remedy is physics once per step, held over the
     RK stages, with the physics in C later.
  2. Moisture adds advected fields (vapour, cloud, rain; later ice).
     Each adds about one Θ-like transport per step.
- float32 (S3b: a further ~2.4× on a stencil) is held in reserve. The C
  file already compiles a float32 set. It is not used until a test shows
  float32 does not change the forecast.

**Status.** S4 (the C core) done and verified. The core is a drop-in
backend; nothing changes for the 12 km production model.

---

## 2026-10-03 — S5: decisions, physics once per step, 40 levels, and a 3 km analysis that fills the gaps

**Decisions (user).**
- The 3 km run gets its own analysis from observations, with lateral
  boundaries from the 12 km forecast of the same cycle (prompt 172).
- Station interpolation to fill the gaps between stations should be part
  of it (prompt 173).

The plan for stage S5 (nine steps) was approved with that addition.

**S5a: physics once per step, and 40 levels** (12 km, 20 campaign cycles, C
backend, 8 threads; predictions in `s5/s5a_predictions.txt`).
- `NHModel(physics_every="step")` / `forecast.py --nh-physics step`
  evaluates drag, mixing, the sponge and the w damping once per large step
  and holds them over the three RK stages, as WRF does. `test_nh3d.py` 6:
  5 evaluations in 5 steps instead of 15.
- `forecast.py --levels N` sets the sigma level count. The sponge depth
  scales with it (10 levels at 40).

| Prediction | Outcome |
|---|---|
| 1. All 40 runs complete | **Held** (20/20 per arm) |
| 2. Physics once per step vs AO within ±0.05 | **Held**: T, u and v all 0.000 to three decimals |
| 3. A 12 km cycle under 2.8 min | **Refuted**: 3.0–4.0 min, though with four runs and test S5b sharing the machine, so not a clean timing |
| 4. 40 vs 20 levels within ±0.15 | **Held**: T +0.116 K (worse, up to +0.25 K at leads 8–10), u −0.100, v −0.091 m/s (better at all 24 leads) |

40 levels moves the lowest level from σ 0.965 (~290 m) to 0.983 (~140 m),
and that changes what the 2 m and 10 m operators extrapolate from. The
better wind is consistent with a lowest level closer to 10 m. The worse
daytime temperature has not been diagnosed. This core has no land surface
yet, so its lowest level has no diurnal cycle. The CAM uses 40 levels
because convection needs them; the temperature difference is recorded and
not tuned.

**S5b: the 3 km analysis and the gap-filling** (8 cycles, 27–28 Sep; rebuilt
from the observations the 12 km cycles archived; scored only on the withheld
ASOS stations, one in five; predictions in `s5/s5b_predictions.txt`).
- `ingest_obs.py --spacing-km 3` produces a 389 × 439 grid. It has new
  options:
  - `--out-root`: a separate tree (`tensors_3km`);
  - `--from-raw --raw-from`: the 12 km cycle's archived observations;
  - `--background-analysis`: the 12 km analysis of the same cycle,
    regridded (`build.regrid_features`), as the first guess.
- The terrain is fetched once at 1 arc-minute: `geo.load_terrain` now picks
  the ETOPO stride from the cell size.
- **Gap-filling** (`--gapfill`): `barnes.barnes_increments_aware` changes
  the surface Barnes analysis three ways. All three constants were chosen
  before the run:
  - weights fall off with the height difference between station and grid
    point, exp(−(Δz/300 m)²), so a valley station does not set a ridge;
  - weights are halved between land and water;
  - a third pass shortens the length scale (120 → 71 → 42 km).

  It is computed in chunks (a 3 km grid against ~2000 stations).
  `test_barnes_aware.py` 4/4.
- `score_withheld` now scores the 10 m wind as well as 2 m T.

Pooled over the 8 cycles (572 withheld temperatures, ~520 per wind
component):

| Analysis | 2 m T RMSE (K) | 10 m u (m/s) | 10 m v (m/s) | Build (s) |
|---|---|---|---|---|
| 12 km, as in production | 1.641 | 1.916 | 2.040 | 9 |
| 3 km, plain Barnes | 1.603 | 1.905 | 2.043 | 53 |
| **3 km, gap-filling** | **1.541** | **1.871** | **1.962** | 124 |
| 3 km, gap-filling, terrain slope limit 0.03 | 1.571 | 1.870 | 1.992 | 120 |

Predictions:
1. Builds under 10 min: **held**.
2. Gap-filling beats plain 3 km and 12 km in T: **held** (−0.10 K against
   12 km; better in 7 of 8 cycles).
3. Plain 3 km within ±0.1 K of 12 km: **held** (−0.04).
4. Gap-filling winds no worse than 12 km by more than 0.1: **held**
   (better: u −0.05, v −0.08).
5. Steeper terrain (slope limit 0.03 instead of 0.0086) at least as good:
   **refuted** (1.571 against 1.541).

**Interpretation.**
- The gain is modest and comes from where the stations are. A 3 km grid
  cannot see structure between stations 30–50 km apart. What the
  gap-filling adds is to stop increments crossing large height
  differences and coastlines, and to let dense clusters of stations
  resolve smaller scales.
- Why steeper terrain did not help is open. One candidate is that
  lapse-correcting a station to steeper grid terrain moves its value
  further with the fixed 6.5 K/km lapse rate, which is wrong under an
  inversion.
- A logistics defect showed up here (P-71). Every Python step sets its
  thread caps to half the machine (`resources.py`). Eight parallel analyses
  therefore took the load to 94 on 104 cores. `tools/cam_cycle.sh` now sets
  `NWP_RESOURCE_FRACTION=0.25`.

**Status.** Physics once per step and 40 levels are used for the 3 km run.
The gap-filling analysis is the 3 km analysis. Next is the first nested
forecast (S5c).

---

## 2026-10-03 (night) — S5c: the first 3 km nested forecast is stable, and four times too slow

**Method.**
- `tools/cam_cycle.sh 2026-09-28T06` with `NWP_SKIP_12KM=1`, so the
  campaign's existing 12 km analysis and forecast were used.
- The 3 km gap-filling analysis (93 s) was followed by the NH core in C at
  26 threads, 40 levels, physics once per step, and edges from the 12 km
  forecast every hour.
- The time step was 25.7 s, set by the new acoustic limit.

Predictions were written before the run (`s5/s5c_predictions.txt`).

| Prediction | Outcome |
|---|---|
| 1. Completes 24 h, wind under 60 m/s | **Refuted** for 24 h: the deadline stopped it at 9.54 h after 84 min. It was stable throughout (max wind 18.9 m/s). |
| 2. Analysis + forecast under 80 min | **Refuted**: 3.7 s per step, ~3.5 h projected for 24 h |
| 3. Peak memory under 30 GB | **Held** (10.2 GB) |
| 4. Edge θ within 1 K of the 12 km driver; interior differs more | **Held**: rms θ difference by band, at 6 h and 9 h: outer 5 cells 0.008 K; rest of the zone 0.35 K; interior 0.54–0.64 K |

**Where the time went** (cProfile of 1 h, 141 steps):

| Part | s per step |
|---|---|
| C dynamics kernels | ~1.1 |
| NumPy physics (mixing 0.85, drag 0.24, the rest 0.3) | 1.44 |
| Relaxation over the whole domain (incl. hydrostatic φ) | 0.47 |
| Boundary-frame time interpolation, whole domain | 0.14 |
| PAV convective adjustment | 0.14 |
| Python/NumPy overhead (array adds of the physics, properties) | ~0.4 |

The frame set-up also took 16 s per hourly frame, about 7 min for
24 h, because the convective adjustment ran over every column.

**Interpretation.** The NH dynamics are no longer the cost; everything
around them is. At 12 km these parts were negligible. At 3 km each one is a
full-domain NumPy pass. This is a defect of the integration, not of the
physics: P-72.

**S5d (the fix, in progress).**
1. Drag and Richardson mixing in C (`phys_col`), transcribed from
   `surface.py` and `turbulence.py`. The sponge and the w damping now act
   only on the levels where they are non-zero.
2. Relaxation over the zone only (`relax_mu`, `relax_col`). The time
   interpolation of the driving frames happens inside the kernel, at zone
   columns only (`BoundaryDriver.bracket`).
3. The physics tendencies are added in place by a threaded kernel.
4. Frames are stabilised in the zone only (`stabilise_frame(columns=)`;
   `test_forecast.py` 16/16).
5. OpenMP threads are bound to consecutive cores
   (`OMP_PLACES=cores OMP_PROC_BIND=close` in `cam_cycle.sh`).

---

## 2026-10-04/05 — S5c2, S5e, S5c3: the 3 km 24 h forecast now fits the budget (53.6 min)

**S5c2** (job ac9ea647; the S5d fixes: physics and zone relaxation in C).
Same case as S5c (28 Sep 06Z, 24 h, 26 threads, `NWP_SKIP_12KM=1`).
Predictions in `s5/s5c2_predictions.txt`.
- Stable (max wind 20.1 m/s). About 1.6 s per step, so the deadline
  stopped it at **21.79 h** after 84.1 min.
- Writing 21 float64 snapshots with `np.savez_compressed` (3.1 GB file,
  one core, zlib level 6) took about 7 min more. The cycle took 92.2 min,
  over the 90 min budget. Peak memory 17.9 GB.
- Identical to S5c to 4e-13 up to 9 h (the C physics transcribes the NumPy
  physics). Edge θ difference from the 12 km driver 0.008 K; interior θ
  difference 0.86 K at 21 h.
- Predictions: (1) completes 24 h — **refuted**; (2) 1.0–1.6 s/step —
  **held**, at the upper edge; (3) frame time — not measured; (4)
  round-off-identical to S5c — **held**.

**Profile 2** (job 9b984695; cProfile, 141 steps = 1 h). Per step, under
the profiler: the dynamics 1.155 s, of which the four kernels that walk one
column at a time took 0.89 s (`ac_col` 0.44, `td_tw` 0.17, `td_col` 0.15,
`as_b` 0.13). Outside the dynamics: PAV 0.135 s, the u/v properties
(division by μ, three times a step) 0.10 s, C physics 0.07 s.

**Why the column kernels were slow.** Arrays are (nz, ny, nx). A kernel
that walks one column at a time reads each level ny·nx elements apart. With
~30 arrays and 40 levels that is ~1200 separate memory streams per column,
more than the hardware prefetcher follows. The point-wise kernels
(`ac_uv`, `td_uv`), which run along rows, cost a fifth as much per array.

**S5e: changes** (all bit-identical to what they replace).
1. **Blocked column kernels** (`ac_col_b`, `td_col_b`, `td_tw_b`, `as_b_b`
   in `nh3d_kernels.c`). The same per-column arithmetic in the same order,
   on blocks of 32 neighbouring columns. The level loop is outside and the
   column loop inside, so each level is a contiguous run. `cnh.py` selects
   them by default (`NWP_NH_LAYOUT=block`; `column` gives the S4 kernels).
2. **PAV adjustment in C** (`pav_col`). It works in place, with the same
   θ/u/v round trip through μ as the NumPy path.
3. **Largest wind for `max_dt` in C** (`uv_maxabs`).
4. **3 km output as float32, zlib level 1.** Grids above 100 000 columns
   (`forecast.LARGE_GRID`) only; 12 km output is unchanged. New options
   `--output-dtype` and `--output-level`. float32 holds θ to 1.5e-5 K and
   π to 0.004 Pa, far below what ASOS verification resolves.

Tests: `test_nh3d_c.py` 13/13 on the server. New tests 10–13: blocked =
column kernels bit for bit (nx = 45, a partial block, both edge modes); C
PAV = NumPy PAV exactly (1327 adjusted columns in both); C max wind = NumPy
exactly. `test_forecast.py` 17/17 (new: output precision rule and
read-back).

**S5e test** (job e4e83f10; predictions in `s5/s5e_predictions.txt`).

| Prediction | Outcome |
|---|---|
| 1. Tests pass; blocked bit-identical to column | **Held** (13/13; differences 0.0) |
| 2. Dynamics bench, block ≤ 0.75 × column | **Held**: 1.148 → 0.736 s/step (0.64×) |
| 3. Real 1 h nested run ≤ 1.15 s/step | **Held**: 0.98 s/step |
| 4. Writing the 1 h output ≤ 2.5 s | **Refuted**, narrowly: 3 s (67 MB), against 10.1 s before |
| 5. Fields at 1 h vs S5c2: max \|Δθ\| ≤ 1e-4 K | **Held**: 1.5e-5 K, the float32 rounding of 354 K |

**S5c3** (job e52a8289; the full cycle through `tools/cam_cycle.sh`,
`NWP_SKIP_12KM=1`, 26 threads at nice 19; predictions in
`s5/s5c3_predictions.txt`).

| Prediction | Outcome |
|---|---|
| 1. Completes 24 h, wind under 60 m/s | **Held**: completed, max wind 19.3 m/s |
| 2. Mean ≤ 1.05 s/step | **Held**: 0.845 s/step (3384 steps in 47.6 min) |
| 3. Cycle ≤ 65 min (≤ 70 with the 12 km run) | **Held**: **53.6 min** (analysis 1.6, forecast set-up ~1, steps 47.6, write 1.3) |
| 4. Peak memory ≤ 16 GB | **Refuted**: 18.2 GB. The snapshots are still kept as float64 in memory and stacked at the end. |
| 5. Bit-identical to S5c2 | **Held**: max \|Δθ\| 1.53e-5 K at every common hour (float32 output rounding only) |
| 6. Thread binding: spread within 10 % of close at 26 threads; 52 threads gain < 15 % | **Held**: spread 0.819 vs close 0.748 s/step (9.5 % slower); unbound 0.730; 52 threads 0.641 (14 % faster) |

- The output was 1613 MB, written in 79 s. `/data5` has 26 TB free.
- The ledger line in `compute_hours.csv`: 53.6 min, 26 threads,
  23.23 core-hours, status 0.
- With the 12 km cycle (3–4 min in S5a, about 5 with its ingest) the pair
  takes about 59 min of the 90 min budget.

**First 3 km verification** (one cycle; 2 m T against ~1040 ASOS reports
an hour). RMSE 1.66 K at 1 h, rising to 3.0–3.1 K at 12–14 h (18–20Z),
2.6–2.7 K at 18–24 h. The bias is near zero for 3 h, then cold, reaching
−1.75 K at 14 h (20Z). A cold afternoon bias is what a core with no land
surface should show (P-59: the NH core has none yet). One cycle shows
nothing about 3 km against 12 km; that is the campaign step.

**Interpretation.**
- The remaining cost is now the dynamics itself: the full nested step
  costs 0.845 s against 0.748 s for the dynamics alone. The physics, PAV,
  relaxation and loop overhead together take ~0.1 s.
- The column kernels were limited by memory access, not by arithmetic.
  Reordering the loops gave 1.56× on the dynamics with no change to a
  single bit of the answer.
- The 1 h test ran at 0.98 s/step and the 24 h run at 0.845. The server
  load was ~24 during both, so this is not contention alone. The likely
  cause is the start-up of the run (compilation cache, first-touch pages),
  amortised over 3384 steps; it was not measured.
- 52 threads would give 14 % more speed for twice the cores. 26 threads
  remain the default, which leaves the second socket to other users.

**Status.** P-72 closed. The plan steps "Cycle script for the pair" and
"Stability and cost trial" are complete. Next is the campaign verification
(3 km against 12 km over the campaign cycles).

---

## 2026-10-05 — S5g: the dry 3 km nest over 20 cycles — better 10 m wind, no better 2 m temperature

**Design** (predictions in `s5/s5g_predictions.txt`, written before the run).
- **Cycles.** The 20 campaign cycles, 26 Sep 00Z – 30 Sep 18Z.
- **How the 3 km runs were made.** `tools/cam_cycle.sh` with
  `NWP_SKIP_12KM=1` produced each 3 km run from:
  - its own gap-filling 3 km analysis;
  - edges from the 12 km forecast in `tensors_3d`, which is test AD's
    chain (hydrostatic, no land surface).
- **Code and machine.** S5e code. 19 new runs; 28 Sep 06Z was reused from
  S5c3, same code. The runs went two at a time, each pinned to one socket
  (taskset even/odd CPUs, 26 threads each), at nice 19.
- **Arms.** All four were verified again here with the same `verify.py`
  (surface operator `auto`) on the same observation archives. They are
  compared only on matches present in all four (478 774 pairs), with
  `tools/cam_compare.py`.

| Arm | What it is |
|---|---|
| cam3 | 3 km NH, 40 levels, dry, no land surface |
| nh12 | 12 km NH, 40 levels, dry (test S5a b2) |
| ad12 | 12 km hydrostatic, 20 levels: the 3 km run's driver (test AD) |
| am12 | 12 km production defaults (test AM): upwind3, land surface, z0 by land/sea, similarity operator |

**Cost.** All 19 runs completed 24 h. The largest wind in the hourly prints
was 18.5–38.3 m/s per run (`cam_<cycle>.log`). Each run's
time steps took 53.3–60.3 min, and each cycle 58.9–66.8 min, with two
cycles sharing the machine. The job took 10.1 h, at most 525 core-hours.

**Result, pooled over 20 cycles and 24 leads** (RMSE; bias in brackets).

| | cam3 | nh12 | ad12 | am12 |
|---|---|---|---|---|
| 2 m T (K), n 165 968 | 3.888 (−0.12) | 3.879 (−0.25) | 3.761 (−0.13) | **3.223** (−0.23) |
| 2 m T, day | 4.197 (−1.48) | 4.241 (−1.61) | 4.033 (−1.51) | 3.335 (+0.22) |
| 2 m T, night | 3.562 (+1.20) | 3.492 (+1.07) | 3.476 (+1.21) | 3.111 (−0.66) |
| 10 m u (m/s), n 156 349 | **2.420** (−0.83) | 2.722 (−1.14) | 2.819 (−1.29) | 2.433 (−0.93) |
| 10 m v (m/s), n 156 457 | 2.202 (−0.47) | 2.279 (−0.69) | 2.359 (−0.88) | **2.141** (−0.20) |

Per cycle, cam3 against nh12:
- T: cam3 better in 10 of 20 cycles (mean +0.005 K);
- u: better in 17 of 20 (−0.262 m/s);
- v: better in 15 of 20 (−0.073 m/s).

By lead (`cam_skill_by_lead.png`), the u gain is 0.08 m/s at 1 h and
0.3–0.4 m/s from 6 h on. So most of it comes from the forecast, not the
analysis.

**By terrain** (relief in a 21 km box of the unlimited ETOPO terrain: 225
flat stations, 119 hilly, 23 mountain), cam3 minus nh12 RMSE:

| | flat | hilly | mountain |
|---|---|---|---|
| 2 m T (K) | +0.041 | −0.022 | −0.133 |
| 10 m u (m/s) | −0.219 | −0.458 | −0.445 |
| 10 m v (m/s) | −0.006 | −0.154 | −0.422 |

**Predictions.**

| Prediction | Outcome |
|---|---|
| 1. All 19 runs complete, wind < 60 m/s | **Held** (largest hourly max 38.3 m/s) |
| 2. Every cycle ≤ 60 min with two pinned streams | **Refuted**: 58.9–66.8 min. Two runs at once cost each one 10–20 % against S5c3's 53.6 min alone. Every cycle is still inside the 90 min budget. |
| 3. T within ±0.15 K of nh12 | **Held** (+0.009 K) |
| 4. T better than nh12 at leads 1–3 h | **Refuted** as stated: better at 1 h (1.748 vs 1.804) and 2 h (2.311 vs 2.326), worse at 3 h (2.863 vs 2.853) |
| 5. Wind better than nh12 by 0–0.2 m/s in u and v | **Refuted for u** (−0.30 m/s, more than predicted); **held for v** (−0.08) |
| 6. T gain larger in mountain than flat terrain | **Held** (mountain −0.13 K, flat +0.04 K), with the class change noted below |
| 7. am12 beats cam3 by > 0.5 K by day, less by night | **Held** (day 0.86 K, night 0.45 K) |
| 8. cam3 beats ad12 in T by ≥ 0.1 K | **Refuted**: cam3 is 0.13 K worse than its own driver |

**Correction, 2026-10-07 (test S6c).** The 10 m wind gain below is a roughness effect, not resolution: the 12 km NH arm
(S5a b2) used z0 0.1 m everywhere, the 3 km run the land/sea map. With matched roughness the two are level in u.
See the S6c entry.

**A change to the instrument after a first look (stated for the record).**
The terrain classes were first computed from the forecast file's terrain.
An interim run on the first 8 cycles put 365 of 367 stations in "flat"
and none in "mountain". The cause is the 3 km terrain's slope limit
(0.0086), which caps the relief in a 21 km box at ~155 m. The classes now
use the unlimited ETOPO terrain on the same grid
(`data/static/terrain_etopo_389x439.npz`, 0–1533 m), with the class edges
unchanged. The interim had already shown the two hilly stations' scores,
so prediction 6 was not judged blind to every number. The by-lead and
pooled results do not depend on the classes.

**Interpretation.**
- **10 m wind: resolution helps, and it helps most over terrain.** The
  dry 3 km run's u error is 11 % below the dry 12 km NH run's, and its
  bias is 0.3 m/s smaller. Over hilly and mountain stations the gain is
  twice the flat-station gain. Without a land surface, the dry 3 km run
  matches production in u (2.420 vs 2.433) and comes within 0.06 m/s in
  v. This is the first CAM result that is a gain rather than a cost.
- **2 m temperature: resolution does nothing until the surface exists.**
  The two dry NH runs are the same to 0.01 K. Both have a day cold bias
  (−1.5 K) and a night warm bias (+1.2 K), the signature of a missing
  diurnal cycle (P-59: the NH core has no land surface). Production's
  land surface is worth 0.67 K over cam3. The dry hydrostatic driver
  being 0.13 K better than either NH run was not expected, and it was
  not diagnosed here. The AO gate measured NH against AL, not AD. A
  candidate is the 40-level daytime T penalty already seen in S5a
  (+0.12 K).
- **What it does not show.** One week, one season, all dry. There is no
  convection to resolve, so this tests terrain and boundary-layer flow,
  not storms. The 3 km analysis differs from the 12 km one (gap-filling),
  so "resolution" here means grid plus analysis. The lead-1 gap shows the
  analysis's share is small for wind.

**After the campaign: snapshots kept at the output precision** (check S5h,
job 738dace6). `run_forecast(snapshot_dtype=)` now holds the 3 km snapshots
as float32, which halves the ~3.9 GB they took in S5c3 (prediction 4 of
S5c3 failed on memory). A 2 h nested run with the change is identical to
S5c3's hours 1–2 (max |Δ| 0.0 in θ, u, v, π; float32 against float32).
Peak memory was 5.2 GB at 2 h. The 24 h peak was not measured.

**Status.** The S5 plan's campaign step is done. The figure is
`docs/cam_skill_by_lead.png`. The next physics for the NH core is clear:
the land surface (force-restore, as in AM) before moisture. It is what
the 2 m temperature lacks, and the CAM's T cannot be judged without it.
That is a change to the roadmap order. **The user chose it (2026-10-05): land surface for the NH core first, then moisture.**

---

## 2026-10-06 — S6a/S6b: the land surface in the NH core; at 12 km it matches production

**Decision (user, 2026-10-05, prompt 175).** The land surface goes into the
NH core before moisture.

**Change.**
- `NHModel._step_land_surface` steps the existing force-restore ground
  model (`land_surface.ForceRestoreSurface`) once per large step, after
  the dynamics and before PAV. This is where `PrimitiveSigma` steps it.
  - The inputs are the hydrostatic model's: lowest-level winds, θ₁, p_s,
    p₁, dp₁ and z₁ (the hypsometric height the drag uses).
  - It adds μ·rate to Th₁, and θ_g reaches the drag with Louis (1979)
    stability.
- `phys_col` in C takes the surface θ and the Louis momentum function, so
  the 3 km physics stays compiled.
- The land surface is on by default for `--core nh`.

**Tests.**
- `test_nh3d.py` 9/9. New:
  - one land step adds exactly ForceRestoreSurface's tendency, to the
    lowest layer only (bit for bit);
  - one hour of sun warms the lowest layer over land by +0.33 K, one hour
    of night cools it by −0.01 K, and water is held;
  - five full steps stay finite.
- `test_nh3d_c.py` 15/15 on the server: C Louis drag = NumPy to 9e-16,
  with 713 stable and 727 unstable columns.
- `test_forecast.py` 17/17.

Predictions were written before the run (`s6/s6ab_predictions.txt`). They
also stated the gate for the 3 km campaign.

**S6a (1 h smoke runs).**
- 3 km, 28 Sep 06Z: land 67.7 % (12 km: 68.1 %); z0 1 m over land and
  0.0002 m over water; ground temperature 276.9–297.1 K at 1 h; finite.
- Cost: 0.98 s/step over the first hour. The dry run's first hour (test
  S5e) was also 0.98; its 24 h mean was 0.845.

**S6b (12 km gate).**
- Configuration: 20 campaign cycles; NH, 40 levels, C, 8 threads, physics
  once per step, upwind3; four pinned streams.
- Two arms:
  - L1: land surface, z0 0.1 m everywhere (one change from the dry S5a b2
    run, "nh12");
  - L2: land surface with the land/sea z0 map (the default, as in AM).
- 40/40 runs completed, at 2.5–2.8 min each.
- Paired matches: 478 786, the same observations as test S5g.

| RMSE (bias) | L2 (NH + land, z0 map) | L1 (NH + land, z0 0.1) | nh12 (NH dry) | am12 (production) |
|---|---|---|---|---|
| 2 m T (K) | 3.328 (−0.43) | 3.332 (−0.54) | 3.879 (−0.25) | **3.223** (−0.23) |
| 2 m T, day | 3.500 (−0.05) | 3.464 (+0.02) | 4.241 (−1.61) | 3.335 (+0.22) |
| 2 m T, night | 3.152 (−0.79) | 3.198 (−1.08) | 3.492 (+1.07) | 3.111 (−0.66) |
| 10 m u (m/s) | **2.383** (−0.77) | 2.685 (−1.08) | 2.722 (−1.14) | 2.433 (−0.93) |
| 10 m v (m/s) | **2.127** (−0.17) | 2.183 (−0.38) | 2.279 (−0.69) | 2.141 (−0.20) |

Per cycle, L2 beats:
- the dry NH run in T in 18 of 20 cycles (−0.54 K on average), in u in 18
  and in v in 16;
- AM in T in 10 of 20 (+0.08 K on average), in u in 9 (−0.03 m/s), in v
  in 10 (−0.01 m/s).

**Predictions.**

| Prediction | Outcome |
|---|---|
| 1. Tests 15/15, 9/9, 17/17 | **Held** |
| 2. 3 km 1 h ≤ 0.95 s/step | **Refuted**, narrowly: 0.98 s/step. The first hour costs the same with and without the land surface (S5e: 0.98). The 24 h mean is measured in S6c. |
| 3. 3 km land fraction within 2 points of 12 km | **Held** (67.7 % vs 68.1 %) |
| 4. Finite; ground 270–320 K | **Held** (276.9–297.1 K) |
| 5. All 40 runs complete | **Held** |
| 6. L2 T within ±0.15 K of AM | **Held** (+0.105 K) |
| 7. L2 day bias within ±0.5 K, night within ±1.0 K | **Held** (−0.05 K, −0.79 K; dry NH −1.61, +1.07) |
| 8. L1 and L2 beat dry NH in T by ≥ 0.4 K | **Held** (−0.55 K both) |
| 9. L2 wind no worse than AM by > 0.1 m/s | **Held**: L2 is better (u −0.05, v −0.01) |

**Gate for S6c: PASS.** The tests passed, all 20 L2 runs completed, and L2's
T (3.328) beat the dry NH run's (3.879). The 3 km campaign started on its
own.

**Interpretation.**
- The surface physics carries across cores. With the same land surface,
  the NH core is within 0.1 K of the hydrostatic production model in
  2 m T. It is slightly better in 10 m wind, with the same day/night
  pattern: the day bias goes from −1.6 K to −0.05 K, and the night bias
  from +1.1 K to −0.8 K.
- The roughness map matters for wind, not temperature. L1 and L2 have the
  same T (3.332 vs 3.328). With uniform z0 0.1 the u error is 0.30 m/s
  larger, and L2 beats L1 in u in all 20 cycles. With a ground that cools
  at night, z0 = 0.1 over land is too smooth: the night u bias is −0.70
  (L1) against −0.44 (L2).
- The night is now too cold (−0.8 K), as it is in AM (−0.66 K). Both
  cores share the surface scheme, so this cold night bias belongs to the
  scheme, not to the NH core. It was not diagnosed here.

---

## 2026-10-07 — S6c: with the land surface, the 3 km CAM has the best 2 m temperature of any NWP1 configuration; its wind gain was roughness, not resolution

**Design** (predictions in `s6/s6c_predictions.txt`, written before any S6
result).
- The campaign started on its own when the S6b gate passed.
- Cycles: the 20 campaign cycles.
- Each run used `tools/cam_cycle.sh` with `NWP_SKIP_12KM=1` and
  `NWP_CAM_BOUNDARY` set to the test AM 12 km forecast of the cycle. That
  is the production edge source. The land surface and the z0 map were on
  (defaults).
- Two cycles ran at a time, one per socket, 26 threads each, nice 19.
- Arms, compared on 478 786 paired matches:
  - cam3L: this run;
  - cam3: S5g, dry, AD edges;
  - nh12L2: S6b, 12 km NH with the land surface;
  - am12: production.

**Cost.**
- All 20 runs completed 24 h. Each run's steps took 55.4–56.1 min, and
  each cycle 61.1–63.3 min, with two sharing the machine.
- Peak memory was 15.1 GB, against 18.2 GB before the float32 snapshots
  (S5h).
- The job ran 11.5 h including the wait for the gate.

**Result** (RMSE; bias in brackets).

| | cam3L | nh12L2 | am12 | cam3 (dry) |
|---|---|---|---|---|
| 2 m T (K) | **2.971** (−0.08) | 3.328 (−0.43) | 3.223 (−0.23) | 3.888 (−0.12) |
| 2 m T, day | **3.116** (+0.30) | 3.500 (−0.05) | 3.335 (+0.22) | 4.197 (−1.48) |
| 2 m T, night | **2.822** (−0.45) | 3.152 (−0.79) | 3.111 (−0.66) | 3.562 (+1.20) |
| 10 m u (m/s) | 2.385 (−0.77) | **2.383** (−0.77) | 2.433 (−0.93) | 2.420 (−0.83) |
| 10 m v (m/s) | 2.181 (−0.19) | **2.127** (−0.17) | 2.141 (−0.20) | 2.202 (−0.47) |

Per cycle, cam3L against each arm:

| | cam3 (dry) | nh12L2 | am12 |
|---|---|---|---|
| T | better in 20 of 20 (−0.87 K) | better in 16 of 20 (−0.33 K) | better in 17 of 20 (−0.25 K) |
| u | — | better in 7 of 20 (+0.01 m/s) | better in 9 of 20 (−0.02 m/s) |
| v | — | better in 6 of 20 (+0.05 m/s) | better in 9 of 20 (+0.04 m/s) |

The temperature gain over nh12L2 grows with lead: −0.07 K at 1 h, −0.30 K
at 6 h, −0.48 K at 12 h, −0.46 K at 18 h and −0.24 K at 24 h. The 3 km
analysis therefore accounts for little of it; it comes from the forecast.
The gain is about the same in every terrain class: flat −0.38 K, hilly
−0.32 K, mountain −0.35 K. Figure: `docs/cam_land_skill.png`.

**Predictions.**

| Prediction | Outcome |
|---|---|
| 1. All 20 complete; each cycle ≤ 70 min | **Held** (61.1–63.3 min) |
| 2. T ≤ nh12L2 + 0.05 K, and ≥ 0.4 K better than dry cam3 | **Held** (−0.36 K; −0.92 K) |
| 3. Wind better than nh12L2 by ≥ 0.15 (u) and ≥ 0.03 (v) | **Refuted**: u +0.002, v +0.054 m/s (worse) |
| 4. u better than am12 | **Held**, narrowly (2.385 vs 2.433; 9 of 20 cycles) |
| 5. T gain over nh12L2 larger in mountain than flat terrain | **Refuted**, narrowly (mountain −0.354, flat −0.379 K) |
| 6. T within ±0.2 K of am12 | **Refuted** on the good side: 0.25 K better |

**Correction to S5g.** S5g reported a 3 km gain in 10 m wind of −0.30 m/s
in u, about twice as large over terrain. That was not resolution:
- The 12 km NH arm there (S5a b2) ran with `--z0 0.1` everywhere. The dry
  3 km run used the default land/sea map (1 m over land).
- S6b measures the map alone: L1 (z0 0.1) to L2 (map), under the same
  land surface, is −0.30 m/s in u, the same size.
- With matched roughness (S6c), cam3L and nh12L2 have the same u (2.385
  vs 2.383), in every terrain class (mountain 2.210 vs 2.195).

So the S5g wind result is withdrawn as a resolution effect. The cause was
a confounded comparison (category E, experiment design): the arm's
configuration was taken from its test name, and its command line was
never checked. The S5g entry carries a note pointing here.

**Interpretation.**
- **2 m T: the 3 km CAM is now the best configuration NWP1 has.** It is
  0.25 K better than production over 20 cycles, in 17 of them. It is
  better by day and by night: the day bias is +0.30 K and the night cold
  bias is −0.45 K, against −0.66 K for AM. The gain over 12 km with
  identical physics is 0.36 K and builds through the day.
- What is not separated: cam3L differs from nh12L2 in grid, analysis
  (3 km gap-filling) and edges (AM's forecast against nh12L2's own). The
  lead structure argues against the analysis. The edges are not tested.
- **10 m wind: no gain from resolution at this stage.** The 3 km and 12 km
  NH runs are level in u, and the 3 km run is 0.05 m/s worse in v. The
  wind skill gain seen so far in NWP1 came from the surface (the land
  surface and the roughness map), not from the grid.
- The CAM's T is now worth judging. That was the reason for putting the
  land surface first.

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

# Problem register

Every problem this project has hit, what it turned out to be, and — for the
closed ones — what was actually done about it. Open problems are first,
because they are the ones that need doing.

**Why this is separate from `RESEARCH_LOG.md`.** The log is chronological: it
records what was tried on a given day, including the attempts that went
nowhere, and it is the right record for the AI-collaboration study. But a
chronological log answers "what happened" and not "what is wrong right now",
and a problem that was diagnosed across four sessions is scattered across four
entries. This file is the by-problem view, one entry per problem, updated in
place when its status changes.

**Every fix entry names the measurement that confirmed it.** "Fixed" without a
number is an assertion. The rule this project runs on is that a fix has to
predict the outcome it was proposed to explain, so each closed entry carries
the before and after.

## Status vocabulary

| status | meaning |
|---|---|
| **OPEN** | reproducible, not fixed |
| **FIXED** | fixed, with the measurement that confirms it |
| **ELIMINATED** | investigated and ruled out as a cause; not a defect |
| **REVERTED** | tried, made things worse, removed |
| **ACCEPTED** | real, understood, deliberately not fixed — with the reason |

## Categories

Reusing the taxonomy in `AI_COLLABORATION.md`: **A** external-interface,
**B** discrete-vs-continuous mathematics, **C** stability / dimensional,
**D** array and language semantics, **E** test design, **F** wrong causal
hypothesis, **G** missing physics, **H** performance, **I** logistics.

---

# OPEN

## P-01 — Tall terrain fails above Nh/U ≈ 1
**Category** G · **Status** ACCEPTED · **Scoped** 2026-09-04

**Symptom.** 4000 m terrain reaches 6/12 forecast hours, with or without
convective adjustment.

**What is known.** The nondimensional mountain height orders every terrain
result:

| terrain | Nh/U | outcome |
|---|---|---|
| 1000 m | 0.38 | 12/12, linear wave |
| 2500 m | 0.96 | 12/12 with convection, 11/12 without |
| 4000 m | 1.19 | 6/12 regardless |

Nh/U ≈ 1 is the classical boundary between a mountain wave that propagates
over the obstacle and one where low-level flow blocks and the wave breaks. The
model reproduces the boundary without being told about it.

**Why ACCEPTED rather than OPEN.** The agreed target is 2 km of terrain, and
the model does **2500 m at 12/12**. The highest point in the Northeast domain
is Mount Washington at 1917 m; on a 12 km grid the cell mean is under 1500 m,
Nh/U ≈ 0.5. 4000 m is a mountain the domain does not contain, and the entry is
kept as a statement of where the physics runs out rather than as work to do.

**Ruled out by measurement** (P-30 to P-41): the sigma coordinate, the
timestep, sponge depth, lid height, the initialization filter, the presence of
convective adjustment, the eddy-diffusivity ceiling.

**Not the answer.** Orographic gravity-wave drag. It parameterizes *subgrid*
orography, and this mountain is 250 km wide on a 12 km grid — resolved by a
factor of twenty. Adding it would double-count the wave the model is already
simulating. Recorded because it is what a literature search suggests first.

**Reopen if** the domain is ever extended west into the Rockies, or the grid
is refined enough that a real ridge reaches Nh/U > 1.

---

## P-02 — The sponge base reflects
**Category** B · **First seen** 2026-09-03 · **Status** OPEN

**Symptom.** Growth over terrain peaks at exactly the level of the sponge's
lower edge, and moves when the edge moves:

| sponge levels | peak growth level (0 = lid) | max\|du\| at 6 h |
|---|---|---|
| 0 | 0 | 60.8 m/s |
| 5 | 5 | 36.4 m/s |
| 8 | 8 | 21.3 m/s |
| 12 | 18 (the surface) | 15.1 m/s |

**Diagnosis.** Partial reflection off the absorbing layer's lower edge. The
amplitude ramp is smooth (raised cosine) but is compressed into ~3 km, and a
mountain wave with a vertical wavelength of several km sees that as abrupt.

**Why it is still open rather than fixed.** Deepening the sponge halves the
amplitude but does **not** change survival — 11/12 at five levels, 11/12 at
eight. It is a real contamination of the upper levels that is not what ends a
run, and the obvious fixes each cost something: a deeper sponge eats the free
troposphere, a higher lid measures neutral-to-worse (P-14).

**This one nearly became a wrong answer.** It was the stated hypothesis for
the 2500 m failure and it is visibly, measurably real. Requiring the fix to
move the survival count is what exposed that it was the wrong cause.

---

## P-03 — 2500 m grows the wind after hour 15
**Category** unknown · **First seen** 2026-09-03 · **Status** OPEN

**Symptom.** With convective adjustment the 2500 m run reaches 16 hours
instead of 11. From hour 15 the wind, which had been pinned at 41–42 m/s for
fourteen hours, starts climbing (47.7 at hour 15, 48.8 at hour 16).

**What is known.** The mode is different from the one convection fixed: Ri
stays near 0.003 rather than going negative, so this is not overturning. Not
yet probed.

---

## P-06 — The observation fetchers have never touched the network
**Category** A · **First seen** 2026-09-01 · **Status** OPEN

**Symptom.** None yet, and that is the concern — the code has never been given
the chance to fail.

`src/verification/fetchers.py` passes 9/9 offline against saved fixtures.
Every single interface defect in this project's history (P-20 to P-24) was
invisible to an offline suite and appeared on first contact with the real
service. There is no reason to expect these to be different.

Two additions are in the same position and are covered by this entry: the
surface-field GRIB search from P-05 (`:(?:HGT|PRES):surface:` and its alias
table), and the whole ASOS path in `src/verify.py` — `fetch_asos` against nine
state networks has been exercised only against a saved payload.

---

## P-07 — The forecast–observation archive has no data in it yet
**Category** I · **First seen** 2026-09-01 · **Status** OPEN ·
**Machinery built** 2026-09-04

**Symptom.** Zero verification pairs on disk. No forecast produced so far can
be scored against what actually happened.

**Still the only time-sensitive item in the register.** Observations remain
downloadable from IEM for years, but the forecast that was valid for them was
never made. Missing a day costs a day of evidence permanently, and no amount
of later effort recovers it.

**What now exists**, all tested offline:

  * `src/verification/sigma_operator.py` — the observation operator for a
    sigma forecast (7/7). The base-class operator returned `field3d[0]` for a
    surface observation, and index 0 is the model LID: a 2 m thermometer would
    have been scored against the 200 hPa field, a **74 K** error that looks
    like a plausible number. Column pressures now follow the terrain, and the
    elevation correction and its size are recorded with every pair.
  * `src/verify.py` — fetch, QC, match, archive (7/7).
  * `tools/daily.sh` — one day of the archive, safe to run from cron: lock
    file, dated log, first-failure exit code, and verification attempted even
    when the forecast step failed, because a forecast that diverged at hour 8
    still produced eight hours worth archiving.

**The design decision that matters.** Raw observations are written verbatim
and compressed **before** any parsing, QC or matching is attempted, with the
forecast copied beside them. Matched pairs are derived data: if the
observation operator changes — and it will, since the elevation correction is
a standard lapse rate that is wrong on exactly the calm clear nights when it
is largest — every match can be recomputed. A failure in parsing or matching
must never cost the raw observations.

**What remains.** A run on the server. The machinery has never met the live
service, which is P-06 and is where this project's defects have always been.

---

# FIXED

## P-47 — A running job was indistinguishable from a frozen one
**Category** C, I · **Status** FIXED · **Fixed** 2026-09-04

**Symptom.** The first real archive run appeared to freeze: no output, no
progress, no way to tell whether anything was happening.

**Diagnosis — three separate causes, all mine.**

1. **Block buffering.** `tools/daily.sh` redirects to a log, and Python
   block-buffers stdout whenever it is not a terminal. The log stayed empty
   for many minutes regardless of what the job was doing.
2. **Progress printed only on the forecast hour.** At dt ≈ 15 s a forecast
   hour is ~240 steps and several minutes of wall clock, so even unbuffered
   there was nothing to see between hours.
3. **A silent network call.** IEM assembles a nine-network, thirteen-hour
   ASOS query on demand and can take minutes. `verify.py` printed nothing
   before or during it.

**Fix.** `python -u` in every step of `daily.sh`; step-level progress in
`forecast.py` every ~0.5% of the run with rate and **ETA**, which is what
turns "it is stuck" into "it has 40 minutes left"; the request URL, the
timeout and the transfer size printed around the ASOS fetch.

Also added, since nothing can be installed on that server and py-spy is
therefore not an option: `faulthandler.register(SIGUSR1)` in both
`forecast.py` and `verify.py`. `kill -USR1 <pid>` prints a traceback of every
thread to stderr and the process carries on — the difference between a slow
step and a genuine hang, with stdlib only.

**Confirmed by.** SIGUSR1 dumps a live traceback and the process continues;
`daily.sh` greps clean for `python -u` on all three steps; suites still green
(`test_verify.py` 7/7, `test_forecast.py` 11/11, `test_fetchers.py` 9/9).

---

## P-48 — The archive would have thrown away the raw observations
**Category** A · **Status** FIXED · **Fixed** 2026-09-04

**Symptom.** None observed — found while investigating P-47, before the first
successful live run.

**Diagnosis.** `verify.py` was written around storing the raw payload verbatim
before anything else, because that is the only irreplaceable part of the
archive (P-07). But it called `fetchers.fetch_asos`, which **parses
internally and returns `Observation` objects**. The "raw" text handed to
`store_raw` would have been a list of objects, and the whole point of the
design was lost. It would have raised on the first live fetch rather than
corrupting anything, but every offline test passed because they all inject a
saved payload and never call the fetcher.

**Fix.** `fetchers.fetch_asos_text` returns the payload verbatim;
`fetch_asos` is now a thin parse over it, so there is one request path and the
archive stores what the service actually sent.

**Confirmed by.** `test_fetchers.py` 9/9 with the split; `test_verify.py` 7/7
including the byte-for-byte round trip.

**Worth noting for P-06.** This is a fetch-path defect that a full offline
suite could not see, found only by reading the call rather than running it.
That is now three of this class in one week.

---


## P-04 — `forecast.py` ran the pressure-coordinate core
**Category** I · **Status** FIXED · **Fixed** 2026-09-04

**Symptom.** A real forecast used the core that P-14 replaced, so it diverged
in 2–3 hours while the sigma core reached 12/12 in the same conditions. The
only core reachable from real data was the broken one; everything measured
since the coordinate change could only be run on idealised states.

**Diagnosis.** Nothing converted isobaric HRRR data on to sigma levels, so
there was no way to build a `PrimitiveSigma` from an analysis.

**Fix.** `src/dynamics/interpolate.py`. Three steps, in order: terrain height
→ surface pressure, by finding the pressure at which the analysis geopotential
height equals the terrain (an interpolation, not a hydrostatic guess, so it
inherits the analysis's own stratification); surface pressure → the target
pressure of each sigma level; analysis columns → those pressures, interpolated
in **log(p)**, since a field is far more nearly linear in log(p) than in p and
the level gaps here run from 25 hPa near the ground to 50 hPa aloft.

Extrapolation was the part that needed care. Theta below the lowest analysis
level follows the **lapse rate of the lowest two levels**, not a constant:
holding theta constant makes the near-surface layer exactly neutral, which the
convective adjustment then reads as marginal everywhere on step one. Wind is
held constant instead — extrapolating a shear downward produces surface winds
the drag scheme then fights.

`forecast.py` now builds a `PrimitiveSigma`, and puts the analysis through the
same **filter → rebalance** sequence measured in P-10, applied identically to
the initial state and to every boundary frame. Surface pressure is prognostic,
so the Davies relaxation drives it at the edges too.

**Confirmed by.** `test_interpolate.py` **7/7** — a field linear in log(p)
reproduced to 0.00e+00, source levels recovered to 3.6e-15, surface pressure
matching a standard atmosphere to **7.6 Pa** across 0–2500 m of terrain, and a
converted analysis integrating 6 h. `test_forecast.py` rewritten for the
sigma path, **11/11**.

**One defect found on the way.** The bracket search in
`surface_pressure_from_heights` had the height ordering backwards, which left
every column above the lowest analysis level pinned at that level's pressure —
a **253 hPa** error over 2500 m terrain. Caught by comparing against a
standard atmosphere, where the right answer is known in closed form. Category
D, and the fourth time an ordering convention has been the defect.

---

## P-05 — `forecast.py` ran over flat ground
**Category** I · **Status** FIXED · **Fixed** 2026-09-04

**Symptom.** `terrain=None` in the driver: every real-data run was over a flat
sea-level plain, in a domain whose defining feature is the Appalachians.

**Fix.** `ingest_hrrr.py` fetches orography once per run — it is static, so
per-hour fetching would be twelve redundant transfers on a shared link — and
writes `terrain.npz` beside the field files. Surface pressure is taken from
the same message set when available; when it is missing the heights supply it
instead, which is why a partial fetch does not stop a forecast.

The driver **refuses to run without terrain** rather than substituting a flat
domain, and names the ingest command in the error. A flat Northeast is not a
degraded forecast, it is a different experiment, and running it silently is
how a result gets misread later. A terrain/field shape mismatch is caught
separately and names `--stride` as the cause.

**Confirmed by.** `test_forecast.py`: surface pressure over an 1800 m ridge
measures **806 hPa** against a hydrostatic estimate of 801, with 996 hPa in
the surrounding lowland; a run directory without terrain is refused; a
mismatched terrain shape is refused.

**Not yet confirmed against the live service.** The GRIB search for surface
fields (`:(?:HGT|PRES):surface:`) and its alias table have never been run
against HRRR — see P-06. Written to fail loudly and specifically.

---

## P-08 — The decisive test integrated a clipped 166 m/s jet
**Category** E · **Status** FIXED · **Fixed** 2026-09-02

**Symptom.** A 2Δx mode growing near the domain boundary, doubling in ~20
minutes against hyperdiffusion tuned for a 3-hour e-folding. Read as a model
instability; survived nine single-candidate patches.

**Diagnosis.** The 6 K meridional temperature contrast implies a **166 m/s**
jet by thermal wind (Ro = 3.2). The test clipped the wind at ±60 m/s, which
destroyed geostrophic balance over **33.6%** of the domain. The clip was the
grid-scale source. The model was never broken.

**Fix.** Contrast reduced to 1.5 K — a realistic 41 m/s Northeast jet — and
the clip removed. `test_primitive_sigma.py` rewritten with the reasoning in
its docstring so the setup cannot quietly drift back.

**Confirmed by.** Flat ground, clean: **12/12 hours, max|u| unchanged at
41.5 m/s** where the clipped version reached 7/12 and 218 m/s.

---

## P-09 — Geostrophic wind from one PGF term over terrain
**Category** B · **Status** FIXED · **Fixed** 2026-09-02

**Symptom.** Terrain test cases initialized with absurd winds; every terrain
row of the boundary-layer baseline was measured against them.

**Diagnosis.** The balanced wind was taken as `-∂φ/∂y / f`. On sigma surfaces
the horizontal force has **two** terms that largely cancel over sloping
ground. Keeping only the first implies an **845 m/s** "balanced" wind over
2500 m terrain.

**Fix.** The initial wind now comes from the full `pressure_gradient_force`.

**Confirmed by.** Initial max|u| over 2500 m terrain: **845 → 41.3 m/s**, with
a physical 14.6 m/s cross-mountain ageostrophic component.

---

## P-10 — Unfiltered white noise in the initial state
**Category** E, G · **Status** FIXED · **Fixed** 2026-09-02

**Symptom.** 1.2 m/s of white noise killed every run within an hour,
regardless of mixing or drag setting.

**Diagnosis.** White noise puts **89%** of its variance at wavelengths the
grid cannot carry, and nonlinear advection amplifies it faster than
hyperdiffusion removes it (measured interior e-folding 10 800 s, 18 400 s next
to a replicate boundary). Real analyses are filtered before integration; this
one was not. The threshold is sharp:

| white noise | survived |
|---|---|
| 0.30 m/s | 12/12 |
| 0.60 m/s | 7/12 |
| 1.20 m/s | 1/12 |

**Fix.** `src/dynamics/initialization.py` — a raised-cosine spectral lowpass,
full response above 8Δx, zero at 4Δx, applied to u, v and the θ deviation from
the level mean. **Order matters and was measured, not assumed:** filtering
changes u, v and θ separately and so reintroduces divergence, which the
rebalance then removes.

| treatment | initial max\|div\| | survived |
|---|---|---|
| none | 3.90e-05 1/s | 1/12 |
| filter only | 9.93e-05 1/s | 11/12 |
| filter, then rebalance | 1.23e-05 1/s | **12/12** |

Note the middle row: filtering *raises* divergence and survives ten hours
longer. Divergence is not the controlling variable — wavenumber content is.

**Confirmed by.** Sub-4Δx wind rms 0.808 → 0.049 m/s; 12/12 at max|u| 42.4.

---

## P-11 — Mountain waves overturn with nothing to remove the instability
**Category** G · **Status** FIXED · **Fixed** 2026-09-03

**Symptom.** 2500 m terrain died at hour 12 with everything else fixed.

**Diagnosis.** Watched hour by hour, the wind never runs away — it sits at
41 m/s from the first hour to the last. What runs away is the stratification:

| hour | 1 | 6 | 8 | 9 | 10 | 11 | 12 |
|---|---|---|---|---|---|---|---|
| min Ri | 11.5 | 0.94 | 0.33 | 0.23 | **−0.05** | **−1.15** | dead |

Ri < 0 is N² < 0 — the mountain wave steepens as it propagates upward and
overturns. Correct physics, missing consequence. `eddy_diffusivity` does treat
Ri ≤ 0 as full-strength mixing, but it is a diffusion capped at 100 m²/s,
relaxing a 600 m layer in dz²/K = 3600 s. The wave steepens faster than an
hour. **The same failure shape as P-10** — a diffusion losing a race — in a
different scheme.

**Fix.** `src/dynamics/convection.py` — dry convective adjustment. Contiguous
unstable segments mixed to their mass-weighted mean, wind mixed over the same
layers so momentum is conserved and convective momentum transport is carried.
Applied as a **post-step adjustment, not a tendency**: an adjustment enforcing
an inequality has no meaningful time derivative, and inside the Runge–Kutta
stages an intermediate state would re-create the instability the final state
must be free of.

**Confirmed by** four predictions written before the run:

| prediction | outcome |
|---|---|
| min Ri floors near 0 | held — 0.09, 0.011, 0.020, 0.028 at hours 10–13 |
| Ri<0 count stops growing | held — 0 for the whole run |
| the run completes 12 hours | held — reached **16** |
| the wind is NOT damped | held — 41–42 m/s through hour 14 |

The fourth mattered most. A scheme buying stability by flattening the flow
would look identical in the first three; that is exactly how P-16 failed.
Conservation measured at 2.8e-16 (heat) and 2.6e-16 (momentum).

---

## P-12 — Convective adjustment converged like a diffusion
**Category** B · **Status** FIXED · **Fixed** 2026-09-03

**Symptom.** The first implementation mixed adjacent unstable pairs. Correct
and conservative, but a fully inverted 20-level column still had **0.26 K** of
spread after 200 sweeps.

**Fix.** Contiguous-segment mixing: a layer joins a segment if the interface
above or below it is unstable, segments are disjoint by construction, each is
mixed to its mass-weighted mean in one operation.

**Confirmed by.** Fully inverted column: 200 sweeps and 0.26 K residual →
**1 sweep and 0.00e+00**.

---

## P-13 — Convective adjustment became the dominant cost
**Category** H · **Status** FIXED · **Fixed** 2026-09-04

**Symptom.** A 12-hour 4000 m run that should take 20 minutes had not finished
in 100. Looked like a hang, not a performance bug.

**Diagnosis.** The sweep touched every column in the domain even though only
0.2–0.4% of interfaces are unstable at any step.

**Fix.** Compact to the columns that actually contain an inversion, run the
scan there, write back.

**Confirmed by.** 1.7 ms on a stable state (early exit), 270 ms with 0.3% of
the domain overturning — cost now proportional to the convection rather than
to the domain.

---

## P-14 — Pressure coordinates could not pose the lower boundary
**Category** B · **Status** FIXED · **Fixed** 2026-09-01

**Symptom.** Divergence within 2–3 forecast hours from real analyses at every
damping setting tried. `max|omega|` roughly quadrupling per hour.

**Diagnosis.** With a rigid flat lower boundary, ω must vanish at both ends of
a column that cannot move. Diagnosing ω from divergence and then forcing both
boundary conditions over-constrains the column; the correction that enforces
ω = 0 at the ground redistributes error through the whole column every step,
feeding a divergence–vertical-velocity feedback with no physical damping.

**Fix.** Terrain-following sigma coordinates with **prognostic surface
pressure**. σ̇ = 0 at lid and ground now falls out of the formulation rather
than being imposed, and the column exchanges mass through a moving surface.

**Confirmed by.** σ̇ at the boundaries verified to 0.00e+00; hydrostatic
consistency 2.1e-03 → **8.2e-15**; a motionless atmosphere over 4000 m terrain
drifts **0.009 m/s in 12 hours**, error linear in slope.

---

## P-15 — Initial analysis divergence
**Category** B · **Status** FIXED · **Fixed** 2026-08-30

Analysis winds are balanced for HRRR's discretisation, not ours. A Helmholtz
split removes the divergent component, solving the Poisson equation in Fourier
space with the eigenvalues of *our discrete* Laplacian so the cancellation is
exact rather than approximate.

**Confirmed by.** max|div| 7.5e-04 → 9.3e-05 1/s; implied ω **60 → 1.08 Pa/s**;
correlation with the rotational flow 0.997.

---

## P-16 — The sponge flattened the jet
**Category** E, G · **Status** FIXED · **Fixed** 2026-09-01

The first absorbing layer relaxed the wind toward the horizontal mean. That
absorbs the waves and also removes a jet, which is legitimate structure.

**Fix.** Relax toward a frozen reference state instead.

**Confirmed by.** The thermal-wind jet persists 24 h at |du|/|u| 0.24% and
spurious |v|/|u| 0.07%, against a jet that was visibly flattened before.

**Caught by.** The thermal-wind test, which is the reason it exists. This is
the failure mode P-11's fourth prediction was written to guard against.

---

## P-17 — Hyperdiffusion six times too weak
**Category** B · **Status** FIXED · **Fixed** 2026-08-31

The coefficient was derived from the continuous k⁴ instead of the eigenvalue
of the *discrete* biharmonic operator, giving an 18-hour e-folding at the grid
scale where 3 hours was intended.

**Fix.** `discrete_biharmonic_eigenvalue()`; the coefficient is now derived
from a requested damping time.

**Confirmed by.** Measured 2Δx e-folding **10 800 s** against the 3 h (10 800 s)
requested — previously 18 h. The measurement also produced the number that
later diagnosed P-10: damping at 3 h cannot hold a mode doubling in 20 min.

---

## P-18 — Hyperdiffusion on surface pressure is a mass source
**Category** B · **Status** FIXED · **Fixed** 2026-09-01

**Symptom.** Surface pressure inflating from 1088 to **1243 hPa** over three
hours.

**Diagnosis.** Hyperdiffusion is only conservative on a periodic domain.
Applied to prognostic surface pressure on a bounded domain it injects mass.

**Fix.** No diffusion on π at all. Grid-scale noise in π has to be controlled
by the wind field that generates it, not by diffusing mass. The reasoning is
in a comment at the call site so it does not get "fixed" back.

**Confirmed by.** Total mass conserved to **0.00e+00** relative drift over 12 h
on a periodic domain (`test_primitive_sigma.py`); surface pressure no longer
inflates.

---

## P-19 — Divergence damping violated its own stability limit
**Category** C · **Status** FIXED · **Fixed** 2026-08-31

Written as a tendency with a coefficient in m²/s. Explicit diffusion needs
ν·dt/dx² ≤ 0.25, and a coefficient chosen without knowing dt violates it.

**Fix.** Rewritten as a dimensionless post-step filter, stable by construction
for any dt. Default **off**.

**Confirmed by.** Stable at every dt tested, and the reason it stays off is
also a measurement: baroclinic growth **1.21×/day at 0.00, 0.33×/day at 0.01,
0.49×/day at 0.10**. Every setting that helps stability suppresses the weather.

---

## P-20 — GRIB search regex matched nothing
**Category** A · **Status** FIXED · **Fixed** 2026-08-29

`^(?:TMP|RH|...)` matched **0 of 708** messages: HRRR inventory entries begin
with a colon. Downloaded nothing and reported success.

**Fix.** `r":(?:TMP|RH|UGRD|VGRD|HGT):\d+ mb:"` — leading colon, no anchor,
with `test_hrrr_search.py` (6/6) pinning it.

**Confirmed by.** Matches **100 of 708** messages (5 variables × 20 levels), and
13 hourly files written where the previous run wrote none.

---

## P-21 — cfgrib renames variables to CF short names
**Category** A · **Status** FIXED · **Fixed** 2026-08-29

`KeyError` on every variable: TMP arrives as `t`, HGT as `gh`.

**Fix.** A `CF_ALIASES` table tried in order.

**Confirmed by.** All five channels extracted from a real GRIB file; ingest
runs end to end.

---

## P-22 — Herbie wrote to a path that did not exist
**Category** A · **Status** FIXED · **Fixed** 2026-08-29

Cache defaulted to `~/data`; the write failed silently and surfaced later as a
misleading `FileNotFoundError`. Reported as "It worked with nothing written".

**Fix.** Explicit cache directory under the run directory, created up front,
with `preflight.py` checking it before anything touches the network.

**Confirmed by.** 13 `.npz` files present on disk after the next run, where the
previous run reported success and wrote zero bytes.

---

## P-23 — `operator.py` shadowed the standard library
**Category** D · **Status** FIXED · **Fixed** 2026-08-30

A module named `operator.py` in `src/verification/` broke `collections`, which
broke `numpy`.

**Fix.** Renamed to `obs_operator.py`.

**Confirmed by.** `import numpy` succeeds from the package directory; the
verification suite runs 9/9.

---

## P-24 — Two `src` directories
**Category** I · **Status** FIXED · **Fixed** 2026-08-29

A `src/src` created during transfer; imports resolved to whichever came first.

**Fix.** Duplicate removed and the transfer path corrected.

**Confirmed by.** A single `src/` on the server; imports resolve to one file.

**Caught by the human, not by the AI or the suite** — the only defect in the
project so far detected that way.

---

## P-25 — `--dry-run` hammered the archive on failure
**Category** C · **Status** FIXED · **Fixed** 2026-08-30

The early return only fired on success, so a failed probe fell through and
attempted all 13 hours.

**Fix.** Return on both paths, alongside `netpolicy.py` (token-bucket limiter,
cache, sequential fetch) after the bandwidth constraint was stated.

**Confirmed by.** `test_netpolicy.py` 9/9; a failing dry run now issues one
request instead of thirteen, and sustained rate stays under the 8 MB/s ceiling.

---

## P-26 — Hourly snapshots landed at 0.86 h and 1.71 h
**Category** C · **Status** FIXED · **Fixed** 2026-08-31

`int(interval / dt)` truncated, so "hourly" output landed at 0.86 h and 1.71 h.

**Fix.** Output is emitted on **target times**, not step counts.

**Confirmed by.** Snapshots at exactly 1.00, 2.00, … h; `test_forecast.py` 7/7.

---

## P-27 — `grid.shift` used axis 1 for x
**Category** D · **Status** FIXED · **Fixed** 2026-08-31

Correct in 2D, where axis 1 is x. In 3D axis 1 is *y*, so every 3D x-derivative
was silently a y-derivative.

**Fix.** Axes mapped onto the last two dimensions regardless of rank.

**Confirmed by.** The 2D suite (8/8) and the 3D suite both pass against the
same operator — the reason the 2D-first order (prompt 42) paid off.

---

## P-28 — The stochastic filter could remove every mode
**Category** C · **Status** FIXED · **Fixed** 2026-09-01

A large enough length scale left the spectral filter with no modes, producing
a constant field with zero variance and no error raised.

**Fix.** Length scale capped; an impossible request raises.

**Confirmed by.** `test_subgrid.py` 7/7, including a case asserting the
perturbation field has non-zero variance at the largest permitted scale.

---

## P-29 — Ekman angle measured against the wrong reference
**Category** E · **Status** FIXED · **Fixed** 2026-09-02

The surface wind was compared to the wind at level 4 aloft, where thermal-wind
turning contaminates the measurement — the drag test read +5.0° with drag and
−10.7° without, which is backwards.

**Fix.** Compare against the **local** geostrophic wind at the same level,
computed from the PGF there.

**Confirmed by.** +26.4° cross-isobar with drag (speed ratio 0.47) against
+12.3° without (0.83).

---

# ELIMINATED

Investigated as causes of the terrain and noise failures, and ruled out. Kept
because a ruled-out candidate is the expensive part of a diagnosis and is
exactly what disappears from a repository.

| # | candidate | measurement that eliminated it |
|---|---|---|
| P-30 | reference-state PGF | no change in survival |
| P-31 | divergence damping | no setting helps without suppressing weather |
| P-32 | sponge strength | survival flat |
| P-33 | hyperdiffusion strength | 3.0 h / 1.0 h / 0.5 h all 2/6 |
| P-34 | level stretching | no change |
| P-35 | level count | no change |
| P-36 | Coriolis energy error | below the growth by orders of magnitude |
| P-37 | terrain smoothing | no change |
| P-38 | balanced surface pressure | no change |
| P-39 | **the timestep, at 4000 m** | dt and dt/2 identical to 4 significant figures for 6 hours |
| P-40 | **eddy-diffusivity ceiling** | K_MAX 100 / 300 / 1000 → 6/12, 6/12, 6/12 |
| P-41 | **lid height** | 200 hPa 11/12, 100 hPa 10–11/12, 50 hPa 9–10/12 — neutral to worse |
| P-46 | **initialization shock over terrain** | the conversion is accurate: geopotential error 0.01 m flat, 4.19 m over 2500 m terrain. The "shock" was a rest-start adjustment plus an inconsistent test analysis |

P-41 is a re-measurement. The original finding was recorded on a state now
known to carry P-08's clipped jet, so it no longer counted as evidence and was
re-run on valid initial states. The conclusion survived.

**P-46 did not survive, and the way it died is worth keeping.** It was opened
on a measurement of 9.1 m/s of "spurious wind" from a converted analysis
started at rest over a 1500 m mountain. Three hypotheses were tested in order:

1. *Geopotential mismatch against the analysis.* Exact hydrostatic inversion
   drove the error to **0.00 m** — and the acceleration did not move,
   2.70 → 2.67 m/s per hour. The hypothesis was wrong, and the inverted
   profile was statically unstable with an 8.5 K sawtooth, exactly as the
   inverse recursion's alternating mode predicts. (`hydrostatic_theta` is kept
   in `interpolate.py` with that warning in its docstring, unused.)
2. *Small-scale structure from the interpolation.* Horizontal filtering and
   one and three passes of vertical smoothing: 9.06 → 9.06 → 9.06 → 9.03 m/s.
   Nothing.
3. *The test.* The acceleration was almost entirely in dv/dt, and the state
   had a meridional temperature gradient and **no wind**. That is not a
   balanced state being corrupted; it is an unbalanced state being correctly
   adjusted. On FLAT ground the same setup drifts 8.98 m/s — terrain
   contributes 2.56 m/s of the 11.54 m/s at 2500 m.

The synthetic analysis was also not hydrostatically self-consistent: it
perturbed temperature by −1.5 K and height by −45 m independently, and those
are not in balance with each other. Rebuilding the heights as the hydrostatic
integral of the temperatures dropped the geopotential error from **140 m to
3.24 m**. Category E, the fourth test-design error of the project, and the
first one caught by the AI rather than by a human noticing an anomaly.

---

# REVERTED

## P-42 — Flux-form potential temperature
**Category** B · **Status** REVERTED

**Symptom.** Unstable within hours where the advective form was clean.

**Diagnosis.** the omega correction breaks discrete continuity. Reverted to the
advective form.

## P-43 — Simmons–Burridge vertical discretisation
**Category** B · **Status** REVERTED

**Symptom.** First attempt differenced half-level Φ and was **300× worse** than the scheme
it replaced. Reverted; the hydrostatic integration stayed as it was.

---

# ACCEPTED

## P-44 — The lowest model level sits at 237 m
**Category** G · **Status** ACCEPTED

Operational models put the lowest level at 10–50 m. With 20 sigma levels and
stretch 1.4 ours is at 237.7 m, which makes the log-law drag coefficient a
coarser approximation than it should be. Accepted for now because fixing it
means more levels, which costs runtime on a shared server, and the drag test
gives a physically correct Ekman spiral at the current spacing.

## P-45 — A single domain-wide roughness length
**Category** G · **Status** ACCEPTED

`z0` is one number. The Northeast domain runs from open sea (0.0002 m) to
forest (1.0 m) — four orders of magnitude — so a single value is poor near the
coast. Accepted until a land-use field is ingested.

---

## Maintaining this register

`python tools/problem.py new "<title>"` appends an OPEN entry with the right
shape. When a problem closes, edit its entry in place: change the status, add
what was done, and **add the measurement that confirms it**. A fix without a
number is an assertion, and this project has already been wrong nine times in
a row while feeling confident.

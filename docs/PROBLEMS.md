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

## P-52 — A dead run is not detected until the end of the forecast hour
**Category** E, H · **Status** OPEN · **First seen** 2026-09-12

**Symptom.** During the P-40 ladder, K_MAX = 110 spent **911 seconds of wall
clock** in forecast hour 7 without the hour finishing. It was not slow
physics. Inside that hour the wind reached **7176 m/s** and the adaptive
timestep had collapsed from 14.8 s to a mean of 0.87 s.

**Diagnosis.** Every ladder and sweep in this project tests for failure at the
HOUR BOUNDARY:

    m.run(3600, dt=dt)
    if not np.isfinite(m.u).all() or np.abs(m.u).max() > 150: break

`run()` itself only breaks on a non-finite surface pressure. A run that is
blowing up but still finite therefore keeps integrating, and because
`max_dt()` ratchets the timestep down as the wind grows — and never lets it
back up — the dying hour gets slower and slower. The state was already
meaningless at 150 m/s; the model then spent a quarter of an hour refining it.

**Why it matters beyond wall clock.**

1. Every sweep in this project is wall-clock bound, and the runs that cost the
   most are the ones whose output will be discarded. The 4000 m ladder spent
   more time on its two failing rungs than on the six that survived.
2. A survival count cannot distinguish "reached hour 8" from "reached hour 8
   at a millisecond timestep". Those are not the same forecast, and nothing in
   the harness currently reports the difference. That is L8.

**What is known.** dt, steps per hour and wall clock per hour are now recorded
by `kmax_binding.py`, and an hour that exceeds a wall-clock budget is reported
as STALLED rather than as survival. That is instrumentation in one script, not
a fix.

**What remains.** The failure test belongs inside `run()`, as a ceiling on
|u| checked with the same cadence as the adaptive dt re-check, so the
integration stops when the state stops being a forecast. It needs a stated
threshold and a switch, because a test that wants to watch a blow-up must be
able to turn it off.

**A caution found while instrumenting this.** The first version of the guard
integrated the hour in ten-minute chunks so it could check between them. That
changes the answer: `run()` truncates its final step to land exactly on the
requested duration, so chunking changes the step sequence, and at K_MAX = 110
the chunked and whole-hour runs diverge — min Ri 0.012 against 0.022 at hour 6
— and then fail differently. In a marginally stable regime, subdividing the
integration differently is not a neutral act. The working guard rides along as
a callback and only reads the clock.

**Update 2026-09-22.** `forecast.run_forecast` now checks for non-finite values and |u| > 150 m/s at the progress cadence (~200 checks per run), not only on the hour. It stopped P-56 at 3.75 h, mid-hour. The same check still does not exist inside `PrimitiveSigma.run()`, which the sweep scripts use.

---

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
convective adjustment, ~~the eddy-diffusivity ceiling~~ — the ceiling is back
on the list of live candidates as of 2026-09-10, because the ladder that
eliminated it was not varying anything (P-40, P-51). It is worth two forecast
hours at 4000 m, which is the only measured movement this problem has had.

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

## P-50 — The radiative lid still dies over tall terrain, at hour 4
**Category** G · **First seen** 2026-09-05 · **Status** OPEN ·
**Three causes found and fixed 2026-09-06, a fourth remains**

**Symptom.** 2500 m terrain, no sponge, radiative lid: **4/12** hours. Was
3/12 before this session's fixes; the sponge gives 12/12.

**What the spatial probe found, and what each fix did.**

| | edge/interior flux ratio | interior \|F\| by hour | dies |
|---|---|---|---|
| original | 1.5 → 2.3 → **42.3** | 0.67, 0.79, 0.77 (flat) | h 4 |
| + edge taper | 0.21, 0.13 | 1.25 → **4.59** | h 3 |
| + high-\|k\| cutoff | 0.15, 0.19 | 1.26 → 2.17 | h 3 |
| + loop opened, raw pi | 0.06 | **25.1** | h 2 |
| + loop opened, transient pi | 0.10, 0.19, 0.25 | 1.52, 2.22, 1.11 (bounded) | h 4 |

**Cause 1 — FFT periodicity at the lateral boundary.** The condition is
evaluated with an FFT, which assumes periodicity; the domain is `replicate`
with a mountain in it, so the field has a step across the wrap and multiplying
by \|k\| turns that step into ringing concentrated at the edges. Interior flux
was flat at 0.67-0.79 Pa/s the whole time — the boundary was amplifying its
own transform error, not radiating a wave. `remove_divergence_spectral` has
the same exposure and survives only because the lateral relaxation overwrites
those cells; the radiation flux goes straight into prognostic surface pressure
with nothing to protect it. **Fixed** by windowing before the transform and
tapering after it.

**Cause 2 — grid-scale feedback.** The transfer is proportional to \|k\|, so the
shortest waves get the largest flux. **Fixed** by a raised-cosine cutoff below
8 grid cells — which is physics, not a fudge: the hydrostatic radiation
condition is only valid where \|k\| << N/U, and that fails long before the grid
scale.

**Cause 3 — a closed loop through surface pressure.** The one that mattered,
and knob-tightening would never have found it:

    F  ->  dpi/dt  ->  pi  ->  phi_top  ->  F

phi_top is the hydrostatic integral from the ground up, so it moves when pi
moves, by about R T / p_s per pascal — roughly 0.94 m²/s² per Pa here. With
rho g / N about 164, the loop gain is **154 \|k\| per second**: an e-folding of
about 130 s for a 120 km wave. And the sign that correctly radiates waves is
the sign that makes this loop grow, so no choice of sign satisfies both.
**Fixed** by subtracting the column's own hydrostatic response — but only the
*transient* part of it. Subtracting the response to raw pi made things worse
(|F| jumped to 25 Pa/s, death an hour earlier), because over 2500 m terrain pi
varies by 27000 Pa from the mountain against a phi' of about 150 m²/s².

**Cause 4 — unknown.** With the flux bounded (1.52, 2.22, 1.11 Pa/s) and the
wind steady at 39 m/s for three hours, the run still ends at hour 4. Nothing
in the recorded diagnostics is running away before it does.

**What improved on the way.** Development is now **stronger**: the suite's
36 h / 18 h ratio went 1.91 → **3.12** against the sponge's 0.75. Every fix
here made the boundary better at its actual job while the terrain case moved
by one hour.

**Off by default.** The sponge remains production.

---

## P-49 — The sponge turns growing weather into decaying weather
**Category** G · **First seen** 2026-09-05 · **Status** OPEN

**Symptom.** Eddy kinetic energy of a growing baroclinic wave, 48 h,
48x48x20, 60 km:

| sponge levels | 6 h | 12 h | 24 h | 36 h | 48 h | |
|---|---|---|---|---|---|---|
| 0 | 6.5e+03 | 4.2e+03 | 3.7e+03 | 4.2e+03 | **1.2e+04** | grows |
| 2 | 3.6e+03 | 2.9e+03 | 1.7e+03 | 1.3e+03 | 1.6e+03 | decays |
| 3 | 2.8e+03 | 2.1e+03 | 1.0e+03 | 8.6e+02 | 7.9e+02 | decays |
| 5 | 2.0e+03 | 1.3e+03 | 8.1e+02 | 7.1e+02 | 6.9e+02 | decays |

**Every** sponge setting turns growth into monotonic decay — including two
levels of twenty. A model that cannot grow a baroclinic wave is not
forecasting weather, it is relaxing toward its initial condition.

**Diagnosis.** The lid is at 200 hPa, so the sponge sits at 301 hPa — in the
upper troposphere, where the upper half of a baroclinic wave lives. Baroclinic
instability is a coupled mode between an upper and a lower wave; damp either
end and it stops growing, whatever it is damped toward. Rayleigh damping
cannot distinguish a mountain wave from a baroclinic wave because in this
configuration they occupy the same levels.

**Four fixes tried, all measured, none works.**

| attempt | development (x/day) | reflection (max\|du\| 6 h) |
|---|---|---|
| no sponge | 3.18 | 60.8 |
| 5 levels, 15 min (the default) | 0.85 | 36.4 |
| rate 1 h / 6 h | 0.78 / 0.82 | 37.8 / 45.5 |
| running reference instead of frozen (6 h / 1 h) | 0.80 / 0.61 | 36.8 / 38.2 |
| divergent component only | 0.22% jet drift, but reflection **55.5** | worse than plain |
| lid 100 / 50 hPa | 0.97 / 1.89 | 53.2 / 60.3 |
| 26 / 30 levels at 50 hPa | nan / 2.45 | 65.7 / 56.4 |

The pattern is the same everywhere: whenever the sponge absorbs (reflection
21-36) development collapses; whenever development survives (1.89-2.45) the
sponge is not absorbing (50-60, which is the no-sponge value).

The divergent-only attempt is worth recording as a failed idea with a good
motivation: gravity waves are divergent and balanced flow is rotational, so
damping only the divergent part should have absorbed the wave and left the jet
alone. It absorbed less than a plain sponge, because a mountain wave is not
purely divergent and its rotational part reflected off the lid untouched.

**Interim decision: nothing changes, and that conclusion took a wrong turn
first.** The default was cut from 5 to 3 on the reasoning that a shallower
layer must do less damage. It does not. By the ratio that matters, three
levels decays FASTER than five:

| sponge | 48 h / 6 h |
|---|---|
| 0 | 1.83 (grows) |
| 2 | 0.44 |
| 3 | **0.28** |
| 5 | 0.34 |

The change also regressed two suites — the decisive noisy case 12/12 → 11/12,
and the Ekman spiral — and improved nothing, so it was reverted the same day.
Depth is not the mechanism, so trading depth buys nothing.

**What would actually fix it — built 2026-09-05, and it does fix this half.**
`src/dynamics/radiation.py`: the hydrostatic Klemp-Durran condition
w(k) = |k| phi'(k) / N, applied as a mass flux through the lid rather than as
damping. Measured, 48 h eddy kinetic energy, ratio of 48 h to 6 h:

| configuration | ratio |
|---|---|
| sponge 5, rigid lid | 0.34 (decays) |
| no sponge, rigid lid | 1.82 |
| no sponge, **radiative** | **2.50** |
| sponge 5, radiative | 0.37 (the sponge dominates) |

Development is not merely preserved, it is better than a rigid lid — an
amplitude e-folding of about 1.8 days against 6.7 for the rigid lid, and 1-3
days is what baroclinic waves actually do.

**It is not yet usable**, because it destabilises tall terrain: see P-50. Off
by default until that is resolved.

**How much this matters right now.** The target is 12-hour forecasts, and
these curves diverge most after 24 h — at 12 h the sponge=3 run holds 2.1e+03
against 4.2e+03. So current forecasts are affected but not invalidated. Any
extension past a day makes this the first thing to fix.

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

## P-53 — Observation-only initial state and frozen lateral boundaries
**Category** G, I · **First seen** 2026-09-22 · **Status** OPEN

**Symptom.** Not a failure yet: a design with a known cost. With no model output allowed and nothing observed after the cycle time, the lateral boundaries can only be held to the initial analysis. At ~20 m/s air crosses ~860 km in 12 h of a ~1300 km domain, so error from the frozen edges should reach much of the interior by hour 12.

**What is known.** Initial state from every reliable observation at the cycle time (`src/analysis/`); upper air from the previous run's forecast where soundings are missing; analysis area extended beyond the domain. Measurement to make: forecast error against distance from the nearest edge, by lead time (research log 2026-09-22, prediction P5).

**Ruled out.** none yet.

---


## P-54 — The radiosonde fetcher's request is rejected by IEM
**Category** A · **First seen** 2026-09-22 · **Status** OPEN

**Symptom.** `fetchers.raob_url()` sends `ts1`/`ts2` and several `station` values in one request. IEM now answers HTTP 422: it requires `sts`/`ets` (ISO, start before end), accepts one 4-character station per request, and wants the `K` prefix (`KOKX`). Every radiosonde fetch would have failed on the server.

**What is known.** Measured on the desktop, 2026-09-22 12Z the day before: with one request per station, 10 of 21 active IDs in and around the domain returned data; Albany, Wallops, CAR and ILN (12Z) and both Canadian sites returned nothing. `NORTHEAST_RAOB` also omits RNK (Blacksburg), which is inside the domain. Fix: new request builder in `src/analysis/sources.py`, one station per request, station list from IEM's RAOB network table.

**Ruled out.** none yet.

---


## P-55 — daily.sh looks for the ingested frames in the wrong directory
**Category** I · **First seen** 2026-09-22 · **Status** OPEN

**Symptom.** Predicted, not yet observed: `tools/daily.sh` gives `forecast.py --run-dir $DATA/tensors/analysis_<stamp>`, but `ingest_hrrr.py` writes to `config.TENSOR_DIR` = `$DATA/tensors_3d/analysis_<stamp>`. The forecast step should fail with 'No live_hrrr_f*.npz', and verify be skipped.

**What is known.** Found by reading the two paths side by side, 2026-09-22, while the first hand run (P-07) was in progress. Confirm or refute from that run's log. The new per-cycle script derives both paths from one variable.

**Ruled out.** none yet.

---




## P-57 — Divergence guard and progress log watch u only
**Category** E/H · **First seen** 2026-09-25 · **Status** FIXED

**Symptom.** Test W width 15, 18.00 h: the hourly log printed "max|u| 17.9 m/s" while the snapshot had max|v| 91.5 m/s. `run_forecast`'s in-loop guard (the 150 m/s ceiling, P-52) and its finiteness checks also looked only at `model.u`.

**What is known.** The first P-56 runaway started in v. A v runaway could therefore pass the ceiling, or go non-finite, without being reported until u followed. Found by comparing the log with `tools/locate_growth.py`, not by a test. Fixed in `src/forecast.py`: the guard uses max(|u|, |v|) and the finiteness of both, and the hourly line prints max|v|. Past divergence times were detected on u and may be late. Locations and ordering come from snapshots and are unaffected.

**Confirmed by.** `src/test_forecast.py` 11/11 after the change. The first server run with the fix must show max|v| in its hourly lines (test O).

**Ruled out.** none.
---


## P-58 — Maps show only F000: forecast snapshots are not on whole hours
**Category** A/E · **First seen** 2026-09-25 · **Status** FIXED

**Symptom.** Found by the user on the first server render (prompt 123): the viewer had only F000. `run_forecast` saves a snapshot at the first step at or after each output time, and the step is 17.1 s, so the snapshots sit at 1.0023 h, 2.0045 h, and so on. `make_maps.py` accepted a snapshot only if it was within 1e-6 h of a whole hour, so it kept none, and only the analysis (hour 0) was drawn.

**What is known.** The synthetic test forecast had exact whole-hour times, so it could not catch this. That is an interface assumption about what `forecast.npz` holds, made without reading `run_forecast`'s output rule (the class of P-54: the offline fixture was the AI's idea of the data, not the data). Fixed with `make_maps.match_hours`, which takes the nearest snapshot within 0.1 h. The synthetic forecast now uses the real step rule.

**Confirmed by.** `src/maps/test_maps.py`, "every forecast hour is found although snapshots land seconds late": snapshot times built by the model's rule give hours 1–16 (16 of 16). The synthetic re-render drew 25 times, 0–24.

**Ruled out.** none.
---


## P-59 — No diurnal cycle: the dry core has no surface heating or radiation
**Category** G · **First seen** 2026-09-25 · **Status** OPEN

**Symptom.** The first real server verification (06Z 2026-09-23, 346–362 surface temperature pairs an hour) has a bias that follows the sun. It is +1.7 to +2.2 °C in the night hours (F001–F005, 07–11Z) and crosses zero at about F006–F007. It reaches **−6.9 °C at F014 (20Z, 4 PM EDT)**, when RMSE is 7.7 °C. Wind speed does the same: +1.6 kt at night, −3.1 kt at F011 (17Z), when daytime mixing is missing.

**What is known.** `src/dynamics` has no surface sensible-heat flux, no radiation and no solar geometry. `radiation.py` is the lid's radiating upper boundary, not physics. So the lowest levels cannot warm by day or cool by night, and a forecast holds roughly the hour-0 temperatures while the real surface warms 6–8 °C into the afternoon. This is the largest error source in the first verified run, larger than anything P-56 has cost so far. Numbers are read from the error maps' header boxes (`maps/err_t_f*.png`); the artifact `verification_20260923_06Z.csv` holds them.

**Candidate responses (none tried).** A surface energy budget with solar geometry and a land/sea surface temperature. Or, as a first step, a prescribed diurnal surface heat flux from solar elevation, which uses no later observations. Either needs a prediction before it is built. Also worth adding: a persistence reference (the hour-0 analysis held fixed) scored the same way, so the model's skill is measured against doing nothing.

**Ruled out.** none.

**The night side (test X, 18Z 2026-09-22 run through the night, 2026-09-26).** With the damping on, the temperature bias rises from +1.3 °C at lead 5 to **+6.2 °C at lead 17 (11Z, around sunrise)**, with RMSE 7.3 °C. It then falls to −1.4 °C by lead 23 (17Z). Without surface cooling the nights stay far too warm, as the days stay too cool (−6.9 °C at 20Z on the 06Z run).

**Step 1, a prescribed diurnal surface heat flux (`--surface-heating`, 2026-10-01; tests AE, AE2, AE20, AG, AH).**
- Over 19 campaign cycles the temperature RMSE falls by 0.25 K on average (better at 21 of 24 leads, in 16 of 19 cycles). The model goes from never beating persistence to beating it at 17 of 24 leads.
- The swing of the bias narrows by about 23 % in the three single cases. The jet case's afternoon cold bias improves by 2 K.
- With momentum mixed by the adjustment, one cycle diverges (P-67).
- It is off by default. Remaining: the evening warm bias the next day (the heat is not removed at night), no clouds, and the wind is unscorable until P-68.

**The user's decision (2026-10-02).** Not the default: "we want a convective allowing model" (`docs/CAM_DESIGN.md`). A CAM takes its surface fluxes from a surface energy budget, not a prescribed curve.

**Step 2, a force-restore ground temperature (`--land-surface`, `src/dynamics/land_surface.py`, 2026-10-02; test AK).** Deardorff (1978) with clear-sky radiation, a bulk sensible heat flux, latent heat as a fixed share of daytime net radiation, and Louis (1979) stability on the drag. It replaces the prescribed flux. It is off by default.

**Test AK (20 cycles).**
- With the similarity operator, the 2 m T RMSE is 3.25 K, against 3.76 for AD and 3.56 for persistence: better than persistence at 18 of 24 leads.
- The day bias is +0.51 K and the night bias −1.33 K, against −2.30 and +1.54 for AD.
- At leads 19–24 the swing is too large (night −2.61 K, day +1.59 K). The model has no clouds, so every night is clear. That is the next limit, not a constant to retune.

---


## P-64 — Cycling never used the previous forecast: snapshots land one step past the hour
**Category** B, A · **First seen** 2026-10-01 · **Status** OPEN (fix in, awaiting test)

**2026-10-02.** The user kept the fix (test AF: no detectable skill difference, 1.578 against 1.599 K). It will close once a committed run shows the forecast first guess used at every chained cycle.

**Symptom.** In the test AD campaign every chained cycle logged a fallback. At 06Z and 18Z the log said "first guess: previous_analysis … (6 h old; previous forecast unusable)". At 00Z and 12Z it said "sounding_mean (previous forecast lacked the lead time)". The ingest documentation says the upper air of a sounding-poor cycle comes from the previous run's 6 h forecast. That has not happened in any cycle.

**What is known.** `run_forecast` emitted a snapshot at the first step at or past each target and stamped it with that step's time. At a 17.13 s step, the 6 h snapshot of a 24 h run is at 21 607.79 s, 7.8 s late (`times_s` of obs_20260926_00: 3614.28, 7212.98, …). `build.first_guess_from_forecast` accepts a lead only within 1e-3 h, 3.6 s. A discrete-time quantity was compared with a continuous one at a tolerance below the discretisation.

**Fix (in the working copy, not yet on the server).**
- (1) When the output interval divides the run, `run_forecast` takes a whole number of steps per interval: 211 steps of 17.0616 s instead of 17.13 s. It stamps each snapshot with its target, and the emission tolerance is half a step. Results therefore change slightly: the step is 0.4 % shorter.
- (2) `first_guess_from_forecast` accepts 0.02 h (72 s), more than one step, so older files qualify.
- Regression test `test_output_times_are_on_the_hour`: 24 snapshots, maximum offset 0.0 s, the 6 h snapshot at 21 600.000 s. `test_forecast.py` 15/15; `test_analysis.py` 21/21.

**Not yet measured.** Whether the forecast first guess gives a better 06/18Z analysis than the 6 h old analysis. The ingest's withheld-station score can decide it (rebuild with `--from-raw`). The user had also been asked to confirm the previous-analysis fallback, which until now has been the only path in use.

---


## P-65 — Surface-heating land mask from smoothed terrain counts half the ocean as land
**Category** G, D · **First seen** 2026-10-01 · **Status** OPEN (fix in, awaiting test)

**Symptom.** Test AE logged "land 85%" for a domain whose southeast part is the Atlantic and the Gulf of Maine.

**What is known.** The mask was `terrain > 0` on the run's `terrain.npz`, which is slope-limited (P-56). The smoothing spreads land heights over the sea: only 15.5 % of it is exactly 0 m (26.1 % ≤ 1 m), against 31.9 % of the raw ETOPO grid in `data/static`, where the ocean is clipped to 0. So about half the ocean was heated by day and cooled by night in test AE.

**Fix (working copy).** The land mask is read from the unsmoothed `data/static/terrain_etopo_<ny>x<nx>.npz`, with a warning fallback to the old mask. Test AE has to be repeated with it, and its first results are reported as contaminated by this.

---


## P-67 — Surface heating destabilises the 28 Sep 00Z cycle at 22.5 h
**Category** C, G · **First seen** 2026-10-01 · **Status** OPEN

**Symptom.** In AE20 (`--surface-heating`, correct land mask, PAV, new defaults), the 00Z 2026-09-28 run diverged at 22.46 h (367 m/s). The same cycle without heating (campaign AD) completed 24 h. max|u| holds at 14–17 m/s until 14 h, then grows through the next day: 22 m/s at 17 h, 24 at 21 h, 39 at 22 h, then it diverges, late afternoon local time (22Z).

**What is known.** `tools/locate_growth.py`:
- The largest hourly changes are in the lowest levels (L17–L19, about 100–300 m above ground) at rows 65–81, columns 41–72: northern New York and New England, terrain up to 515 m (mean 196 m), 15–30 cells from the edge.
- From the first hours, at night, hourly changes reach 4–8 m/s with 10–30 points over 5 m/s.
- From 12 h they reach 11–17 m/s, and at 21–22 h the change rises to 29–31 m/s at L15 (about 500 m).
- The unheated run has the same noisy patch (5–11 m/s hourly changes at L18, rows 65–80, columns 41–64), but it stays bounded (max|u| near 14 m/s throughout).

So the heating amplifies something already present rather than creating it.

**Candidates (test AG, one change each, heating otherwise on):**
- (a) the adjustment's instant mixing of momentum through the daytime mixed layer;
- (b) the steady night cooling over sloping terrain (drainage, or the sigma pressure-gradient error made worse by strong surface inversions);
- (c) the daytime heating itself.

**Test AG (one change each).**
- AG1, no momentum mixing: completes (max|u| ≤ 15.6 m/s).
- AG2, no night cooling: completes (≤ 15.4).
- AG3, night cooling only (no daytime heating; momentum still mixed): diverges at 17.48 h.

So it is the steady night cooling of the lowest layer together with the adjustment's instant momentum mixing; neither alone causes it. The candidate configuration is heating with `--no-conv-momentum` (test AH, 20 cycles). The cold-pool dynamics are not yet explained.

**Consequence.** `--surface-heating` stays off by default until this is understood.


**2026-10-02: not specific to the heating (tests AJ2 and the location check).** The 28 Sep 06Z cycle with `--z0 1.0` (no heating) diverged at 18.99 h, and 5 of 8 `--no-mixing` runs diverged.
- The z0 run and the P-67 heated run start in the same box. It is the second or third level above the ground at 45.2–45.4 N, 72.2–72.8 W (northern Vermont, about 20 cells from the northern edge, 170–310 m terrain), in the afternoon and evening of 28 Sep.
- Growth: 17 → 41 m/s between 15 and 18 h in the z0 run; 22 → 31 m/s at 17 h and 33 → 51 m/s at 21–22 h in the P-67 run.
- The control has its fastest low-level wind in the same box (16–21 m/s) and survives.
- ~~Next: a tendency budget at that point (pressure gradient, advection, drag, relaxation), not another variant.~~ done the same day:
- **The budget (`tools/budget28.py`, which closes exactly).** At a strict local speed maximum, the second-order centred horizontal advection supplies +7 to +66 m/s per hour along the wind, where the continuous term is zero. Pressure gradient and mixing oppose it, and hyperdiffusion removes only 2–7 m/s per hour. So the growth is the advection scheme's dispersion error. Test AL tries third-order upwind-biased advection (`--advection upwind3`).
- **Test AL (2026-10-02).** With `--advection upwind3`, both failing configurations (28 Sep 06Z with z0 = 1.0, and 28 Sep 00Z heated with momentum mixed) complete 24 h. The fastest wind anywhere stays at 15.6 and 14.9 m/s, against 41.5 and 51.1 m/s before divergence. On 20 cycles skill is unchanged (T +0.002 K, u +0.005, v 0.000 m/s), and runs are about 18 % slower. **Fix available.** It closes when upwind3 becomes the default (the user's decision).
---


## P-68 — Wind verified at the lowest model level (~300 m) against 10 m observations
**Category** A, E · **First seen** 2026-10-01 · **Status** OPEN

**Symptom.** Over the 20 campaign cycles, persistence at lead 1 has a pooled wind bias of −2 m/s in both u and v and an RMSE of about 3.5 m/s. Persistence at lead 1 is the analysis, which assimilated these observations.

**What is known.** By cycle, the persistence wind speed at the stations is 1.7–2.3 times the observed 10 m speed (11.1 against 4.8 m/s at 27 Sep 12Z; 2.0 against 0.9 at 30 Sep 06Z). The operator takes u and v of the lowest model level, about 300 m above the ground (`model_level_height_m` in the matches), with no reduction to 10 m. So every wind score, model or persistence, carries a speed-ratio error that has nothing to do with the forecast. The temperature operator has a related known weakness: a standard lapse rate from 300 m to the station, wrong in a mixed layer by day and in an inversion by night.

**Step 1 (2026-10-01).** `verify.py --wind-operator log10m`: the model's neutral log law from z1 above the model ground to 10 m (factor about 0.59), with z0 = 0.1 m as in the model's drag. Default unchanged.

**Test AI, the 20 campaign cycles re-scored from the archived observations:**
- Persistence lead-1 speed ratio 1.86 → 1.10.
- Persistence lead-1 RMSE: u 3.49 → 2.28, v 3.11 → 1.94.
- Temperature identical.
- The model still loses to persistence in wind (u +0.34, v +0.15 m/s pooled), and its speed ratio grows from 1.17 to 1.51 over 24 h, against 1.10 to 1.32 for persistence.

~~**Open.** Whether to make `log10m` the default (the user's decision). A stability-dependent reduction, since the neutral log law over-reduces above about 50 m. The model's wind growth with lead time is a model error, not an operator one.~~ (superseded 2026-10-02, below)

**`log10m` is the default (2026-10-02, the user's decision).** Archives verified before that date used `lowest` unless their name ends in `_w`.

**Step 2, a similarity operator (`--surface-operator similarity`, `src/verification/surface_similarity.py`).** Monin–Obukhov between the model's ground temperature and the lowest level, for both the 10 m wind and the 2 m temperature. It needs a `--land-surface` forecast. The time-of-day split showed that the "growth with lead time" was mostly a night-time excess (P-69), so the neutral operator was wrong at night in the way the old one was wrong at all hours.

---

## P-69 — The 10 m wind is 1.6–1.8 times the observed at night: no stable surface layer
**Category** G · **First seen** 2026-10-02 · **Status** OPEN

**Symptom.** In the 20 campaign cycles under the 10 m operator, the dry model's 10 m wind is about right by day: speed ratio 1.10, against 0.95 for persistence. At night it is 1.63 (1.37 at leads 1–6, 1.81 at 13–18), against persistence's 1.41. The model's mean speed is about 3.9 m/s at all hours, while the observed falls from 3.6 by day to 2.2–2.7 at night. Test AI's "growth with lead time" (1.17 → 1.51) is mostly the share of night hours in each lead.

**What is known (tests AJ, AJ2).**
- Against the analysis valid at the same time, the model's lowest-level wind exceeds it by 1.4 m/s at 6 h and 2.3 m/s at 24 h.
- The excess is confined to the lowest two levels. It is largest over the sea (ratio 1.55 at 24 h) and smallest in the relaxed edge zone (1.19).
- The sea result rules out one too-small roughness length as the cause.
- The 2 m temperature shows the same day/night split: bias −2.30 K by day and +1.54 K at night. The prescribed heating fixed part of the day and none of the night.
- Common cause, to be tested: the ground has no temperature. The surface layer never becomes stable, the drag is always neutral, and both 10 m operators assume a neutral profile.

**Test AJ2 (8 cycles).**
- Richardson mixing and convective momentum mixing are not the source: removing either leaves the excess unchanged.
- Doubling the drag (z0 = 1.0) removes 66 % of the lowest-level excess at 6 h and 33 % at 24 h, and lowers wind RMSE by 0.56 (u) and 0.26 (v) m/s.
- By 24 h the whole lower troposphere is too fast (+2.5 m/s at levels 13–18), which points at the resolved dynamics.

**Test AK (20 cycles).** With the force-restore ground temperature and the similarity operator, the night ratio falls from 1.63 to 0.92 and the day ratio is 1.22. Wind RMSE falls against AD (u −0.107, v −0.163 m/s), but u is still 0.23 m/s worse than persistence. Most of the night excess was the neutral 10 m diagnosis of a surface layer that should be stable. **Status:** largely explained; it stays open until the land surface is adopted or rejected (the user's decision) and the remaining u error is understood (too little turning across the isobars, test AJ2).

---


# FIXED

## P-66 — Verification takes 10-14 min per run: all-pairs buddy search in QC
**Category** H · **First seen** 2026-10-01 · **Status** FIXED · **Fixed** 2026-10-01

**Symptom.** On the server each `verify.py` run spent 10–14 min after the ASOS fetch, which itself took about 1 s. The fetch-free `--report-only` runs were just as slow. A campaign of 20 cycles × (model + persistence + heated) would have taken about 13 h of single-threaded verification. The production cycle verifies once a day, inside a 1.5 h budget.

**What is known.** `run_qc` (`src/verification/observations.py`) finds each observation's buddies by scanning every other observation in Python (same variable, same level, within 30 min and 150 km). That is O(n²) over the 10⁴–10⁵ observations of a 24 h window. The same function runs in the ingest (64–195 s per cycle) and in `make_maps.py`.

**Fix (working copy).** `buddy_neighbours` groups by variable, sorts by time, finds the 30-minute window by binary search, and applies the level and distance tests with NumPy, using the same formulas in the same order. The old search stays as `buddy_neighbours_reference`. Test `test_fast_buddy_search_matches_reference` (3 672 observations, 44 200 neighbour links, 56 rejections): identical neighbours, in the same order, and identical verdicts; 0.07 s against 3.47 s. `test_verification.py` 10/10; `test_analysis.py` 21/21.

**Confirmed by.** A real 24 h ASOS window on the server (`asos_2026092601_2026092700`, 50 170 observations):
- fast search 7.6 s, all-pairs reference 849.1 s;
- identical verdicts (0 differ; 645 buddy rejections in each).

---

## P-63 — Calm cases run 3.8x slower: the convective adjustment hits its 20-sweep cap every step
**Category** H, C · **First seen** 2026-09-26 · **Status** FIXED · **Fixed** 2026-10-01

**Symptom.** At torch × 8, calm cases (06Z 2026-09-23, 18Z 2026-09-22) run about 3 steps/s; the jet case (12Z 2026-09-25) runs 11 steps/s. That holds with each run alone on a quiet server (test Y: 2.9 against 11.1) and with denormals flushed (3.0 against 11.2).

**What is known.** Test Z counted and timed `dry_convective_adjustment` over 1 h:
- **Calm case:** 20 sweeps, the cap, on all 210 steps, taking 51.6 s of about 72 s of stepping. About 0.8 % of interfaces are unstable at the start of every call (mean 8.19e-3, max 9.16e-3), nearly constant through the hour.
- **Jet case:** 0 sweeps on all 229 steps, 0.2 s in total.

So the calm case never finishes adjusting within a step, and the instability is back at the next call. A sweep is two Python loops over 20 levels, several hundred small array operations. On a desktop test with 3 K noise the scheme does converge, but it needs 69 sweeps: mixing uncovers new unstable interfaces at segment edges. The cap of 20 is therefore not a guard but binding.

**Candidates (test Z2, `tools/convection_check.py`).** (a) The relaxation re-imposes an unstable frozen boundary state after every step, so the adjustment fights it indefinitely. The fix would be to adjust the driving state once. (b) The interior instability needs more than 20 sweeps. The fix would be a converging algorithm, and a cheaper sweep.

**Z2 result: (a).** In the calm 1 h state, all 1 436 unstable interfaces (7.08e-3) lie within 14 cells of the edge: 747, 470 and 219 at distances 0–4, 5–9 and 10–14, and none deeper. All of them are in the lowest two interfaces, and the largest violation is 0.471 K. The jet state has none. The run has one driving frame, the 18Z analysis, held for 24 h. Its afternoon superadiabatic surface layer is pulled back into the edge columns by the relaxation after every step. With a higher cap the adjustment converges in 30 sweeps, so the cap of 20 is too low for this state, but the state comes back every step.

**Fix (awaiting test AA).** `stabilise_frame` in `src/forecast.py` adjusts each boundary frame once, to convergence, and returns a stable frame unchanged. `--raw-boundaries` gives the old behaviour. Regression tests: `test_stable_frame_is_unchanged` and `test_relaxation_toward_stabilised_frame_stays_stable` (`src/test_forecast.py`, 14/14). It is confirmed if AA1 runs at 8 steps/s or faster and its 1 h state has 0 unstable interfaces.

**Test AA: the frame fix is right, but it is not the cost.** The frame was unstable (9.14e-3), the fix cleared it, and skill moved by at most 0.02. The jet case stayed bit-identical. The speed stayed at 3.2 steps/s (06Z case: 2.9).

**Test AB (class split at every call): the algorithm is the cost.** Every call in both calm cases runs to the 20-sweep cap and leaves instability behind. The violations at call start are 1e-8–1e-3 K (96 % and 79 %), not round-off. They sit on layers the adjustment itself mixed to neutral, which each step nudges by about 1e-4 K. Clearing them needs 23–24 sweeps. The jet case has no neutral layers, so it costs nothing. **Fix candidate:** pool-adjacent-violators (`dry_convective_adjustment_pav`, `--conv-scheme pav`), exact in one pass. It mixes less than the sweep scheme (the least-change stable state), so it changes the answer, and test AC judges it on skill.

**Ruled out.** Machine load and denormals (test Y).

**Fix.** Two parts. (1) The boundary frame is adjusted once, to convergence (`stabilise_frame`, test AA), so the relaxation no longer re-imposes an unstable layer. (2) The convective adjustment is pool-adjacent-violators (`dry_convective_adjustment_pav`, now the default; `--conv-scheme sweep` gives the old scheme). PAV is exact in one pass, the least-change stable state.

**Confirmed by.** Test AC (2026-10-01), on the server against the sweep scheme:
- Calm 18Z 2026-09-22, 24 h: 8.1 min (26.6 with the sweep scheme); 0 unstable interfaces in all 24 hourly states (the sweep run: in 24 of 24).
- Calm 06Z 2026-09-23: 7.8 against 27.2 min.
- Skill identical to the 0.01 printed at every lead in both.
- The jet case is bit-identical.
- `test_convection.py` 12/12 where torch is installed (the desktop torch environment and the server, 2026-10-01), and 11/11 without torch, where the torch-agreement test is skipped (the run after PAV became the default). The suite calls both schemes directly, so the default does not enter it. That includes the P-63 near-neutral chain the 20-sweep cap leaves behind.

---

## P-62 — Every forecast lost its final output time (float sum of steps)
**Category** D · **First seen** 2026-09-26 · **Status** FIXED · **Fixed** 2026-09-26

**Symptom.** 24 h runs with 15-minute output wrote 95 snapshots ending at 23.75 h, not 96 ending at 24.00 h, and 8 h runs ended at 7.75 h. The re-verified test X therefore had no lead 24.

**What is known.** `run_forecast` emits a snapshot when `model.time >= target - 1e-9`. `model.time` is a sum of n_steps floats (5 040 steps of 17.1 s) and ends a few nanoseconds short of the duration. So the last target, which no later step can cross, was never reached. Intermediate targets were unaffected.

**Fix.** The last step always writes the final pending target.

**Confirmed by.** A stub model with 17.1 s steps reproduced it: 95 snapshots, the last at 23.753 h. With the fix, `test_final_output_time_is_written` in `src/test_forecast.py` gets 96, the last at 24.000000 h. The suite passes 12/12.
---

## P-61 — Verification drops every hour after 00Z (date-only ASOS request)
**Category** A · **First seen** 2026-09-26 · **Status** FIXED · **Fixed** 2026-09-26

**Symptom.** Test X verified two 24 h forecasts from 18Z 2026-09-22, but the scores stop at lead 6.25 h, which is 00:15Z on 2026-09-23. Lead 6 has 281 temperature pairs against about 357 at leads 1–5. Test W's comparisons were complete: the Q case to 8 h and the 06Z case to 16 h. Both windows stayed inside one UTC day.

**What is known.** `src/verification/fetchers.asos_url` sent only `year1/month1/day1` and `year2/month2/day2`. IEM then ends the request at 00Z of the end date. So any verification window that crosses midnight UTC loses every hour after 00Z. Almost every 24 h cycle crosses it. The analysis ingest (`src/analysis/sources.asos_url`) was never affected: it sends exact `sts`/`ets` timestamps.

**Fix.** `fetchers.asos_url` now sends `sts`/`ets` to the minute, as the ingest does. `test_fetchers.py` checks a window crossing midnight (17:30Z to 18:30Z the next day) and fails on a date-only request; the suite passes 9/9. The old archives hold the truncated observations under the same window name, so they cannot be reused.

**Ruled out.** none.

**Confirmed by.** Test X re-verified on the server into fresh archives (x0b, x1b). X1 is scored at every lead from 1 to 23 h, with 342–362 temperature pairs per lead; before the fix it stopped at 6.25 h. X0 is scored to 14 h, its last snapshot. Lead 24 was missing because the forecast itself ended at 23.75 h (P-62), not because of the request.
---

## P-56 — The first observation-built forecast diverges at 3.75 h
**Category** F?, G? · **First seen** 2026-09-22 · **Status** FIXED · **Fixed** 2026-09-26

**Symptom.** 2026-09-21 12Z from observations only, 24 h requested: max|u| 46 m/s for 3 h, then 372 m/s at 3.75 h, stopped inside the hour by the new guard. Desktop, 12 km grid, ETOPO terrain 0–1161 m.

**What is known.** Probe (`src/analysis/probe_obs_blowup.py`, 5-min snapshots): the runaway is the meridional wind at level 17 of 20 (near the ground), first growing by more than 20 % between 2.50 and 2.59 h, 14 cells from the edge, over 812 m of terrain in northern Maine. Only 2 points grew by more than 5 m/s, in a ~5.5 Δx pattern. The initial divergence did not reach the filter's target (4.7e-4 → 9.2e-5 1/s). Next: record every field around (row 82, col 89) from 2.0 h, check the static stability of that column in the initial state (the surface blend adds increments to the lowest 1000 m), and compare with an HRRR-seeded run of the same cycle (`--source hrrr`), which has survived 12 h before.

**Ruled out.** the edges: max|u| was pinned at the lid by the frozen boundary, and the runaway began 14 cells inside it.

**Second case, server, 2026-09-22 18Z.** Standard-atmosphere first guess (no soundings, no previous run), initial max|u| 6.2 m/s; the wind grew steadily from hour 3 (10 → 18 → 15 → 38 m/s) and reached 435 m/s at 6.31 h. A near-calm start dying rules out jet strength. Shared by both cases, and not by the HRRR runs that survived 12 h: the observation-built lower atmosphere (surface blend), ETOPO block-averaged terrain, and frozen single-frame edges.

**Ruled out, test A (2026-09-22).** The surface blend: with it switched off the run still diverged, at 4.10 h, first growing at 2.84–2.92 h in v at level 17 over 712 m of terrain in central NY, 47 cells from every edge.

**Terrain implicated (2026-09-22).** The same state over flat ground held max|u| at 46 m/s to 12 h, then went non-finite at 13.0 h (a sudden negative pressure or column depth, a different signature). Next suspect: below-ground pressure-level values in the observation analysis, which HRRR supplies smoothly and the analysis does not. Test and prediction in the research log.

**Ruled out, below-ground values (2026-09-22).** Extrapolating the 7507 pressure-level values under the terrain (up to 7.55 K) left the run dying at exactly 3.75 h. Next: test B, HRRR terrain under the observation state (server).

**Reframed, 2026-09-23.** Test C, an HRRR start with frozen edges, died at 3.16 h, and the record shows no real-data forecast of the sigma core had ever run before 2026-09-22. So this is not an observation problem: no real state over real terrain has survived. The earlier text's comparison with "HRRR runs that reached 12 h" was wrong; those runs never existed. Leading mechanism: slope. Real terrain at 12 km reaches 0.0316, against 0.0086 for the steepest idealised terrain that survived 12/12. Test S (slope-limited ETOPO) in the research log.

**Early failure explained (test S, 2026-09-23).** ETOPO at 12 km reaches a slope of 0.0536 (model measure), against 0.0086 for the steepest idealised terrain that survived. Smoothed to 0.0083 (11 passes, peak 1161 → 881 m), the same state survived to 13.7 h instead of 3.75 h. Ingest now limits slope to 0.0086 by default. **What remains** is a second failure at hours 13–14, seen on both flat and smoothed terrain: theta minimum falling and sigma_dot growing from hour ~9. Frozen edges are the first suspect.

**Second failure located (2026-09-23).** The first growth came at 11.59–11.67 h, 10 cells from the western (inflow) edge: the inner boundary of the 10-cell relaxation zone. It is violent (theta ±6.7 K and pi 6.2 hPa in 5 min) and near grid scale. Mechanism: the frozen edge pulls the inflow boundary back to hour 0 while the interior has moved on. Candidates: a gentler/wider relaxation for a frozen driver, or inflow-only relaxation. Not yet tested.

**Second case (server 06Z 2026-09-23, 16.31 h).** The largest change stayed 10–13 cells from the edge, almost always at the south-west corner of the zone boundary over 540–880 m terrain (2 of 30 early, small maxima were further north on the same column), from hour 1, in near-windless flow (standard-atmosphere first guess). So an inflow mismatch is not required. Open: the zone boundary itself (H1) or the terrain at that place (H2). Test W, relaxation width 6 vs 15, separates them.

**Test W (2026-09-25).** The growth moves with the relaxation width: edge distance mostly 6–8 at width 6 and 15–17 at width 15. Width 6 ran away over flat coastal ground, so terrain is not required. The failure is made at the inner boundary of the zone, a band pinned to hour 0 beside a free interior. Width 15 delays it by about 2 h. Next: test O, no relaxation (O1) and alpha 0.1 (O2).

**Test O (2026-09-25).** With no relaxation the run dies at 6.67 h at the physical southern edge (edge distance 0), so some zone is needed. With alpha 0.1 at width 10 it lasts to 21.67 h (from 16.31 h), still failing about 10 cells in at the south-west corner. Strength of the pull toward the frozen state is the strongest control found (+5.4 h); width 15 gave about +2 h. Next: test P (width 15 + alpha 0.1; alpha 0.03), then a 12Z case with soundings before any default changes.

**Test P (2026-09-26).** Width 15 with alpha 0.1 completed 24 h on the 06Z case with no 15-minute change above 5 m/s anywhere: the first clean 24 h real-data forecast. Width 10 with alpha 0.03 also reached 24 h, but the zone-boundary growth was rising at the end. Next: test Q, a 12Z case with soundings, before changing the default.

**Test Q (2026-09-26).** On a 12Z case with soundings, width 15 with alpha 0.1 did **not** survive: it diverged at 12.65 h against 7.45 h for the default. Both failed from an interior jet-level disturbance at 4 h that the zone does not control (P-60). The default is unchanged: width 15 with alpha 0.1 removed zone-boundary growth on one case and did nothing for the other failure.

**Seen in P-60 test S (2026-09-26).** With vertical mixing off, the south-east zone-boundary point r10–12 c97–99 (edge 10, L03–L05, over the sea) ran away first: v changed by 36.5 m/s in 15 min at 3.75–4.00 h. With default mixing it peaks at 3.3 m/s, so mixing is holding back a zone-boundary instability at that corner.

**Test W (2026-09-26).** With divergence damping, the 06Z 2026-09-23 case runs 24 h clean with both the default zone (W3, which was predicted to fail) and the wide weak zone (W4). The Q case also completes (W1, W2). Width 15, alpha 0.1 and C = 0.0064 are now the forecast defaults. **Still OPEN** until 18Z 2026-09-22, the first server cycle, rebuilt from raw, completes with them (test X). The 12Z 2026-09-21 case cannot be rerun, because its raw files were in the wiped desktop scratch data.

**Fix.** The forecast defaults since 2026-09-26: divergence damping (Skamarock and Klemp 1992) at C = 0.0064 (ν about 5.4–5.9e4 m²/s), and the relaxation zone widened and weakened to width 15, alpha 0.1. Damping alone was enough on the 06Z case (W3). The wider, weaker zone removed the remaining zone-boundary activity on the Q case (W1 against W2).

**Confirmed by.** Every real case that can still be run, each run 24 h, old settings against new:

| Case | Old settings (width 10, alpha 1, no damping) | New defaults |
|---|---|---|
| 06Z 2026-09-23 (calm) | diverged 16.31 h (south-west zone corner) | W4: 24 h, largest 15-min change 0.5 m/s |
| 12Z 2026-09-25 (jet) | diverged 7.45 h (P-60) | W2: 24 h, no point changing > 5 m/s |
| 18Z 2026-09-22 (calm; rebuilt from raw, test X) | X0: diverged 14.05 h (south-west zone corner, 9–12 cells in, from about 9.75 h) | X1: 24 h, no point changing > 5 m/s |

Skill before the old run fails is unchanged within ± 0.2 (K, m/s) over hours 1–4 (Q), 1–16 (06Z) and 1–6 (18Z). The largest degradation is +0.11 (v, Q case, lead 1). The largest change, −0.17 m/s (u, 06Z), is an improvement. Not rerun: 12Z 2026-09-21, the case that opened this entry, whose raw observations went with the wiped desktop scratch data. Reopen if a live cycle diverges with the new defaults.
---

## P-60 — Jet-level instability 17–21 cells inside the domain kills the 12Z case at 4 h
**Category** C · **First seen** 2026-09-26 · **Status** FIXED · **Fixed** 2026-09-26

**Symptom.** Test Q started 2026-09-25 12Z from real soundings (first guess `sounding_mean`, max|u| 27.6 and max|v| 33.2 m/s). With the default zone it diverged at **7.45 h**; with width 15 and alpha 0.1 it diverged at **12.65 h**. In both runs the change is small (≤ 3.3 m/s per 15 min) until **4.00–4.25 h**. Then it grows at the same place: levels **L03–L05 (about 270–330 hPa, jet level)**, rows 76–78, columns 86–90 (**45.3–45.5 N, 69.4–68.8 W, central Maine**), **17–21 cells from the nearest edge**. From 5 h on, 5 100–11 500 points change by more than 5 m/s every 15 minutes, so neither forecast is usable after about 4.5 h.

**What is known.** The onset does not depend on the relaxation zone. The same time, place and levels appear with the default zone (10 cells) and with width 15 and alpha 0.1. The onset point lies 7–11 cells inside the default zone's inner boundary and 2–6 inside the wide zone's, so this is not P-56's zone-boundary growth. The numpy and torch runs of the default case (Q0n, Q0t) differ by 1.6e-13 at 0.25 h. That difference grows with an **e-folding time of 24 min** from the start (2.0e-9 at 4 h) and 33 min afterwards: an unstable mode is present from hour 0. On the 06Z 2026-09-23 case the same pair grew with an e-folding time of 263 min for 8 h, then 74 min. An e-folding time of 24 min is faster than inertial instability can grow (its growth rate is at most about f, an e-folding of about 3 h). Shear instability (Ri < 0.25; Miles 1961; Howard 1961), static instability, or a numerical mode can grow that fast. Levels L00–L04 are the wind sponge, so the onset band L03–L05 straddles the sponge base, as P-56's growth sat at the lateral zone's inner edge.

**Test R (2026-09-26).** The mode is in the Maine box from hour 1. It spans the lid (L00) to L04, where the box's strongest wind (32.6 m/s) is at the lid itself, and its e-folding shortens from 46 to 17 min over hours 0.5–4. It is still in the linear range, so the flow under it becomes more unstable over those hours.

**Measurement R2 (1–4 h).** Round the mode, the minimum Ri at L03/L04 falls from 0.91 to 0.44, 0.30 and 0.28. N2 stays positive and eta/f above 0. Domain points with Ri < 0.25 at L03/L04 go from 0 to 19 to 579 by 4 h, just before the runaway. Mixing is exactly zero for Ri ≥ 0.25, so the sharpening layer has no vertical dissipation.

**Candidates.** (c) The rigid 200 hPa lid cutting through the jet. It cannot be tested cleanly, because the analysis stops at 200 hPa. **Test T (2026-09-26).** Halving the timestep leaves the onset at 4.00–4.25 h. Hyperdiffusion ×4 delays it to 5.75–6.00 h, and the mode's roughness is 0.99 (2–4Δx). P-60 is a grid-scale (about 2–3Δx) spatial-discretization mode, independent of the timestep. **Measurement U (2026-09-26).** The one steady energy source is the mode's own vertical motion acting on the jet's vertical shear, d(sigma_dot) dU/dsigma, at 1.5–2.5e-3 s⁻¹ at 2, 3 and 4 h. The measured energy growth is 1.5–1.9e-3 s⁻¹. The pressure-gradient term swings in sign, as an oscillating mode's does. This happens while Ri ≥ 0.28, so the discrete mode feels less restoring force than Ri implies. Candidate (f): the Lorenz-grid vertical computational mode (Arakawa and Konor 1996). **U2 (2026-09-26).** The time mean holds up: growth 1.08e-3 s⁻¹, and the shear term is 71 % of the total. Theta's adjacent-level correlation is 1.00 at every level, so the mode is smooth and deep, not a Lorenz zigzag: (f) is refuted. Candidate (g): a divergent grid-scale gravity-wave mode. Test V: divergence damping (Skamarock and Klemp 1992), `--div-damp` 0.01 and 0.003.

**Test V (2026-09-26): (g) holds.** With divergence damping at ν = 5.9e4 and 1.8e4 m²/s, no point changes by more than 5 m/s in 0–6.25 h, and both runs complete 8 h. The resolved maxima (27.6 and 33.2 m/s) are unchanged. The mode's divergence/vorticity ratio is 2–52 at L00–L05. (ν was 0.64 of the intended value because of a bug, since fixed; the logs' ν is the one quoted.) Next, test W: 24 h on both cases, with and without the wide weak zone, plus verification. **Candidate treatment, not yet a default.**

**Ruled out.**
- The lateral relaxation settings as the cause of the onset (Q0 vs Q1: same onset).
- The torch backend. Q0n and Q0t diverge at the same 7.45 h, and their difference is round-off amplified by the mode.
- (a) An unstable flow at the start. At t+0.25 h round the mode, eta/f ≥ 0.66, Ri ≥ 1.59 and N2 > 0 at every level, and no point in the domain has Ri < 0.25 at L00–L18.
- (b) The sponge base. With 8 and 3 sponge levels the onset stays at 4.0–4.5 h in the same place. With 8 levels it moves up to the lid, not down with the sponge base.
- Time discretization. With half the timestep the onset is unchanged (test T1).
- (f) The Lorenz-grid computational mode. Theta in the mode has an adjacent-level correlation of 1.00, with no zigzag (U2).
- (e) Missing vertical dissipation, and vertical mixing generally. With Ri_c 1.0, and with mixing off, the onset is unchanged (4.00–4.25 h at r77 c88; test S).

**Fix.** Divergence damping (Skamarock and Klemp 1992), `subgrid.divergence_damping`, as a forecast default: C = 0.0064, ν = C dx dy / dt ≈ 5.9e4 m²/s. It damps the divergent wind only: 2Δx in about 10 min, 10Δx in about 105 min. The wider, weaker relaxation zone (width 15, alpha 0.1) became the default at the same time, for P-56.

**Confirmed by.** Q case, 12Z 2026-09-25:
- V1 (ν 5.88e4) and V2 (ν 1.77e4): no point changes by more than 5 m/s in 0–6.25 h, where the default blew up from 4.00 h; both complete 8 h.
- W1 (default zone): 24 h with no interior onset; the only fast points, at most 6, are at the zone boundary.
- W2 (new defaults): 24 h with no point changing by more than 5 m/s.
- Skill before the undamped run's onset (hours 1–4) is unchanged within 0.11 (K, m/s).
One case: a second jet case would strengthen it.
---

## P-40 — The eddy-diffusivity ceiling was binding, and was never tested
**Category** E, G · **Status** FIXED · **Fixed** 2026-09-12 ·
**Eliminated in error 2026-09-04, reopened 2026-09-08**

**Symptom.** A K_MAX ladder at 4000 m returned identical survival at every
setting — 6/12 at 100, 300 and 1000 m²/s — and was filed as an elimination.
The ladder was not varying anything (P-51): it assigned `turbulence.K_MAX` at
runtime, and the value the mixing scheme uses is bound into
`vertical_mixing`'s signature at import.

**The ladder, re-run through the constructor.** 4000 m terrain, 8-level
sponge, clean and filtered, 12-hour ceiling (2026-09-11/12):

| K_MAX (m²/s) | 100 | 110 | 125 | 150 | 200 | 250 | 300 | 1000 |
|---|---|---|---|---|---|---|---|---|
| survived | 6/12 | 6/12 | **7/12** | **8/12** | 8/12 | 8/12 | 8/12 | 8/12 |

A monotone ramp between 100 and 150, and flat above it. **The earlier note
that it "saturates somewhere between 300 and 1000" was an artifact of having
only three rungs**; the whole effect is bought by the first 50 m²/s.

**Dissipation or suppression — the part that decides whether the hours are
worth anything.** A scheme can buy survival by dissipating what kills the run
or by flattening the flow until nothing is left to break, and a survival count
cannot tell them apart (L2; the sponge failed exactly this way in P-16 and
P-49). Four predictions were written into `kmax_binding.py` and committed
before any run finished. Compared at a COMMON hour (h6), not at each run's own
last hour:

| K_MAX | max\|u\| | jet | overturning | N² mid | clip% |
|---|---|---|---|---|---|
| 100 | 54.8 | 48.1 | 0.369% | 2.214e-04 | 0.11 |
| 150 | 54.7 | 48.4 | 0.375% | 2.214e-04 | 0.02 |
| 200 | 54.7 | 48.5 | 0.375% | 2.214e-04 | 0.01 |
| 300 | 54.7 | 48.5 | 0.379% | 2.214e-04 | 0.00 |
| 1000 | 54.7 | 48.5 | 0.378% | 2.214e-04 | 0.00 |

- **P2 held, and it was the one written to be able to fail.** max|u| does not
  fall as the ceiling rises — 54.7 at every setting — and the jet is if
  anything marginally *stronger* with more mixing available, 48.1 → 48.5.
  Nothing is being flattened. This is not suppression.
- **P4 held exactly.** Mid-level stratification is 2.214e-04 at every setting,
  to four significant figures. The column is not being mixed out.
- **P3 FAILED.** The prediction was that more available mixing would reduce
  the overturning fraction, since it is available precisely where Ri ≤ 0. It
  does not: 0.369% → 0.379%, very slightly the wrong way, and the convective
  adjustment fires at the same rate at every setting. **So the two hours are
  not bought by suppressing the overturning, and what they ARE bought by is
  not established.** See "what is not explained" below.
- **P5 held.** At K_MAX = 1000 the realized diffusivity peaks at 605 m²/s and
  the clip fraction is 0.00% throughout, so above roughly 600 the parameter is
  inert by construction. That explains 1000 ≡ 300. It does not explain
  150 ≡ 300, where the ceiling still binds.

**What the clip column says.** At 100 the ceiling truncates the diffusivity
the scheme itself asked for on 0.11% of interfaces; at 150 that falls to
0.02%, and above ~600 the formula never asks for more. A cap biting on about
one interface in a thousand, in the breaking region, was worth two forecast
hours.

**Fix.** The default ceiling is raised from 100 to 200 m²/s. 200 rather than
150 for margin: the ramp is complete by 150, 200 measures identically to 150
and 300 on every discriminator, and observed diffusivities in a breaking
mountain wave are 10²–10³ m²/s, so 100 was low on physical grounds as well.

**Confirmed by.** The ladder above, and by the production case being
untouched — 2500 m terrain with a 5-level sponge, which is the terrain this
project is actually for (P-01):

| K_MAX | survived | max\|u\| | jet | min Ri | N² mid |
|---|---|---|---|---|---|
| 100 | 12/12 | 44.0 | 23.5 | 0.085 | 2.013e-04 |
| 150 | 12/12 | 43.9 | 23.5 | 0.068 | 2.013e-04 |
| 300 | 12/12 | 43.9 | 23.5 | 0.063 | 2.013e-04 |

Suites re-run at the new default 2026-09-12, all green: `test_primitive_sigma` 7/7, `test_surface` 6/6, `test_initialization` 5/5, `test_convection` 5/5, `test_radiation` 7/7, `test_subgrid` 7/7, `test_sigma` 7/7, `test_boundaries` 6/6, `test_shallow_water` 8/8.

**What is NOT explained, and is the honest residue of this entry.** P3's
failure means the mechanism is open. More available diffusivity does not
reduce overturning, does not change the stratification, does not change the
wind, and yet moves survival by two hours. The most likely remaining
explanation is that the extra diffusivity acts on MOMENTUM in the breaking
layer — removing the shear that would otherwise concentrate into the runaway
that ends the run — rather than on the buoyancy the overturning fraction
measures. That is a hypothesis and it has not been tested. Testing it means
looking at the momentum budget in the breaking layer, not at another ladder.

**What this does not fix.** 4000 m still fails, at hour 9 rather than hour 7,
and the failure looks the same at every ceiling. P-01 stands: Nh/U ≈ 1 is
where this model's physics runs out, and the ceiling moved the boundary
without removing it.

---

## P-51 — A sensitivity ladder set a module global that nothing read
**Category** D · **Status** FIXED · **Fixed** 2026-09-10

**Symptom.** A K_MAX ladder returned peak |v| of 44.1 at both 100 and 400 —
identical to one decimal place. Earlier, the same mechanism produced a
recorded negative result: P-40's "K_MAX 100 / 300 / 1000 → 6/12, 6/12, 6/12".

**Diagnosis.** Python binds default arguments once, when the `def` is
executed:

    def vertical_mixing(..., k_max=K_MAX, ...):

`K_MAX` is read at import and frozen into the signature. `primitive_sigma`
then called `vertical_mixing(u, v, theta, pi, lev)` with no `k_max`, so the
model always used the value from the moment `turbulence` was first imported.
`kmax_ladder.py` assigned `turbulence.K_MAX = kmax` on each pass, which
changes the module attribute and nothing else. Every rung of the ladder ran
the same experiment.

Measured directly, 2026-09-10, on a neutral column with 12 m/s of shear per
level:

| how the ceiling was set | max K |
|---|---|
| `turbulence.K_MAX = 400`, default call | 100.0 |
| `k_max=400` passed explicitly | 400.0 |

`ri_crit` and `mixing_length` had exactly the same exposure and were fixed
with it, though no experiment is known to have varied them this way.

**Fix.** The three mixing parameters are now instance state on
`PrimitiveSigma` (`k_max`, `ri_crit`, `mixing_length`), defaulted from the
module values at construction and passed explicitly into `vertical_mixing` on
every call. `lid_test.build_on` forwards `**model_kw` to the constructor, and
`kmax_ladder.py` varies the ceiling that way instead of by assignment. The
comment at the point where someone would reach for the global says what
happens if they do, and quotes P-40.

**Confirmed by.** `test_primitive_sigma.test_mixing_knobs_are_connected`: max
K is 100.0 at `k_max=100` and 400.0 at `k_max=400`, and the constructed
default equals `turbulence.K_MAX`. It is a guard rail, not a physics test —
its whole job is to fail if the path is ever re-frozen, rather than let a
ladder return the same number three times and be believed. Sigma core suite
7/7 with it added.

**Why category D and not E.** The experiment was designed correctly; the
language semantics silently defeated it. What makes it expensive is that the
failure mode of this defect is *plausible output* — three identical numbers
look like a flat sensitivity, which is a publishable-shaped result.

---

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
| P-41 | **lid height** | 200 hPa 11/12, 100 hPa 10–11/12, 50 hPa 9–10/12 — neutral to worse |
| P-46 | **initialization shock over terrain** | the conversion is accurate: geopotential error 0.01 m flat, 4.19 m over 2500 m terrain. The "shock" was a rest-start adjustment plus an inconsistent test analysis |

**P-40 has been WITHDRAWN from this table and reopened.** It read
"K_MAX 100 / 300 / 1000 → 6/12, 6/12, 6/12", and three byte-identical survival
counts are not an elimination, they are a broken experiment: the ladder set a
module global that nothing read (P-51). Re-run through the constructor it is
6/12, **8/12**, **8/12**. Two forecast hours were available the whole time and
the register said the candidate was dead.

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

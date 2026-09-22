# Learning log — lessons, and where they came back

The research log records what happened on a day. The problem register records
what was wrong. Neither answers the question this file exists for: **does a
lesson learned once get applied again?**

That is the AI-collaboration study's most interesting measurable. Anyone can
learn a lesson in the moment. The evidence that it was actually learned is a
later episode where it fires *without being prompted*, on a different problem,
in a different part of the code.

## How to read an entry

Each lesson has an **origin** — where it was paid for, usually expensively —
and a list of **callbacks**: later moments where it was applied. A callback is
only listed if it changed what happened next. Recognising a lesson in
hindsight is not a callback.

Callbacks are marked:

| mark | meaning |
|---|---|
| **H** | the human prompted it |
| **A** | the AI applied it unprompted |
| **T** | a tool or test enforced it, with no one remembering |

The **T** category matters most in the long run. A lesson that survives only
in someone's memory decays; one wired into a suite or a script does not.

---

## L1 — Probe the error, don't guess-check

**Origin.** 2026-09-02, prompt 53: *"I think taking a step back and probing the
error is a better idea then guess checking."* Said after nine single-candidate
patches, all measured, all negative. It is the most consequential instruction
in the project and it is eighteen words long.

**Callbacks.**

- **A** · 2026-09-02 — probed the noisy-state failure instead of patching it:
  step-by-step recording located a 2Δx mode near the boundary, hyperdiffusion
  at a 3 h e-folding against a mode doubling in 20 min. Three defects found in
  one session, all in the test setup.
- **A** · 2026-09-03 — the tall-terrain failure: rather than adding damping,
  moved the sponge base and watched whether the growth peak followed it (it
  did: 0 → 5 → 8 → 18).
- **A** · 2026-09-03 — hour-by-hour watch of the 2500 m run, which showed the
  wind flat at 41 m/s and the Richardson number going negative at hour 10. The
  failure was mountain-wave overturning, invisible in any summary statistic.
- **A** · 2026-09-04 — P-46 opened, three hypotheses tested in order, two
  eliminated by measurement, the third turned out to be the test itself.
- **A** · 2026-09-06 — **the clearest one.** After two failed adjustments to
  the radiative boundary, stopped tuning and wrote the feedback loop down:

      F -> dpi/dt -> pi -> phi_top -> F

  with a gain of 154·|k| per second — a 130-second e-folding. That number
  explained the observed growth and nothing else did, and no setting could
  have fixed it, because the sign that radiates waves correctly is the sign
  that makes the loop grow. Two requirements with opposite signs is not a
  tuning problem, and only the arithmetic showed that.

**Transfer verdict.** Prompted once, applied five times unprompted, including
on a problem class (a control-loop instability) unlike the one it was learned
on. This is the strongest transfer in the project.

---

## L2 — A fix must predict the outcome it was proposed to explain

**Origin.** 2026-09-03. Deepening the sponge halved the reflection amplitude,
36.4 → 21.3 m/s, visibly and measurably. It changed survival by nothing:
11/12 either way. A patch-and-check loop would have accepted it as the answer.

**Callbacks.**

- **A** · 2026-09-03 — convective adjustment: four predictions written before
  the run, including one that could only fail — *the wind must NOT be damped*.
  A scheme buying stability by flattening the flow would have looked identical
  on the other three, and that is exactly how the first sponge failed (P-16).
- **T** · 2026-09-04 — `tools/problem.py check` flags any FIXED entry with no
  "Confirmed by" measurement. On its first run it found **thirteen** fixes
  asserted with no number attached. All were real; writing the numbers out
  forced them to be found again.
- **A** · 2026-09-05 — P-49: required each candidate sponge fix to move the
  survival count, not just the reflection. Six settings, none passed.
- **A** · 2026-09-10 — reopening P-40 (below), the "what remains" list leads
  with the prediction that must be written before the re-run: if a higher
  ceiling buys hours by flattening the jet rather than by dissipating the
  breaking wave, max|u| falls, and the survival count alone cannot tell those
  apart.
- **A** · 2026-09-12 — **the lesson paid out, and not in the way expected.**
  Four predictions for the K_MAX ladder were written into the script's
  docstring and committed to git *while the runs were still going*, so they
  could not drift toward the data. P2 (max|u| must not fall) held and licensed
  the change. **P3 failed** — more available mixing does not reduce the
  overturning fraction — and that failure is the only thing the session
  learned: it turned "the ceiling dissipates the overturning" from a
  conclusion into an open question, and stopped a mechanism being narrated
  into the register alongside a real measurement. The predictions that held
  permitted a change; the one that failed produced knowledge.

---

## L3 — Suspect the test before the model

**Origin.** Paid for four times. A `tanh` jet on a periodic domain. An eddy
energy ratio divided by ~0, reported as "×1.6e29". A test jet at Ro ≈ 7 read
as a flat-ground instability. An Ekman angle measured against the wind aloft,
which came out backwards.

**Callbacks.**

- **A** · 2026-09-02 — the decisive test was integrating a **166 m/s** jet
  clipped at 60, destroying geostrophic balance over 34% of the domain. The
  clip was the grid-scale source the model had been blamed for.
- **A** · 2026-09-04 — P-46's third hypothesis was "the test", and it was
  right: the state had a temperature gradient and no wind, so it *had* to
  accelerate. On flat ground the same setup drifted 8.98 m/s of the 11.54.
- **A** · 2026-09-05 — caught a metric artifact before believing it: at
  sponge = 2 the growth ratio collapsed while max|v| was unchanged at 12.8.
  Those cannot both mean suppression. The full curve settled it.
- **A** · 2026-09-05 — caught it again in the opposite direction: the first
  radiative test asserted the disturbance **at the lid** should be smaller, and
  it was larger. That is correct behaviour — a transparent lid lets the wave
  through, so the top level is more active, not less.

- **A** · 2026-09-12 — a new shape of the same lesson: not a test that reads
  the model wrongly, but an instrument that **changes what it measures**. A
  stall guard was written to integrate each forecast hour in ten-minute chunks
  so it could check the clock between them. It works, and it is not neutral:
  `run()` truncates its final step to land exactly on the requested duration,
  so chunking changes the step sequence, and at K_MAX = 110 the chunked and
  whole-hour runs diverged (min Ri 0.012 against 0.022 at hour 6) and then
  failed differently. Caught because the number moved when only the harness
  had changed. Replaced with a callback that reads the clock and touches
  nothing, then checked against the un-chunked run before being trusted.

**Transfer verdict.** Now the third hypothesis by habit rather than the last
resort. Four of the five callbacks are the AI catching its own test, and the
newest one is the harder case: the test was not misreading the model, it was
quietly running a different model.

---

## L12 — Identical results at different settings are a bug report

**Origin.** 2026-09-08. A K_MAX ladder returned peak |v| of 44.1 at both 100
and 400 — identical to one decimal place. The knob was not connected:
`vertical_mixing(..., k_max=K_MAX)` binds the module global when the module is
imported, so assigning `turbulence.K_MAX` at runtime does nothing.

The same bug had already produced a **recorded negative result**: P-40, "K_MAX
100 / 300 / 1000 → 6/12, 6/12, 6/12". Three settings giving byte-identical
survival counts, read as a clean elimination rather than as a broken
experiment. Re-run properly it gives 6/12, **8/12**, **8/12** — two forecast
hours that were available the whole time.

**Callbacks.**

- **T** · 2026-09-10 — the lesson stopped being a habit and became a test.
  `test_primitive_sigma.test_mixing_knobs_are_connected` puts a neutral,
  strongly sheared column through the model at two ceilings and asserts the
  peak diffusivity is 100.0 at `k_max=100` and 400.0 at `k_max=400`. It
  measures no physics at all. Its only job is to fail the day the path is
  re-frozen, instead of letting a ladder return the same number three times
  and be believed.

**Where it lives.** Instance state instead of module globals (P-51), with the
reason written at the point where someone would otherwise reach for the
global, and the guard rail above behind it. That is a **T**, but a narrow one:
it protects this knob. The general form — *identical output at different
settings is a bug report* — is still only a habit, and the register now says
so in P-40 rather than leaving it to memory.

**A note on how this one was found.** Not by review, and not by a test. By
noticing that two numbers agreed too well. That is worth saying because it is
the only defect in the project detected by a result being *too clean*, and
nothing in the toolchain looks for that.

---

## L4 — Check against a known answer, not against yourself

**Origin.** 2026-09-04. The height-bracket search in
`surface_pressure_from_heights` had the ordering backwards — after sorting by
descending pressure, index 0 is the *lowest* height. Every column above the
lowest analysis level stayed pinned at that level's pressure: a **253 hPa**
error over 2500 m terrain. A self-consistency check passes this. Comparing
against a standard atmosphere, where the answer is known in closed form, does
not.

**Callbacks.**

- **A** · 2026-09-04 — when the conversion looked wrong, built a synthetic
  analysis whose heights are the hydrostatic integral of its own temperatures,
  rather than trusting the one already written. The geopotential error fell
  from 140 m to 3.24 m: the original test data had not been in hydrostatic
  balance with itself.
- **T** · standing — `test_interpolate.py` and `test_sigma_operator.py` both
  compare against closed-form atmospheres rather than against the model.

---

## L5 — A diffusion loses a race against a growing mode

**Origin.** 2026-09-02. Hyperdiffusion tuned for a 3-hour e-folding at the
grid scale, against white noise amplified by advection with a 20-minute
doubling time. Damping lost by a factor of nine. The fix was not more
damping — it was removing the variance before the run started.

**Callbacks.**

- **A** · 2026-09-03 — **pattern recognised, not rediscovered.** The mountain
  wave overturns; `eddy_diffusivity` does treat Ri ≤ 0 as full-strength mixing,
  but it is a diffusion capped at 100 m²/s, relaxing a 600 m layer in
  dz²/K = 3600 s, and the wave steepens faster than an hour. Same shape, a
  different scheme. Recognising it pointed straight at convective *adjustment*
  — a rearrangement, not a faster diffusion — which took 2500 m terrain from
  11/12 to 12/12 and out to 16 hours.

**Transfer verdict.** The single clearest case of a lesson generalising across
schemes — with one correction now attached to it. The reason to prefer
rearrangement over a faster diffusion stands. The evidence cited alongside it,
that raising K_MAX "does nothing (6/12 at 100, 300 and 1000)", was P-40, and
P-40 was a broken experiment. Re-run, the ceiling is worth two hours. **The
lesson was right and one of its supporting numbers was not**, which is its own
small warning: a good pattern attracts confirming evidence and does not check
it.

---

## L6 — Order is a measurement, not a convention

**Origin.** 2026-09-02. Filtering and rebalancing an initial state:

| treatment | initial max\|div\| | survived |
|---|---|---|
| none | 3.90e-05 1/s | 1/12 |
| filter only | 9.93e-05 1/s | 11/12 |
| filter, then rebalance | 1.23e-05 1/s | **12/12** |

Filtering *raises* divergence and survives ten hours longer. Divergence was
not the controlling variable; wavenumber content was.

**Callbacks.**

- **A** · 2026-09-03 — convective adjustment applied as a **post-step**
  operator, not inside the Runge-Kutta stages: an intermediate stage would
  re-create the instability the final state is meant to be free of.
- **A** · 2026-09-04 — the pressure-to-sigma conversion is three ordered steps,
  and the driver applies filter-then-rebalance identically to the initial state
  and every boundary frame, so the edges are never prepared differently from
  the interior.
- **A** · 2026-09-04 — `verify.py` writes raw observations **before** parsing,
  QC or matching is attempted, because everything downstream is recomputable
  and the raw payload is not.

---

## L7 — A constraint can expire silently

**Origin.** 2026-09-05. "A deep sponge flattens the jet" was true when
recorded (P-16, a sponge that relaxed toward the horizontal mean). The sponge
was then changed to relax toward a frozen reference, which for a steady jet
*is* the jet. The constraint was false from that moment and was carried for
four more days, ruling out an experiment that turned out to be free.

**Callbacks.**

- **A** · 2026-09-03 — re-measured the lid-height decision on corrected initial
  states, because the original measurement had been taken on a state now known
  to carry a clipped 166 m/s jet. The conclusion survived; the point is that it
  was re-run rather than cited.
- **T** · standing — measurements are now written into the code beside the
  setting they justify, with dates, so the next person to reach for a knob can
  see what has already been tried and when.
- **T** · 2026-09-10 — `tools/stale.py` flags a dated measurement whose file
  has been modified well after the number was taken: the code moved and
  nobody re-ran it. It cannot know whether a claim is still true; it can know
  that nothing has re-checked it. A false positive costs one re-run; the false
  negative it is built against cost four days.

**What building it found, and a correction.** The interesting output is not
the staleness list. It is `--undated`. On 2026-09-10, against this repository:

| | count |
|---|---|
| measurements carrying a date | 432 |
| measurements with no date within twelve lines | 476 |
| dated measurements in files changed >14 days later (git dates) | **0** |
| the same, using filesystem mtimes on a fresh clone | 12 |

A number with no date cannot expire, because nothing can tell when the code
moved past it — so slightly more than half of what this project reports is
outside the reach of any staleness check, and dating those is what would make
the check worth anything.

The mtime row is the tool's honest weakness and is why it prefers git: a fresh
copy stamps every file at once, so everything looks as though it changed
today. All twelve of those hits were files nobody had touched.

**And the counts above already moved once, on 2026-09-12**, which is the joke
this lesson keeps making at its own expense. Appending a research-log entry
that day updated `RESEARCH_LOG.md`'s file date, and the checker promptly
flagged twelve August measurements inside it as stale — every one of them an
entry that is frozen by design and that nobody has any business re-running. A
check with a 100% false-positive rate gets switched off, so append-only
records are now out of scope, and the counts it reports are:

| | 2026-09-10 | 2026-09-12 |
|---|---|---|
| dated measurements | 432 | 169 |
| undated | 476 | 466 |
| stale by git dates | 0 | 0 |

The drop in "dated" is entirely the research and prompt logs leaving scope,
not measurements disappearing. Recording both columns rather than overwriting
the first, because a number that changed when the *instrument* changed and not
the code is exactly the thing this file exists to keep visible.

**The correction.** A note written on 2026-09-08 recorded these counts as
"5 dated, 19 undated" and drew the conclusion that the project's measurements
were almost entirely undated. Those figures came from an earlier version of
the checker and **do not reproduce** — the tool as committed reads comments,
docstrings and prose rather than every line, and treats a dated research-log
heading as dating its section. The direction of the finding survives (there
are more undated numbers than dated ones); the numbers did not. A measurement
about the project's measurements, recorded without the script that produced
it, went stale in two days. There is no better demonstration of L7 available.

---

## L8 — A metric has a blind spot; look at the curve

**Origin.** 2026-09-05. Eddy energy on day 2 over day 1 cannot tell a wave
that never grew from one that grew fast and saturated. At sponge = 2 the ratio
said 0.93 while max|v| said nothing had changed. The full time series settled
it — monotonic decay with any sponge, growth without — and also showed the
no-sponge case oscillating by a factor of two between samples, meaning the
ratio had been partly sampling luck for four days.

**Callbacks.**

- **A** · 2026-09-05 — the radiative boundary's development test first
  compared 24 h to 12 h and reported 0.61 against 0.73: no discrimination,
  because both configurations are still shedding the initial transient there.
  Moved to 36 h over 18 h **chosen from the measured curve**, which separates
  cleanly (0.75 against 3.12). An assertion window picked for speed rather
  than from the data is not a test.
- **A** · 2026-09-12 — the survival count has a blind spot of its own, and it
  cost fifteen minutes of wall clock before anyone looked. A K_MAX = 110 run
  sat in forecast hour 7 for **911 seconds**: not slow physics, but a wind of
  **7176 m/s** with the timestep collapsed from 14.8 s to 0.87 s. Every sweep
  in this project tests for failure at the hour boundary, so a run that is
  blowing up and still finite keeps integrating, and gets slower as it does.
  "Reached hour 8" and "reached hour 8 at a millisecond timestep" are not the
  same forecast and nothing reported the difference. Now P-52, with dt, steps
  per hour and wall clock recorded beside the survival count.

---

## L9 — An offline suite cannot see an interface

**Origin.** 2026-09-01 and before. Five interface defects — a GRIB regex
matching 0 of 708 messages, cfgrib's CF renaming, Herbie writing to a path
that did not exist, a module shadowing the standard library, a duplicated
`src` — every one of them passed a full offline suite and appeared on first
contact with the real service.

**Callbacks.**

- **A** · 2026-09-04 — found P-48 by *reading the call rather than running
  it*: `verify.py` was built around storing the raw payload verbatim, but
  called `fetch_asos`, which parses internally and returns objects. Every
  offline test passed because they all inject a saved payload and never call
  the fetcher.
- **A** · 2026-09-04 — the surface-field GRIB search was written and
  immediately labelled untested against the live service, in the code and in
  the register, rather than described as done.
- **T** · standing — P-06 exists specifically to hold the list of code paths
  that have never touched the network.

---

## L10 — Protect the irreplaceable thing first

**Origin.** 2026-09-04, designing the verification archive. Observations stay
downloadable for years; the forecast that was valid for them is only makeable
on the day. So raw payloads are written verbatim and compressed before
anything that could throw, with the forecast copied beside them.

**Callbacks.**

- **A** · 2026-09-04 — `tools/pull.sh` and every sync since exclude `data/`
  explicitly, and say so in the code, so an update can never overwrite the
  archive.
- **A** · 2026-09-04 — `tools/daily.sh` attempts verification **even when the
  forecast step failed**, because a run that diverged at hour 8 still produced
  eight hours worth archiving.

---

## L11 — Guard rails catch what self-review does not

**Origin.** Continuous, but sharpest on 2026-09-05: the sponge default was cut
from 5 to 3 on a reading contradicted by the table written in the same commit
— three levels decays *faster* than five. Self-review passed it. The
regression suite did not: the decisive noisy case went 12/12 → 11/12 and the
Ekman test broke. Reverted the same day, both suites recovered.

**Callbacks.**

- **T** · 2026-09-04 — `tools/problem.py check` found 20 gaps in the first
  draft of the register.
- **T** · 2026-09-04 — `tools/manifest.py --check` answers the question
  `pull.sh` and `checklayout.py` cannot: is this copy byte-for-byte what it
  should be. It is what found that thirty files had never been pushed.
- **T** · 2026-09-04 — `tools/checklayout.py` detects the `src/src` and
  nested-package shapes that a human caught once and no test ever had.
- **T** · 2026-09-10 — the K_MAX guard rail (L12). Note what it is *not*: it
  does not check that the mixing scheme is right, only that its knob is
  attached. Guard rails of this kind are cheap and specific, and this project
  has now been burned twice by a parameter that was not doing anything.

---

## Scorecard

| lesson | origin | callbacks | H | A | T |
|---|---|---|---|---|---|
| L1 probe, don't guess-check | human, prompt 53 | 5 | 1 | 5 | 0 |
| L2 a fix must predict | measurement | 5 | 0 | 4 | 1 |
| L3 suspect the test | four test-design defects | 5 | 0 | 5 | 0 |
| L4 check against a known answer | a 253 hPa bug | 2 | 0 | 1 | 1 |
| L5 diffusion loses a race | hyperdiffusion vs noise | 1 | 0 | 1 | 0 |
| L6 order is a measurement | filter/balance | 3 | 0 | 3 | 0 |
| L7 constraints expire | a 4-day stale belief | 3 | 0 | 1 | 2 |
| L8 metrics have blind spots | a saturation artifact | 2 | 0 | 2 | 0 |
| L9 offline suites miss interfaces | five interface defects | 3 | 0 | 2 | 1 |
| L10 protect the irreplaceable | archive design | 2 | 0 | 2 | 0 |
| L11 guard rails beat review | a default changed wrongly | 4 | 0 | 0 | 4 |
| L12 identical results = a bug | a frozen default | 1 | 0 | 0 | 1 |

**What the counts say, with the caveat that they are self-reported and n = 1.**
Every lesson but one originated in a *measurement that surprised someone* —
not in reasoning, and not in a code review. Only L1 came from the human
directly, and it is also the one that transferred furthest.

**Where the lessons live now** is the more useful column. Five of twelve have
no enforcement behind them: they exist as habits, and habits are exactly what
a fresh session does not have. The seven with a **T** are the ones that will
still be working in a month.

**Two of the entries above now correct themselves** (L5's supporting number,
L7's own counts), and both corrections came from re-running something rather
than re-reading it. That is the pattern this file is for.

**2026-09-12 adds a third kind of entry**: a prediction that failed on
purpose. L2 has always said a fix must predict the outcome it was proposed to
explain. What the K_MAX session shows is the sharper version — the prediction
that FAILS is the one that pays. Two predictions held and merely licensed a
default change; the third failed and is the reason the register now says the
mechanism is unknown instead of asserting one. A session where every
prediction holds has probably not asked anything.

---

## Maintaining this file

Add a callback when a lesson visibly changes what happens next — not when it
is merely recalled. If a lesson fires three times unprompted, consider whether
it can be moved from habit to a **T**: a test, a check in `tools/problem.py`,
or a comment written at the point where someone would otherwise do the wrong
thing.

Numbers in this file are measurements like any other, and `tools/stale.py`
reads this file too. Date them.

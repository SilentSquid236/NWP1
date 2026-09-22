---
name: nwp-debug
description: Use when an NWP1 model run fails, diverges, gives a surprising number, or a proposed fix needs testing. The project's probe-first method for diagnosing numerical and physical failures before changing code.
---

# Diagnosing a failure in NWP1

This procedure exists because nine single-candidate patches were once applied
to a model that was not broken — the test setup was. Every step below was
learned from a specific failure; the origins are in `docs/LEARNING_LOG.md`.

## 1. Before touching code, suspect the setup

Check, in this order, and write down what you find:

- Is the initial state **balanced**? A temperature gradient with no wind must
  accelerate — that is geostrophic adjustment, not a bug. Compare against a
  flat-ground control before blaming terrain.
- Is the test jet **realistic**? Compute Ro = U/(fL). A 1.5 K meridional
  contrast gives ~41 m/s; 6 K gives 166 m/s.
- Is any synthetic analysis **hydrostatically self-consistent**? Heights must
  be the integral of the temperatures, not perturbed independently.
- Is the thing being varied **actually varying**? Identical results at
  different settings mean a disconnected knob. Module constants bound as
  default arguments do not change when the module global does — pass
  `k_max`, `ri_crit` through `PrimitiveSigma(...)`.

## 2. Probe: locate the failure

Record every prognostic field, per level, on a short interval, and locate:

- **where in space** — edges vs interior (an FFT on a non-periodic domain puts
  its error at the boundary)
- **which level** — index 0 is the lid, index -1 the ground
- **which scale** — dominant wavelength; 2Δx means numerics, domain-scale means
  physics or imbalance
- **which field first** — the first non-finite value, and its neighbours one
  step earlier
- **which quantity runs away** — often not the one being watched (P-50's was
  the meridional wind)

If nothing recorded runs away before the failure, the failing quantity is not
being recorded. Widen the net.

## 3. After two failed adjustments, stop tuning

Write the mechanism down. If it is a loop, write the loop and estimate its gain
from the numbers you have. A gain with an e-folding of seconds explains a
failure no parameter setting can fix — and two requirements with opposite signs
is not a tuning problem.

## 4. Predict before running

In `docs/RESEARCH_LOG.md`, before the run, state what the fix should do —
including **at least one prediction that can only fail** if the fix is working
for the wrong reason. For damping-type fixes that is always: *the balanced flow
must not be weakened.* Commit the predictions before the results exist.

## 5. Measure the right thing

- A fix must move the outcome it was proposed for, not just a proxy.
- Look at the **whole curve**, not two points: a ratio cannot tell saturation
  from suppression.
- Pick assertion windows from the measured curve, not for speed.
- Compare against a **known closed-form answer** (standard atmosphere,
  analytic balance), not against the model's own output.

## 6. Respect the server

50% of cores maximum; no installs; long runs go in the background with
unbuffered output (`python -u`) and progress lines, because a silent run is
indistinguishable from a hung one. `kill -USR1 <pid>` dumps a traceback
without stopping it.

## 7. Record it

Every elimination is a result. Add it to the ELIMINATED table in
`docs/PROBLEMS.md` with the measurement that eliminated it. Date every number.

# Using an AI Collaborator to Build a Weather Model

A study document. The research question is not whether the model works — it
is what an AI collaborator gets right and wrong when building scientific
software, and what catches the errors.

This project is the primary source. The AI (Claude) wrote essentially all the
code; the human set direction, ran everything on real hardware, and supplied
the ground truth. Both roles matter to the result.

## Methodological caveat, stated first

**n = 1, uncontrolled, and partly self-reported.** The AI is both subject and
observer here, which is a real weakness: errors it never noticed cannot appear
in its own tally. The bug counts below are therefore a *lower bound*.

Mitigations used: every claim is tied to a test or a measurement in the log;
failures are recorded with the same weight as successes; and the human ran the
code independently, which is how several errors surfaced at all.

Anyone extending this should treat the taxonomy as a hypothesis-generator, not
a result.

## What was built

A dry hydrostatic limited-area weather model, from nothing, in roughly four
days of intermittent work: C-grid dynamical core (shallow water, then 3D
primitive equations), Davies boundary relaxation, subgrid dissipation and
stochastic perturbation, HRRR ingestion, an observation-based verification
harness, and adaptive bias correction. ~70 tests.

It is not a usable forecast system — see `docs/CAPABILITIES.md` and
`docs/STABILITY.md`. That is part of the finding.

## Error taxonomy

Every defect the AI introduced, classified. Detail and numbers in
`docs/RESEARCH_LOG.md`.

### A. External-interface assumptions — 5 errors, the largest category

Conventions of systems the AI could not observe directly.

| error | consequence |
|---|---|
| GRIB search regex anchored with `^` | matched 0 of 708 messages; downloaded nothing |
| cfgrib renames variables to CF short names | `KeyError` on every variable |
| Herbie's cache defaults to `~/data` | silent write failure → misleading `FileNotFoundError` |
| `operator.py` shadowed the stdlib | broke `collections`, then `numpy` |
| assumed home-directory quota was the cause | wrong hypothesis; wasted a cycle |

**Every one survived a full offline test suite.** None could have been caught
without either real data or real environment inspection.

### B. Discrete-vs-continuous mathematics — 2 errors

Correct continuum reasoning applied to a discrete operator.

- Hyperdiffusion coefficient derived from continuous `k⁴` rather than the
  discrete Laplacian's 2Δx response: damping **6× weaker** than intended.
- Flux-form theta transport assumed exact discrete continuity, which the
  omega boundary correction breaks: unstable, reverted.

### C. Stability and dimensional analysis — 2 errors

- Divergence damping written as a tendency with a coefficient in m²/s,
  violating `nu·dt/dx² ≤ 0.25`. Blew up the model it was meant to stabilise.
- Sponge layer relaxing wind toward the horizontal mean — absorbs waves, but
  also flattens a jet, which is legitimate structure.

### D. Array and language semantics — 2 errors

- `axis=1` meaning x in 2D but y in 3D: the core differenced north–south when
  it meant east–west.
- `parse_raob_csv` shadowed its own parameter, silently ignoring caller input.

### E. Test-design errors — 3 errors

Failures that looked like model bugs but were bugs in the question:

- Test configurations violating their own CFL limit.
- A `tanh` jet used as a balanced state on a *periodic* domain.
- Eddy-energy growth measured as a ratio to a near-zero baseline (×1.6e29).

### F. Wrong causal hypotheses — 3

Stated confidently, then falsified by measurement:

- Vector-invariant momentum would fix energy drift. (It was time truncation;
  the change improved enstrophy conservation 14× instead.)
- Aliasing explained the excessive initial divergence. (Block averaging
  changed nothing.)
- Herbie's home-directory cache caused the download failure. (235 GB free, no
  quota.)

## What caught the errors

| detection mechanism | errors caught | share |
|---|---|---|
| AI-written tests, run offline | 7 | ~41% |
| First contact with real data | 5 | ~29% |
| Measurement designed to test a stated hypothesis | 3 | ~18% |
| AI self-review before delivery | 1 | ~6% |
| Human noticing an anomaly (`src/src`) | 1 | ~6% |

**The strongest single tool was the analytic test** — asserting an exact
answer rather than a tolerance. Six such tests hold to machine precision
(0.00e+00, 1.08e-18, 8.9e-16), and each would break instantly on an indexing
or sign error.

**The second strongest was the convergence test.** Where no exact answer
exists, refining the grid and checking the order (3.97, 3.99 against a
theoretical 4.0) distinguishes truncation error from a bug. That distinction
is not otherwise available.

## Observations

**1. The AI tested everything it could test, and that was not enough.** Every
offline-testable component worked on first contact with real data. All five
category-A failures were at interfaces with systems whose behaviour had to be
assumed. For scientific software, the boundary with external data formats is
where AI-written code fails, not the algorithms.

**2. Predictions were wrong roughly as often as they were right** on
mechanism (category F). The AI reliably produced plausible causal stories; the
stories were unreliable. What made this tractable was measuring the *predicted
effect specifically* rather than observing that things improved.

**3. Domain knowledge was mostly correct; discretisation was where it slipped.**
The AI knew the primitive equations, C-grid staggering, Davies relaxation,
Sadourny's scheme, and SPPT. Its errors clustered where continuum mathematics
meets a finite grid (category B) — knowledge that is real but shallower than
it appears.

**4. Negative results required deliberate effort to preserve.** The natural
pull is toward reporting what worked. `docs/STABILITY.md` exists because the
alternative was a reader asking "have you tried more damping?" — a question
the measurement table answers and prose would not.

**5. The human's contribution was direction and reality.** Three decisions
changed the project's course and none came from the AI: abandoning the neural
emulator for physics, insisting on observations over model output, and
requiring shared-resource discipline. The AI executed each well and proposed
none of them.

## What to measure going forward

For each work session, record in `docs/RESEARCH_LOG.md`:

- **defects introduced**, by taxonomy category above
- **detection mechanism** for each
- **stated hypotheses** and whether measurement confirmed them
- **time to first real-data contact** — category-A errors are invisible before
  it, so this is the dominant latency in the loop
- **human interventions** that changed direction, separated from those that
  corrected execution

Open questions worth designing for:

- Does the category-A share fall as the AI accumulates context about a
  specific external system, or stay flat?
- Would requiring an explicit stability analysis before writing any damping
  term have prevented both category-C errors? (Cheap to test: adopt the rule
  and count.)
- Does asking for a hypothesis *in writing before* each change reduce
  category-F errors, or only make them legible?

## Honest summary

An AI collaborator produced a verified dry dynamical core, a complete
observation-verification harness, and ~70 meaningful tests in days rather than
months. It also introduced at least 17 defects, cost several cycles on
confident-but-wrong diagnoses, and did not identify the project's most
important decision.

The resulting model does not yet forecast. Whether that is a failure depends
on what was being tested — as an exercise in AI-assisted scientific software,
the failure mode is the data.

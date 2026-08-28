# Methodology

How this project decides whether something works. Written down because the
practices below caught real bugs, and because a research project should be
able to state its standards of evidence.

## 1. Test against analytic answers, not tolerances

A tolerance chosen so the test passes measures nothing. Where an exact answer
exists, assert it:

| test | exact answer | measured |
|---|---|---|
| resting atmosphere stays at rest | 0 | 0.00e+00 |
| mass conservation, flux form | 0 | 0.00e+00 |
| hydrostatic integral, isothermal | R·T·ln(p_s/p) | 6e-12 m error |
| bilinear interpolation of a linear field | exact | 8.9e-16 |
| divergence of a rotational flow | 0 | 1.08e-18 |
| vertical interpolation, linear in log(p) | exact | 0.00e+00 |

These are not approximations that happen to be close. They are properties of
the discretisation, and any indexing or sign error breaks them immediately.

## 2. When no exact answer exists, test the convergence order

Some quantities are only correct in the limit. A balanced state is balanced
to discretisation accuracy, so a threshold on its residual is arbitrary.
Refine the grid instead:

```
geostrophic imbalance:  6.58e-04 → 1.65e-04 → 4.14e-05
ratios:                 3.97, 3.99          (2nd order = 4.0)
```

**Truncation error converges; bugs do not.** This is stronger evidence than
any single-resolution tolerance.

## 3. Write negative tests

A method that only ever helps in its own test suite is not being tested. The
bias corrector is asserted to be **neutral** when there is no bias to remove
(measured −3.6%). A method showing large gains there would be fitting noise,
and would destroy signal on real data.

Similarly: the CFL limit test asserts the model **blows up** above the limit.
A stability bound nothing violates is not a bound.

## 4. State the hypothesis before the result

Twice this project predicted a mechanism and was wrong:

- Vector-invariant momentum was expected to fix energy drift. It did not —
  the drift was time truncation. It improved enstrophy conservation 14×
  instead.
- Aliasing was expected to explain excessive initial divergence. Block
  averaging changed nothing; the cause was elsewhere.

Both were caught by measuring the predicted effect specifically rather than
observing that "things improved". A change kept for the wrong reason is a
change that will be removed for the wrong reason later.

## 5. Record what failed, with numbers

`docs/STABILITY.md` exists because a reader needs to know that damping was
tried across a grid of configurations and none worked — otherwise the obvious
first suggestion is "have you tried more damping?". The measurement table is
the answer.

Failures documented so far: flux-form theta transport (unstable, with the
mechanism), divergence damping (works, suppresses baroclinic growth, with the
numbers), stronger hyperdiffusion (no effect on survival time).

## 6. Distinguish the model's error from the harness's error

Several "failures" were bugs in the test, not the code: CFL-violating test
configs, a non-periodic initial condition on a periodic domain, an eddy-energy
ratio dividing by near-zero. Each looked like a model bug at first.

Before concluding the code is wrong, check that the test is asking a
well-posed question.

## 7. Verify against observations, never against another model

The project's central methodological claim. A model trained on — or scored
against — another model's output is bounded by that model and inherits its
biases as if they were physics.

Consequences, enforced in code:
- Verification uses ASOS, mesonets, and radiosondes only.
- `test_no_model_sources` fails if a model-derived source enters the
  observation stream.
- SPC mesoanalysis is excluded as truth: it is a gridded model analysis.
- HRRR supplies initial and boundary conditions only, which is a different
  relationship — the interior evolves under its own physics.

## 8. Interfaces with external systems are where bugs live

Every offline-testable component worked on first contact with real data.
Every failure was at a boundary with an external system whose conventions had
been assumed: GRIB index formatting, cfgrib's CF renaming, Herbie's cache
behaviour, a module name colliding with the standard library.

Practice: capture real samples of external formats and test parsers against
them offline. The GRIB regex bug is now caught by index lines copied verbatim
from a live file.

## 9. Shared-resource discipline is part of the method

The compute runs on a machine shared with ~30 researchers. CPU is capped at
50% and adapts to observed load; network is rate-limited to ~8 MB/s with
sequential requests and a local cache. Both are documented, defaulted
conservatively, and overridable.

This is not politeness for its own sake — a job that degrades everyone else's
work will be killed, and results that cannot be reproduced because the run was
terminated are not results.

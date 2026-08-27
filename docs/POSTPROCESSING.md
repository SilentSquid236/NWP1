# Post-Processing: Correcting Physics Output with Learned Patterns

Status: adaptive bias correction implemented and tested (`src/postproc/`).
Learned (NN) post-processing designed, not built — it is blocked on data that
does not exist yet.

## The idea

A physics core produces a forecast. That forecast has systematic errors: a
grid cell whose elevation does not match the station in it, a boundary layer
scheme that is always too shallow at night, a coastline the grid cannot
resolve. Physics will not remove these, because they are artefacts of
discretisation and missing processes rather than of the equations.

Post-processing learns those errors from past forecasts and removes them.
Operationally this is Model Output Statistics, and it is one of the largest
sources of forecast skill improvement in the last fifty years — often worth
more than a resolution increase.

## The one rule that makes this legitimate

**Corrections are learned against OBSERVATIONS, never against another model.**

Training a correction network toward HRRR would make this an HRRR emulator
with extra steps — the exact failure that motivated abandoning the Conv3d
approach. Every error signal must be `(our forecast − observation)`, where the
observation is a radiosonde, METAR, or other real measurement.

This has a consequence worth stating plainly: post-processing can only improve
what we verify against. Corrections exist for the stations and variables we
have observations for, and nowhere else.

## Stage 1 — adaptive bias correction (built)

`src/postproc/bias_correction.py`. A Kalman filter tracks a slowly varying
bias per `(station, lead_time, variable)`:

```
innovation = (forecast − observation) − bias_estimate
bias      <- bias + K · innovation
K          = P / (P + R)
```

The gain K adapts on its own: high when the estimate is uncertain, falling as
evidence accumulates (measured: 1.000 → 0.145 → 0.041 over 30 updates). This
is why it works on **days** of data rather than the years a network needs, and
why it tracks a bias that changes sign — seasonal drift, an instrument swap, a
model change.

Guardrails, all deliberate and tested:

- **Refuses to act** below `min_samples`. One outlier is not a bias.
- **Hard cap** on the applied correction. A pathological error stream produces
  a bounded correction, not an unbounded one.
- **Output only.** Corrections are never fed back into the model state, which
  would violate the conservation properties of the core.

### Measured behaviour

| test | result |
|---|---|
| learns a constant 2.5 K warm bias | estimated 2.50 K, RMSE −98% |
| no correction before evidence | 0.00 at 5 samples, 4.79 at 15 |
| correction capped | 100-unit error → 3.0 applied |
| **neutral when no bias exists** | **skill −3.6% (near zero)** |
| tracks a sign flip | +3.00 → −2.00 K |
| keys independent | A/f06 +4.88, B/f06 −2.93, A/f24 +0.00 |

The fourth row is the important one. With unbiased random error there is
nothing systematic to remove, and the corrector is correctly near-neutral. A
method that showed large gains there would be fitting noise — and would
destroy signal on real data.

## Stage 2 — learned post-processing (designed)

Where a network genuinely beats the Kalman filter is in **conditional** bias:
error that depends on the situation rather than being constant. The model may
be 3 K too warm on clear calm nights and unbiased in overcast wind — a single
running bias averages those into uselessness, while a network conditioned on
the right features separates them.

### Features

Everything must be available at forecast time:

- forecast value, and forecast values of related variables at the same point
- lead time, pressure level
- time of day, day of year (both as sin/cos pairs — midnight is adjacent to
  23:00, and a raw hour number does not say so)
- local static fields: elevation difference between station and grid cell,
  land/sea fraction, terrain slope
- recent verified errors at this and neighbouring stations — the "recent
  patterns" signal
- forecast stability indicators: static stability, wind speed, cloud proxy

### Target

`observation − forecast` at the station. Predicting the *correction* rather
than the value keeps the physics forecast as the baseline, so the worst case
is predicting zero and falling back to raw physics — a much safer failure mode
than predicting the field outright.

### Model

Start with gradient-boosted trees, not a neural network. With a few thousand
samples per station they beat neural networks on tabular problems, train in
seconds, and expose feature importances that tell you *why* the model is
biased — which feeds back into fixing the physics. Move to a network only when
data volume justifies it and spatial structure needs to be exploited, at which
point a small CNN over the forecast field is the natural form.

### The blocker

**There is no training data.** Not a technical gap — an archival one. Every
forecast the model produces has to be stored alongside the observations that
verify it, from the first run onward. A year of daily forecasts is a few
thousand samples per station; that is the minimum for conditional correction,
and none of it can be recovered retroactively.

So archiving must start with the first real forecast, before any of this can
be built.

## Evaluation, and how this fails

Post-processing that is not verified honestly will silently make forecasts
worse. Required, every time:

1. **Held-out data.** Never score on what was fitted. Split by time, not
   randomly — random splits leak, because adjacent hours are nearly identical.
2. **Beat the raw physics.** The skill score against uncorrected output is the
   entire point. Negative means stop.
3. **Beat the trivial baselines.** Persistence, and climatology. A corrector
   that loses to "tomorrow equals today" is not adding anything.
4. **Check the regime changes specifically.** This is where post-processing
   does damage: it drags forecasts toward recent conditions, so it is worst
   exactly when the weather changes — which is when the forecast matters most.
   Score cold-front passages and rapid changes separately from quiet days.

A corrector that improves mean scores while degrading the events people care
about is a net loss, and aggregate statistics will hide it.

## Where this sits

```
HRRR ---> initial + boundary conditions
              |
        physics core (dynamics/)
              |
        raw forecast ------------------> verification vs obs
              |                                  |
        bias correction <--- recent errors <-----+
              |
        corrected forecast
```

The verification harness feeds the corrector, so it has to be built first —
which is also stage 1 of the data assimilation plan in
`docs/DATA_ASSIMILATION.md`. One piece of work unlocks both.

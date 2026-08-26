# Data Assimilation Layer — Design

Status: design draft. No code written yet.

## Purpose

The emulator learns `state(T) -> state(T+1)` from HRRR analysis fields. HRRR is
itself a model product, so an emulator trained only on HRRR can never be better
than HRRR — it inherits every bias. The assimilation layer exists to pull the
initial state toward **real observations** before the forecast steps forward,
and to give us an independent yardstick for whether the model is any good.

Two distinct uses, often conflated:

1. **Analysis correction** — adjust the F00 initial state using obs valid at
   that time, then forecast from the corrected state.
2. **Verification** — hold obs back and score forecasts against them. This is
   the only honest measure of skill, since scoring against HRRR just measures
   agreement with the thing we copied.

Build (2) first. Without it, (1) cannot be evaluated.

## Core concepts

For a background state `x_b` (the HRRR field) and observations `y`:

```
innovation   d = y - H(x_b)          how much obs disagree with the background
analysis   x_a = x_b + K d           corrected state
```

`H` is the **observation operator**: it maps the gridded model state into
observation space, so a grid field can be compared with a point measurement.
`K` is the **gain**: how much to trust the obs relative to the background.

Full DA (3D-Var, EnKF) solves for `K` from error covariances. That is a large
undertaking. A successive-corrections scheme (Cressman or Barnes) approximates
it with a distance-weighted influence radius, needs no matrix inversion, and is
a reasonable first implementation:

```
increment(grid_point) = sum_i w_i * d_i / sum_i w_i
w_i = exp(-r_i^2 / R^2)          r = distance from ob to grid point
```

`R` is the influence radius: too small and obs affect nothing, too large and a
single station smears across the domain. Start near the correlation length
scale of the field — roughly 100 km for surface temperature, less for moisture.

## Observation operator by source

| source | type | cadence | H implementation |
|---|---|---|---|
| Radiosonde (RAOB) | point, upper air | 00Z / 12Z | bilinear horizontal + interpolate in log(pressure) to model levels |
| METAR / ASOS | point, surface | hourly | bilinear horizontal at lowest level; correct for station-vs-grid elevation |
| MRMS | already gridded | ~2 min | regrid to model grid; no scattered-point problem |
| NEXRAD Level II | volumetric, per-radar | ~5 min | needs a compositor — defer, use MRMS instead |

Vertical interpolation must be done in **log pressure**, not height or linear
pressure — atmospheric variables are far closer to linear in log(p), and the
difference is large enough to matter at the levels we carry.

Elevation correction for surface obs is not optional. A station in a valley and
its grid cell average can differ by hundreds of metres, which is several
degrees of temperature. Uncorrected, that difference is read as an innovation
and injected as a spurious correction.

## Observation record schema

One flat structure for every source, so quality control and assimilation stay
source-agnostic:

```python
{
  "time":      datetime,     # UTC, observation valid time
  "lat":       float,
  "lon":       float,
  "elevation": float,        # metres MSL, station height
  "pressure":  float | None, # hPa; None for surface obs
  "variable":  str,          # one of config.CHANNELS
  "value":     float,        # SI units, converted at ingest
  "error_std": float,        # assumed obs error, per source and variable
  "source":    str,          # "raob" | "metar" | "mrms"
  "station":   str,
}
```

`error_std` is what makes obs comparable. A radiosonde temperature is good to
roughly 0.5 K; a METAR temperature to about 1 K; both carry representativeness
error, because a point measurement is being compared with a ~3 km cell average.
That representativeness term usually dominates instrument error and should be
included in `error_std`, not ignored.

## Quality control

Bad obs are worse than no obs — a single mis-decoded value will inject a large
false increment. Minimum gates, in order:

1. **Range check** — physically impossible values (RH > 100%, T < 180 K).
2. **Gross innovation check** — reject when `|d| > 5 * error_std`. Catches
   decode errors and misplaced stations.
3. **Buddy check** — compare against nearby obs of the same variable; an ob
   disagreeing sharply with all its neighbours is suspect.
4. **Blacklist** — persistent offenders, kept in a config file.

Log every rejection with its reason. Silent dropping makes it impossible to
tell "no obs available" from "all obs rejected", and those need different fixes.

## Where it plugs in

```
Herbie -> HRRR GRIB2 -> background state x_b (F00)
                              |
obs fetch -> QC -> H(x_b) -> innovations d
                              |
                    Cressman/Barnes -> increments
                              |
                     analysis state x_a  ---> emulator ---> forecast T+1
                                                              |
                        withheld obs -----> verification -----+
```

Two training strategies, worth deciding deliberately:

- **Train on corrected states.** Assimilate first, train the emulator on
  analysis states. Simple, but obs errors get baked into the learned dynamics.
- **Train on HRRR, assimilate at inference.** Keep learned dynamics clean and
  apply corrections only when forecasting. Cleaner separation, and it lets the
  same trained model be evaluated with and without assimilation.

The second is recommended. It preserves an ablation: the same weights, scored
with and without obs, isolates exactly what assimilation contributes.

## A tension worth resolving early

**Regional domains and radiosondes conflict.** There are ~92 CONUS radiosonde
sites launching twice daily. A regional subset of a few hundred grid points may
contain **one site, or none**. Upper-air assimilation with a single sounding
twice a day is close to meaningless, and the 20-level state cannot be
meaningfully constrained by surface obs alone.

Three ways out:

1. Size the domain around radiosonde coverage — big enough for 4–6 sites.
2. Accept a surface-only assimilation regionally, and treat upper levels as
   unconstrained.
3. Use aircraft obs (AMDAR/ACARS), which are far denser in the vertical near
   airports, though access is more restricted than the public NWS feeds.

This decision should be made before the domain is fixed, because it determines
the domain.

## Staging

1. Verification harness — fetch obs, QC, score HRRR against them. Establishes
   the baseline and exercises `H` with no assimilation involved.
2. Surface assimilation — METAR/ASOS, lowest model level, Cressman increments.
3. Upper-air assimilation — radiosondes, log-p interpolation.
4. Radar — MRMS, gridded, for moisture and convective structure.

Each stage is independently testable, and each should show a measurable change
in verification scores. A stage that changes nothing is a bug, not a no-op.

## Open questions

- Time window: how far from valid time can an ob be used? ±30 min is typical
  for hourly cycling.
- Bias correction: radiosonde types have known systematic biases, some
  solar-angle dependent.
- Do we assimilate all five channels, or only those obs constrain well?
  Geopotential height is poorly observed directly.
- Normalisation interacts with this — increments must be applied in physical
  units, before whatever scaling the network sees.

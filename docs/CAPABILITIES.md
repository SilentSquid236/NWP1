# What This System Can and Cannot Do

An honest assessment. Written because a model that runs and produces
plausible-looking fields is easy to mistake for a model that forecasts.

## What works today

Given HRRR initial and boundary conditions, the system integrates the dry
hydrostatic primitive equations over the Northeast and produces a physically
plausible evolution of temperature and wind. Baroclinic waves develop and
amplify from smooth initial states.

Verified properties, all from the test suites rather than from inspection:

| property | evidence |
|---|---|
| mass conservation | exact (0.00e+00) over 12 h |
| resting atmosphere stays at rest | exact |
| gravity wave speed | 0.7% of sqrt(gH) |
| geostrophic balance | 0.00% drift over 24 h |
| thermal-wind balance | 3.9% drift over 24 h |
| numerical convergence | 2nd order (ratios 3.97, 3.99) |
| open boundaries | 0.5% reflection vs a rigid edge |
| baroclinic growth | eddies to 23 m/s, +21%/day |

61 tests across dynamics, boundaries, verification, fetchers, and
post-processing. Physics tests assert analytic answers or convergence order,
not tolerances chosen to pass.

## What is missing, and what it costs

The core is **dry, flat, hydrostatic, and unforced**. Each of those is a
first-order omission for surface weather.

| missing | consequence |
|---|---|
| **moisture / latent heat** | no clouds, no precipitation at all. Condensation drives mid-latitude cyclones; storms will be systematically too weak |
| **radiation** | no diurnal cycle. Nothing heats the ground by day or cools it at night — surface temperature goes wrong within hours |
| **boundary layer / surface fluxes** | no friction, no mixing. Low-level winds too strong, no daytime mixing, no nocturnal inversion |
| **terrain** | pressure coordinates assume a flat lower boundary. The Appalachians and the coastline both shape Northeast weather heavily |
| **nonhydrostatic dynamics** | no convection. Thunderstorms cannot form by construction |

## Stability on real data — read this first

**The core is stable on idealised balanced states and unstable on real HRRR
analyses**, diverging after 2–3 forecast hours. All 29 idealised physics tests
pass; the operators are verified to machine precision. The failure appears
only with realistic sheared, noisy initial states, and no damping
configuration fixes it — see `docs/STABILITY.md` for the full measurements.

The likely cause is the pressure vertical coordinate's ill-posed lower
boundary. The fix is sigma coordinates, already needed for terrain.

Everything below describes what the model would deliver once integrable.

## What to expect in practice

**Reasonable for 12-24 h:** mid-tropospheric flow — 500 hPa heights,
upper-level winds, the large-scale wave pattern. This is mostly dry adiabatic
dynamics, which is what the core does correctly.

**Poor:** anything near the surface. Temperature, wind, humidity at 2 m and
10 m all depend on the boundary layer and radiation, neither of which exists.

**Absent:** precipitation, cloud, convection, fog, snow. There is no water in
the model.

The fair comparison is not HRRR. It is roughly an operational model of the
early 1960s — Charney-era dry baroclinic prediction. That was a landmark, and
it also could not forecast an afternoon.

**Do not use this for forecasts anyone relies on.** A plausible-looking wind
field conceals the absence of precipitation and the diurnal cycle.

## Performance, measured

Timings for a 12-hour forecast, 20 levels, single core:

| resolution | grid | ms/step | dt | steps | 12 h forecast |
|---|---|---|---|---|---|
| 12 km | 110x97 | 138 | 46.3 s | 933 | **2.1 min** |
| 6 km | 219x194 | 708 | 23.1 s | 1867 | **22 min** |
| 3 km | 438x388 | 4040 | 11.6 s | 3734 | **4.2 hours** |

Cost scales as roughly dx^-3: halving the grid spacing quadruples the cell
count and halves the timestep.

### The 104-core box does not help

Measured: `OMP_NUM_THREADS=1` gives 596 ms/step, `=2` gives 585 ms/step. No
scaling. NumPy's elementwise arithmetic and `np.roll` are single-threaded --
only BLAS calls (matrix multiplication) use multiple cores, and this code
contains none. The dynamics is memory-bandwidth bound on one core.

The resource governor in `resources.py` still matters (it constrains anything
that *does* thread, and documents intent), but it is not buying speed here.

Making the box useful would require explicit parallelism: domain decomposition
with `multiprocessing`, or `numba`/`numexpr` on the tendency kernels. Neither
is a small change, and neither is worth doing before the physics is worth
running at scale.

### Resolution should match the physics

**Run at 6 km, or even 12 km.** The core is hydrostatic and has no convection
scheme, so 3 km resolves nothing the equations can represent -- it is false
precision at 12x the cost. The hydrostatic approximation is itself marginal at
3 km and comfortable at 6-12 km. Coarser is both faster and more defensible.

At 12 km a 12-hour forecast finishes in two minutes, which makes iteration
practical. That matters more right now than resolution.

## Where the value actually is

1. **The core is correct and proven correct.** Errors in advection, Coriolis
   staggering, or the pressure gradient become invisible once moisture and
   convection are added -- a numerical instability is indistinguishable from
   real convection. Getting this verified first is the part that cannot be
   retrofitted.
2. **The verification harness is model-independent.** It scores anything
   against real observations, including HRRR. That has standalone value.
3. **Every line is understood.** Which is the point of building rather than
   downloading WRF.

## Realistic next steps, in order

1. **Run on real HRRR data.** Ingestion and Herbie have never been exercised
   against live data.
2. **Verify against ASOS and start the archive.** Forecast-observation pairs
   are not recoverable retroactively; every day without archiving is lost
   training data for post-processing.
3. **Then choose:** boundary layer + radiation (largest gain for surface
   forecasts) or moisture (largest gain for storms). Moisture is the bigger
   job -- vapour transport, condensation, latent heating, then microphysics --
   and is realistically a months-long project.

Even with all of that, the result is a research model, not a competitive one.
Operational centres staff these with teams over years. That is worth knowing
at the outset, and it does not make the exercise less worthwhile.

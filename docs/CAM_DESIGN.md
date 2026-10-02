# Toward a convection-allowing NWP1 — design and roadmap (draft for review)

**Status.** A draft, written 2026-10-02 (prompt 163). The user set the goal
in answer to the question of whether to switch surface heating on by default:
*"no we want a convective allowing model"*. In the follow-up question they
chose *"plan toward a 3 km non-hydrostatic CAM"*. Nothing in this file has
been built yet. Every cost below is an estimate scaled from a measured run,
and stage S3 exists to replace the estimates with measurements before any
choice of grid is made.

---

## 1. What "convection-allowing" requires, and where NWP1 is

A convection-allowing model (CAM) resolves deep convective storms on the grid
rather than parameterising them. Weisman et al. (1997) found that about 4 km
grid spacing is the coarsest at which a squall line's mesoscale structure and
evolution are reproduced without a convective parameterisation. Operational
CAMs run at about 3 km: the HRRR is the US example (Benjamin et al. 2016;
Dowell et al. 2022). At these scales vertical accelerations are no longer
negligible, so a CAM needs non-hydrostatic equations.

| requirement | NWP1 today | gap |
|---|---|---|
| grid spacing ≤ 4 km | 12 km, 97 × 110 | 3–4× finer in each direction |
| non-hydrostatic dynamics (prognostic w) | hydrostatic sigma primitive equations | new dynamical core |
| no convection parameterisation | dry convective adjustment (PAV, `convection.py`) | drop it; overturning is resolved |
| moisture with explicit microphysics | dry | water vapour, cloud, rain, then ice |
| surface fluxes of heat and moisture (convection initiation) | drag only; a prescribed heating flux exists as an option, off by default | surface layer, land surface, radiation |
| a boundary layer, plus 3D turbulence at grid scale | Richardson-number vertical mixing | a PBL scheme and a Smagorinsky (1963) style 3D closure |
| storm-scale initial conditions | observation-only analysis (ASOS, raobs, buoys, MRMS radar available in `src/analysis/sources.py`) | radar and moisture in the analysis; spin-up accepted at first |

**Why surface heating matters for a CAM.** In a CAM, convection starts where
the boundary layer is warmed and moistened from below. The prescribed flux
(`diurnal.py`) was a first stand-in for that, and the user's decision keeps it
off by default. A CAM still needs surface sensible and latent heat fluxes.
They will come from radiation and a land-surface balance (stage S5), not from
a prescribed sine curve.

## 2. What is kept, replaced and dropped

**Kept.**
- The observation pipeline, QC and analysis (`src/analysis/`, `src/verification/`).
- Verification, including the 10 m operator.
- The maps and viewer.
- The NumPy/torch backend.
- The C-grid horizontal staggering and RK3 time stepping (Wicker and Skamarock 2002).
- The terrain pipeline (ETOPO1).
- Davies boundary relaxation and the test discipline: predictions first, class-splitting, a problem register.

**Replaced.**
- The hydrostatic core (`primitive_sigma.py`). The proposed replacement uses
  hydrostatic pressure as the vertical coordinate (Laprise 1992), with w and
  geopotential prognostic. It is split-explicit: acoustic modes are
  sub-stepped inside each RK3 stage and vertically implicit
  (Klemp et al. 2007; Skamarock and Klemp 2008). This is the WRF-ARW design, and it keeps
  the sigma-pressure coordinate NWP1 already uses. It is chosen over a height
  coordinate because the existing vertical discretisation, analysis conversion
  and diagnostics carry over.
- Richardson mixing, by a PBL scheme plus a 3D deformation-based closure.

**Dropped.**
- The dry convective adjustment. A non-hydrostatic model overturns unstable
  layers itself.
- The convective momentum mixing option.
- The prescribed heating flux, once S5 exists.

**The 12 km hydrostatic model stays** as the outer model. It is the only
source of lateral boundaries for a nested CAM in an observation-only system,
and it costs 7–9 min per 24 h.

## 3. Cost against the 1.5 h budget (estimates)

**Measured base.**
- 24 h at 12 km, 97 × 110 × 20 = 213,400 cells, dt 17.06 s (5,065 steps).
- 7.3–9.3 min on the server at torch × 8 threads, float64; 8.3 min is used below.
- That is about 4.6 × 10⁻⁷ s per cell-step.

**Scaling.** T ≈ 8.3 min × (cells / 213,400) × (17.06 s / dt) × k.
- dt ≈ 6 s per km of grid spacing, the usual split-explicit choice.
- k is the extra cost per step of the non-hydrostatic core: acoustic
  sub-steps, the w equation and the vertical implicit solve. k = 2.5 is
  assumed, and is unmeasured.
- Physics (microphysics, PBL, radiation) is not included and will add to every row.

| configuration | columns × levels | cells / now | dt | T (k = 2.5, today's settings) | T with 24 threads + float32 (assumed 2.5–4× faster) |
|---|---|---|---|---|---|
| 3 km, full domain, 50 levels | 388 × 440 × 50 | 40.0 | 18 s | ~13 h | ~3.3–5.2 h |
| 3 km, full domain, 40 levels | 388 × 440 × 40 | 32.0 | 18 s | ~10.5 h | ~2.6–4.2 h |
| 4 km, full domain, 40 levels | 291 × 330 × 40 | 18.0 | 24 s | ~4.4 h | ~1.1–1.8 h |
| 3 km nest, 600 × 600 km, 40 levels (+ 12 km parent) | 200 × 200 × 40 | 7.5 | 18 s | ~2.6 h | ~45–70 min |
| 3 km nest, 450 × 450 km, 40 levels (+ 12 km parent) | 150 × 150 × 40 | 4.2 | 18 s | ~1.5 h | ~30–45 min |

**Reading of the table.**
- A 3 km model over the whole present domain does not fit 1.5 h on the
  server's CPUs under the ≤ 26-thread rule, even with optimistic speed-ups.
- Two options are plausible:
  1. a 3 km nest over part of the domain, for example the I-95 corridor and
     the Hudson valley, inside the 12 km model;
  2. 4 km over the whole domain, if float32 and more threads deliver what
     they are assumed to.
- Both depend on k and on the speed-ups, which are unmeasured. Stage S3
  measures them before the choice.
- **Decision for the user later:** nest or full domain, and whether the 1.5 h
  budget or the thread cap may change for the CAM.

## 4. Roadmap

Each stage ends at a gate test with predictions written before it, as for
the dry model. No stage starts until the previous gate holds.

| stage | content | gate |
|---|---|---|
| S0 (now) | finish the dry 12 km fixes: near-surface wind excess (P-69), the heating instability (P-67), the night warm bias | P-69 and P-67 closed or accepted with a measurement |
| S1 | 2D (x–z) non-hydrostatic core: Laprise coordinate, split-explicit acoustic steps, vertically implicit w–geopotential solve | standard benchmarks: the Straka et al. (1993) density current within the published spread at 100–200 m; a linear hydrostatic mountain wave against its analytic solution; a rising warm bubble that conserves mass to round-off |
| S2 | 3D non-hydrostatic core at 12 km on the present domain, same boundaries and analysis | over the 20-cycle campaign, skill within ± 0.1 K and ± 0.1 m/s of the hydrostatic model (non-hydrostatic effects are small at 12 km); no instability |
| S3 | benchmark: the S2 core at 4 km and 3 km on sized arrays, float32 vs float64, 8/16/24 threads | measured k and speed-ups; pick the grid (section 3) with the user |
| S4 | moisture: water vapour with positive-definite advection, saturation adjustment, warm rain (Kessler 1969) | Bryan and Fritsch (2002) moist benchmark; dewpoint verified against ASOS; vapour mass conserved |
| S5 | surface layer with heat and moisture fluxes; a simple land-surface balance; shortwave and longwave radiation; PBL scheme | a diurnal 2 m temperature cycle with no prescribed flux; the evening warm bias smaller than with the prescribed flux |
| S6 | ice microphysics, single-moment (Lin et al. 1983) | simulated reflectivity verified against MRMS; precipitation verified with the fractions skill score (Roberts and Lean 2008) |
| S7 | the CAM in the daily cycle, nested in or replacing the 12 km model | 24 h inside 1.5 h on 20 consecutive cycles |

**Known risks, stated now.**
- **Observation-only initial conditions.** Without radar data assimilation, a
  CAM spends its first 6–12 h spinning up convection. MRMS is available, so
  latent heating nudging from radar is a candidate later. It is listed as a
  risk rather than a stage.
- **Time step.** The S1 acoustic step and the S2 large step are the two
  stability limits to measure. The present model's external-wave limit
  (dt 17 s) is set by a different wave.
- **Cost.** Section 3 may come out worse once physics is added. The nest is
  the fall-back.

## 5. Classic references for the storm-scale side

Klemp and Wilhelmson (1978) is the first three-dimensional cloud model with
the time-split treatment of sound waves that the S1 core inherits.

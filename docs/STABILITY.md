# Stability on Real Analysis Data — Findings

Status: **the dry core is stable on idealised states and unstable on real
HRRR analyses.** It survives 2–3 forecast hours before diverging. This is a
structural limitation, not a tuning problem, and it is documented here rather
than hidden behind damping that appears to help.

## What happens

Initialised from HRRR at 12 km over the Northeast:

```
initial state  max|div| 9.3e-05 1/s,  max|omega| 1.08 Pa/s      <- healthy
+  1.0 h       max|u|  43.0 m/s,      max|omega| 2.37 Pa/s
+  2.0 h       max|u|  51.8 m/s,      max|omega| 9.48 Pa/s
+  3.0 h       NaN
```

Vertical velocity roughly quadruples per hour until the run dies. The initial
state is fine — the growth happens during integration.

## What was tried, and measured

**Initial-state balancing (kept).** Analysis winds are balanced for HRRR's
discretisation, not ours. A Helmholtz split removes the divergent component,
solving the Poisson equation in Fourier space with the eigenvalues of our
*discrete* Laplacian, so cancellation is exact rather than approximate.

| | before | after |
|---|---|---|
| max\|div\| | 7.5e-04 1/s | 9.3e-05 1/s |
| implied omega | 60 Pa/s | 1.08 Pa/s |
| correlation with rotational flow | — | 0.997 |

This works and is kept. It is necessary but not sufficient.

**Divergence damping (available, default OFF).** First written as a tendency
with a coefficient in m²/s. That was itself unstable: explicit diffusion needs
`nu·dt/dx² ≤ 0.25`, and a coefficient chosen without knowing dt violates it.
Rewritten as a dimensionless post-step filter, stable by construction for any
dt.

It extended survival from 1 h to 3 h. But it damps the weather:

| div_damp | baroclinic growth/day | verdict |
|---|---|---|
| 0.00 | 1.21x | healthy |
| 0.01 | 0.33x | suppressed |
| 0.10 | 0.49x | suppressed |

Any level that helps stability also suppresses baroclinic development. Off by
default.

**Sponge layer (kept).** Amplified divergence damping near the lid, where a
rigid top reflects vertically propagating waves. The first version relaxed
wind toward the horizontal mean and was caught by the thermal-wind test —
it flattened the jet, which is legitimate structure, not wave noise.

**Stronger hyperdiffusion (does not help).** Damping times from 3 h down to
0.5 h at the grid scale, crossed with divergence damping:

| hyper damping | div_damp | hours survived |
|---|---|---|
| 3.0 h | 0.10 | 2 / 6 |
| 1.0 h | 0.10 | 2 / 6 |
| 0.5 h | 0.10 | 2 / 6 |
| 0.5 h | 0.20 | 3 / 6, then omega ~1e31 |

No damping configuration reaches 6 hours. The strongest setting is worse than
a moderate one, which is characteristic of a structural problem rather than
insufficient dissipation.

## Diagnosis

The numerics themselves are verified: the divergence operator returns
**1.08e-18** on a discretely consistent rotational flow, geostrophic imbalance
converges at exactly second order (ratios 3.97, 3.99), and all 29 idealised
tests pass. The instability appears only with realistic — sheared, stratified,
noisy — initial states.

The most likely cause is the **vertical coordinate**. In pure pressure
coordinates with a rigid flat lower boundary, `omega` must vanish at both the
lid and the ground. Diagnosing it from divergence and then forcing both
boundary conditions over-constrains the column: the linear correction that
enforces `omega = 0` at the ground redistributes error through the whole
column every step. That correction feeds vertical advection, which changes
divergence, which changes omega — a tight feedback with no physical damping.

This is precisely why operational models use terrain-following (sigma or
hybrid) coordinates. It is not only about representing mountains; the lower
boundary condition is ill-posed in pressure coordinates over any real surface.

## What this means

The core is a correct **dry hydrostatic solver for smooth balanced states**.
It is not yet usable for forecasts from real analyses.

The fix is sigma coordinates, which was already on the roadmap for terrain and
now turns out to be required for stability as well. That is a substantial
change — every vertical operator, the hydrostatic integration, and the omega
diagnosis all move to the new coordinate — but it addresses the root cause
rather than suppressing the symptom.

Interim options, in decreasing honesty:

1. **Sigma coordinates.** The real fix.
2. **Shorter forecasts.** 2 hours is genuinely integrable today.
3. Heavier damping — rejected: it suppresses the weather and still fails.

## Reproducing

```bash
python src/ingest_hrrr.py --start <date>T00 --hours 13 --stride 4
python src/forecast.py --run-dir <run> --hours 12
```

The forecast prints `max|div|` and `max|omega|` before stepping and after each
hour, so the growth is visible rather than inferred from a crash.

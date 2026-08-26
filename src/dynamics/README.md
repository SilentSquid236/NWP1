# Dynamics

A real dynamical core: the equations are integrated, not learned.

## Stage 1 — shallow water (this directory)

The 2D shallow-water equations on an Arakawa C-grid, beta-plane, doubly
periodic, stepped with the Wicker–Skamarock RK3 scheme.

Default form is **vector-invariant** (Sadourny 1975), with potential vorticity
q = (zeta+f)/h and Bernoulli function B = K + g*h:

```
du/dt = + q_bar * V_bar - dB/dx
dv/dt = - q_bar * U_bar - dB/dy
dh/dt = -d(hu)/dx - d(hv)/dy
```

Coriolis and momentum advection merge into a single vorticity flux; the
pressure gradient and kinetic-energy gradient merge into the Bernoulli term.
The textbook advective form is kept as `form="advective"` for comparison.

Why start here: these equations contain advection, Coriolis, the pressure
gradient and gravity-wave propagation — every hard part of atmospheric
dynamics except vertical structure and moisture — and they have **analytic
solutions to test against**. A bug found here takes an afternoon; the same bug
in a 3D moist model is nearly invisible, because you cannot tell a numerical
instability from real convection.

## Limited-area boundaries

`boundaries.py` implements Davies (1976) relaxation: in a band around the
perimeter the state is nudged toward externally supplied values after every
step, with alpha ~= 1 at the outermost row tapering smoothly to 0 inside.
Without this, a wave reaching the edge reflects and contaminates the interior.

`BoundaryDriver` interpolates hourly driving data (HRRR) in time -- the model
steps every few seconds, and a boundary that jumps once an hour radiates a
spurious wave inward on every update.

```python
gr    = CGrid(nx, ny, dx, dy, edge_mode="replicate")
relax = DaviesRelaxation(gr, width=12)
driver = BoundaryDriver(times, external_states)
run_limited_area(model, driver, relax, duration=24*3600)
```

This makes the model depend on HRRR at its edges. Every operational regional
model works this way, and it is categorically different from *learning* HRRR:
only the perimeter is imposed, and the interior evolves under its own physics.

### Validation

```bash
cd src/dynamics && python test_shallow_water.py   # 8 physics tests
cd src/dynamics && python test_boundaries.py      # 6 boundary tests
```

| test | what it proves | result |
|---|---|---|
| rest stays at rest | no spurious forcing | exact (0.00e+00) |
| mass conservation | flux-form continuity correct | exact (0.00e+00) |
| gravity wave speed | pressure gradient + divergence | 0.7% vs sqrt(gH) |
| geostrophic balance | Coriolis sign and staggering | 0.00% drift over 24 h |
| CFL limit | timestep bound is real | stable at 0.5x, blows up at 4x |
| energy drift | scheme not creating energy | −4.3% over 24 h at dt_max |
| drift vs dt | loss is time truncation, not scheme | 7.2x reduction when dt halved |
| enstrophy | vector-invariant vs advective | 14x better conservation |

Boundary tests:

| test | what it proves | result |
|---|---|---|
| weight taper | alpha smooth, zero in interior | exact |
| steady state | relaxation sign correct | 0.00e+00 deviation |
| zone tracks external | forcing actually applied | edge u = 5.000 (target 5.0) |
| **wave exits, no reflection** | **the point of the scheme** | **0.5% of rigid-edge case** |
| time interpolation | no hourly boundary jumps | exact, clamped outside range |
| interior mass stable | open boundary not leaking | 4.9e-05 over 12 h |

### What the vector-invariant rewrite actually changed

It was expected to fix the energy drift. **It did not**, and measuring
carefully is what revealed why: the drift is dominated by RK3 *time*
truncation, not by the spatial scheme. Both forms lose the same energy at the
same dt, and halving dt cuts the loss ~7x:

| dt | advective | vector-invariant |
|---|---|---|
| dt_max | −4.325% | −4.320% |
| /2 | −0.611% | −0.604% |
| /4 | −0.084% | −0.077% |
| /8 | −0.017% | −0.010% |

The real benefit shows up in **potential enstrophy**, the second invariant,
which governs the turbulent cascade. In a vorticity-rich shear flow the
vector-invariant form conserves it ~14x better (0.0015% vs 0.0218% over 48 h).
That matters for long integrations, and the form generalises cleanly to 3D --
which is the stronger reason to keep it.

Practical consequence: **if energy conservation matters for a run, reduce dt.**
Do not expect the spatial scheme to rescue a marginal timestep.

### Known limitations

- **Boundary driving data is not yet wired to HRRR.** `BoundaryDriver` takes
  arrays; converting ingested HRRR fields into shallow-water boundary states
  is the connecting piece still missing.
- **Second-order centred advection.** Fine here; a 5th-order upwind scheme
  (as WRF uses) is less diffusive for sharp features.
- **Enstrophy is conserved well but not exactly.** The Arakawa-Lamb scheme
  conserves both energy and enstrophy by construction; Sadourny conserves one
  at a time depending on the variant chosen.

### The CFL constraint, concretely

Explicit schemes cannot step faster than the fastest signal crossing a cell.
For H = 10 km the external gravity wave runs at ~313 m/s:

| dx | max dt |
|---|---|
| 20 km | 36 s |
| 10 km | 18 s |
| 3 km | 5.4 s |

At 3 km, a 24-hour forecast is ~16,000 steps. This is the central cost of
explicit atmospheric models, and the reason operational cores use
semi-implicit or split-explicit time stepping.

## Next stages

2. Dry hydrostatic primitive equations — 3D on pressure levels
3. Moisture — vapour transport, condensation, latent heating
4. Parameterizations — boundary layer, radiation, microphysics

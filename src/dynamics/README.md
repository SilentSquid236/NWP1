# Dynamics

A real dynamical core: the equations are integrated, not learned.

## Stage 1 — shallow water (this directory)

The 2D shallow-water equations on an Arakawa C-grid, beta-plane, doubly
periodic, stepped with the Wicker–Skamarock RK3 scheme.

```
du/dt = -u du/dx - v du/dy + f v - g dh/dx
dv/dt = -u dv/dx - v dv/dy - f u - g dh/dy
dh/dt = -d(hu)/dx - d(hv)/dy
```

Why start here: these equations contain advection, Coriolis, the pressure
gradient and gravity-wave propagation — every hard part of atmospheric
dynamics except vertical structure and moisture — and they have **analytic
solutions to test against**. A bug found here takes an afternoon; the same bug
in a 3D moist model is nearly invisible, because you cannot tell a numerical
instability from real convection.

### Validation

```bash
cd src/dynamics && python test_shallow_water.py
```

| test | what it proves | result |
|---|---|---|
| rest stays at rest | no spurious forcing | exact (0.00e+00) |
| mass conservation | flux-form continuity correct | exact (0.00e+00) |
| gravity wave speed | pressure gradient + divergence | 0.7% vs sqrt(gH) |
| geostrophic balance | Coriolis sign and staggering | 0.00% drift over 24 h |
| CFL limit | timestep bound is real | stable at 0.5x, blows up at 4x |
| energy drift | scheme not creating energy | −4.3% over 24 h |

### Known limitations

- **Momentum is in advective form**, so energy is conserved only
  approximately (the −4.3% above). A flux-form or vector-invariant momentum
  equation conserves it far better and is the right upgrade before 3D.
- **Doubly periodic only.** A limited-area domain needs a Davies relaxation
  zone driven by external boundary data.
- **Second-order centred advection.** Fine here; a 5th-order upwind scheme
  (as WRF uses) is less diffusive for sharp features.

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

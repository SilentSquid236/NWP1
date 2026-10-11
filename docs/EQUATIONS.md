# Equations of NWP1

This document lists every equation the model uses, grouped by component. For
each one it gives the equation, where it lives in the code, the parameter
values the code sets, why it was chosen, and where it came from. The reasons
are taken from the project's own record: code docstrings and comments,
`docs/CAM_DESIGN.md`, `docs/PROBLEMS.md` and `docs/RESEARCH_LOG.md`. Test and
problem identifiers (test AN, P-63, S5e) point into those files. A published
source is named only when the code or the docs name it, or when it is the
standard origin of a textbook equation. Where the record does not say,
this document says "origin not recorded" or "reason not recorded". Full
citations are in the References section at the end and in
`docs/REFERENCES.md`.

Written 2026-10-07 from the code as it stood that day. The code is the
authority. If this file and the code disagree, the code is right and this file
is out of date.

## Contents

1. Coordinates and grids
2. The hydrostatic production core (12 km)
3. The non-hydrostatic core (3 km CAM; also 12 km NH runs)
4. Physics: surface drag, vertical mixing, convective adjustment, land surface, sun
5. Moisture (warm rain, CAM stage S7)
6. Lateral boundaries and nesting
7. Initialisation and the observation analysis
8. Diagnostics and verification
9. Physical constants
10. Gaps in the record
11. Equation sources
12. References

## Conventions used throughout

- **Vertical index.** Index 0 is the model lid and index N-1 is the lowest
  full level, everywhere in NWP1 (`sigma.py`). Half levels (layer interfaces,
  w levels) run from h = 0 (lid) to h = N (ground). Full level k sits between
  half levels k and k+1.
- **Horizontal index.** Arrays are (ny, nx), row j = 0 at the southern edge,
  column i = 0 at the western edge.
- **Units.** SI throughout: Pa, K, m, s. Mixing ratios in kg/kg.
- **Two configurations are in daily use.**
  - 12 km production: hydrostatic core, 97 x 110 columns, 20 sigma levels,
    `tools/daily.sh` calling `src/forecast.py` with its defaults.
  - 3 km CAM: non-hydrostatic core with C kernels, 389 x 439 columns, 40
    levels, physics once per step, nested in the 12 km forecast
    (`tools/cam_cycle.sh`: `--core nh --nh-backend c --nh-physics step
    --levels 40 --boundary-forecast`).
- **Defaults** quoted below are those of the code. Where `forecast.py`
  overrides a module default, the forecast value is given and the module
  value is noted.
- Equation numbers are local to this document. "Eq. 2.7" means equation
  2.7; "section 2.7" means the subsection headed 2.7. They are numbered
  separately.

---

## 1. Coordinates and grids

### 1.1 Terrain-following sigma coordinate

$$
\sigma = \frac{p - p_t}{p_s - p_t}, \qquad \pi \equiv p_s - p_t, \qquad p(\sigma) = p_t + \sigma\,\pi
\tag{1.1}
$$

- $p$ pressure, $p_s$ surface pressure, $p_t$ lid pressure, $\pi$ the
  column's pressure depth (prognostic). $\sigma = 0$ at the lid and 1 at the
  ground.
- **Code.** `src/dynamics/sigma.py`: `SigmaLevels.pressure`,
  `SigmaLevels.pressure_half`.
- **Parameters.** $p_t$ = 20 000 Pa (200 hPa).
- **Why.** The pressure-coordinate core could not pose its lower boundary:
  omega had to vanish at both ends of a column that could not move, and the
  correction that enforced it fed a divergence/vertical-velocity loop. Runs
  from real analyses diverged within 2-3 hours at every damping setting
  (P-14, `docs/STABILITY.md`). In sigma the ground is $\sigma = 1$ whatever
  the terrain, $\dot\sigma = 0$ at both ends is a real boundary condition
  because $\pi$ is prognostic, and terrain comes with the coordinate. The
  200 hPa lid was measured against 100 and 50 hPa over 2500 m terrain:
  11/12, 10-11/12 and 9-10/12 forecast hours survived (docstring of
  `SigmaLevels.__init__`, re-measured 2026-09-03). "A dry model with no
  stratospheric physics gains little from a 50 hPa lid."
- **Source.** Phillips (1957).
- **Limits.** The pressure-gradient force splits into two large terms that
  nearly cancel over steep terrain (section 2.4).

### 1.2 Level placement

$$
\sigma^h_h = \left(\frac{h}{N}\right)^{s}, \quad h = 0,\dots,N; \qquad
\sigma_k = \tfrac12\left(\sigma^h_k + \sigma^h_{k+1}\right); \qquad
\Delta\sigma_k = \sigma^h_{k+1} - \sigma^h_k
\tag{1.2}
$$

- **Code.** `sigma.py`, `SigmaLevels.__init__`.
- **Parameters.** stretch $s = 1.4$. $N = 20$ at 12 km (`config.N_LEVELS`,
  the analysis's pressure-level count), $N = 40$ at 3 km (`--levels 40`).
  Lowest full level: $\sigma_{N-1}$ = 0.9654 at 20 levels (about 237-290 m
  above ground), 0.9826 at 40 levels (about 140 m).
- **Why.** The stretch concentrates resolution near the surface "where the
  boundary layer lives" (comment). Twenty levels match the analysis so that
  a different count would not silently change the vertical resolution of
  every earlier result (`forecast.py` comment). Forty levels for the CAM
  because convection needs them (S5a: T +0.116 K worse, u -0.100 and
  v -0.091 m/s better than 20 levels; recorded and not tuned).
- **Source.** Origin of the power-law stretch not recorded.
- **Limits.** P-44 (ACCEPTED): the lowest level at 237 m makes the log-law
  drag a coarse approximation.

### 1.3 Exner function and potential temperature

$$
\Pi = \left(\frac{p}{P_0}\right)^{\kappa}, \qquad T = \theta\,\Pi, \qquad \theta = T\left(\frac{P_0}{p}\right)^{\kappa}, \qquad \kappa = R_d/c_p
\tag{1.3}
$$

- **Code.** `sigma.py` (`SigmaLevels.exner`, `temperature`),
  `interpolate.theta_from_T`, and inline in most modules.
- **Parameters.** $P_0$ = 100 000 Pa. $R_d$ and $c_p$ differ between modules
  (section 9).
- **Why.** Theta is the model's thermodynamic variable because it is
  conserved by dry adiabatic motion, so interpolating it "does not invent
  heating the way interpolating T through a deep layer does"
  (`forecast.hrrr_to_sigma_state`).
- **Source.** Textbook (Poisson's equation); no source cited, by the rule in
  `docs/REFERENCES.md`.

### 1.4 Arakawa C grid and its operators

Placement: $h$ (scalars, $\pi$, $\theta$) at cell centres; $u_{j,i}$ on the
western face of cell $i$, at $x = (i - \tfrac12)\Delta x$; $v_{j,i}$ on the
southern face of cell $j$.

$$
\delta^{+}_x a = \frac{a_{i+1}-a_i}{\Delta x}, \quad
\delta^{-}_x a = \frac{a_i-a_{i-1}}{\Delta x}, \quad
\overline{a}^{u} = \tfrac12(a_i + a_{i-1}), \quad
\overline{a}^{v} = \tfrac12(a_j + a_{j-1})
\tag{1.4}
$$

$$
\overline{v}^{u}_{j,i} = \tfrac14\left(v_{j,i}+v_{j+1,i}+v_{j,i-1}+v_{j+1,i-1}\right), \qquad
\overline{u}^{v}_{j,i} = \tfrac14\left(u_{j,i}+u_{j,i+1}+u_{j-1,i}+u_{j-1,i+1}\right)
\tag{1.5}
$$

- **Code.** `src/dynamics/grid.py`, class `CGrid`: `dx_forward`,
  `dx_backward`, `dy_forward`, `dy_backward`, `h_to_u`, `h_to_v`, `v_to_u`,
  `u_to_v`, `shift`.
- **Parameters.** `edge_mode="replicate"` for the limited area (an index
  beyond the edge repeats the edge value); `"periodic"` for idealised tests.
- **Why.** "The pressure gradient and divergence become adjacent differences
  rather than 2*dx differences, which kills the checkerboard null mode that
  plagues the A-grid. It is what WRF, MPAS, and most operational models use"
  (`grid.py` docstring). P-27: `shift` maps its axis onto the last two
  dimensions so the same operator serves 2-D and 3-D fields.
- **Source.** Arakawa and Lamb (1977).

### 1.5 Domain and grid spacing

$$
\Delta y = \frac{(\varphi_{max}-\varphi_{min})\cdot 111\,132}{n_y}, \qquad
\Delta x = \frac{(\lambda_{max}-\lambda_{min})\cdot 111\,320\,\cos\varphi_0}{n_x}, \qquad
\varphi_0 = \tfrac12(\varphi_{min}+\varphi_{max})
\tag{1.6}
$$

$$
(n_y, n_x) = \Big(\mathrm{round}(L_y/d),\ \mathrm{round}(L_x/d)\Big)
\tag{1.7}
$$

- $\varphi$ latitude, $\lambda$ longitude (degrees), $d$ the target spacing.
- **Code.** `src/forecast.py`, `build_grid`; `src/analysis/geo.py`,
  `grid_shape`, `spacing_m`, `cell_centres`.
- **Parameters.** Domain 37.0-47.5 N, 82.0-66.0 W (`config.DOMAIN`).
  12 km: 97 x 110, $\Delta y$ = 12 029.8 m, $\Delta x$ = 11 985.6 m.
  3 km: 389 x 439, $\Delta y$ = 2999.7 m, $\Delta x$ = 3003.2 m.
  Cells are uniform in latitude and longitude; the model treats them as a
  Cartesian plane with $\Delta x$ evaluated at the central latitude.
- **Why.** Spacing is derived from the domain and array shape, "so a change
  to the subset bounds cannot silently desynchronise the dynamics from the
  data" (`build_grid`). The analysis uses the same formula so it places
  values on exactly the forecast's cell centres; a half-cell shift "is also
  a spurious pressure gradient everywhere" (`geo.py`).
- **Source.** Origin of the 111 132 and 111 320 m per degree constants not
  recorded.

### 1.6 Coriolis parameter on a beta plane

$$
f(y) = f_0 + \beta\,(y - y_0), \qquad f_0 = 2\Omega\sin\varphi_0, \qquad \beta = \frac{2\Omega\cos\varphi_0}{a}, \qquad y_0 = \tfrac12 L_y
\tag{1.8}
$$

- **Code.** `forecast.build_grid` (values), `grid.CGrid.__init__`
  (`f_h`, `f_u`, `f_v`).
- **Parameters.** $\Omega$ = 7.2921e-5 s^-1, $a$ = 6 371 000 m. For this
  domain $f_0$ = 9.8059e-5 s^-1, $\beta$ = 1.6945e-11 m^-1 s^-1.
  $f_u$ uses the centre rows, $f_v$ the v rows.
- **Why.** Reason for a beta plane rather than $f$ at each latitude not
  recorded.
- **Source.** Textbook.

### 1.7 Terrain

Terrain is ETOPO1 block-averaged onto the grid, ocean clipped to 0 m, then
smoothed by repeated 1-2-1 passes until its maximum slope is at most a limit:

$$
S = \max_{j,i}\sqrt{(\delta^-_x h)^2 + (\delta^-_y h)^2}
\tag{1.9}
$$

$$
\tilde h_i = \tfrac14\left(h_{i-1} + 2h_i + h_{i+1}\right)\ \text{in }x\text{, then the same in }y
\tag{1.10}
$$

$$
w_{forced} = U\,S
\tag{1.11}
$$

- **Code.** `geo.load_terrain`, `geo.block_average`, `geo.limit_slope`;
  `sigma.terrain_slope`, `sigma.smooth_terrain`, `sigma.forced_ascent`.
- **Parameters.** Slope limit 0.0086 (`ingest_obs.py --max-slope`), up to
  200 passes. ETOPO stride 2 arc-minutes at 12 km, 1 arc-minute when a cell
  is smaller than 4 arc-minutes (S5).
- **Why.** 0.0086 is the steepest idealised terrain the core survived 12/12 h
  over (2500 m ridge). Raw ETOPO at 12 km reaches 0.0536 and every real-terrain
  run died within 4 h; the same state smoothed to 0.0083 lived to 13.7 h
  (P-56, `geo.limit_slope`). Slope, not height, governs the forcing: "a 6000 m
  mountain with a gentle slope integrates LONGER than a 3000 m mountain with
  a steep one" (`sigma.forced_ascent`). At 3 km a steeper limit (0.03) gave a
  worse analysis, 1.571 K against 1.541 K (S5b).
- **Source.** Terrain data: NOAA National Geophysical Data Center (2009).
  Origin of the 1-2-1 filter not recorded.

### 1.8 Map projection (display only)

Maps are drawn on a spherical Lambert conformal conic. The model itself does
not use it.

$$
n = \frac{\ln(\cos\varphi_1/\cos\varphi_2)}{\ln\left[\tan(\tfrac\pi4+\tfrac{\varphi_2}2)/\tan(\tfrac\pi4+\tfrac{\varphi_1}2)\right]}, \quad
F = \frac{\cos\varphi_1\,\tan^n(\tfrac\pi4+\tfrac{\varphi_1}2)}{n}, \quad
\rho = \frac{aF}{\tan^n(\tfrac\pi4+\tfrac\varphi2)}
\tag{1.12}
$$

$$
x = \rho\sin\left[n(\lambda-\lambda_0)\right], \qquad y = \rho_0 - \rho\cos\left[n(\lambda-\lambda_0)\right], \qquad k = \frac{n\rho}{a\cos\varphi}
\tag{1.13}
$$

- **Code.** `src/maps/geography.py`, class `LambertConformal`.
- **Parameters.** Standard parallels 39 and 45 N, centre 42 N 74 W,
  $a$ = 6 371 000 m.
- **Why.** "The projection US regional products use, with scale error under
  1 % across 37-47.5 N" (docstring).
- **Source.** Snyder (1987).

---

## 2. The hydrostatic production core (12 km)

`src/dynamics/primitive_sigma.py`, class `PrimitiveSigma`, built by
`forecast.py`. Prognostic: $u$, $v$, $\theta$, $\pi$. Diagnosed each call:
$\Phi$ (hydrostatic), $\dot\sigma$ (continuity).

### 2.1 The prognostic system

$$
\frac{\partial \pi}{\partial t} = -\int_0^1 \nabla\cdot(\pi\mathbf V)\,d\sigma
\tag{2.1}
$$

$$
\frac{\partial u}{\partial t} = -A(u) - \mathcal V(u) + f_u\,\overline{v}^{u} + F_x + \mathcal D_u, \qquad
\frac{\partial v}{\partial t} = -A(v) - \mathcal V(v) - f_v\,\overline{u}^{v} + F_y + \mathcal D_v
\tag{2.2}
$$

$$
\frac{\partial \theta}{\partial t} = -A(\theta) - \mathcal V(\theta) + \mathcal D_\theta
\tag{2.3}
$$

- $A$ horizontal advection (section 2.6), $\mathcal V$ vertical advection
  (section 2.5), $F$ pressure-gradient force (section 2.4), $\mathcal D$ the sum of hyperdiffusion,
  divergence damping, drag, vertical mixing, sponge and (optionally) the
  stochastic factor. For $\theta$ the advecting wind is
  $\tfrac12(u_i + u_{i+1})$, $\tfrac12(v_j + v_{j+1})$.
- **Code.** `PrimitiveSigma.tendencies`.
- **Why.** Prognostic surface pressure is "the whole point": $\dot\sigma = 0$
  at lid and ground falls out of the formulation, verified to 0.00e+00 in
  `test_sigma.py` (module docstring; P-14).
- **Source.** Phillips (1957) for the coordinate. The primitive equations are
  textbook.
- **Limits.** P-42 (REVERTED): a flux-form theta was unstable "where the
  advective form was clean", so theta stays in advective form here.

### 2.2 Continuity and the vertical velocity $\dot\sigma$

$$
D_k = \delta^{+}_x\!\left(\overline{\pi}^{u}u_k\right) + \delta^{+}_y\!\left(\overline{\pi}^{v}v_k\right), \qquad
\frac{\partial\pi}{\partial t} = -\sum_{k=0}^{N-1} D_k\,\Delta\sigma_k + F_{top}
\tag{2.4}
$$

$$
\pi\,\dot\sigma_h = -\sum_{k=0}^{h-1} D_k\,\Delta\sigma_k - \sigma^h_h\,\frac{\partial \pi}{\partial t} + F_{top}, \qquad h = 0,\dots,N
\tag{2.5}
$$

- **Code.** `sigma.continuity`.
- **Parameters.** $F_{top} = 0$ (rigid lid) in production; it is the lid
  mass flux of the optional radiative boundary (Eq. 2.19).
- **Why.** With $F_{top} = 0$, $\dot\sigma$ vanishes at $\sigma = 0$ and
  $\sigma = 1$ by construction, "no correction is applied afterwards, which is
  the structural difference from the pressure-coordinate version". Writing
  the constant as $F$ and not $F(1-s)$ "is what makes sigma_dot come out as F
  at the top and exactly 0 at the ground" (docstring).
- **Source.** Textbook sigma-coordinate continuity; no source cited.

### 2.3 Hydrostatic geopotential

$$
\Phi_{N-1} = \Phi_s + R_d\,T_{N-1}\ln\frac{p_s}{p_{N-1}}, \qquad
\Phi_k = \Phi_{k+1} + R_d\,\tfrac12(T_k+T_{k+1})\ln\frac{p_{k+1}}{p_k}
\tag{2.6}
$$

- $\Phi_s = g\,h$, $h$ terrain height.
- **Code.** `sigma.hydrostatic_geopotential`.
- **Why.** Integrated upward from the surface, $d\Phi/d\ln p = -R_d T$
  (docstring). A Simmons and Burridge (1981) half-level form exists in the
  same file (`alpha_factors`, `geopotential_sb`, `pressure_gradient_sb`)
  but is not called by the model: its first use was "300x worse than the
  scheme it replaced" and was reverted (P-43). Hydrostatic consistency was
  reached instead through the form of the pressure gradient (Eq. 2.7).
- **Source.** Textbook (hypsometric equation).

### 2.4 Pressure-gradient force

Plain form, with $\ln p$ differenced directly:

$$
F_x = -\left[\delta^-_x\Phi + \overline{R_d T}^{u}\,\delta^-_x\ln p\right], \qquad
F_y = -\left[\delta^-_y\Phi + \overline{R_d T}^{v}\,\delta^-_y\ln p\right]
\tag{2.7}
$$

Reference-state form (the default, `ref_pgf=True`):

$$
T_0(\sigma_k) = \langle T\rangle_k, \qquad \psi = \Phi + R_d\,T_0\ln p, \qquad
F_x = -\left[\delta^-_x\psi + \overline{R_d(T - T_0)}^{u}\,\delta^-_x\ln p\right]
\tag{2.8}
$$

- $\langle\cdot\rangle_k$ the horizontal mean over the domain on sigma level
  $k$, recomputed every tendency call.
- **Code.** `sigma.pressure_gradient_force`; the reference profile is built
  in `PrimitiveSigma.tendencies`.
- **Why.**
  - Both signs of the sigma form "were derived wrong twice before being
    measured". The check: for an isothermal atmosphere at rest over a ridge
    the two terms must cancel (residual 2.9e-04 against terms of 6e-02).
  - Differencing $\ln p$ directly instead of expanding it as
    $\sigma\nabla\pi/p$ at cell centres took the residual for
    $T = T(p)$ over 3000 m terrain from 2.1e-03 to 8.2e-15. That residual was
    a spurious force proportional to slope, largest in the boundary layer,
    and "why forecast survival fell monotonically with terrain height"
    (docstring; P-14 "hydrostatic consistency 2.1e-03 -> 8.2e-15").
  - The reference-state form does the large cancellation analytically. The
    record is mixed on its value: the constructor comment says "the plain
    form is stable on flat ground but not over terrain", while the research
    log entry of 2026-08-28 found it had "no effect on survival time" at the
    time. It is the default in the code.
  - P-09: initialising winds from only the first term implied an 845 m/s
    "balanced" wind over 2500 m terrain; the full force gives 41.3 m/s.
- **Source.** The sigma-coordinate transformation is textbook. The
  reference-state form is called a "standard remedy" in the log; its origin
  is not recorded.

### 2.5 Vertical advection

$$
\mathcal V(a)_k = \frac{\dot\sigma_{k+1}a^h_{k+1} - \dot\sigma_k a^h_k}{\Delta\sigma_k} - a_k\,\frac{\dot\sigma_{k+1}-\dot\sigma_k}{\Delta\sigma_k}, \qquad
a^h_h = \tfrac12(a_{h-1}+a_h),\ \ a^h_0 = a_0,\ \ a^h_N = a_{N-1}
\tag{2.9}
$$

- **Code.** `sigma.vertical_advection`.
- **Why.** "Flux-form so that a constant field is transported exactly: the two
  half level fluxes cancel identically. An advective form would leave a
  residual proportional to the divergence of sigma_dot, which is precisely the
  term that misbehaves near the boundaries" (docstring).
- **Source.** Origin not recorded.

### 2.6 Horizontal advection

Second-order centred (`"centred2"`):

$$
A(a) = u\,\frac{a_{i+1}-a_{i-1}}{2\Delta x} + v\,\frac{a_{j+1}-a_{j-1}}{2\Delta y}
\tag{2.10}
$$

Third-order upwind-biased (`"upwind3"`, the default since 2026-10-03), per
direction:

$$
u\,\partial_x a \approx u\,\frac{-a_{i+2}+8a_{i+1}-8a_{i-1}+a_{i-2}}{12\Delta x} + \lvert u\rvert\,\frac{a_{i+2}-4a_{i+1}+6a_i-4a_{i-1}+a_{i-2}}{12\Delta x}
\tag{2.11}
$$

- **Code.** `PrimitiveSigma._horiz_adv`; `forecast.py --advection`
  (default `upwind3`; the class default is still `centred2`).
- **Why.** centred2 is non-dissipative, and its dispersion error can raise a
  local maximum: at the 28 Sep hot spot it was +7 to +66 m/s per hour along
  the wind at a strict local speed maximum, where the continuous term is
  zero (P-67 budget). Test AL: upwind3 completed both P-67 failures with the
  fastest wind 15.6 and 14.9 m/s (against 41.5 and 51.1 before divergence),
  with skill unchanged (T +0.002 K, u +0.005 m/s). centred2 failed on four
  real cases in all. Adopted as the default by the user after the holdout,
  test AN (2026-10-03). It is also the family the non-hydrostatic core uses.
- **Source.** Wicker and Skamarock (2002).

### 2.7 Hyperdiffusion

$$
\mathcal D^{(4)}(a) = -K_4\,\mathcal L\!\left(\mathcal L(a)\right), \qquad
\mathcal L(a) = \frac{a_{i+1}-2a_i+a_{i-1}}{\Delta x^2} + \frac{a_{j+1}-2a_j+a_{j-1}}{\Delta y^2}
\tag{2.12}
$$

$$
K_4 = \frac{1}{\tau_4\left(4/\Delta x^2 + 4/\Delta y^2\right)^2}
\tag{2.13}
$$

- Applied to $u$, $v$ and $\theta - \langle\theta\rangle_k$. **Not** applied
  to $\pi$.
- **Code.** `src/dynamics/subgrid.py`: `hyperdiffusion`,
  `discrete_biharmonic_eigenvalue`, `recommended_hyper_coeff`; call in
  `PrimitiveSigma.tendencies`.
- **Parameters.** $\tau_4$ = 3 h. $K_4$ = 3.008e10 m^4/s at 12 km
  (`--hyper-factor` 1). Stability limit $\Delta t \le 0.5\,\Delta x^4/(16K_4)$.
- **Why.** Fourth-order damping goes as $k^4$, "so halving the wavelength
  multiplies the damping by 16, leaving synoptic features essentially
  untouched". $K_4$ uses the eigenvalue of the **discrete** operator at
  2$\Delta x$; the continuous $(\pi/\Delta x)^4$ made it about six times too
  weak (18 h e-folding where 3 h was intended; P-17). No diffusion on $\pi$:
  on a bounded domain it is a mass source (p_s inflated from 1088 to
  1243 hPa in three hours; P-18).
- **Source.** "Three hours at the grid scale is a common operational choice"
  (docstring); origin not recorded.

### 2.8 Divergence damping

$$
D = \delta^+_x u + \delta^+_y v, \qquad
\mathcal D_u = \nu\,\delta^-_x D, \qquad \mathcal D_v = \nu\,\delta^-_y D, \qquad
\nu = C\,\frac{\Delta x\,\Delta y}{\Delta t_0}
\tag{2.14}
$$

- A one-dimensional wave of $n\Delta x$ decays at
  $(2\sin(\pi/n))^2\,\nu/\Delta x^2$; stability limit
  $\nu\,(4/\Delta x^2+4/\Delta y^2)\,\Delta t \le 2.5\times 0.4$.
- **Code.** `subgrid.divergence_damping`,
  `subgrid.divergence_damping_stability_dt`; $\nu$ is set in `forecast.main`.
- **Parameters.** $C$ = 0.0064 (`--div-damp`). $\Delta t_0$ is `max_dt()` of
  the loaded initial state. $\nu \approx$ 5.4-5.9e4 m^2/s at 12 km; 2$\Delta x$
  e-folds in about 10 min, 10$\Delta x$ in about 105 min (P-60).
- **Why.** P-60 was a divergent 2-3$\Delta x$ mode at jet level that killed
  the 12Z 2026-09-25 case at 4 h. Test V: with damping no point changed by
  more than 5 m/s in 0-6.25 h, the resolved maxima unchanged. Test W: every
  real case then ran 24 h. Default since 2026-09-26 (P-56, P-60). The damping
  acts on the divergent wind only, so balanced flow is untouched.
- **History.** P-19 records an earlier version as a dimensionless post-step
  filter, default off, because every setting that helped stability
  suppressed baroclinic growth (1.21x/day at 0, 0.33x/day at 0.01). The
  present tendency form replaced it for P-60.
- **Source.** Skamarock and Klemp (1992).

### 2.9 Upper sponge

$$
R_k = r\cdot\tfrac12\left[1-\cos\!\left(\pi\,\frac{K_s - k}{K_s}\right)\right],\quad k < K_s; \qquad R_k = 0,\ k \ge K_s
\tag{2.15}
$$

$$
\left.\frac{\partial u}{\partial t}\right|_{sponge} = -R_k\,(u - u_{ref}), \qquad \left.\frac{\partial v}{\partial t}\right|_{sponge} = -R_k\,(v - v_{ref})
\tag{2.16}
$$

- $u_{ref}$, $v_{ref}$ frozen at the first step (the prepared initial state).
- **Code.** `PrimitiveSigma.__init__` (profile), `set_reference`,
  `tendencies`.
- **Parameters.** $r$ = 1/900 s^-1. $K_s$ = 5 at 20 levels;
  `forecast.py` sets $K_s = \max(1, \mathrm{round}(5N/20))$, so 10 at 40
  levels, keeping the sponge base at $\sigma^h$ = 0.1436. At 5 levels the
  weights are 1.0, 0.905, 0.655, 0.345, 0.095 of $r$.
- **Why.** The rigid lid reflects vertically propagating gravity waves from
  terrain. Relaxing to a frozen reference rather than the horizontal mean,
  because relaxing to the mean flattened a jet (P-16; thermal-wind jet then
  held to 0.24 % over 24 h). Depth 5: cutting it to 3 regressed two suites
  and improved nothing (constructor comment; P-49). The rate "was chosen, not
  derived" (`sponge_rate_ladder.py`). Reason for the raised-cosine profile
  not recorded.
- **Source.** Origin not recorded.
- **Limits.** P-49 (OPEN): every sponge setting turns a growing baroclinic
  wave into a decaying one, because the sponge base (about 301 hPa) is in the
  upper troposphere. P-02 (OPEN): the sponge base partially reflects.

### 2.10 Time stepping

Three-stage Runge-Kutta:

$$
X^{*} = X^n + \frac{\Delta t}{3}\,\mathcal F(X^n), \qquad
X^{**} = X^n + \frac{\Delta t}{2}\,\mathcal F(X^{*}), \qquad
X^{n+1} = X^n + \Delta t\,\mathcal F(X^{**})
\tag{2.17}
$$

Then, on the completed step and in this order: the prescribed surface heat
flux if on (section 4.11); the land surface (section 4.8); the convective
adjustment (section 4.7); the radiative-lid reference update if on. The
lateral relaxation (section 6.1) follows in `forecast.run_forecast`.

- **Code.** `PrimitiveSigma.step`.
- **Why.** Adjustments "enforce an inequality on the state and have no
  meaningful time derivative", so they act on the completed step; "an
  intermediate stage would otherwise re-create the instability the final
  state is meant to be free of". The heat flux goes first so the adjustment
  mixes a layer the sun has just made unstable in the same step (comments).
- **Source.** Wicker and Skamarock (2002).

### 2.11 Time step

$$
\Delta t = \min\left(\frac{0.6\,\min(\Delta x,\Delta y)}{\sqrt2\,(c_L + \max\lvert\mathbf u\rvert)},\ \frac{0.6\,\min\Delta\sigma}{\max\lvert\dot\sigma\rvert},\ \Delta t_{K_4},\ \Delta t_{\nu}\right), \qquad c_L = \sqrt{R_d\,T_{max}}
\tag{2.18}
$$

- $T_{max}$ is clipped to 150-350 K. The step is re-checked every 50 steps
  (`adaptive=True` in `run`). `forecast.run_forecast` shortens it so a whole
  number of steps fits each output hour (P-64), giving about 17 s at 12 km.
  A run stops if $\max(\lvert u\rvert,\lvert v\rvert) >$ 150 m/s (`U_CEILING`,
  P-52, P-57).
- **Code.** `PrimitiveSigma.max_dt`, `PrimitiveSigma.run`,
  `forecast.run_forecast`.
- **Why.** With prognostic surface pressure the fastest signal is the
  external (Lamb) wave at about 300 m/s; carrying over the
  pressure-coordinate core's 100 m/s "makes dt about 3x too large", and
  surface pressure goes negative within about 20 steps (docstring). Reason
  for the 0.6 safety factor not recorded.
- **Source.** Textbook CFL condition.

### 2.12 Options present in the code but off in production

**Radiative upper boundary** (`radiation.py`, `radiative_top=False`):

$$
\hat\omega_{top}(\mathbf k) = -\rho_{top}\,g\,\frac{\lvert\mathbf k\rvert}{N_{top}}\,\hat\phi'_{top}(\mathbf k)
\tag{2.19}
$$

with $\phi'$ the lid geopotential minus a running low-pass ($\tau$ = 3 h),
$N_{top}$ one value from the top four layers, a cosine edge taper and a
spectral roll-off (8-16$\Delta x$). It enters as $F_{top}$ in Eqs. 2.4-2.5.
The sign is a parameter (`radiation_sign=-1`) measured by `test_radiation.py`.
Sources: Klemp and Durran (1983); Bougeault (1983). Off because it "still
dies over tall terrain, at hour 4" (P-50) and the sponge stays "until
something better is demonstrated" (constructor comment, P-49).

**Stochastic tendency perturbations** (`subgrid.StochasticPerturbation`,
`--stochastic`):

$$
\mathcal F \to \mathcal F\,(1 + r), \qquad r \leftarrow a\,r + \sqrt{1-a^2}\,r_{new}, \quad a = e^{-\Delta t/\tau}, \qquad \hat r_{new} \propto e^{-k^2L^2/2}\,\hat\xi
\tag{2.20}
$$

$\xi$ white Gaussian noise, the filtered field normalised to unit variance,
times the amplitude, clipped to $\pm$0.9. Forecast settings: amplitude 0.3,
$\tau$ = 6 h, $L$ = 300 km (capped at a quarter of the domain, P-28). Source:
Buizza et al. (1999).

---

## 3. The non-hydrostatic core

`src/dynamics/nh3d.py`: class `NH3D` (the dynamics) and `NHModel` (the adapter
that looks like `PrimitiveSigma` to `forecast.py`). Compiled counterparts are in
`nh3d_kernels.c` (kernels `td_*`, `ac_*`, `as_*`, `adv_q`, `q_stage`,
`phys_col`, `pav_col`, `moist_col`, `relax_col`, `relax_mu`). The C file states
that every kernel "is a transcription of a NumPy expression in nh3d.py, in
the same order of operations, and is tested against it (test_nh3d_c.py). The
NumPy code stays the reference." The equations below are therefore those of
the NumPy code.

Design choice (`docs/CAM_DESIGN.md` section 2): hydrostatic pressure as the
vertical coordinate, with w and geopotential prognostic, split-explicit with
vertically implicit acoustic substeps, "the WRF-ARW design", chosen over a
height coordinate "because the existing vertical discretisation, analysis
conversion and diagnostics carry over". Built in stages with benchmark gates:
S1a 2-D height coordinate (`nh2d.py`, Straka et al. 1993 density current),
S1b 2-D mass coordinate (`nh2d_mass.py`, Smith 1979 mountain wave), S2 3-D
(`nh3d.py`, test AO).

### 3.1 Variables

$$
\mu = \pi = p_s - p_t, \quad U = \overline{\mu}^{u}u, \quad V = \overline{\mu}^{v}v, \quad W = \mu\,w, \quad \Theta = \mu\,\theta\ (\mu\,\theta_\rho\ \text{with moisture}), \quad \phi = g z
\tag{3.1}
$$

- $U$, $V$, $\Theta$ on full (mass) levels; $W$ and $\phi$ on the N+1 w levels;
  $\phi_N = g h$ fixed.
- **Code.** `NH3D.__init__`, `NHModel.set_state`.
- **Source.** Laprise (1992); Skamarock and Klemp (2008).

### 3.2 Inverse density and equation of state

$$
\alpha_k = \frac{\phi_k - \phi_{k+1}}{\mu\,\Delta\sigma_k}, \qquad
p_k = P_0\left(\frac{R_d\,\Theta_k}{\mu\,P_0\,\alpha_k}\right)^{\gamma}, \qquad \gamma = \frac{c_p}{c_v},\ c_v = c_p - R_d
\tag{3.2}
$$

- **Code.** `NH3D.diagnose`.
- **Parameters.** The core's constants come from `nh2d.py`: $R_d$ = 287.04,
  $c_p$ = 1004.5, $c_v$ = 717.46 J/(kg K), $\gamma$ = 1.40008.
- **Why.** Pressure comes from the equation of state, so it is
  non-hydrostatic; with $\theta_\rho$ in $\Theta$ the same equation is exact
  for moist air (section 5.1).
- **Source.** Laprise (1992); Skamarock and Klemp (2008).

### 3.3 Continuity and $\Omega$

$$
\frac{\partial\mu}{\partial t} = -\sum_{k} (\nabla\!\cdot\!\mathbf V)_k\,\Delta\sigma_k, \qquad
(\nabla\!\cdot\!\mathbf V)_k = \delta^+_x U_k + \delta^+_y V_k
\tag{3.3}
$$

$$
\Omega_0 = 0, \qquad \Omega_{h+1} = \Omega_h - \left(\frac{\partial\mu}{\partial t} + (\nabla\!\cdot\!\mathbf V)_h\right)\Delta\sigma_h, \qquad \Omega_N = 0, \qquad \dot\sigma = \Omega/\mu
\tag{3.4}
$$

- $\Omega = \mu\dot\sigma$, positive downward (toward larger $\sigma$).
- **Code.** `NH3D.tendencies`, `NH3D.omega`.
- **Why.** Column mass in flux form "so mass is conserved to round-off"
  (`nh2d_mass.py`; S1b mass change -1.3e-16).

### 3.4 Momentum and the perturbation pressure gradient

Reference state: $\bar\mu$, $\bar p_k = p_t + \sigma_k\bar\mu$, $\bar\phi$,
$\bar\alpha$ (section 3.7). Perturbations $p' = p - \bar p$,
$\alpha' = \alpha - \bar\alpha$, $\mu' = \mu - \bar\mu$, $\phi' = \phi - \bar\phi$.
$\overline{\phi}^{\sigma}$ is the average of the two w levels around a mass
level.

$$
\mathrm{PGF}_x = \overline{\mu\alpha}^{u}\,\delta^-_x p' + \overline{\mu\alpha' + \mu'\bar\alpha}^{u}\,\delta^-_x\bar p + \overline{\bar\mu}^{u}\,\delta^-_x\overline{\phi'}^{\sigma} + \overline{\partial_\sigma p'}^{u}\,\delta^-_x\overline{\phi}^{\sigma}
\tag{3.5}
$$

$$
\frac{\partial U}{\partial t} = -\mathrm{PGF}_x + f_u\,\overline{\overline{V}^{h}}^{u} - \overline{\mu}^{u}\left[A_3(u; u) + A_3(u; \overline v) + \overline{\dot\sigma}\,\partial_\sigma u\right] + \mathcal P_U
\tag{3.6}
$$

and the same for $V$ with $-f_v\,\overline{\overline{U}^{h}}^{v}$.

- $\partial_\sigma p'$ at mass levels from w-level values: $p'_w = 0$ at the
  lid, the two-level average inside, linear extrapolation at the ground.
  $A_3$ is Eq. 2.11. $\partial_\sigma u$ is centred, one-sided at the ends.
  $\mathcal P_U$ is the physics (section 3.10).
- **Code.** `NH3D.tendencies` (variables `pp`, `alp`, `mup`, `pgx`, `pgy`,
  `FU`, `FV`), `NH3D._u3`.
- **Why.** "The horizontal pressure gradient uses the perturbation form
  against a hydrostatic reference sounding, with the reference balance
  removed analytically", so the large opposing terms over terrain cancel
  exactly: the reference balance
  $\bar\mu\bar\alpha\,\partial_x\bar p + \bar\mu\,\partial_x\bar\phi = 0$ is
  subtracted (`nh2d_mass.py` docstring). Test: an atmosphere at rest over a
  2500 m mountain, different from the reference profile, shows no growth in
  24 h (S2, `docs/nh3d_terrain_rest.png`), where the first hydrostatic core
  diverged at 3 h.
- **Source.** Laprise (1992); Skamarock and Klemp (2008).

### 3.5 Vertical momentum and geopotential

$$
\frac{\partial W_h}{\partial t} = g\left(\frac{p'_h - p'_{h-1}}{\Delta\sigma^w_h} - \mu'\right) - \mu\left(u_w\frac{w_{i+1}-w_{i-1}}{2\Delta x} + v_w\frac{w_{j+1}-w_{j-1}}{2\Delta y} + \dot\sigma\,\partial_\sigma w\right) + \mathcal P_W, \qquad h < N
\tag{3.7}
$$

$$
\frac{\partial\phi_h}{\partial t} = -\left(u_w\frac{\phi_{i+1}-\phi_{i-1}}{2\Delta x} + v_w\frac{\phi_{j+1}-\phi_{j-1}}{2\Delta y} + \dot\sigma\,\partial_\sigma\phi\right) + g\,w, \qquad \frac{\partial\phi_N}{\partial t} = 0,\ \ W_N = 0
\tag{3.8}
$$

- $\Delta\sigma^w_0 = \sigma_0$, $\Delta\sigma^w_h = \sigma_h - \sigma_{h-1}$,
  $\Delta\sigma^w_N = 1-\sigma_{N-1}$; $p'_{-1} = 0$ at the lid. $u_w$, $v_w$
  are cell-centre winds averaged to w levels.
- **Code.** `NH3D.tendencies` (`FW`, `Fphi`).
- **Source.** Laprise (1992); Skamarock and Klemp (2008).

### 3.6 Potential temperature, flux form, third order

$$
\frac{\partial\Theta}{\partial t} = -\left[\delta^+_x\mathcal F_3(U,\theta) + \delta^+_y\mathcal F_3(V,\theta)\right] - \frac{\Omega_{k+1}\theta^w_{k+1} - \Omega_k\theta^w_k}{\Delta\sigma_k}
\tag{3.9}
$$

$$
\mathcal F_3(U,q)_{i-1/2} = U\,\frac{7(q_i+q_{i-1}) - (q_{i+1}+q_{i-2})}{12} - \lvert U\rvert\,\frac{3(q_i-q_{i-1}) - (q_{i+1}-q_{i-2})}{12}
\tag{3.10}
$$

$$
\theta^w_h = \frac{7(\theta_{h-1}+\theta_h) - (\theta_{h-2}+\theta_{h+1})}{12} - \mathrm{sgn}(\Omega_h)\,\frac{3(\theta_h-\theta_{h-1}) - (\theta_{h+1}-\theta_{h-2})}{12}, \quad 2 \le h \le N-2
\tag{3.11}
$$

with the two-level average at $h$ = 1 and N-1 and the end values at
$h$ = 0 and N.

- **Code.** `NH3D._flux3`, `NH3D._vface3`, `NH3D.tendencies` (`FTh`).
- **Why.** S1b: a second-order centred theta flux undershot at the sharp cold
  front of the density current (min theta' -18.2 K against an initial
  -15 K). With the third-order flux in x and in the column interior: front
  15.25 km, min theta' -9.99 K at 100 m, against about 15.5 km and -9.8 K in
  the reference (research log 2026-10-02, "the S1b fix").
- **Source.** Wicker and Skamarock (2002) for the third-order upwind-biased
  family (named in `nh2d.py` and `primitive_sigma.py` as `upwind3`).

### 3.7 Reference sounding

$$
\bar\theta(z) = 250\,\exp\!\left(\frac{g z}{c_p\cdot 250}\right), \qquad
\bar\Pi_{i} = \bar\Pi_{i-1} - \frac{g\,\Delta z}{c_p\,\tfrac12(\bar\theta_i + \bar\theta_{i-1})}, \quad \bar\Pi_0 = 1\ (p = P_0\ \text{at}\ z = 0)
\tag{3.12}
$$

$$
\bar\mu = \bar p(h) - p_t, \quad \bar\phi_h = g\,\bar z(p_t + \sigma^h_h\bar\mu), \quad
\bar\alpha_k = \frac{\bar\phi_k - \bar\phi_{k+1}}{\bar\mu\Delta\sigma_k}, \quad
\bar\theta_k = \frac{P_0\,\bar\alpha_k}{R_d}\left(\frac{\bar p_k}{P_0}\right)^{1/\gamma}
\tag{3.13}
$$

- $\bar z(p)$ by interpolation in $\ln p$ of the fine-grid sounding.
- **Code.** `nh2d_mass.ReferenceSounding`; `NH3D.__init__`.
- **Parameters.** Isothermal 250 K; $\Delta z$ = 5 m; integrated to
  $z_{top}$ + 2000 m with $z_{top}$ = 20 000 m.
- **Why.** The reference "only has to be hydrostatic" (comment). $\bar\theta$
  is set "consistent with the discrete EOS so that p' = 0 at rest"
  (`nh2d_mass.py`). Reason for 20 000 m not recorded.
- **Source.** Origin not recorded.

### 3.8 Split-explicit time integration

RK3 as Eq. 2.17, but each stage advances the stage state $X^*$ by $n$
acoustic substeps of $\Delta\tau = \Delta t/n_s$, with $n = n_s/3,\ n_s/2,\ n_s$
in the three stages. In each stage the slow tendency $\mathcal F(X^*)$ is held
fixed and the perturbations $X'' = X - X^*$ start from $X^n - X^*$.

- **Code.** `NH3D.step`, `NH3D._acoustic`.
- **Parameters.** $n_s$ = 6 (must be a multiple of 6), off-centring
  $\beta$ = 0.1, divergence damping coefficient $\beta_d$ = 0.1.
- **Why.** The design document's choice (introduction to section 3). In S1a
  the density current at $n_s$ = 12 matched $n_s$ = 6 (front 15.35 km in
  both, min theta' -10.20 against -10.22 K at 100 m). Reason for the specific values 0.1 and 0.1 not recorded.
- **Source.** Wicker and Skamarock (2002); Klemp et al. (2007).

### 3.9 The acoustic substep

Linearised about $X^*$:

$$
c^2 = \gamma\,p^*\alpha^*, \qquad
p''_k = \frac{c^2}{\alpha^*\Theta^*}\,\Theta''_k - \frac{c^2}{\alpha^{*2}\mu^*\Delta\sigma_k}\left(\phi''_k - \phi''_{k+1}\right)
\tag{3.14}
$$

$$
p_d = p''^{\tau} + \beta_d\left(p''^{\tau} - p''^{\tau-\Delta\tau}\right)
\tag{3.15}
$$

$$
U''^{\tau+\Delta\tau} = U''^{\tau} + \Delta\tau\left[F_U - \overline{\mu^*}^{u}\overline{\alpha^*}^{u}\,\delta^-_x p_d - \overline{\mu^*}^{u}\,\delta^-_x\overline{\phi''}^{\sigma}\right]
\tag{3.16}
$$

then $\mu''$ and $\Omega''$ from Eqs. 3.3-3.4 with $\mathbf V''$, and

$$
\Theta''^{\tau+\Delta\tau} = \Theta''^{\tau} + \Delta\tau\left[F_\Theta - \nabla\!\cdot\!\left(\mathbf V''\,\overline{\theta^*}\right) - \frac{\Omega''_{k+1}\overline{\theta^*}^w_{k+1} - \Omega''_k\overline{\theta^*}^w_k}{\Delta\sigma_k}\right]
\tag{3.17}
$$

$$
\Phi = \phi''^{\tau} + \Delta\tau\left[F_\phi - \frac{U''_w\partial_x\phi^* + V''_w\partial_y\phi^* + \Omega''\partial_\sigma\phi^*}{\mu^*} + a_-\,\frac{g\,W''^{\tau}}{\mu^*}\right], \quad
W_e = W''^{\tau} + \Delta\tau\left[F_W + g\left(a_-\,\delta_\sigma p_d - \mu''\right)\right]
\tag{3.18}
$$

Vertically implicit $W''$, per column, $h = 0,\dots,N-1$:

$$
\left(1 + G_hB\,(E_h + E_{h-1})\right)W''_h - G_hB\,E_h\,W''_{h+1} - G_hB\,E_{h-1}\,W''_{h-1} = W_{e,h} + G_h\,S_h
\tag{3.19}
$$

$$
G_h = \frac{\Delta\tau\,g\,a_+}{\Delta\sigma^w_h}, \quad B = \frac{\Delta\tau\,a_+\,g}{\mu^*}, \quad E_k = \frac{c^2}{\alpha^{*2}\mu^*\Delta\sigma_k}, \quad a_\pm = \tfrac12(1\pm\beta), \quad \phi''^{\tau+\Delta\tau} = \Phi + B\,W''^{\tau+\Delta\tau}
\tag{3.20}
$$

- $S_h = p''_h(\Theta'',\Phi) - p''_{h-1}(\Theta'',\Phi)$ (with $p''_{-1} = 0$);
  $E_{-1} = 0$; $W''_N = 0$; $\phi''_N = 0$. Solved by the Thomas algorithm.
  $\overline{\theta^*}^w$ in Eq. 3.17 is the plain two-level average.
- **Code.** `NH3D._acoustic`, `thomas`.
- **Why.** The time-split treatment of sound waves (CAM_DESIGN section 5);
  "the vertically implicit w-pressure coupling" checked first on flat ground
  in S1a, where the sound speed came out within 0.15 %.
- **Source.** Klemp et al. (2007); Klemp and Wilhelmson (1978) for time-split
  sound waves; Skamarock and Klemp (1992) for the pressure-based divergence
  damping (named in `nh2d.py`). The Thomas algorithm is textbook.

### 3.10 Physics coupling and the large step

- Physics tendencies are mass-weighted: $\mathcal P_U = \overline{\mu}^u\,\Delta u_{phys}$,
  $\mathcal P_\Theta = \mu\,\Delta\theta_{mix}$ (times the $\theta_\rho$
  factor with moisture), from the drag (Eq. 4.6) and mixing (Eqs. 4.7-4.9).
- Sponge on the mass-weighted wind: $\mathcal P_U \mathrel{-}= R_k\,(U - \overline{\mu}^u u_{ref})$,
  $u_{ref}$ the state given to `set_state`.
- **w damping:**

$$
\mathcal P_{W,h} = -2R_h\,W_h, \qquad h = 0,\dots,N-1
\tag{3.21}
$$

- `physics_every="step"` (the 3 km setting) evaluates all of this once from
  the step-start state and holds it over the three stages.
- **Code.** `NHModel._physics`, `NHModel._physics_c`, `NHModel.step`.
- **Why.** "Physics once per step ... as WRF does": S5a matched stage physics
  to 0.000 in T, u and v over 20 cycles and cut the evaluations from 15 to 5
  per 5 steps. Reason for the factor 2 in the w damping not recorded.

$$
\Delta t_{NH} = \min\left(60\ \mathrm{s},\ \frac{0.7\,\min(\Delta x,\Delta y)}{\max(\lvert u\rvert,\lvert v\rvert,10)},\ \frac{n_s\cdot 0.5\,\min(\Delta x,\Delta y)}{350}\right)
\tag{3.22}
$$

- **Code.** `NHModel.max_dt` (`C_SOUND` = 350 m/s, `ACOUSTIC_COURANT` = 0.5).
- **Why.** The acoustic limit is 103 s at 12 km and never binds; at 3 km it is
  25.7 s, and without it the advective limit allowed about 50 s, "a substep
  Courant number near 1 for sound" (docstring, S5). Reasons for the 60 s cap
  and the 0.7 factor not recorded.

### 3.11 The two-dimensional benchmark cores

These are not used in forecasts. They are the gates S1a and S1b.

`nh2d.py` (height coordinate, Exner form):

$$
\frac{du}{dt} = -c_p\theta\,\partial_x\pi', \quad
\frac{dw}{dt} = -c_p\theta\,\partial_z\pi' + g\frac{\theta'}{\bar\theta}, \quad
\frac{d\theta'}{dt} = -w\,\frac{d\bar\theta}{dz}, \quad
\frac{\partial\pi'}{\partial t} = -\frac{\bar c^2}{c_p\bar\rho\bar\theta^2}\,\nabla\!\cdot(\bar\rho\bar\theta\mathbf V)
\tag{3.23}
$$

with $\bar c^2 = (c_p/c_v)R_d\bar\pi\bar\theta$ and a constant-K diffusion
($K$ = 75 m^2/s for the density current). Sources: Klemp and Wilhelmson
(1978); Skamarock and Klemp (1992). `nh2d_mass.py` is the 2-D form of
sections 3.2-3.9, plus a $\sin^2$ Rayleigh sponge used only in its tests.
Gates: Straka et al. (1993) density current; Smith (1979) linear mountain
wave (correlation 0.92 and slope 0.86 against the analytic solution at 120
levels).

---

## 4. Physics

### 4.1 Lowest-level height and layer depth

$$
z_1 = \frac{R_d T_{N-1}}{g}\ln\frac{p_s}{p_{N-1}}, \qquad \Delta z_1 = \frac{R_d T_{N-1}}{g}\ln\frac{p_s}{p^h_{N-1}}
\tag{4.1}
$$

- **Code.** `src/dynamics/surface.py`: `lowest_level_height`, `surface_drag`.
- **Why.** $z_1$ "must come from the model's own layer thickness, not a
  hardcoded constant, or the drag silently changes with the vertical grid"
  (docstring).
- **Source.** Textbook (hypsometric equation).

### 4.2 Neutral drag coefficient

$$
C_{d,n} = \left(\frac{\kappa_{vK}}{\ln(z_1/z_0)}\right)^2, \qquad z_1 \ge 2z_0
\tag{4.2}
$$

- **Code.** `surface.neutral_drag_coefficient`.
- **Parameters.** $\kappa_{vK}$ = 0.4. Roughness map since 2026-10-03 (test
  AM): $z_0$ = 1.0 m over land, 0.0002 m over water; land is raw ETOPO
  terrain above 0 m (P-65). `--z0 0.1` is the old single value.
- **Why.** Mixing redistributes momentum within a column; "only drag REMOVES
  it". The roughness map: test AM, u RMSE -0.274 m/s against AK (24 of 24
  leads); the values are the forest and open-sea entries of `surface.py`'s
  table, "not fitted". P-45 (one roughness length) was accepted until this.
- **Source.** Log law (textbook). The roughness values are called
  "conventional" in `surface.py`; their source is not named.

### 4.3 Bulk Richardson number

$$
Ri_b = \frac{g\,z_1\,(\theta_1 - \theta_s)}{\theta_1\,\max(u_1^2+v_1^2,\ 0.01)}
\tag{4.3}
$$

- **Code.** `surface.bulk_richardson`.
- $\theta_s$ is the ground potential temperature from the land surface
  (section 4.8). Without a land surface $\theta_s$ is `None` and the drag is
  neutral.

### 4.4 Stability functions for the drag

Cut-off form (used when a surface temperature is given without the land
surface; not in production):

$$
F(Ri_b) = \begin{cases}\min(1 + 10\lvert Ri_b\rvert,\ 4) & Ri_b < 0\\ (1 - Ri_b/0.2)^2 & 0 \le Ri_b < 0.2\\ 0 & Ri_b \ge 0.2\end{cases}
\tag{4.4}
$$

Louis (1979) momentum function (used with the land surface, the production
setting):

$$
F_m = \begin{cases}\dfrac{1}{1 + 2b\,Ri\,/\sqrt{1 + d\,Ri}} & Ri \ge 0\\[2ex] 1 - \dfrac{2b\,Ri}{1 + 3b\,c_m\,a^2\sqrt{-Ri\,z/z_0}} & Ri < 0\end{cases}, \qquad a^2 = C_{d,n}
\tag{4.5}
$$

- **Code.** `surface.stability_function`; `land_surface.louis_momentum`.
- **Parameters.** $b$ = 5, $d$ = 5, $c_m$ = 7.4; $z/z_0 = \max(z_1,2z_0)/z_0$.
- **Why.** "Neglecting F entirely would apply daytime drag to a calm, clear
  night." The cut-off form drops drag to zero at $Ri_b$ = 0.2, which with a
  cooling ground "decouples the surface completely. The ground would then
  cool without limit." Louis's long stable tail weakens the drag at night
  without switching it off (`land_surface.py`).
- **Source.** Louis (1979) for Eq. 4.5. Eq. 4.4 is called "Louis-type"; the
  constants 10, 4 and 0.2 have no recorded origin.

### 4.5 Surface drag tendency

$$
\left.\frac{\partial u}{\partial t}\right|_{N-1} = -\frac{C_d\,\lvert\mathbf V_1\rvert\,u_1}{\Delta z_1}, \qquad
\left.\frac{\partial v}{\partial t}\right|_{N-1} = -\frac{C_d\,\lvert\mathbf V_1\rvert\,v_1}{\Delta z_1}, \qquad C_d = C_{d,n}F, \qquad u_* = \sqrt{C_d}\,\lvert\mathbf V_1\rvert
\tag{4.6}
$$

- $\lvert\mathbf V_1\rvert = \sqrt{u_1^2+v_1^2}$ with $u_1$, $v_1$ taken at the
  same array index (their own C-grid points).
- **Code.** `surface.surface_drag`.
- **Why.** Bulk aerodynamic stress $\tau = \rho C_d\lvert\mathbf V\rvert\mathbf V$
  on the lowest level only; "the mixing scheme spreads the effect upward from
  there, which is how a boundary layer actually deepens". Checked by the
  Ekman spiral (docstring; P-29).
- **Source.** Textbook bulk formula.

### 4.6 Richardson-number vertical mixing

At each interior interface $k+\tfrac12$ (between full levels $k$ and $k+1$):

$$
\Delta z = \max\left(\left\lvert\frac{R_d\,\bar T}{g}\ln\frac{p_{k+1}}{p_k}\right\rvert,\ 1\right), \quad
N^2 = \frac{g}{\bar\theta}\,\frac{\theta_k - \theta_{k+1}}{\Delta z}, \quad
S^2 = \left(\frac{\Delta u}{\Delta z}\right)^2 + \left(\frac{\Delta v}{\Delta z}\right)^2, \quad
Ri = \frac{N^2}{S^2}
\tag{4.7}
$$

$$
K = \min\left(l^2\,S\,f(Ri),\ K_{max}\right), \qquad
f(Ri) = \begin{cases}1 & Ri \le 0\\ (1 - Ri/Ri_c)^2 & 0 < Ri < Ri_c\\ 0 & Ri \ge Ri_c\end{cases}
\tag{4.8}
$$

$$
\left.\frac{\partial a}{\partial t}\right|_k = \frac{\mathcal F_{k+1/2} - \mathcal F_{k-1/2}}{\tfrac12(\Delta z_{k+1/2} + \Delta z_{k-1/2})}, \qquad
\mathcal F_{k+1/2} = K_{k+1/2}\,\frac{a_{k+1} - a_k}{\Delta z_{k+1/2}}, \qquad \mathcal F_{-1/2} = \mathcal F_{N-1/2} = 0
\tag{4.9}
$$

(the top and bottom levels divide by the single interface's $\Delta z$).
$Ri = \infty$ where $S^2 \le$ 1e-12.

- Applied to $u$, $v$, $\theta$ and, with moisture, $q_v$, $q_c$, $q_r$ with
  the same $K$.
- **Code.** `src/dynamics/turbulence.py`: `richardson`, `eddy_diffusivity`,
  `vertical_mixing`.
- **Parameters.** $Ri_c$ = 0.25, $K_{max}$ = 200 m^2/s, $l$ = 150 m. Passed
  through the `PrimitiveSigma` constructor; module globals set at run time do
  nothing (P-51).
- **Why.** The model develops real shear instability and had nothing to
  dissipate it: survival fell from 12/12 h to 2/12 h as the per-level shear
  rose from 3.2 to 8.4 m/s (docstring). $K_{max}$ 200, not 100: at 100 the
  ceiling truncated the diffusivity on 0.11 % of interfaces and cost two
  forecast hours over 4000 m terrain (6/12 against 8/12); above about 600 the
  parameter is inert; 200 is chosen over 150 for margin (P-40, measured
  2026-09-11/12). Reason for $l$ = 150 m not recorded.
- **Source.** Louis (1979) for the shape of the stability function; Miles
  (1961) and Howard (1961) for $Ri_c$ = 0.25.
- **Limits.** P-40: a higher ceiling does not reduce the overturning fraction,
  so the mechanism behind the two hours is open.

### 4.7 Dry convective adjustment

Wherever $\theta$ decreases with height, adjacent layers are mixed to their
mass-weighted mean. Pool-adjacent-violators, per column from the ground up,
with blocks $B$:

$$
w_k = \Delta\sigma_k\,\pi, \qquad \bar a_B = \frac{\sum_{k\in B} w_k a_k}{\sum_{k\in B} w_k}, \qquad
\text{merge } B_{upper} \text{ into } B_{lower} \text{ while } \bar\theta_{B_{upper}} < \bar\theta_{B_{lower}} - \epsilon
\tag{4.10}
$$

- $a \in \{\theta, u, v\}$ (and $q_v$, $q_c$, $q_r$ with moisture). In the NH
  core $\theta$ is $\theta_\rho$ (the core's own variable).
- **Code.** `src/dynamics/convection.py`: `dry_convective_adjustment_pav`
  (the model), `dry_convective_adjustment` (segment sweep, used for boundary
  frames, section 6.3). C kernel `pav_col` for the NH core (S5e).
- **Parameters.** $\epsilon$ = 1e-10 K; `mix_momentum=True`. Applied after
  each step (section 2.10).
- **Why.** Over 2500 m terrain the minimum Ri fell from 11.5 at hour 1 to
  -0.05 at hour 10: the mountain wave overturned. Diffusive mixing with K
  capped at 100 m^2/s relaxes a 600 m layer in an hour, slower than the wave
  steepens; "overturning has to be removed by rearrangement" (P-11).
  Pairwise mixing converged like a diffusion (P-12); segment mixing hit its
  20-sweep cap every step in calm cases and cost about 70 % of the step
  (P-63). PAV gives "the stable state closest to the input (mass-weighted
  least squares)" in one pass; default since test AC (P-63). Mixing $q$ over
  the same blocks was added on 2026-10-07 after the 12 km moist run's day
  dewpoint bias (P-74, S7b).
- **Source.** Manabe et al. (1965) for dry convective adjustment; Ayer et al.
  (1955) for pool-adjacent-violators.
- **Limits.** "This is the dry, hard-adjustment form: no entrainment, no
  cloud, no downdrafts, no mass flux" (docstring). The CAM design drops it
  eventually, since a non-hydrostatic model overturns unstable layers itself
  (CAM_DESIGN section 2); it is still applied in the 3 km runs.

### 4.8 Force-restore ground temperature

$$
\frac{dT_g}{dt} = C_G\,(R_n - H - LE) - \frac{2\pi}{\tau_d}\,(T_g - T_2), \qquad
C_G = \frac{2\sqrt\pi}{\rho_c\,d_1}, \qquad d_1 = \sqrt{k_s\tau_d}
\tag{4.11}
$$

$$
R_n = (1-A)\,SW_\downarrow + \varepsilon_g\left(\varepsilon_a\sigma_{SB}T_a^4 - \sigma_{SB}T_g^4\right), \qquad
SW_\downarrow = S_0\,\tau\,\max(\sin e, 0), \qquad T_a = \theta_1\Pi_1
\tag{4.12}
$$

$$
H = \rho\,c_p\,C_h\,U\,(\theta_g - \theta_1)\,\Pi_s, \quad
C_h = a^2 F_h, \quad U = \max(\lvert\mathbf V_1\rvert, U_{min}), \quad
\rho = \frac{p_s}{R_d T_g}, \quad \theta_g = T_g/\Pi_s, \quad LE = EF\,\max(R_n, 0)
\tag{4.13}
$$

Louis (1979) heat function, with the momentum function's tail on the stable
side (`long_tail=True`):

$$
F_h = \begin{cases}F_m(Ri) & Ri > 0\ \ (\text{long tail})\\[1ex] 1 - \dfrac{3b\,Ri}{1 + 3b\,c_h\,a^2\sqrt{-Ri\,z_1/z_0}} & Ri \le 0\end{cases}, \qquad
Ri = \frac{g z_1(\theta_1-\theta_g)}{\theta_1 U^2}
\tag{4.14}
$$

(Louis's own stable form, $1/(1 + 3b\,Ri\sqrt{1+d\,Ri})$, is in
`louis_heat` and is used when `long_tail=False`.) Heating of the lowest
layer:

$$
\left.\frac{\partial\theta}{\partial t}\right|_{N-1} = \frac{g\,H}{c_p\,\Delta p_1}\left(\frac{P_0}{p_1}\right)^{\kappa}, \qquad \Delta p_1 = \Delta\sigma_{N-1}\pi
\tag{4.15}
$$

- Over water $T_g = T_2$ (held). Initial $T_g = T_2 = T_1 + \Gamma z_1$ with
  $\Gamma$ = 0.0065 K/m, from the cycle-time lowest-level air temperature.
  Shortwave at mid-step. $\theta_g$ goes to the drag (Eqs. 4.3-4.6).
- **Code.** `src/dynamics/land_surface.py`: `ForceRestoreSurface.step`,
  `force_restore_coefficient`, `louis_heat`, `louis_momentum`;
  `PrimitiveSigma._step_land_surface`; `NHModel._step_land_surface`.
- **Parameters.** $A$ = 0.18, $\varepsilon_g$ = 0.95, $\varepsilon_a$ = 0.79,
  $\rho_c$ = 2.0e6 J m^-3 K^-1, $k_s$ = 5.0e-7 m^2/s, $\tau_d$ = 86 400 s
  ($d_1$ = 0.208 m, $C_G$ = 8.53e-6 K m^2/J), $EF$ = 0.5, $U_{min}$ = 1.0 m/s,
  $b = d = 5$, $c_h$ = 5.3, $\sigma_{SB}$ = 5.670374e-8 W m^-2 K^-4, $S_0$ and
  $\tau$ from section 4.10.
- **Why.** P-59 (no diurnal cycle), P-69 (night 10 m wind 1.6-1.8 times
  observed) and the night warm bias share one cause: "the ground has no
  temperature of its own". The user declined the prescribed flux as the
  default because every CAM takes its fluxes from a surface energy budget;
  this is "that budget, in its simplest established form". Test AK: with the
  similarity operator the 2 m T RMSE fell to 3.25 K (persistence 3.56) and
  the night wind ratio from 1.63 to 0.92. On by default since 2026-10-03
  (test AN) for the hydrostatic core and since S6 for the NH core. The long
  stable tail for the heat flux "is a choice made from that measurement. It
  is not taken from a reference": with Louis's heat function the night flux
  collapsed to -1.1 W/m^2 at $Ri_b$ 2.0, a runaway decoupling, where clear
  nights have tens of W/m^2.
- **Source.** Deardorff (1978) for force-restore; Louis (1979) for $F_m$,
  $F_h$ and their constants. The table of constants in the docstring gives a
  reason for each (mixed forest and cropland, moist loam, a Bowen ratio of 1)
  but names no published source except Brutsaert (1975) for $\varepsilon_a$.
- **Limits** (docstring): no clouds, one soil type, one albedo, the Great Lakes
  counted as land by the terrain mask, $EF$ a constant. At leads 19-24 the
  diurnal swing is too large with no clouds (test AM).

### 4.9 Clear-sky air emissivity

$$
\varepsilon_a = 1.24\left(\frac{e}{T}\right)^{1/7}, \qquad e \text{ in hPa}
\tag{4.16}
$$

- **Code.** `land_surface.brutsaert_emissivity`. The model does not call it;
  it uses the constant `EPS_A` = 0.79, which is this formula at $e$ = 12 hPa
  and $T$ = 285 K (1.24 (12/285)^(1/7) = 0.789).
- **Why.** "A late-September Northeast value; there are no clouds"
  (docstring).
- **Source.** Brutsaert (1975).

### 4.10 Solar geometry

$$
\delta = 23.45^\circ\,\sin\!\left(2\pi\,\frac{284 + n}{365}\right), \qquad
\omega_h = 15^\circ\left(t_{UTC} + \frac{\lambda}{15} - 12\right), \qquad
\sin e = \sin\varphi\sin\delta + \cos\varphi\cos\delta\cos\omega_h
\tag{4.17}
$$

- $n$ day of the year, $t_{UTC}$ in hours, $\lambda$ longitude in degrees.
- **Code.** `src/dynamics/diurnal.py`: `solar_declination`,
  `sin_solar_elevation`; used by the land surface and the prescribed flux.
- **Parameters.** $S_0$ = 1361 W/m^2, clear-sky transmission $\tau$ = 0.75.
- **Why.** The smallest step that puts the sun in "using no information from
  after the cycle time and no one else's model output" (P-59). $\tau$ = 0.75
  "because the dry model has no clouds".
- **Source.** Cooper (1969) for the declination; Kopp and Lean (2011) for
  $S_0$. Origin of $\tau$ = 0.75 not recorded.
- **Limits.** The equation of time is ignored (solar noon up to 16 min off).

### 4.11 Prescribed surface heat flux (option, off by default)

$$
H = f\,S_0\,\tau\,\max(\sin e, 0) + H_{night}\ \ \text{over land}, \qquad H = 0\ \ \text{over water}
\tag{4.18}
$$

applied through Eq. 4.15 before the convective adjustment.

- **Code.** `diurnal.surface_heat_flux`, `DiurnalHeating.theta_tendency`;
  `forecast.py --surface-heating` (hydrostatic core only).
- **Parameters.** $f$ = 0.2, $H_{night}$ = -30 W/m^2 (peaks near 120 W/m^2 at
  a late-September noon).
- **Why.** P-59 step 1: over 19 cycles it lowered 2 m T RMSE by 0.25 K. It
  stays off because of P-67 (with the instant momentum mixing, the night
  cooling destabilised a cycle) and because the user chose a surface energy
  budget instead (2026-10-02). Origin of $f$ and $H_{night}$ not recorded.

---

## 5. Moisture (warm rain, CAM stage S7)

`src/dynamics/moisture.py` and the moist parts of `nh3d.py`. Non-hydrostatic
core only (`forecast.py --moist` refuses the hydrostatic core).

### 5.1 Density potential temperature

$$
\theta_\rho = \theta\,\frac{1 + q_v/\varepsilon}{1 + q_v + q_c + q_r}, \qquad \Theta = \mu\,\theta_\rho, \qquad \varepsilon = R_d/R_v
\tag{5.1}
$$

- **Code.** `moisture.rho_factor`; `NHModel.enable_moisture`, `set_state`,
  `theta`, `_rho_f`. $\theta$ itself is never stored:
  $\theta = \theta_\rho/(\text{factor})$.
- **Parameters.** $\varepsilon$ = 0.622.
- **Why.** With $\theta_\rho$ in place of $\theta$ the core's equation of
  state (Eq. 3.2) is exact for moist air, since $p\alpha = R_d T_\rho$ with
  $\alpha$ the inverse density of the whole air; vapour buoyancy and
  condensate loading both enter through it, and "not one line of the dry
  solver changed" (S7 log entry, `moisture.py`).
- **Source.** Emanuel (1994), section 4.3 (cited in `moisture.py`).
- **Limits.** P-74: $\mu$ does not change when rain falls out or the surface
  evaporates (at most about 0.3 hPa in a day of heavy rain); mixing ratios
  are relative to dry air in the Kessler formulas but weighted by the
  whole-air $\mu$ (about 1 %); constant $c_p$ and $L_v$; no ice.

### 5.2 Saturation vapour pressure and mixing ratio

$$
e_s(T) = 611.2\,\exp\!\left(\frac{17.67\,(T - 273.15)}{T - 29.65}\right)\ \mathrm{Pa}, \qquad
q_s = \frac{\varepsilon\,e_s}{\max(p - e_s,\ 1)}
\tag{5.2}
$$

- **Code.** `moisture.es_bolton`, `moisture.qsat`; the same formula in
  `analysis/build.saturation_vapour_pressure`.
- **Why.** Part of the Kessler transcription from WRF (section 5.7). Test:
  $e_s$(20 C) = 2336.9 Pa.
- **Source.** Bolton (1980).

### 5.3 Vapour from relative humidity

$$
e = \frac{\mathrm{RH}}{100}\,e_s(T), \qquad q_v = \frac{\varepsilon\,e}{\max(p - e,\ 1)}
\tag{5.3}
$$

- RH clipped to 0-100 %. The analysis RH is interpolated to the sigma levels
  linearly in $\ln p$, held beyond the end levels (`forecast.rh_on_sigma`),
  and $T$ is the prepared initial (or frame) temperature.
- **Code.** `moisture.qv_from_rh`; `forecast.qv_from_rh_sigma`. Initial state
  and every analysis boundary frame. $q_c = q_r = 0$ at the start.

### 5.4 Dewpoint

$$
e = \frac{q_v\,p}{\varepsilon + q_v}\ (\ge 10^{-3}\ \mathrm{Pa}), \qquad
T_d = \frac{29.65\,\ln(e/611.2) - 17.67\cdot 273.15}{\ln(e/611.2) - 17.67}
\tag{5.4}
$$

- **Code.** `moisture.dewpoint` (the exact inverse of Eq. 5.2).

### 5.5 Transport with substep-mean fluxes

In each RK stage, with the acoustic substeps' mean fluxes

$$
\bar U = U^* + \frac1n\sum_{m=1}^{n} U''_m, \qquad \bar V = V^* + \frac1n\sum_{m=1}^{n} V''_m
\tag{5.5}
$$

$$
Q^{stage} = Q^n + n\Delta\tau\left[-\delta^+_x\mathcal F_3(\bar U, q^*) - \delta^+_y\mathcal F_3(\bar V, q^*) - \frac{\bar\Omega_{k+1}q^{w}_{k+1} - \bar\Omega_k q^{w}_k}{\Delta\sigma_k} + \mathcal P_Q\right], \qquad q^* = Q^*/\mu^*
\tag{5.6}
$$

- $Q = \mu q$ for each species; $\bar\Omega$ from Eqs. 3.3-3.4 with $\bar U$,
  $\bar V$; $\mathcal F_3$ and $q^w$ as Eqs. 3.10-3.11; $\mathcal P_Q$ the
  mixing tendency.
- **Code.** `NH3D._acoustic` (flux sums), `NH3D._adv_q`, `NH3D._q_stage`,
  `NH3D.step`; C kernels `adv_q`, `q_stage`.
- **Why.** "The fluxes are the mean of the acoustic substeps' fluxes. With
  those, the mass the Q update implies equals the mass the acoustic update
  gave mu, so a uniform q stays uniform to round-off. With the stage-start
  fluxes it would not" (S7). Test: uniform $q_v$ stays uniform to 1.3e-15
  over terrain in a noisy 3-D flow.
- **Source.** The third-order flux is Eq. 3.10. Advecting scalars with the
  mass fluxes averaged over the acoustic substeps is what WRF ARW does, for the
  same consistency reason (Skamarock and Klemp 2008). The S7
  entry in the research log does not cite it; the choice was made on the
  uniform-q requirement.

### 5.6 Hole filling

$$
\text{owe} = -\sum_k w_k\min(q_k,0), \quad \text{have} = \sum_k w_k\max(q_k,0), \quad
q_k \leftarrow \max(q_k,0)\left(1 - \frac{\text{owe}}{\text{have}}\right)\ \text{if have} > \text{owe, else } 0
\tag{5.7}
$$

- $w_k = \Delta\sigma_k$. Mass that cannot be repaid is counted in
  `NHModel.q_lost`.
- **Code.** `moisture.fill_negative`; called in `NHModel._step_moisture`
  before Kessler; C `moist_col`.
- **Why.** Removes negative mixing ratios "conserving the column mass".
- **Source.** Origin not recorded.

### 5.7 Kessler warm rain

One step per large step, after the dynamics and the land surface, on dry
$\theta$ and the core's non-hydrostatic $p$ (Eq. 3.2), $\rho = 1/\alpha$ and
$\Delta z = (\phi_k - \phi_{k+1})/g$.

Sedimentation (first, as in WRF):

$$
v_t = 36.34\,\left(10^{-3}\rho\,q_r\right)^{0.1364}\sqrt{\rho_{b}/\rho}, \qquad
n_f = \max\!\left(1,\ \left\lfloor 0.5 + \frac{\max_k(v_t\Delta t/\Delta z)}{0.75}\right\rfloor\right)
\tag{5.8}
$$

$$
q_{r,k} \leftarrow q_{r,k} - \frac{\Delta t}{n_f}\,\frac{\mathcal R_k - \mathcal R_{k-1}}{\rho_k\Delta z_k}, \qquad \mathcal R = \rho\,q_r\,v_t, \qquad \text{rain} \mathrel{+}= \mathcal R_{N-1}\frac{\Delta t}{n_f}\ \ (\mathrm{mm})
\tag{5.9}
$$

Autoconversion and accretion (accretion implicit in $q_r$):

$$
f_n = \frac{1}{1 + 2.2\,\Delta t\,q_r^{0.875}}, \qquad
\Delta q_{r} = q_c(1 - f_n) + f_n\,10^{-3}\,\Delta t\,\max(q_c - 10^{-3},\ 0)
\tag{5.10}
$$

Saturation adjustment (one Newton step) and rain evaporation:

$$
\Pi = (p/P_0)^\kappa,\ \ T = \Pi\theta,\ \ \gamma_L = \frac{L_v}{c_p\Pi},\ \ f_5 = \frac{17.67\,(273.15 - 29.65)\,L_v}{c_p},\ \ 
P = \frac{q_v - q_s}{1 + \dfrac{p}{p - e_s}\,\dfrac{q_s f_5}{(T - 29.65)^2}}
\tag{5.11}
$$

$$
E_r = \min\!\left[\Delta t\,\frac{(1.6 + 124.9\,(\rho_c q_r)^{0.2046})(\rho_c q_r)^{0.525}}{2.55\times10^{8}/(p\,q_s) + 5.4\times10^{5}}\,\frac{\max(q_s - q_v, 0)}{\rho_c\,q_s},\ \max(-P - q_c, 0),\ q_r\right], \qquad \rho_c = 10^{-3}\rho
\tag{5.12}
$$

$$
P' = \max(P, -q_c):\quad \theta \mathrel{+}= \gamma_L(P' - E_r),\quad q_v \leftarrow \max(q_v - P' + E_r, 0),\quad q_c \mathrel{+}= P',\quad q_r \mathrel{-}= E_r
\tag{5.13}
$$

- $\rho_b$ is the density of the lowest layer. $\Theta$ is rebuilt as
  $\mu\theta\,\cdot$(factor) from the result.
- **Code.** `moisture.kessler`, `moisture._vt`; `NHModel._step_moisture`;
  C `moist_col`.
- **Parameters.** $L_v$ = 2.5e6 J/kg; KC1 = 1e-3 s^-1, KC2 = 1e-3 kg/kg,
  KC3 = 2.2, KC4 = 0.875; fall Courant number at most 0.75 per substep;
  $c_p$ and $\kappa$ from `nh2d.py`. `warm_rain=False` keeps only the
  saturation adjustment (the reversible case).
- **Why.** Roadmap stage S4 (now S7): "warm rain (Kessler 1969)", with the
  Bryan and Fritsch (2002) benchmark as the gate. It is "transcribed from
  WRF's module_mp_kessler.F (the reference warm-rain scheme)". Tests: a 5 %
  supersaturated column adjusts to within 5e-4 of saturation in one step,
  conserving water and $c_pT + L_vq_v$ to 1e-13; 20 min of rain through 30
  layers closes the column water budget exactly.
- **Source.** Kessler (1969) as in Klemp and Wilhelmson (1978), via WRF's
  `module_mp_kessler.F` (WRF version not recorded). Terminal velocity: Klemp
  and Wilhelmson (1978). Rain evaporation: Klemp and Wilhelmson (1978),
  eq. 2.14 (as cited in `moisture.py`). Bolton (1980) for $e_s$.

**Gate (test only, `test_moist_bubble.py`).** The Bryan and Fritsch (2002)
moist/dry bubble pair, 20 x 10 km, $\Delta x$ = 200 m, 50 levels, 1000 s,
saturated moist base state with $q_t$ = 0.02 built with the model's own
saturation adjustment. Moist minus dry: max w +7.2 %, height of max
$\theta_\rho'$ +7.0 %, against a 10 % gate (S7).

### 5.8 Surface vapour flux

$$
E = \begin{cases}LE/L_v & \text{land}\\ \rho\,C_h\,U\,\left(q_s(T_w, p_s) - q_{v,1}\right) & \text{water}\end{cases}, \qquad
\Delta(\mu q_v)_{N-1} = \Delta t\,\frac{g\,E}{\Delta\sigma_{N-1}}
\tag{5.14}
$$

- $C_h$, $U$, $\rho$ as in Eq. 4.13. Negative $E$ over water (dew) is allowed.
- **Code.** `ForceRestoreSurface.step` (with `qv1`);
  `NHModel._step_land_surface`.
- **Why.** "Over land the latent heat LE that the energy budget already
  spends becomes vapour, E = LE / Lv, so the ground's energy balance is the
  one of S6 unchanged" (comment).
- **Limits.** S7b: the 12 km moist run's 2 m dewpoint lost to persistence in
  every cycle (RMSE 3.113 against 2.278 K), with a day bias of +1.14 K over
  hilly and +1.55 K over mountain stations. The record points at vertical
  transport of vapour out of the lowest layers; the PAV now mixes $q$
  (section 4.7), tested in arm nh12M2 (P-74).

### 5.9 Physics and boundaries with moisture

- Physics sees dry $\theta = \theta_\rho/f$ and the step-start $q$. Its
  $\theta$ tendency becomes a $\theta_\rho$ tendency at fixed $q$:
  $\mathcal P_\Theta \leftarrow f\,\mathcal P_\Theta$ (`NHModel._physics`).
- Edge relaxation acts on $\theta_\rho$, toward the frame's $\theta$ times
  the frame's factor $f(q_{frame})$, using the model's own $q$ for any
  species a frame lacks; $q$ is relaxed toward the frame's $q$
  (`NHModel.relax_with`, `_frame_theta_rho`, `_relax_q_zone`). Weights as
  in section 6.1.

---

## 6. Lateral boundaries and nesting

### 6.1 Davies relaxation

$$
d_{j,i} = \min\left(i,\ n_x - 1 - i,\ j,\ n_y - 1 - j\right), \quad
s = \mathrm{clip}\!\left(1 - \frac{d}{W},\ 0,\ 1\right), \quad
\alpha = \alpha_{max}\cdot\tfrac12\left(1 - \cos\pi s\right)
\tag{6.1}
$$

$$
\psi \leftarrow \psi + \alpha\,(\psi_{ext} - \psi), \qquad \psi \in \{u, v, \theta, \pi\}
\tag{6.2}
$$

applied to the state after every step, the same 2-D weights at every level.

- **Code.** `src/dynamics/boundaries.py`: `relaxation_weights`,
  `DaviesRelaxation`; `src/forecast.py`: `Relaxation3D.apply`,
  `run_forecast`.
- **Parameters.** $W$ = 15 cells, $\alpha_{max}$ = 0.1 (forecast defaults
  since 2026-09-26; module defaults 10 and 1.0). Profile cosine. Surface
  pressure: `pi_mode="full"` (point by point) by default; `"mean"` relaxes only
  the zone's weighted mean of $\pi$ and `"none"` leaves it free
  (`--relax-pi`; P-76 in `docs/PROBLEMS.md`, under test as E1).
- **Why.** "The smooth taper is the whole point: a hard switch from prescribed
  to free is itself a discontinuity, and reflects" (docstring). $\pi$ is
  relaxed because it is prognostic: leaving it free while relaxing the wind
  "drives the boundary column toward a mass field the incoming flow does not
  support" (comment). Width 15 and alpha 0.1: tests O, P and W (P-56). Pull
  strength toward the frozen state was the strongest control found (+5.4 h);
  with divergence damping every real case then ran 24 h. Applied to the state
  rather than as a tendency, "which is how Davies forcing is normally
  implemented and keeps it independent of the time scheme".
- **Source.** Davies (1976).
- **Limits.** P-70: per-step weights depend on the step length. The NH
  adapter corrects for it (section 6.2); the hydrostatic core does not, because its
  weights "were tuned at the step it runs at".

### 6.2 Relaxation per unit time (non-hydrostatic core)

$$
\alpha_{eff} = 1 - (1 - \alpha)^{\Delta t/\Delta t_{ref}}
\tag{6.3}
$$

$$
\phi \leftarrow \phi + \alpha_{eff}\,(\phi_{hyd} - \phi), \qquad W \leftarrow (1 - \alpha_{eff})\,W, \qquad
\phi_{hyd,k} = \phi_{hyd,k+1} + \frac{R_d\,\theta_k\Pi_k}{p_k}\,\mu\,\Delta\sigma_k,\ \ \phi_{hyd,N} = g h
\tag{6.4}
$$

- $u$, $v$, $\theta$ (or $\theta_\rho$), $q$ and $\mu$ are relaxed with
  $\alpha_{eff}$ as in Eq. 6.2, then $U$, $V$, $\Theta$, $Q$ rebuilt with the
  new $\mu$. $\Delta t_{ref}$ is the hydrostatic core's `max_dt()` for the
  run (16.7 s on the P-70 case). The same $\phi_{hyd}$ integral loads the
  initial state with $W = 0$ (`NHModel.set_state`, `hydrostatic_phi`).
- **Code.** `NHModel.relax_with`, `NHModel.relax_with_driver` (C kernels
  `relax_mu`, `relax_col`, zone columns only, S5d).
- **Why.** The Davies weights were set with the hydrostatic core's ~17 s
  step; at the NH core's 60 s step they relaxed about 3.5 times more weakly
  per hour. After the fix the NH-minus-hydrostatic difference fell from 0.38
  to 0.025 m/s rms in u at 1 h (P-70).

### 6.3 Boundary frames: preparation, time interpolation, stabilisation

$$
\psi_{ext}(t) = (1 - w)\,\psi_{a} + w\,\psi_{b}, \qquad w = \frac{t - t_a}{t_b - t_a}
\tag{6.5}
$$

- Each analysis frame goes through the same conversion, filter and balance
  as the initial state (sections 7.1-7.4), then the convective adjustment run to
  convergence (segment sweep, cap 1000 sweeps; zone columns only for nest
  frames). A run from observations has one frame, held for the whole run.
- **Code.** `boundaries.BoundaryDriver.at`, `bracket`;
  `forecast.stabilise_frame`, `forecast.main`.
- **Why.** Without time interpolation "the boundary jumps once an hour and
  radiates a spurious wave inward". Same preparation, because otherwise "the
  relaxation would drive one toward the other every step". Stabilisation:
  the relaxation re-imposed an unstable afternoon surface layer in the edge
  columns after every step, and the adjustment ran to its cap on every call
  (P-63, test Z2). One frame: "nothing observed after the cycle time may
  enter" (P-53).

### 6.4 Nest: horizontal remapping

$$
x_{src} = \frac{i + \tfrac12\,[\text{c}]}{n_{dst}}\,n_{src} - \tfrac12\,[\text{c}], \qquad
a_{dst} = (1-w_x)\left[(1-w_y)a_{j_0,i_0} + w_y a_{j_1,i_0}\right] + w_x\left[(1-w_y)a_{j_0,i_1} + w_y a_{j_1,i_1}\right]
\tag{6.6}
$$

- The bracketed $\tfrac12$ applies to cell-centred ("c") points, not to face
  ("f") points; indices clipped to the source grid.
- **Code.** `src/nest.py`: `fractional_index`, `regrid`.
- **Why.** Both grids are regular in latitude and longitude over the same
  domain, so a 3 km point's position in the 12 km grid is a fractional index,
  separately for each C-grid staggering (docstring).

### 6.5 Nest: surface pressure and vertical remapping

$$
p_{s,dst} = p_{s,src}\,\exp\!\left(-\frac{g\,\Delta h}{R_d\,\bar T}\right), \qquad \bar T = T_{low} - \tfrac12\Gamma\,\Delta h, \qquad \Delta h = h_{dst} - h_{src}
\tag{6.7}
$$

$$
\ln\theta,\ u,\ v,\ q\ \text{linear in}\ \ln p\ \text{between}\ p_{src,k} = p_t + \sigma^{src}_k\pi_{src}\ \text{and}\ p_{dst,m} = p_t + \sigma^{dst}_m\pi_{dst}
\tag{6.8}
$$

- Held beyond the end levels, except $\ln\theta$ extrapolated above the top
  level. $q \ge 0$ enforced. Terrain at u and v points is the two-cell mean.
- **Code.** `nest.surface_pressure_to`, `nest.remap_column`,
  `nest.frame_from_coarse`, `nest.scalar_from_coarse`,
  `nest.frames_from_forecast`.
- **Parameters.** $\Gamma$ = 0.0065 K/m; $T_{low}$ the 12 km lowest-level
  temperature; constants from `nest.py` ($R_d$ = 287.04).
- **Why.** $\ln\theta$ "is exactly linear in ln p for an isothermal layer".
  Extrapolating theta at the top because "the 3 km top level can sit above
  the 12 km one, and a constant theta there would be a neutral layer under
  the lid". Held $q$ "so a non-negative field stays non-negative". Edges
  follow the 12 km forecast hour by hour by the user's choice (prompt 172).
- **Source.** Hypsometric equation (textbook).

---

## 7. Initialisation and the observation analysis

### 7.1 Surface pressure from analysis heights

$$
\ln p_s = \ln p_k + \frac{h - z_k}{z_{k+1} - z_k}\,(\ln p_{k+1} - \ln p_k)\quad (z_k \le h < z_{k+1})
\tag{7.1}
$$

$$
\ln p_s = \ln p_0 + \frac{g\,(z_0 - h)}{R_d\,T_{layer}}, \qquad T_{layer} = \frac{g\,(z_1 - z_0)}{R_d\ln(p_0/p_1)}\quad (h < z_0)
\tag{7.2}
$$

- $z_k$ analysis geopotential height on pressure level $p_k$ (index 0 the
  highest pressure).
- **Code.** `src/dynamics/interpolate.py`, `surface_pressure_from_heights`.
- **Why.** Interpolating against height "inherits the analysis's own
  stratification" rather than assuming a lapse rate. 1000 hPa sits about
  100 m above sea level, so most of a coastal domain needs the hydrostatic
  extension. A reversed bracket once gave a 253 hPa error over 2500 m terrain,
  caught against a standard atmosphere (comment).

### 7.2 Pressure levels to sigma levels

$$
a(p) = a_{k} + \frac{\ln p - \ln p_{k}}{\ln p_{k+1} - \ln p_{k}}\,(a_{k+1} - a_{k})
\tag{7.3}
$$

- Targets $p = p_t + \sigma_k\pi$ with $\pi$ from Eqs. 7.1-7.2. Winds held
  outside the source range; $\theta$ (from $T$ by Eq. 1.3) extrapolated linearly
  in $\ln p$.
- **Code.** `interpolate.interp_log_p`, `interpolate.pressure_to_sigma`.
- **Why.** Fields are far more nearly linear in $\ln p$ than in $p$; with
  25 hPa gaps near the ground and 50-100 hPa aloft "the difference is not
  cosmetic". Holding $\theta$ constant below the lowest level makes the
  near-surface layer exactly neutral, which the convective adjustment reads
  as on the edge of instability; extrapolating the wind's shear downward
  makes surface winds the drag then fights (docstring).
- **Note.** `interpolate.hydrostatic_theta` (an exact inverse of Eq. 2.6,
  written for P-46) is defined but not called by the model.

### 7.3 Spectral low-pass filter

$$
r(\lambda) = \begin{cases}0 & \lambda \le 4\Delta\\ \tfrac12\left[1 - \cos\!\left(\pi\,\dfrac{\lambda - 4}{8 - 4}\right)\right] & 4\Delta < \lambda < 8\Delta\\ 1 & \lambda \ge 8\Delta\end{cases}, \qquad
\hat a \to r\,\hat a
\tag{7.4}
$$

- $\lambda = 1/\sqrt{k_x^2 + k_y^2}$ in grid cells. Applied to $u$, $v$ and
  $\theta - \langle\theta\rangle_k$.
- **Code.** `src/dynamics/initialization.py`: `spectral_lowpass`,
  `filter_initial_state`; `forecast.main`.
- **Why.** P-10: white noise put 89 % of its variance at scales the grid
  cannot carry, and advection amplified the 2$\Delta x$ amplitude (doubling in
  about 20 min) faster than the 3 h hyperdiffusion removed it; 0.30 m/s of
  noise survived 12 h, 0.60 m/s did not. A raised cosine because "a sharp
  spectral truncation rings in physical space". Theta filtered with the wind
  so the pressure-gradient force does not put the noise back. Order measured:
  none 1/12 h, filter only 11/12 h, filter then rebalance 12/12 h. Reason for
  the exact 4$\Delta x$ and 8$\Delta x$ cut-offs not recorded.

### 7.4 Removing initial divergence

$$
\nabla^2_d\chi = D, \qquad \hat\chi = \frac{\hat D}{-4\sin^2(k_x/2)/\Delta x^2 - 4\sin^2(k_y/2)/\Delta y^2}, \qquad
u \leftarrow u - \delta^-_x\chi, \quad v \leftarrow v - \delta^-_y\chi
\tag{7.5}
$$

then 60 passes of

$$
u \leftarrow u + c\,\delta^-_x D, \qquad v \leftarrow v + c\,\delta^-_y D, \qquad c = 0.2\,\min(\Delta x,\Delta y)^2
\tag{7.6}
$$

- Per level; the mean mode is untouched. $D = \delta^+_x u + \delta^+_y v$.
- **Code.** `subgrid.remove_divergence_spectral`,
  `subgrid.balance_initial_state` (`method="hybrid"`).
- **Why.** "Analysis winds are not in balance with OUR discretisation": they
  carry divergence of order 1e-3 s^-1 against about 1e-5 in the atmosphere,
  tens of Pa/s of vertical motion. Using the eigenvalues of the discrete
  Laplacian makes the cancellation exact (P-15: implied omega 60 to
  1.08 Pa/s). The FFT assumes periodicity, so the 60-pass polish cleans the
  outermost cells. Reasons for $c$ and the pass count not recorded.
- **Source.** Helmholtz decomposition (textbook).

### 7.5 Observation analysis: successive correction of innovations

$$
\mathrm{inc}(\mathbf x) = \sum_{p=1}^{P}\frac{\sum_i w_i^{(p)}(\mathbf x)\,\big(d_i - \mathrm{inc}^{(p-1)}_i\big)}{\sum_i w_i^{(p)}(\mathbf x) + \lambda}, \qquad
w_i^{(p)} = \exp\!\left(-\frac{r_i^2}{L_p^2}\right)\left(\frac{s_{ref}}{s_i}\right)^2, \qquad L_{p+1} = L_p\sqrt{\gamma}
\tag{7.7}
$$

$$
r^2 = \left[(\varphi_a - \varphi_b)\,111.132\right]^2 + \left[(\lambda_a - \lambda_b)\,111.320\cos\tfrac{\varphi_a+\varphi_b}{2}\right]^2\ \ (\mathrm{km}^2)
\tag{7.8}
$$

Gap-filling variant for the 3 km analysis:

$$
w_i \leftarrow w_i\,\exp\!\left[-\left(\frac{z_{\mathbf x} - z_i}{H}\right)^2\right]\times\begin{cases}c_{coast} & \text{land/water mismatch}\\ 1 & \text{otherwise}\end{cases}
\tag{7.9}
$$

- $d_i = y_i - H(x_b)$ the innovation against the first guess, $s_i$ the
  observation error, $s_{ref}$ their median, $\mathrm{inc}^{(p-1)}_i$ the
  increment already analysed at observation $i$.
- **Code.** `src/analysis/barnes.py`: `barnes_increments`,
  `barnes_increments_aware`, `sq_distance_km`; called from `build.py`.
- **Parameters.**

| analysis | L (km) | lambda | passes | where |
|---|---|---|---|---|
| upper air, per level (T, RH, u, v) | 450 | 0.25 | 2 | `build.analyse_upper` |
| sea-level pressure (background: median of the obs) | 300 | 0.1 | 3 | `build.analyse_pmsl` |
| surface T, RH, u, v | 120 | 0.25 | 2 (3 with gap-filling) | `build.blend_surface` |

  $\gamma$ = 0.35 (so 120, 71, 42 km over three passes). Gap-filling:
  $H$ = 300 m, $c_{coast}$ = 0.5 (`ingest_obs.py --gapfill`, 3 km only).
- **Why.** Analysing innovations with $\lambda$ in the denominator keeps the
  first guess far from every observation; classical normalised Barnes would
  give a point 800 km away the nearest sounding's value at full strength.
  "lam is the ratio of observation to first-guess error variance in
  disguise". The error factor makes a buoy (1 K) count more than an ASOS
  thermometer (1.5 K) (docstring). Gap-filling (S5b, constants chosen before
  the run): "a valley station does not set a ridge"; 3 km gap-filling scored
  1.541 K on withheld stations against 1.641 K for the 12 km analysis.
  Reasons for the length scales, the $\lambda$ values, the pass counts and
  $\gamma$ = 0.35 not recorded.
- **Source.** Barnes (1964); Koch et al. (1983) for the multi-pass
  convergence parameter $\gamma$.

### 7.6 Analysis steps

Soundings to analysis levels (superobs): linear in $\ln p$ between reported
levels, skipped where the bracketing levels are more than 150 hPa apart.
First guess: the previous forecast at the lead (tolerance 0.02 h, P-64);
else the previous analysis when there are no soundings; else the mean
sounding per level; else a standard atmosphere:

$$
T_{std}(p) = \max\left(288.15\left(\frac{p}{101\,325}\right)^{R_d\Gamma/g},\ 216.65\right), \quad \mathrm{RH} = 50\ \%, \quad u = v = 0
\tag{7.10}
$$

Heights, hydrostatically from the analysed sea-level pressure:

$$
T_v = T(1 + 0.608\,q), \qquad q = \frac{\varepsilon\,e}{p - (1-\varepsilon)e}, \qquad e = \frac{\mathrm{RH}}{100}e_s(T)
\tag{7.11}
$$

$$
z_0 = \frac{R_d\,\bar T_v}{g}\ln\frac{p_{msl}}{p_0},\ \ \bar T_v = T_{v,0} + \tfrac12\Gamma z_0\ (\text{two passes}); \qquad
z_k = z_{k-1} + \frac{R_d}{g}\,\tfrac12(T_{v,k-1} + T_{v,k})\ln\frac{p_{k-1}}{p_k}
\tag{7.12}
$$

Surface blend:

$$
y_T \leftarrow y_T + \Gamma\,(z_{station} - z_{grid}), \qquad
a_k \leftarrow a_k + \mathrm{clip}\!\left(1 - \frac{Z_k - h}{1000\ \mathrm{m}},\ 0,\ 1\right)\mathrm{inc}_{sfc}
\tag{7.13}
$$

- The first guess is evaluated at 2 m (T, RH) or 10 m (u, v) above the grid
  terrain by `value_at_height` (linear in z; below the lowest level
  $T = T_0 + \Gamma(Z_0 - z)$, others held). Stations more than 600 m from
  the grid terrain are not used. RH clipped to 1-100 %. Heights are computed
  again after the blend. One ASOS station in five is withheld by
  $\mathrm{crc32}(\text{id}) \bmod 5 = 0$. Forecast first guesses below the
  ground use $T = T_{N-1}(p/p_{N-1})^{R_d\Gamma/g}$. The 3 km first guess is
  the 12 km analysis regridded bilinearly (`regrid_features`).
- **Code.** `src/analysis/build.py`: `superob_soundings`,
  `first_guess_from_*`, `standard_temperature`, `virtual_temperature`,
  `hydrostatic_heights`, `value_at_height`, `blend_surface`,
  `sigma_to_pressure_levels`, `is_withheld`, `build_analysis`.
- **Parameters.** $R_d$ = 287.05, $\Gamma$ = 0.0065 K/m, $\varepsilon$ = 0.622.
- **Why.** Heights "never analysed independently: heights perturbed
  separately from temperature were one of the four worst test-setup defects
  in this project". Sea-level pressure is "the anchor for every height".
  Superob gap rule: "a profile with a hole is not evidence about what is in
  the hole". Withheld stations: "any other score of a cycle-time analysis is
  a score against its own inputs" (docstring). Reasons for the 1000 m blend
  depth, the 600 m mismatch limit and the 150 hPa gap not recorded.
- **Source.** Bolton (1980) for $e_s$; virtual temperature and the
  hypsometric equation are textbook; the US standard atmosphere is named in
  the docstring without a citation.
- **Option, off.** `extrapolate_below_ground` replaces below-ground
  pressure-level values (T along $\Gamma$ from the first level above
  ground). Tested as the cause of P-56 and refuted, so it is off.

### 7.7 Relative humidity from radiosonde dewpoint

$$
\mathrm{RH} = 100\,\frac{\exp\!\left(a\,T_d/(b + T_d)\right)}{\exp\!\left(a\,T/(b + T)\right)}, \qquad a = 17.625,\ b = 243.04\ {}^\circ\mathrm{C}
\tag{7.14}
$$

- **Code.** `src/verification/fetchers.py`, `rh_from_dewpoint`.
- **Source.** Alduchov and Eskridge (1996).
- **Note.** Observations use this Magnus form while the model and analysis
  use Bolton (Eq. 5.2). Reason for the two forms not recorded.

---

## 8. Diagnostics and verification

### 8.1 Map diagnostics

Column state uses Eqs. 1.1 and 1.3 and the hydrostatic integral Eq. 2.6. Then:

$$
T_g = T_{N-1} + \Gamma\,(Z_{N-1} - h), \qquad p_{msl} = p_s\left(\frac{T_g + \Gamma h}{T_g}\right)^{g/(R_d\Gamma)}
\tag{8.1}
$$

$$
T(p) = T_{N-1}\left(\frac{p}{p_{N-1}}\right)^{R_d\Gamma/g}, \qquad
Z(p) = Z_{N-1} + \frac{T_{N-1}}{\Gamma}\left[1 - \left(\frac{p}{p_{N-1}}\right)^{R_d\Gamma/g}\right]\qquad (p > p_{N-1})
\tag{8.2}
$$

$$
\eta = \frac{\partial v}{\partial x} - \frac{\partial u}{\partial y} + 2\Omega\sin\varphi, \qquad
\omega(p) = -\int_{200\ \mathrm{hPa}}^{p}\nabla\!\cdot\!\mathbf V\,dp', \quad \omega(200\ \mathrm{hPa}) = 0
\tag{8.3}
$$

- Between model levels fields are linear in $\ln p$. Winds are de-staggered
  to centres, $u_c = \tfrac12(u_i + u_{i+1})$, $v_c = \tfrac12(v_j + v_{j+1})$.
  $\eta$ by centred differences (`np.gradient`); $\omega$ by the trapezoid
  rule at 25 hPa spacing. Points below the ground are flagged.
- **Code.** `src/maps/derive.py`: `column_state`, `ground_temperature`,
  `mslp`, `to_pressure`, `relative_vorticity`, `kinematic_omega`,
  `destagger`, `snapshot`.
- **Parameters.** $\Gamma$ = 0.0065 K/m; $\Omega$ = 7.292e-5 s^-1;
  constants from `sigma.py`.
- **Why.** "Everything a map product needs is derived here from the model's
  own variables with the model's own constants." 1000 hPa is extrapolated
  over high ground "as on every operational chart". Near-surface map fields
  are the lowest model level, and "its maps say so instead of pretending"
  (docstring).
- **Source.** Standard-lapse-rate reduction and the kinematic method are
  textbook; origin not recorded.

### 8.2 Observation operator

Horizontal: bilinear in the cell-centre positions
$f_x = x/\Delta x - \tfrac12$, $f_y = y/\Delta y - \tfrac12$, with
$x = (\lambda - \lambda_{min})\cdot 111\,320\cos\varphi_0$ and
$y = (\varphi - \varphi_{min})\cdot 111\,132$; no extrapolation outside the
domain. Vertical: linear in $\ln p$ in the column's own pressures (Eq. 1.1);
`None` outside the column. Surface observations use index N-1, never index 0.

- **Code.** `src/verification/obs_operator.py`, `GridInterpolator`;
  `sigma_operator.SigmaInterpolator`.
- **Why.** In sigma there is no single 1-D array of level pressures. The base
  class's surface fallback `field3d[0]` is the model lid here: "a 2 m
  thermometer would have been scored against the 200 hPa field ... off by
  seventy kelvin" (docstring). Extrapolated matches are "a fabricated number".

### 8.3 Station temperature, standard operator

$$
z_1^{asl} = h + \frac{R_d T_1}{g}\ln\frac{p_s}{p_1}, \qquad T_{station} = T_1 - \Gamma\left(z_{station} - z_1^{asl}\right)
\tag{8.4}
$$

- **Code.** `SigmaInterpolator.lowest_level_height`, `station_temperature`;
  `obs_operator.elevation_correct_temperature`.
- **Parameters.** $\Gamma$ = 6.5e-3 K/m; $R_d$ = 287.05.
- **Why.** A valley station and its grid cell can differ by hundreds of
  metres, "which is several degrees" (docstring).
- **Limits.** "The standard lapse rate is wrong exactly when it matters most
  -- a clear calm night with a surface inversion" (docstring; P-68).
  Production uses the similarity operator (section 8.5) whenever the forecast carries
  a ground temperature.

### 8.4 10 m wind, neutral log law

$$
f_{10} = \frac{\ln(10/z_0)}{\ln(z_{agl}/z_0)}\ \ (z_{agl} > 10\ \mathrm m), \qquad f_{10} = 1\ \ (z_{agl} \le 10\ \mathrm m), \qquad z_{agl} = z_1^{asl} - h
\tag{8.5}
$$

- **Code.** `SigmaInterpolator.wind_10m_factor`; `verify.py
  --wind-operator log10m` (the default since 2026-10-02, the user's decision).
- **Parameters.** $z_0$ = 0.1 m.
- **Why.** The lowest level is about 300 m up; without a reduction "every
  wind score is mostly the speed ratio". It is "the same profile the model's
  surface drag assumes" (docstring). Test AI: persistence lead-1 speed ratio
  1.86 to 1.10 (P-68).
- **Note.** The operator's $z_0$ = 0.1 m predates the land/sea roughness map
  of test AM (1.0 m over land). Reason for keeping 0.1 in the operator not
  recorded.

### 8.5 Similarity operator (10 m wind and 2 m temperature)

$$
Ri_b = \zeta\,\frac{\ln(z_1/z_{0h}) - \psi_h(\zeta) + \psi_h(\zeta z_{0h}/z_1)}{\left[\ln(z_1/z_0) - \psi_m(\zeta) + \psi_m(\zeta z_0/z_1)\right]^2}, \qquad
Ri_b = \frac{g z_1(\theta_1 - \theta_g)}{\theta_1 U_1^2}
\tag{8.6}
$$

$$
\frac{U_{10}}{U_1} = \frac{\ln(10/z_0) - \psi_m(10\zeta/z_1) + \psi_m(\zeta z_0/z_1)}{\ln(z_1/z_0) - \psi_m(\zeta) + \psi_m(\zeta z_0/z_1)}, \qquad
\theta_2 = \theta_g + (\theta_1 - \theta_g)\,\frac{\ln(2/z_{0h}) - \psi_h(2\zeta/z_1) + \psi_h(\zeta z_{0h}/z_1)}{\ln(z_1/z_{0h}) - \psi_h(\zeta) + \psi_h(\zeta z_{0h}/z_1)}
\tag{8.7}
$$

Unstable ($\zeta < 0$), $x = (1 - 16\zeta)^{1/4}$:

$$
\psi_m = 2\ln\frac{1+x}{2} + \ln\frac{1+x^2}{2} - 2\arctan x + \frac\pi2, \qquad \psi_h = 2\ln\frac{1+x^2}{2}
\tag{8.8}
$$

Stable ($\zeta \ge 0$):

$$
\psi_m = -\left[a\zeta + b\left(\zeta - \frac cd\right)e^{-d\zeta} + \frac{bc}{d}\right], \qquad
\psi_h = -\left[\left(1 + \frac{2a\zeta}{3}\right)^{3/2} + b\left(\zeta - \frac cd\right)e^{-d\zeta} + \frac{bc}{d} - 1\right]
\tag{8.9}
$$

- $\zeta = z_1/L$ found by bisection on Eq. 8.6 in $[-50, 500]$, 80 iterations.
  $T_2 = \theta_2\Pi_s$ at the model ground, then $-\Gamma(z_{station} - h)$
  to the station.
- **Code.** `src/verification/surface_similarity.py`: `psi_m`, `psi_h`,
  `solve_zeta`, `surface_diagnostics`; `SigmaInterpolator.surface_similarity`,
  `station_temperature_similarity`; `verify.resolve_surface_operator`.
- **Parameters.** $a$ = 1, $b$ = 2/3, $c$ = 5, $d$ = 0.35;
  $z_0$ = 0.1 m, $z_{0h} = z_0/10$; $U_1 \ge$ 0.5 m/s; $z_1 \ge$ 20 m;
  $\kappa_{vK}$ = 0.4. Used for every forecast with a ground temperature, not
  for persistence (`verify.py --surface-operator auto`, default since
  2026-10-03).
- **Why.** A clear night is not neutral: under the neutral operator the night
  10 m wind was 1.6-1.8 times the observed (P-69). With similarity the night
  ratio fell from 1.63 to 0.92 (test AK). Persistence keeps the standard
  operator, because holding the cycle-time stability all day reduces its
  daytime wind as if it were still night (`resolve_surface_operator`).
  Reasons for $z_{0h} = z_0/10$ and the floors not recorded; $z_0$ as in 8.4.
- **Source.** Businger-Dyer forms integrated by Paulson (1970) for Eq. 8.8
  (the Businger-Dyer origin itself is not cited further in the repo);
  Beljaars and Holtslag (1991) for Eq. 8.9.
- **Limits.** Similarity holds in the surface layer; the lowest level at
  about 230 m "is often above the stable boundary layer" at night, so the
  diagnostic is an extrapolation (docstring).

### 8.6 Station dewpoint

$$
p_{station} = p_s\exp\!\left(-\frac{g\,(z_{station} - h)}{R_d\,T_1}\right), \qquad T_{d,station} = T_d\!\left(\max(q_{v,N-1}, 10^{-6}),\ p_{station}\right)
\tag{8.10}
$$

- $T_d(q, p)$ is Eq. 5.4. $q_v$ of the lowest level is taken as constant
  through the surface layer.
- **Code.** `SigmaInterpolator.station_dewpoint`; `verify.py` variable
  `DPT` (ASOS `dwpf`).
- **Parameters.** $R_d$ = 287.04 and $g$ = 9.80665 written as literals.

### 8.7 Scores

$$
\mathrm{bias} = \frac1n\sum(f - o), \quad
\mathrm{RMSE} = \sqrt{\frac1n\sum(f - o)^2}, \quad
\mathrm{MAE} = \frac1n\sum\lvert f - o\rvert, \quad
r = \mathrm{corr}(f, o), \quad
\mathrm{skill} = 1 - \frac{\mathrm{RMSE}}{\mathrm{RMSE}_{ref}}
\tag{8.11}
$$

- Non-finite pairs dropped; $r$ only when $n > 1$ and both standard
  deviations are non-zero. The reference is persistence: the prepared
  initial state written at every output time (`forecast.py --persistence`,
  P-59).
- **Code.** `src/verification/scoring.py`: `scores`, `scores_by`,
  `skill_vs_reference`.
- **Why.** "Beating nothing is not a result" (docstring).

### 8.8 Precipitation

The moist core writes accumulated rain (mm) from the sedimentation flux (Eq. 5.9)
at each output time. No precipitation scoring equation was found in `src/` or
`tools/` when this file was written: MRMS products are fetched
(`fetchers.mrms_url`) but not matched to the forecast. The CAM roadmap names
the fractions skill score (Roberts and Lean 2008) for a later stage. This
section should be updated when the scoring is written.

---

## 9. Physical constants

Values exactly as in the code.

| Quantity | Value | Where defined |
|---|---|---|
| Gas constant of dry air $R_d$ | 287.058 J/(kg K) | `dynamics/sigma.py`, `dynamics/vertical.py`: hydrostatic core, physics, land surface, maps, initial conversion |
| | 287.04 | `dynamics/nh2d.py` (imported by `nh2d_mass.py`, `nh3d.py`: the NH equation of state, reference state, hydrostatic phi, Kessler); `nest.py`; literal in `sigma_operator.station_dewpoint` |
| | 287.05 | `analysis/build.py`; `verification/sigma_operator.py` |
| Specific heat $c_p$ | 1004.6 J/(kg K) | `sigma.py`, `vertical.py`; literal in `build.first_guess_from_forecast` and in `sigma_operator.py`'s KAPPA |
| | 1004.5 | `nh2d.py`, `nest.py` |
| $c_v = c_p - R_d$ | 717.46 | `nh2d.py` |
| $\gamma = c_p/c_v$ | 1.40008 | `nh2d_mass.py`, `nh3d.py` |
| $\kappa = R_d/c_p$ | 0.285744 | `sigma.py` |
| | 0.285754 | `nh2d.py`, `nest.py` |
| | 0.285736 | `sigma_operator.py` (287.05/1004.6) |
| Gravity $g$ | 9.80665 m/s^2 | every module (`G0` or `G`) |
| Reference pressure $P_0$ | 100 000 Pa | every module |
| Lid pressure $p_t$ | 20 000 Pa | `sigma.SigmaLevels` |
| $\varepsilon = R_d/R_v$ | 0.622 | `moisture.py`, `build.py` |
| Latent heat $L_v$ | 2.5e6 J/kg | `moisture.py` |
| Bolton constants | 611.2 Pa, 17.67, 29.65 K, 273.15 K | `moisture.py`, `build.py` |
| Kessler constants | 1e-3 s^-1, 1e-3 kg/kg, 2.2, 0.875; 36.34, 0.1364; 1.6, 124.9, 0.2046, 0.525, 2.55e8, 5.4e5 | `moisture.py` (and `nh3d_kernels.c`) |
| Virtual-temperature factor | 0.608 | `build.virtual_temperature` |
| Lapse rate $\Gamma$ | 0.0065 K/m | `build.py`, `land_surface.py`, `nest.py`, `maps/derive.py`; 6.5e-3 in `obs_operator.py` |
| von Karman $\kappa_{vK}$ | 0.4 | `surface.py`, `surface_similarity.py`; literal in `land_surface.py` |
| Stefan-Boltzmann $\sigma_{SB}$ | 5.670374e-8 W m^-2 K^-4 | `land_surface.py` |
| Solar constant $S_0$ | 1361 W/m^2 | `diurnal.py` |
| Earth rotation $\Omega$ | 7.2921e-5 s^-1 | `forecast.build_grid` |
| | 7.292e-5 s^-1 | `maps/derive.py` |
| Earth radius $a$ | 6 371 000 m | `forecast.build_grid` (beta), `maps/geography.py` |
| Metres per degree | 111 132 (lat), 111 320 (lon, times cos) | `analysis/geo.py`, `verification/obs_operator.py`, `forecast.build_grid`; km values 111.132, 111.320 in `analysis/barnes.py` |
| Sound speed bound | 350 m/s | `nh3d.NHModel.C_SOUND` |
| Magnus constants | 17.625, 243.04 C | `verification/fetchers.py` |
| US standard atmosphere | 288.15 K, 101 325 Pa, 216.65 K | `build.standard_temperature` |

**Where they differ.** Three values of $R_d$ (287.058, 287.05, 287.04) and two
of $c_p$ (1004.6, 1004.5) are in use. The hydrostatic core and the physics use
287.058 and 1004.6; the non-hydrostatic dynamics, Kessler and the nest use
287.04 and 1004.5; the analysis and the verification operator use 287.05.
The NH core is loaded with a state converted using the `sigma.py` values and
then evaluated with the `nh2d.py` values. The differences are below 0.01 %
in $\kappa$. The reason for the different values is not recorded, and no
entry in `docs/PROBLEMS.md` covers it.

---

## 10. Gaps in the record

Every place above that says the origin or reason is not recorded. Numbers in
parentheses are section numbers.

1. Power-law stretch of the sigma levels (1.2): origin not recorded.
2. Metres-per-degree constants 111 132 and 111 320 (1.5): origin not recorded.
3. Beta plane rather than latitude-dependent f (1.6): reason not recorded.
4. 1-2-1 terrain filter (1.7): origin not recorded.
5. Flux-form vertical advection with the $a\,\partial_\sigma\dot\sigma$
   correction (2.5): origin not recorded.
6. Reference-state pressure-gradient form (2.4): origin not recorded; the
   record is mixed on whether it helps.
7. Hyperdiffusion damping time of 3 h (2.7): origin not recorded ("a common
   operational choice").
8. Sponge (2.9): origin not recorded; rate 1/900 s "chosen, not derived";
   reason for the raised-cosine profile not recorded.
9. CFL safety factor 0.6 (2.11): reason not recorded.
10. NH reference sounding (3.7): origin not recorded; reason for
    $z_{top}$ = 20 000 m not recorded.
11. NH off-centring 0.1 and acoustic divergence damping 0.1 (3.8): reason
    not recorded.
12. NH w-damping factor 2 (3.10): reason not recorded.
13. NH 60 s step cap and 0.7 advective factor (3.10): reason not recorded.
14. Cut-off drag stability constants 10, 4 and 0.2 (4.4): origin not
    recorded.
15. Roughness lengths (4.2) and land-surface constants other than
    $\varepsilon_a$ (4.8): reasons given, published source not named.
16. Mixing length 150 m (4.6): reason not recorded.
17. Clear-sky transmission 0.75 (4.10): origin not recorded.
18. Prescribed-flux constants f = 0.2 and H_night = -30 W/m^2 (4.11): origin
    not recorded.
19. Substep-mean flux transport of moisture (5.5): origin not recorded.
20. Hole filling (5.6): origin not recorded.
21. WRF version of `module_mp_kessler.F` (5.7): not recorded.
22. Filter cut-offs 4 and 8 grid lengths (7.3): reason not recorded.
23. Divergence polish coefficient 0.2 min(dx,dy)^2 and 60 passes (7.4):
    reason not recorded.
24. Barnes length scales, lambda values, pass counts and gamma = 0.35 (7.5):
    reason not recorded.
25. Surface blend depth 1000 m, station mismatch limit 600 m and superob gap
    150 hPa (7.6): reason not recorded.
26. Magnus form for observations against Bolton in the model (7.7): reason
    not recorded.
27. MSLP reduction and kinematic omega (8.1): origin not recorded (textbook).
28. $z_0$ = 0.1 m in the verification operators against the 1.0/0.0002 m map
    in the model (8.4, 8.5): reason not recorded.
29. $z_{0h} = z_0/10$ and the floors $U_1 \ge$ 0.5 m/s, $z_1 \ge$ 20 m (8.5):
    reason not recorded.
30. Three values of $R_d$ and two of $c_p$ (9): reason not recorded.

---

## 11. Equation sources

One row per numbered equation, in order. **Source(s)** lists the works named
for that equation in its section's **Source** line, or for a quantity it takes
from another numbered equation (noted in the last column). "textbook (none
cited)" means the section calls the equation textbook physics; "none stated"
means the section names no source. **How the source was used** is one of: as
published; adapted (what changed); constants from a source; project's own
(derived here); origin not recorded. The reasons are in the equation's own
section.

| Eq. | What it is | Source(s) | How the source was used |
|---|---|---|---|
| 1.1 | Sigma coordinate and p(sigma) | Phillips (1957) | as published |
| 1.2 | Stretched sigma levels | none stated | origin not recorded (power-law stretch); s = 1.4, N and p_t chosen by measurement here |
| 1.3 | Exner function, T from theta | textbook (none cited) | as published (textbook, Poisson's equation) |
| 1.4 | C-grid differences and two-point averages | Arakawa and Lamb (1977) | as published (C-grid staggering); replicate edges are this code's limited-area choice |
| 1.5 | Four-point v-to-u and u-to-v averages | Arakawa and Lamb (1977) | as published |
| 1.6 | Grid spacing from the domain | none stated | project's own (derived here); the metres-per-degree constants: origin not recorded |
| 1.7 | Grid shape for a target spacing | none stated | project's own (derived here) |
| 1.8 | Beta-plane Coriolis parameter | textbook (none cited) | as published (textbook); reason for a beta plane not recorded |
| 1.9 | Maximum terrain slope | NOAA National Geophysical Data Center (2009) | project's own (derived here) on ETOPO1 terrain (data source) |
| 1.10 | 1-2-1 terrain smoothing | NOAA National Geophysical Data Center (2009) | origin not recorded (the filter); applied to ETOPO1 until the slope is at most 0.0086 (P-56) |
| 1.11 | Forced ascent w = U S | none stated | project's own (derived here) |
| 1.12 | Lambert conformal constants n, F, rho | Snyder (1987) | as published (spherical form) |
| 1.13 | Lambert conformal x, y and scale factor | Snyder (1987) | as published (spherical form) |
| 2.1 | Surface-pressure (pi) tendency | Phillips (1957) | as published (sigma-coordinate primitive equations; the equations are textbook) |
| 2.2 | Horizontal momentum equations | Phillips (1957) | as published (textbook primitive equations in sigma) |
| 2.3 | Potential temperature equation | Phillips (1957) | as published (advective form; a flux form was reverted, P-42) |
| 2.4 | Discrete continuity, pi tendency | textbook (none cited) | as published (textbook); the lid flux F_top is the project's own addition for the radiating lid |
| 2.5 | Diagnostic sigma-dot | textbook (none cited) | adapted: textbook continuity with the constant written as F_top so sigma-dot is F_top at the lid and 0 at the ground (derived here) |
| 2.6 | Hydrostatic geopotential | textbook (none cited) | as published (hypsometric equation); the Simmons and Burridge (1981) form in the same file is not used (P-43) |
| 2.7 | Pressure-gradient force, plain form | textbook (none cited) | adapted: textbook sigma transformation with ln p differenced directly on the C grid instead of expanded (residual 2.1e-03 to 8.2e-15, P-14) |
| 2.8 | Pressure-gradient force, reference-state form | none stated | origin not recorded (called a "standard remedy") |
| 2.9 | Vertical advection, flux form with correction | none stated | origin not recorded |
| 2.10 | Horizontal advection, centred2 | none stated | origin not recorded (the Source line of section 2.6 names a paper for upwind3 only) |
| 2.11 | Horizontal advection, upwind3 | Wicker and Skamarock (2002) | as published (written as the fourth-order centred difference plus a \|u\|-scaled fourth-derivative term) |
| 2.12 | Biharmonic hyperdiffusion operator | none stated | origin not recorded |
| 2.13 | Hyperdiffusion coefficient | none stated | project's own (derived here) (discrete 2-dx eigenvalue, P-17); the 3 h damping time: origin not recorded |
| 2.14 | Divergence damping | Skamarock and Klemp (1992) | adapted: explicit tendency with nu = C dx dy / dt0, C = 0.0064 set by P-60 tests V and W |
| 2.15 | Sponge rate profile | none stated | origin not recorded (rate "chosen, not derived") |
| 2.16 | Sponge relaxation of u, v | none stated | origin not recorded; relaxing to a frozen reference rather than the mean is the project's choice (P-16) |
| 2.17 | Three-stage Runge-Kutta step | Wicker and Skamarock (2002) | as published |
| 2.18 | CFL time step | textbook (none cited) | adapted: textbook CFL with the Lamb-wave speed sqrt(R T) as the signal speed (docstring); safety 0.6 reason not recorded |
| 2.19 | Radiating upper boundary (option, off) | Klemp and Durran (1983); Bougeault (1983) | adapted: condition turned into a lid mass flux for Eq. 2.4; running low-pass reference, edge taper and spectral roll-off added here (P-50) |
| 2.20 | SPPT-style tendency perturbation (option, off) | Buizza et al. (1999) | adapted ("SPPT-style"): the Gaussian filter, AR(1) time correlation and the 0.9 clip are as coded here |
| 3.1 | NH prognostic variables | Laprise (1992); Skamarock and Klemp (2008) | as published |
| 3.2 | Inverse density and equation of state | Laprise (1992); Skamarock and Klemp (2008) | as published |
| 3.3 | Column-mass tendency | none stated | origin not recorded (no Source line; the coordinate is Laprise 1992, section 3.1) |
| 3.4 | Omega from continuity | none stated | origin not recorded (no Source line; the coordinate is Laprise 1992, section 3.1) |
| 3.5 | Perturbation pressure-gradient force | Laprise (1992); Skamarock and Klemp (2008) | as published; the reference balance is removed against this code's sounding (Eq. 3.12) |
| 3.6 | U and V momentum equations | Laprise (1992); Skamarock and Klemp (2008); Wicker and Skamarock (2002) | as published; horizontal momentum advection is upwind3 (Eq. 2.11) |
| 3.7 | W equation | Laprise (1992); Skamarock and Klemp (2008) | as published |
| 3.8 | Geopotential equation | Laprise (1992); Skamarock and Klemp (2008) | as published |
| 3.9 | Theta, flux form | Wicker and Skamarock (2002) | adapted: third-order flux in x, y and the column interior, adopted after the S1b density-current failure |
| 3.10 | Third-order upwind-biased flux | Wicker and Skamarock (2002) | as published (flux form of upwind3) |
| 3.11 | Third-order vertical face values | Wicker and Skamarock (2002) | adapted: third order on faces 2 to N-2 only; two-level average at h = 1 and N-1, end values at the boundaries |
| 3.12 | Reference sounding, isothermal 250 K | none stated | origin not recorded (the reference "only has to be hydrostatic") |
| 3.13 | Reference mu, phi, alpha, theta | none stated | project's own (derived here) (theta-bar set so p' = 0 at rest with the discrete equation of state) |
| 3.14 | Linearised acoustic pressure | Klemp et al. (2007); Klemp and Wilhelmson (1978) | as published |
| 3.15 | Acoustic divergence damping | Skamarock and Klemp (1992) | as published (pressure-based form); coefficient 0.1 reason not recorded |
| 3.16 | Acoustic U update | Klemp et al. (2007); Klemp and Wilhelmson (1978) | as published |
| 3.17 | Acoustic Theta update | Klemp et al. (2007); Klemp and Wilhelmson (1978) | as published |
| 3.18 | Explicit parts of phi and W | Klemp et al. (2007); Klemp and Wilhelmson (1978) | as published; off-centring 0.1 reason not recorded |
| 3.19 | Vertically implicit W tridiagonal | Klemp et al. (2007); Klemp and Wilhelmson (1978) | as published (Thomas algorithm, textbook) |
| 3.20 | Implicit coefficients and phi update | Klemp et al. (2007); Klemp and Wilhelmson (1978) | as published |
| 3.21 | Rayleigh damping of W | none stated | origin not recorded (factor 2 on the sponge rate not explained) |
| 3.22 | NH large time step | none stated | project's own (derived here) (acoustic limit added in S5); 60 s cap and 0.7 reason not recorded |
| 3.23 | 2-D Exner-form benchmark equations | Klemp and Wilhelmson (1978); Skamarock and Klemp (1992) | as published (benchmark core only) |
| 4.1 | Lowest-level height and layer depth | textbook (none cited) | as published (hypsometric equation) |
| 4.2 | Neutral drag coefficient | textbook (none cited) | as published (log law); roughness values "conventional", source not named |
| 4.3 | Bulk Richardson number | none stated | origin not recorded |
| 4.4 | Cut-off drag stability function | Louis (1979) | adapted: "Louis-type" shape; the constants 10, 4 and 0.2 have no recorded origin |
| 4.5 | Louis momentum stability function | Louis (1979) | as published (b = d = 5, c_m = 7.4) |
| 4.6 | Surface drag tendency | textbook (none cited) | as published (bulk aerodynamic formula) |
| 4.7 | Gradient Richardson number, N2, S2 | none stated | origin not recorded |
| 4.8 | Eddy diffusivity K(Ri) | Louis (1979); Miles (1961); Howard (1961) | adapted: the Louis (1979) shape "with the constants exposed"; constants from Miles (1961) and Howard (1961) for Ri_c = 0.25; K_max 200 measured (P-40); l = 150 m reason not recorded |
| 4.9 | Vertical diffusion flux divergence | none stated | origin not recorded |
| 4.10 | Pool-adjacent-violators adjustment | Manabe et al. (1965); Ayer et al. (1955) | adapted: Manabe et al. dry adjustment solved exactly by the Ayer et al. algorithm; u, v and q mixed over the same blocks (P-67, P-74) |
| 4.11 | Force-restore ground temperature | Deardorff (1978) | as published; soil constants given with reasons, source not named |
| 4.12 | Net radiation and clear-sky sunshine | Brutsaert (1975); Kopp and Lean (2011); Cooper (1969) | project's own (derived here); constants from Brutsaert (1975) (epsilon_a) and Kopp and Lean (2011) (S0); sin e from Eq. 4.17 |
| 4.13 | Sensible and latent heat fluxes | Louis (1979) | project's own (derived here) (bulk heat flux with the Louis 1979 transfer function); LE = EF max(Rn, 0) is the dry model's stand-in |
| 4.14 | Louis heat stability function with long tail | Louis (1979) | adapted: Louis form on the unstable side; stable side uses the momentum function's tail, "a choice made from that measurement ... not taken from a reference" |
| 4.15 | Heating of the lowest layer by H | none stated | origin not recorded |
| 4.16 | Brutsaert clear-sky emissivity | Brutsaert (1975) | as published (not called; the constant 0.79 is its value at e = 12 hPa, T = 285 K) |
| 4.17 | Solar declination, hour angle, elevation | Cooper (1969) | as published for the declination (Cooper 1969); the elevation formula is not cited; equation of time ignored |
| 4.18 | Prescribed surface heat flux (option, off) | Kopp and Lean (2011) | project's own (derived here); S0 from Kopp and Lean (2011); f, tau and H_night: origin not recorded |
| 5.1 | Density potential temperature | Emanuel (1994) | as published (Emanuel 1994, sec. 4.3) |
| 5.2 | Bolton saturation vapour pressure and q_s | Bolton (1980) | as published |
| 5.3 | Vapour mixing ratio from RH | Bolton (1980) | origin not recorded (the conversion); e_s from Bolton (1980) |
| 5.4 | Dewpoint from q_v and p | Bolton (1980) | adapted: the exact algebraic inverse of Bolton's formula (Eq. 5.2) |
| 5.5 | Substep-mean mass fluxes | Skamarock and Klemp (2008) | project's own (derived here) (chosen on the uniform-q requirement); the same as WRF ARW per Skamarock and Klemp (2008), which the S7 record does not cite |
| 5.6 | Moisture stage update | Wicker and Skamarock (2002); Skamarock and Klemp (2008) | adapted: Theta's third-order flux (Eq. 3.10) applied to mu q with the substep-mean fluxes of Eq. 5.5 |
| 5.7 | Column hole filling | none stated | origin not recorded |
| 5.8 | Rain terminal velocity and fall substeps | Klemp and Wilhelmson (1978); WRF `module_mp_kessler.F`; Kessler (1969) | as published (transcribed from WRF module_mp_kessler.F; terminal velocity of Klemp and Wilhelmson 1978) |
| 5.9 | Rain sedimentation and surface rain | Kessler (1969); Klemp and Wilhelmson (1978); WRF `module_mp_kessler.F` | as published (WRF transcription; index 0 is the lid here, the ground in WRF) |
| 5.10 | Autoconversion and accretion | Kessler (1969); Klemp and Wilhelmson (1978); WRF `module_mp_kessler.F` | as published (constants from the WRF code) |
| 5.11 | Saturation adjustment, one Newton step | Kessler (1969); Klemp and Wilhelmson (1978); WRF `module_mp_kessler.F`; Bolton (1980) | as published (WRF transcription; e_s from Bolton 1980) |
| 5.12 | Rain evaporation | Klemp and Wilhelmson (1978); WRF `module_mp_kessler.F` | as published (Klemp and Wilhelmson 1978, eq. 2.14, via WRF) |
| 5.13 | Kessler state update | Kessler (1969); Klemp and Wilhelmson (1978); WRF `module_mp_kessler.F` | as published (WRF transcription) |
| 5.14 | Surface vapour flux | Louis (1979); Bolton (1980) | project's own (derived here) (E = LE/Lv over land keeps the S6 energy balance; bulk formula over water); C_h from Eq. 4.13, q_s from Eq. 5.2 |
| 6.1 | Davies relaxation weights | Davies (1976) | as published (cosine taper, width 15, alpha_max 0.1) |
| 6.2 | Davies relaxation update | Davies (1976) | as published (applied to the state after each step) |
| 6.3 | Relaxation weight per unit time | none stated | project's own (derived here) (P-70) |
| 6.4 | Hydrostatic phi relaxation and W damping in the zone | none stated | project's own (derived here) |
| 6.5 | Linear time interpolation of frames | none stated | project's own (derived here) |
| 6.6 | Nest bilinear regridding with staggering | none stated | project's own (derived here) |
| 6.7 | Nest surface pressure moved to new terrain | textbook (none cited) | as published (hypsometric equation, layer-mean T from the standard lapse rate) |
| 6.8 | Nest vertical remapping in ln p | none stated | project's own (derived here) (ln theta linear in ln p) |
| 7.1 | Surface pressure by interpolation in height | none stated | project's own (derived here) |
| 7.2 | Hydrostatic extension below the lowest level | none stated | project's own (derived here) |
| 7.3 | Interpolation linear in ln p | none stated | project's own (derived here) |
| 7.4 | Raised-cosine spectral low-pass | none stated | project's own (derived here); cut-offs 4 and 8 grid lengths reason not recorded |
| 7.5 | Helmholtz divergence removal | textbook (none cited) | adapted: textbook Helmholtz split solved with the eigenvalues of the discrete Laplacian (P-15) |
| 7.6 | Divergence polish passes | none stated | project's own (derived here); coefficient and pass count reason not recorded |
| 7.7 | Successive-correction increments | Barnes (1964); Koch et al. (1983) | adapted: Barnes analysis of innovations with lambda in the denominator and error weights; multi-pass gamma from Koch et al. (1983); length scales, lambda and gamma values reason not recorded |
| 7.8 | Tangent-plane distance | none stated | project's own (derived here) |
| 7.9 | Terrain- and coast-aware weights | none stated | project's own (derived here) (S5b; constants chosen before the run) |
| 7.10 | Standard-atmosphere first guess | textbook (none cited) | as published (US standard atmosphere, named without a citation) |
| 7.11 | Virtual temperature | textbook (none cited); Bolton (1980) | as published (textbook); e_s from Bolton (1980) |
| 7.12 | Hydrostatic heights from sea-level pressure | textbook (none cited) | as published (hypsometric equation, two fixed-point passes at the bottom) |
| 7.13 | Surface blend and lapse correction | none stated | project's own (derived here); 1000 m depth reason not recorded |
| 7.14 | RH from dewpoint (Magnus) | Alduchov and Eskridge (1996) | constants from Alduchov and Eskridge (1996) |
| 8.1 | Ground temperature and MSLP reduction | none stated | origin not recorded (textbook standard-lapse reduction) |
| 8.2 | Below-ground T and Z on pressure levels | none stated | origin not recorded (textbook) |
| 8.3 | Absolute vorticity and kinematic omega | none stated | origin not recorded (textbook) |
| 8.4 | Station temperature, lapse correction | none stated | project's own (derived here) |
| 8.5 | Neutral 10 m wind factor | none stated | project's own (derived here) (P-68; the same log law as Eq. 4.2) |
| 8.6 | Bulk Richardson number from zeta | Paulson (1970); Beljaars and Holtslag (1991) | as published (Monin-Obukhov similarity, psi from Eqs. 8.8-8.9) |
| 8.7 | 10 m wind and 2 m theta profiles | Paulson (1970); Beljaars and Holtslag (1991) | as published (Monin-Obukhov similarity); z0h = z0/10 reason not recorded |
| 8.8 | Unstable psi functions | Paulson (1970) | as published (Businger-Dyer forms integrated by Paulson 1970) |
| 8.9 | Stable psi functions | Beljaars and Holtslag (1991) | as published (a = 1, b = 2/3, c = 5, d = 0.35) |
| 8.10 | Station dewpoint | Bolton (1980) | project's own (derived here); T_d from Eq. 5.4 |
| 8.11 | Bias, RMSE, MAE, correlation, skill | none stated | origin not recorded (standard scores) |

---

## 12. References

In the style of `docs/REFERENCES.md` (AMS). Entries marked * are named in the
code but are not in `docs/REFERENCES.md`; their bibliographic details beyond
what the code gives were not checked against the repository.

Alduchov, O. A., and R. E. Eskridge, 1996: Improved Magnus form approximation of saturation vapor pressure. *J. Appl. Meteor.*, **35**, 601–609, https://doi.org/10.1175/1520-0450(1996)035<0601:IMFAOS>2.0.CO;2. Used in: Eq. 7.14.

Arakawa, A., and V. R. Lamb, 1977: Computational design of the basic dynamical processes of the UCLA general circulation model. *Methods in Computational Physics: Advances in Research and Applications*, Vol. 17, Academic Press, 173–265, https://doi.org/10.1016/B978-0-12-460817-7.50009-4. Used in: Eqs. 1.4, 1.5.

Ayer, M., H. D. Brunk, G. M. Ewing, W. T. Reid, and E. Silverman, 1955: An empirical distribution function for sampling with incomplete information. *Ann. Math. Stat.*, **26**, 641–647, https://doi.org/10.1214/aoms/1177728423. Used in: Eq. 4.10.

Barnes, S. L., 1964: A technique for maximizing details in numerical weather map analysis. *J. Appl. Meteor.*, **3**, 396–409, https://doi.org/10.1175/1520-0450(1964)003<0396:ATFMDI>2.0.CO;2. Used in: Eq. 7.7.

Beljaars, A. C. M., and A. A. M. Holtslag, 1991: Flux parameterization over land surfaces for atmospheric models. *J. Appl. Meteor.*, **30**, 327–341, https://doi.org/10.1175/1520-0450(1991)030<0327:FPOLSF>2.0.CO;2. Used in: Eqs. 8.6, 8.7, 8.9.

Bolton, D., 1980: The computation of equivalent potential temperature. *Mon. Wea. Rev.*, **108**, 1046–1053, https://doi.org/10.1175/1520-0493(1980)108<1046:TCOEPT>2.0.CO;2. Used in: Eqs. 5.2, 5.3, 5.4, 5.11, 5.14, 7.11, 8.10.

Bougeault, P., 1983: A non-reflective upper boundary condition for limited-height hydrostatic models. *Mon. Wea. Rev.*, **111**, 420–429, https://doi.org/10.1175/1520-0493(1983)111<0420:ANRUBC>2.0.CO;2. Used in: Eq. 2.19.

Brutsaert, W., 1975: On a derivable formula for long-wave radiation from clear skies. *Water Resour. Res.*, **11**, 742–744, https://doi.org/10.1029/WR011i005p00742. Used in: Eqs. 4.12, 4.16.

Bryan, G. H., and J. M. Fritsch, 2002: A benchmark simulation for moist nonhydrostatic numerical models. *Mon. Wea. Rev.*, **130**, 2917–2928, https://doi.org/10.1175/1520-0493(2002)130<2917:ABSFMN>2.0.CO;2. Used in: background only (cited in prose, section 5.7).

Buizza, R., M. Miller, and T. N. Palmer, 1999: Stochastic representation of model uncertainties in the ECMWF ensemble prediction system. *Quart. J. Roy. Meteor. Soc.*, **125**, 2887–2908, https://doi.org/10.1002/qj.49712556006. Used in: Eq. 2.20.

Cooper, P. I., 1969: The absorption of radiation in solar stills. *Sol. Energy*, **12**, 333–346, https://doi.org/10.1016/0038-092X(69)90047-4. Used in: Eqs. 4.12, 4.17.

Davies, H. C., 1976: A lateral boundary formulation for multi-level prediction models. *Quart. J. Roy. Meteor. Soc.*, **102**, 405–418, https://doi.org/10.1002/qj.49710243210. Used in: Eqs. 6.1, 6.2.

Deardorff, J. W., 1978: Efficient prediction of ground surface temperature and moisture, with inclusion of a layer of vegetation. *J. Geophys. Res.*, **83**, 1889–1903, https://doi.org/10.1029/JC083iC04p01889. Used in: Eq. 4.11.

\*Emanuel, K. A., 1994: *Atmospheric Convection*. Oxford University Press. Cited in `moisture.py` as "Emanuel 1994, sec. 4.3"; title, publisher and author initials not confirmed from the repository. Used in: Eq. 5.1.

Howard, L. N., 1961: Note on a paper of John W. Miles. *J. Fluid Mech.*, **10**, 509–512, https://doi.org/10.1017/S0022112061000317. Used in: Eq. 4.8.

Kessler, E., 1969: *On the Distribution and Continuity of Water Substance in Atmospheric Circulations*. Meteor. Monogr., No. 32, Amer. Meteor. Soc., 84 pp., https://doi.org/10.1007/978-1-935704-36-2_1. Used in: Eqs. 5.8, 5.9, 5.10, 5.11, 5.13.

Klemp, J. B., and D. R. Durran, 1983: An upper boundary condition permitting internal gravity wave radiation in numerical mesoscale models. *Mon. Wea. Rev.*, **111**, 430–444, https://doi.org/10.1175/1520-0493(1983)111<0430:AUBCPI>2.0.CO;2. Used in: Eq. 2.19.

Klemp, J. B., and R. B. Wilhelmson, 1978: The simulation of three-dimensional convective storm dynamics. *J. Atmos. Sci.*, **35**, 1070–1096, https://doi.org/10.1175/1520-0469(1978)035<1070:TSOTDC>2.0.CO;2. Used in: Eqs. 3.14, 3.16, 3.17, 3.18, 3.19, 3.20, 3.23, 5.8, 5.9, 5.10, 5.11, 5.12, 5.13.

Klemp, J. B., W. C. Skamarock, and J. Dudhia, 2007: Conservative split-explicit time integration methods for the compressible nonhydrostatic equations. *Mon. Wea. Rev.*, **135**, 2897–2913, https://doi.org/10.1175/MWR3440.1. Used in: Eqs. 3.14, 3.16, 3.17, 3.18, 3.19, 3.20.

Koch, S. E., M. desJardins, and P. J. Kocin, 1983: An interactive Barnes objective map analysis scheme for use with satellite and conventional data. *J. Climate Appl. Meteor.*, **22**, 1487–1503, https://doi.org/10.1175/1520-0450(1983)022<1487:AIBOMA>2.0.CO;2. Used in: Eq. 7.7.

Kopp, G., and J. L. Lean, 2011: A new, lower value of total solar irradiance: Evidence and climate significance. *Geophys. Res. Lett.*, **38**, L01706, https://doi.org/10.1029/2010GL045777. Used in: Eqs. 4.12, 4.18.

Laprise, R., 1992: The Euler equations of motion with hydrostatic pressure as an independent variable. *Mon. Wea. Rev.*, **120**, 197–207, https://doi.org/10.1175/1520-0493(1992)120<0197:TEEOMW>2.0.CO;2. Used in: Eqs. 3.1, 3.2, 3.5, 3.6, 3.7, 3.8.

Louis, J.-F., 1979: A parametric model of vertical eddy fluxes in the atmosphere. *Bound.-Layer Meteor.*, **17**, 187–202, https://doi.org/10.1007/BF00117978. Used in: Eqs. 4.4, 4.5, 4.8, 4.13, 4.14, 5.14.

Manabe, S., J. Smagorinsky, and R. F. Strickler, 1965: Simulated climatology of a general circulation model with a hydrologic cycle. *Mon. Wea. Rev.*, **93**, 769–798, https://doi.org/10.1175/1520-0493(1965)093<0769:SCOAGC>2.3.CO;2. Used in: Eq. 4.10.

Miles, J. W., 1961: On the stability of heterogeneous shear flows. *J. Fluid Mech.*, **10**, 496–508, https://doi.org/10.1017/S0022112061000305. Used in: Eq. 4.8.

NOAA National Geophysical Data Center, 2009: ETOPO1 1 arc-minute global relief model. NOAA National Centers for Environmental Information, accessed 26 September 2026 through the NOAA CoastWatch ERDDAP server (dataset `etopo180`), https://doi.org/10.7289/V5C8276M. Used in: Eqs. 1.9, 1.10.

Paulson, C. A., 1970: The mathematical representation of wind speed and temperature profiles in the unstable atmospheric surface layer. *J. Appl. Meteor.*, **9**, 857–861, https://doi.org/10.1175/1520-0450(1970)009<0857:TMROWS>2.0.CO;2. Used in: Eqs. 8.6, 8.7, 8.8.

Phillips, N. A., 1957: A coordinate system having some special advantages for numerical forecasting. *J. Meteor.*, **14**, 184–185, https://doi.org/10.1175/1520-0469(1957)014<0184:ACSHSS>2.0.CO;2. Used in: Eqs. 1.1, 2.1, 2.2, 2.3.

Roberts, N. M., and H. W. Lean, 2008: Scale-selective verification of rainfall accumulations from high-resolution forecasts of convective events. *Mon. Wea. Rev.*, **136**, 78–97, https://doi.org/10.1175/2007MWR2123.1. Used in: background only (cited in prose, section 8.8).

Simmons, A. J., and D. M. Burridge, 1981: An energy and angular-momentum conserving vertical finite-difference scheme and hybrid vertical coordinates. *Mon. Wea. Rev.*, **109**, 758–766, https://doi.org/10.1175/1520-0493(1981)109<0758:AEAAMC>2.0.CO;2. Used in: background only (cited in prose, section 2.3).

Skamarock, W. C., and J. B. Klemp, 1992: The stability of time-split numerical methods for the hydrostatic and the nonhydrostatic elastic equations. *Mon. Wea. Rev.*, **120**, 2109–2127, https://doi.org/10.1175/1520-0493(1992)120<2109:TSOTSN>2.0.CO;2. Used in: Eqs. 2.14, 3.15, 3.23.

Skamarock, W. C., and J. B. Klemp, 2008: A time-split nonhydrostatic atmospheric model for weather research and forecasting applications. *J. Comput. Phys.*, **227**, 3465–3485, https://doi.org/10.1016/j.jcp.2007.01.037. Used in: Eqs. 3.1, 3.2, 3.5, 3.6, 3.7, 3.8, 5.5, 5.6.

Smith, R. B., 1979: The influence of mountains on the atmosphere. *Adv. Geophys.*, **21**, 87–230, https://doi.org/10.1016/S0065-2687(08)60262-9. Used in: background only (cited in prose, introduction to section 3 and section 3.11).

Snyder, J. P., 1987: *Map Projections—A Working Manual*. U.S. Geological Survey Professional Paper 1395, https://doi.org/10.3133/pp1395. Used in: Eqs. 1.12, 1.13.

Straka, J. M., R. B. Wilhelmson, L. J. Wicker, J. R. Anderson, and K. K. Droegemeier, 1993: Numerical solutions of a non-linear density current: A benchmark solution and comparisons. *Int. J. Numer. Methods Fluids*, **17**, 1–22, https://doi.org/10.1002/fld.1650170103. Used in: background only (cited in prose, introduction to section 3 and section 3.11).

Wicker, L. J., and W. C. Skamarock, 2002: Time-splitting methods for elastic models using forward time schemes. *Mon. Wea. Rev.*, **130**, 2088–2097, https://doi.org/10.1175/1520-0493(2002)130<2088:TSMFEM>2.0.CO;2. Used in: Eqs. 2.11, 2.17, 3.6, 3.9, 3.10, 3.11, 5.6.

\*WRF (Weather Research and Forecasting model) source code, `module_mp_kessler.F`. Named in `moisture.py` as the source of the transcription; the WRF version is not recorded. Used in: Eqs. 5.8, 5.9, 5.10, 5.11, 5.12, 5.13.

Named in the code without a separate reference: the Businger-Dyer profile
functions (through Paulson 1970) and the US standard atmosphere.

"""
Dissipation and stochastic variance -- what a model needs to run on real
atmospheres rather than analytic test cases.

Two separate problems.

DISSIPATION. Nonlinear flow cascades energy toward small scales. A grid can
only represent scales down to 2*dx, so without a sink energy accumulates
there as grid-scale noise until the field is dominated by it. Real models
apply scale-selective damping: strong at the grid scale, negligible on
resolved features. Fourth-order hyperdiffusion is the standard choice --
its damping goes as k^4, so halving the wavelength multiplies the damping
by 16, leaving synoptic features essentially untouched.

STOCHASTIC VARIANCE. The equations describe resolved motion. Everything
unresolved -- turbulence, convection below the grid scale, terrain drag
detail -- is absent, and its effect on the resolved flow is not zero and not
deterministic. Operational centres add stochastic perturbations (ECMWF's
SPPT, SKEB) because a single deterministic integration is systematically
overconfident: it produces one trajectory when the atmosphere admits a
distribution of them.

The perturbations here are SPPT-style: tendencies are multiplied by (1 + r),
where r is a random field that is SMOOTH in space and CORRELATED in time.
Both properties matter. White noise would be scrubbed out by diffusion and
would inject grid-scale energy; a field that changes discontinuously each
step would act like noise rather than like a persistent unresolved process.
"""

import numpy as np


# ---------------------------------------------------------------------------
# Dissipation
# ---------------------------------------------------------------------------

def hyperdiffusion(a, grid, coeff):
    """
    -K * del^4 a, applied as the negative biharmonic so it damps.

    coeff has units m^4/s. See recommended_hyper_coeff() for scaling.
    """
    lap = ((grid.shift(a, 1, 1) - 2 * a + grid.shift(a, -1, 1)) / grid.dx**2 +
           (grid.shift(a, 1, 0) - 2 * a + grid.shift(a, -1, 0)) / grid.dy**2)
    lap2 = ((grid.shift(lap, 1, 1) - 2 * lap + grid.shift(lap, -1, 1)) / grid.dx**2 +
            (grid.shift(lap, 1, 0) - 2 * lap + grid.shift(lap, -1, 0)) / grid.dy**2)
    return -coeff * lap2


def discrete_biharmonic_eigenvalue(grid):
    """
    Response of the DISCRETE biharmonic operator to the 2dx checkerboard --
    the shortest wave the grid can hold, and the one noise lives on.

    For a checkerboard the discrete Laplacian gives -(4/dx^2 + 4/dy^2), so the
    biharmonic gives that squared. Note this is NOT (pi/dx)^4: the discrete
    Laplacian's response at 2dx is 4/dx^2, while the continuous operator would
    give (pi/dx)^2 = 9.87/dx^2. Using the continuous value makes the damping
    ~6x weaker than intended -- which is exactly the bug this replaced.
    """
    return (4.0 / grid.dx**2 + 4.0 / grid.dy**2) ** 2


def recommended_hyper_coeff(grid, damping_time=3 * 3600.0):
    """
    Coefficient that damps the 2dx checkerboard on `damping_time`.

    Three hours at the grid scale is a common operational choice: fast enough
    to control noise, slow enough to leave resolved weather alone. Waves at
    4dx are damped ~16x more slowly, at 8dx ~256x, and so on.
    """
    return 1.0 / (damping_time * discrete_biharmonic_eigenvalue(grid))


def hyper_stability_dt(grid, coeff, safety=0.5):
    """
    Explicit biharmonic diffusion has its own stability limit, which can be
    tighter than the advective CFL if the coefficient is large. Check it.
    """
    dx = min(grid.dx, grid.dy)
    return safety * dx**4 / (16.0 * coeff) if coeff > 0 else np.inf


# ---------------------------------------------------------------------------
# Stochastic variance
# ---------------------------------------------------------------------------

class StochasticPerturbation:
    """
    SPPT-style multiplicative tendency perturbation.

        tendency -> tendency * (1 + r)

    r is generated as spatially smoothed Gaussian noise, evolved in time as a
    first-order autoregressive process so it decorrelates over `tau` rather
    than flickering every step.

    Parameters
    ----------
    amplitude : std dev of r. 0.2-0.5 is the operational range; larger values
                make individual runs unstable rather than merely uncertain.
    tau       : temporal decorrelation time (s). ~6 h is typical.
    length_scale : spatial smoothing scale (m). Should be well above the grid
                scale, or the perturbations are just noise for the diffusion
                to remove -- and well below the domain size, or every
                resolvable mode is filtered away. Capped at a quarter of the
                domain for that reason.
    clip      : r is truncated to +-clip to prevent a tail event from
                reversing the sign of a tendency, which is unphysical.
    """

    def __init__(self, grid, amplitude=0.3, tau=6 * 3600.0,
                 length_scale=500e3, clip=0.9, seed=None, nz=1):
        self.grid = grid
        self.amplitude = float(amplitude)
        self.tau = float(tau)
        self.length_scale = float(length_scale)
        self.clip = float(clip)
        self.rng = np.random.default_rng(seed)
        self.nz = nz

        self.r = self._draw()

    def _smooth(self, field):
        """
        Spectral smoothing: keep only wavenumbers below the length scale.
        Implemented with an FFT so the filter is exact and isotropic.
        """
        gr = self.grid
        ny, nx = gr.ny, gr.nx
        kx = np.fft.fftfreq(nx, d=gr.dx) * 2 * np.pi
        ky = np.fft.fftfreq(ny, d=gr.dy) * 2 * np.pi
        KX, KY = np.meshgrid(kx, ky)
        k2 = KX**2 + KY**2

        # A length scale approaching the domain size filters out every
        # resolvable mode and leaves a constant field with zero variance.
        # Cap it so at least a few wavelengths fit.
        L = min(self.length_scale, 0.25 * min(gr.Lx, gr.Ly))

        filt = np.exp(-0.5 * k2 * L**2)
        filt[0, 0] = 0.0        # drop the mean: a perturbation has zero mean

        out = np.real(np.fft.ifft2(np.fft.fft2(field, axes=(-2, -1)) * filt,
                                   axes=(-2, -1)))
        std = out.std()
        if std < 1e-30:
            # Should not happen now, but never silently return a dead field.
            raise RuntimeError(
                f"stochastic smoothing produced zero variance "
                f"(length_scale={self.length_scale:.0f} m, "
                f"domain {gr.Lx:.0f}x{gr.Ly:.0f} m)")
        return out / std

    def _draw(self):
        shape = ((self.nz, self.grid.ny, self.grid.nx) if self.nz > 1
                 else (self.grid.ny, self.grid.nx))
        raw = self.rng.normal(0.0, 1.0, shape)
        return np.clip(self._smooth(raw) * self.amplitude,
                       -self.clip, self.clip)

    def advance(self, dt):
        """
        AR(1) update:  r <- a r + sqrt(1 - a^2) * new,   a = exp(-dt/tau)

        The sqrt keeps the variance stationary; without it the field would
        decay toward zero over a long run.
        """
        a = np.exp(-dt / self.tau)
        self.r = np.clip(a * self.r + np.sqrt(max(0.0, 1 - a * a)) * self._draw(),
                         -self.clip, self.clip)
        return self.r

    def apply(self, tendency):
        """Multiply a tendency by (1 + r), broadcasting over levels."""
        r = self.r
        if tendency.ndim == 3 and r.ndim == 2:
            r = r[None, :, :]
        return tendency * (1.0 + r)

    def __repr__(self):
        return (f"StochasticPerturbation(amp={self.amplitude}, "
                f"tau={self.tau/3600:.1f}h, L={self.length_scale/1000:.0f}km)")


def perturb_initial_state(field, rng, amplitude, grid, length_scale=300e3):
    """
    Add a smooth random perturbation to an initial field.

    This is how ensemble members are generated: the atmosphere's state is
    never known exactly, and small initial differences grow. Perturbing the
    initial condition samples that uncertainty, and is a different mechanism
    from perturbing the tendencies (model error).
    """
    p = StochasticPerturbation(grid, amplitude=1.0, length_scale=length_scale,
                               seed=int(rng.integers(2**31)),
                               nz=field.shape[0] if field.ndim == 3 else 1)
    r = p.r
    if field.ndim == 3 and r.ndim == 2:
        r = r[None, :, :]
    return field + amplitude * r


# ---------------------------------------------------------------------------
# Initialisation: removing spurious divergence
# ---------------------------------------------------------------------------

def divergence_damping(u, v, grid, coeff):
    """
    Tendency that damps the DIVERGENT part of the flow, leaving the rotational
    part untouched:

        du/dt += nu * d(div)/dx ,   dv/dt += nu * d(div)/dy

    Standard in operational models. It targets exactly the component that
    drives spurious vertical motion, without smearing the vorticity that
    carries the weather.
    """
    div = grid.dx_forward(u) + grid.dy_forward(v)
    return coeff * grid.dx_backward(div), coeff * grid.dy_backward(div)


def remove_divergence_spectral(u, v, grid):
    """
    Remove the divergent part of a wind field EXACTLY, in one step.

    Helmholtz: any flow splits into rotational + divergent parts. The
    divergent part is the gradient of a velocity potential chi satisfying

        laplacian(chi) = div(u, v)

    Subtracting grad(chi) leaves a non-divergent field. Solving in Fourier
    space with the eigenvalues of OUR DISCRETE Laplacian -- not the continuous
    -k^2 -- means the discrete divergence cancels to machine precision rather
    than approximately.

    The FFT assumes periodicity, which this domain does not have, so expect
    some error in the outermost cells. The relaxation zone overwrites those
    anyway, and the interior gain is worth far more.
    """
    div = grid.dx_forward(u) + grid.dy_forward(v)

    ny, nx = div.shape[-2:]
    kx = 2 * np.pi * np.fft.fftfreq(nx)
    ky = 2 * np.pi * np.fft.fftfreq(ny)
    KX, KY = np.meshgrid(kx, ky)

    # Eigenvalue of the 5-point Laplacian built from our forward/backward pair.
    lam = (-4 * np.sin(KX / 2) ** 2 / grid.dx ** 2
           - 4 * np.sin(KY / 2) ** 2 / grid.dy ** 2)
    lam[0, 0] = 1.0                     # mean mode: no correction

    chi_hat = np.fft.fft2(div, axes=(-2, -1)) / lam
    chi_hat[..., 0, 0] = 0.0
    chi = np.real(np.fft.ifft2(chi_hat, axes=(-2, -1)))

    return u - grid.dx_backward(chi), v - grid.dy_backward(chi)


def balance_initial_state(u, v, grid, target_div=2e-5, max_iter=400,
                          verbose=True, method="hybrid", edge=6):
    """
    Iteratively remove grid-scale divergence from an initial wind field.

    WHY THIS IS NEEDED

    Analysis winds are not in balance with OUR discretisation. Interpolated,
    coarsened, and differenced with our operators, they carry divergence of
    order 1e-3 s^-1 where the real atmosphere has ~1e-5. Integrated over an
    800 hPa column that is tens of Pa/s of vertical motion -- a hundred times
    reality -- which drives a violent gravity-wave adjustment and destroys the
    forecast in the first hour.

    Real models solve this with digital-filter or variational initialisation.
    This is the cheap version: gradient descent on divergence squared, which
    is what divergence damping does when iterated to convergence. The
    rotational flow -- the part carrying the weather -- is preserved.

    Returns (u, v, info).
    """
    u, v = u.copy(), v.copy()

    def maxdiv(a, b):
        return float(np.abs(grid.dx_forward(a) + grid.dy_forward(b)).max())

    d0 = maxdiv(u, v)
    speed0 = float(np.abs(u).max())

    if method in ("spectral", "hybrid"):
        # Exact in one shot. Applied per level for a 3D field.
        if u.ndim == 3:
            for k in range(u.shape[0]):
                u[k], v[k] = remove_divergence_spectral(u[k], v[k], grid)
        else:
            u, v = remove_divergence_spectral(u, v, grid)
        it = 1
        if method == "hybrid":
            # The FFT assumes periodicity, so the outermost cells retain
            # error. A short iterative polish cleans them; the interior is
            # already exact and is barely touched.
            coeff = 0.2 * min(grid.dx, grid.dy) ** 2
            for _ in range(60):
                du, dv = divergence_damping(u, v, grid, coeff)
                u += du
                v += dv
            it = 61
    else:
        coeff = 0.2 * min(grid.dx, grid.dy) ** 2
        it = 0
        for it in range(1, max_iter + 1):
            du, dv = divergence_damping(u, v, grid, coeff)
            u += du
            v += dv
            if it % 10 == 0 and maxdiv(u, v) < target_div:
                break

    d = maxdiv(u, v)
    e = edge
    sl = (Ellipsis, slice(e, -e), slice(e, -e))
    d_interior = float(np.abs((grid.dx_forward(u) + grid.dy_forward(v))[sl]).max())
    info = {
        "div_interior": d_interior,
        "iterations": it,
        "div_before": d0,
        "div_after": d,
        "omega_before_Pa_s": d0 * 80000,
        "omega_after_Pa_s": d * 80000,
        "speed_before": speed0,
        "speed_after": float(np.abs(u).max()),
        "converged": d < target_div,
    }
    if verbose:
        print(f"  initialisation : max|div| {d0:.2e} -> {d:.2e} 1/s "
              f"({method}, {it} pass{'es' if it > 1 else ''})"
              + ("" if info["converged"] else "  (did not reach target)"))
        print(f"                   interior max|div| {d_interior:.2e} 1/s "
              f"(edges are relaxed to the driver anyway)")
        print(f"                   implied omega {info['omega_before_Pa_s']:.1f}"
              f" -> {d_interior * 80000:.2f} Pa/s interior; "
              f"max|u| {speed0:.1f} -> {info['speed_after']:.1f} m/s")
    return u, v, info

"""
Radiative (non-reflecting) upper boundary.

WHY A DAMPING LAYER CANNOT DO THIS JOB (P-49)

A Rayleigh sponge removes wave AMPLITUDE. Measured over 48 hours, every
setting of it -- any depth from two levels to twelve, any rate from fifteen
minutes to six hours, frozen or running reference, lid at 200, 100 or 50 hPa,
20, 26 or 30 levels -- turns a growing baroclinic wave into a decaying one:

    absorbing well (reflection 21-36 m/s)  ->  development 0.85-0.96 x/day
    development alive (1.89-2.45)          ->  reflection 50-60, not absorbing

The reason is that amplitude is not what distinguishes the two. A mountain
wave and a baroclinic wave sit in the same levels at similar amplitudes, and a
sponge sees only "deviation from a reference". It cannot tell them apart
because it is not measuring anything that differs.

WHAT DOES DIFFER: VERTICAL ENERGY FLUX

A gravity wave PROPAGATES vertically -- it carries energy upward out of the
domain. Balanced flow does not; its vertical group velocity is essentially
zero. A boundary that passes vertical wave flux and nothing else is therefore
selective by construction, and it damps nothing at all: it just stops the
energy coming back.

THE CONDITION

For hydrostatic internal gravity waves, linearised about a resting reference,
the upward-radiating solution satisfies, for each horizontal wavenumber k:

    w_hat(k) = (|k| / N) phi_hat'(k)

(Klemp & Durran 1983; Bougeault 1983 for the hydrostatic case.) phi' is the
geopotential PERTURBATION -- deviation from the horizontal mean -- so k = 0 is
excluded automatically and the mean state is untouched. |k| = 0 gives no flux,
which is correct: a horizontally uniform column has nothing to radiate.

Converting to the mass flux this model's continuity equation wants:

    p = p_top + sigma * pi,  so at sigma = 0:  omega_top = pi * sigmadot_top
    omega = -rho g w

    F(k) = omega_top(k) = -rho_top g (|k| / N) phi_hat'(k)

which is in Pa/s and goes straight into `sigma.continuity(top_flux=F)`.

THE SIGN IS NOT DECIDED BY ARGUMENT

Get it backwards and the boundary becomes a wave SOURCE rather than a sink --
it pumps energy in at exactly the rate it should be letting it out, and the
run looks unstable rather than obviously mis-signed. The convention depends on
the sign of the assumed time dependence, which is a choice, so the sign here
is a parameter and `test_radiation.py` MEASURES which value absorbs by
integrating the same mountain wave both ways.
"""

import numpy as np

from sigma import RD, G0, P0, KAPPA


def horizontal_wavenumber(grid, shape):
    """|k| in 1/m for each Fourier mode of a (ny, nx) field."""
    ny, nx = shape
    kx = 2 * np.pi * np.fft.fftfreq(nx, d=grid.dx)
    ky = 2 * np.pi * np.fft.fftfreq(ny, d=grid.dy)
    KX, KY = np.meshgrid(kx, ky)
    return np.sqrt(KX ** 2 + KY ** 2)


def brunt_vaisala(theta, pi, lev, k_top=4):
    """
    A single representative N for the radiation condition, from the layers
    just below the lid.

    One number, not a field: the condition is applied in Fourier space, where
    a spatially varying N would not factor out of the transform. Averaging the
    top few layers is the usual compromise and is where the outgoing wave
    actually leaves.
    """
    p = lev.pressure(pi)
    T = theta * (p / P0) ** KAPPA
    k_top = min(k_top, lev.nz - 1)

    th = theta[:k_top + 1]
    dth = th[1:] - th[:-1]
    th_half = 0.5 * (th[:-1] + th[1:])
    T_half = 0.5 * (T[:k_top + 1][:-1] + T[:k_top + 1][1:])
    dz = RD * T_half / G0 * np.log(p[1:k_top + 1] / p[:k_top])
    dz = np.maximum(np.abs(dz), 1.0)

    # z increases as the index decreases, so d(theta)/dz flips sign.
    N2 = -(G0 / th_half) * dth / dz
    return float(np.sqrt(max(np.nanmean(N2), 1e-8)))


def edge_taper(shape, width):
    """
    Cosine taper to zero over the outermost `width` cells.

    The radiation condition is evaluated with an FFT, which assumes the domain
    is periodic. This one is not: `edge_mode="replicate"`, with a mountain in
    it, so the field has a step across the wrap and multiplying by |k| turns
    that step into Gibbs ringing concentrated at the boundary.

    `remove_divergence_spectral` has the same problem and gets away with it
    because the lateral relaxation zone overwrites the outermost cells. The
    radiation flux has no such protection -- it goes straight into prognostic
    surface pressure -- and the measurement is unambiguous (P-50):

        hour              1      2      3
        |F| edge/inner  1.49   2.28  42.32   -> non-finite at hour 4

    The interior flux was flat at 0.67-0.79 Pa/s the whole time. The boundary
    was not radiating a wave; it was amplifying its own transform error.
    """
    ny, nx = shape
    if width <= 0:
        return np.ones(shape)

    def ramp(n):
        w = np.ones(n)
        k = min(width, n // 2)
        edge = 0.5 * (1 - np.cos(np.pi * (np.arange(k) + 0.5) / k))
        w[:k] = edge
        w[-k:] = edge[::-1]
        return w

    return ramp(ny)[:, None] * ramp(nx)[None, :]


def spectral_rolloff(kmag, grid, cutoff_dx=8.0, rolloff_dx=16.0):
    """
    Raised-cosine response: 0 at or below `cutoff_dx` grid cells, 1 at or
    above `rolloff_dx`, smooth in between. A hard cut rings in physical space
    and puts back some of what it removes.
    """
    dmin = min(grid.dx, grid.dy)
    with np.errstate(divide="ignore"):
        wl = np.where(kmag > 0, 2 * np.pi / np.maximum(kmag, 1e-30), np.inf)
    wl_cells = wl / dmin

    resp = np.ones_like(wl_cells)
    resp[wl_cells <= cutoff_dx] = 0.0
    band = (wl_cells > cutoff_dx) & (wl_cells < rolloff_dx)
    x = (wl_cells[band] - cutoff_dx) / (rolloff_dx - cutoff_dx)
    resp[band] = 0.5 * (1 - np.cos(np.pi * x))
    return resp


def radiative_top_flux(phi, theta, pi, lev, grid, sign=-1.0, n_min=1e-3,
                       reference=None, pi_reference=None, taper=10,
                       cutoff_dx=8.0, rolloff_dx=16.0):
    """
    Mass flux through the lid, in Pa/s, for a non-reflecting upper boundary.

    Returns a (ny, nx) field suitable for `continuity(top_flux=...)`.

    `sign` is +1 or -1 and selects which of the two vertically propagating
    solutions the boundary admits. One of them radiates energy out; the other
    pumps it in. Which is which depends on a convention, so it is measured
    rather than asserted -- see `test_radiation.py`.

    `reference` is a slowly-varying geopotential field to radiate the
    DEVIATION from, and over terrain it is not optional.

    WHY. Applied to the raw perturbation, this condition took a 2500 m terrain
    run from 9/12 hours to 3/12, and halving the timestep changed nothing, so
    it was not a CFL problem. The cause is that over a mountain the top-level
    geopotential perturbation is dominated by the terrain's own steady
    hydrostatic imprint -- a balanced anomaly that sits there and does not
    propagate. Radiating it pumps mass out of the columns over the mountain,
    continuously, for the whole run.

    A steady anomaly has zero intrinsic frequency and should radiate nothing.
    Subtracting a running low-pass leaves the transient part, which is the
    part that actually propagates -- the same distinction the Rayleigh sponge
    could not make, made here on the quantity that carries the flux.
    """
    phi_top = np.asarray(phi[0], dtype=float)
    perturbation = phi_top - phi_top.mean()
    if reference is not None:
        ref = np.asarray(reference, dtype=float)
        perturbation = perturbation - (ref - ref.mean())

    # REMOVE THE COLUMN'S OWN HYDROSTATIC RESPONSE TO SURFACE PRESSURE.
    #
    # This is the closed loop that killed the terrain case, and it is worth
    # writing out because tightening knobs would never have found it:
    #
    #     F  ->  dpi/dt  ->  pi  ->  phi_top  ->  F
    #
    # phi_top is the hydrostatic integral from the ground up, so it moves when
    # pi moves, by about R T / p_s per pascal -- roughly 0.94 m^2/s^2 per Pa
    # here. With F = rho g |k| phi' / N and rho g / N about 164, the loop gain
    # is 154 |k| per second: an e-folding of about 130 s for a 120 km wave.
    # And the sign that correctly radiates waves is the sign that makes THIS
    # loop grow, so the two cannot be satisfied by choosing a sign.
    #
    # Measured: interior |F| growing 1.26 -> 2.17 -> non-finite, with the edge
    # artifact already fixed and a grid-scale cutoff already applied. Neither
    # touched it, because the loop is not a boundary artifact and not
    # grid-scale.
    #
    # The condition is meant to act on the WAVE's geopotential perturbation.
    # The part of phi'_top that merely follows pi is the column's balanced
    # hydrostatic response, not a wave, so subtracting it both restores the
    # intended meaning and opens the loop.
    # BOTH TERMS MUST BE TRANSIENT. Subtracting the response to the raw pi
    # perturbation makes it worse, not better: over 2500 m terrain pi varies
    # by 27000 Pa because of the MOUNTAIN, so pi' is +/-13000 Pa against a
    # phi' of about 150 m^2/s^2. That subtraction does not open the loop, it
    # injects a terrain-shaped signal two orders of magnitude larger than the
    # wave. Measured: |F| jumped to 25 Pa/s and the run died an hour sooner.
    #
    # The loop runs through CHANGES in pi, so it is the change that has to be
    # removed -- pi against its own running low-pass, matching the treatment
    # phi_top already gets.
    if pi_reference is not None:
        p_s = lev.p_top + np.asarray(pi, dtype=float)
        T_col = float(np.mean(theta * (lev.pressure(pi) / P0) ** KAPPA))
        dphi_dpi = RD * T_col / p_s
        pi_pert = np.asarray(pi, dtype=float) - np.asarray(pi_reference,
                                                           dtype=float)
        pi_pert = pi_pert - float(np.mean(pi_pert))
        perturbation = perturbation - dphi_dpi * pi_pert

    N = max(brunt_vaisala(theta, pi, lev), n_min)
    kmag = horizontal_wavenumber(grid, perturbation.shape)

    p_top = lev.p_top
    T_top = float(np.mean(theta[0] * (lev.pressure(pi)[0] / P0) ** KAPPA))
    rho_top = p_top / (RD * max(T_top, 100.0))

    # Window BEFORE the transform, so the step across the wrap is smaller,
    # and taper the result AFTER it, so whatever ringing survives does not
    # reach the prognostic surface pressure. Both are needed: windowing alone
    # leaves a smaller step, and tapering alone lets the ringing form.
    win = edge_taper(perturbation.shape, taper)

    # HIGH-WAVENUMBER CUTOFF, and it is physics rather than a fudge.
    #
    # The transfer function is proportional to |k|, so the shortest waves get
    # the largest flux -- and the flux changes surface pressure, which changes
    # the geopotential, which changes the flux. That loop has gain rising with
    # |k|, so it goes unstable at the grid scale first. Measured with the edge
    # artifact already fixed (P-50): interior |F| grew 1.25 -> 4.59 Pa/s in one
    # hour and the run died at hour 3.
    #
    # The hydrostatic radiation condition is only valid where the wave is
    # hydrostatic, |k| << N/U, which fails long before the grid scale. So the
    # cutoff removes exactly the wavenumbers the condition was never derived
    # for -- and which this model cannot represent as waves anyway.
    response = spectral_rolloff(kmag, grid, cutoff_dx, rolloff_dx)

    hat = np.fft.fft2(perturbation * win)
    w_hat = (kmag / N) * hat * response
    w = np.real(np.fft.ifft2(w_hat)) * win

    return sign * rho_top * G0 * w


def top_flux_stability_dt(F, pi, lev, safety=0.25):
    """
    The flux is explicit, so it has a timestep limit of its own: it must not
    move more than a fraction of the top layer's mass in one step.
    """
    F = np.abs(np.asarray(F))
    fmax = float(F.max())
    if fmax <= 0:
        return np.inf
    layer_mass = float(lev.dsigma[0] * np.asarray(pi).min())
    return safety * layer_mass / fmax

"""
Initialization filtering: remove grid-scale variance before the first step.

WHY THIS EXISTS

Measured on flat ground with a properly balanced 41 m/s jet, the model
integrates 12 hours cleanly with up to 0.30 m/s of white noise added to the
wind, and fails within 1-7 hours at 0.60 m/s and above. The threshold is
sharp, and it is not a boundary-layer problem: neither Richardson mixing nor
surface drag moves it (both schemes, all four combinations, 1/12 h at 1.2 m/s
noise).

The reason is a rate mismatch that can be stated in numbers. Hyperdiffusion is
tuned to damp a 2dx mode with a 3-hour e-folding (10 800 s in the interior,
18 400 s next to a replicate boundary, measured). White noise puts most of its
variance AT that scale, and nonlinear advection amplifies it faster than
3 hours -- the probe shows the 2dx amplitude doubling in about 20 minutes.
Damping loses the race.

Operational models do not integrate raw analyses for this reason. The analysis
is balanced AND filtered first (digital-filter initialization, incremental
analysis update). `balance_initial_state` already removes the divergent part;
this module removes the part that lives at scales the grid cannot carry.

    balance_initial_state  ->  no spurious gravity waves  (divergence)
    filter_grid_scale      ->  no unresolved variance     (wavenumber)

The two are independent and both are needed. Filtering is done in wavenumber
space with a smooth (raised-cosine) rolloff rather than a hard cut, because a
sharp spectral truncation rings in physical space and puts back some of what
it removed.
"""

import numpy as np


def spectral_lowpass(field, grid, cutoff_dx=4.0, rolloff_dx=8.0):
    """
    Smoothly remove variance at wavelengths shorter than `cutoff_dx` grid
    cells, leaving everything longer than `rolloff_dx` untouched.

    Response is 1 for wavelengths >= rolloff_dx, 0 for <= cutoff_dx, and a
    raised cosine in between. Works on 2D or 3D (level, y, x) fields.
    """
    if cutoff_dx >= rolloff_dx:
        raise ValueError("cutoff_dx must be shorter than rolloff_dx")

    a = np.asarray(field, dtype=float)
    ny, nx = a.shape[-2:]

    # Wavelength in grid cells for each Fourier mode.
    kx = np.fft.fftfreq(nx)            # cycles per grid cell
    ky = np.fft.fftfreq(ny)
    kk = np.sqrt(kx[None, :] ** 2 + ky[:, None] ** 2)
    with np.errstate(divide="ignore"):
        wl = np.where(kk > 0, 1.0 / np.maximum(kk, 1e-30), np.inf)

    resp = np.ones_like(wl)
    resp[wl <= cutoff_dx] = 0.0
    band = (wl > cutoff_dx) & (wl < rolloff_dx)
    x = (wl[band] - cutoff_dx) / (rolloff_dx - cutoff_dx)
    resp[band] = 0.5 * (1 - np.cos(np.pi * x))

    out = np.real(np.fft.ifft2(np.fft.fft2(a, axes=(-2, -1)) * resp,
                               axes=(-2, -1)))
    return out


def grid_scale_energy(field, grid, cutoff_dx=4.0):
    """
    Fraction of the field's variance living at wavelengths <= cutoff_dx.

    This is the number to look at before deciding whether an initial state
    needs filtering: a balanced analysis carries a few percent, white noise
    carries most of its variance there.
    """
    a = np.asarray(field, dtype=float)
    a = a - a.mean()
    total = float((a ** 2).sum())
    if total <= 0:
        return 0.0
    smooth = spectral_lowpass(a, grid, cutoff_dx, cutoff_dx * 2.0)
    return float(((a - smooth) ** 2).sum() / total)


def filter_initial_state(u, v, theta, grid, cutoff_dx=4.0, rolloff_dx=8.0,
                         filter_theta=True):
    """
    Apply the lowpass to an initial state, returning filtered copies.

    Temperature is filtered with the wind by default: leaving grid-scale
    structure in theta while removing it from the wind leaves the state
    unbalanced at exactly the scales that were just cleaned, and the
    pressure-gradient force puts the noise straight back into the wind.
    """
    uf = spectral_lowpass(u, grid, cutoff_dx, rolloff_dx)
    vf = spectral_lowpass(v, grid, cutoff_dx, rolloff_dx)
    if not filter_theta:
        return uf, vf, theta
    # Filter the DEVIATION from the horizontal mean so the mean profile,
    # which carries the stratification, is untouched.
    ref = theta.mean(axis=(-2, -1), keepdims=True)
    tf = ref + spectral_lowpass(theta - ref, grid, cutoff_dx, rolloff_dx)
    return uf, vf, tf

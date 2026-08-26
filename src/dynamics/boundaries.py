"""
Davies relaxation for a limited-area domain.

A regional model has open edges. Something must tell the solution what is
arriving from outside, or waves reflect off the boundary and contaminate the
interior within hours. The standard treatment (Davies 1976) is a relaxation
zone: in a band of cells around the perimeter, the model state is nudged
toward externally supplied values after every step,

    phi <- phi + alpha * (phi_ext - phi)

with alpha ~= 1 at the outermost row (essentially prescribed) falling smoothly
to 0 at the inner edge of the zone. The smooth taper is the whole point: a
hard switch from prescribed to free is itself a discontinuity, and reflects.

The external values come from a coarser or larger-domain model -- for this
project, HRRR. That dependency is not a weakness; every operational regional
model on earth is driven this way. It is fundamentally different from LEARNING
the driving model, because only the boundary is imposed and the interior
evolves under its own physics.
"""

import numpy as np


def relaxation_weights(nx, ny, width=10, alpha_max=1.0, profile="cosine"):
    """
    Build the alpha field: alpha_max at the perimeter, 0 inside the zone.

    width   : zone thickness in cells. 8-15 is typical; too thin reflects,
              too thick wastes domain.
    profile : "cosine" (smooth, recommended) or "quadratic" or "linear".
    """
    if width < 1:
        return np.zeros((ny, nx))
    if 2 * width >= min(nx, ny):
        raise ValueError(f"relaxation zone {width} too wide for {nx}x{ny} grid")

    ix = np.arange(nx)
    iy = np.arange(ny)
    # Distance (in cells) from the nearest edge.
    dist_x = np.minimum(ix, nx - 1 - ix)
    dist_y = np.minimum(iy, ny - 1 - iy)
    dist = np.minimum(dist_x[None, :], dist_y[:, None]).astype(float)

    s = np.clip(1.0 - dist / width, 0.0, 1.0)      # 1 at edge -> 0 at zone edge

    if profile == "cosine":
        w = 0.5 * (1.0 - np.cos(np.pi * s))        # smooth at both ends
    elif profile == "quadratic":
        w = s ** 2
    elif profile == "linear":
        w = s
    else:
        raise ValueError(f"unknown profile {profile!r}")

    return alpha_max * w


class DaviesRelaxation:
    """Applies the relaxation after each model step."""

    def __init__(self, grid, width=10, alpha_max=1.0, profile="cosine"):
        self.grid = grid
        self.width = width
        self.alpha = relaxation_weights(grid.nx, grid.ny, width,
                                        alpha_max, profile)

    def apply(self, model, ext):
        """
        ext : dict with any of "u", "v", "h" -- the external state to relax
              toward. Missing keys are left alone.
        """
        a = self.alpha
        if "u" in ext:
            model.u += a * (ext["u"] - model.u)
        if "v" in ext:
            model.v += a * (ext["v"] - model.v)
        if "h" in ext:
            model.h += a * (ext["h"] - model.h)

    @property
    def interior_fraction(self):
        """Fraction of the domain not touched by relaxation."""
        return float((self.alpha == 0).sum()) / self.alpha.size

    def __repr__(self):
        return (f"DaviesRelaxation(width={self.width}, "
                f"interior={self.interior_fraction:.1%})")


class BoundaryDriver:
    """
    Holds a time series of external states and interpolates between them.

    Driving data arrives hourly (HRRR); the model steps every few seconds, so
    boundary values must be interpolated in time. Without interpolation the
    boundary jumps once an hour and radiates a spurious wave inward on every
    update.
    """

    def __init__(self, times, states):
        """
        times  : increasing seconds from model t=0
        states : list of dicts with "u", "v", "h" arrays
        """
        if len(times) != len(states):
            raise ValueError("times and states must be the same length")
        if len(times) < 1:
            raise ValueError("need at least one external state")
        order = np.argsort(times)
        self.times = np.asarray(times, dtype=float)[order]
        self.states = [states[i] for i in order]

    def at(self, t):
        """Linearly interpolated external state at time t (seconds)."""
        if len(self.times) == 1 or t <= self.times[0]:
            return self.states[0]
        if t >= self.times[-1]:
            return self.states[-1]

        j = int(np.searchsorted(self.times, t))
        t0, t1 = self.times[j - 1], self.times[j]
        w = (t - t0) / (t1 - t0)
        a, b = self.states[j - 1], self.states[j]
        return {k: (1.0 - w) * a[k] + w * b[k] for k in a if k in b}

    @property
    def span_hours(self):
        return (self.times[-1] - self.times[0]) / 3600.0

    def __repr__(self):
        return (f"BoundaryDriver({len(self.times)} frames, "
                f"{self.span_hours:.1f} h)")


def run_limited_area(model, driver, relax, duration, dt=None, callback=None,
                     every=0):
    """
    Integrate with boundary relaxation applied after every step.

    Relaxation is applied to the state, not added as a tendency: it is a
    correction imposed after the dynamics have acted, which is how Davies
    forcing is normally implemented and keeps it independent of the time
    scheme.
    """
    dt = dt or model.max_dt()
    n = int(np.ceil(duration / dt))
    dt = duration / n

    for k in range(n):
        model.step(dt)
        relax.apply(model, driver.at(model.time))
        if callback and every and (k + 1) % every == 0:
            callback(model)
    return n

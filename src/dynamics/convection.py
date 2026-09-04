"""
Dry convective adjustment.

WHY, WITH THE MEASUREMENT THAT MOTIVATED IT

Over 2500 m terrain with a clean, filtered, balanced 41 m/s flow, the model
survives 11 of 12 hours. The wind never runs away -- it sits at 41 m/s to the
end. What runs away is the STRATIFICATION:

    hour     1     2     3     4     6     7     8     9    10    11    12
    min Ri  11.5   3.8   2.1   1.6  0.94  0.47  0.33  0.23 -0.05 -1.15  dead
    Ri<0.25    0     0     0     0     0     0     0     1     7    57

A negative Richardson number means N^2 < 0: potential temperature decreasing
with height. The mountain wave steepens as it propagates upward, overturns,
and the model has nothing that removes the resulting statically unstable
layer. This is real physics -- mountain waves do break -- and the failure is
the absence of the process that follows breaking, not a numerical defect.

WHY THE MIXING SCHEME DOES NOT COVER IT

`turbulence.eddy_diffusivity` does treat Ri <= 0 as full-strength mixing, but
it is a diffusion with K capped at 100 m^2/s. Across a 600 m layer that is a
relaxation time of dz^2/K = 3600 s. The wave steepens faster than an hour, so
diffusion loses the race exactly as hyperdiffusion lost the race against
grid-scale noise. Convective overturning is not a slow diffusion; it is a
rearrangement that happens on the eddy turnover time, which at these scales is
effectively instant.

WHAT THIS DOES

Wherever theta decreases with height, the unstable layers are mixed to a
common mass-weighted mean -- neutral stratification, enthalpy conserved. Wind
is mixed over the same layers with the same weights, so momentum is conserved
too and the scheme carries convective momentum transport rather than leaving
the wind untouched in a column that has just overturned.

It is applied as a POST-STEP ADJUSTMENT, not as a tendency. An adjustment that
enforces an inequality has no meaningful time derivative, and putting it
inside the Runge-Kutta stages would let the intermediate states re-create the
instability the final state is supposed to be free of.

LIMITS, STATED PLAINLY

This is the dry, hard-adjustment form: no entrainment, no cloud, no
downdrafts, no mass flux. It removes static instability and nothing else. When
moisture arrives it will need to be replaced by something that carries latent
heat, not extended.
"""

import numpy as np


def unstable_fraction(theta, tol=1e-10):
    """Fraction of interfaces where theta decreases with height (index 0 = top)."""
    return float(np.mean(theta[:-1] < theta[1:] - tol))


def dry_convective_adjustment(theta, u, v, pi, lev, max_sweeps=20,
                              mix_momentum=True, tol=1e-10):
    """
    Remove static instability by mass-weighted mixing of adjacent layers.

    Arrays are (nz, ny, nx) with index 0 at the model top. A column is stable
    when theta DECREASES with index. Returns (theta, u, v, info) with new
    arrays; the inputs are not modified.

    The sweep is repeated because mixing one pair can destabilise the pair
    above or below it. Convergence is monotone -- each mix strictly reduces
    the number of unstable interfaces or leaves it unchanged -- so the sweep
    cap is a guard, not a tuning knob.
    """
    theta = np.array(theta, dtype=float, copy=True)
    u = np.array(u, dtype=float, copy=True)
    v = np.array(v, dtype=float, copy=True)

    # Layer mass per unit area: dp/g, and dp = dsigma * pi.
    dm = lev.dsigma.reshape(-1, *([1] * np.ndim(pi))) * np.asarray(pi)
    dm = np.broadcast_to(dm, theta.shape).astype(float)

    before = unstable_fraction(theta, tol)
    nz = theta.shape[0]
    fields = [theta] + ([u, v] if mix_momentum else [])
    sweeps = 0

    for sweeps in range(1, max_sweeps + 1):
        bad = theta[:-1] < theta[1:] - tol
        if not bad.any():
            sweeps -= 1
            break

        # CONTIGUOUS-SEGMENT MIXING, not pairwise.
        #
        # A layer belongs to a mixing segment if the interface above it or the
        # interface below it is unstable. Contiguous segments are disjoint by
        # construction, so mixing each one to its mass-weighted mean is exact
        # and conserves the column integral to round-off.
        #
        # Pairwise mixing was tried first and is conservative but converges
        # like a diffusion: a fully inverted 20-level column still had 0.26 K
        # of spread after 200 sweeps. Segment mixing settles it in one.
        member = np.zeros(theta.shape, dtype=bool)
        member[:-1] |= bad
        member[1:] |= bad

        # Segment id: increments whenever a new segment starts, going down.
        starts = member.copy()
        starts[1:] &= ~member[:-1]
        seg = np.cumsum(starts, axis=0)

        # Forward pass: accumulate mass and mass-weighted sums per segment.
        acc_m = np.zeros(theta.shape[1:])
        acc = [np.zeros(theta.shape[1:]) for _ in fields]
        sums_m = np.zeros(theta.shape)
        sums = [np.zeros(theta.shape) for _ in fields]
        for k in range(nz):
            reset = starts[k]
            acc_m = np.where(reset, dm[k], acc_m + dm[k] * member[k])
            for j, a in enumerate(fields):
                acc[j] = np.where(reset, dm[k] * a[k],
                                  acc[j] + dm[k] * a[k] * member[k])
            sums_m[k] = acc_m
            for j in range(len(fields)):
                sums[j][k] = acc[j]

        # Backward pass: the last layer of a segment holds the totals; carry
        # them back up to every layer in the same segment.
        tot_m = np.zeros(theta.shape)
        tot = [np.zeros(theta.shape) for _ in fields]
        carry_m = np.zeros(theta.shape[1:])
        carry = [np.zeros(theta.shape[1:]) for _ in fields]
        last = member.copy()
        last[:-1] &= ~member[1:]
        for k in range(nz - 1, -1, -1):
            take = last[k]
            carry_m = np.where(take, sums_m[k], carry_m)
            for j in range(len(fields)):
                carry[j] = np.where(take, sums[j][k], carry[j])
            tot_m[k] = carry_m
            for j in range(len(fields)):
                tot[j][k] = carry[j]

        with np.errstate(invalid="ignore", divide="ignore"):
            for j, a in enumerate(fields):
                mean = np.where(tot_m > 0, tot[j] / np.maximum(tot_m, 1e-30),
                                a)
                a[...] = np.where(member, mean, a)

    info = {
        "unstable_before": before,
        "unstable_after": unstable_fraction(theta, tol),
        "sweeps": sweeps,
    }
    return theta, u, v, info

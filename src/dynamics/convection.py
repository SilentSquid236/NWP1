"""
Dry convective adjustment (after Manabe et al. 1965).

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

from backend import xp_of, to_numpy


def unstable_fraction(theta, tol=1e-10):
    """Fraction of interfaces where theta decreases with height (index 0 = top)."""
    return float(xp_of(theta).mean(theta[:-1] < theta[1:] - tol))


def dry_convective_adjustment(theta, u, v, pi, lev, max_sweeps=20,
                              mix_momentum=True, tol=1e-10):
    """
    Remove static instability by mass-weighted mixing of adjacent layers.

    Arrays are (nz, ny, nx) with index 0 at the model top. A column is stable
    when theta DECREASES with index. Returns (theta, u, v, info) with new
    arrays; the inputs are not modified.

    The sweep is repeated because mixing one pair can destabilise the pair
    above or below it. Convergence is monotone -- each mix strictly reduces
    the number of unstable interfaces or leaves it unchanged.

    THE CAP BINDS (P-63, test AB). This was written as "the sweep cap is a
    guard, not a tuning knob", and it is not true. Once a column has been
    mixed to neutral, every later step nudges the neutral layers by ~1e-4 K,
    and clearing that chain takes 23-34 sweeps. At the cap of 20 every call
    in the calm cases ran to the cap, left instability behind, and cost ~70 %
    of the step. `dry_convective_adjustment_pav` below reaches the converged
    answer in one pass.
    """
    xp = xp_of(theta, u, v, pi)
    if xp.name == "numpy":
        theta = np.array(theta, dtype=float, copy=True)
        u = np.array(u, dtype=float, copy=True)
        v = np.array(v, dtype=float, copy=True)
    else:
        theta, u, v = xp.copy(theta), xp.copy(u), xp.copy(v)

    # Layer mass per unit area: dp/g, and dp = dsigma * pi.
    pi_a = xp.asarray(pi)
    dm = xp.asarray(lev.dsigma).reshape(-1, *([1] * xp.ndim(pi_a))) * pi_a
    dm = xp.astype(xp.broadcast_to(dm, tuple(theta.shape)), float)

    before = unstable_fraction(theta, tol)
    nz = theta.shape[0]
    fields = [theta] + ([u, v] if mix_momentum else [])
    sweeps = 0

    for sweeps in range(1, max_sweeps + 1):
        bad = theta[:-1] < theta[1:] - tol
        if not bool(bad.any()):
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
        member = xp.zeros(tuple(theta.shape), dtype=bool)
        member[:-1] |= bad
        member[1:] |= bad

        # Segment id: increments whenever a new segment starts, going down.
        starts = xp.copy(member)
        starts[1:] &= ~member[:-1]

        # Forward pass: accumulate mass and mass-weighted sums per segment.
        shape2, shape3 = tuple(theta.shape[1:]), tuple(theta.shape)
        acc_m = xp.zeros(shape2)
        acc = [xp.zeros(shape2) for _ in fields]
        sums_m = xp.zeros(shape3)
        sums = [xp.zeros(shape3) for _ in fields]
        for k in range(nz):
            reset = starts[k]
            acc_m = xp.where(reset, dm[k], acc_m + dm[k] * member[k])
            for j, a in enumerate(fields):
                acc[j] = xp.where(reset, dm[k] * a[k],
                                  acc[j] + dm[k] * a[k] * member[k])
            sums_m[k] = acc_m
            for j in range(len(fields)):
                sums[j][k] = acc[j]

        # Backward pass: the last layer of a segment holds the totals; carry
        # them back up to every layer in the same segment.
        tot_m = xp.zeros(shape3)
        tot = [xp.zeros(shape3) for _ in fields]
        carry_m = xp.zeros(shape2)
        carry = [xp.zeros(shape2) for _ in fields]
        last = xp.copy(member)
        last[:-1] &= ~member[1:]
        for k in range(nz - 1, -1, -1):
            take = last[k]
            carry_m = xp.where(take, sums_m[k], carry_m)
            for j in range(len(fields)):
                carry[j] = xp.where(take, sums[j][k], carry[j])
            tot_m[k] = carry_m
            for j in range(len(fields)):
                tot[j][k] = carry[j]

        with xp.errstate(invalid="ignore", divide="ignore"):
            for j, a in enumerate(fields):
                mean = xp.where(tot_m > 0, tot[j] / xp.maximum(tot_m, 1e-30),
                                a)
                a[...] = xp.where(member, mean, a)

    info = {
        "unstable_before": before,
        "unstable_after": unstable_fraction(theta, tol),
        "sweeps": sweeps,
    }
    return theta, u, v, info


def dry_convective_adjustment_pav(theta, u, v, pi, lev, mix_momentum=True,
                                  tol=1e-10):
    """
    The same adjustment as `dry_convective_adjustment`, exact in one pass.

    POOL-ADJACENT-VIOLATORS (Ayer et al. 1955). Walk each column from the
    ground up, keeping a stack of mixed blocks. Push the next layer as a block
    of its own; while the block on top is colder than the block below it (by
    more than `tol`), merge the two into their mass-weighted mean. At the lid
    the blocks are stable by construction, so there is nothing left to sweep.

    The answer is the stable state closest to the input (mass-weighted least
    squares): the least mixing that removes the instability. It is NOT where
    the sweep scheme ends. The sweep scheme merges every contiguous run of
    layers next to an unstable interface, across the stable interfaces between
    them, so it mixes more: 16 K apart on a column with 3 K noise, 2e-4 K on a
    near-neutral one (measured 2026-09-29; test_convection.py checks the
    least-change property). u and v are mixed over the same blocks, and mass-
    weighted theta, u and v are conserved as in the sweep scheme.

    Only columns with an unstable interface are touched, and the work is on
    NumPy arrays (for torch CPU tensors the conversion is a copy of those few
    columns). A stable state comes back as unchanged copies, value for value.
    Returns (theta, u, v, info) with new arrays; the inputs are not modified.
    """
    xp = xp_of(theta, u, v, pi)
    torch_in = xp.name != "numpy"
    th = np.array(to_numpy(theta), dtype=float, copy=True)
    uu = np.array(to_numpy(u), dtype=float, copy=True)
    vv = np.array(to_numpy(v), dtype=float, copy=True)
    pa = np.asarray(to_numpy(pi), dtype=float)

    nz = th.shape[0]
    t2, u2, v2 = th.reshape(nz, -1), uu.reshape(nz, -1), vv.reshape(nz, -1)
    bad = t2[:-1] < t2[1:] - tol
    before = float(bad.mean())
    cb = np.nonzero(bad.any(axis=0))[0]          # columns needing work

    if cb.size:
        n = cb.size
        ar = np.arange(n)
        w = np.asarray(lev.dsigma, float)[:, None] * pa.reshape(1, -1)[:, cb]
        fields = [t2] + ([u2, v2] if mix_momentum else [])
        S_w = np.zeros((nz, n))
        S = [np.zeros((nz, n)) for _ in fields]  # mass-weighted sums
        start = np.zeros((nz, n), dtype=np.int64)
        top = np.full(n, -1, dtype=np.int64)

        for j in range(nz):                       # j = 0 is the lowest layer
            k = nz - 1 - j
            top += 1
            wk = w[k]
            S_w[top, ar] = wk
            for s, f in zip(S, fields):
                s[top, ar] = wk * f[k, cb]
            start[top, ar] = j
            act = ar
            while True:
                act = act[top[act] >= 1]
                if act.size == 0:
                    break
                t = top[act]
                upper = S[0][t, act] / S_w[t, act]
                lower = S[0][t - 1, act] / S_w[t - 1, act]
                viol = upper < lower - tol
                act, t = act[viol], t[viol]
                if act.size == 0:
                    break
                S_w[t - 1, act] += S_w[t, act]
                for s in S:
                    s[t - 1, act] += s[t, act]
                top[act] -= 1

        # Write the block means back, walking up again: a column moves to its
        # next block at the level where that block starts.
        b = np.zeros(n, dtype=np.int64)
        for j in range(nz):
            k = nz - 1 - j
            nxt = np.minimum(b + 1, nz - 1)
            b = np.where((b + 1 <= top) & (start[nxt, ar] == j), b + 1, b)
            for s, f in zip(S, fields):
                f[k, cb] = s[b, ar] / S_w[b, ar]

    after = float((t2[:-1] < t2[1:] - tol).mean())
    info = {"unstable_before": before, "unstable_after": after,
            "sweeps": 1 if cb.size else 0, "columns": int(cb.size),
            "scheme": "pav"}
    if torch_in:
        import torch
        th, uu, vv = (torch.from_numpy(a) for a in (th, uu, vv))
    return th, uu, vv, info

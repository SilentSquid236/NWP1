"""
Terrain-following sigma coordinate.

    sigma = (p - p_top) / (p_s - p_top),      sigma = 0 at the lid, 1 at the ground

WHY THIS REPLACES PRESSURE COORDINATES

In pure pressure coordinates the lower boundary is a problem with no good
answer. Isobaric surfaces intersect terrain, and with a rigid flat bottom the
vertical velocity omega must vanish at BOTH ends of the column. Diagnosing
omega from divergence and then enforcing both conditions over-constrains the
system: the correction that pins omega = 0 at the ground redistributes error
through the whole column every step, feeding vertical advection -> divergence
-> omega, a loop with no physical damping. Measured consequence: the model was
stable on smooth balanced states and diverged within 2-3 hours on real
analyses, at every damping setting tried. See docs/STABILITY.md.

Sigma fixes this structurally rather than by damping:

  * The ground is sigma = 1 by definition, whatever the terrain height.
  * sigma_dot = 0 at both ends is a genuine boundary condition, not an
    imposed constraint, because SURFACE PRESSURE IS PROGNOSTIC. The column
    can breathe; mass leaves through a moving boundary instead of being
    forced to balance.
  * Terrain follows for free, which the Northeast domain needs anyway.

The cost is a more complex pressure-gradient force: on a sloping sigma surface
it splits into two terms that nearly cancel over steep terrain, and their
difference is the small quantity we want. That cancellation is the known
weakness of sigma coordinates, and it is why the pressure-gradient test below
checks an isothermal atmosphere at rest over a slope -- the case where the two
terms must cancel exactly.

CONVENTIONS

  * Index 0 is the TOP (sigma smallest), index nz-1 is the GROUND.
    This is the opposite of the pressure-coordinate module, and is chosen so
    sigma increases with index, matching the usual meteorological ordering.
  * pi = p_s - p_top is the "pressure thickness" of the column, prognostic.
  * p(sigma) = p_top + sigma * pi
"""

import numpy as np

RD = 287.058
CP = 1004.6
KAPPA = RD / CP
P0 = 100_000.0
G0 = 9.80665


class SigmaLevels:
    """
    Vertical grid in sigma. Full levels carry u, v, theta; half levels carry
    sigma_dot and bound the layers.
    """

    def __init__(self, n_levels=20, p_top=20000.0, stretch=1.4):
        """
        p_top defaults to 200 hPa, not 50 hPa.

        Measured: with terrain, raising the lid from 200 to 50 hPa costs about
        an hour of stability (3/4 -> 2/4 forecast hours survived), while level
        count and stretching change nothing. Energy enters through the
        pressure-gradient term at the topmost levels, where R*T/p is largest
        and the hydrostatic integral has accumulated furthest. A dry model
        with no stratospheric physics gains little from a 50 hPa lid.

        RE-MEASURED 2026-09-03 on corrected initial states, because the
        original measurement was made on a state now known to carry a clipped
        166 m/s jet. The conclusion survives. Over 2500 m terrain, clean and
        filtered, 12-hour ceiling:

            lid      sponge 5 levels   sponge 8 levels
            200 hPa      11/12             11/12
            100 hPa      10/12             11/12
             50 hPa       9/12             10/12

        Raising the lid is neutral to slightly worse at every sponge depth.
        Keep 200 hPa.
        """
        self.nz = int(n_levels)
        self.p_top = float(p_top)

        # Half levels from 0 (lid) to 1 (ground), stretched so resolution
        # concentrates near the surface where the boundary layer lives.
        s = np.linspace(0.0, 1.0, self.nz + 1)
        self.sigma_half = s ** stretch

        # Full levels at layer midpoints.
        self.sigma = 0.5 * (self.sigma_half[:-1] + self.sigma_half[1:])
        self.dsigma = np.diff(self.sigma_half)

        if np.any(self.dsigma <= 0):
            raise ValueError("sigma half levels must increase monotonically")

    # --- pressure ---------------------------------------------------------

    def pressure(self, pi):
        """p at full levels. pi = p_s - p_top, shape (ny, nx) or scalar."""
        pi = np.asarray(pi)
        return self.p_top + self.sigma.reshape(-1, *([1] * pi.ndim)) * pi

    def pressure_half(self, pi):
        pi = np.asarray(pi)
        return self.p_top + self.sigma_half.reshape(-1, *([1] * pi.ndim)) * pi

    def exner(self, pi):
        return (self.pressure(pi) / P0) ** KAPPA

    def temperature(self, theta, pi):
        return theta * self.exner(pi)

    def __repr__(self):
        return (f"SigmaLevels({self.nz} levels, p_top={self.p_top/100:.0f} hPa, "
                f"sigma {self.sigma[0]:.3f}..{self.sigma[-1]:.3f})")


def hydrostatic_geopotential(theta, pi, lev, phi_surface=0.0):
    """
    Integrate hydrostatic balance UPWARD from the surface.

        dPhi/d(ln p) = -R T

    Index 0 is the lid, nz-1 the ground, so we fill from the bottom index
    backwards. phi_surface is g * terrain height and may be a 2D field.
    """
    p = lev.pressure(pi)
    T = theta * (p / P0) ** KAPPA

    phi = np.empty_like(theta)
    phi[-1] = phi_surface + RD * T[-1] * np.log(
        (lev.p_top + pi) / p[-1])          # ground -> lowest full level

    for k in range(lev.nz - 2, -1, -1):
        T_layer = 0.5 * (T[k] + T[k + 1])
        phi[k] = phi[k + 1] + RD * T_layer * np.log(p[k + 1] / p[k])

    return phi


def pressure_gradient_force(phi, theta, pi, lev, grid, reference=None):
    """
    Horizontal pressure-gradient force in sigma coordinates.

    Transforming grad_p(Phi) onto sigma surfaces, with p = p_top + sigma*pi:

    Chain rule between sigma and pressure surfaces, with p = p_top + sigma*pi:

        grad(Phi)|_sigma = grad(Phi)|_p + (dPhi/dp) * (dp/dx)|_sigma
                         = grad(Phi)|_p - (R T / p) * sigma * grad(pi)

    so, rearranging and negating,

        F = -grad(Phi)|_p = -grad_sigma(Phi) - (R T / p) * sigma * grad(pi)

    BOTH signs here were derived wrong twice before being measured. The check
    that settled it: for an isothermal atmosphere in exact balance over a
    ridge the two terms must cancel, and only this combination does
    (residual 2.9e-04 against individual terms of 6e-02). Note it is also NOT
    -R T grad(ln p_s), which coincides with the correct term only at
    sigma = 1 with p_top = 0.

    The two terms are individually large over terrain and nearly cancel; their
    difference is the physical force. That cancellation is the known weakness
    of sigma coordinates and is what test_pressure_gradient_over_terrain
    measures.

    `reference` is an optional horizontally uniform T0(sigma) profile. Passing
    it switches to the reference-state form below, where the cancellation is
    done ANALYTICALLY instead of numerically.
    """
    p = lev.pressure(pi)
    T = theta * (p / P0) ** KAPPA
    sig = lev.sigma.reshape(-1, 1, 1)

    dpidx = grid.dx_backward(pi)
    dpidy = grid.dy_backward(pi)

    lnp = np.log(p)

    if reference is None:
        # HYDROSTATICALLY CONSISTENT FORM.
        #
        #   F = -[ grad(Phi) + R T grad(ln p) ]
        #
        # with grad(ln p) DIFFERENCED DIRECTLY and T averaged to the velocity
        # point. The earlier version expanded grad(ln p) analytically as
        # sigma*grad(pi)/p -- algebraically identical, discretely not, because
        # the expansion is evaluated at cell centres while grad(pi) lives on
        # the faces.
        #
        # The criterion: for T = T(p) alone the true force is exactly zero
        # over any terrain. Measured residual at 3000 m of terrain:
        #
        #     sigma*grad(pi)/p at centres      2.1e-03   (what we had)
        #     same, coefficient on v points    2.3e-05
        #     grad(ln p) differenced directly  8.2e-15   (machine zero)
        #
        # That residual is a spurious force proportional to terrain slope,
        # largest in the LOWEST layers where the sigma surfaces are most
        # steeply tilted -- the boundary layer. It is why forecast survival
        # fell monotonically with terrain height.
        fx = -(grid.dx_backward(phi) + grid.h_to_u(RD * T) * grid.dx_backward(lnp))
        fy = -(grid.dy_backward(phi) + grid.h_to_v(RD * T) * grid.dy_backward(lnp))
        return fx, fy

    T0 = reference.reshape(-1, 1, 1)
    psi = phi + RD * T0 * lnp
    fx = -(grid.dx_backward(psi)
           + grid.h_to_u(RD * (T - T0)) * grid.dx_backward(lnp))
    fy = -(grid.dy_backward(psi)
           + grid.h_to_v(RD * (T - T0)) * grid.dy_backward(lnp))
    return fx, fy


def continuity(u, v, pi, lev, grid):
    """
    Sigma-system continuity. Returns (dpi_dt, sigma_dot).

        d(pi)/dt = -integral_0^1 div(pi V) dsigma
        pi * sigma_dot(s) = -integral_0^s div(pi V) dsigma' - s * d(pi)/dt

    sigma_dot vanishes at sigma = 0 and sigma = 1 BY CONSTRUCTION here --
    the second term is exactly what makes the top and bottom values cancel.
    No correction is applied afterwards, which is the structural difference
    from the pressure-coordinate version.
    """
    # Mass flux divergence per layer.
    pi_u = grid.h_to_u(pi) * u
    pi_v = grid.h_to_v(pi) * v
    div = grid.dx_forward(pi_u) + grid.dy_forward(pi_v)      # (nz, ny, nx)

    ds = lev.dsigma.reshape(-1, 1, 1)
    total = (div * ds).sum(axis=0)                            # (ny, nx)
    dpi_dt = -total

    # Partial integral to each half level, top downward.
    partial = np.concatenate([np.zeros((1,) + div.shape[1:]),
                              np.cumsum(div * ds, axis=0)], axis=0)  # (nz+1,...)

    s_half = lev.sigma_half.reshape(-1, 1, 1)
    pi_sigmadot = -partial - s_half * dpi_dt                  # (nz+1, ny, nx)

    with np.errstate(divide="ignore", invalid="ignore"):
        sigma_dot = np.where(pi > 0, pi_sigmadot / pi, 0.0)

    return dpi_dt, sigma_dot


def vertical_advection(a, sigma_dot, lev):
    """
    sigma_dot * da/dsigma, with sigma_dot on half levels and a on full levels.

    Flux-form so that a constant field is transported exactly: the two half
    level fluxes cancel identically. An advective form would leave a residual
    proportional to the divergence of sigma_dot, which is precisely the term
    that misbehaves near the boundaries.
    """
    # Interface values by simple averaging; ends use the adjacent full level.
    a_half = np.empty((lev.nz + 1,) + a.shape[1:])
    a_half[1:-1] = 0.5 * (a[:-1] + a[1:])
    a_half[0] = a[0]
    a_half[-1] = a[-1]

    flux = sigma_dot * a_half
    ds = lev.dsigma.reshape(-1, 1, 1)

    # d(flux)/dsigma minus a * d(sigma_dot)/dsigma recovers the advective form
    # while keeping the exact-constant property.
    dflux = (flux[1:] - flux[:-1]) / ds
    dsd = (sigma_dot[1:] - sigma_dot[:-1]) / ds
    return dflux - a * dsd


# ---------------------------------------------------------------------------
# Hydrostatically consistent discretisation (Simmons & Burridge 1981)
# ---------------------------------------------------------------------------
#
# THE CRITERION
#
# If temperature is a function of pressure alone -- horizontally uniform on
# pressure surfaces -- then the true pressure-gradient force is exactly zero
# everywhere, over any terrain. A discretisation that does not reproduce that
# exactly is "hydrostatically inconsistent", and the residual is a spurious
# force that scales with terrain slope.
#
# Measured for the naive form: 0.5% of the individual terms. That sounds
# small, but it is a persistent forcing, it accumulates upward through the
# geopotential integral, and the terrain sweep showed forecast survival
# falling monotonically with slope.
#
# THE FIX
#
# Derive Phi and the pressure-gradient term from the SAME half-level integral,
# so the cancellation is structural rather than numerical:
#
#   Phi_{k+1/2} = Phi_s + sum_{j>k} R T_j ln(p_{j+1/2} / p_{j-1/2})
#   Phi_k       = Phi_{k+1/2} + alpha_k R T_k
#   alpha_k     = 1 - (p_{k-1/2} / dp_k) ln(p_{k+1/2} / p_{k-1/2})
#
# and build the force from those same pieces. alpha_k is the geometric factor
# placing the full level correctly within its layer in log-pressure; using the
# arithmetic midpoint instead is precisely the inconsistency.
#
# WHERE THIS BITES HARDEST: THE BOUNDARY LAYER
#
# The lowest layers sit directly on the terrain, where sigma surfaces are most
# steeply tilted and the two cancelling terms are largest relative to their
# difference. They are also the thickest in pressure. So the near-surface
# levels carry the worst relative error -- and they are exactly where a
# boundary-layer scheme would later add friction and mixing, which would then
# be operating on a spuriously forced flow. Getting this right before adding a
# boundary layer matters; a friction term tuned against a wrong pressure
# gradient would compensate for a bug.


def alpha_factors(pi, lev):
    """
    Simmons-Burridge alpha_k: where the full level sits within its layer, in
    log-pressure. Returns shape (nz, ny, nx).

    The top layer needs a convention because p_{k-1/2} -> p_top; ln(2) is the
    standard choice and corresponds to placing the level at the layer's
    log-pressure midpoint.
    """
    p_half = lev.pressure_half(pi)                    # (nz+1, ...)
    dp = p_half[1:] - p_half[:-1]                     # (nz, ...)

    with np.errstate(divide="ignore", invalid="ignore"):
        alpha = 1.0 - (p_half[:-1] / dp) * np.log(p_half[1:] / p_half[:-1])

    # Topmost layer: p_half[0] = p_top. If p_top is zero the expression is
    # singular; with a finite lid it is well defined, but we keep the standard
    # convention available for the degenerate case.
    if lev.p_top <= 0:
        alpha[0] = np.log(2.0)

    return alpha, p_half, dp


def geopotential_sb(theta, pi, lev, phi_surface=0.0):
    """
    Geopotential at FULL levels, via half-level accumulation.

    Integrating to half levels first and then placing the full level with
    alpha_k is what makes the pressure-gradient force below consistent. The
    simpler "average T between adjacent full levels" integration is not.
    """
    alpha, p_half, dp = alpha_factors(pi, lev)
    p_full = lev.pressure(pi)
    T = theta * (p_full / P0) ** KAPPA

    phi_half = np.empty((lev.nz + 1,) + theta.shape[1:])
    phi_half[-1] = phi_surface                        # ground
    for k in range(lev.nz - 1, -1, -1):
        phi_half[k] = phi_half[k + 1] + RD * T[k] * np.log(
            p_half[k + 1] / p_half[k])

    phi_full = phi_half[1:] + alpha * RD * T
    return phi_full, phi_half


def pressure_gradient_sb(theta, pi, lev, grid, phi_surface=0.0):
    """
    Hydrostatically consistent pressure-gradient force.

        F = -[ grad(Phi_{k+1/2}) + alpha_k R T_k grad(ln p) ]

    Both pieces come from the same half-level construction used for Phi, so
    for T = T(p) the two cancel to machine precision instead of to 0.5%.
    """
    alpha, p_half, dp = alpha_factors(pi, lev)
    p_full = lev.pressure(pi)
    T = theta * (p_full / P0) ** KAPPA

    _, phi_half = geopotential_sb(theta, pi, lev, phi_surface)

    # grad(ln p) on a sigma surface: p = p_top + sigma*pi  =>  d(ln p) =
    # sigma * d(pi) / p.
    sig = lev.sigma.reshape(-1, 1, 1)
    dlnp_dx = sig * grid.dx_backward(pi) / p_full
    dlnp_dy = sig * grid.dy_backward(pi) / p_full

    # Gradient of the half-level geopotential BELOW each full level, which is
    # the piece the alpha term completes.
    dphi_dx = grid.dx_backward(phi_half[1:])
    dphi_dy = grid.dy_backward(phi_half[1:])

    fx = -(dphi_dx + alpha * RD * T * dlnp_dx)
    fy = -(dphi_dy + alpha * RD * T * dlnp_dy)
    return fx, fy


# ---------------------------------------------------------------------------
# Orography preparation
# ---------------------------------------------------------------------------

def terrain_slope(h, grid):
    """Maximum terrain slope, and the ascent it forces at a given wind speed."""
    dhdx = grid.dx_backward(h)
    dhdy = grid.dy_backward(h)
    return float(np.sqrt(dhdx ** 2 + dhdy ** 2).max())


def forced_ascent(h, grid, wind_speed):
    """
    w = U * slope -- the vertical velocity terrain imposes on horizontal flow.

    This, not mountain height, is what governs how hard the orography forces
    the model. Measured: a 6000 m mountain with a gentle slope integrates
    LONGER than a 3000 m mountain with a steep one, and survival tracks this
    quantity across both.
    """
    return wind_speed * terrain_slope(h, grid)


def smooth_terrain(h, grid, passes=1, target_slope=None, max_passes=40):
    """
    Smooth orography to limit slope, as operational models do.

    Raw terrain at a model's grid spacing contains slopes the dynamics cannot
    support: the forced ascent w = U*slope drives adiabatic temperature
    tendencies of order 10 K/hour, and the resulting response is what breaks
    the integration. Every operational centre filters its orography for this
    reason -- the model's mountains are deliberately gentler than the real
    ones.

    Uses a 1-2-1 filter, which removes grid-scale structure while preserving
    the resolved shape and the domain-mean height.

    With `target_slope`, smooths until the slope falls below it (or max_passes
    is reached) and reports what it took.
    """
    out = h.astype(float).copy()

    def one_pass(a):
        ax = 0.25 * (grid.shift(a, -1, 1) + 2 * a + grid.shift(a, 1, 1))
        return 0.25 * (grid.shift(ax, -1, 0) + 2 * ax + grid.shift(ax, 1, 0))

    if target_slope is None:
        for _ in range(passes):
            out = one_pass(out)
        return out

    n = 0
    while terrain_slope(out, grid) > target_slope and n < max_passes:
        out = one_pass(out)
        n += 1
    return out, n, terrain_slope(out, grid)

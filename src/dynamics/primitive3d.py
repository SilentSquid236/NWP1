"""
Dry hydrostatic primitive equations on pressure levels.

Prognostic variables: u, v, theta (potential temperature).
Diagnosed each step: geopotential Phi (hydrostatic), omega (continuity).

    du/dt = -(u du/dx + v du/dy + omega du/dp) + f v - dPhi/dx
    dv/dt = -(u dv/dx + v dv/dy + omega dv/dp) - f u - dPhi/dy
    dth/dt = -(u dth/dx + v dth/dy + omega dth/dp)
    dPhi/dp = -R T / p                                  (hydrostatic)
    du/dx + dv/dy + domega/dp = 0                       (continuity)

DRY means theta is materially conserved -- no condensation, no latent heat, no
radiation. The thermodynamic equation has no source term at all. That is the
point: every term here is transport and pressure work, and a bug in any of
them is visible against an analytic balance. Adding moisture to a core whose
dry dynamics are unverified makes it impossible to tell a physics error from a
numerical one.

The hydrostatic approximation replaces the vertical momentum equation with
hydrostatic balance. Valid when horizontal scales greatly exceed vertical
ones -- true at 20 km grid spacing, marginal at 3 km, and wrong for
convection. At 3 km a nonhydrostatic core is the correct choice; this one is
honest at coarser resolution.
"""

import numpy as np

from grid import CGrid
from vertical import (PressureLevels, hydrostatic_geopotential, diagnose_omega,
                      T_from_theta, RD, CP, KAPPA, G0)
from subgrid import (hyperdiffusion, recommended_hyper_coeff,
                     hyper_stability_dt, divergence_damping)


class Primitive3D:
    def __init__(self, grid: CGrid, levels: PressureLevels, nu=0.0,
                 hyper=None, stochastic=None, div_damp=None,
                 sponge_levels=4, sponge_strength=8.0):
        """
        nu         : harmonic diffusion (m^2/s), usually 0 -- prefer hyper.
        hyper      : biharmonic coefficient (m^4/s). None = auto (damps the
                     2dx wave in ~3 h). 0 disables it, which is appropriate
                     ONLY for idealised tests -- a real run without a
                     grid-scale sink accumulates noise until it is useless.
        stochastic : a StochasticPerturbation, or None for a deterministic run.
        """
        self.grid = grid
        self.lev = levels
        self.nu = float(nu)
        self.hyper = (recommended_hyper_coeff(grid) if hyper is None
                      else float(hyper))
        self.stochastic = stochastic

        # DIVERGENCE DAMPING, as a post-step filter.
        #
        # Balancing the initial state removes the divergence present at t=0,
        # but integration keeps generating more -- from boundary relaxation,
        # from imperfectly balanced driving frames, and from the flow's own
        # adjustment. With no sink that gravity-wave energy accumulates and
        # the run dies.
        #
        # This is DIMENSIONLESS and applied after the RK3 update, not as a
        # tendency. Written as a tendency with a coefficient in m^2/s it must
        # satisfy nu*dt/dx^2 <= 0.25, and a coefficient chosen without knowing
        # dt violates that and blows the model up -- which is exactly what the
        # first version of this did. As a post-step filter the increment is
        # beta * dx^2 * grad(div), where beta IS the stability number, so it
        # is stable by construction for any dt.
        # DEFAULT OFF. Measured: any value >= 0.01 suppresses baroclinic
        # growth (eddy energy 1.21x/day -> 0.33x/day at 0.01), i.e. it damps
        # the weather along with the waves. It also did NOT rescue a run from
        # real analysis data -- see docs/STABILITY.md. Available for
        # experiments; not a fix.
        self.div_damp = 0.0 if div_damp is None else float(div_damp)
        if not 0.0 <= self.div_damp <= 0.25:
            raise ValueError(f"div_damp is a stability number and must be in "
                             f"[0, 0.25]; got {self.div_damp}")

        # SPONGE LAYER. The rigid lid reflects vertically propagating gravity
        # waves back into the domain, where they interfere and grow.
        #
        # The sponge amplifies DIVERGENCE damping near the top rather than
        # damping the wind itself. Gravity waves are divergent; balanced
        # large-scale flow is rotational. Relaxing the full wind toward its
        # horizontal mean would absorb the waves but also flatten a jet --
        # legitimate structure -- which is exactly what the thermal-wind test
        # caught when this was first written that way.
        self.sponge_levels = int(sponge_levels)
        self.sponge_strength = float(sponge_strength)
        self._sponge = np.ones((levels.nz, 1, 1))
        for k in range(self.sponge_levels):
            frac = (self.sponge_levels - k) / self.sponge_levels
            self._sponge[levels.nz - 1 - k, 0, 0] = \
                1.0 + self.sponge_strength * 0.5 * (1 - np.cos(np.pi * frac))

        shape = (levels.nz, grid.ny, grid.nx)
        self.u = np.zeros(shape)
        self.v = np.zeros(shape)
        self.theta = np.zeros(shape)

        self.time = 0.0
        self.step_count = 0

    # --- diagnostics -------------------------------------------------------

    def geopotential(self, theta=None):
        return hydrostatic_geopotential(
            self.theta if theta is None else theta, self.lev)

    def divergence(self, u, v):
        gr = self.grid
        return (gr.dx_forward(u) + gr.dy_forward(v))

    def omega(self, u=None, v=None):
        u = self.u if u is None else u
        v = self.v if v is None else v
        return diagnose_omega(self.divergence(u, v), self.lev)

    def temperature(self):
        return T_from_theta(self.theta, self.lev.p.reshape(-1, 1, 1))

    # --- operators ---------------------------------------------------------

    def _ddp(self, a):
        """d/dp with non-uniform level spacing (second-order interior)."""
        return np.gradient(a, self.lev.p, axis=0)

    def _horiz_adv(self, a, u_at_a, v_at_a):
        gr = self.grid
        dadx = 0.5 * (gr.dx_forward(a) + gr.dx_backward(a))
        dady = 0.5 * (gr.dy_forward(a) + gr.dy_backward(a))
        return u_at_a * dadx + v_at_a * dady

    def tendencies(self, u, v, theta):
        gr = self.grid

        phi = hydrostatic_geopotential(theta, self.lev)
        om = diagnose_omega(self.divergence(u, v), self.lev)

        # Interpolate cross-terms onto each variable's own grid point.
        v_at_u = gr.v_to_u(v)
        u_at_v = gr.u_to_v(u)
        om_at_u = gr.h_to_u(om)
        om_at_v = gr.h_to_v(om)

        du = -(self._horiz_adv(u, u, v_at_u) + om_at_u * self._ddp(u)) \
             + gr.f_u * v_at_u - gr.dx_backward(phi)

        dv = -(self._horiz_adv(v, u_at_v, v) + om_at_v * self._ddp(v)) \
             - gr.f_v * u_at_v - gr.dy_backward(phi)

        # theta lives at cell centres; average u and v onto centres.
        #
        # NOTE: this is ADVECTIVE form. Flux form would conserve the domain
        # integral of theta exactly, but only if the discrete continuity
        # relation div(v) + domega/dp = 0 held exactly -- and it does not,
        # because diagnose_omega() applies a linear correction to force omega
        # to zero at the lid and ground. Multiplying that residual by theta
        # (~300 K) is a large spurious heating, and flux form is unstable
        # here as a result. Doing it properly means constructing omega from
        # the same discrete operators used in the flux divergence, so that
        # continuity is satisfied to machine precision. That is the correct
        # fix and is not yet implemented; see README.
        u_at_h = 0.5 * (u + gr.shift(u, 1, 1))
        v_at_h = 0.5 * (v + gr.shift(v, 1, 0))
        dth = -(self._horiz_adv(theta, u_at_h, v_at_h) + om * self._ddp(theta))

        if self.nu > 0:
            du = du + self.nu * self._laplacian(u)
            dv = dv + self.nu * self._laplacian(v)
            dth = dth + self.nu * self._laplacian(theta)

        if self.hyper > 0:
            gr = self.grid
            du = du + hyperdiffusion(u, gr, self.hyper)
            dv = dv + hyperdiffusion(v, gr, self.hyper)
            # theta is damped about its level mean, so the dissipation cannot
            # erode the background stratification -- only the perturbations.
            th_ref = theta.mean(axis=(1, 2)).reshape(-1, 1, 1)
            dth = dth + hyperdiffusion(theta - th_ref, gr, self.hyper)

        if self.stochastic is not None:
            du = self.stochastic.apply(du)
            dv = self.stochastic.apply(dv)
            dth = self.stochastic.apply(dth)

        return du, dv, dth

    def _laplacian(self, a):
        gr = self.grid
        return ((gr.shift(a, 1, 1) - 2 * a + gr.shift(a, -1, 1)) / gr.dx**2 +
                (gr.shift(a, 1, 0) - 2 * a + gr.shift(a, -1, 0)) / gr.dy**2)

    # --- time stepping -----------------------------------------------------

    def max_dt(self, safety=0.6, wave_speed=None):
        """
        CFL for the hydrostatic system.

        With a rigid lid the external (Lamb) mode is suppressed, so the fastest
        signal is the gravest internal gravity wave -- of order 100 m/s for a
        troposphere-deep mode. We take that plus the wind speed. The vertical
        CFL on omega is checked too and usually is not limiting.
        """
        gr = self.grid
        c = 100.0 if wave_speed is None else wave_speed
        speed = c + max(np.abs(self.u).max(), np.abs(self.v).max(), 1e-9)
        dt_h = safety * min(gr.dx, gr.dy) / (speed * np.sqrt(2.0))

        om = np.abs(self.omega()).max()
        dt_v = np.inf
        if om > 0:
            dt_v = safety * self.lev.dp.min() / om

        # Explicit biharmonic diffusion has its own stability limit.
        dt_hyper = hyper_stability_dt(gr, self.hyper)

        return float(min(dt_h, dt_v, dt_hyper))

    def apply_divergence_filter(self):
        """
        Post-step divergence damping. Absorbs gravity waves; leaves the
        rotational (weather-carrying) flow alone. The sponge scales it up near
        the lid, where a rigid top would otherwise reflect waves back down.
        """
        if self.div_damp <= 0:
            return
        gr = self.grid
        beta = np.minimum(self.div_damp * self._sponge, 0.25)
        div = gr.dx_forward(self.u) + gr.dy_forward(self.v)
        self.u += beta * gr.dx ** 2 * gr.dx_backward(div)
        self.v += beta * gr.dy ** 2 * gr.dy_backward(div)

    def step(self, dt):
        if self.stochastic is not None:
            self.stochastic.advance(dt)
        u0, v0, t0 = self.u, self.v, self.theta

        du, dv, dth = self.tendencies(u0, v0, t0)
        u1, v1, t1 = u0 + dt / 3 * du, v0 + dt / 3 * dv, t0 + dt / 3 * dth

        du, dv, dth = self.tendencies(u1, v1, t1)
        u2, v2, t2 = u0 + dt / 2 * du, v0 + dt / 2 * dv, t0 + dt / 2 * dth

        du, dv, dth = self.tendencies(u2, v2, t2)
        self.u = u0 + dt * du
        self.v = v0 + dt * dv
        self.theta = t0 + dt * dth

        self.apply_divergence_filter()

        self.time += dt
        self.step_count += 1

    def run(self, duration, dt=None, callback=None, every=0):
        dt = dt or self.max_dt()
        n = int(np.ceil(duration / dt))
        dt = duration / n
        for k in range(n):
            self.step(dt)
            if callback and every and (k + 1) % every == 0:
                callback(self)
        return n

    # --- integral quantities ----------------------------------------------

    def total_theta(self):
        """Mass-weighted potential temperature. Conserved by dry advection."""
        w = np.empty(self.lev.nz)
        w[1:-1] = 0.5 * (self.lev.p[:-2] - self.lev.p[2:])
        w[0] = 0.5 * (self.lev.p[0] - self.lev.p[1])
        w[-1] = 0.5 * (self.lev.p[-2] - self.lev.p[-1])
        return float((self.theta * w.reshape(-1, 1, 1)).sum())

    def total_energy(self):
        """Kinetic + available potential (cp*T), mass weighted."""
        w = np.gradient(-self.lev.p) / G0
        T = self.temperature()
        ke = 0.5 * (self.u**2 + self.v**2)
        return float(((ke + CP * T) * w.reshape(-1, 1, 1)).sum())

    def diagnostics(self):
        om = self.omega()
        return {
            "time_h": self.time / 3600,
            "steps": self.step_count,
            "theta_int": self.total_theta(),
            "energy": self.total_energy(),
            "max|u|": float(np.abs(self.u).max()),
            "max|omega| Pa/s": float(np.abs(om).max()),
            "theta_min": float(self.theta.min()),
            "theta_max": float(self.theta.max()),
        }

    def __repr__(self):
        return (f"Primitive3D({self.lev.nz}x{self.grid.ny}x{self.grid.nx}, "
                f"hyper={self.hyper:.2e}, divdamp={self.div_damp:.2f}, "
                f"sponge={self.sponge_levels}L, "
                f"stoch={'on' if self.stochastic else 'off'}, "
                f"t={self.time/3600:.2f}h)")

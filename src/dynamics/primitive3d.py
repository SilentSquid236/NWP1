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
from subgrid import hyperdiffusion, recommended_hyper_coeff, hyper_stability_dt


class Primitive3D:
    def __init__(self, grid: CGrid, levels: PressureLevels, nu=0.0,
                 hyper=None, stochastic=None):
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
                f"hyper={self.hyper:.2e}, "
                f"stoch={'on' if self.stochastic else 'off'}, "
                f"t={self.time/3600:.2f}h)")

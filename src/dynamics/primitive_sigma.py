"""
Dry hydrostatic primitive equations on terrain-following sigma levels.

Prognostic: u, v, theta, and pi = p_s - p_top (the column's pressure depth).
Diagnosed each step: Phi (hydrostatic), sigma_dot (continuity).

    dpi/dt = -integral_0^1 div(pi V) dsigma
    du/dt  = -adv_h(u) - sigma_dot du/dsigma + f v + F_x
    dv/dt  = -adv_h(v) - sigma_dot dv/dsigma - f u + F_y
    dth/dt = -adv_h(th) - sigma_dot dth/dsigma
    F      = -grad_sigma(Phi) - (R T / p) sigma grad(pi)

WHAT CHANGED FROM THE PRESSURE-COORDINATE VERSION

Surface pressure is now PROGNOSTIC. That is the whole point: in pressure
coordinates omega had to vanish at both ends of a column that could not move,
which over-constrained the system and required a correction whose residual fed
a divergence/vertical-velocity feedback. Measured consequence: divergence
within 2-3 hours on real analyses at every damping setting tried.

Here sigma_dot = 0 at lid and ground falls out of the formulation, verified
to 0.00e+00 in test_sigma.py, and the column exchanges mass through a moving
surface instead.

Terrain comes with the coordinate rather than being a separate feature.
"""

import numpy as np

from grid import CGrid
from sigma import (SigmaLevels, hydrostatic_geopotential, continuity,
                   vertical_advection, pressure_gradient_force,
                   RD, CP, KAPPA, P0, G0)
from subgrid import hyperdiffusion, recommended_hyper_coeff, hyper_stability_dt
from turbulence import vertical_mixing, richardson, mixing_stability_dt


class PrimitiveSigma:
    def __init__(self, grid: CGrid, levels: SigmaLevels, terrain=None,
                 hyper=None, stochastic=None, ref_pgf=True,
                 sponge_levels=5, sponge_rate=1.0 / 900.0, mixing=True):
        self.grid = grid
        self.lev = levels

        shape = (levels.nz, grid.ny, grid.nx)
        self.u = np.zeros(shape)
        self.v = np.zeros(shape)
        self.theta = np.zeros(shape)
        self.pi = np.full((grid.ny, grid.nx), 95_000.0)

        # Surface geopotential. Zero terrain is flat ground at sea level.
        self.terrain = (np.zeros((grid.ny, grid.nx)) if terrain is None
                        else np.asarray(terrain, dtype=float))
        self.phi_s = G0 * self.terrain

        self.hyper = (recommended_hyper_coeff(grid) if hyper is None
                      else float(hyper))
        self.stochastic = stochastic

        # Reference-state pressure-gradient force. The plain form is stable on
        # flat ground but not over terrain -- see docs/RESEARCH_LOG.md.
        self.ref_pgf = bool(ref_pgf)

        # ABSORBING LAYER at the model top.
        #
        # Flow over terrain excites vertically propagating gravity waves. The
        # rigid lid reflects them straight back down, where they interfere
        # with the upgoing waves and grow. Cross-sections show exactly this:
        # over flat ground the vertical velocity is horizontally uniform and
        # harmless, while over a mountain cellular structure appears on the
        # flanks and extends upward until the run dies.
        #
        # The sponge relaxes the wind toward a REFERENCE state rather than
        # toward the horizontal mean. Relaxing to the mean also absorbs the
        # waves, but it flattens a jet -- legitimate structure -- which the
        # thermal-wind test caught when the pressure-coordinate version was
        # written that way.
        self.sponge_levels = int(sponge_levels)
        self.sponge_rate = float(sponge_rate)
        self._sponge = np.zeros((levels.nz, 1, 1))
        for k in range(self.sponge_levels):
            frac = (self.sponge_levels - k) / self.sponge_levels
            self._sponge[k, 0, 0] = self.sponge_rate * 0.5 * (
                1 - np.cos(np.pi * frac))
        self._u_ref = None
        self._v_ref = None

        # Richardson-number vertical mixing. The model develops real shear
        # instability and without this has nothing to dissipate it.
        self.mixing = bool(mixing)
        self._K_last = None

        self.time = 0.0
        self.step_count = 0

    # --- diagnostics -------------------------------------------------------

    @property
    def surface_pressure(self):
        return self.lev.p_top + self.pi

    def pressure(self):
        return self.lev.pressure(self.pi)

    def temperature(self):
        return self.theta * (self.pressure() / P0) ** KAPPA

    def geopotential(self, theta=None, pi=None):
        return hydrostatic_geopotential(
            self.theta if theta is None else theta,
            self.pi if pi is None else pi,
            self.lev, phi_surface=self.phi_s)

    def sigma_dot(self, u=None, v=None, pi=None):
        _, sd = continuity(self.u if u is None else u,
                           self.v if v is None else v,
                           self.pi if pi is None else pi,
                           self.lev, self.grid)
        return sd

    # --- operators ---------------------------------------------------------

    def _horiz_adv(self, a, u_at_a, v_at_a):
        gr = self.grid
        dadx = 0.5 * (gr.dx_forward(a) + gr.dx_backward(a))
        dady = 0.5 * (gr.dy_forward(a) + gr.dy_backward(a))
        return u_at_a * dadx + v_at_a * dady

    def _laplacian(self, a):
        gr = self.grid
        return ((gr.shift(a, 1, 1) - 2 * a + gr.shift(a, -1, 1)) / gr.dx**2 +
                (gr.shift(a, 1, 0) - 2 * a + gr.shift(a, -1, 0)) / gr.dy**2)

    def set_reference(self, u=None, v=None):
        """
        Freeze the state the sponge relaxes toward. Called automatically on the
        first step if not set explicitly.
        """
        self._u_ref = (self.u if u is None else u).copy()
        self._v_ref = (self.v if v is None else v).copy()

    def tendencies(self, u, v, theta, pi):
        gr, lev = self.grid, self.lev

        phi = hydrostatic_geopotential(theta, pi, lev, phi_surface=self.phi_s)
        dpi_dt, sd = continuity(u, v, pi, lev, gr)
        # Reference profile: the horizontal-mean temperature on each sigma
        # surface. Recomputed each call so it tracks the evolving state.
        T_ref = (theta * (lev.pressure(pi) / P0) ** KAPPA).mean(axis=(1, 2))
        fx, fy = pressure_gradient_force(phi, theta, pi, lev, gr,
                                         reference=T_ref if self.ref_pgf else None)

        v_at_u = gr.v_to_u(v)
        u_at_v = gr.u_to_v(u)

        du = (-self._horiz_adv(u, u, v_at_u)
              - vertical_advection(u, sd, lev)
              + gr.f_u * v_at_u + fx)

        dv = (-self._horiz_adv(v, u_at_v, v)
              - vertical_advection(v, sd, lev)
              - gr.f_v * u_at_v + fy)

        u_at_h = 0.5 * (u + gr.shift(u, 1, 1))
        v_at_h = 0.5 * (v + gr.shift(v, 1, 0))
        dth = (-self._horiz_adv(theta, u_at_h, v_at_h)
               - vertical_advection(theta, sd, lev))

        if self.hyper > 0:
            du = du + hyperdiffusion(u, gr, self.hyper)
            dv = dv + hyperdiffusion(v, gr, self.hyper)
            th_ref = theta.mean(axis=(1, 2), keepdims=True)
            dth = dth + hyperdiffusion(theta - th_ref, gr, self.hyper)
            # NOTE: no diffusion on pi. Hyperdiffusion is only conservative on
            # a periodic domain; applied to the prognostic surface pressure on
            # a bounded domain it acts as a MASS SOURCE. Measured: p_s
            # inflating from 1088 to 1243 hPa over three hours. Grid-scale
            # noise in pi has to be controlled by the wind field that
            # generates it, not by diffusing mass.

        if self.mixing:
            mu, mv, mth, K = vertical_mixing(u, v, theta, pi, lev)
            du = du + mu
            dv = dv + mv
            dth = dth + mth
            self._K_last = K

        if self.sponge_levels > 0 and self._u_ref is not None:
            du = du - self._sponge * (u - self._u_ref)
            dv = dv - self._sponge * (v - self._v_ref)

        if self.stochastic is not None:
            du = self.stochastic.apply(du)
            dv = self.stochastic.apply(dv)
            dth = self.stochastic.apply(dth)

        return du, dv, dth, dpi_dt

    # --- time stepping -----------------------------------------------------

    def max_dt(self, safety=0.6, wave_speed=None):
        """
        CFL limit.

        With PROGNOSTIC surface pressure the fastest signal is the EXTERNAL
        (Lamb) wave at sqrt(R*T) ~ 300 m/s, not the internal gravity wave.
        Pressure coordinates with a rigid lid suppress that mode, so the
        pressure-coordinate core could use ~100 m/s -- carrying that number
        over here makes dt about 3x too large, and the model blows up with
        surface pressure going NEGATIVE within ~20 steps.

        This is the cost of a prognostic free surface, and it is why
        operational models split the fast external mode off into a
        sub-stepped or semi-implicit solver rather than resolving it
        explicitly.
        """
        gr = self.grid
        if wave_speed is None:
            T = np.clip(self.temperature(), 150.0, 350.0)
            wave_speed = float(np.sqrt(RD * T.max()))
        speed = wave_speed + max(np.abs(self.u).max(), np.abs(self.v).max(), 1e-9)
        dt_h = safety * min(gr.dx, gr.dy) / (speed * np.sqrt(2.0))

        sd = np.abs(self.sigma_dot()).max()
        dt_v = np.inf
        if sd > 0:
            dt_v = safety * self.lev.dsigma.min() / sd

        return float(min(dt_h, dt_v, hyper_stability_dt(gr, self.hyper)))

    def step(self, dt):
        if self.sponge_levels > 0 and self._u_ref is None:
            self.set_reference()
        if self.stochastic is not None:
            self.stochastic.advance(dt)

        u0, v0, t0, p0 = self.u, self.v, self.theta, self.pi

        du, dv, dth, dp = self.tendencies(u0, v0, t0, p0)
        u1, v1, t1, p1 = (u0 + dt / 3 * du, v0 + dt / 3 * dv,
                          t0 + dt / 3 * dth, p0 + dt / 3 * dp)

        du, dv, dth, dp = self.tendencies(u1, v1, t1, p1)
        u2, v2, t2, p2 = (u0 + dt / 2 * du, v0 + dt / 2 * dv,
                          t0 + dt / 2 * dth, p0 + dt / 2 * dp)

        du, dv, dth, dp = self.tendencies(u2, v2, t2, p2)
        self.u = u0 + dt * du
        self.v = v0 + dt * dv
        self.theta = t0 + dt * dth
        self.pi = p0 + dt * dp

        self.time += dt
        self.step_count += 1

    def run(self, duration, dt=None, callback=None, every=0, adaptive=True,
            recheck_steps=50):
        """
        Integrate forward.

        ADAPTIVE TIMESTEP. The CFL limit depends on the wind, which grows
        during a run -- a dt chosen from the initial state can be violated an
        hour later, and the failure looks like a physics instability rather
        than a stale timestep. Re-checking periodically and subdividing costs
        one diagnostic evaluation per check.
        """
        dt = dt or self.max_dt()
        elapsed = 0.0
        k = 0
        while elapsed < duration - 1e-9:
            if adaptive and k % recheck_steps == 0:
                limit = self.max_dt()
                if dt > limit:
                    dt = limit
            step_dt = min(dt, duration - elapsed)
            self.step(step_dt)
            elapsed += step_dt
            k += 1
            if callback and every and k % every == 0:
                callback(self)
            if not np.isfinite(self.pi).all():
                break
        return k

    # --- integrals ---------------------------------------------------------

    def total_mass(self):
        """Column mass ~ integral of pi over the domain."""
        return float(self.pi.sum() * self.grid.dx * self.grid.dy / G0)

    def total_theta(self):
        ds = self.lev.dsigma.reshape(-1, 1, 1)
        return float((self.theta * self.pi * ds).sum())

    @property
    def lamb_wave_speed(self):
        """External gravity wave speed, sqrt(R T) -- sets the timestep."""
        return float(np.sqrt(RD * np.clip(self.temperature(), 150, 350).max()))

    def diagnostics(self):
        sd = self.sigma_dot()
        return {
            "time_h": self.time / 3600,
            "steps": self.step_count,
            "mass": self.total_mass(),
            "theta_int": self.total_theta(),
            "max|u|": float(np.abs(self.u).max()),
            "max|sigma_dot|": float(np.abs(sd).max()),
            "p_s_min_hPa": float(self.surface_pressure.min() / 100),
            "p_s_max_hPa": float(self.surface_pressure.max() / 100),
            "theta_min": float(self.theta.min()),
            "theta_max": float(self.theta.max()),
        }

    def __repr__(self):
        return (f"PrimitiveSigma({self.lev.nz}x{self.grid.ny}x{self.grid.nx}, "
                f"terrain {self.terrain.min():.0f}-{self.terrain.max():.0f} m, "
                f"t={self.time/3600:.2f}h)")

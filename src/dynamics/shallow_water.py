"""
Shallow-water dynamical core on an Arakawa C-grid, beta-plane, periodic.

Equations (advective form, single layer of depth h):

    du/dt = -u du/dx - v du/dy + f v - g dh/dx
    dv/dt = -u dv/dx - v dv/dy - f u - g dh/dy
    dh/dt = -d(h u)/dx - d(h v)/dy

The continuity equation is written in FLUX form so total mass is conserved to
machine precision -- it becomes a telescoping sum over the periodic domain.
Momentum is in advective form, which is simpler and adequate here; a flux-form
or vector-invariant momentum equation conserves energy better and is the right
upgrade before moving to 3D.

Time integration is the three-stage Runge-Kutta of Wicker & Skamarock (2002),
the scheme WRF uses. It is stable for the CFL numbers we need and damps the
computational mode that plagues leapfrog.
"""

import numpy as np

from grid import CGrid

G = 9.80665      # m s^-2


class ShallowWaterModel:
    def __init__(self, grid: CGrid, H=10_000.0, nu=0.0, g=G):
        """
        grid : CGrid
        H    : mean layer depth (m). Gravity wave speed is sqrt(gH).
        nu   : optional harmonic diffusion (m^2/s), 0 = off.
        """
        self.grid = grid
        self.H = float(H)
        self.nu = float(nu)
        self.g = float(g)

        self.u = np.zeros((grid.ny, grid.nx))
        self.v = np.zeros((grid.ny, grid.nx))
        self.h = np.full((grid.ny, grid.nx), self.H)

        self.time = 0.0
        self.step_count = 0

    # --- physics -----------------------------------------------------------

    @property
    def gravity_wave_speed(self):
        return np.sqrt(self.g * self.H)

    def max_dt(self, safety=0.8):
        """
        CFL limit. Explicit schemes cannot step faster than the fastest signal
        crossing a cell -- here the external gravity wave, ~313 m/s for
        H = 10 km. This is why explicit models take small steps.
        """
        gr = self.grid
        speed = self.gravity_wave_speed + max(np.abs(self.u).max(),
                                              np.abs(self.v).max())
        return safety * min(gr.dx, gr.dy) / (speed * np.sqrt(2.0))

    def tendencies(self, u, v, h):
        gr = self.grid
        dx, dy = gr.dx, gr.dy

        # --- continuity, flux form ----------------------------------------
        # Mass flux through each face uses the depth interpolated to that face.
        U = gr.h_to_u(h) * u                       # at u points
        V = gr.h_to_v(h) * v                       # at v points
        dh = -(gr.dx_forward(U, dx) + gr.dy_forward(V, dy))

        # --- u momentum ----------------------------------------------------
        v_at_u = gr.v_to_u(v)
        dudx = 0.5 * (gr.dx_forward(u, dx) + gr.dx_backward(u, dx))   # centred
        dudy = 0.5 * (gr.dy_forward(u, dy) + gr.dy_backward(u, dy))
        du = (-u * dudx - v_at_u * dudy
              + gr.f_u * v_at_u
              - self.g * gr.dx_backward(h, dx))

        # --- v momentum ----------------------------------------------------
        u_at_v = gr.u_to_v(u)
        dvdx = 0.5 * (gr.dx_forward(v, dx) + gr.dx_backward(v, dx))
        dvdy = 0.5 * (gr.dy_forward(v, dy) + gr.dy_backward(v, dy))
        dv = (-u_at_v * dvdx - v * dvdy
              - gr.f_v * u_at_v
              - self.g * gr.dy_backward(h, dy))

        if self.nu > 0:
            du += self.nu * self._laplacian(u)
            dv += self.nu * self._laplacian(v)

        return du, dv, dh

    def _laplacian(self, a):
        gr = self.grid
        return ((np.roll(a, -1, axis=1) - 2 * a + np.roll(a, 1, axis=1)) / gr.dx**2 +
                (np.roll(a, -1, axis=0) - 2 * a + np.roll(a, 1, axis=0)) / gr.dy**2)

    # --- time stepping -----------------------------------------------------

    def step(self, dt):
        """One RK3 step (Wicker & Skamarock 2002)."""
        u0, v0, h0 = self.u, self.v, self.h

        du, dv, dh = self.tendencies(u0, v0, h0)
        u1, v1, h1 = u0 + dt / 3 * du, v0 + dt / 3 * dv, h0 + dt / 3 * dh

        du, dv, dh = self.tendencies(u1, v1, h1)
        u2, v2, h2 = u0 + dt / 2 * du, v0 + dt / 2 * dv, h0 + dt / 2 * dh

        du, dv, dh = self.tendencies(u2, v2, h2)
        self.u, self.v, self.h = u0 + dt * du, v0 + dt * dv, h0 + dt * dh

        self.time += dt
        self.step_count += 1

    def run(self, duration, dt=None, callback=None, every=0):
        """Integrate forward. Returns the number of steps taken."""
        dt = dt or self.max_dt()
        n = int(np.ceil(duration / dt))
        dt = duration / n                        # land exactly on the target
        for k in range(n):
            self.step(dt)
            if callback and every and (k + 1) % every == 0:
                callback(self)
        return n

    # --- diagnostics -------------------------------------------------------

    def total_mass(self):
        return float(self.h.sum() * self.grid.dx * self.grid.dy)

    def total_energy(self):
        """Kinetic + available potential energy."""
        gr = self.grid
        h_u, h_v = gr.h_to_u(self.h), gr.h_to_v(self.h)
        ke = 0.5 * (h_u * self.u**2 + h_v * self.v**2).sum()
        pe = 0.5 * self.g * ((self.h - self.H)**2).sum()
        return float((ke + pe) * gr.dx * gr.dy)

    def potential_vorticity(self):
        """(zeta + f) / h, at cell corners."""
        gr = self.grid
        zeta = (gr.dx_forward(self.v, gr.dx) - gr.dy_forward(self.u, gr.dy))
        return (zeta + gr.f_h) / self.h

    def diagnostics(self):
        return {
            "time_h": self.time / 3600,
            "steps": self.step_count,
            "mass": self.total_mass(),
            "energy": self.total_energy(),
            "max|u|": float(np.abs(self.u).max()),
            "max|v|": float(np.abs(self.v).max()),
            "h_min": float(self.h.min()),
            "h_max": float(self.h.max()),
        }

    def __repr__(self):
        return (f"ShallowWaterModel({self.grid.nx}x{self.grid.ny}, H={self.H:.0f}m, "
                f"c={self.gravity_wave_speed:.0f}m/s, t={self.time/3600:.2f}h)")

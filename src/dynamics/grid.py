"""
Arakawa C-grid for a doubly-periodic beta-plane.

Variable placement (the whole point of the C-grid):

        v[j+1]
          ^
          |
  u[i] -> h[i,j] -> u[i+1]
          |
          ^
        v[j]

  h  cell centre        (ny, nx)
  u  western cell face  (ny, nx)   -- u[j,i] sits at x = (i - 1/2) dx
  v  southern cell face (ny, nx)   -- v[j,i] sits at y = (j - 1/2) dy

Why C-grid: the pressure gradient and divergence become adjacent differences
rather than 2*dx differences, which kills the checkerboard null mode that
plagues the A-grid. It is what WRF, MPAS, and most operational models use.

Indices are periodic, implemented with np.roll, so there are no boundary
special cases in the interior operators. Lateral boundary conditions for a
limited-area domain come later, as a relaxation zone applied after each step.
"""

import numpy as np


class CGrid:
    def __init__(self, nx, ny, dx, dy, f0=1.0e-4, beta=1.6e-11):
        self.nx, self.ny = int(nx), int(ny)
        self.dx, self.dy = float(dx), float(dy)
        self.f0, self.beta = float(f0), float(beta)

        self.Lx = self.nx * self.dx
        self.Ly = self.ny * self.dy

        # Cell-centre coordinates
        self.xc = (np.arange(self.nx) + 0.5) * self.dx
        self.yc = (np.arange(self.ny) + 0.5) * self.dy
        self.Xc, self.Yc = np.meshgrid(self.xc, self.yc)

        # Face coordinates
        self.xu = np.arange(self.nx) * self.dx            # u points
        self.yv = np.arange(self.ny) * self.dy            # v points

        # Coriolis on a beta-plane, centred on the domain midpoint.
        y0 = 0.5 * self.Ly
        self.f_h = self.f0 + self.beta * (self.Yc - y0)                 # centres
        self.f_u = self.f0 + self.beta * (self.Yc - y0)                 # same rows
        _, Yv = np.meshgrid(self.xc, self.yv)
        self.f_v = self.f0 + self.beta * (Yv - y0)

    # --- differencing (periodic) -------------------------------------------

    @staticmethod
    def dx_forward(a, dx):
        """(a[i+1] - a[i]) / dx  -- centre -> face+1 or face -> centre."""
        return (np.roll(a, -1, axis=1) - a) / dx

    @staticmethod
    def dx_backward(a, dx):
        """(a[i] - a[i-1]) / dx"""
        return (a - np.roll(a, 1, axis=1)) / dx

    @staticmethod
    def dy_forward(a, dy):
        return (np.roll(a, -1, axis=0) - a) / dy

    @staticmethod
    def dy_backward(a, dy):
        return (a - np.roll(a, 1, axis=0)) / dy

    # --- interpolation between staggered points ----------------------------

    @staticmethod
    def h_to_u(h):
        """Cell centre -> western face: average of the two straddling cells."""
        return 0.5 * (h + np.roll(h, 1, axis=1))

    @staticmethod
    def h_to_v(h):
        return 0.5 * (h + np.roll(h, 1, axis=0))

    @staticmethod
    def v_to_u(v):
        """Four-point average of v onto u points. Needed for Coriolis."""
        return 0.25 * (v + np.roll(v, -1, axis=0)
                       + np.roll(v, 1, axis=1)
                       + np.roll(np.roll(v, -1, axis=0), 1, axis=1))

    @staticmethod
    def u_to_v(u):
        return 0.25 * (u + np.roll(u, -1, axis=1)
                       + np.roll(u, 1, axis=0)
                       + np.roll(np.roll(u, -1, axis=1), 1, axis=0))

    def __repr__(self):
        return (f"CGrid({self.nx}x{self.ny}, dx={self.dx/1000:.0f}km, "
                f"dy={self.dy/1000:.0f}km, Lx={self.Lx/1000:.0f}km, "
                f"f0={self.f0:.2e}, beta={self.beta:.2e})")

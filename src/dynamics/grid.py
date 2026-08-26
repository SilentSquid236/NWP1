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

Interior operators are index-shifts; edge_mode decides whether those shifts
wrap (periodic) or replicate the edge value (limited-area). Lateral
boundary forcing is applied as a Davies relaxation zone -- see boundaries.py.
"""

import numpy as np


class CGrid:
    def __init__(self, nx, ny, dx, dy, f0=1.0e-4, beta=1.6e-11,
                 edge_mode="periodic"):
        """
        edge_mode : "periodic" wraps (idealised tests);
                    "replicate" repeats the edge value, so the domain has real
                    boundaries. A limited-area run uses "replicate" plus a
                    Davies relaxation zone driven by external data.
        """
        if edge_mode not in ("periodic", "replicate"):
            raise ValueError(f"unknown edge_mode {edge_mode!r}")
        self.edge_mode = edge_mode
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
        # Corner points (j-1/2, i-1/2) share the v rows in y.
        self.f_corner = self.f_v.copy()

    # --- shifting, honouring edge_mode -------------------------------------

    def shift(self, a, n, axis):
        """
        a[i+n] with the grid's edge treatment.

        periodic  : wraps around (np.roll)
        replicate : values beyond the edge repeat the edge value, i.e. a
                    zero-gradient (Neumann) condition. Combined with a
                    relaxation zone this stops the domain wrapping while
                    keeping the interior operators unchanged.
        """
        if self.edge_mode == "periodic":
            return np.roll(a, -n, axis=axis)

        out = np.roll(a, -n, axis=axis)
        if n > 0:      # pulled from beyond the far edge
            if axis == 1:
                out[:, -n:] = a[:, -1][:, None]
            else:
                out[-n:, :] = a[-1, :][None, :]
        elif n < 0:    # pulled from before the near edge
            k = -n
            if axis == 1:
                out[:, :k] = a[:, 0][:, None]
            else:
                out[:k, :] = a[0, :][None, :]
        return out

    # --- differencing ------------------------------------------------------

    def dx_forward(self, a, dx=None):
        """(a[i+1] - a[i]) / dx"""
        return (self.shift(a, 1, 1) - a) / (dx or self.dx)

    def dx_backward(self, a, dx=None):
        """(a[i] - a[i-1]) / dx"""
        return (a - self.shift(a, -1, 1)) / (dx or self.dx)

    def dy_forward(self, a, dy=None):
        return (self.shift(a, 1, 0) - a) / (dy or self.dy)

    def dy_backward(self, a, dy=None):
        return (a - self.shift(a, -1, 0)) / (dy or self.dy)

    # --- interpolation between staggered points ----------------------------

    def h_to_u(self, h):
        """Cell centre -> western face."""
        return 0.5 * (h + self.shift(h, -1, 1))

    def h_to_v(self, h):
        return 0.5 * (h + self.shift(h, -1, 0))

    def v_to_u(self, v):
        """Four-point average of v onto u points. Needed for Coriolis."""
        return 0.25 * (v + self.shift(v, 1, 0)
                       + self.shift(v, -1, 1)
                       + self.shift(self.shift(v, 1, 0), -1, 1))

    def u_to_v(self, u):
        return 0.25 * (u + self.shift(u, 1, 1)
                       + self.shift(u, -1, 0)
                       + self.shift(self.shift(u, 1, 1), -1, 0))

    def __repr__(self):
        return (f"CGrid({self.nx}x{self.ny}, dx={self.dx/1000:.0f}km, "
                f"dy={self.dy/1000:.0f}km, Lx={self.Lx/1000:.0f}km, "
                f"f0={self.f0:.2e}, beta={self.beta:.2e}, "
                f"edges={self.edge_mode})")

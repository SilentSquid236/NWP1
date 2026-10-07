"""
A 2-D (x-z) compressible non-hydrostatic core: CAM roadmap stage S1a.

WHY THIS EXISTS

The roadmap (docs/CAM_DESIGN.md) replaces the hydrostatic sigma core with a
split-explicit non-hydrostatic one. Stage S1 builds the time-splitting and the
vertically implicit acoustic solver first, in two dimensions, and holds them
to published benchmarks before any of it touches the forecast model.

S1 is done in two steps:
- S1a (this file): height coordinate, flat ground, periodic in x. This is the
  configuration of the classic split-explicit benchmark papers, so the
  time-splitting and the acoustic solver are checked against them alone.
- S1b: the same solver in the hydrostatic-pressure (mass) coordinate of
  Laprise (1992), with terrain, checked against a mountain wave. The
  vertically implicit w-pressure coupling below carries over unchanged in
  structure.

The split into S1a and S1b is a deviation from the roadmap as first written,
recorded in the research log. A failure in a first benchmark could otherwise
come from the solver or from the coordinate; this way it can only come from
one of them.

EQUATIONS (dry; Exner-pressure form, as in Klemp and Wilhelmson 1978 and
Skamarock and Klemp 1992)

    du/dt   = -(u du/dx + w du/dz) - cp theta dpi'/dx                     + D(u)
    dw/dt   = -(u dw/dx + w dw/dz) - cp theta dpi'/dz + g theta'/theta_bar + D(w)
    dth'/dt = -(u dth'/dx + w dth'/dz) - w dtheta_bar/dz                  + D(theta')
    dpi'/dt = -(cbar^2 / (cp rho_bar theta_bar^2)) div(rho_bar theta_bar V)

- theta = theta_bar(z) + theta' is the potential temperature, and
  pi = pi_bar(z) + pi' is the Exner pressure. The reference state is
  hydrostatic: d pi_bar/dz = -g / (cp theta_bar).
- cbar^2 = (cp/cv) R pi_bar theta_bar is the squared speed of sound.
- D is a constant-K second-order diffusion (the density-current benchmark
  prescribes K = 75 m^2/s).

GRID
- Arakawa C: u on x faces (u[:, i] is the face west of cell i), w on z faces
  (w[k] at z = k dz), theta' and pi' at cell centres.
- Periodic in x. Rigid lid and floor (w = 0).

TIME INTEGRATION (Wicker and Skamarock 2002; Klemp et al. 2007)
- Three-stage Runge-Kutta for the slow terms (advection, buoyancy,
  dtheta_bar/dz, diffusion), evaluated at the stage state.
- Inside each stage, ns/3, ns/2 and ns acoustic steps of dtau = dt/ns start
  again from the state at the beginning of the large step, with the slow
  forcing held fixed.
  - u is forward.
  - w and pi' are implicit in the vertical, with off-centring beta. That is
    a tridiagonal solve per column.
  - A pressure-based divergence damping (Skamarock and Klemp 1992) keeps the
    acoustic modes in check.
- Horizontal advection is third-order upwind-biased (the `upwind3` of the
  forecast model). Vertical advection is second-order centred.
"""
import numpy as np

G = 9.80665
RD = 287.04
CP = 1004.5
CV = CP - RD
P0 = 1.0e5
KAPPA = RD / CP


def reference_state(z_c, z_w, theta_c, p_surface=P0):
    """
    Hydrostatic Exner pressure on faces (integrated up from the floor with the
    centre thetas) and at centres (the face average), with the densities.
    """
    nz = len(z_c)
    pi_w = np.empty(nz + 1)
    pi_w[0] = (p_surface / P0) ** KAPPA
    for k in range(nz):
        pi_w[k + 1] = pi_w[k] - G * (z_w[k + 1] - z_w[k]) / (CP * theta_c[k])
    pi_c = 0.5 * (pi_w[1:] + pi_w[:-1])
    return pi_c, pi_w


def _thomas(a, b, c, d):
    """Solve tridiagonal systems along axis 0 (a: sub, b: diag, c: super)."""
    n = b.shape[0]
    cp_ = np.empty_like(b)
    dp_ = np.empty_like(d)
    cp_[0] = c[0] / b[0]
    dp_[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp_[i - 1]
        cp_[i] = c[i] / m
        dp_[i] = (d[i] - a[i] * dp_[i - 1]) / m
    x = np.empty_like(d)
    x[-1] = dp_[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp_[i] - cp_[i] * x[i + 1]
    return x


class NH2D:
    def __init__(self, nx, nz, dx, dz, theta_bar, K=0.0, ns=6, beta=0.1,
                 div_damp=0.1, rayleigh=None):
        self.nx, self.nz, self.dx, self.dz = nx, nz, float(dx), float(dz)
        self.z_c = (np.arange(nz) + 0.5) * dz
        self.z_w = np.arange(nz + 1) * dz
        self.tb_c = np.asarray(theta_bar(self.z_c), dtype=float) * np.ones(nz)
        self.tb_w = np.asarray(theta_bar(self.z_w), dtype=float) * np.ones(nz + 1)
        self.pi_c, self.pi_w = reference_state(self.z_c, self.z_w, self.tb_c)
        self.rho_c = P0 * self.pi_c ** (CV / RD) / (RD * self.tb_c)
        self.rho_w = P0 * self.pi_w ** (CV / RD) / (RD * self.tb_w)
        self.c2_c = (CP / CV) * RD * self.pi_c * self.tb_c
        self.dtb_dz_c = np.gradient(self.tb_c, self.z_c) if nz > 2 else np.zeros(nz)
        self.K = float(K)
        self.ns = int(ns)
        if self.ns % 6:
            raise ValueError("ns must be a multiple of 6 (the stages use ns/3, ns/2 and ns)")
        self.beta = float(beta)
        self.div_damp = float(div_damp)
        self.rayleigh = None if rayleigh is None else np.asarray(rayleigh, dtype=float)
        self.u = np.zeros((nz, nx))
        self.w = np.zeros((nz + 1, nx))
        self.th = np.zeros((nz, nx))
        self.pp = np.zeros((nz, nx))
        self.time = 0.0

    # --- operators -------------------------------------------------------
    @staticmethod
    def _xadv_u3(a, vel, dx):
        """vel * da/dx, third-order upwind-biased, periodic in x."""
        p1, m1 = np.roll(a, -1, -1), np.roll(a, 1, -1)
        p2, m2 = np.roll(a, -2, -1), np.roll(a, 2, -1)
        d4 = (-p2 + 8 * p1 - 8 * m1 + m2) / (12 * dx)
        diss = (p2 - 4 * p1 + 6 * a - 4 * m1 + m2) / (12 * dx)
        return vel * d4 + np.abs(vel) * diss

    def _lap(self, a, rigid=False):
        """Second-order Laplacian, zero-gradient at the floor and lid."""
        lx = (np.roll(a, -1, -1) - 2 * a + np.roll(a, 1, -1)) / self.dx ** 2
        up = np.concatenate([a[1:], a[-1:]], axis=0)
        dn = np.concatenate([a[:1], a[:-1]], axis=0)
        lz = (up - 2 * a + dn) / self.dz ** 2
        if rigid:
            lz[0] = 0.0
            lz[-1] = 0.0
        return lx + lz

    @staticmethod
    def _zgrad_c(a, dz):
        g = np.empty_like(a)
        g[1:-1] = (a[2:] - a[:-2]) / (2 * dz)
        g[0] = (a[1] - a[0]) / dz
        g[-1] = (a[-1] - a[-2]) / dz
        return g

    def slow(self, u, w, th):
        """Slow tendencies: advection, buoyancy, background theta, diffusion, damping."""
        dx, dz = self.dx, self.dz
        w_c = 0.5 * (w[1:] + w[:-1])
        u_c = 0.5 * (u + np.roll(u, -1, -1))
        w_u = 0.5 * (w_c + np.roll(w_c, 1, -1))
        u_w = np.empty_like(w)
        u_w[1:-1] = 0.5 * (u_c[1:] + u_c[:-1])
        u_w[0] = u_c[0]
        u_w[-1] = u_c[-1]
        Ru = -(self._xadv_u3(u, u, dx) + w_u * self._zgrad_c(u, dz))
        Rth = (-(self._xadv_u3(th, u_c, dx) + w_c * self._zgrad_c(th, dz))
               - w_c * self.dtb_dz_c[:, None])
        Rw = np.zeros_like(w)
        Rw[1:-1] = -(self._xadv_u3(w, u_w, dx)[1:-1]
                     + w[1:-1] * (w[2:] - w[:-2]) / (2 * dz))
        Rw[1:-1] += G * 0.5 * (th[1:] + th[:-1]) / self.tb_w[1:-1, None]
        if self.K > 0:
            Ru += self.K * self._lap(u)
            Rth += self.K * self._lap(th)
            Rw += self.K * self._lap(w, rigid=True)
        if self.rayleigh is not None:
            Rw -= self.rayleigh[:, None] * w
        Rw[0] = 0.0
        Rw[-1] = 0.0
        return Ru, Rw, Rth

    def _acoustic(self, u, w, th, pp, Ru, Rw, Rth, nsteps, dtau):
        """nsteps acoustic steps with the slow forcing held fixed."""
        dx, dz, nz = self.dx, self.dz, self.nz
        ap, am = 0.5 * (1 + self.beta), 0.5 * (1 - self.beta)
        th_full = self.tb_c[:, None] + th
        th_u = 0.5 * (th_full + np.roll(th_full, 1, -1))
        th_w = np.empty((nz + 1, self.nx))
        th_w[1:-1] = 0.5 * (th_full[1:] + th_full[:-1])
        th_w[0], th_w[-1] = th_full[0], th_full[-1]
        C = (self.c2_c / (CP * self.rho_c * self.tb_c ** 2))[:, None]
        rt_c = (self.rho_c * self.tb_c)[:, None]
        rt_w = np.broadcast_to((self.rho_w * self.tb_w)[:, None], w.shape)
        A = dtau * CP * th_w[1:-1] * ap / dz
        Bc = dtau * C * ap / dz
        lower = -A * Bc[:-1] * rt_w[:-2]
        diag = 1.0 + A * (Bc[1:] + Bc[:-1]) * rt_w[1:-1]
        upper = -A * Bc[1:] * rt_w[2:]
        pp_old = pp
        for _ in range(nsteps):
            ppd = pp + self.div_damp * (pp - pp_old)
            pp_old = pp
            u = u + dtau * (Ru - CP * th_u * (ppd - np.roll(ppd, 1, -1)) / dx)
            hdiv = (np.roll(rt_c * u, -1, -1) - rt_c * u) / dx
            vdiv_old = (rt_w[1:] * w[1:] - rt_w[:-1] * w[:-1]) / dz
            pp_star = pp - dtau * C * (hdiv + am * vdiv_old)
            rhs = w[1:-1] + dtau * (Rw[1:-1] - CP * th_w[1:-1] * am * (ppd[1:] - ppd[:-1]) / dz)
            r = rhs - A * (pp_star[1:] - pp_star[:-1])
            w = np.zeros_like(w)
            w[1:-1] = _thomas(lower, diag, upper, r)
            vdiv_new = (rt_w[1:] * w[1:] - rt_w[:-1] * w[:-1]) / dz
            pp = pp_star - dtau * C * ap * vdiv_new
            th = th + dtau * Rth
        return u, w, th, pp

    def step(self, dt):
        """One large step: RK3 with acoustic sub-steps."""
        u0, w0, th0, pp0 = self.u, self.w, self.th, self.pp
        us, ws, ths, pps = u0, w0, th0, pp0
        dtau = dt / self.ns
        for n in (self.ns // 3, self.ns // 2, self.ns):
            Ru, Rw, Rth = self.slow(us, ws, ths)
            us, ws, ths, pps = self._acoustic(u0, w0, th0, pp0, Ru, Rw, Rth, n, dtau)
        self.u, self.w, self.th, self.pp = us, ws, ths, pps
        self.time += dt

    def run(self, t_end, dt, callback=None, every=0):
        n = int(round(t_end / dt))
        for i in range(n):
            self.step(dt)
            if callback is not None and every and (i + 1) % every == 0:
                callback(self)
        return self

    def mass_deviation(self):
        """Domain integral of the linearised density perturbation (diagnostic)."""
        return float(np.sum(self.pp)) * self.dx * self.dz

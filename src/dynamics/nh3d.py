"""
The 3-D non-hydrostatic core in NWP1's sigma (mass) coordinate: CAM roadmap stage S2.

WHAT IT IS

This is the S1b solver (nh2d_mass.py, checked against the mountain-wave and
density-current benchmarks) extended to three dimensions. It is laid out on
NWP1's own grids, so that it can replace PrimitiveSigma inside the existing
forecast pipeline: the same analysis, boundaries, physics, output and
verification.
- Horizontal: CGrid (Arakawa C; replicate edges for a limited area, periodic
  for idealised tests). Coriolis from the grid's f_u and f_v.
- Vertical: SigmaLevels, with index 0 the LID, as everywhere in NWP1. Full
  (mass) levels carry U, V and Theta. Half (w) levels carry W, phi and Omega.
  The ground is half level nz.
- mu = p_s - p_top is NWP1's pi.

VARIABLES
- mu (ny, nx).
- U = mu u and V = mu v, on the u and v points of the mass levels.
- W = mu w (nz+1, ny, nx).
- Theta = mu theta.
- phi = g z on the w levels; phi[nz] = g h is fixed.

EQUATIONS
As in nh2d_mass.py, with the y direction and rotation added:

    dU/dt = -PGF_x + f V(at u) - mu (u . grad u + sigdot du/dsigma) + D
    dV/dt = -PGF_y - f U(at v) - mu (u . grad v + sigdot dv/dsigma) + D

The horizontal pressure gradient uses the perturbation form against a
hydrostatic reference sounding, with the reference balance removed
analytically.

TIME INTEGRATION
RK3 (Wicker and Skamarock 2002) with acoustic sub-steps linearised about
each stage state (Klemp et al. 2007):
- U and V forward;
- then mu, Omega and Theta;
- then a vertically implicit tridiagonal solve for W'' per column, with
  phi'' following;
- a pressure-based divergence damping.

Horizontal advection is upwind3: the momentum in advective form, Theta as a
third-order flux. Vertical advection is second-order centred for momentum
and third-order upwind-biased for the Theta flux.
"""
import numpy as np

from backend import xp_of
from nh2d import G, RD, CP, CV, P0, KAPPA
from nh2d_mass import ReferenceSounding

GAMMA = CP / CV


def thomas(a, b, c, d):
    """Tridiagonal solve along axis 0 (a sub, b diag, c super), any backend."""
    xp = xp_of(d)
    n = b.shape[0]
    cp_ = xp.empty_like(b)
    dp_ = xp.empty_like(d)
    cp_[0] = c[0] / b[0]
    dp_[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp_[i - 1]
        cp_[i] = c[i] / m
        dp_[i] = (d[i] - a[i] * dp_[i - 1]) / m
    x = xp.empty_like(d)
    x[-1] = dp_[-1]
    for i in range(n - 2, -1, -1):
        x[i] = dp_[i] - cp_[i] * x[i + 1]
    return x


class NH3D:
    def __init__(self, grid, levels, terrain=None, theta_ref=None, K=0.0, ns=6,
                 beta=0.1, div_damp=0.1, z_top_ref=None):
        self.grid, self.lev = grid, levels
        gr, lev = grid, levels
        ny, nx, nz = gr.ny, gr.nx, lev.nz
        self.nx, self.ny, self.nz = nx, ny, nz
        self.p_top = float(lev.p_top)
        h = np.zeros((ny, nx)) if terrain is None else np.asarray(terrain, dtype=float)
        self.h = h
        if theta_ref is None:          # isothermal 250 K reference: only has to be hydrostatic
            theta_ref = lambda z: 250.0 * np.exp(G * z / (CP * 250.0))
        zt = 20000.0 if z_top_ref is None else z_top_ref
        self.ref = ReferenceSounding(theta_ref, zt)
        sh = np.asarray(lev.sigma_half, dtype=float)         # 0 (lid) .. 1 (ground)
        sf = np.asarray(lev.sigma, dtype=float)
        self.sh, self.sf = sh, sf
        self.ds = np.diff(sh)                                # > 0
        self.dsw = np.empty(nz + 1)                          # sigma spacing around each w level
        self.dsw[0] = sf[0]
        self.dsw[1:nz] = sf[1:] - sf[:-1]
        self.dsw[nz] = 1.0 - sf[-1]
        # reference state
        self.mub = self.ref.p_of_z(h) - self.p_top
        pw = self.p_top + sh[:, None, None] * self.mub[None]
        self.phib = G * self.ref.z_of_p(pw)
        self.phib[-1] = G * h
        self.alb = (self.phib[:-1] - self.phib[1:]) / (self.mub[None] * self.ds[:, None, None])
        self.pb = self.p_top + sf[:, None, None] * self.mub[None]
        self.thb = (P0 * self.alb / RD) * (self.pb / P0) ** (1.0 / GAMMA)
        self.K, self.ns, self.beta, self.div_damp = float(K), int(ns), float(beta), float(div_damp)
        if self.ns % 6:
            raise ValueError("ns must be a multiple of 6")
        self.f_u = np.asarray(gr.f_u, dtype=float)
        self.f_v = np.asarray(gr.f_v, dtype=float)
        # state (rest, reference)
        self.mu = self.mub.copy()
        self.U = np.zeros((nz, ny, nx))
        self.V = np.zeros((nz, ny, nx))
        self.W = np.zeros((nz + 1, ny, nx))
        self.Th = self.mub[None] * self.thb
        self.phi = self.phib.copy()
        self.time = 0.0
        self.extra_tendency = None       # callable(state) -> dict of additions (physics hook)

    # --- grid helpers (axis 1 = x, 0 = y in CGrid convention) -------------
    def _sh(self, a, n, axis):
        return self.grid.shift(a, n, axis)

    def h_to_u(self, a):
        return 0.5 * (a + self._sh(a, -1, 1))

    def h_to_v(self, a):
        return 0.5 * (a + self._sh(a, -1, 0))

    def u_to_h(self, a):
        return 0.5 * (a + self._sh(a, 1, 1))

    def v_to_h(self, a):
        return 0.5 * (a + self._sh(a, 1, 0))

    def ddx_u(self, a):          # centre field -> d/dx at u points
        return (a - self._sh(a, -1, 1)) / self.grid.dx

    def ddy_v(self, a):
        return (a - self._sh(a, -1, 0)) / self.grid.dy

    def div_h(self, U, V):       # face fluxes -> divergence at centres
        return ((self._sh(U, 1, 1) - U) / self.grid.dx + (self._sh(V, 1, 0) - V) / self.grid.dy)

    def _u3(self, a, vel, axis, d):
        p1, m1 = self._sh(a, 1, axis), self._sh(a, -1, axis)
        p2, m2 = self._sh(a, 2, axis), self._sh(a, -2, axis)
        xp = xp_of(a)
        return (vel * (-p2 + 8 * p1 - 8 * m1 + m2) / (12 * d)
                + xp.abs(vel) * (p2 - 4 * p1 + 6 * a - 4 * m1 + m2) / (12 * d))

    def _flux3(self, F, q, axis):
        """Third-order upwind-biased flux F q at the 'west/south' faces."""
        xp = xp_of(q)
        q0, qm1 = q, self._sh(q, -1, axis)
        qp1, qm2 = self._sh(q, 1, axis), self._sh(q, -2, axis)
        return F * (7.0 * (q0 + qm1) - (qp1 + qm2)) / 12.0 - xp.abs(F) * (3.0 * (q0 - qm1) - (qp1 - qm2)) / 12.0

    def _vface3(self, q, Om):
        """theta on w levels for the vertical flux (index 0 = lid)."""
        xp = xp_of(q)
        nz = self.nz
        qf = xp.empty((nz + 1,) + tuple(q.shape[1:]), dtype=q.dtype) if xp.__name__ != "numpy" else np.empty((nz + 1,) + q.shape[1:])
        qf[0], qf[-1] = q[0], q[-1]
        qf[1:-1] = 0.5 * (q[1:] + q[:-1])
        if nz >= 4:
            a, b, c, d = q[:-3], q[1:-2], q[2:-1], q[3:]   # for face k = 2..nz-2: k-2, k-1 (above), k (below), k+1
            avg = (7.0 * (b + c) - (a + d)) / 12.0
            dif = (3.0 * (c - b) - (d - a)) / 12.0
            # Omega > 0 is downward (sigma increasing): upwind side is above (b): q = avg - dif
            qf[2:-2] = avg - xp.sign(Om[2:-2]) * dif
        return qf

    # --- diagnostics -------------------------------------------------------
    def diagnose(self, mu, Th, phi):
        al = (phi[:-1] - phi[1:]) / (mu[None] * self.ds[:, None, None])
        p = P0 * (RD * Th / (mu[None] * P0 * al)) ** GAMMA
        return al, p

    def omega(self, divV, dmu):
        xp = xp_of(divV)
        Om = xp.zeros((self.nz + 1,) + tuple(divV.shape[1:])) if xp_of(divV).__name__ != "numpy" else np.zeros((self.nz + 1,) + divV.shape[1:])
        for k in range(self.nz):
            Om[k + 1] = Om[k] - (dmu + divV[k]) * self.ds[k]
        Om[-1] = 0.0
        return Om

    def winds(self):
        """u, v at their C-grid points and w on w levels (NumPy views)."""
        return (self.U / self.h_to_u(self.mu)[None], self.V / self.h_to_v(self.mu)[None],
                self.W / self.mu[None])

    def theta(self):
        return self.Th / self.mu[None]

    def total_mass(self):
        return float(np.sum(self.mu)) * self.grid.dx * self.grid.dy

    # --- tendencies --------------------------------------------------------
    def tendencies(self, mu, U, V, W, Th, phi):
        nz = self.nz
        gr = self.grid
        al, p = self.diagnose(mu, Th, phi)
        mu_u, mu_v = self.h_to_u(mu), self.h_to_v(mu)
        u = U / mu_u[None]
        v = V / mu_v[None]
        w = W / mu[None]
        th = Th / mu[None]
        divV = self.div_h(U, V)
        dmu = -(divV * self.ds[:, None, None]).sum(axis=0)
        Om = self.omega(divV, dmu)
        sd = Om / mu[None]
        # pressure gradient, perturbation form
        pp = p - self.pb
        alp = al - self.alb
        mup = mu - self.mub
        phip = phi - self.phib
        phim = 0.5 * (phi[:-1] + phi[1:])
        phipm = 0.5 * (phip[:-1] + phip[1:])
        pw = np.empty((nz + 1,) + mu.shape)
        pw[0] = 0.0
        pw[1:-1] = 0.5 * (pp[1:] + pp[:-1])
        pw[-1] = pp[-1] + (pp[-1] - pp[-2]) * (1.0 - self.sf[-1]) / (self.sf[-1] - self.sf[-2])
        dpp_ds = (pw[1:] - pw[:-1]) / self.ds[:, None, None]
        A1, A2 = mu[None] * al, mu[None] * alp + mup[None] * self.alb
        pgx = (self.h_to_u(A1) * self.ddx_u(pp) + self.h_to_u(A2) * self.ddx_u(self.pb)
               + self.h_to_u(self.mub)[None] * self.ddx_u(phipm) + self.h_to_u(dpp_ds) * self.ddx_u(phim))
        pgy = (self.h_to_v(A1) * self.ddy_v(pp) + self.h_to_v(A2) * self.ddy_v(self.pb)
               + self.h_to_v(self.mub)[None] * self.ddy_v(phipm) + self.h_to_v(dpp_ds) * self.ddy_v(phim))
        # advection of momentum (advective form)
        sd_m = 0.5 * (sd[:-1] + sd[1:])
        def dds_m(a):
            g = np.empty_like(a)
            g[1:-1] = (a[2:] - a[:-2]) / (self.sf[2:] - self.sf[:-2])[:, None, None]
            g[0] = (a[1] - a[0]) / (self.sf[1] - self.sf[0])
            g[-1] = (a[-1] - a[-2]) / (self.sf[-1] - self.sf[-2])
            return g
        v_u = self.h_to_u(self.v_to_h(v))
        u_v = self.h_to_v(self.u_to_h(u))
        FU = (-pgx + self.f_u[None] * self.h_to_u(self.v_to_h(V))
              - mu_u[None] * (self._u3(u, u, 1, gr.dx) + self._u3(u, v_u, 0, gr.dy) + self.h_to_u(sd_m) * dds_m(u)))
        FV = (-pgy - self.f_v[None] * self.h_to_v(self.u_to_h(U))
              - mu_v[None] * (self._u3(v, u_v, 1, gr.dx) + self._u3(v, v, 0, gr.dy) + self.h_to_v(sd_m) * dds_m(v)))
        # Theta flux form
        Fz = Om * self._vface3(th, Om)
        FTh = (-((self._sh(self._flux3(U, th, 1), 1, 1) - self._flux3(U, th, 1)) / gr.dx
                 + (self._sh(self._flux3(V, th, 0), 1, 0) - self._flux3(V, th, 0)) / gr.dy)
               - (Fz[1:] - Fz[:-1]) / self.ds[:, None, None])
        # W and phi
        FW = np.zeros_like(W)
        FW[0] = G * (pp[0] / self.dsw[0] - mup)
        FW[1:-1] = G * ((pp[1:] - pp[:-1]) / self.dsw[1:-1, None, None] - mup[None])
        uc, vc = self.u_to_h(u), self.v_to_h(v)
        def to_w(a):
            o = np.empty((nz + 1,) + a.shape[1:]); o[1:-1] = 0.5 * (a[1:] + a[:-1]); o[0] = a[0]; o[-1] = a[-1]; return o
        u_w, v_w = to_w(uc), to_w(vc)
        dw_ds = np.zeros_like(W)
        dw_ds[1:-1] = (w[2:] - w[:-2]) / (self.sh[2:] - self.sh[:-2])[:, None, None]
        advw = (u_w * (self._sh(w, 1, 1) - self._sh(w, -1, 1)) / (2 * gr.dx)
                + v_w * (self._sh(w, 1, 0) - self._sh(w, -1, 0)) / (2 * gr.dy) + sd * dw_ds)
        FW[:-1] -= (mu[None] * advw)[:-1]
        dphi_ds = np.zeros_like(phi)
        dphi_ds[1:-1] = (phi[2:] - phi[:-2]) / (self.sh[2:] - self.sh[:-2])[:, None, None]
        dphi_ds[0] = (phi[1] - phi[0]) / (self.sh[1] - self.sh[0])
        Fphi = (-(u_w * (self._sh(phi, 1, 1) - self._sh(phi, -1, 1)) / (2 * gr.dx)
                  + v_w * (self._sh(phi, 1, 0) - self._sh(phi, -1, 0)) / (2 * gr.dy)
                  + sd * dphi_ds) + G * w)
        Fphi[-1] = 0.0
        FW[-1] = 0.0
        F = [dmu, FU, FV, FW, FTh, Fphi]
        if self.extra_tendency is not None:
            extra = self.extra_tendency(mu, U, V, W, Th, phi)
            for i, key in enumerate(("mu", "U", "V", "W", "Th", "phi")):
                if key in extra:
                    F[i] = F[i] + extra[key]
        return F

    def _acoustic(self, X0, Xs, F, nsteps, dtau):
        nz, gr = self.nz, self.grid
        mus, Us, Vs, Ws, Ths, phis = Xs
        dmuF, FU, FV, FW, FTh, Fphi = F
        ap, am = 0.5 * (1 + self.beta), 0.5 * (1 - self.beta)
        al, p = self.diagnose(mus, Ths, phis)
        c2 = GAMMA * p * al
        mu_u, mu_v = self.h_to_u(mus), self.h_to_v(mus)
        ths = Ths / mus[None]
        Qc = c2 / al / Ths
        E = c2 / (al ** 2 * mus[None] * self.ds[:, None, None])
        al_u, al_v = self.h_to_u(al), self.h_to_v(al)
        dphi_ds = np.zeros_like(phis)
        dphi_ds[1:-1] = (phis[2:] - phis[:-2]) / (self.sh[2:] - self.sh[:-2])[:, None, None]
        dphi_ds[0] = (phis[1] - phis[0]) / (self.sh[1] - self.sh[0])
        dphidx = (self._sh(phis, 1, 1) - self._sh(phis, -1, 1)) / (2 * gr.dx)
        dphidy = (self._sh(phis, 1, 0) - self._sh(phis, -1, 0)) / (2 * gr.dy)
        Bv = dtau * ap * G / mus
        Gk = dtau * G * ap / self.dsw[:nz, None, None]
        # tridiagonal coefficients (rows k = 0..nz-1)
        lower = np.zeros((nz,) + mus.shape); upper = np.zeros((nz,) + mus.shape)
        diag = 1.0 + Gk * Bv[None] * E
        diag[1:] += Gk[1:] * Bv[None] * E[:-1]
        lower[1:] = -Gk[1:] * Bv[None] * E[:-1]
        upper[:-1] = -Gk[:-1] * Bv[None] * E[:-1 + 0] if False else -Gk[:-1] * Bv[None] * E[:-1]
        upper[-1] = 0.0
        m2 = X0[0] - mus; U2 = X0[1] - Us; V2 = X0[2] - Vs; W2 = X0[3] - Ws; T2 = X0[4] - Ths; f2 = X0[5] - phis
        def pdd(T2, f2):
            return Qc * T2 - E * (f2[:-1] - f2[1:])
        p_old = pdd(T2, f2)
        for _ in range(nsteps):
            p2 = pdd(T2, f2)
            pd = p2 + self.div_damp * (p2 - p_old)
            p_old = p2
            f2m = 0.5 * (f2[:-1] + f2[1:])
            U2 = U2 + dtau * (FU - (mu_u[None] * al_u * self.ddx_u(pd) + mu_u[None] * self.ddx_u(f2m)))
            V2 = V2 + dtau * (FV - (mu_v[None] * al_v * self.ddy_v(pd) + mu_v[None] * self.ddy_v(f2m)))
            divV = self.div_h(U2, V2)
            dmu2 = -(divV * self.ds[:, None, None]).sum(axis=0)
            m2 = m2 + dtau * (dmuF + dmu2)
            Om2 = self.omega(divV, dmu2)
            Fz = Om2 * self._vface_lin(ths)
            T2 = T2 + dtau * (FTh - ((self._sh(U2 * self.h_to_u(ths), 1, 1) - U2 * self.h_to_u(ths)) / gr.dx
                                     + (self._sh(V2 * self.h_to_v(ths), 1, 0) - V2 * self.h_to_v(ths)) / gr.dy)
                              - (Fz[1:] - Fz[:-1]) / self.ds[:, None, None])
            U2c, V2c = self.u_to_h(U2), self.v_to_h(V2)
            def to_w(a):
                o = np.empty((nz + 1,) + a.shape[1:]); o[1:-1] = 0.5 * (a[1:] + a[:-1]); o[0] = a[0]; o[-1] = a[-1]; return o
            Phi = f2 + dtau * (Fphi - (to_w(U2c) * dphidx + to_w(V2c) * dphidy + Om2 * dphi_ds) / mus[None]
                               + am * G * W2 / mus[None])
            Phi[-1] = 0.0
            dpo = np.zeros_like(W2)
            dpo[0] = pd[0] / self.dsw[0]
            dpo[1:-1] = (pd[1:] - pd[:-1]) / self.dsw[1:-1, None, None]
            Wexp = W2 + dtau * (FW + G * (am * dpo - m2[None]))
            Q = Qc * T2
            S = np.empty((nz,) + mus.shape)
            S[0] = Q[0] - E[0] * (Phi[0] - Phi[1])
            S[1:] = (Q[1:] - E[1:] * (Phi[1:-1] - Phi[2:])) - (Q[:-1] - E[:-1] * (Phi[:-2] - Phi[1:-1]))
            r = Wexp[:-1] + Gk * S
            Wn = thomas(lower, diag, upper, r)
            W2 = np.zeros_like(W2); W2[:-1] = Wn
            f2 = Phi + Bv[None] * W2
            f2[-1] = 0.0
        return (mus + m2, Us + U2, Vs + V2, Ws + W2, Ths + T2, phis + f2)

    def _vface_lin(self, th):
        o = np.empty((self.nz + 1,) + th.shape[1:])
        o[1:-1] = 0.5 * (th[1:] + th[:-1]); o[0] = th[0]; o[-1] = th[-1]
        return o

    def step(self, dt):
        X0 = (self.mu, self.U, self.V, self.W, self.Th, self.phi)
        Xs = X0
        dtau = dt / self.ns
        for n in (self.ns // 3, self.ns // 2, self.ns):
            F = self.tendencies(*Xs)
            Xs = self._acoustic(X0, Xs, F, n, dtau)
        self.mu, self.U, self.V, self.W, self.Th, self.phi = Xs
        self.time += dt


class NHModel:
    """
    The non-hydrostatic core behind PrimitiveSigma's interface (stage S2).

    forecast.run_forecast drives a model through step(dt), time, u/v/theta/pi,
    lev, max_dt() and sigma_dot(), plus a relaxation hook. This class provides
    exactly that, so the same analysis, boundaries, output and verification
    are used for both cores. Physics comes from the hydrostatic model it was
    built from:
    - surface drag and Richardson mixing as tendencies (multiplied by mu);
    - the pool-adjacent-violators adjustment after each step;
    - the same top sponge on u and v, and a Rayleigh damping of w in the
      sponge levels.
    Prescribed heating and the land surface are not wired in yet.
    """

    def __init__(self, hydro, theta_ref=None, ns=6, dt_max=60.0, div_damp=0.1):
        self.h = hydro
        self.grid, self.lev = hydro.grid, hydro.lev
        self.core = NH3D(hydro.grid, hydro.lev, terrain=hydro.terrain, theta_ref=theta_ref,
                         ns=ns, div_damp=div_damp)
        self.core.extra_tendency = self._physics
        self.dt_max = float(dt_max)
        self.land_surface = None
        self.surface_heating = None
        self._u_ref = self._v_ref = None
        self.time = 0.0
        self.step_count = 0
        self._last_dt = None
        # RELAXATION PER UNIT TIME, NOT PER STEP. The Davies weights in
        # forecast.py are applied once per step, and their values were set
        # with the hydrostatic core's ~17 s step. At 60 s the same weights
        # would relax ~3.5x more weakly per hour. The NH adapter rescales them
        # to the same e-folding time: a_eff = 1 - (1 - a)**(dt / dt_ref).
        try:
            self.relax_dt_ref = float(hydro.max_dt())
        except Exception:
            self.relax_dt_ref = None
        nz = self.lev.nz
        self._wdamp = np.zeros((nz + 1, 1, 1))
        self._wdamp[:nz] = np.asarray(hydro._sponge, dtype=float) * 2.0

    # --- state as the hydrostatic interface sees it ------------------------
    @property
    def pi(self):
        return self.core.mu

    @property
    def surface_pressure(self):
        return self.lev.p_top + self.core.mu

    @property
    def u(self):
        return self.core.U / self.core.h_to_u(self.core.mu)[None]

    @property
    def v(self):
        return self.core.V / self.core.h_to_v(self.core.mu)[None]

    @property
    def theta(self):
        return self.core.Th / self.core.mu[None]

    def set_state(self, u, v, theta, pi):
        """Load a hydrostatic state; W = 0 and phi integrated hydrostatically."""
        c = self.core
        c.mu = np.array(pi, dtype=float, copy=True)
        c.U = np.asarray(u, dtype=float) * c.h_to_u(c.mu)[None]
        c.V = np.asarray(v, dtype=float) * c.h_to_v(c.mu)[None]
        c.Th = np.asarray(theta, dtype=float) * c.mu[None]
        c.W = np.zeros_like(c.W)
        c.phi = self.hydrostatic_phi(c.mu, np.asarray(theta, dtype=float))
        self._u_ref, self._v_ref = np.array(u, copy=True), np.array(v, copy=True)

    def hydrostatic_phi(self, mu, theta):
        c = self.core
        p = c.p_top + c.sf[:, None, None] * mu[None]
        alpha = RD * theta * (p / P0) ** KAPPA / p
        phi = np.empty((self.lev.nz + 1,) + mu.shape)
        phi[-1] = G * np.asarray(self.h.terrain, dtype=float)
        for k in range(self.lev.nz - 1, -1, -1):
            phi[k] = phi[k + 1] + alpha[k] * mu * c.ds[k]
        return phi

    def relax_with(self, relax, ext):
        """Davies relaxation of the NH state toward the driving (hydrostatic) frame."""
        c = self.core
        a2 = np.asarray(relax.alpha2d, dtype=float)
        if self.relax_dt_ref and self._last_dt:
            a2 = 1.0 - (1.0 - a2) ** (self._last_dt / self.relax_dt_ref)
        a3 = a2[None]
        u, v, th = self.u, self.v, self.theta
        if "u" in ext:
            u = u + a3 * (np.asarray(ext["u"]) - u)
        if "v" in ext:
            v = v + a3 * (np.asarray(ext["v"]) - v)
        if "theta" in ext:
            th = th + a3 * (np.asarray(ext["theta"]) - th)
        mu = c.mu
        if "pi" in ext:
            mu = mu + a2 * (np.asarray(ext["pi"]) - mu)
        c.mu = mu
        c.U = u * c.h_to_u(mu)[None]
        c.V = v * c.h_to_v(mu)[None]
        c.Th = th * mu[None]
        # geopotential toward hydrostatic balance and w toward zero in the zone
        phi_h = self.hydrostatic_phi(mu, th)
        c.phi = c.phi + a3 * (phi_h - c.phi)
        c.W = c.W * (1.0 - a3)

    # --- physics as tendencies ----------------------------------------------
    def _physics(self, mu, U, V, W, Th, phi):
        from surface import surface_drag
        from turbulence import vertical_mixing
        c, hy = self.core, self.h
        mu_u, mu_v = c.h_to_u(mu), c.h_to_v(mu)
        u, v, th = U / mu_u[None], V / mu_v[None], Th / mu[None]
        FU = np.zeros_like(U); FV = np.zeros_like(V); FT = np.zeros_like(Th)
        if hy.drag:
            du, dv, info = surface_drag(u, v, th, mu, self.lev, z0=hy.z0, theta_s=hy.theta_surface)
            FU += mu_u[None] * du
            FV += mu_v[None] * dv
        if hy.mixing:
            mu_, mv_, mth, K = vertical_mixing(u, v, th, mu, self.lev, ri_crit=hy.ri_crit,
                                               k_max=hy.k_max, mixing_length=hy.mixing_length)
            FU += mu_u[None] * mu_
            FV += mu_v[None] * mv_
            FT += mu[None] * mth
        if self._u_ref is not None and hy.sponge_levels > 0:
            sp = np.asarray(hy._sponge, dtype=float)
            FU -= sp * (U - mu_u[None] * self._u_ref)
            FV -= sp * (V - mu_v[None] * self._v_ref)
        return {"U": FU, "V": FV, "Th": FT, "W": -self._wdamp * W}

    def step(self, dt):
        self.core.step(dt)
        self._last_dt = float(dt)
        hy = self.h
        if hy.convection:
            from convection import dry_convective_adjustment_pav
            th, u, v, info = dry_convective_adjustment_pav(self.theta, self.u, self.v, self.core.mu,
                                                           self.lev, mix_momentum=hy.conv_mix_momentum)
            c = self.core
            c.Th = th * c.mu[None]
            c.U = u * c.h_to_u(c.mu)[None]
            c.V = v * c.h_to_v(c.mu)[None]
        self.time += dt
        self.step_count += 1

    def max_dt(self, safety=0.7):
        umax = float(max(np.abs(self.u).max(), np.abs(self.v).max(), 10.0))
        return min(self.dt_max, safety * min(self.grid.dx, self.grid.dy) / umax)

    def sigma_dot(self):
        c = self.core
        divV = c.div_h(c.U, c.V)
        dmu = -(divV * c.ds[:, None, None]).sum(axis=0)
        return c.omega(divV, dmu) / c.mu[None]

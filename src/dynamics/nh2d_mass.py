"""
The 2-D non-hydrostatic core in the hydrostatic-pressure (mass) coordinate: CAM roadmap stage S1b.

WHY

This is the coordinate the roadmap chose (docs/CAM_DESIGN.md): Laprise
(1992), as in WRF-ARW (Skamarock and Klemp 2008). It is NWP1's sigma
(p_h = p_top + sigma mu) made non-hydrostatic. It follows the terrain, keeps
the existing analysis conversion, and carries the column mass mu in flux
form, so mass is conserved to round-off.

Stage S1a (nh2d.py) checked the time-splitting and the vertically implicit
acoustic solver on flat ground in height coordinates. This file puts the
same structure on the mass coordinate.

VARIABLES (x-sigma; periodic in x; sigma = 1 at the ground, 0 at the lid)
- mu(x): the column mass p_s - p_top, at cell centres.
- U = mu u on x faces and mass levels.
- W = mu w at cell centres on w levels.
- Theta = mu theta at centres and mass levels.
- phi = g z, the geopotential, at centres on w levels.
- Diagnosed: the inverse density alpha_k = (phi_{k+1} - phi_k) / (mu dsigma_k),
  and the pressure p = p0 (R Theta / (mu p0 alpha))^(cp/cv). The pressure
  is p_top at the lid.

EQUATIONS (dry, no rotation)

    dU/dt     = -mu alpha dp/dx - (dp/dsigma) dphi/dx + advection + D
    dW/dt     =  g (dp/dsigma - mu)                   + advection + D
    dTheta/dt = -d(U theta)/dx - d(Omega theta)/dsigma            + D
    dmu/dt    = -sum_k d(U_k)/dx dsigma_k
    dphi/dt   = -(U dphi/dx + Omega dphi/dsigma)/mu + g W/mu

- Omega = mu dsigma/dt comes from the continuity equation.
- The horizontal pressure gradient uses a hydrostatic reference state
  (p_bar = p_top + sigma mu_bar, phi_bar from the reference sounding). Its
  balance mu_bar alpha_bar dp_bar/dx + mu_bar dphi_bar/dx = 0 is removed
  analytically, so the large opposing terms over terrain cancel exactly:

    PGF = mu alpha dp'/dx + (mu alpha' + mu' alpha_bar) dp_bar/dx
          + mu_bar dphi'/dx + (dp'/dsigma) dphi/dx

TIME INTEGRATION
RK3 (Wicker and Skamarock 2002) with acoustic sub-steps (Klemp et al. 2007).
- In each stage the perturbations X'' = X - X* from the stage state X*
  start at X^t - X* and are integrated with the full tendency F(X*) plus the
  acoustic operator linearised about X*.
- U is forward, then mu, Omega and Theta.
- W and phi are implicit in the vertical, with off-centring beta: a
  tridiagonal solve for W'' per column, with the pressure linearised as
  p''_k = (c^2/alpha) Theta''/Theta - c^2/(alpha^2 mu dsigma_k) (phi''_{k+1} - phi''_k).
- A pressure-based divergence damping is applied, as in S1a.
"""
import numpy as np

from nh2d import G, RD, CP, CV, P0, KAPPA, _thomas

GAMMA = CP / CV


class ReferenceSounding:
    """theta_bar(z) integrated hydrostatically on a fine grid; z(p) by interpolation in ln p."""

    def __init__(self, theta_of_z, z_top, p_surface=P0, dz=5.0):
        z = np.arange(0.0, z_top + 2000.0 + dz, dz)
        th = theta_of_z(z) * np.ones_like(z)
        pi = np.empty_like(z)
        pi[0] = (p_surface / P0) ** KAPPA
        for i in range(1, len(z)):
            pi[i] = pi[i - 1] - G * dz / (CP * 0.5 * (th[i] + th[i - 1]))
        self.z, self.p = z, P0 * pi ** (1.0 / KAPPA)
        self.theta_of_z = theta_of_z

    def p_of_z(self, z):
        return np.exp(np.interp(z, self.z, np.log(self.p)))

    def z_of_p(self, p):
        return np.interp(-np.log(p), -np.log(self.p), self.z)


class NH2DMass:
    def __init__(self, nx, nz, dx, z_top, theta_of_z, terrain=None, K=0.0, ns=6,
                 beta=0.1, div_damp=0.1, sponge_depth=0.0, sponge_rate=0.0, u0=0.0):
        self.nx, self.nz, self.dx = nx, nz, float(dx)
        self.ref = ReferenceSounding(theta_of_z, z_top)
        self.p_top = float(self.ref.p_of_z(z_top))
        h = np.zeros(nx) if terrain is None else np.asarray(terrain, dtype=float)
        self.h = h
        # sigma levels from the reference pressure at equally spaced heights over flat ground
        z_w = np.linspace(0.0, z_top, nz + 1)
        p_w0 = self.ref.p_of_z(z_w)
        self.sw = (p_w0 - self.p_top) / (p_w0[0] - self.p_top)          # 1 ... 0
        self.sw[0], self.sw[-1] = 1.0, 0.0
        self.ds = self.sw[:-1] - self.sw[1:]                             # > 0
        self.sm = 0.5 * (self.sw[:-1] + self.sw[1:])
        self.dsw = np.empty(nz + 1)                                      # sigma distance between mass levels around w level k
        self.dsw[1:-1] = self.sm[:-1] - self.sm[1:]
        self.dsw[0] = 1.0 - self.sm[0]
        self.dsw[-1] = self.sm[-1]
        # reference state
        self.mub = self.ref.p_of_z(h) - self.p_top                      # (nx,)
        pw = self.p_top + self.sw[:, None] * self.mub[None, :]
        self.phib = G * self.ref.z_of_p(pw)                              # (nz+1, nx)
        self.phib[0] = G * h
        self.alb = (self.phib[1:] - self.phib[:-1]) / (self.mub[None, :] * self.ds[:, None])
        self.pb = self.p_top + self.sm[:, None] * self.mub[None, :]
        # theta consistent with the discrete EOS so that p' = 0 at rest
        self.thb = (P0 * self.alb / RD) * (self.pb / P0) ** (1.0 / GAMMA)
        self.K, self.ns, self.beta, self.div_damp = float(K), int(ns), float(beta), float(div_damp)
        if self.ns % 6:
            raise ValueError("ns must be a multiple of 6")
        # sponge: Rayleigh damping toward the reference/background in the top layers
        zc = 0.5 * (self.phib[1:] + self.phib[:-1]) / G
        zb = z_top - sponge_depth
        self.rd_m = np.where(zc > zb, sponge_rate * np.sin(0.5 * np.pi * (zc - zb) / max(sponge_depth, 1.0)) ** 2, 0.0)
        zw = self.phib / G
        self.rd_w = np.where(zw > zb, sponge_rate * np.sin(0.5 * np.pi * (zw - zb) / max(sponge_depth, 1.0)) ** 2, 0.0)
        self.u0 = float(u0)
        # state
        self.mu = self.mub.copy()
        self.U = np.full((nz, nx), u0) * self._xf(self.mu)[None, :]
        self.W = np.zeros((nz + 1, nx))
        self.Th = self.mub[None, :] * self.thb
        self.phi = self.phib.copy()
        self.time = 0.0

    # --- helpers ---------------------------------------------------------
    @staticmethod
    def _xf(a):
        """centre -> west face average (periodic)."""
        return 0.5 * (a + np.roll(a, 1, -1))

    @staticmethod
    def _xc(a):
        """west face -> centre average."""
        return 0.5 * (a + np.roll(a, -1, -1))

    def _dxc(self, a):
        """d/dx of a face field, at centres."""
        return (np.roll(a, -1, -1) - a) / self.dx

    def _dxf(self, a):
        """d/dx of a centre field, at west faces."""
        return (a - np.roll(a, 1, -1)) / self.dx

    def diagnose(self, mu, Th, phi):
        al = (phi[1:] - phi[:-1]) / (mu[None, :] * self.ds[:, None])
        p = P0 * (RD * Th / (mu[None, :] * P0 * al)) ** GAMMA
        return al, p

    def omega(self, U, dmu_dt):
        """Omega on w levels from continuity (zero at the ground and the lid)."""
        div = self._dxc(U)                                               # (nz, nx)
        Om = np.zeros((self.nz + 1, self.nx))
        for k in range(self.nz - 1, -1, -1):                             # from the lid down
            Om[k] = Om[k + 1] - (dmu_dt + div[k]) * self.ds[k]
        Om[0] = 0.0
        return Om

    @staticmethod
    def _u3(a, vel, dx):
        p1, m1 = np.roll(a, -1, -1), np.roll(a, 1, -1)
        p2, m2 = np.roll(a, -2, -1), np.roll(a, 2, -1)
        return vel * (-p2 + 8 * p1 - 8 * m1 + m2) / (12 * dx) + np.abs(vel) * (p2 - 4 * p1 + 6 * a - 4 * m1 + m2) / (12 * dx)

    @staticmethod
    def _flux3(U, q):
        """Third-order upwind-biased flux U q at the west faces (periodic)."""
        q0, qm1 = q, np.roll(q, 1, -1)
        qp1, qm2 = np.roll(q, -1, -1), np.roll(q, 2, -1)
        avg = (7.0 * (q0 + qm1) - (qp1 + qm2)) / 12.0
        dif = (3.0 * (q0 - qm1) - (qp1 - qm2)) / 12.0
        return U * avg - np.abs(U) * dif

    def _vface3(self, q, Om):
        """theta on w levels for the vertical flux: third-order upwind-biased in the interior,
        centred next to the boundaries (Omega is zero at the ground and the lid)."""
        nz = self.nz
        qf = np.empty((nz + 1, q.shape[1]))
        qf[0], qf[-1] = q[0], q[-1]
        qf[1:-1] = 0.5 * (q[1:] + q[:-1])
        if nz >= 4:
            # face k (2..nz-2) between mass levels k-1 (below) and k (above);
            # Omega > 0 is downward (sigma increasing), so the upwind side is above.
            a, b, c, d = q[:-3], q[1:-2], q[2:-1], q[3:]               # k-2, k-1, k, k+1 for k = 2..nz-2
            avg = (7.0 * (b + c) - (a + d)) / 12.0
            dif = (3.0 * (c - b) - (d - a)) / 12.0
            # for upward motion (Omega < 0) the upwind cell is below (b): q = avg - dif;
            # for downward motion (Omega > 0) it is above (c): q = avg + dif
            s_ = np.sign(Om[2:-2])
            qf[2:-2] = avg + s_ * dif
        return qf

    def tendencies(self, mu, U, W, Th, phi):
        """Full tendencies F(X) for (mu, U, W, Th, phi)."""
        nz, dx = self.nz, self.dx
        al, p = self.diagnose(mu, Th, phi)
        mu_f = self._xf(mu)
        u = U / mu_f[None, :]
        w = W / mu[None, :]
        th = Th / mu[None, :]
        dmu = -np.sum(self._dxc(U) * self.ds[:, None], axis=0)
        Om = self.omega(U, dmu)
        sd = Om / mu[None, :]                                             # sigma-dot on w levels
        # --- pressure gradient (perturbation form)
        pp = p - (self.p_top + self.sm[:, None] * self.mub[None, :])    # p'
        alp = al - self.alb
        mup = mu - self.mub
        phip = phi - self.phib
        pb = self.pb
        phim = 0.5 * (phi[1:] + phi[:-1])
        phipm = 0.5 * (phip[1:] + phip[:-1])
        # dp'/dsigma at mass levels (centred, one-sided at ends; p' = 0 at the lid face)
        pw = np.empty((nz + 1, self.nx))
        pw[1:-1] = 0.5 * (pp[1:] + pp[:-1])
        pw[0] = pp[0] + (pp[0] - pp[1]) * (1.0 - self.sm[0]) / (self.sm[0] - self.sm[1]) if nz > 1 else pp[0]
        pw[-1] = 0.0
        dpp_ds = (pw[:-1] - pw[1:]) / self.ds[:, None]
        pgf = (self._xf(mu * al) * self._dxf(pp)
               + self._xf(mu * alp + mup * self.alb) * self._dxf(pb)
               + self._xf(self.mub)[None, :] * self._dxf(phipm)
               + self._xf(dpp_ds) * self._dxf(phim))
        # --- U: advective form, mass weighted
        sd_m = 0.5 * (sd[1:] + sd[:-1])
        def dds_m(a):
            g = np.empty_like(a)
            g[1:-1] = (a[:-2] - a[2:]) / (self.sm[:-2] - self.sm[2:])[:, None]
            g[0] = (a[0] - a[1]) / (self.sm[0] - self.sm[1])
            g[-1] = (a[-2] - a[-1]) / (self.sm[-2] - self.sm[-1])
            return g
        w_m = 0.5 * (w[1:] + w[:-1])
        FU = -pgf - mu_f[None, :] * (self._u3(u, u, dx) + self._xf(sd_m) * dds_m(u))
        # --- Theta: flux form, centred
        thf = self._xf(th)
        thw = np.empty((nz + 1, self.nx))
        thw[1:-1] = 0.5 * (th[1:] + th[:-1])
        thw[0], thw[-1] = th[0], th[-1]
        # -d(Omega theta)/dsigma: sigma decreases upward, so it is (F_{k+1} - F_k) / ds_k
        FTh = -self._dxc(self._flux3(U, th)) + ((Om * self._vface3(th, Om))[1:] - (Om * self._vface3(th, Om))[:-1]) / self.ds[:, None]
        # --- W and phi on w levels
        FW = np.zeros_like(W)
        dpw = (pp[:-1] - pp[1:]) / self.dsw[1:-1, None]                  # dp'/dsigma at interior w levels
        FW[1:-1] = G * (dpw - mup[None, :])
        FW[-1] = G * ((pp[-1] - 0.0) / self.dsw[-1] - mup)
        u_w = np.empty_like(W)
        uc = self._xc(u)
        u_w[1:-1] = 0.5 * (uc[1:] + uc[:-1]); u_w[0] = uc[0]; u_w[-1] = uc[-1]
        dw_ds = np.zeros_like(W)
        dw_ds[1:-1] = (w[:-2] - w[2:]) / (self.sw[:-2] - self.sw[2:])[:, None]
        FW[1:] -= (mu[None, :] * (self._u3(w, u_w, dx) + sd * dw_ds))[1:]
        dphi_ds = np.zeros_like(phi)
        dphi_ds[1:-1] = (phi[:-2] - phi[2:]) / (self.sw[:-2] - self.sw[2:])[:, None]
        dphi_ds[-1] = (phi[-2] - phi[-1]) / (self.sw[-2] - self.sw[-1])
        dphi_dx = (np.roll(phi, -1, -1) - np.roll(phi, 1, -1)) / (2 * dx)
        Fphi = -(u_w * dphi_dx + sd * dphi_ds) + G * w
        Fphi[0] = 0.0
        # --- diffusion (on coordinate surfaces, physical vertical spacing)
        if self.K > 0:
            lap = lambda a: (np.roll(a, -1, -1) - 2 * a + np.roll(a, 1, -1)) / dx ** 2
            dzm = (phi[1:] - phi[:-1]) / G
            def lapz(a, dz):
                up = np.concatenate([a[1:], a[-1:]], 0); dn = np.concatenate([a[:1], a[:-1]], 0)
                return (up - 2 * a + dn) / dz ** 2
            FU += mu_f[None, :] * self.K * (lap(u) + lapz(u, self._xf(dzm)))
            FTh += mu[None, :] * self.K * (lap(th - self.thb) + lapz(th - self.thb, dzm))
            dzw = np.concatenate([dzm[:1], 0.5 * (dzm[1:] + dzm[:-1]), dzm[-1:]], 0)
            lw = lap(w) + lapz(w, dzw); lw[0] = 0; lw[-1] = 0
            FW += mu[None, :] * self.K * lw
        # --- sponge
        if np.any(self.rd_m > 0):
            FU -= self._xf(self.rd_m) * (U - self.u0 * mu_f[None, :])
            FTh -= self.rd_m * (Th - mu[None, :] * self.thb)
            FW -= self.rd_w * W
        FW[0] = 0.0
        return dmu, FU, FW, FTh, Fphi

    def step(self, dt):
        X0 = (self.mu, self.U, self.W, self.Th, self.phi)
        Xs = X0
        dtau = dt / self.ns
        for n in (self.ns // 3, self.ns // 2, self.ns):
            F = self.tendencies(*Xs)
            Xs = self._acoustic(X0, Xs, F, n, dtau)
        self.mu, self.U, self.W, self.Th, self.phi = Xs
        self.time += dt

    def _acoustic(self, X0, Xs, F, nsteps, dtau):
        nz, dx = self.nz, self.dx
        mus, Us, Ws, Ths, phis = Xs
        dmuF, FU, FW, FTh, Fphi = F
        ap, am = 0.5 * (1 + self.beta), 0.5 * (1 - self.beta)
        al, p = self.diagnose(mus, Ths, phis)
        c2 = GAMMA * p * al
        mu_f = self._xf(mus)
        ths = Ths / mus[None, :]
        thf = self._xf(ths)
        thw = np.empty((nz + 1, self.nx)); thw[1:-1] = 0.5 * (ths[1:] + ths[:-1]); thw[0], thw[-1] = ths[0], ths[-1]
        Qc = c2 / al / Ths                                              # p'' per Theta''
        E = c2 / (al ** 2 * mus[None, :] * self.ds[:, None])             # p'' per phi'' difference
        phim = 0.5 * (phis[1:] + phis[:-1])
        dphim_dx = self._dxf(phim)
        dphi_ds = np.zeros_like(phis)
        dphi_ds[1:-1] = (phis[:-2] - phis[2:]) / (self.sw[:-2] - self.sw[2:])[:, None]
        dphi_ds[-1] = (phis[-2] - phis[-1]) / (self.sw[-2] - self.sw[-1])
        Bv = dtau * ap * G / mus                                         # (nx,)
        Gk = dtau * G * ap / self.dsw[:, None]                           # (nz+1, 1)
        # perturbations relative to the stage state, starting at X^t - X*
        m2 = X0[0] - mus; U2 = X0[1] - Us; W2 = X0[2] - Ws; T2 = X0[3] - Ths; f2 = X0[4] - phis
        def pdd(T2, f2):
            return Qc * T2 - E * (f2[1:] - f2[:-1])
        p_old = pdd(T2, f2)
        for _ in range(nsteps):
            p2 = pdd(T2, f2)
            pd = p2 + self.div_damp * (p2 - p_old)
            p_old = p2
            # U
            U2 = U2 + dtau * (FU - (mu_f[None, :] * self._xf(al) * self._dxf(pd) + self._xf(mus)[None, :] * self._dxf(0.5 * (f2[1:] + f2[:-1]))))
            # mu, Omega, Theta
            divU = self._dxc(U2)
            dmu2 = -np.sum(divU * self.ds[:, None], axis=0)
            m2 = m2 + dtau * (dmuF + dmu2)
            Om2 = np.zeros((nz + 1, self.nx))
            for k in range(nz - 1, -1, -1):
                Om2[k] = Om2[k + 1] - (dmu2 + divU[k]) * self.ds[k]
            Om2[0] = 0.0
            T2 = T2 + dtau * (FTh - self._dxc(U2 * thf) + ((Om2 * thw)[1:] - (Om2 * thw)[:-1]) / self.ds[:, None])
            # explicit parts of phi'' and W''
            dphi_dx_w = (np.roll(phis, -1, -1) - np.roll(phis, 1, -1)) / (2 * dx)
            U2w = np.empty_like(W2); U2c = self._xc(U2)
            U2w[1:-1] = 0.5 * (U2c[1:] + U2c[:-1]); U2w[0] = U2c[0]; U2w[-1] = U2c[-1]
            Phi = f2 + dtau * (Fphi - (U2w * dphi_dx_w + Om2 * dphi_ds) / mus[None, :] + am * G * W2 / mus[None, :])
            Phi[0] = 0.0
            dp_old = np.zeros_like(W2)
            dp_old[1:-1] = (pd[:-1] - pd[1:]) / self.dsw[1:-1, None]
            dp_old[-1] = pd[-1] / self.dsw[-1]
            Wexp = W2 + dtau * (FW + G * (am * dp_old - m2[None, :]))
            # tridiagonal for W''_k, k = 1..nz
            Q = Qc * T2
            S = np.zeros_like(W2)
            S[1:-1] = (Q[:-1] - E[:-1] * (Phi[1:-1] - Phi[:-2])) - (Q[1:] - E[1:] * (Phi[2:] - Phi[1:-1]))
            S[-1] = Q[-1] - E[-1] * (Phi[-1] - Phi[-2])
            lower = np.zeros((nz, self.nx)); diag = np.ones((nz, self.nx)); upper = np.zeros((nz, self.nx))
            # interior k = 1..nz-1 -> rows 0..nz-2
            g_ = Gk[1:-1]
            diag[:-1] = 1.0 + g_ * Bv[None, :] * (E[:-1] + E[1:])
            lower[:-1] = -g_ * Bv[None, :] * E[:-1]
            upper[:-1] = -g_ * Bv[None, :] * E[1:]
            lower[0] = 0.0                                                 # W''_0 = 0
            # top k = nz -> row nz-1
            diag[-1] = 1.0 + Gk[-1] * Bv * E[-1]
            lower[-1] = -Gk[-1] * Bv * E[-1]
            r = Wexp[1:] + Gk[1:] * S[1:]
            Wn = _thomas(lower, diag, upper, r)
            W2 = np.zeros_like(W2); W2[1:] = Wn
            f2 = Phi + Bv[None, :] * W2
            f2[0] = 0.0
        return (mus + m2, Us + U2, Ws + W2, Ths + T2, phis + f2)

    # --- diagnostics -------------------------------------------------------
    def total_mass(self):
        return float(np.sum(self.mu)) * self.dx

    def w_field(self):
        return self.W / self.mu[None, :]

    def u_centres(self):
        return self._xc(self.U / self._xf(self.mu)[None, :])

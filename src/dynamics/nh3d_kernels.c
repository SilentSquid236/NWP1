/*
 * Compiled kernels for the non-hydrostatic core (nh3d.py) -- CAM stage S4a.
 *
 * WHY C. Test S3a: the NumPy core needs 42 h for a 24 h forecast at 3 km x 40
 * levels on one core; the budget is 1.5 h. Test S3b: on one stencil of this
 * kind, a fused C loop with OpenMP is 216x faster than NumPy (float32, 26
 * threads). NumPy makes one pass over memory per operation; a fused loop reads
 * each input once. The user chose this route (prompt 170).
 *
 * RULES.
 *  - Every kernel here is a transcription of a NumPy expression in nh3d.py, in
 *    the same order of operations, and is tested against it
 *    (test_nh3d_c.py). The NumPy code stays the reference.
 *  - Nothing is installed: the server's gcc compiles this file at first use
 *    (cnh.py), and Python calls it through ctypes.
 *  - Arrays are C-contiguous (nz, ny, nx), level 0 = lid, as in nh3d.py.
 *    w-level arrays have nz+1 levels. Horizontal neighbours come from index
 *    tables (xm1, xp1, ym1, yp1) built from the grid's own shift(), so the
 *    periodic and replicate edge modes are both exact.
 *  - REAL is double by default; compiling with -DREAL=float -DSFX=_f32 gives
 *    the float32 set.
 */
#include <math.h>
#include <stddef.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#ifndef REAL
#define REAL double
#define SFX _f64
#endif
#define CAT_(a, b) a##b
#define CAT(a, b) CAT_(a, b)
#define FN(name) CAT(name, SFX)

#if defined(_WIN32)
#define EXPORT __declspec(dllexport)
#else
#define EXPORT
#endif

EXPORT int FN(nwp_set_threads)(int n)
{
#ifdef _OPENMP
    if (n > 0) omp_set_num_threads(n);
    return omp_get_max_threads();
#else
    (void)n;
    return 1;
#endif
}

/*
 * Acoustic substep, part 1 (nh3d._acoustic, pdd and divergence damping):
 *   p2 = Qc*T2 - E*(f2[k] - f2[k+1]);  pd = p2 + damp*(p2 - p_old);  p_old = p2
 */
EXPORT void FN(ac_pd)(int nz, int ny, int nx,
                      const REAL *Qc, const REAL *T2, const REAL *E, const REAL *f2,
                      REAL *p_old, REAL damp, REAL *pd)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel for schedule(static)
    for (int k = 0; k < nz; k++) {
        const size_t o = (size_t)k * N;
        for (size_t c = o; c < o + N; c++) {
            REAL p2 = Qc[c] * T2[c] - E[c] * (f2[c] - f2[c + N]);
            pd[c] = p2 + damp * (p2 - p_old[c]);
            p_old[c] = p2;
        }
    }
}

/*
 * Acoustic substep, part 2: the horizontal momentum update.
 *   U2 += dtau*(FU - (cu*ddx_u(pd) + mu_u*ddx_u(f2m)))   with cu = mu_u*al_u
 *   V2 += dtau*(FV - (cv*ddy_v(pd) + mu_v*ddy_v(f2m)))   f2m = (f2[k]+f2[k+1])/2
 */
EXPORT void FN(ac_uv)(int nz, int ny, int nx, const int *xm1, const int *ym1,
                      REAL dtau, REAL dx, REAL dy,
                      const REAL *FU, const REAL *FV,
                      const REAL *cu, const REAL *mu_u, const REAL *cv, const REAL *mu_v,
                      const REAL *pd, const REAL *f2, REAL *U2, REAL *V2)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++) {
            const size_t row = (size_t)k * N + (size_t)j * nx;
            const size_t rows = (size_t)k * N + (size_t)ym1[j] * nx;
            const size_t h = (size_t)j * nx;
            for (int i = 0; i < nx; i++) {
                const size_t c = row + i, cw = row + xm1[i], cs = rows + i;
                REAL f2m = (REAL)0.5 * (f2[c] + f2[c + N]);
                REAL f2w = (REAL)0.5 * (f2[cw] + f2[cw + N]);
                REAL f2s = (REAL)0.5 * (f2[cs] + f2[cs + N]);
                U2[c] = U2[c] + dtau * (FU[c] - (cu[c] * ((pd[c] - pd[cw]) / dx)
                                                 + mu_u[h + i] * ((f2m - f2w) / dx)));
                V2[c] = V2[c] + dtau * (FV[c] - (cv[c] * ((pd[c] - pd[cs]) / dy)
                                                 + mu_v[h + i] * ((f2m - f2s) / dy)));
            }
        }
}

/*
 * Acoustic substep, part 3: everything else, one column at a time --
 * divergence, mu, omega, the theta flux, phi (explicit part), the
 * vertically implicit w solve (Thomas), and phi from the new w.
 * Reads U2, V2 at horizontal neighbours (final after part 2); writes only its
 * own column of m2, T2, W2, f2.
 */
EXPORT void FN(ac_col)(int nz, int ny, int nx, const int *xp1, const int *yp1,
                       REAL dtau, REAL dx, REAL dy, REAL G, REAL am,
                       const REAL *ds, const REAL *dsw, const REAL *Gk,
                       const REAL *dmuF, const REAL *FTh, const REAL *Fphi, const REAL *FW,
                       const REAL *thw, const REAL *thu, const REAL *thv,
                       const REAL *Qc, const REAL *E,
                       const REAL *dphidx, const REAL *dphidy, const REAL *dphi_ds,
                       const REAL *mus, const REAL *Bv,
                       const REAL *lower, const REAL *diag, const REAL *upper,
                       const REAL *pd, const REAL *U2, const REAL *V2,
                       REAL *m2, REAL *T2, REAL *W2, REAL *f2)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel
    {
        REAL divV[nz], Om[nz + 1], Phi[nz + 1], r[nz], cp[nz], dp[nz], Uc[nz], Vc[nz];
        #pragma omp for collapse(2) schedule(static)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const size_t c2 = (size_t)j * nx + i;
                const size_t ce = (size_t)j * nx + xp1[i];
                const size_t cn = (size_t)yp1[j] * nx + i;
                /* divergence and its column integral (sequential sum, as NumPy) */
                REAL dsum = 0;
                for (int k = 0; k < nz; k++) {
                    const size_t o = (size_t)k * N;
                    REAL dv = (U2[o + ce] - U2[o + c2]) / dx + (V2[o + cn] - V2[o + c2]) / dy;
                    divV[k] = dv;
                    dsum = dsum + dv * ds[k];
                    Uc[k] = (REAL)0.5 * (U2[o + c2] + U2[o + ce]);
                    Vc[k] = (REAL)0.5 * (V2[o + c2] + V2[o + cn]);
                }
                const REAL dmu2 = -dsum;
                m2[c2] = m2[c2] + dtau * (dmuF[c2] + dmu2);
                const REAL m2n = m2[c2];
                Om[0] = 0;
                for (int k = 0; k < nz; k++) Om[k + 1] = Om[k] - (dmu2 + divV[k]) * ds[k];
                Om[nz] = 0;
                /* theta: flux form, horizontal from face fluxes, vertical linear faces */
                for (int k = 0; k < nz; k++) {
                    const size_t o = (size_t)k * N, c = o + c2;
                    REAL fx = (U2[o + ce] * thu[o + ce] - U2[c] * thu[c]) / dx;
                    REAL fy = (V2[o + cn] * thv[o + cn] - V2[c] * thv[c]) / dy;
                    REAL Fz1 = Om[k + 1] * thw[(size_t)(k + 1) * N + c2];
                    REAL Fz0 = Om[k] * thw[o + c2];
                    T2[c] = T2[c] + dtau * (FTh[c] - (fx + fy) - (Fz1 - Fz0) / ds[k]);
                }
                /* phi, explicit part, on w levels */
                const REAL mu = mus[c2];
                for (int k = 0; k <= nz; k++) {
                    const size_t cw = (size_t)k * N + c2;
                    REAL uw = (k == 0) ? Uc[0] : (k == nz) ? Uc[nz - 1] : (REAL)0.5 * (Uc[k] + Uc[k - 1]);
                    REAL vw = (k == 0) ? Vc[0] : (k == nz) ? Vc[nz - 1] : (REAL)0.5 * (Vc[k] + Vc[k - 1]);
                    Phi[k] = f2[cw] + dtau * (Fphi[cw] - (uw * dphidx[cw] + vw * dphidy[cw] + Om[k] * dphi_ds[cw]) / mu
                                              + am * G * W2[cw] / mu);
                }
                Phi[nz] = 0;
                /* right-hand side of the w equation */
                REAL Qprev = 0, Eprev = 0;
                for (int k = 0; k < nz; k++) {
                    const size_t c = (size_t)k * N + c2;
                    REAL dpo = (k == 0) ? pd[c2] / dsw[0] : (pd[c] - pd[c - N]) / dsw[k];
                    REAL Wexp = W2[c] + dtau * (FW[c] + G * (am * dpo - m2n));
                    REAL Q = Qc[c] * T2[c];
                    REAL S = Q - E[c] * (Phi[k] - Phi[k + 1]);
                    if (k > 0) S = S - (Qprev - Eprev * (Phi[k - 1] - Phi[k]));
                    Qprev = Q; Eprev = E[c];
                    r[k] = Wexp + Gk[k] * S;
                }
                /* Thomas, as nh3d.thomas */
                cp[0] = upper[c2] / diag[c2];
                dp[0] = r[0] / diag[c2];
                for (int k = 1; k < nz; k++) {
                    const size_t c = (size_t)k * N + c2;
                    REAL m = diag[c] - lower[c] * cp[k - 1];
                    cp[k] = upper[c] / m;
                    dp[k] = (r[k] - lower[c] * dp[k - 1]) / m;
                }
                REAL x = dp[nz - 1];
                W2[(size_t)(nz - 1) * N + c2] = x;
                for (int k = nz - 2; k >= 0; k--) {
                    x = dp[k] - cp[k] * x;
                    W2[(size_t)k * N + c2] = x;
                }
                W2[(size_t)nz * N + c2] = 0;
                const REAL bv = Bv[c2];
                for (int k = 0; k < nz; k++) {
                    const size_t c = (size_t)k * N + c2;
                    f2[c] = Phi[k] + bv * W2[c];
                }
                f2[(size_t)nz * N + c2] = 0;
            }
    }
}


/* ======================================================================
 * TENDENCIES (nh3d.tendencies), in three passes. Each pass reads the
 * previous pass's arrays at horizontal neighbours, so a pass must finish
 * before the next starts.
 * ==================================================================== */

/*
 * Pass 1, one column at a time: diagnose, the face winds, the divergence
 * integral, Omega, and the column quantities the pressure gradient needs.
 *   out 3-D (nz):   pp, dpp, phim, phipm, A1, A2, u, v, th, sdm
 *   out 3-D (nz+1): Om, sd, w
 *   out 2-D:        dmu, mup
 */
EXPORT void FN(td_col)(int nz, int ny, int nx, const int *xm1, const int *xp1,
                       const int *ym1, const int *yp1, REAL dx, REAL dy,
                       REAL P0_, REAL RD_, REAL GAMMA_,
                       const REAL *ds, const REAL *sf,
                       const REAL *pb, const REAL *alb, const REAL *mub, const REAL *phib,
                       const REAL *mu, const REAL *U, const REAL *V, const REAL *W,
                       const REAL *Th, const REAL *phi,
                       REAL *pp, REAL *dpp, REAL *phim, REAL *phipm, REAL *A1, REAL *A2,
                       REAL *u, REAL *v, REAL *th, REAL *sdm,
                       REAL *Om, REAL *sd, REAL *w, REAL *dmu, REAL *mup)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel
    {
        REAL divV[nz], pw[nz + 1], Omc[nz + 1];
        #pragma omp for collapse(2) schedule(static)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const size_t c2 = (size_t)j * nx + i;
                const size_t cw = (size_t)j * nx + xm1[i], ce = (size_t)j * nx + xp1[i];
                const size_t cs = (size_t)ym1[j] * nx + i, cn = (size_t)yp1[j] * nx + i;
                const REAL m = mu[c2];
                const REAL mu_u = (REAL)0.5 * (m + mu[cw]);
                const REAL mu_v = (REAL)0.5 * (m + mu[cs]);
                const REAL mp = m - mub[c2];
                mup[c2] = mp;
                REAL dsum = 0;
                for (int k = 0; k < nz; k++) {
                    const size_t o = (size_t)k * N, c = o + c2;
                    REAL al = (phi[c] - phi[c + N]) / (m * ds[k]);
                    REAL p = P0_ * pow(RD_ * Th[c] / (m * P0_ * al), GAMMA_);
                    pp[c] = p - pb[c];
                    REAL alp = al - alb[c];
                    A1[c] = m * al;
                    A2[c] = m * alp + mp * alb[c];
                    phim[c] = (REAL)0.5 * (phi[c] + phi[c + N]);
                    phipm[c] = (REAL)0.5 * ((phi[c] - phib[c]) + (phi[c + N] - phib[c + N]));
                    u[c] = U[c] / mu_u;
                    v[c] = V[c] / mu_v;
                    th[c] = Th[c] / m;
                    REAL dv = (U[o + ce] - U[c]) / dx + (V[o + cn] - V[c]) / dy;
                    divV[k] = dv;
                    dsum = dsum + dv * ds[k];
                }
                const REAL dm = -dsum;
                dmu[c2] = dm;
                Omc[0] = 0;
                for (int k = 0; k < nz; k++) Omc[k + 1] = Omc[k] - (dm + divV[k]) * ds[k];
                Omc[nz] = 0;
                for (int k = 0; k <= nz; k++) {
                    const size_t c = (size_t)k * N + c2;
                    Om[c] = Omc[k];
                    sd[c] = Omc[k] / m;
                    w[c] = W[c] / m;
                }
                for (int k = 0; k < nz; k++)
                    sdm[(size_t)k * N + c2] = (REAL)0.5 * (Omc[k] / m + Omc[k + 1] / m);
                /* pressure on w levels and its sigma derivative */
                pw[0] = 0;
                for (int k = 1; k < nz; k++)
                    pw[k] = (REAL)0.5 * (pp[(size_t)k * N + c2] + pp[(size_t)(k - 1) * N + c2]);
                {
                    REAL a = pp[(size_t)(nz - 1) * N + c2], b = pp[(size_t)(nz - 2) * N + c2];
                    pw[nz] = a + (a - b) * ((REAL)1.0 - sf[nz - 1]) / (sf[nz - 1] - sf[nz - 2]);
                }
                for (int k = 0; k < nz; k++)
                    dpp[(size_t)k * N + c2] = (pw[k + 1] - pw[k]) / ds[k];
            }
    }
}

/* vertical derivative on mass levels, as dds_m in nh3d.tendencies */
#define DDS_M(a, k, c)                                                              \
    ((k) == 0 ? ((a)[(c) + N] - (a)[(c)]) / (sf[1] - sf[0])                         \
     : (k) == nz - 1 ? ((a)[(c)] - (a)[(c) - N]) / (sf[nz - 1] - sf[nz - 2])        \
     : ((a)[(c) + N] - (a)[(c) - N]) / (sf[(k) + 1] - sf[(k) - 1]))

/* fourth-order centred + upwind dissipation, as nh3d._u3 */
#define U3(vel, a0, p1, m1, p2, m2, d)                                              \
    ((vel) * (((-(p2) + 8 * (p1)) - 8 * (m1)) + (m2)) / (12 * (d))                  \
     + fabs(vel) * ((((p2) - 4 * (p1)) + 6 * (a0)) - 4 * (m1) + (m2)) / (12 * (d)))

/*
 * Pass 2: the momentum tendencies FU, FV (pressure gradient, Coriolis,
 * advection), at every u and v point.
 */
EXPORT void FN(td_uv)(int nz, int ny, int nx,
                      const int *xm1, const int *xp1, const int *xm2, const int *xp2,
                      const int *ym1, const int *yp1, const int *ym2, const int *yp2,
                      REAL dx, REAL dy, const REAL *sf,
                      const REAL *f_u, const REAL *f_v, const REAL *mu, const REAL *mub,
                      const REAL *pb, const REAL *pp, const REAL *dpp, const REAL *phim,
                      const REAL *phipm, const REAL *A1, const REAL *A2,
                      const REAL *U, const REAL *V, const REAL *u, const REAL *v,
                      const REAL *sdm, REAL *FU, REAL *FV)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++) {
            const size_t o = (size_t)k * N;
            const int jm = ym1[j], jp = yp1[j], jm2 = ym2[j], jp2 = yp2[j];
            for (int i = 0; i < nx; i++) {
                const int im = xm1[i], ip = xp1[i], im2 = xm2[i], ip2 = xp2[i];
                const size_t h = (size_t)j * nx + i, hw = (size_t)j * nx + im, hs = (size_t)jm * nx + i;
                const size_t c = o + h, cw = o + hw, cs = o + hs;
                /* ---- u point (j, i-1/2) ---- */
                REAL mu_u = (REAL)0.5 * (mu[h] + mu[hw]);
                REAL pgx = (REAL)0.5 * (A1[c] + A1[cw]) * ((pp[c] - pp[cw]) / dx)
                         + (REAL)0.5 * (A2[c] + A2[cw]) * ((pb[c] - pb[cw]) / dx)
                         + (REAL)0.5 * (mub[h] + mub[hw]) * ((phipm[c] - phipm[cw]) / dx)
                         + (REAL)0.5 * (dpp[c] + dpp[cw]) * ((phim[c] - phim[cw]) / dx);
                /* f * V at the u point: h_to_u(v_to_h(V)) */
                REAL Vh = (REAL)0.5 * (V[c] + V[o + (size_t)jp * nx + i]);
                REAL Vhw = (REAL)0.5 * (V[cw] + V[o + (size_t)jp * nx + im]);
                REAL VatU = (REAL)0.5 * (Vh + Vhw);
                REAL vh = (REAL)0.5 * (v[c] + v[o + (size_t)jp * nx + i]);
                REAL vhw = (REAL)0.5 * (v[cw] + v[o + (size_t)jp * nx + im]);
                REAL v_u = (REAL)0.5 * (vh + vhw);
                REAL ua = u[c];
                REAL advx = U3(ua, ua, u[o + (size_t)j * nx + ip], u[cw], u[o + (size_t)j * nx + ip2],
                               u[o + (size_t)j * nx + im2], dx);
                REAL advy = U3(v_u, ua, u[o + (size_t)jp * nx + i], u[cs], u[o + (size_t)jp2 * nx + i],
                               u[o + (size_t)jm2 * nx + i], dy);
                REAL sdu = (REAL)0.5 * (sdm[c] + sdm[cw]);
                REAL advz = sdu * DDS_M(u, k, c);
                FU[c] = (-pgx + f_u[h] * VatU) - mu_u * ((advx + advy) + advz);
                /* ---- v point (j-1/2, i) ---- */
                REAL mu_v = (REAL)0.5 * (mu[h] + mu[hs]);
                REAL pgy = (REAL)0.5 * (A1[c] + A1[cs]) * ((pp[c] - pp[cs]) / dy)
                         + (REAL)0.5 * (A2[c] + A2[cs]) * ((pb[c] - pb[cs]) / dy)
                         + (REAL)0.5 * (mub[h] + mub[hs]) * ((phipm[c] - phipm[cs]) / dy)
                         + (REAL)0.5 * (dpp[c] + dpp[cs]) * ((phim[c] - phim[cs]) / dy);
                REAL Uh = (REAL)0.5 * (U[c] + U[o + (size_t)j * nx + ip]);
                REAL Uhs = (REAL)0.5 * (U[cs] + U[o + (size_t)jm * nx + ip]);
                REAL UatV = (REAL)0.5 * (Uh + Uhs);
                REAL uh = (REAL)0.5 * (u[c] + u[o + (size_t)j * nx + ip]);
                REAL uhs = (REAL)0.5 * (u[cs] + u[o + (size_t)jm * nx + ip]);
                REAL u_v = (REAL)0.5 * (uh + uhs);
                REAL va = v[c];
                REAL bdx = U3(u_v, va, v[o + (size_t)j * nx + ip], v[cw], v[o + (size_t)j * nx + ip2],
                              v[o + (size_t)j * nx + im2], dx);
                REAL bdy = U3(va, va, v[o + (size_t)jp * nx + i], v[cs], v[o + (size_t)jp2 * nx + i],
                              v[o + (size_t)jm2 * nx + i], dy);
                REAL sdv = (REAL)0.5 * (sdm[c] + sdm[cs]);
                REAL bdz = sdv * DDS_M(v, k, c);
                FV[c] = (-pgy - f_v[h] * UatV) - mu_v * ((bdx + bdy) + bdz);
            }
        }
}

/* third-order upwind-biased face flux F q at the west/south face, as nh3d._flux3 */
#define FLUX3(F, q0, qm1, qp1, qm2)                                                 \
    ((F) * ((REAL)7.0 * ((q0) + (qm1)) - ((qp1) + (qm2))) / (REAL)12.0              \
     - fabs(F) * ((REAL)3.0 * ((q0) - (qm1)) - ((qp1) - (qm2))) / (REAL)12.0)

/*
 * Pass 3, one column at a time: the theta flux divergence (third-order
 * upwind-biased faces), and the W and phi tendencies.
 */
EXPORT void FN(td_tw)(int nz, int ny, int nx,
                      const int *xm1, const int *xp1, const int *xm2, const int *xp2,
                      const int *ym1, const int *yp1, const int *ym2, const int *yp2,
                      REAL dx, REAL dy, REAL G_, const REAL *ds, const REAL *dsw, const REAL *sh,
                      const REAL *mu, const REAL *U, const REAL *V, const REAL *th,
                      const REAL *Om, const REAL *pp, const REAL *mup,
                      const REAL *u, const REAL *v, const REAL *w, const REAL *sd, const REAL *phi,
                      REAL *FTh, REAL *FW, REAL *Fphi)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel
    {
        REAL Fz[nz + 1], uc[nz], vc[nz];
        #pragma omp for collapse(2) schedule(static)
        for (int j = 0; j < ny; j++)
            for (int i = 0; i < nx; i++) {
                const int im = xm1[i], ip = xp1[i];
                const int jm = ym1[j], jp = yp1[j];
                const size_t c2 = (size_t)j * nx + i;
                const REAL m = mu[c2];
                /* vertical theta flux on w levels (_vface3) */
                Fz[0] = Om[c2] * th[c2];
                Fz[nz] = Om[(size_t)nz * N + c2] * th[(size_t)(nz - 1) * N + c2];
                for (int k = 1; k < nz; k++) {
                    const size_t c = (size_t)k * N + c2;
                    REAL qf;
                    if (nz >= 4 && k >= 2 && k <= nz - 2) {
                        REAL a = th[c - 2 * N], b = th[c - N], cc = th[c], d = th[c + N];
                        REAL avg = ((REAL)7.0 * (b + cc) - (a + d)) / (REAL)12.0;
                        REAL dif = ((REAL)3.0 * (cc - b) - (d - a)) / (REAL)12.0;
                        REAL o = Om[c];
                        REAL sg = (o > 0) ? (REAL)1 : (o < 0) ? (REAL)-1 : (REAL)0;
                        qf = avg - sg * dif;
                    } else {
                        qf = (REAL)0.5 * (th[c] + th[c - N]);
                    }
                    Fz[k] = Om[c] * qf;
                }
                for (int k = 0; k < nz; k++) {
                    const size_t o = (size_t)k * N, c = o + c2;
                    const size_t rowj = o + (size_t)j * nx, rowjp = o + (size_t)jp * nx;
                    /* x faces at i (west) and ip (east, = shifted flux) */
                    REAL fw = FLUX3(U[c], th[c], th[rowj + im], th[rowj + ip], th[rowj + xm2[i]]);
                    const int e = ip;
                    REAL fe = FLUX3(U[rowj + e], th[rowj + e], th[rowj + xm1[e]], th[rowj + xp1[e]],
                                    th[rowj + xm2[e]]);
                    /* y faces at j (south) and jp (north) */
                    const size_t rjm = o + (size_t)jm * nx, rjm2 = o + (size_t)ym2[j] * nx;
                    REAL fs = FLUX3(V[c], th[c], th[rjm + i], th[rowjp + i], th[rjm2 + i]);
                    const int nj = jp;
                    const size_t rn = o + (size_t)nj * nx, rnm = o + (size_t)ym1[nj] * nx;
                    const size_t rnp = o + (size_t)yp1[nj] * nx, rnm2 = o + (size_t)ym2[nj] * nx;
                    REAL fn = FLUX3(V[rn + i], th[rn + i], th[rnm + i], th[rnp + i], th[rnm2 + i]);
                    FTh[c] = -((fe - fw) / dx + (fn - fs) / dy) - (Fz[k + 1] - Fz[k]) / ds[k];
                    uc[k] = (REAL)0.5 * (u[c] + u[rowj + ip]);
                    vc[k] = (REAL)0.5 * (v[c] + v[rowjp + i]);
                }
                /* W and phi on w levels */
                const REAL mp = mup[c2];
                for (int k = 0; k <= nz; k++) {
                    const size_t c = (size_t)k * N + c2;
                    const size_t rowj = (size_t)k * N + (size_t)j * nx;
                    REAL uw = (k == 0) ? uc[0] : (k == nz) ? uc[nz - 1] : (REAL)0.5 * (uc[k] + uc[k - 1]);
                    REAL vw = (k == 0) ? vc[0] : (k == nz) ? vc[nz - 1] : (REAL)0.5 * (vc[k] + vc[k - 1]);
                    REAL dwds = (k == 0 || k == nz) ? (REAL)0
                              : (w[c + N] - w[c - N]) / (sh[k + 1] - sh[k - 1]);
                    REAL dphids = (k == nz) ? (REAL)0
                                : (k == 0) ? (phi[c + N] - phi[c]) / (sh[1] - sh[0])
                                : (phi[c + N] - phi[c - N]) / (sh[k + 1] - sh[k - 1]);
                    const size_t ce = rowj + ip, cw = rowj + im;
                    const size_t cn = (size_t)k * N + (size_t)jp * nx + i, cs = (size_t)k * N + (size_t)jm * nx + i;
                    if (k < nz) {
                        REAL fw0 = (k == 0) ? G_ * (pp[c2] / dsw[0] - mp)
                                            : G_ * ((pp[c] - pp[c - N]) / dsw[k] - mp);
                        REAL advw = uw * (w[ce] - w[cw]) / (2 * dx)
                                  + vw * (w[cn] - w[cs]) / (2 * dy) + sd[c] * dwds;
                        FW[c] = fw0 - m * advw;
                        Fphi[c] = -(uw * (phi[ce] - phi[cw]) / (2 * dx)
                                    + vw * (phi[cn] - phi[cs]) / (2 * dy) + sd[c] * dphids) + G_ * w[c];
                    } else {
                        FW[c] = 0;
                        Fphi[c] = 0;
                    }
                }
            }
    }
}


/* ======================================================================
 * SET-UP OF THE ACOUSTIC SUBSTEPS (the part of nh3d._acoustic before its
 * loop), in two passes, and the final sum after it.
 * ==================================================================== */

/* Pass A, pointwise: al, Qc, E and theta at the stage state. */
EXPORT void FN(as_a)(int nz, int ny, int nx, REAL P0_, REAL RD_, REAL GAMMA_,
                     const REAL *ds, const REAL *mus, const REAL *Ths, const REAL *phis,
                     REAL *al, REAL *Qc, REAL *E, REAL *ths)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int k = 0; k < nz; k++)
        for (int j = 0; j < ny; j++) {
            const size_t o = (size_t)k * N + (size_t)j * nx;
            const size_t h = (size_t)j * nx;
            for (int i = 0; i < nx; i++) {
                const size_t c = o + i;
                const REAL m = mus[h + i];
                REAL a = (phis[c] - phis[c + N]) / (m * ds[k]);
                REAL p = P0_ * pow(RD_ * Ths[c] / (m * P0_ * a), GAMMA_);
                REAL c2 = GAMMA_ * p * a;
                al[c] = a;
                ths[c] = Ths[c] / m;
                Qc[c] = c2 / a / Ths[c];
                E[c] = c2 / (a * a * m * ds[k]);
            }
        }
}

/* Pass B, one column at a time: every coefficient that is fixed during the
 * substeps, the tridiagonal matrix, and the starting perturbations. */
EXPORT void FN(as_b)(int nz, int ny, int nx, const int *xm1, const int *xp1,
                     const int *ym1, const int *yp1, REAL dx, REAL dy, REAL bvc,
                     const REAL *sh, const REAL *Gk,
                     const REAL *mus, const REAL *al, const REAL *ths, const REAL *phis,
                     const REAL *Us, const REAL *Vs, const REAL *Ws, const REAL *Ths,
                     const REAL *mu0, const REAL *U0, const REAL *V0, const REAL *W0,
                     const REAL *Th0, const REAL *phi0, const REAL *Qc, const REAL *E,
                     REAL *cu, REAL *mu_u, REAL *cv, REAL *mu_v, REAL *thu, REAL *thv, REAL *thw,
                     REAL *dphidx, REAL *dphidy, REAL *dphi_ds, REAL *Bv,
                     REAL *lower, REAL *diag, REAL *upper,
                     REAL *m2, REAL *U2, REAL *V2, REAL *W2, REAL *T2, REAL *f2, REAL *p_old)
{
    const size_t N = (size_t)ny * nx;
    #pragma omp parallel for collapse(2) schedule(static)
    for (int j = 0; j < ny; j++)
        for (int i = 0; i < nx; i++) {
            const size_t c2 = (size_t)j * nx + i;
            const size_t hw = (size_t)j * nx + xm1[i], he = (size_t)j * nx + xp1[i];
            const size_t hs = (size_t)ym1[j] * nx + i, hn = (size_t)yp1[j] * nx + i;
            const REAL m = mus[c2];
            const REAL mu_u_ = (REAL)0.5 * (m + mus[hw]);
            const REAL mu_v_ = (REAL)0.5 * (m + mus[hs]);
            mu_u[c2] = mu_u_;
            mu_v[c2] = mu_v_;
            const REAL bv = bvc / m;
            Bv[c2] = bv;
            m2[c2] = mu0[c2] - m;
            for (int k = 0; k < nz; k++) {
                const size_t o = (size_t)k * N, c = o + c2;
                cu[c] = mu_u_ * ((REAL)0.5 * (al[c] + al[o + hw]));
                cv[c] = mu_v_ * ((REAL)0.5 * (al[c] + al[o + hs]));
                thu[c] = (REAL)0.5 * (ths[c] + ths[o + hw]);
                thv[c] = (REAL)0.5 * (ths[c] + ths[o + hs]);
                U2[c] = U0[c] - Us[c];
                V2[c] = V0[c] - Vs[c];
                T2[c] = Th0[c] - Ths[c];
                REAL d = (REAL)1.0 + Gk[k] * bv * E[c];
                if (k > 0) {
                    d = d + Gk[k] * bv * E[c - N];
                    lower[c] = -Gk[k] * bv * E[c - N];
                } else {
                    lower[c] = 0;
                }
                diag[c] = d;
                upper[c] = (k < nz - 1) ? -Gk[k] * bv * E[c] : (REAL)0;
            }
            for (int k = 0; k <= nz; k++) {
                const size_t o = (size_t)k * N, c = o + c2;
                thw[c] = (k == 0) ? ths[c2] : (k == nz) ? ths[(size_t)(nz - 1) * N + c2]
                       : (REAL)0.5 * (ths[c] + ths[c - N]);
                dphidx[c] = (phis[o + he] - phis[o + hw]) / (2 * dx);
                dphidy[c] = (phis[o + hn] - phis[o + hs]) / (2 * dy);
                dphi_ds[c] = (k == nz) ? (REAL)0
                           : (k == 0) ? (phis[c + N] - phis[c]) / (sh[1] - sh[0])
                           : (phis[c + N] - phis[c - N]) / (sh[k + 1] - sh[k - 1]);
                W2[c] = W0[c] - Ws[c];
                f2[c] = phi0[c] - phis[c];
            }
            for (int k = 0; k < nz; k++) {
                const size_t c = (size_t)k * N + c2;
                p_old[c] = Qc[c] * T2[c] - E[c] * (f2[c] - f2[c + N]);
            }
        }
}

/* out = a + b over n elements (the stage state plus its perturbation) */
EXPORT void FN(add_into)(size_t n, const REAL *a, const REAL *b, REAL *out)
{
    #pragma omp parallel for schedule(static)
    for (size_t c = 0; c < n; c++) out[c] = a[c] + b[c];
}

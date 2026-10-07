
"""Tendency budget at the 28 Sep northern-Vermont hot spot (P-67), from saved hourly states.
usage: python budget28.py FORECAST_NPZ Z0   (run from the dev tree root)"""
import sys, numpy as np
sys.path[:0] = ["src", "src/dynamics", "src/verification", "src/analysis", "."]
import forecast, primitive_sigma as PS
from backend import to_numpy
from surface import surface_drag
fc_path, z0 = sys.argv[1], sys.argv[2]
D = "/data5/pierce/AINWP/data/tensors_3d/obs_20260928_06"
J0, J1, I0, I1 = 70, 80, 58, 72          # box around 45.2-45.4N, 72.2-72.8W
def budget(model, driver, relax, duration, **kw):
    z = np.load(fc_path); hrs = np.rint(z["times_s"] / 3600).astype(int)
    gr, lev = model.grid, model.lev
    for n, h in enumerate(hrs):
        if h < 12: continue
        u, v, th, pi = (z[k][n].astype(float) for k in ("u", "v", "theta", "pi"))
        full = model.tendencies(u, v, th, pi)
        du, dv = full[0], full[1]
        phi = PS.hydrostatic_geopotential(th, pi, lev, phi_surface=model.phi_s)
        T_ref = np.mean(th * (lev.pressure(pi) / PS.P0) ** PS.KAPPA, axis=(1, 2))
        fx, fy = PS.pressure_gradient_force(phi, th, pi, lev, gr, reference=T_ref if model.ref_pgf else None)
        v_at_u = gr.v_to_u(v); u_at_v = gr.u_to_v(u)
        cu, cv = np.asarray(gr.f_u) * v_at_u, -np.asarray(gr.f_v) * u_at_v
        hu = PS.hyperdiffusion(u, gr, model.hyper) if model.hyper > 0 else 0 * u
        hv = PS.hyperdiffusion(v, gr, model.hyper) if model.hyper > 0 else 0 * v
        ddu, ddv = PS.divergence_damping(u, v, gr, model.div_damp) if model.div_damp > 0 else (0 * u, 0 * v)
        gu, gv, _ = surface_drag(u, v, th, pi, lev, z0=model.z0)
        mu, mv, _, K = PS.vertical_mixing(u, v, th, pi, lev, ri_crit=model.ri_crit, k_max=model.k_max, mixing_length=model.mixing_length)
        dpi_dt, sd = PS.continuity(u, v, pi, lev, gr)
        hau = -model._horiz_adv(u, u, v_at_u); hav = -model._horiz_adv(v, u_at_v, v)
        vau = -PS.vertical_advection(u, sd, lev); vav = -PS.vertical_advection(v, sd, lev)
        spd = np.hypot(u, v)
        sub = spd[:, J0:J1, I0:I1]; k, j, i = np.unravel_index(np.argmax(sub), sub.shape); j += J0; i += I0
        terms = {"total": (du, dv), "pgf": (fx, fy), "coriolis": (cu, cv), "hyper": (hu, hv), "divdamp": (ddu, ddv),
                 "drag": (gu, gv), "mixing": (mu, mv), "hadv": (hau, hav), "vadv": (vau, vav)}
        ax, ay = u[k, j, i] / spd[k, j, i], v[k, j, i] / spd[k, j, i]
        out = {nm: (a[k, j, i] * ax + b[k, j, i] * ay) * 3600 for nm, (a, b) in terms.items()}
        out["residual"] = out["total"] - sum(out[x] for x in ("pgf", "coriolis", "hyper", "divdamp", "drag", "mixing", "hadv", "vadv"))
        # speed of the 3x3x3 neighbourhood: is (k,j,i) a strict local maximum of speed?
        nb = spd[max(k-1,0):k+2, j-1:j+2, i-1:i+2]
        out["nbr_max-V"] = float(nb.max() - spd[k, j, i]) * 3600 / 3600
        # pgf + coriolis = ageostrophic acceleration
        print(f"{h:2d} h  level {k} j {j} i {i}  |V| {spd[k,j,i]:5.1f} m/s  along-wind tendency (m/s per h): "
              + "  ".join(f"{nm} {val:+6.2f}" for nm, val in out.items())
              + f"  | K {float(K[min(k,K.shape[0]-1), j, i]) if hasattr(K,'shape') else K}")
        # cross-wind (to the left) PGF+Coriolis
        lx, ly = -ay, ax
        print(f"      cross-wind (left +): pgf {(fx[k,j,i]*lx+fy[k,j,i]*ly)*3600:+6.2f}  coriolis {(cu[k,j,i]*lx+cv[k,j,i]*ly)*3600:+6.2f}"
              f"   terrain {model.phi_s[j,i]/9.81 if model.phi_s is not None else 0:.0f} m, slope to E {((model.phi_s[j,i+1]-model.phi_s[j,i-1])/9.81/(2*gr.dx)) if model.phi_s is not None else 0:+.4f}, to N {((model.phi_s[j+1,i]-model.phi_s[j-1,i])/9.81/(2*gr.dy)) if model.phi_s is not None else 0:+.4f}")
    info = kw.get("info") if kw.get("info") is not None else {}
    info["stopped"] = "budget"
    return []
forecast.run_forecast = budget
sys.argv = ["forecast.py", "--run-dir", D, "--hours", "24", "--z0", z0, "--backend", "numpy", "--out", "/tmp/_budget_unused.npz"]
forecast.main()

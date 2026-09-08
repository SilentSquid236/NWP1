"""
P-02: can the sponge absorb gravity waves without damping the flow?

WHAT IS MEASURED, AND WHY THE OBVIOUS FIX IS NOT AVAILABLE

Growth over terrain peaks at exactly the sponge's lower edge and moves when
the edge moves -- partial reflection off the absorbing layer:

    sponge levels |  0   |  5   |  8   | 12
    peak growth k |  0   |  5   |  8   | 18
    max|du| (6 h) | 60.8 | 36.4 | 21.3 | 15.1

Deepening the sponge halves the reflection. But a 12-level sponge on a
20-level model damps the free troposphere, and an earlier version that relaxed
toward the horizontal mean was caught by the thermal-wind test for flattening
a jet (P-16). So the sponge has to stay shallow, and the reflection stays.

THE HYPOTHESIS, STATED BEFORE THE RUN

Gravity waves are DIVERGENT. Balanced flow -- the jet, the thing that must not
be damped -- is ROTATIONAL. Helmholtz splits them exactly, and `subgrid.
remove_divergence_spectral` already does that split with the eigenvalues of
our discrete Laplacian.

So: damp only the divergent component in the sponge. The prediction is that
this lets the layer be DEEP without flattening anything, because there is
almost nothing divergent in a balanced jet to damp.

    1. reflection amplitude falls at least as much as a deep plain sponge
    2. the thermal-wind jet survives, which a deep plain sponge should not
    3. survival over terrain is no worse

Prediction 2 is the one worth stating. A scheme that quietly damps the jet
would look identical on 1 and 3, and that is exactly how P-16 failed.
"""
import numpy as np
np.seterr(all="ignore")

from grid import CGrid
from sigma import SigmaLevels, RD, G0, P0, KAPPA
from primitive_sigma import PrimitiveSigma
from subgrid import remove_divergence_spectral
from initialization import filter_initial_state
from subgrid import balance_initial_state
from lid_test import build_on


def make_divergent_sponge(m):
    """
    Replace the model's sponge with one that damps only the divergent wind.

    Applied as a post-step operator: the rotational part is left untouched,
    and only the divergent residual inside the layer is relaxed away.
    """
    prof = m._sponge[:, 0, 0].copy()          # per-level rate, 0 below the layer

    orig_step = m.step

    def step(dt):
        # Turn the built-in sponge off for this step and do it ourselves.
        saved = m._sponge.copy()
        m._sponge = np.zeros_like(m._sponge)
        orig_step(dt)
        m._sponge = saved

        active = prof > 0
        if not active.any():
            return
        u_rot = np.empty_like(m.u)
        v_rot = np.empty_like(m.v)
        for k in np.nonzero(active)[0]:
            u_rot[k], v_rot[k] = remove_divergence_spectral(m.u[k], m.v[k],
                                                            m.grid)
            a = prof[k] * dt
            # Relax the DIVERGENT residual toward zero; the rotational part
            # passes through untouched.
            m.u[k] = m.u[k] - a * (m.u[k] - u_rot[k])
            m.v[k] = m.v[k] - a * (m.v[k] - v_rot[k])

    m.step = step
    return m


def thermal_wind_case(sponge_levels, divergent=False):
    """The test that caught P-16: a balanced jet that must persist."""
    gr = CGrid(48, 48, 50e3, 50e3, f0=1.0e-4, beta=0.0)
    lev = SigmaLevels(20)
    m = PrimitiveSigma(gr, lev, sponge_levels=sponge_levels)
    T0, dT = 260.0, 4.0
    k_y = 2 * np.pi / gr.Ly
    m.pi = np.full((gr.ny, gr.nx), 101325.0) - lev.p_top
    p = lev.pressure(m.pi)
    T = T0 + dT * np.cos(k_y * gr.Yc)
    m.theta = T / (p / P0) ** KAPPA
    phi = m.geopotential()
    for k in range(lev.nz):
        dphidy = 0.5 * (gr.dy_forward(phi[k]) + gr.dy_backward(phi[k]))
        m.u[k] = -dphidy / gr.f0
    if divergent:
        make_divergent_sponge(m)
    u0 = m.u.copy()
    m.run(24 * 3600, dt=m.max_dt())
    ok = np.isfinite(m.u).all()
    drift = (float(np.abs(m.u - u0).max() / np.abs(u0).max()) if ok
             else float("nan"))
    return drift, float(np.abs(u0).max())


def terrain_case(sponge_levels, divergent=False, hours=6):
    m = build_on(SigmaLevels(20), 2500.0, sponge_levels)
    if divergent:
        make_divergent_sponge(m)
    u0 = m.u.copy()
    dt = m.max_dt()
    done = 0
    for _ in range(hours):
        m.run(3600, dt=dt)
        if not np.isfinite(m.u).all():
            break
        done += 1
    if not np.isfinite(m.u).all():
        return float("nan"), -1, done
    d = m.u - u0
    amp = np.sqrt((d ** 2).mean(axis=(1, 2)))
    return float(np.abs(d).max()), int(np.argmax(amp)), done


if __name__ == "__main__":
    print("P-02: sponge reflection vs jet preservation\n")
    print(f"{'sponge':>7} {'kind':>11} {'jet drift 24h':>14} "
          f"{'refl max|du| 6h':>16} {'peak k':>7}")
    for nsp, kind in ((5, "plain"), (8, "plain"), (12, "plain"),
                      (8, "divergent"), (12, "divergent")):
        div = kind == "divergent"
        drift, jet = thermal_wind_case(nsp, div)
        refl, peak, done = terrain_case(nsp, div)
        print(f"{nsp:7d} {kind:>11} {drift*100:13.2f}% "
              f"{refl:16.2f} {peak:7d}", flush=True)

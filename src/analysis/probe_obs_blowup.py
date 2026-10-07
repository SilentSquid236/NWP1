"""
P-56: WHERE does the first observation-built forecast die?

    python src/analysis/probe_obs_blowup.py <run_dir> [hours]

2026-09-21 12Z, built from observations only, ran 3 h at max|u| 46 m/s and
then reached 372 m/s before hour 4 (stopped by the in-hour guard). Before
anything is changed, locate it (nwp-debug, step 2): record every prognostic
field every 5 minutes and report, for the first 5-minute interval in which
max|u| or max|v| grows by more than 20 %,

  * where: grid index, distance from the nearest edge in cells
  * which level (index 0 is the LID)
  * which field ran away first: u, v, theta, pi
  * the local terrain height and the dominant scale of the anomaly

Same setup as forecast.main(), step for step -- a probe that prepares the
state differently from the forecast would be probing a different run.
"""

import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
SRC = HERE.parent
for p in (SRC, SRC.parent, SRC / "dynamics", SRC / "verification"):
    sys.path.insert(0, str(p))

import numpy as np

import forecast as F
import config
from vertical import PressureLevels
from sigma import SigmaLevels
from primitive_sigma import PrimitiveSigma
from boundaries import BoundaryDriver
from initialization import filter_initial_state
from subgrid import balance_initial_state


def main():
    run_dir = Path(sys.argv[1])
    hours = float(sys.argv[2]) if len(sys.argv) > 2 else 4.0
    files, source = F.driving_frames(run_dir)
    lev = SigmaLevels(config.N_LEVELS)
    fields0, meta0 = F.load_state(files[0])
    grid = F.build_grid(fields0, PressureLevels(config.PRESSURE_LEVELS))
    terrain, p_sfc = F.load_terrain(run_dir, fields0.shape[-2:])
    pi0, u0, v0, th0 = F.hrrr_to_sigma_state(fields0, lev, terrain, p_surface=p_sfc)
    pb, ub, vb, tb = pi0.copy(), u0.copy(), v0.copy(), th0.copy()
    ub, vb, tb = filter_initial_state(ub, vb, tb, grid)
    ub, vb, _ = balance_initial_state(ub, vb, grid, verbose=False)
    driver = BoundaryDriver([0.0], [F.state_to_boundary(ub, vb, tb, pb)])
    model = PrimitiveSigma(grid, lev, terrain=terrain)
    u0, v0, th0 = filter_initial_state(u0, v0, th0, grid)
    u0, v0, _ = balance_initial_state(u0, v0, grid)
    model.pi = pi0.copy()
    model.u, model.v, model.theta = u0.copy(), v0.copy(), th0.copy()
    relax = F.Relaxation3D(grid, width=10)

    info = {}
    snaps = F.run_forecast(model, driver, relax, hours * 3600, output_every=300,
                           progress=False, info=info)
    print(f"\n{len(snaps)} snapshots, {info['stopped']}")
    t = np.array([s[0] for s in snaps]) / 3600
    umax = np.array([np.nanmax(np.abs(s[1])) for s in snaps])
    vmax = np.array([np.nanmax(np.abs(s[2])) for s in snaps])
    print("  t(h)   max|u|  max|v|  theta range      pi range (hPa)")
    for i in range(0, len(snaps), max(1, len(snaps) // 30)):
        th, pi = snaps[i][3], snaps[i][4]
        print(f"  {t[i]:5.2f} {umax[i]:7.1f} {vmax[i]:7.1f}  "
              f"{np.nanmin(th):6.1f}-{np.nanmax(th):6.1f}  "
              f"{np.nanmin(pi)/100:6.1f}-{np.nanmax(pi)/100:6.1f}")

    grow = np.where((umax[1:] > 1.2 * umax[:-1]) | (vmax[1:] > 1.2 * vmax[:-1]))[0]
    if grow.size == 0:
        print("\nno 5-minute interval grew by more than 20 %: the runaway is "
              "faster than 5 minutes or not in u, v -- widen the net")
        return 0
    i = int(grow[0]) + 1
    print(f"\nfirst >20 % growth between t={t[i-1]:.2f} h and t={t[i]:.2f} h")
    ny, nx = terrain.shape
    for name, idx in (("u", 1), ("v", 2)):
        a = np.abs(snaps[i][idx])
        k, j, ii = np.unravel_index(np.nanargmax(a), a.shape)
        edge = min(j, ii, ny - 1 - j, nx - 1 - ii)
        prev = np.abs(snaps[i - 1][idx])
        print(f"  |{name}| max {a[k, j, ii]:.1f} (was {prev[k, j, ii]:.1f}) at level {k} "
              f"(0 = lid), row {j}, col {ii}, {edge} cells from the edge, "
              f"terrain {terrain[j, ii]:.0f} m")
        # Where is the growth spread? Count points that grew > 5 m/s.
        d = a - prev
        big = d > 5.0
        if big.any():
            kk, jj, xx = np.where(big)
            e = np.minimum.reduce([jj, xx, ny - 1 - jj, nx - 1 - xx])
            print(f"      {big.sum()} points grew >5 m/s; levels {sorted(set(kk.tolist()))[:10]}; "
                  f"edge distance min {e.min()} median {int(np.median(e))} cells")
            row = d[k, j]
            s = np.sign(row - row.mean())
            crossings = int((np.diff(s) != 0).sum())
            print(f"      along row {j}: {crossings} sign changes over {nx} cells "
                  f"(~{2 * nx / max(crossings, 1):.1f} dx wavelength)")
    dth = snaps[i][3] - snaps[i - 1][3]
    dpi = snaps[i][4] - snaps[i - 1][4]
    print(f"  theta change max {np.nanmax(np.abs(dth)):.2f} K; "
          f"pi change max {np.nanmax(np.abs(dpi))/100:.2f} hPa")
    np.savez_compressed(run_dir / "probe_blowup.npz", t=t, umax=umax, vmax=vmax,
                        u_last=snaps[i][1], u_prev=snaps[i - 1][1],
                        v_last=snaps[i][2], v_prev=snaps[i - 1][2], terrain=terrain)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

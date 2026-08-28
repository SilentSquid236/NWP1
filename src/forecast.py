"""
End-to-end forecast driver: initialise from HRRR, integrate our own physics,
verify against observations, archive the result.

    python src/forecast.py --start 2026-08-01T00 --hours 24

WHAT DEPENDS ON HRRR, AND WHAT DOES NOT

  initial conditions   HRRR, at cold start. Once DA cycling exists this
                       becomes our own previous forecast corrected by
                       observations, and HRRR is only a first guess.
  boundary conditions  HRRR, permanently. A bounded domain must be told what
                       arrives at its edges every timestep. Unavoidable
                       without going global; every operational regional model
                       works this way.
  verification truth   NEVER HRRR. Observations only -- ASOS, mesonets,
                       radiosondes.

The forecast in between is ours: our equations, our numerics, our errors.
"""

import argparse
import sys
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import resources
RESOURCE_PLAN = resources.apply()

import numpy as np

import config
sys.path.insert(0, str(Path(__file__).resolve().parent / "dynamics"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "verification"))

from grid import CGrid
from vertical import PressureLevels, theta_from_T, T_from_theta
from primitive3d import Primitive3D
from boundaries import DaviesRelaxation, BoundaryDriver
from subgrid import StochasticPerturbation


# ---------------------------------------------------------------------------
# HRRR state -> model state
# ---------------------------------------------------------------------------

def load_state(path):
    """Load one ingested HRRR field set (.npz)."""
    z = np.load(path, allow_pickle=False)
    return z["hrrr_features"], {k: z[k] for k in z.files if k != "hrrr_features"}


def hrrr_to_model_state(fields, levels, channels=None):
    """
    Convert an ingested HRRR array [C, L, Y, X] to (u, v, theta).

    HRRR gives temperature; the core is formulated in potential temperature,
    which is the variable conserved by dry adiabatic motion. Converting on the
    way in means the conversion happens once, not every timestep.
    """
    channels = list(channels or config.CHANNELS)
    idx = {c: i for i, c in enumerate(channels)}
    for need in ("TMP", "UGRD", "VGRD"):
        if need not in idx:
            raise KeyError(f"channel {need} missing from {channels}")

    T = fields[idx["TMP"]].astype(float)
    u = fields[idx["UGRD"]].astype(float)
    v = fields[idx["VGRD"]].astype(float)

    p = levels.p.reshape(-1, 1, 1)
    theta = theta_from_T(T, p)
    return u, v, theta


def build_grid(fields, levels, domain=None):
    """
    Construct the model grid to match the ingested field dimensions.

    Grid spacing is derived from the domain extent and array shape rather than
    assumed, so a change to the subset bounds cannot silently desynchronise
    the dynamics from the data.
    """
    domain = domain or config.DOMAIN
    _, _, ny, nx = fields.shape

    lat0 = 0.5 * (domain["lat_min"] + domain["lat_max"])
    dy = (domain["lat_max"] - domain["lat_min"]) * 111_132.0 / ny
    dx = (domain["lon_max"] - domain["lon_min"]) * 111_320.0 * \
        np.cos(np.radians(lat0)) / nx

    # Beta-plane centred on the domain.
    omega = 7.2921e-5
    f0 = 2 * omega * np.sin(np.radians(lat0))
    beta = 2 * omega * np.cos(np.radians(lat0)) / 6_371_000.0

    return CGrid(nx, ny, dx, dy, f0=f0, beta=beta, edge_mode="replicate")


def state_to_boundary(u, v, theta):
    """Package a model state as boundary-driver input."""
    return {"u": u, "v": v, "theta": theta}


# ---------------------------------------------------------------------------
# Boundary relaxation for the 3D core
# ---------------------------------------------------------------------------

class Relaxation3D:
    """
    Davies relaxation applied to a 3D state.

    The 2D weight field is reused at every level: the lateral boundary is a
    vertical wall, so the taper depends only on horizontal distance from the
    edge.
    """

    def __init__(self, grid, width=10, alpha_max=1.0, profile="cosine"):
        self.inner = DaviesRelaxation(grid, width, alpha_max, profile)
        self.alpha = self.inner.alpha[None, :, :]
        self.width = width

    def apply(self, model, ext):
        a = self.alpha
        if "u" in ext:
            model.u += a * (ext["u"] - model.u)
        if "v" in ext:
            model.v += a * (ext["v"] - model.v)
        if "theta" in ext:
            model.theta += a * (ext["theta"] - model.theta)

    @property
    def interior_fraction(self):
        return self.inner.interior_fraction


def run_forecast(model, driver, relax, duration, dt=None, output_every=None,
                 progress=True):
    """
    Integrate with boundary relaxation, collecting output states.

    Returns a list of (valid_seconds, u, v, theta) snapshots.
    """
    dt = dt or model.max_dt()
    n_steps = int(np.ceil(duration / dt))
    dt = duration / n_steps
    interval = output_every or duration

    # Emit on TARGET TIMES, not on a step count. Deriving a stride as
    # int(interval / dt) truncates, so snapshots drift steadily earlier than
    # requested -- an hourly output at dt=771 s would land at 0.86 h, 1.71 h,
    # 2.57 h. Crossing a target time is exact regardless of dt.
    targets = list(np.arange(interval, duration + 1e-9, interval))
    if not targets or targets[-1] < duration - 1e-9:
        targets.append(duration)

    snapshots = []
    next_i = 0
    for k in range(n_steps):
        model.step(dt)
        relax.apply(model, driver.at(model.time))

        if next_i < len(targets) and model.time >= targets[next_i] - 1e-9:
            next_i += 1
            snapshots.append((model.time, model.u.copy(), model.v.copy(),
                              model.theta.copy()))
            if progress:
                d = model.diagnostics()
                print(f"  +{model.time/3600:5.1f} h  "
                      f"max|u| {d['max|u|']:6.1f} m/s  "
                      f"theta {d['theta_min']:.1f}-{d['theta_max']:.1f} K  "
                      f"max|omega| {d['max|omega| Pa/s']:.3f} Pa/s")
            if not np.isfinite(model.u).all():
                print("  FORECAST DIVERGED -- stopping")
                break

    return snapshots


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--run-dir", required=True,
                   help="Directory of ingested live_hrrr_f*.npz files")
    p.add_argument("--hours", type=int, default=12)
    p.add_argument("--output-every", type=float, default=1.0,
                   help="Snapshot interval in hours")
    p.add_argument("--relax-width", type=int, default=10)
    p.add_argument("--stochastic", action="store_true",
                   help="Enable SPPT-style tendency perturbations")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--out", default=None, help="Where to write forecast .npz")
    args = p.parse_args()

    run_dir = Path(args.run_dir)
    files = sorted(run_dir.glob("live_hrrr_f*.npz"),
                   key=lambda q: int(q.stem.split("_f")[-1]))
    if not files:
        print(f"No live_hrrr_f*.npz in {run_dir}. Run ingest_hrrr.py first.")
        return 1

    print("NWP forecast")
    print(config.describe())
    print(resources.describe(RESOURCE_PLAN))
    print(f"  driving frames : {len(files)} from {run_dir}\n")

    levels = PressureLevels(config.PRESSURE_LEVELS)

    fields0, meta0 = load_state(files[0])
    grid = build_grid(fields0, levels)
    print(f"  grid           : {grid}")

    u0, v0, th0 = hrrr_to_model_state(fields0, levels)

    # Boundary frames from the remaining files, hourly.
    times, states = [], []
    for i, f in enumerate(files[:args.hours + 1]):
        fl, _ = load_state(f)
        u, v, th = hrrr_to_model_state(fl, levels)
        times.append(i * 3600.0)
        states.append(state_to_boundary(u, v, th))
    driver = BoundaryDriver(times, states)
    print(f"  boundaries     : {driver}")

    stoch = None
    if args.stochastic:
        stoch = StochasticPerturbation(grid, amplitude=0.3, tau=6 * 3600,
                                       length_scale=300e3, seed=args.seed)
        print(f"  stochastic     : {stoch}")

    model = Primitive3D(grid, levels, stochastic=stoch)
    model.u, model.v, model.theta = u0, v0, th0

    relax = Relaxation3D(grid, width=args.relax_width)
    print(f"  relaxation     : width {args.relax_width}, "
          f"interior {relax.interior_fraction:.0%}")
    print(f"  timestep       : {model.max_dt():.1f} s\n")

    snaps = run_forecast(model, driver, relax, args.hours * 3600,
                         output_every=args.output_every * 3600)

    out = Path(args.out or (run_dir / "forecast.npz"))
    np.savez_compressed(
        out,
        times_s=np.array([s[0] for s in snaps]),
        u=np.stack([s[1] for s in snaps]),
        v=np.stack([s[2] for s in snaps]),
        theta=np.stack([s[3] for s in snaps]),
        levels_hPa=np.array(config.PRESSURE_LEVELS),
        lat=meta0.get("lat"), lon=meta0.get("lon"),
    )
    print(f"\nWrote {len(snaps)} snapshots -> {out}")
    print("Next: verify against observations "
          "(src/verification) and archive the matches.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

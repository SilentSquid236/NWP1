"""
End-to-end forecast driver: initialise from HRRR, integrate our own physics,
verify against observations, archive the result.

    python src/forecast.py --start 2026-08-01T00 --hours 24

WHICH CORE THIS RUNS

The SIGMA core, on terrain-following levels. It used to build a `Primitive3D`
on pressure levels, which was the core measured to diverge in two to three
hours from real analyses (P-14 in docs/PROBLEMS.md) -- so the only core a real
forecast could reach was the broken one, while the one that reaches 12/12 h
could only be run on idealised states. `src/dynamics/interpolate.py` closes
that gap: HRRR arrives isobaric, the model works in sigma, and the conversion
now happens on the way in.

Three things happen to the analysis before the first step, and the order was
measured rather than assumed (see docs/RESEARCH_LOG.md, 2026-09-02):

  1. interpolate on to sigma levels over the real terrain
  2. FILTER -- remove variance at wavelengths the grid cannot carry
  3. REBALANCE -- filtering u, v and theta separately reintroduces divergence

Measured on the idealised equivalent: no filter 1/12 h, filter only 11/12 h,
filter then rebalance 12/12 h.

WHAT THE FORECAST IS BUILT FROM (since 2026-09-22)

  initial conditions   OBSERVATIONS valid at the cycle time (src/analysis/,
                       written by src/ingest_obs.py as obs_analysis_f00.npz).
                       Where soundings are missing, the upper air is this
                       model's previous forecast. No HRRR.
  boundary conditions  the initial analysis, HELD FIXED for the whole run.
                       Nothing observed after the cycle time may enter, so
                       there is nothing else to relax toward. A single
                       driving frame is exactly that: BoundaryDriver returns
                       it at every time.
  verification truth   observations only, after the forecast window closes.

  --source hrrr        the old path (live_hrrr_f*.npz from ingest_hrrr.py,
                       hourly HRRR analyses at the edges) is still readable,
                       as a labelled baseline and nothing more.

The forecast in between is ours: our equations, our numerics, our errors.
"""

import argparse
import faulthandler
import signal
import sys
import time
from datetime import datetime, timedelta
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

import resources
RESOURCE_PLAN = resources.apply()

# See src/verify.py for why: `kill -USR1 <pid>` dumps a traceback without
# stopping the run, and needs nothing installed.
faulthandler.enable()
if hasattr(signal, "SIGUSR1"):
    faulthandler.register(signal.SIGUSR1)

import numpy as np

import config
sys.path.insert(0, str(Path(__file__).resolve().parent / "dynamics"))
sys.path.insert(0, str(Path(__file__).resolve().parent / "verification"))

from grid import CGrid
from vertical import PressureLevels, theta_from_T, T_from_theta
from sigma import SigmaLevels
from primitive_sigma import PrimitiveSigma
from interpolate import pressure_to_sigma, surface_pressure_from_heights
from initialization import filter_initial_state
from boundaries import DaviesRelaxation, BoundaryDriver
from subgrid import StochasticPerturbation, balance_initial_state


# ---------------------------------------------------------------------------
# HRRR state -> model state
# ---------------------------------------------------------------------------

FEATURE_KEYS = ("features", "hrrr_features")


def load_state(path):
    """
    Load one ingested field set (.npz): an observation analysis written by
    src/ingest_obs.py (key "features") or an HRRR frame (key "hrrr_features").
    """
    z = np.load(path, allow_pickle=False)
    key = next((k for k in FEATURE_KEYS if k in z.files), None)
    if key is None:
        raise KeyError(f"{path} has none of {FEATURE_KEYS}; it has {z.files}")
    return z[key], {k: z[k] for k in z.files if k not in FEATURE_KEYS}


def driving_frames(run_dir):
    """
    The frames for a run, observation analysis first.

    An observation run has ONE frame, obs_analysis_f00.npz, and so frozen
    boundaries; an HRRR run has live_hrrr_f00..fNN. A directory holding both
    is refused -- which source drove a forecast must never be a guess.
    """
    run_dir = Path(run_dir)
    key = lambda q: int(q.stem.split("_f")[-1])
    obs = sorted(run_dir.glob("obs_analysis_f*.npz"), key=key)
    hrrr = sorted(run_dir.glob("live_hrrr_f*.npz"), key=key)
    if obs and hrrr:
        raise SystemExit(f"{run_dir} holds both observation and HRRR frames; "
                         f"refusing to guess which drives the run")
    return (obs, "observations") if obs else (hrrr, "hrrr")


def hrrr_channels(fields, channels=None):
    """Pull the named channels out of an ingested [C, L, Y, X] array."""
    channels = list(channels or config.CHANNELS)
    idx = {c: i for i, c in enumerate(channels)}
    for need in ("TMP", "UGRD", "VGRD", "HGT"):
        if need not in idx:
            raise KeyError(
                f"channel {need} missing from {channels}. HGT is required: "
                f"the sigma conversion locates the surface by finding the "
                f"pressure at which the analysis height equals the terrain.")
    return (fields[idx["TMP"]].astype(float),
            fields[idx["UGRD"]].astype(float),
            fields[idx["VGRD"]].astype(float),
            fields[idx["HGT"]].astype(float))


def hrrr_to_sigma_state(fields, lev, terrain, p_surface=None, channels=None):
    """
    Convert an ingested HRRR field set to a sigma-coordinate model state.

    Returns (pi, u, v, theta). Temperature becomes potential temperature on
    the way in: theta is the model's variable and is conserved by dry
    adiabatic motion, so interpolating it does not invent heating the way
    interpolating T through a deep layer does.
    """
    T, u, v, z = hrrr_channels(fields, channels)
    p_pa = np.asarray(config.PRESSURE_LEVELS, dtype=float) * 100.0
    return pressure_to_sigma(u, v, T, z, p_pa, terrain, lev,
                             p_surface=p_surface)


def load_terrain(run_dir, shape):
    """
    Terrain and, if present, surface pressure for the run.

    Refuses to substitute flat ground. A forecast over a flat Northeast is not
    a degraded forecast, it is a different experiment, and silently running it
    is how a result gets misread later.
    """
    path = Path(run_dir) / "terrain.npz"
    if not path.exists():
        raise SystemExit(
            f"No terrain.npz in {run_dir}.\n"
            f"The sigma core is a terrain-following model; running it over "
            f"flat ground in this domain would not be a meaningful forecast.\n"
            f"Re-run: python src/ingest_hrrr.py --start ... "
            f"(it fetches terrain once per run).")
    z = np.load(path, allow_pickle=False)
    terrain = z["terrain"].astype(float)
    p_sfc = z["p_surface"].astype(float) if "p_surface" in z.files else None
    if terrain.shape != shape:
        raise SystemExit(
            f"terrain.npz is {terrain.shape} but the fields are {shape}. "
            f"They were ingested with different --stride settings.")
    return terrain, p_sfc


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


def state_to_boundary(u, v, theta, pi=None):
    """Package a model state as boundary-driver input."""
    out = {"u": u, "v": v, "theta": theta}
    if pi is not None:
        out["pi"] = pi
    return out


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
        # Surface pressure is PROGNOSTIC in the sigma core, so it has to be
        # relaxed at the edges too. Leaving it free while relaxing the wind
        # drives the boundary column toward a mass field the incoming flow
        # does not support.
        if "pi" in ext:
            model.pi += self.alpha[0] * (ext["pi"] - model.pi)

    @property
    def interior_fraction(self):
        return self.inner.interior_fraction


U_CEILING = 150.0     # m/s -- the range limit on observed wind; beyond it
                      # the state is not weather and the run is over (P-52)


def run_forecast(model, driver, relax, duration, dt=None, output_every=None,
                 progress=True, deadline_s=None, info=None):
    """
    Integrate with boundary relaxation, collecting output states.

    Returns a list of (valid_seconds, u, v, theta, pi) snapshots.

    Two ways to stop early, both recorded in `info["stopped"]`:

      deadline  wall clock since the call exceeded `deadline_s`. The hours
                already reached are kept, so a cut-short run can still be
                archived and verified (the 1.5 h budget, prompt 98).
      diverged  checked at the progress cadence, not only on the hour: a
                non-finite value or |u| > U_CEILING. P-52 was a run that kept
                refining a meaningless state for a quarter of an hour; with a
                fixed dt it would not slow down, but it would still spend the
                rest of its budget on nothing.
    """
    info = {} if info is None else info
    info["stopped"] = "completed"
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

    # PROGRESS EVERY FEW SECONDS, NOT EVERY FORECAST HOUR.
    #
    # At dt ~ 15 s a forecast hour is 240 steps and several minutes of wall
    # clock. Printing only on the hour means minutes of silence, which is
    # indistinguishable from a hang -- and that is exactly how the first real
    # run was reported. The ETA is what turns "it is stuck" into "it has 40
    # minutes to go".
    t_start = time.time()
    every_n = max(1, n_steps // 200)

    for k in range(n_steps):
        model.step(dt)
        relax.apply(model, driver.at(model.time))

        if k % every_n == 0:
            # Both components: the P-56 runaways are in v first, and a guard
            # on u alone let one reach 91 m/s unreported (P-57).
            umax = float(max(np.nanmax(np.abs(model.u)),
                             np.nanmax(np.abs(model.v))))
            finite = np.isfinite(model.u).all() and np.isfinite(model.v).all()
            if not finite or umax > U_CEILING:
                print(f"\n  DIVERGED at t+{model.time/3600:.2f} h "
                      f"(max|u,v| {umax:.0f} m/s) -- stopping", flush=True)
                info["stopped"] = f"diverged at {model.time/3600:.2f} h"
                break
            if deadline_s is not None and time.time() - t_start > deadline_s:
                print(f"\n  DEADLINE reached at t+{model.time/3600:.2f} h "
                      f"after {(time.time()-t_start)/60:.1f} min -- stopping",
                      flush=True)
                info["stopped"] = f"deadline at {model.time/3600:.2f} h"
                break

        if progress and k and k % every_n == 0:
            el = time.time() - t_start
            rate = (k + 1) / el
            eta = (n_steps - k - 1) / rate
            print(f"    step {k+1}/{n_steps}  "
                  f"t+{model.time/3600:4.2f} h  "
                  f"{rate:.1f} steps/s  "
                  f"elapsed {el/60:.1f} min  ETA {eta/60:.1f} min",
                  end="\r", flush=True)

        if next_i < len(targets) and model.time >= targets[next_i] - 1e-9:
            next_i += 1
            snapshots.append((model.time, model.u.copy(), model.v.copy(),
                              model.theta.copy(), model.pi.copy()))
            if progress:
                print(" " * 96, end="\r")      # clear the progress line
                ps = model.surface_pressure
                print(f"  +{model.time/3600:5.1f} h  "
                      f"max|u| {np.abs(model.u).max():6.1f} m/s  "
                      f"max|v| {np.abs(model.v).max():6.1f}  "
                      f"theta {model.theta.min():.1f}-{model.theta.max():.1f} K  "
                      f"p_s {ps.min()/100:.0f}-{ps.max()/100:.0f} hPa  "
                      f"max|sigma_dot| {np.abs(model.sigma_dot()).max():.2e}")
            if not (np.isfinite(model.u).all() and np.isfinite(model.v).all()):
                print("  FORECAST DIVERGED -- stopping")
                info["stopped"] = f"diverged at {model.time/3600:.2f} h"
                break

    info["wall_s"] = time.time() - t_start
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
    p.add_argument("--relax-alpha", type=float, default=1.0,
                   help="Relaxation weight per step at the outer edge "
                        "(the cosine ramp scales from it); 0 < alpha <= 1")
    p.add_argument("--stochastic", action="store_true",
                   help="Enable SPPT-style tendency perturbations")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--no-balance", action="store_true",
                   help="Skip initial divergence removal. The forecast will "
                        "almost certainly blow up; useful only for showing "
                        "why the balancing step exists.")
    p.add_argument("--out", default=None, help="Where to write forecast .npz")
    p.add_argument("--deadline-min", type=float, default=None,
                   help="Stop integrating after this many minutes of wall "
                        "clock and write the hours reached (cycle budget)")
    args = p.parse_args()

    run_dir = Path(args.run_dir)
    files, source = driving_frames(run_dir)
    if not files:
        print(f"No obs_analysis_f*.npz or live_hrrr_f*.npz in {run_dir}. "
              f"Run src/ingest_obs.py (or ingest_hrrr.py) first.")
        return 1

    print("NWP forecast")
    print(config.describe())
    print(resources.describe(RESOURCE_PLAN))
    print(f"  driving frames : {len(files)} from {run_dir} ({source})")
    if len(files) == 1:
        print("  boundaries     : held at the initial state for the whole run "
              "(nothing observed later may enter)")
    print()

    levels = PressureLevels(config.PRESSURE_LEVELS)
    # Number of SIGMA levels, deliberately the same count as the analysis has
    # pressure levels -- not because they must match, but because a different
    # count would silently change the vertical resolution of every result
    # measured so far.
    lev = SigmaLevels(config.N_LEVELS)

    fields0, meta0 = load_state(files[0])
    grid = build_grid(fields0, levels)
    print(f"  grid           : {grid}")
    print(f"  vertical       : {lev}")

    terrain, p_sfc = load_terrain(run_dir, fields0.shape[-2:])
    print(f"  terrain        : {terrain.min():.0f}-{terrain.max():.0f} m"
          + ("" if p_sfc is None else "   (surface pressure from HRRR)"))

    pi0, u0, v0, th0 = hrrr_to_sigma_state(fields0, lev, terrain,
                                           p_surface=p_sfc)

    # Boundary frames. Each is put through the SAME conversion as the initial
    # state -- if the edges were prepared differently from the interior, the
    # relaxation would drive one toward the other every step.
    times, states = [], []
    for i, f in enumerate(files[:args.hours + 1]):
        fl, _ = load_state(f)
        pi_b, u, v, th = hrrr_to_sigma_state(fl, lev, terrain,
                                             p_surface=p_sfc)
        if not args.no_balance:
            u, v, th = filter_initial_state(u, v, th, grid)
            u, v, _ = balance_initial_state(u, v, grid, verbose=False)
        times.append(i * 3600.0)
        states.append(state_to_boundary(u, v, th, pi_b))
    driver = BoundaryDriver(times, states)
    print(f"  boundaries     : {driver}")

    stoch = None
    if args.stochastic:
        stoch = StochasticPerturbation(grid, amplitude=0.3, tau=6 * 3600,
                                       length_scale=300e3, seed=args.seed)
        print(f"  stochastic     : {stoch}")

    model = PrimitiveSigma(grid, lev, terrain=terrain, stochastic=stoch)

    # PREPARE THE INITIAL STATE. Order measured, not assumed.
    #
    #   filter    -- white grid-scale variance is amplified by advection
    #                faster than hyperdiffusion removes it. Measured
    #                threshold: 0.30 m/s survives 12 h, 0.60 m/s does not.
    #   rebalance -- filtering u, v and theta separately puts divergence back
    #                into a balanced state.
    #
    # Measured on the idealised equivalent: none 1/12 h, filter only 11/12 h,
    # filter then rebalance 12/12 h. Note that the filter-only case has HIGHER
    # divergence than the unfiltered one and still survives ten hours longer:
    # wavenumber content is the controlling variable, not divergence.
    if not args.no_balance:
        u0, v0, th0 = filter_initial_state(u0, v0, th0, grid)
        u0, v0, binfo = balance_initial_state(u0, v0, grid)
        if binfo.get("omega_after_Pa_s", 0.0) > 5.0:
            print("  WARNING: initial divergence is still large; expect a "
                  "noisy first hour.")

    # Assign COPIES. The relaxation updates the model state in place, and the
    # same arrays are still referenced by the boundary frames; sharing them
    # would let the first relaxation step quietly rewrite the driving data.
    model.pi = pi0.copy()
    model.u, model.v, model.theta = u0.copy(), v0.copy(), th0.copy()

    ps = model.surface_pressure
    print(f"  initial state  : max|u| {np.abs(model.u).max():.1f} m/s, "
          f"p_s {ps.min()/100:.0f}-{ps.max()/100:.0f} hPa, "
          f"max|sigma_dot| {np.abs(model.sigma_dot()).max():.2e} 1/s")

    if not 0.0 < args.relax_alpha <= 1.0:
        raise SystemExit(f"--relax-alpha {args.relax_alpha} is outside (0, 1]")
    relax = Relaxation3D(grid, width=args.relax_width,
                         alpha_max=args.relax_alpha)
    print(f"  relaxation     : width {args.relax_width}, "
          f"alpha {args.relax_alpha:g} at the edge, "
          f"interior {relax.interior_fraction:.0%}")
    print(f"  timestep       : {model.max_dt():.1f} s "
          f"(external wave ~290 m/s sets this)\n")

    info = {}
    snaps = run_forecast(model, driver, relax, args.hours * 3600,
                         output_every=args.output_every * 3600,
                         deadline_s=(None if args.deadline_min is None
                                     else 60.0 * args.deadline_min),
                         info=info)
    if not snaps:
        print(f"\nNo forecast hour completed ({info.get('stopped')}); "
              f"nothing written.")
        return 1

    out = Path(args.out or (run_dir / "forecast.npz"))
    np.savez_compressed(
        out,
        times_s=np.array([s[0] for s in snaps]),
        u=np.stack([s[1] for s in snaps]),
        v=np.stack([s[2] for s in snaps]),
        theta=np.stack([s[3] for s in snaps]),
        pi=np.stack([s[4] for s in snaps]),
        sigma=lev.sigma,
        p_top=lev.p_top,
        terrain=terrain,
        lat=meta0.get("lat"), lon=meta0.get("lon"),
        run_time=meta0.get("run_time", np.array("")),
        source=np.array(source),
        hours_requested=args.hours,
        stopped=np.array(info.get("stopped", "")),
        wall_s=info.get("wall_s", np.nan),
    )
    print(f"\nWrote {len(snaps)} snapshots -> {out}  ({info.get('stopped')}, "
          f"{info.get('wall_s', float('nan'))/60:.1f} min)")
    print("Next: verify once the forecast window has closed "
          "(src/verify_pending.py).")
    stopped = info.get("stopped", "")
    # Output is written in every case; the code says how the run ended, so the
    # cycle script's log shows a cut-short or diverged run as such.
    return 0 if stopped == "completed" else (2 if stopped.startswith("deadline") else 3)


if __name__ == "__main__":
    raise SystemExit(main())

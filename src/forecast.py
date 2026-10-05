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
from backend import xp_of, to_numpy
from interpolate import pressure_to_sigma, surface_pressure_from_heights
from initialization import filter_initial_state
from boundaries import DaviesRelaxation, BoundaryDriver
from subgrid import StochasticPerturbation, balance_initial_state
from convection import dry_convective_adjustment, unstable_fraction


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


FRAME_MAX_SWEEPS = 1000   # a one-off per frame, so run to convergence


def stabilise_frame(theta, u, v, pi, lev, columns=None):
    """
    Remove static instability from a boundary frame, once (P-63).

    Uses the same mass-weighted adjustment the model applies after every
    step, but run to convergence. Returns (theta, u, v, info); a frame that is
    already stable comes back unchanged. Raises if it does not converge,
    because an edge that stays unstable would bring P-63 straight back.

    columns (CAM stage S5d): a 2-D boolean mask of the columns that matter
    -- the relaxation zone. The adjustment is column by column, so only those
    columns are adjusted; the others are returned as they were. At 3 km a
    whole-frame adjustment took 16 s per hourly frame.
    """
    if columns is not None:
        cols = np.flatnonzero(np.asarray(columns).ravel())
        nz = theta.shape[0]
        sub = lambda a: np.ascontiguousarray(a.reshape(a.shape[0], -1)[:, cols])[:, None, :]
        th_s, u_s, v_s = sub(theta), sub(u), sub(v)
        pi_s = np.ascontiguousarray(np.asarray(pi).ravel()[cols])[None, :]
        th_s, u_s, v_s, info = stabilise_frame(th_s, u_s, v_s, pi_s, lev)
        out = []
        for full, part in ((theta, th_s), (u, u_s), (v, v_s)):
            a = np.array(full, copy=True)
            a.reshape(nz, -1)[:, cols] = part[:, 0, :]
            out.append(a)
        return out[0], out[1], out[2], info
    info = {"unstable_before": unstable_fraction(theta), "sweeps": 0,
            "unstable_after": 0.0}
    if info["unstable_before"] == 0:
        return theta, u, v, info
    theta, u, v, info = dry_convective_adjustment(
        theta, u, v, pi, lev, max_sweeps=FRAME_MAX_SWEEPS)
    if info["unstable_after"] > 0:
        raise RuntimeError(f"boundary frame still {info['unstable_after']:.2e} "
                           f"unstable after {info['sweeps']} sweeps")
    return theta, u, v, info

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
    Davies (1976) relaxation applied to a 3D state.

    The 2D weight field is reused at every level: the lateral boundary is a
    vertical wall, so the taper depends only on horizontal distance from the
    edge.
    """

    def __init__(self, grid, width=10, alpha_max=1.0, profile="cosine"):
        self.inner = DaviesRelaxation(grid, width, alpha_max, profile)
        self.alpha = self.inner.alpha[None, :, :]
        self.alpha2d = self.alpha[0]          # one array, so torch caches it once
        self.width = width

    def apply(self, model, ext):
        # Backend-neutral: for torch the weights and the driving state are
        # converted (and cached) so a NumPy array never meets a tensor.
        xp = xp_of(model.u)
        a = xp.asarray(self.alpha)
        if "u" in ext:
            model.u += a * (xp.asarray(ext["u"]) - model.u)
        if "v" in ext:
            model.v += a * (xp.asarray(ext["v"]) - model.v)
        if "theta" in ext:
            model.theta += a * (xp.asarray(ext["theta"]) - model.theta)
        # Surface pressure is PROGNOSTIC in the sigma core, so it has to be
        # relaxed at the edges too. Leaving it free while relaxing the wind
        # drives the boundary column toward a mass field the incoming flow
        # does not support.
        if "pi" in ext:
            model.pi += xp.asarray(self.alpha2d) * (xp.asarray(ext["pi"]) - model.pi)

    @property
    def interior_fraction(self):
        return self.inner.interior_fraction


Z0_LAND_DEFAULT = 1.0     # m, land roughness since 2026-10-03 (test AM)
Z0_SEA_DEFAULT = 0.0002   # m, water roughness since 2026-10-03 (test AM)

U_CEILING = 150.0     # m/s -- the range limit on observed wind; beyond it
                      # the state is not weather and the run is over (P-52)


def _land_mask(terrain):
    """Land = raw ETOPO terrain above 0 m; falls back to the run terrain."""
    # LAND FROM THE UNSMOOTHED TERRAIN (P-65). The run's terrain.npz is
    # slope-limited, and the smoothing spreads land heights out over the
    # sea: 15.5 % of it is exactly 0 m against 31.9 % of the raw ETOPO grid
    # (ocean clipped to 0). "terrain > 0" on it called half the ocean land.
    raw = Path(config.DATA_ROOT) / "static" / f"terrain_etopo_{terrain.shape[0]}x{terrain.shape[1]}.npz"
    if raw.exists():
        return np.load(raw)["terrain"] > 0.0, raw.name
    land_src = "the smoothed run terrain (no raw ETOPO file; coastal sea counted as land)"
    print(f"  WARNING: {raw} not found; land mask from {land_src}")
    return terrain > 0.0, land_src


def run_forecast(model, driver, relax, duration, dt=None, output_every=None,
                 progress=True, deadline_s=None, info=None, snapshot_dtype=None):
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
    interval = output_every or duration

    # OUTPUT TIMES ON THE HOUR (P-64). When the output interval divides the
    # run, take a whole number of (slightly shorter) steps per interval, so
    # every snapshot falls exactly on its target. Before this the 6 h
    # snapshot sat 7.8 s past 6 h (one step's overshoot), the next cycle's
    # ingest asks for 6 h to within 3.6 s, and so no cycle ever used the
    # previous forecast as its first guess.
    n_out = duration / interval
    if round(n_out) >= 1 and abs(n_out - round(n_out)) < 1e-9:
        per_out = int(np.ceil(interval / dt - 1e-9))
        n_steps = per_out * int(round(n_out))
    else:
        n_steps = int(np.ceil(duration / dt))
    dt = duration / n_steps

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
    xp = xp_of(model.u)
    torch_run = xp.name == "torch"

    for k in range(n_steps):
        model.step(dt)
        if hasattr(model, "relax_with_driver"):
            model.relax_with_driver(relax, driver, model.time)
        elif hasattr(model, "relax_with"):
            model.relax_with(relax, driver.at(model.time))
        else:
            relax.apply(model, driver.at(model.time))

        if k % every_n == 0:
            # Both components: the P-56 runaways are in v first, and a guard
            # on u alone let one reach 91 m/s unreported (P-57).
            if torch_run:
                finite = bool(xp.isfinite(model.u).all()) and bool(xp.isfinite(model.v).all())
                umax = float(max(xp.abs(model.u).max(), xp.abs(model.v).max())) \
                    if finite else float("inf")
            else:
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

        # The LAST step always writes the final target. model.time is a sum of
        # n_steps floats and can end a few ns short of `duration`, so a 1e-9 s
        # tolerance never fired there: every run lost its last snapshot (24 h
        # runs ended at 23.75 h; found 2026-09-26 in the test X verification).
        #
        # Half a step of tolerance: with the step count above, a target is a
        # whole number of steps away and model.time reaches it to round-off
        # (P-64). The snapshot is then stamped with the target itself.
        if next_i < len(targets) and (model.time >= targets[next_i] - 0.5 * dt
                                      or k == n_steps - 1):
            target = targets[next_i]
            next_i += 1
            snap = [to_numpy(a) if torch_run else a.copy()
                    for a in (model.u, model.v, model.theta, model.pi)]
            if snapshot_dtype is not None:
                # S5e: kept at the output precision (float32 for the 3 km
                # run), so the 24 snapshots take half the memory. The file
                # is the same: _write_forecast casts to this dtype anyway.
                snap = [np.asarray(a, dtype=snapshot_dtype) for a in snap]
            stamp = target if abs(model.time - target) < 1e-6 * dt else model.time
            snapshots.append((stamp, *snap))
            if getattr(model, "land_surface", None) is not None:
                info.setdefault("tg", []).append(np.array(model.land_surface.Tg, copy=True))
            if progress:
                print(" " * 96, end="\r")      # clear the progress line
                su, sv, sth, spi = snap
                ps = model.lev.p_top + spi
                print(f"  +{model.time/3600:5.1f} h  "
                      f"max|u| {np.abs(su).max():6.1f} m/s  "
                      f"max|v| {np.abs(sv).max():6.1f}  "
                      f"theta {sth.min():.1f}-{sth.max():.1f} K  "
                      f"p_s {ps.min()/100:.0f}-{ps.max()/100:.0f} hPa  "
                      f"max|sigma_dot| {float(xp.abs(model.sigma_dot()).max()):.2e}")
            if not (np.isfinite(snap[0]).all() and np.isfinite(snap[1]).all()):
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
    # DEFAULTS CHANGED 2026-09-26 (P-56 tests P and W, P-60 tests V and W):
    # width 15, alpha 0.1 and divergence damping C = 0.0064. Both real cases
    # (06Z 2026-09-23 calm, 12Z 2026-09-25 jet) then run 24 h clean. The old
    # behaviour is --relax-width 10 --relax-alpha 1 --div-damp 0.
    p.add_argument("--relax-width", type=int, default=15)
    p.add_argument("--backend", choices=("numpy", "torch"), default="numpy",
                   help="array backend for the core: numpy (one core) or "
                        "torch (multi-threaded CPU, same float64 physics)")
    p.add_argument("--threads", type=int, default=8,
                   help="torch threads (the Xeon's measured sweet spot is ~8)")
    p.add_argument("--relax-alpha", type=float, default=0.1,
                   help="Relaxation weight per step at the outer edge "
                        "(the cosine ramp scales from it); 0 < alpha <= 1")
    p.add_argument("--dt-factor", type=float, default=1.0,
                   help="Multiply the CFL timestep (P-60 test T; must be in (0, 1])")
    p.add_argument("--hyper-factor", type=float, default=1.0,
                   help="Multiply the recommended hyperdiffusion coefficient "
                        "(P-60 test T)")
    p.add_argument("--div-damp", type=float, default=0.0064,
                   help="Divergence damping, as C in nu = C dx dy / dt "
                        "(Skamarock and Klemp 1992); 0 is off (P-60 tests V, W)")
    p.add_argument("--ri-crit", type=float, default=None,
                   help="Richardson number below which vertical mixing acts "
                        "(default: turbulence.RI_CRIT, 0.25; P-60 test S)")
    p.add_argument("--no-mixing", action="store_true",
                   help="Turn off turbulent vertical mixing (diagnosis only)")
    p.add_argument("--core", choices=("hydrostatic", "nh"), default="hydrostatic",
                   help="Dynamical core: the hydrostatic sigma model (default) or the "
                        "non-hydrostatic mass-coordinate core (nh3d.py; CAM stage S2; NumPy only)")
    p.add_argument("--nh-backend", choices=("numpy", "c"), default="numpy",
                   help="With --core nh: NumPy (the reference) or compiled C kernels with "
                        "OpenMP (nh3d_kernels.c, built with gcc at first use; CAM stage S4). "
                        "--threads sets the OpenMP thread count")
    p.add_argument("--nh-physics", choices=("stage", "step"), default="stage",
                   help="With --core nh: evaluate the physics in every RK stage (default) or "
                        "once per step, held over the stages (CAM stage S5)")
    p.add_argument("--boundary-forecast", default=None,
                   help="Nest: lateral boundaries from this coarser forecast.npz of the "
                        "same domain and cycle, interpolated to this grid every hour "
                        "(nest.py; CAM stage S5). Default: the analysis held fixed")
    p.add_argument("--advection", choices=("centred2", "upwind3"), default="upwind3",
                   help="Horizontal advection: third-order upwind-biased (default since "
                        "2026-10-03; Wicker and Skamarock 2002; P-67 tests AL, AN) or "
                        "second-order centred (the old default)")
    p.add_argument("--land-surface", dest="land_surface", action="store_true", default=None,
                   help="Force-restore ground temperature with a surface energy "
                        "budget (land_surface.py; P-59, tests AK, AM, AN). ON by default "
                        "since 2026-10-03 for the hydrostatic core; off with "
                        "--surface-heating or --core nh")
    p.add_argument("--no-land-surface", dest="land_surface", action="store_false",
                   help="Turn the land surface off (the pre-2026-10-03 configuration)")
    p.add_argument("--z0", type=float, default=None,
                   help="One roughness length everywhere, m (test AJ). Without it the "
                        "default since 2026-10-03 is a land/sea map: --z0-land 1.0, "
                        "--z0-sea 0.0002 (test AM). --z0 0.1 is the old default")
    p.add_argument("--z0-land", type=float, default=None,
                   help="Roughness over land, m (default 1.0 unless --z0 is given; "
                        "land from the raw ETOPO mask)")
    p.add_argument("--z0-sea", type=float, default=None,
                   help="Roughness over water, m (default 0.0002 unless --z0 is given)")
    p.add_argument("--no-drag", action="store_true",
                   help="Turn off surface drag (diagnosis only; test AJ)")
    p.add_argument("--sponge-levels", type=int, default=None,
                   help="Levels below the lid in the wind sponge (P-60 test R); "
                        "default 5 at 20 levels (the measured choice), scaled with "
                        "--levels to keep the same depth in sigma")
    p.add_argument("--levels", type=int, default=None,
                   help="Number of sigma levels (default: the analysis's pressure-level "
                        "count, 20). CAM stage S5 uses 40")
    p.add_argument("--stochastic", action="store_true",
                   help="Enable SPPT-style tendency perturbations (Buizza et al. 1999)")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument("--no-balance", action="store_true",
                   help="Skip initial divergence removal. The forecast will "
                        "almost certainly blow up; useful only for showing "
                        "why the balancing step exists.")
    p.add_argument("--conv-scheme", choices=("sweep", "pav"), default="pav",
                   help="Convective adjustment: pav (pool-adjacent-violators, "
                        "exact in one pass; the default since test AC, P-63) "
                        "or sweep (segment mixing, cap 20; the old model)")
    p.add_argument("--surface-heating", action="store_true",
                   help="Prescribed diurnal surface heat flux from solar "
                        "elevation over land (P-59; src/dynamics/diurnal.py)")
    p.add_argument("--sh-fraction", type=float, default=0.2,
                   help="Share of clear-sky sunshine into sensible heat")
    p.add_argument("--sh-tau", type=float, default=0.75,
                   help="Clear-sky transmission")
    p.add_argument("--sh-night", type=float, default=-30.0,
                   help="Steady ground cooling flux, W/m2 (negative)")
    p.add_argument("--persistence", action="store_true",
                   help="Write the prepared initial state at every output "
                        "time instead of integrating: the do-nothing "
                        "reference forecast for verification (P-59)")
    p.add_argument("--no-conv-momentum", action="store_true",
                   help="The convective adjustment mixes theta only, not u and v "
                        "(diagnosis, P-67)")
    p.add_argument("--raw-boundaries", action="store_true",
                   help="Do not convectively adjust the boundary frames "
                        "(the behaviour before P-63; comparison only)")
    p.add_argument("--out", default=None, help="Where to write forecast .npz")
    p.add_argument("--output-dtype", choices=("auto", "f64", "f32"), default="auto",
                   help="precision of the written fields (default auto: float32 for grids "
                        f"above {LARGE_GRID} columns, i.e. the 3 km run; float64 otherwise)")
    p.add_argument("--output-level", type=int, default=None, choices=range(0, 10),
                   metavar="0-9", help="zlib level for the output (0 = none); default: "
                   "np.savez_compressed for float64, level 1 for float32")
    p.add_argument("--deadline-min", type=float, default=None,
                   help="Stop integrating after this many minutes of wall "
                        "clock and write the hours reached (cycle budget)")
    args = p.parse_args()

    # DEFAULTS CHANGED 2026-10-03 (user decision after the holdout, test AN):
    # upwind3 advection, the force-restore land surface, and a land/sea
    # roughness map. The old configuration is
    #   --advection centred2 --no-land-surface --z0 0.1
    if args.z0 is None:
        if args.z0_land is None and args.z0_sea is None:
            args.z0_land, args.z0_sea = Z0_LAND_DEFAULT, Z0_SEA_DEFAULT
        args.z0 = 0.1
    if args.land_surface is None:
        args.land_surface = not args.surface_heating and args.core == "hydrostatic"
        if args.core != "hydrostatic" and not args.surface_heating:
            print("  NOTE: --core nh has no land surface yet; running without it")

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
    if len(files) == 1 and not args.boundary_forecast:
        print("  boundaries     : held at the initial state for the whole run "
              "(nothing observed later may enter)")
    print()

    levels = PressureLevels(config.PRESSURE_LEVELS)
    # Number of SIGMA levels, deliberately the same count as the analysis has
    # pressure levels -- not because they must match, but because a different
    # count would silently change the vertical resolution of every result
    # measured so far.
    nlev = config.N_LEVELS if args.levels is None else int(args.levels)
    if nlev < 4:
        raise SystemExit(f"--levels {nlev}: need at least 4")
    lev = SigmaLevels(nlev)
    if args.sponge_levels is None:
        # Same depth in sigma as the measured 5 of 20: the top quarter of the
        # half levels, which with the 1.4 stretch is sigma < 0.14 either way.
        args.sponge_levels = max(1, int(round(5 * nlev / 20)))

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
    n_unstable_frames, frame_sweeps, frame_frac = 0, 0, 0.0
    for i, f in enumerate(files[:args.hours + 1]):
        fl, _ = load_state(f)
        pi_b, u, v, th = hrrr_to_sigma_state(fl, lev, terrain,
                                             p_surface=p_sfc)
        if not args.no_balance:
            u, v, th = filter_initial_state(u, v, th, grid)
            u, v, _ = balance_initial_state(u, v, grid, verbose=False)
        # ADJUST THE DRIVING STATE ONCE (P-63). The model removes static
        # instability after every step, but the relaxation then pulls the
        # edge columns back toward the frame. A frame with an unstable layer
        # (a daytime superadiabatic surface layer, test Z2: all of it within
        # 14 cells of the edge, the lowest two interfaces) re-creates the
        # instability every step, and the adjustment ran to its sweep cap on
        # every call. The frame gets the same adjustment the model applies,
        # run to convergence, so the edges drive toward a state the model
        # itself can hold.
        if not args.raw_boundaries:
            th, u, v, cinfo = stabilise_frame(th, u, v, pi_b, lev)
            if cinfo["unstable_before"] > 0:
                n_unstable_frames += 1
                frame_sweeps = max(frame_sweeps, cinfo["sweeps"])
                frame_frac = max(frame_frac, cinfo["unstable_before"])
        times.append(i * 3600.0)
        states.append(state_to_boundary(u, v, th, pi_b))
    if args.boundary_forecast:
        # NESTED RUN (CAM stage S5): the edges follow a coarser forecast of the
        # same cycle, hour by hour, instead of the analysis held fixed. Hour 0
        # is this run's own initial state; later hours come from nest.py.
        import nest
        nf = nest.frames_from_forecast(args.boundary_forecast, terrain, lev.sigma, args.hours)
        times, states = times[:1], states[:1]
        # Only the relaxation zone (plus one cell) ever reads a frame.
        zone = Relaxation3D(grid, width=args.relax_width, alpha_max=args.relax_alpha).alpha2d > 0
        zone_d = zone.copy()
        zone_d[1:] |= zone[:-1]; zone_d[:-1] |= zone[1:]; zone_d[:, 1:] |= zone[:, :-1]; zone_d[:, :-1] |= zone[:, 1:]
        for t, u, v, th, pi_b in nf:
            if t <= 0.0:
                continue
            if not args.raw_boundaries:
                th, u, v, cinfo = stabilise_frame(th, u, v, pi_b, lev, columns=zone_d)
                if cinfo["unstable_before"] > 0:
                    n_unstable_frames += 1
                    frame_sweeps = max(frame_sweeps, cinfo["sweeps"])
                    frame_frac = max(frame_frac, cinfo["unstable_before"])
            times.append(t)
            states.append(state_to_boundary(u, v, th, pi_b))
        print(f"  nested in      : {args.boundary_forecast} ({len(nf)} hourly frames)")
        if times[-1] < args.hours * 3600.0 - 1e-6:
            print(f"  NOTE: the coarse forecast ends at {times[-1] / 3600:.2f} h; the edges are "
                  f"held at its last frame after that")
    driver = BoundaryDriver(times, states)
    print(f"  boundaries     : {driver}")
    if args.raw_boundaries:
        print("  frame adjust   : off (--raw-boundaries)")
    else:
        print(f"  frame adjust   : {n_unstable_frames} of {len(states)} frames "
              f"unstable (max {frame_frac:.2e} of interfaces), "
              f"converged in <= {frame_sweeps} sweeps")

    stoch = None
    if args.stochastic:
        stoch = StochasticPerturbation(grid, amplitude=0.3, tau=6 * 3600,
                                       length_scale=300e3, seed=args.seed)
        print(f"  stochastic     : {stoch}")

    if not 0 <= args.sponge_levels < lev.nz:
        raise SystemExit(f"--sponge-levels {args.sponge_levels} is outside 0..{lev.nz - 1}")
    if args.ri_crit is not None and not args.ri_crit > 0:
        raise SystemExit(f"--ri-crit {args.ri_crit} must be positive")
    if not 0.0 < args.dt_factor <= 1.0:
        raise SystemExit(f"--dt-factor {args.dt_factor} is outside (0, 1]")
    if not args.hyper_factor >= 0.0:
        raise SystemExit(f"--hyper-factor {args.hyper_factor} is negative")
    hyper = None
    if args.hyper_factor != 1.0:
        from subgrid import recommended_hyper_coeff
        hyper = args.hyper_factor * recommended_hyper_coeff(grid)
    z0 = args.z0
    z0_note = f"z0 = {args.z0:g} m"
    if args.z0_land is not None or args.z0_sea is not None:
        zl = args.z0 if args.z0_land is None else args.z0_land
        zs = args.z0 if args.z0_sea is None else args.z0_sea
        if not (zl > 0 and zs > 0):
            raise SystemExit("--z0-land and --z0-sea must be positive")
        land_z0, src_z0 = _land_mask(terrain)
        z0 = np.where(land_z0, zl, zs)
        z0_note = f"z0 = {zl:g} m over land, {zs:g} m over water ({src_z0})"
    model = PrimitiveSigma(grid, lev, terrain=terrain, stochastic=stoch,
                           sponge_levels=args.sponge_levels, hyper=hyper,
                           ri_crit=args.ri_crit, mixing=not args.no_mixing,
                           drag=not args.no_drag, z0=z0)
    if args.hyper_factor != 1.0:
        print(f"  hyperdiffusion : x{args.hyper_factor:g} ({model.hyper:.3g})")
    if not args.z0 > 0:
        raise SystemExit(f"--z0 {args.z0} must be positive")
    if args.div_damp < 0:
        raise SystemExit(f"--div-damp {args.div_damp} is negative")
    print(f"  mixing         : "
          + ("off" if args.no_mixing else f"on below Ri {model.ri_crit:g}"))
    print(f"  drag           : "
          + ("off" if args.no_drag else f"on, {z0_note}"))
    print(f"  sponge         : {args.sponge_levels} levels below the lid")
    if args.surface_heating:
        from diurnal import DiurnalHeating
        if not all(k in meta0 for k in ("lat", "lon", "run_time")):
            raise SystemExit("--surface-heating needs lat, lon and run_time in the "
                             f"driving frame; {files[0].name} has {sorted(meta0)}")
        start = np.datetime64(str(np.asarray(meta0["run_time"])), "s").astype(datetime)
        land, land_src = _land_mask(terrain)
        model.surface_heating = DiurnalHeating(
            np.asarray(meta0["lat"]), np.asarray(meta0["lon"]), land, start,
            f_sensible=args.sh_fraction, tau=args.sh_tau, h_night=args.sh_night)
        print(f"  surface heat   : {model.surface_heating}")
    model.conv_scheme = args.conv_scheme
    model.advection = args.advection
    if args.advection != "centred2":
        print(f"  advection      : {args.advection} (horizontal)")
    model.conv_mix_momentum = not args.no_conv_momentum
    if args.no_conv_momentum:
        print("  convection     : momentum NOT mixed (--no-conv-momentum)")
    print(f"  convection     : {args.conv_scheme}"
          + (" (cap 20 sweeps)" if args.conv_scheme == "sweep"
             else " (pool-adjacent-violators, one pass)"))

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
    if args.core == "nh":
        # STAGE S2: the non-hydrostatic core behind the same interface. The
        # hydrostatic model built above supplies the grid, terrain and physics
        # settings; its prepared state is loaded with w = 0 and a hydrostatic phi.
        from nh3d import NHModel
        if args.land_surface or args.surface_heating:
            raise SystemExit("--core nh does not yet support --land-surface or --surface-heating")
        if args.backend == "torch":
            print("  NOTE: --core nh runs on NumPy; --backend torch ignored")
            args.backend = "numpy"
        hydro = model
        model = NHModel(hydro, backend=args.nh_backend, threads=args.threads,
                        physics_every=args.nh_physics)
        model.set_state(u0, v0, th0, pi0)
        if args.nh_backend == "c":
            print(f"  NH kernels     : compiled C (nh3d_kernels.c), {model.core.threads} OpenMP threads")
        if args.nh_physics == "step":
            print("  NH physics     : once per step, held over the RK stages")
        print(f"  core           : non-hydrostatic (nh3d.py, ns {model.core.ns}, dt <= {model.dt_max:g} s; "
              f"edge relaxation rescaled to the hydrostatic step {model.relax_dt_ref:.1f} s)")

    if args.land_surface:
        if args.surface_heating:
            raise SystemExit("--land-surface replaces --surface-heating; use one")
        if not all(k in meta0 for k in ("lat", "lon", "run_time")):
            raise SystemExit("--land-surface needs lat, lon and run_time in the "
                             f"driving frame; {files[0].name} has {sorted(meta0)}")
        from land_surface import ForceRestoreSurface
        from surface import lowest_level_height
        from sigma import P0, KAPPA
        start = np.datetime64(str(np.asarray(meta0["run_time"])), "s").astype(datetime)
        land, land_src = _land_mask(terrain)
        p1 = lev.p_top + float(np.asarray(lev.sigma)[-1]) * pi0
        z1 = lowest_level_height(th0, pi0, lev)
        model.land_surface = ForceRestoreSurface(
            np.asarray(meta0["lat"]), np.asarray(meta0["lon"]), land, start,
            th0[-1] * (p1 / P0) ** KAPPA, z1, z0=z0)
        model.theta_surface = model.land_surface.Tg / ((lev.p_top + pi0) / P0) ** KAPPA
        print(f"  land surface   : {model.land_surface} (land from {land_src}); "
              f"drag stability Louis (1979)")

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
    if args.div_damp > 0:
        # nu from the CFL step of the LOADED initial state. (Before 2026-09-26
        # this ran before the state was assigned, when u = v = 0 and theta was
        # unset, so dt came out 24.5 s instead of 15.8 s on the Q case and nu
        # was 0.64 of the intended C dx dy / dt. Test V's logs print the nu
        # actually used, which is what its record quotes.)
        dt0 = model.max_dt()
        model.div_damp = args.div_damp * grid.dx * grid.dy / dt0
        rate = 4.0 * model.div_damp / grid.dx ** 2
        print(f"  div. damping   : C {args.div_damp:g}, nu {model.div_damp:.3g} m2/s "
              f"(2dx e-folds in {1.0 / rate / 60.0:.1f} min, 10dx in "
              f"{1.0 / (rate * np.sin(np.pi / 10) ** 2) / 60.0:.0f} min)")
    print(f"  timestep       : {model.max_dt():.1f} s "
          f"(external wave ~290 m/s sets this)\n")

    if args.persistence:
        # THE DO-NOTHING REFERENCE (P-59). The prepared initial state --
        # same conversion, filter and balance as the model starts from -- is
        # written at every output time, so verification scores it exactly
        # the way it scores the model. Skill means beating this.
        interval = args.output_every * 3600
        duration = args.hours * 3600
        times = list(np.arange(interval, duration + 1e-9, interval))
        if not times or times[-1] < duration - 1e-9:
            times.append(duration)
        snaps = [(float(t), to_numpy(model.u).copy(), to_numpy(model.v).copy(),
                  to_numpy(model.theta).copy(), to_numpy(model.pi).copy())
                 for t in times]
        info = {"stopped": "persistence", "wall_s": 0.0}
        if model.land_surface is not None:
            # Held like everything else: the ground temperature at the cycle
            # time, so the similarity operator can score persistence too.
            info["tg"] = [np.array(model.land_surface.Tg, copy=True) for _ in times]
        print(f"  persistence    : initial state written at {len(snaps)} times, "
              f"no integration\n")
        return _write_forecast(args, run_dir, snaps, info, lev, terrain, meta0, source)

    if args.backend == "torch":
        n = model.to_backend("torch", threads=args.threads)
        print(f"  backend        : torch, {n} threads (float64)\n")
    else:
        print("  backend        : numpy (one core)\n")

    info = {}
    dt_run = None
    if args.dt_factor != 1.0:
        dt_run = args.dt_factor * model.max_dt()
        print(f"  timestep used  : {dt_run:.1f} s (x{args.dt_factor:g})\n")
    out_dtype = _output_mode(args, np.shape(to_numpy(model.pi)))[0]
    snaps = run_forecast(model, driver, relax, args.hours * 3600, dt=dt_run,
                         output_every=args.output_every * 3600,
                         deadline_s=(None if args.deadline_min is None
                                     else 60.0 * args.deadline_min),
                         info=info,
                         snapshot_dtype=(None if out_dtype == np.float64 else out_dtype))
    if not snaps:
        print(f"\nNo forecast hour completed ({info.get('stopped')}); "
              f"nothing written.")
        return 1
    return _write_forecast(args, run_dir, snaps, info, lev, terrain, meta0, source)


# OUTPUT SIZE (CAM stage S5e). At 3 km x 40 levels a 24 h run is 3.4 GB of
# float64 fields; np.savez_compressed (zlib level 6, one core) took ~7 min of
# the 90 min cycle in test S5c2. For grids above LARGE_GRID columns, "auto"
# writes float32 (theta to ~2e-5 K, pi to ~0.01 Pa: far below what ASOS
# verification resolves) with zlib level 1. Smaller grids (the 12 km run)
# keep float64 and np.savez_compressed, unchanged.
LARGE_GRID = 100_000


def _output_mode(args, shape):
    """(dtype, zlib level or None for np.savez_compressed) for the output."""
    want = getattr(args, "output_dtype", "auto") or "auto"
    large = int(np.prod(shape[-2:])) > LARGE_GRID
    if want == "auto":
        want = "f32" if large else "f64"
    level = getattr(args, "output_level", None)
    if level is None and want == "f32":
        level = 1
    return (np.float32 if want == "f32" else np.float64), level


def _savez_level(path, level, **arrays):
    """np.savez with a chosen zlib level (0 = stored, no compression)."""
    import zipfile
    comp = zipfile.ZIP_STORED if level == 0 else zipfile.ZIP_DEFLATED
    kw = {} if level == 0 else {"compresslevel": int(level)}
    with zipfile.ZipFile(path, "w", compression=comp, allowZip64=True, **kw) as zf:
        for name, a in arrays.items():
            with zf.open(name + ".npy", "w", force_zip64=True) as f:
                np.lib.format.write_array(f, np.asanyarray(a), allow_pickle=True)


def _write_forecast(args, run_dir, snaps, info, lev, terrain, meta0, source):
    """Write the snapshots in the format verify.py and make_maps.py read."""
    out = Path(args.out or (run_dir / "forecast.npz"))
    dtype, level = _output_mode(args, np.shape(snaps[0][1]))
    t_w = time.time()
    save = np.savez_compressed if level is None else \
        (lambda path, **kw: _savez_level(path, level, **kw))
    save(
        out,
        times_s=np.array([s[0] for s in snaps]),
        u=np.stack([s[1] for s in snaps]).astype(dtype, copy=False),
        v=np.stack([s[2] for s in snaps]).astype(dtype, copy=False),
        theta=np.stack([s[3] for s in snaps]).astype(dtype, copy=False),
        pi=np.stack([s[4] for s in snaps]).astype(dtype, copy=False),
        output_dtype=np.array(np.dtype(dtype).name),
        sigma=lev.sigma,
        p_top=lev.p_top,
        terrain=terrain,
        lat=meta0.get("lat"), lon=meta0.get("lon"),
        run_time=meta0.get("run_time", np.array("")),
        source=np.array(source),
        hours_requested=args.hours,
        stopped=np.array(info.get("stopped", "")),
        wall_s=info.get("wall_s", np.nan),
        **({"tg": np.stack(info["tg"])} if info.get("tg") else {}),
    )
    print(f"\nWrote {len(snaps)} snapshots -> {out}  ({info.get('stopped')}, "
          f"{info.get('wall_s', float('nan'))/60:.1f} min)")
    print(f"  output         : {np.dtype(dtype).name}, "
          f"{'zlib 6 (savez_compressed)' if level is None else f'zlib {level}'}, "
          f"{out.stat().st_size / 1e6:.0f} MB in {time.time() - t_w:.0f} s")
    print("Next: verify once the forecast window has closed "
          "(src/verify_pending.py).")
    stopped = info.get("stopped", "")
    # Output is written in every case; the code says how the run ended, so the
    # cycle script's log shows a cut-short or diverged run as such.
    return 0 if stopped in ("completed", "persistence") else \
        (2 if stopped.startswith("deadline") else 3)


if __name__ == "__main__":
    raise SystemExit(main())

#!/usr/bin/env python3
"""
Generate docs/STRUCTURE.md -- an annotated tree of the project.

The tree SHAPE is read from disk so it cannot drift out of date. The one-line
annotations live in the table below and are the part a human maintains; a file
on disk with no annotation is printed with a "(unannotated)" marker so new
files are visible rather than silently blending in.

Run from the project root:  python tools/tree.py
"""

import os
import sys

SKIP_DIRS = {"__pycache__", ".git", ".idea", ".vscode", "herbie_cache"}
SKIP_EXT = {".pyc", ".npy", ".npz", ".log"}

NOTES = {
    # --- top level ---------------------------------------------------------
    "CLAUDE.md": "READ FIRST -- brief for a fresh session: constraints, method, state",
    "CLAUDE_SCIENCE.md": "how to bring the project into Claude Science",
    "config.py": "domain, channels, pressure levels, env-driven paths",
    "resources.py": "CPU governor: 50% ceiling, adapts to other users' load",
    "netpolicy.py": "token-bucket rate limiter, download cache, polite fetcher",
    "preflight.py": "GO/NO-GO checks before a run touches the network",
    "diagnose_herbie.py": "five-section diagnostic for HRRR fetch failures",
    "test_env.py": "environment probe: what is installed on the server",
    "test_netpolicy.py": "bandwidth policy suite (9/9)",
    "requirements.txt": "reference only -- nothing may be installed on the server",

    # --- src ---------------------------------------------------------------
    "src/ingest_hrrr.py": "Herbie fetch, domain cut, stride coarsening -> .npz",
    "src/forecast.py": "end-to-end driver, sigma core over real terrain",
    "src/verify.py": "verify a forecast against observations and archive the pairs",
    "src/test_forecast.py": "driver suite (11/11)",
    "src/test_hrrr_search.py": "GRIB search-string suite (6/6)",
    "src/nwp_emulator_3d.py": "Conv3d state-to-state emulator  [SUPERSEDED -- bounded by its teacher]",
    "src/autoregressive_dataset.py": "training pairs for the emulator; numeric f-hour sort  [SUPERSEDED]",
    "src/train_autoregressive.py": "emulator training loop  [SUPERSEDED by the physics core]",

    # --- dynamics ----------------------------------------------------------
    "src/dynamics/grid.py": "Arakawa C-grid, staggering, adjoint averaging, edge modes",
    "src/dynamics/shallow_water.py": "2D vector-invariant and advective forms (8/8)",
    "src/dynamics/boundaries.py": "Davies relaxation, limited-area driver (6/6)",
    "src/dynamics/vertical.py": "pressure-coordinate vertical operators",
    "src/dynamics/primitive3d.py": "pressure-coordinate 3D core  [SUPERSEDED by sigma]",
    "src/dynamics/sigma.py": "terrain-following coordinate: continuity, PGF, slopes (7/7)",
    "src/dynamics/primitive_sigma.py": "THE CORE -- sigma primitive equations, prognostic p_s",
    "src/dynamics/subgrid.py": "hyperdiffusion, SPPT, Helmholtz balancing (7/7)",
    "src/dynamics/turbulence.py": "Richardson-number vertical mixing",
    "src/dynamics/surface.py": "bulk aerodynamic drag, log law, Louis stability",
    "src/dynamics/convection.py": "dry convective adjustment (post-step, conservative)",
    "src/dynamics/radiation.py": "radiative upper boundary: wave flux through the lid",
    "src/dynamics/radiation_probe.py": "where a reflected wave shows up, rigid vs radiative",
    "src/dynamics/radiation_vs_sponge.py": "development and terrain survival, both boundaries",
    "src/dynamics/radiation_final.py": "the same, with transient-only radiation",
    "src/dynamics/radiation_dt.py": "is the terrain failure a timestep problem? (no)",
    "src/dynamics/sponge_divergent.py": "damp only the divergent wind? (worse)",
    "src/dynamics/sponge_depth_weather.py": "development against sponge depth",
    "src/dynamics/sponge_rate_ladder.py": "development against sponge rate",
    "src/dynamics/sponge_running_ref.py": "frozen vs running sponge reference",
    "src/dynamics/sponge_lid_weather.py": "development against lid height",
    "src/dynamics/sponge_levels_lid.py": "levels x lid: can more levels buy both?",
    "src/dynamics/sponge_curve.py": "eddy energy curves -- suppression or saturation?",
    "src/dynamics/sponge_decision.py": "12 h survival against sponge depth",
    "src/dynamics/sponge_minimum.py": "smallest sponge that survives both cases",
    "src/dynamics/sponge_three.py": "development for shallow sponges",
    "src/dynamics/sponge_edge_test.py": "does growth track the sponge base?",
    "src/dynamics/endgame_probe.py": "hour-by-hour watch of the tall-terrain failure",
    "src/dynamics/endgame_convection.py": "the same, with and without convection",
    "src/dynamics/sponge_edge_test.py": "does the growth peak track the sponge base?",
    "src/dynamics/lid_test.py": "lid height and sponge depth against survival",
    "src/dynamics/rest_terrain_test.py": "motionless atmosphere over tall terrain",
    "src/dynamics/terrain_probe.py": "terrain baseline with and without the filter",
    "src/dynamics/terrain_matrix2.py": "terrain rows re-measured with convection",
    "src/dynamics/probe_4000.py": "hour-by-hour watch at 2500 m against 4000 m",
    "src/dynamics/balance_check.py": "initial-state tendencies and Nh/U by terrain height",
    "src/dynamics/dt_check_4000.py": "timestep sensitivity at 4000 m",
    "src/dynamics/dt_probe_4000.py": "the same, instrumented hour by hour",
    "src/dynamics/kmax_ladder.py": "eddy-diffusivity ceiling sensitivity",
    "src/dynamics/initialization.py": "spectral lowpass -- removes unresolved initial variance",
    "src/dynamics/interpolate.py": "pressure levels -> sigma levels over terrain",
    "src/dynamics/diagnose_growth.py": "energy budget by term / level / wavenumber",
    "src/dynamics/probe_failure.py": "step-by-step failure probe; locates the growing mode",
    "src/dynamics/probe_p50.py": "locates where the radiative lid fails over terrain",
    "src/dynamics/kmax_binding.py": "P-40: does the eddy-diffusivity ceiling bind, and how",
    "src/dynamics/kmax_production.py": "P-40: the same ceiling on the 2500 m production case",
    "src/dynamics/probe_shock.py": "measures the geopotential error the conversion introduces",
    "src/dynamics/probe_shock_consistent.py": "the same, on a hydrostatically self-consistent analysis",
    "src/dynamics/probe_4000.py": "hour-by-hour watch of the 4000 m failure",
    "src/dynamics/balance_check.py": "initial-state balance and Nh/U by terrain height",
    "src/dynamics/dt_check_4000.py": "timestep against survival at 4000 m",
    "src/dynamics/sweep_boundary_layer.py": "mixing x drag x terrain x noise matrix",
    "src/dynamics/noise_ladder.py": "initial-noise amplitude threshold",
    "src/dynamics/order_test.py": "filter/balance ordering experiment",
    "src/dynamics/visualize_instability.py": "cross-sections and growth curves",
    "src/dynamics/instability_cross_section.png": "figure: w and theta' through the mountain",
    "src/dynamics/instability_growth.png": "figure: max|u| against forecast hour",

    # --- verification / postproc -------------------------------------------
    "src/verification/observations.py": "ASOS / mesonet / raob record types and QC",
    "src/verification/fetchers.py": "IEM and MRMS clients (observations only, never HRRR)",
    "src/verification/obs_operator.py": "model state -> observation space (pressure levels)",
    "src/verification/sigma_operator.py": "observation operator for a sigma forecast",
    "src/verification/scoring.py": "bias, RMSE, skill scores",
    "src/postproc/bias_correction.py": "Kalman-filter bias correction, MOS (7/7)",

    # --- docs --------------------------------------------------------------
    "docs/RESEARCH_LOG.md": "dated entries: hypothesis, method, result, interpretation",
    "docs/METHODOLOGY.md": "how claims are established in this project",
    "docs/AI_COLLABORATION.md": "defect taxonomy A-F for the AI-to-build study",
    "docs/PROMPT_LOG.md": "every human prompt, classified -- the study's input record",
    "docs/TOKEN_COST.md": "what the project costs -- billed vs API-equivalent",
    "docs/PROBLEMS.md": "problem register: what is wrong, what fixed it, what ruled it out",
    "docs/REFERENCES.md": "outside work used, in AMS format, and where each is used",
    "docs/CAPABILITIES.md": "what the model can and cannot do, stated up front",
    "docs/STABILITY.md": "the stability investigation  [conclusion superseded]",
    "docs/DATA_ASSIMILATION.md": "observation ingest and analysis design",
    "docs/POSTPROCESSING.md": "neural post-processing design",
    "docs/STRUCTURE.md": "this file, generated by tools/tree.py",
    "docs/LEARNING_LOG.md": "lessons, and the later moments where each one fired again",
    "docs/MANIFEST.txt": "sha256 and size of every file, generated by tools/manifest.py",
    "docs/instability_cross_section.png": "figure: w and theta' through the mountain",
    "docs/instability_growth.png": "figure: max|u| against forecast hour",

    "tools/newlog.py": "append a dated research-log entry from the template",
    "tools/tree.py": "generates docs/STRUCTURE.md",
    "tools/problem.py": "adds to and audits docs/PROBLEMS.md",
    "tools/daily.sh": "one forecast cycle from cron (obs -> analysis -> forecast); `verify` mode scores closed windows",
    "tools/bench_threads.py": "numpy vs torch thread scaling on model-sized arrays, before any port",
    "src/ingest_obs.py": "one cycle's initial state from observations at or before the cycle time",
    "src/make_maps.py": "renders a run's product maps and its HTML viewer (<rundir>/maps/index.html)",
    "src/maps/__init__.py": "forecast maps package",
    "src/maps/geography.py": "Lambert conformal projection in NumPy; Natural Earth state/coast lines, fetched once and cached",
    "src/maps/derive.py": "map diagnostics from the sigma state or the analysis: heights, pressure levels, MSLP, vorticity",
    "src/maps/render.py": "the product maps (surface, upper air, analysis with reports, forecast-minus-observed)",
    "src/maps/viewer.py": "self-contained Pivotal-style HTML viewer: hover readout, click-for-sounding (skew-T, hodograph)",
    "src/maps/test_maps.py": "map tests: projection, clipping, standard-atmosphere diagnostics, rendering, viewer",
    "tools/fetch_boundaries.py": "fetches the Natural Earth lines ahead of time, or writes the bundled copy",
    "src/dynamics/backend.py": "array backends for the core: NumPy (default) or PyTorch (multi-threaded CPU, float64)",
    "src/dynamics/backend_reference.py": "the realistic reference integration used to compare backends",
    "src/dynamics/test_backend.py": "backend tests: torch reproduces numpy to round-off; no cache leakage",
    "tools/mode_structure.py": "where the fastest-growing mode lives (difference of two round-off-different runs) and the initial stability there",
    "tools/mode_budget.py": "which tendency term feeds a growing mode (kinetic-energy budget of the round-off difference)",
    "tools/check_refs.py": "every author-year citation has an AMS-format entry in docs/REFERENCES.md",
    "src/dynamics/test_div_damping.py": "divergence damping: rotational flow untouched, 2dx decay rate, torch = numpy, default bit-identical",
    "tools/score_by_lead.py": "verification RMSE and bias by variable and forecast hour, two archives side by side",
    "tools/compare_forecasts.py": "hour-by-hour difference between two forecast files (e.g. numpy vs torch)",
    "tools/check_backend.py": "on a new machine: torch vs numpy speed and agreement at several thread counts",
    "src/verify_pending.py": "verifies every archived forecast whose window has closed, once",
    "src/analysis/sources.py": "one adapter per observation source; missing sources skipped and logged",
    "src/analysis/build.py": "first guess, sounding superobs, Barnes increments, hydrostatic heights",
    "src/analysis/barnes.py": "successive-correction analysis of increments against a first guess",
    "src/analysis/geo.py": "the forecast's grid, bilinear sampling, ETOPO terrain via ERDDAP",
    "src/analysis/probe_obs_blowup.py": "P-56: where the first observation-built forecast dies",
    "src/analysis/probe_terrain_b.py": "P-56 test B: same observations over HRRR terrain (diagnostic only)",
    "src/analysis/testdata/asos_2026092112_sample.csv": "live IEM ASOS payload, every 10th row",
    "src/analysis/testdata/ndbc_41025_5day_sample.txt": "live NDBC 5-day file, head",
    "src/analysis/testdata/ndbc_active_sample.xml": "live NDBC station list, 25 stations",
    "src/analysis/testdata/raob_KIAD_2026092112.csv": "live IEM sounding, Sterling VA",
    "src/analysis/testdata/raob_KOKX_2026092112.csv": "live IEM sounding, Upton NY",
    "src/analysis/testdata/raob_network.geojson": "live IEM RAOB station table, analysis box",
    "tools/checklayout.py": "checks for src/src nesting, missing and duplicate modules",
    "tools/pull.sh": "update from GitHub over curl -- no git needed on the server",
    "tools/locate_growth.py": "where a saved forecast starts to run away: largest change per snapshot, edge distance",
    "tools/manifest.py": "writes and checks docs/MANIFEST.txt, file by file",
    "tools/tokens.py": "token ledger: billed cost vs API-equivalent shadow price",
    "tools/stale.py": "flags measurements whose file moved after the number was taken",
    "tools/apply_sync.py": "applies a sync archive safely -- no nesting, never touches data/",
    "docs/token_ledger.csv": "per-session token counts; measured and estimated kept apart",

    "skills/nwp-debug/SKILL.md": "skill: diagnose before patching",
    "skills/nwp-record-session/SKILL.md": "skill: research log, register, prompts, tokens",
    "skills/nwp-sync/SKILL.md": "skill: patch route, pull.sh, apply_sync",
}

DIR_NOTES = {
    "skills": "recurring procedures as SKILL.md -- how habits survive a new session",
    "src/dynamics": "the model itself",
    "src/verification": "observations and scoring -- never model output",
    "src/postproc": "learned correction of a finished forecast",
    "docs": "the research record",
    "tools": "maintenance scripts",
    "data": "run outputs and caches (not in version control)",
}


def walk(root, prefix="", rel=""):
    lines = []
    try:
        names = sorted(os.listdir(os.path.join(root, rel)))
    except OSError:
        return lines
    entries = []
    for n in names:
        r = os.path.join(rel, n) if rel else n
        full = os.path.join(root, r)
        if os.path.isdir(full):
            if n in SKIP_DIRS:
                continue
            entries.append((n, r, True))
        else:
            if os.path.splitext(n)[1] in SKIP_EXT or n.startswith("."):
                continue
            if n == "__init__.py":
                continue
            entries.append((n, r, False))
    entries.sort(key=lambda e: (not e[2], e[0]))

    # data/ is deliberately not in version control, so a fresh checkout does
    # not have it and the walk would silently drop it -- and then the layout
    # doc would stop mentioning the one directory in this project that must
    # never be overwritten (P-07). It is part of the documented layout whether
    # or not this particular copy has run anything yet.
    if rel == "" and not os.path.isdir(os.path.join(root, "data")):
        lines.append("|-- data/" + " " * 30 + "# "
                     + DIR_NOTES["data"] + " [not in this copy]")

    for i, (n, r, isdir) in enumerate(entries):
        last = i == len(entries) - 1
        branch = "`-- " if last else "|-- "
        key = r.replace(os.sep, "/")
        if isdir:
            note = DIR_NOTES.get(key, "")
            lines.append(f"{prefix}{branch}{n}/"
                         + (f"{' ' * max(1, 34 - len(prefix) - len(n))}# {note}"
                            if note else ""))
            lines += walk(root, prefix + ("    " if last else "|   "), r)
        else:
            note = NOTES.get(key)
            if note is None and n.startswith("test_") and n.endswith(".py"):
                note = f"suite for {n[5:]}"
            if note is None and n == "README.md":
                note = ("orientation for this directory" if rel
                        else "project overview and how to run it")
            if note is None:
                note = "(unannotated)"
            pad = max(1, 38 - len(prefix) - len(n))
            lines.append(f"{prefix}{branch}{n}{' ' * pad}# {note}")
    return lines


HEADER = """# Project structure

Generated by `tools/tree.py` -- re-run it after adding a file. The tree shape
comes from disk; the annotations are maintained in that script. A file showing
`(unannotated)` has been added without being described.

READING ORDER, for someone new to the code:

0. `docs/PROBLEMS.md` -- what is currently wrong, and what has been ruled out
1. `src/dynamics/grid.py` -- the staggering everything else assumes
2. `src/dynamics/shallow_water.py` -- the same operators in 2D, where they are
   easy to check
3. `src/dynamics/sigma.py` -- the vertical coordinate
4. `src/dynamics/primitive_sigma.py` -- the model
5. `docs/RESEARCH_LOG.md` -- why it looks the way it does
6. `docs/LEARNING_LOG.md` -- which of those lessons came back, and where

Every module has a matching `test_*.py`. The suites are the specification:
when a scheme's behaviour is in question, the test that pins it down is more
informative than the code.

```
"""

FOOTER = """```

## Test suites

These counts are maintained by hand in `tools/tree.py`, which means they are
exactly the kind of number that goes stale without anything noticing. The
"last run" column starts on 2026-09-10 and records when a suite was actually
executed. A dash means nobody has re-run it since the column existed: the
count beside it is what the last person to run it wrote down, and how long ago
that was is not recoverable, which is the whole argument for the column.

| suite | covers | status | last run |
|---|---|---|---|
| `test_shallow_water.py` | 2D dynamics | 8/8 | 2026-09-12 |
| `test_boundaries.py` | Davies relaxation, limited area | 6/6 | 2026-09-12 |
| `test_sigma.py` | coordinate, PGF, continuity | 7/7 | 2026-09-12 |
| `test_subgrid.py` | hyperdiffusion, SPPT, balancing | 7/7 | 2026-09-12 |
| `test_surface.py` | drag, log law, Ekman spiral | 6/6 | 2026-09-12 |
| `test_initialization.py` | spectral filter, noise threshold | 5/5 | 2026-09-12 |
| `test_convection.py` | dry convective adjustment | 5/5 | 2026-09-12 |
| `test_interpolate.py` | pressure -> sigma conversion | 8/8 | — |
| `test_radiation.py` | radiative upper boundary | 7/7 | 2026-09-12 |
| `test_primitive_sigma.py` | the 3D core | 7/7 | 2026-09-12 |
| `test_primitive3d.py` | superseded pressure core | 8/8 | — |
| `test_forecast.py` | end-to-end driver | 11/11 | — |
| `test_hrrr_search.py` | GRIB interface | 6/6 | — |
| `test_verification.py` | observation handling | 9/9 | 2026-09-10 |
| `test_sigma_operator.py` | sigma observation operator | 7/7 | — |
| `test_verify.py` | verification archiver | 7/7 | — |
| `test_fetchers.py` | IEM / MRMS clients | 9/9 | 2026-09-10 |
| `test_bias_correction.py` | post-processing | 7/7 | — |
| `test_netpolicy.py` | bandwidth policy | 9/9 | 2026-09-10 |

## Data flow

```
observations  ->  verification/fetchers  ->  observations  ->  sigma_operator
(ASOS, mesonet,        |                       (parse, QC)             |
 raob, MRMS)           v                                              v
                  stored VERBATIM                                  verify.py
                  (the irreplaceable part;                    (match, archive)
                   everything else is                               |
                   recomputable from it)                            v
                                                        scoring / postproc
HRRR analysis ->  ingest_hrrr  ->  interpolate  ->  initialization  ->
(initial state only,  (+ terrain)    (p -> sigma)     (filter, then
 never verification)                                   rebalance)
                                                            |
                                                            v
                                                     primitive_sigma
                                                (+ turbulence, surface,
                                                 convection, subgrid)
                                                                    |
                                                                    v
                                                          postproc/bias_correction
```

The separation on the left is deliberate and is a constraint of the project:
HRRR may seed a forecast, but it may never verify one. Verification comes from
instruments.
"""


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    body = "\n".join(walk(root))
    out = HEADER + "NWP_Deployment_Package/\n" + body + "\n" + FOOTER
    path = os.path.join(root, "docs", "STRUCTURE.md")
    with open(path, "w", encoding="utf-8", newline="\n") as f:  # LF on Windows too
        f.write(out)
    print(out)
    # Count the tree, not the header sentence that explains the marker --
    # this used to report "1 unannotated file(s)" on a fully annotated tree,
    # which is a check that cries wolf every single run.
    n_un = body.count("(unannotated)")
    print(f"\nwrote {path}"
          + (f"  -- {n_un} unannotated file(s)" if n_un else ""), file=sys.stderr)


if __name__ == "__main__":
    main()

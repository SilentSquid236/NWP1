"""
Pre-flight check. Run this BEFORE a live test run.

    python preflight.py                     # check with defaults
    python preflight.py --hours 12 --stride 4

Every check is cheap and read-only. The point is to fail in ten seconds on a
missing package or an unwritable directory rather than twenty minutes into a
download.

Exit code 0 means GO, 1 means something needs fixing first.
"""

import argparse
import importlib
import os
import shutil
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))
sys.path.insert(0, str(HERE / "src"))
sys.path.insert(0, str(HERE / "src" / "dynamics"))

OK, WARN, FAIL = "OK  ", "WARN", "FAIL"
rows = []


def check(name, status, detail):
    rows.append((status, name, detail))
    print(f"  [{status}] {name}\n        {detail}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--hours", type=int, default=12)
    ap.add_argument("--stride", type=int, default=4)
    ap.add_argument("--skip-network", action="store_true")
    ap.add_argument("--skip-tests", action="store_true")
    args = ap.parse_args()

    print("\nNWP pre-flight\n" + "=" * 66)

    # --- 1. Python and packages -------------------------------------------
    v = sys.version_info
    check("python version", OK if v >= (3, 8) else FAIL,
          f"{v.major}.{v.minor}.{v.micro} at {sys.executable}")

    required = ["numpy"]
    ingest_only = ["xarray", "cfgrib", "eccodes", "herbie"]
    missing_req, missing_ing = [], []
    versions = {}
    for m in required + ingest_only:
        try:
            mod = importlib.import_module(m)
            versions[m] = getattr(mod, "__version__", "?")
        except Exception:
            (missing_req if m in required else missing_ing).append(m)

    check("core packages", FAIL if missing_req else OK,
          f"numpy {versions.get('numpy', 'MISSING')}" +
          (f"; MISSING: {missing_req}" if missing_req else ""))
    check("ingestion packages", WARN if missing_ing else OK,
          (f"missing {missing_ing} -- ingest will fail, forecast from existing "
           f".npz still works" if missing_ing else
           ", ".join(f"{m} {versions[m]}" for m in ingest_only)))

    # --- 2. Configuration --------------------------------------------------
    try:
        import config
        d = config.DOMAIN
        span_lat = d["lat_max"] - d["lat_min"]
        span_lon = d["lon_max"] - d["lon_min"]
        sane = (config.N_LEVELS > 3 and 0 < span_lat < 40 and 0 < span_lon < 60)
        check("configuration", OK if sane else FAIL,
              f"{config.N_CHANNELS} vars x {config.N_LEVELS} levels, "
              f"{d['name']} {span_lat:.1f}x{span_lon:.1f} deg, "
              f"levels {config.PRESSURE_LEVELS[0]}-{config.PRESSURE_LEVELS[-1]} hPa")
    except Exception as e:
        check("configuration", FAIL, f"config.py did not import: {e}")
        config = None

    # --- 3. Data directory and disk ---------------------------------------
    if config:
        root = config.DATA_ROOT
        env_set = "NWP_DATA_ROOT" in os.environ
        try:
            root.mkdir(parents=True, exist_ok=True)
            probe = root / ".preflight_write_test"
            probe.write_text("x")
            probe.unlink()
            writable = True
        except Exception as e:
            writable = False

        check("data directory", OK if writable else FAIL,
              f"{root} {'(writable)' if writable else 'NOT WRITABLE'}"
              + ("" if env_set else "  -- NWP_DATA_ROOT unset, using default"))

        if writable:
            free_gb = shutil.disk_usage(root).free / 1e9
            ny, nx = 388 // args.stride, 438 // args.stride
            per_hour_mb = (config.N_CHANNELS * config.N_LEVELS * ny * nx * 4) / 1e6
            need_gb = per_hour_mb * (args.hours + 1) / 1000
            ratio = free_gb / max(need_gb, 1e-9)
            check("disk space", OK if ratio > 20 else (WARN if ratio > 3 else FAIL),
                  f"{free_gb:.1f} GB free; this run needs ~{need_gb*1000:.0f} MB "
                  f"({per_hour_mb:.1f} MB/hour at stride {args.stride})")

    # --- 4. Herbie cache location -----------------------------------------
    # Herbie defaults to ~/data. On a shared box that is a quota'd home
    # directory, and a failed write there does NOT raise -- the file simply
    # never appears and xarray later reports a baffling FileNotFoundError.
    # ingest_hrrr.py overrides save_dir; verify that target is usable.
    if config:
        cache = Path(os.environ.get("NWP_HERBIE_DIR",
                                    config.DATA_ROOT / "herbie_cache"))
        try:
            cache.mkdir(parents=True, exist_ok=True)
            probe = cache / ".preflight_probe"
            probe.write_bytes(b"x" * 4096)
            probe.unlink()
            cache_ok, why = True, ""
        except Exception as e:
            cache_ok, why = False, f"{type(e).__name__}: {e}"

        free_gb = shutil.disk_usage(cache if cache.exists()
                                    else cache.parent).free / 1e9
        home_free = shutil.disk_usage(Path.home()).free / 1e9
        check("Herbie cache", OK if cache_ok else FAIL,
              f"{cache} "
              + (f"writable, {free_gb:.1f} GB free" if cache_ok
                 else f"NOT WRITABLE -- {why}")
              + f"  (home has {home_free:.1f} GB; Herbie's default ~/data is "
                f"deliberately overridden)")

    # --- 5. Compute budget -------------------------------------------------
    try:
        import resources
        plan = resources.plan()
        load = resources.load_average()
        others = max(0.0, load - 1)
        check("compute budget", OK,
              f"{plan['total_cores']} cores, 1-min load {load:.1f} "
              f"(~{others:.0f} used by others); ceiling "
              f"{plan['ceiling_cores']} cores. NOTE: dynamics is "
              f"single-threaded numpy -- extra cores do not speed it up")
    except Exception as e:
        check("compute budget", WARN, f"resources.py: {e}")

    # --- 6. Estimated runtime ---------------------------------------------
    ms_per_step = {1: 4040.0, 2: 708.0, 4: 138.0}.get(args.stride)
    if ms_per_step:
        dt = {1: 11.6, 2: 23.1, 4: 46.3}[args.stride]
        steps = args.hours * 3600 / dt
        minutes = ms_per_step * steps / 1000 / 60
        check("estimated runtime", OK if minutes < 90 else WARN,
              f"{args.hours} h forecast at ~{3*args.stride} km: "
              f"{steps:.0f} steps x {ms_per_step:.0f} ms = "
              f"{minutes:.1f} min (single core)"
              + ("  -- consider a larger --stride" if minutes > 90 else ""))

    # --- 7. Network --------------------------------------------------------
    if not args.skip_network:
        try:
            from netpolicy import PoliteFetcher, estimate_ingest_mb, max_mbps
            f = PoliteFetcher()
            targets = [
                ("AWS HRRR", "https://noaa-hrrr-bdp-pds.s3.amazonaws.com/?list-type=2&max-keys=1"),
                ("IEM observations", "https://mesonet.agron.iastate.edu/"),
            ]
            for label, url in targets:
                t0 = time.time()
                try:
                    f.get(url, timeout=20, use_cache=False)
                    check(f"reachable: {label}", OK,
                          f"responded in {time.time()-t0:.2f}s")
                except Exception as e:
                    check(f"reachable: {label}", FAIL,
                          f"{type(e).__name__}: {str(e)[:90]}")

            est = estimate_ingest_mb(args.hours + 1, config.N_LEVELS if config else 20)
            check("bandwidth plan", OK,
                  f"~{est['total_download_MB']:.0f} MB total, capped at "
                  f"{max_mbps():.0f} MB/s (NWP_MAX_MBPS to raise)")
        except Exception as e:
            check("network checks", WARN, f"skipped: {e}")

    # --- 8. Fast test suites ----------------------------------------------
    if not args.skip_tests:
        suites = [("dynamics", "src/dynamics/test_primitive3d.py"),
                  ("verification", "src/verification/test_verification.py"),
                  ("forecast", "src/test_forecast.py")]
        for label, rel in suites:
            path = HERE / rel
            if not path.exists():
                check(f"tests: {label}", WARN, f"{rel} not found")
                continue
            t0 = time.time()
            r = subprocess.run([sys.executable, path.name], cwd=path.parent,
                               capture_output=True, text=True, timeout=600)
            tail = [l for l in r.stdout.strip().splitlines() if "passed" in l]
            check(f"tests: {label}", OK if r.returncode == 0 else FAIL,
                  f"{tail[-1] if tail else 'no result'} "
                  f"({time.time()-t0:.1f}s)")

    # --- verdict -----------------------------------------------------------
    fails = [r for r in rows if r[0] == FAIL]
    warns = [r for r in rows if r[0] == WARN]
    print("=" * 66)
    if fails:
        print(f"NO-GO -- {len(fails)} blocking issue(s):")
        for _, name, _ in fails:
            print(f"    - {name}")
        return 1

    print(f"GO" + (f" -- with {len(warns)} warning(s)" if warns else ""))
    for _, name, _ in warns:
        print(f"    ! {name}")

    if config:
        run = f"analysis_{time.strftime('%Y%m%d')}_00"
        print(f"\nSuggested first run:")
        print(f"    python src/ingest_hrrr.py --start <YYYY-MM-DD>T00 "
              f"--hours {args.hours + 1} --stride {args.stride}")
        print(f"    python src/forecast.py --run-dir "
              f"{config.TENSOR_DIR}/{run} --hours {args.hours}")
        print(f"\nIngest {args.hours + 1} hours for a {args.hours} h forecast: "
              f"the extra frame drives the boundaries at the final hour.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

"""
Environment smoke test. Run this first on any new machine.

    python test_env.py

Reports which parts of the stack are present. A missing GRIB library only
matters on the machine that ingests data; a training-only box needs torch.
"""

import importlib
import sys

CORE = ["torch", "numpy"]
DATA = ["xarray", "cfgrib", "eccodes", "herbie", "metpy"]


def check(name):
    try:
        mod = importlib.import_module(name)
        return True, getattr(mod, "__version__", "?")
    except Exception as e:
        return False, f"{type(e).__name__}: {e}"


def main():
    print("NWP environment check")
    print(f"  python   : {sys.version.split()[0]}  ({sys.platform})")
    print(f"  executable: {sys.executable}\n")

    ok_core = True
    print("Training stack (required to train):")
    for m in CORE:
        ok, info = check(m)
        print(f"  {'OK  ' if ok else 'FAIL'} {m:10s} {info}")
        ok_core &= ok

    print("\nIngestion stack (required only to build tensors from GRIB):")
    ok_data = True
    for m in DATA:
        ok, info = check(m)
        print(f"  {'OK  ' if ok else '--  '} {m:10s} {info}")
        ok_data &= ok

    print()
    try:
        import config
        print("Configuration:")
        print(config.describe())
    except Exception as e:
        print(f"config.py not importable from here: {e}")

    print()
    if ok_core and ok_data:
        print("Full pipeline available: this machine can ingest AND train.")
    elif ok_core:
        print("Training only: ingest elsewhere, copy tensors here, then train.")
    else:
        print("Cannot train here -- torch is missing or broken.")
    return 0 if ok_core else 1


if __name__ == "__main__":
    raise SystemExit(main())

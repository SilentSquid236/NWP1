"""
Diagnose a Herbie download failure. Read-only apart from one small test file.

    python diagnose_herbie.py
"""
import os, shutil, subprocess, sys, traceback
from datetime import datetime, timedelta
from pathlib import Path

print("=" * 66)
print("1. Where does Herbie want to write, and can it?")
print("=" * 66)
try:
    from herbie import Herbie, config as hconfig
    print(f"  herbie version : {__import__('herbie').__version__}")
    save_dir = Path(hconfig.get('default', {}).get('save_dir', Path.home() / 'data'))
    print(f"  save_dir       : {save_dir}")
    print(f"  config file    : {Path.home() / '.config' / 'herbie' / 'config.toml'}")
except Exception as e:
    print(f"  herbie import FAILED: {e}")
    raise SystemExit(1)

home = Path.home()
du = shutil.disk_usage(home)
print(f"  home free      : {du.free/1e9:.2f} GB of {du.total/1e9:.1f} GB")

try:
    save_dir.mkdir(parents=True, exist_ok=True)
    probe = save_dir / ".write_probe"
    probe.write_bytes(b"x" * 1024)
    probe.unlink()
    print(f"  writable       : YES")
except Exception as e:
    print(f"  writable       : NO -- {type(e).__name__}: {e}")

print()
r = subprocess.run(["quota", "-s"], capture_output=True, text=True)
if r.returncode == 0 and r.stdout.strip():
    print("  quota:")
    for line in r.stdout.strip().splitlines():
        print(f"    {line}")
else:
    print("  quota         : no per-user quota reported")

print()
print("=" * 66)
print("2. Can Herbie find the file remotely?")
print("=" * 66)
when = datetime.utcnow() - timedelta(hours=12)
when = when.replace(minute=0, second=0, microsecond=0)
print(f"  testing        : {when:%Y-%m-%d %H}Z  (12 h ago, definitely archived)")

try:
    H = Herbie(when.strftime("%Y-%m-%d %H:%M"), model="hrrr", product="prs",
               fxx=0, verbose=True)
    print(f"  found source   : {getattr(H, 'grib_source', None)}")
    print(f"  remote URL     : {getattr(H, 'grib', None)}")
    print(f"  index URL      : {getattr(H, 'idx', None)}")
    print(f"  local path     : {getattr(H, 'get_localFilePath', lambda *a: '?')()}")
except Exception as e:
    print(f"  Herbie() FAILED: {type(e).__name__}: {e}")
    traceback.print_exc()
    raise SystemExit(1)

print()
print("=" * 66)
print("3. Does the variable search match anything?")
print("=" * 66)
search = r"^(?:TMP|RH|UGRD|VGRD|HGT):\d+ mb:"
try:
    idx = H.inventory(search)
    print(f"  regex          : {search}")
    print(f"  messages match : {len(idx)}")
    if len(idx):
        print(f"  example rows:")
        for _, row in idx.head(3).iterrows():
            print(f"    {row.get('search_this', row.to_dict())}")
    else:
        print("  NOTHING MATCHED -- the regex does not fit this file's naming.")
        allidx = H.inventory()
        print(f"  file has {len(allidx)} messages; first few:")
        for _, row in allidx.head(8).iterrows():
            print(f"    {row.get('search_this', '')}")
except Exception as e:
    print(f"  inventory FAILED: {type(e).__name__}: {e}")
    traceback.print_exc()

print()
print("=" * 66)
print("4. Can it actually download a small subset?")
print("=" * 66)
try:
    ds = H.xarray(":TMP:500 mb:")
    print(f"  SUCCESS -- got {type(ds).__name__}")
    if hasattr(ds, "dims"):
        print(f"  dims           : {dict(ds.dims)}")
except TypeError:
    try:
        ds = H.xarray(searchString=":TMP:500 mb:")
        print("  SUCCESS via legacy searchString= argument")
        print(f"  dims           : {dict(ds.dims)}")
    except Exception as e:
        print(f"  FAILED: {type(e).__name__}: {e}")
        traceback.print_exc()
except Exception as e:
    print(f"  FAILED: {type(e).__name__}: {e}")
    traceback.print_exc()

print()
print("=" * 66)
print("5. Tools Herbie shells out to")
print("=" * 66)
for tool in ("curl", "wget"):
    print(f"  {tool:6s}: {shutil.which(tool) or 'NOT FOUND'}")

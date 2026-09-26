#!/usr/bin/env python3
"""Fetch the Natural Earth state and coast lines once and cache them.

    python tools/fetch_boundaries.py            # -> <data>/static/boundaries_ne50m.npz
    python tools/fetch_boundaries.py --bundle   # -> src/maps/boundaries_ne50m.npz

make_maps.py fetches on first use anyway. This is for fetching ahead of time,
or for writing the bundled copy that machines without network use.
"""
import argparse
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))
sys.path.insert(0, str(ROOT / "src"))

import config                                            # noqa: E402
from netpolicy import PoliteFetcher                      # noqa: E402
from maps.geography import (CACHE_NAME, build_boundaries,  # noqa: E402
                            save_boundaries)


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--bundle", action="store_true",
                   help="write src/maps/boundaries_ne50m.npz instead of the cache")
    a = p.parse_args()
    out = (ROOT / "src" / "maps" / CACHE_NAME if a.bundle
           else Path(config.DATA_ROOT) / "static" / CACHE_NAME)
    f = PoliteFetcher()
    layers = build_boundaries(config.DOMAIN, lambda u: f.get_text(u, timeout=300))
    save_boundaries(out, layers)
    for k, (x, _) in layers.items():
        print(f"  {k:10s} {int((~(x != x)).sum()):7d} points")
    print(f"wrote {out}  ({out.stat().st_size / 1e3:.0f} kB)")


if __name__ == "__main__":
    main()

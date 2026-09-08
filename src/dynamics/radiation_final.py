"""Transient-only radiation: the two measurements that decide it."""
import numpy as np
np.seterr(all="ignore")
from radiation_vs_sponge import curve, terrain

hrs = list(range(6, 49, 6))
configs = [("sponge 5, rigid", 5, False),
           ("no sponge, rigid", 0, False),
           ("no sponge, RADIATIVE", 0, True)]

print("baroclinic development\n")
print("config                | " + " ".join(f"{h:>8}h" for h in hrs) + " | 48h/6h")
for name, sp, rad in configs:
    c = curve(sp, rad)
    c += [float("nan")] * (len(hrs) - len(c))
    r = c[-1] / c[0] if np.isfinite(c[-1]) and c[0] else float("nan")
    print(f"{name:21} | " + " ".join(f"{v:9.2e}" for v in c) + f" | {r:6.2f}",
          flush=True)

print("\n2500 m terrain, 12 h\n")
print(f"{'config':21} {'survived':>9} {'max|u|':>8}")
for name, sp, rad in configs:
    d, u = terrain(sp, rad)
    print(f"{name:21} {d:6d}/12 {u:8.1f}", flush=True)

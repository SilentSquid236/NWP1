"""What does a 3-level sponge cost in development? It survives 12/12."""
import numpy as np
np.seterr(all="ignore")
from sponge_depth_weather import growth
print("baroclinic development, shallow sponges\n")
print(f"{'sponge':>7} {'eddy x/day':>11} {'max|v|':>8}")
for nsp in (0, 2, 3, 4, 5):
    g, vmax, jet = growth(nsp)
    print(f"{nsp:7d} {g:11.2f} {vmax:8.1f}", flush=True)

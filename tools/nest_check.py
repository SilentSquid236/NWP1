"""
How a 3 km nested forecast sits against the 12 km forecast that drives it
(CAM stage S5).

    python tools/nest_check.py <3 km forecast.npz> <12 km forecast.npz>

For each output hour: the 12 km state interpolated to the 3 km grid exactly as
the boundary driver does it (nest.frame_from_coarse), and the rms difference
of theta, u and v and of surface pressure in three bands -- the outer 5 cells
(where the relaxation is strongest), the rest of the relaxation zone, and the
interior. Also the largest wind speed and the theta range, to see whether the
3 km run stays physical.
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
import nest  # noqa: E402


def bands(ny, nx, width=15):
    j, i = np.meshgrid(np.arange(ny), np.arange(nx), indexing="ij")
    d = np.minimum(np.minimum(j, ny - 1 - j), np.minimum(i, nx - 1 - i))
    return {"edge<5": d < 5, "zone5-15": (d >= 5) & (d < width), "interior": d >= width}


def main(fine, coarse):
    a = np.load(fine, allow_pickle=False)
    b = np.load(coarse, allow_pickle=False)
    h3 = np.asarray(a["terrain"], float)
    tb = np.asarray(b["times_s"], float)
    bm = bands(*h3.shape)
    print(f"3 km {h3.shape}, {len(a['times_s'])} snapshots; 12 km {b['u'].shape[-2:]}")
    print("hour  band        d_theta(K)  d_u(m/s)  d_v(m/s)  d_ps(hPa)   max|U|3km  theta 3km")
    for n, t in enumerate(np.asarray(a["times_s"], float)):
        if t % 21600 and n != len(a["times_s"]) - 1:
            continue
        m = int(np.argmin(np.abs(tb - t)))
        if abs(tb[m] - t) > 1.0:
            continue
        uu, vv, th, pi = nest.frame_from_coarse(b["u"][m], b["v"][m], b["theta"][m], b["pi"][m],
                                                b["sigma"], float(b["p_top"]),
                                                np.asarray(b["terrain"], float), h3, a["sigma"])
        spd = float(np.hypot(a["u"][n], a["v"][n]).max())
        for k, msk in bm.items():
            r = lambda x, y: float(np.sqrt(np.mean((x - y)[:, msk] ** 2)))
            dps = float(np.sqrt(np.mean((a["pi"][n] - pi)[msk] ** 2))) / 100.0
            print(f"{t / 3600:4.0f}  {k:10s}  {r(a['theta'][n], th):9.3f}  {r(a['u'][n], uu):8.3f}  "
                  f"{r(a['v'][n], vv):8.3f}  {dps:9.3f}   {spd:8.1f}  "
                  f"{float(a['theta'][n].min()):.1f}-{float(a['theta'][n].max()):.1f}")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

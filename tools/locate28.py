
import numpy as np
V = "/data5/pierce/AINWP/data/verification_tests"
for name, f in (("z0x 06Z", f"{V}/aj2/z0x/obs_20260928_06.npz"), ("ctl 06Z", f"{V}/aj2/ctl/obs_20260928_06.npz"),
                ("AE20 00Z (P-67)", f"{V}/ae20/20260928_00.npz")):
    try: z = np.load(f)
    except Exception as e: print(name, "missing", e); continue
    lat, lon = z["lat"], z["lon"]; hrs = np.rint(z["times_s"] / 3600).astype(int)
    print("==", name, "hours", hrs.min(), "-", hrs.max())
    for n in range(len(hrs)):
        if hrs[n] < 12: continue
        s = np.hypot(z["u"][n], z["v"][n]); k, j, i = np.unravel_index(np.argmax(s), s.shape)
        th = z["theta"][n]
        print(f"  {hrs[n]:2d} h max speed {s.max():5.1f} at level {k} j {j} i {i} ({lat[j,i]:.2f}N {lon[j,i]:.2f}) "
              f"edge dist {min(j, i, s.shape[1]-1-j, s.shape[2]-1-i)}; terrain {z['terrain'][j,i]:.0f} m; "
              f"dtheta(k,k+1) there {th[min(k+1,19),j,i]-th[k,j,i]:+.2f}")

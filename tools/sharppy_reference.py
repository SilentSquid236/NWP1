"""SHARPpy's indices for model columns: the reference for src/maps/indices.py.

    python tools/sharppy_reference.py COLUMNS.npz OUT.json      (an environment with SHARPpy)

COLUMNS.npz holds (levels, N) arrays p (Pa), T (K), Td (K), z (m MSL), u, v (m/s) with
level 0 = the surface, and terrain (N,). The output has the same keys as
indices.column_indices, in the same units, for each column (NaN where SHARPpy has none).
SHARPpy: Blumberg et al. (2017), sharppy.sharptab (no GUI needed).
"""
import json
import sys
import warnings

import numpy as np
import sharppy.sharptab.interp as interp
import sharppy.sharptab.params as params
import sharppy.sharptab.profile as profile
import sharppy.sharptab.utils as utils
import sharppy.sharptab.winds as winds

warnings.filterwarnings("ignore")
KT = 1.0 / 0.514444


def f(x):
    """A float, or NaN for SHARPpy's masked or missing values."""
    try:
        x = float(x)
    except (TypeError, ValueError, np.ma.MaskError):
        return float("nan")
    return x if np.isfinite(x) and abs(x) < 1e5 else float("nan")


def one(p, T, Td, z, u, v):
    spd, wdir = np.hypot(u, v) * KT, (np.degrees(np.arctan2(-u, -v)) + 360.0) % 360.0
    prof = profile.create_profile(profile="default", pres=p / 100.0, hght=z, tmpc=T - 273.15,
                                  dwpc=Td - 273.15, wspd=spd, wdir=wdir, missing=-9999, strictQC=False)
    sb, ml, mu = params.parcelx(prof, flag=1), params.parcelx(prof, flag=4), params.parcelx(prof, flag=3)
    agl = lambda h: interp.to_agl(prof, h)
    sfc = prof.pres[prof.sfc]
    p6 = interp.pres(prof, interp.to_msl(prof, 6000.0)); p1 = interp.pres(prof, interp.to_msl(prof, 1000.0))
    s1, s6 = winds.wind_shear(prof, pbot=sfc, ptop=p1), winds.wind_shear(prof, pbot=sfc, ptop=p6)
    rstu, rstv, _, _ = winds.non_parcel_bunkers_motion(prof)   # the method indices.py uses
    return {
        "sbcape": f(sb.bplus), "sbcin": f(sb.bminus), "sblcl": f(sb.lclhght), "sblfc": f(sb.lfchght),
        "sbel": f(sb.elpres), "li": f(sb.li5),
        "mlcape": f(ml.bplus), "mlcin": f(ml.bminus), "mllcl": f(ml.lclhght),
        "mucape": f(mu.bplus), "mucin": f(mu.bminus), "mup": f(mu.pres),
        "pw": f(params.precip_water(prof)) * 25.4,
        "lr03": f(params.lapse_rate(prof, 0.0, 3000.0, pres=False)),
        "lr75": f(params.lapse_rate(prof, 700.0, 500.0, pres=True)),
        "shr1": f(utils.mag(*s1)), "shr6": f(utils.mag(*s6)),
        "rm_u": f(rstu), "rm_v": f(rstv),
        "srh1": f(winds.helicity(prof, 0.0, 1000.0, stu=rstu, stv=rstv)[0]),
        "srh3": f(winds.helicity(prof, 0.0, 3000.0, stu=rstu, stv=rstv)[0]),
    }


def main(src, dst):
    z = np.load(src)
    n = z["p"].shape[1]
    rows = []
    for i in range(n):
        try:
            rows.append(one(z["p"][:, i], z["T"][:, i], z["Td"][:, i], z["z"][:, i], z["u"][:, i], z["v"][:, i]))
        except Exception as e:                      # noqa: BLE001
            rows.append({"error": f"{type(e).__name__}: {e}"})
    keys = sorted({k for r in rows for k in r if k != "error"})
    out = {k: [r.get(k, float("nan")) for r in rows] for k in keys}
    out["errors"] = [r.get("error") for r in rows]
    json.dump(out, open(dst, "w"))
    print(f"{n} columns, {sum(1 for r in rows if 'error' in r)} errors")


if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

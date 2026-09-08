"""
The decisive one: can more levels buy both?

    p_top   sponge base  eddy x/day  reflection
    200 hPa  none            3.18       60.8
    200 hPa  301 hPa         0.85       36.4     <- sponge sits ON the wave
    100 hPa  214 hPa         0.97       53.2
     50 hPa  170 hPa         1.89       60.3     <- wave free, wave not absorbed

Raising the lid moves the sponge off the baroclinic wave and development
recovers, 0.85 -> 1.89. But reflection climbs back to the no-sponge value,
because twenty levels stretched over a deeper domain leave the sponge covering
too little of the wave to absorb it.

Both symptoms have one cause: TWENTY LEVELS IS NOT ENOUGH TO HOLD A
TROPOSPHERE AND A SPONGE AT THE SAME TIME. With p_top = 200 hPa there is no
room above the weather; with p_top = 50 hPa there is room but no resolution
left to fill it.

PREDICTION: a 50 hPa lid with 26-30 levels keeps tropospheric resolution AND
gives the sponge somewhere stratospheric to sit -- development near 3, and
reflection back down toward the 36 that a well-placed sponge achieved.

If reflection stays at 60, the sponge depth is not the issue and the honest
conclusion is that a rigid lid needs a radiative boundary condition, not a
thicker blanket.
"""
import numpy as np
np.seterr(all="ignore")

from sigma import SigmaLevels
import sponge_lid_weather as S

if __name__ == "__main__":
    print("levels x lid, 5-level sponge unless noted\n")
    print(f"{'nz':>4} {'p_top':>8} {'sponge':>7} {'base':>7} "
          f"{'eddy x/day':>11} {'refl':>8} {'dt':>7} {'survived':>9}")
    for nz, p_top, nsp in ((20, 5000.0, 5),
                           (26, 5000.0, 5),
                           (26, 5000.0, 8),
                           (30, 5000.0, 8)):
        lev = SigmaLevels(nz, p_top=p_top)

        import primitive_sigma, lid_test
        from grid import CGrid
        from sigma import P0, KAPPA

        # development
        gr = CGrid(48, 48, 60e3, 60e3, f0=1.0e-4, beta=1.6e-11)
        m = primitive_sigma.PrimitiveSigma(gr, lev, sponge_levels=nsp)
        k_y = 2 * np.pi / gr.Ly
        m.pi = np.full((gr.ny, gr.nx), 101325.0) - lev.p_top
        p = lev.pressure(m.pi)
        m.theta = (258.0 + 6.0 * np.cos(k_y * gr.Yc)) / (p / P0) ** KAPPA
        phi = m.geopotential()
        for k in range(lev.nz):
            m.u[k] = -0.5 * (gr.dy_forward(phi[k])
                             + gr.dy_backward(phi[k])) / gr.f0
        m.theta += 0.5 * np.sin(4 * np.pi * gr.Xc / gr.Lx) * np.sin(k_y * gr.Yc)

        def eddy(mm):
            up = mm.u - mm.u.mean(axis=2, keepdims=True)
            vp = mm.v - mm.v.mean(axis=2, keepdims=True)
            return float((up ** 2 + vp ** 2).sum())

        dt = m.max_dt()
        m.run(24 * 3600, dt=dt); e1 = eddy(m)
        m.run(24 * 3600, dt=dt); e2 = eddy(m)
        g = e2 / e1 if (np.isfinite(m.u).all() and e1 > 0) else float("nan")

        refl, done = S.terrain(nsp, p_top) if nz == 20 else (None, None)
        if refl is None:
            mm = lid_test.build_on(lev, 2500.0, nsp)
            u0 = mm.u.copy()
            d2 = mm.max_dt()
            done = 0
            for _ in range(6):
                mm.run(3600, dt=d2)
                if not np.isfinite(mm.u).all():
                    break
                done += 1
            refl = (float(np.abs(mm.u - u0).max())
                    if np.isfinite(mm.u).all() else float("nan"))

        print(f"{nz:4d} {p_top/100:7.0f}h {nsp:7d} {p[nsp-1,0,0]/100:6.0f}h "
              f"{g:11.2f} {refl:8.2f} {dt:7.1f} {done:6d}/6", flush=True)

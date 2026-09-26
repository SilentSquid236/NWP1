#!/usr/bin/env python3
"""
Tests for divergence damping (subgrid.divergence_damping).
Method: Skamarock and Klemp (1992). Run from the repo root:

    PYTHONPATH=src/dynamics python src/dynamics/test_div_damping.py
"""
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parent))

import backend                                          # noqa: E402
from grid import CGrid                                  # noqa: E402
from subgrid import divergence_damping                  # noqa: E402

NX, NY, DX, DY, NU = 40, 36, 12_000.0, 11_000.0, 5.0e4
IN = (slice(None), slice(3, -3), slice(3, -3))           # interior, away from edges


def grid():
    return CGrid(NX, NY, DX, DY, edge_mode="replicate")


def test_rotational_flow_is_untouched():
    """A streamfunction flow has zero discrete divergence, so zero damping."""
    rng = np.random.default_rng(3)
    psi = rng.standard_normal((4, NY + 1, NX + 1)) * 1e6     # corners
    u = -(psi[:, 1:, :-1] - psi[:, :-1, :-1]) / DY          # west faces
    v = (psi[:, :-1, 1:] - psi[:, :-1, :-1]) / DX           # south faces
    du, dv = divergence_damping(u, v, grid(), NU)
    scale = NU * np.abs(u).max() / DX ** 2
    worst = max(np.abs(du[IN]).max(), np.abs(dv[IN]).max()) / scale
    assert worst < 1e-10, worst
    return f"rotational flow: interior tendency {worst:.1e} of the 2dx scale"


def test_2dx_wave_decays_at_4nu_over_dx2():
    i = np.arange(NX)
    u = np.broadcast_to((-1.0) ** i, (4, NY, NX)).copy()
    v = np.zeros_like(u)
    du, dv = divergence_damping(u, v, grid(), NU)
    rate = -du[IN] / u[IN]
    expect = 4.0 * NU / DX ** 2
    assert np.allclose(rate, expect, rtol=1e-12), (rate.min(), rate.max(), expect)
    assert np.abs(dv[IN]).max() == 0.0
    return f"2dx wave: decay rate {rate.mean():.4e} 1/s, expected {expect:.4e}"


def test_torch_matches_numpy():
    if backend.TORCH is None:
        return "torch not installed here: skipped"
    rng = np.random.default_rng(5)
    u, v = rng.standard_normal((2, 4, NY, NX))
    ref = divergence_damping(u, v, grid(), NU)
    T = backend.TORCH
    got = divergence_damping(T.asarray(u), T.asarray(v), grid(), NU)
    worst = max(float(np.abs(backend.to_numpy(g) - r).max() / np.abs(r).max())
                for g, r in zip(got, ref))
    assert worst < 1e-12, worst
    return f"torch against numpy: {worst:.1e} relative"


def test_off_by_default_and_bit_identical():
    """The core's default has no damping, and tendencies are unchanged by it."""
    from sigma import SigmaLevels
    from primitive_sigma import PrimitiveSigma
    g = CGrid(NX, NY, DX, DY, f0=1e-4, beta=1.6e-11, edge_mode="replicate")
    lev = SigmaLevels(8)
    m = PrimitiveSigma(g, lev)
    assert m.div_damp == 0.0
    rng = np.random.default_rng(7)
    m.u = 10.0 + rng.standard_normal((8, NY, NX))
    m.v = rng.standard_normal((8, NY, NX))
    # Stable: theta falls from 340 K at the lid (index 0) to 300 K at the ground.
    m.theta = (300.0 + np.linspace(40.0, 0.0, 8)[:, None, None]
               + 0.1 * rng.standard_normal((8, NY, NX)))
    m.pi = np.full((NY, NX), 80_000.0)
    m.set_reference()
    a = m.tendencies(m.u, m.v, m.theta, m.pi)
    m.div_damp = NU
    b = m.tendencies(m.u, m.v, m.theta, m.pi)
    ddu, ddv = divergence_damping(m.u, m.v, g, NU)
    assert np.array_equal(b[0] - a[0], ddu) or np.allclose(b[0] - a[0], ddu, rtol=0, atol=1e-18)
    assert np.allclose(b[1] - a[1], ddv, rtol=0, atol=1e-18)
    assert np.array_equal(a[2], b[2]) and np.array_equal(a[3], b[3])
    return "core: off by default; when on, adds exactly the damping to u and v only"


if __name__ == "__main__":
    tests = [test_rotational_flow_is_untouched, test_2dx_wave_decays_at_4nu_over_dx2,
             test_torch_matches_numpy, test_off_by_default_and_bit_identical]
    ok = 0
    for t in tests:
        try:
            print(f"  PASS  {t.__name__}: {t()}")
            ok += 1
        except Exception as e:                            # noqa: BLE001
            print(f"  FAIL  {t.__name__}: {type(e).__name__}: {e}")
    print(f"{ok}/{len(tests)} passed")
    sys.exit(0 if ok == len(tests) else 1)

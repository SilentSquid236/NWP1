"""
The torch backend must reproduce the NumPy core.

Both run the same realistic case (backend_reference.build: 800 m ridge,
baroclinic jet, noise, mixing, drag and convective adjustment on) for 150
steps. Float64 arithmetic in a different order differs only by round-off,
so the states must agree to 1e-9 relative. If PyTorch is not installed the
test says so and passes (the NumPy path is then the only one in use).

Run:  python src/dynamics/test_backend.py
"""
import sys
import time
from pathlib import Path

import numpy as np

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE))

import backend                                    # noqa: E402
import backend_reference as br                    # noqa: E402

results = []


def report(name, ok, detail):
    results.append(bool(ok))
    print(f"  [{'PASS' if ok else 'FAIL'}] {name}\n        {detail}")


def run(backend_name, steps=150, threads=None, nx=44, ny=40):
    m = br.build(nx, ny)
    if backend_name != "numpy":
        m.to_backend(backend_name, threads=threads)
    t0 = time.time()
    for _ in range(steps):
        m.step(20.0)
    return {k: np.asarray(backend.to_numpy(getattr(m, k))) for k in ("u", "v", "theta", "pi")}, time.time() - t0


def test_torch_matches_numpy():
    if backend.TORCH is None:
        report("torch backend matches numpy (skipped: PyTorch not installed)", True, "numpy only")
        return
    a, ta = run("numpy")
    b, tb = run("torch", threads=4)
    rel = {k: float(np.abs(a[k] - b[k]).max() / np.abs(a[k]).max()) for k in a}
    report("torch backend reproduces numpy to round-off over 150 steps",
           max(rel.values()) < 1e-9 and all(np.isfinite(b[k]).all() for k in b),
           "max relative difference " + ", ".join(f"{k} {v:.1e}" for k, v in rel.items())
           + f"; numpy {ta:.1f} s, torch {tb:.1f} s")


def test_numpy_path_untouched_by_torch_state():
    """A NumPy model after a torch model has run: the cache must not leak tensors in."""
    if backend.TORCH is None:
        report("numpy path unaffected by the torch cache (skipped)", True, "numpy only")
        return
    run("torch", steps=5, threads=2)
    a, _ = run("numpy", steps=20)
    b, _ = run("numpy", steps=20)
    report("numpy path unaffected by a torch run in the same process",
           all(isinstance(a[k], np.ndarray) and np.array_equal(a[k], b[k]) for k in a),
           "two numpy runs after a torch run are identical")


def test_round_trip():
    if backend.TORCH is None:
        report("state round-trip numpy -> torch -> numpy (skipped)", True, "numpy only")
        return
    m = br.build(12, 10)
    u0 = m.u.copy()
    m.to_backend("torch", threads=1)
    ok_t = m.backend == "torch"
    m.to_backend("numpy")
    report("state round-trip numpy -> torch -> numpy is exact",
           ok_t and m.backend == "numpy" and np.array_equal(m.u, u0), f"backend after: {m.backend}")


if __name__ == "__main__":
    print("=" * 62)
    print("Backends")
    print("=" * 62)
    for fn in (test_torch_matches_numpy, test_numpy_path_untouched_by_torch_state, test_round_trip):
        try:
            fn()
        except Exception as e:
            report(fn.__name__, False, f"raised {type(e).__name__}: {e}")
    print("=" * 62)
    n = sum(results)
    print(f"{n}/{len(results)} passed\n")
    raise SystemExit(0 if n == len(results) else 1)

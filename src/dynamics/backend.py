"""
Array backends for the dynamical core: NumPy (default; Harris et al. 2020) or
PyTorch (CPU; Paszke et al. 2019).

WHY. The core is element-wise array code. NumPy runs it on one core, so a
24 h forecast took ~40 min whatever the thread settings said. PyTorch runs
the same operations multi-threaded: measured on the Xeon at 8 threads,
6-14x on model-sized arrays (tools/bench_threads.py, prompt 115).

HOW. The physics is written ONCE. Every hot function asks `xp_of(array)`
for a namespace with NumPy's names, and gets NumPy itself for NumPy input
(so the NumPy path is unchanged, bit for bit) or `TORCH` for tensors.
`TORCH` maps each name onto the torch equivalent, in float64, on the CPU.

Constants held as NumPy arrays (sigma levels, Coriolis, terrain, the sponge,
relaxation weights) must not meet a tensor directly: `ndarray * tensor`
silently converts the tensor back to NumPy. Hot code therefore passes them
through `xp.asarray`, which is free for NumPy and cached for torch.
"""

import contextlib

import numpy as np

try:
    import torch
except ImportError:                     # the NumPy path needs nothing else
    torch = None


class _TorchNamespace:
    """The NumPy names the core uses, implemented with torch (float64, CPU)."""

    pi = np.pi
    inf = np.inf
    name = "torch"
    CACHE_MAX = 256

    def __init__(self):
        self._cache = {}

    # --- conversion -------------------------------------------------------
    def asarray(self, a, dtype=None):
        if isinstance(a, torch.Tensor):
            return a if dtype is None else a.to(self._dt(dtype))
        if isinstance(a, np.ndarray):
            key = id(a)
            hit = self._cache.get(key)
            if hit is not None and hit[0] is a:
                return hit[1]
            t = torch.from_numpy(np.ascontiguousarray(a))
            if t.dtype == torch.float32:
                t = t.to(torch.float64)
            # Cache the long-lived constants (grid, levels, weights, a frozen
            # driving frame). Holding `a` keeps its id unique. Bounded, so an
            # array made fresh every step (an interpolated driving state)
            # cannot grow the cache without limit.
            if len(self._cache) >= self.CACHE_MAX:
                self._cache.pop(next(iter(self._cache)))
            self._cache[key] = (a, t)
            return t
        return torch.as_tensor(a, dtype=self._dt(dtype) if dtype else torch.float64)

    @staticmethod
    def _dt(dtype):
        if dtype in (bool, np.bool_, "bool"):
            return torch.bool
        if dtype in (int, np.int64, "int64"):
            return torch.int64
        return torch.float64

    def _t(self, x):
        return x if isinstance(x, torch.Tensor) else self.asarray(x)

    # --- construction -----------------------------------------------------
    def zeros(self, shape, dtype=float):
        return torch.zeros(shape, dtype=self._dt(dtype))

    def ones(self, shape, dtype=float):
        return torch.ones(shape, dtype=self._dt(dtype))

    def empty(self, shape, dtype=float):
        return torch.empty(shape, dtype=self._dt(dtype))

    def full(self, shape, value, dtype=float):
        return torch.full(shape, value, dtype=self._dt(dtype))

    def zeros_like(self, a, dtype=None):
        return torch.zeros_like(a) if dtype is None else torch.zeros_like(a, dtype=self._dt(dtype))

    def ones_like(self, a):
        return torch.ones_like(a)

    def empty_like(self, a):
        return torch.empty_like(a)

    def copy(self, a):
        return a.clone()

    def astype(self, a, dtype):
        return a.to(self._dt(dtype))

    # --- element-wise -----------------------------------------------------
    def abs(self, a):
        return torch.abs(a)

    def sqrt(self, a):
        return torch.sqrt(a)

    def log(self, a):
        return torch.log(a)

    def exp(self, a):
        return torch.exp(a)

    def cos(self, a):
        return torch.cos(a)

    def sin(self, a):
        return torch.sin(a)

    def isfinite(self, a):
        return torch.isfinite(a)

    def where(self, c, a, b):
        like = a if isinstance(a, torch.Tensor) else b
        if not isinstance(a, torch.Tensor):
            a = torch.as_tensor(a, dtype=like.dtype if isinstance(like, torch.Tensor) else torch.float64)
        if not isinstance(b, torch.Tensor):
            b = torch.as_tensor(b, dtype=a.dtype)
        return torch.where(c, a, b)

    def maximum(self, a, b):
        if not isinstance(b, torch.Tensor):
            return torch.clamp(a, min=b)
        if not isinstance(a, torch.Tensor):
            return torch.clamp(b, min=a)
        return torch.maximum(a, b)

    def minimum(self, a, b):
        if not isinstance(b, torch.Tensor):
            return torch.clamp(a, max=b)
        if not isinstance(a, torch.Tensor):
            return torch.clamp(b, max=a)
        return torch.minimum(a, b)

    def clip(self, a, lo, hi):
        return torch.clamp(a, lo, hi)

    # --- reductions and shape ---------------------------------------------
    def sum(self, a, axis=None):
        return a.sum() if axis is None else a.sum(dim=axis)

    def mean(self, a, axis=None, keepdims=False):
        if a.dtype == torch.bool:
            a = a.to(torch.float64)
        if axis is None:
            return a.mean()
        return a.mean(dim=axis, keepdim=keepdims)

    def max(self, a):
        return a.max()

    def min(self, a):
        return a.min()

    def any(self, a):
        return a.any()

    def cumsum(self, a, axis=0):
        return torch.cumsum(a, dim=axis)

    def concatenate(self, arrays, axis=0):
        return torch.cat([self._t(x) for x in arrays], dim=axis)

    def roll(self, a, shift, axis):
        return torch.roll(a, shifts=shift, dims=axis)

    def expand_dims(self, a, axis):
        return a.unsqueeze(axis)

    def broadcast_to(self, a, shape):
        return a.expand(shape)

    def ndim(self, a):
        return a.dim()

    def errstate(self, **kw):
        return contextlib.nullcontext()


TORCH = _TorchNamespace() if torch is not None else None


class _NumpyNamespace:
    """NumPy with the few extra names the core uses (asarray is np.asarray)."""
    name = "numpy"

    def __getattr__(self, k):
        return getattr(np, k)

    @staticmethod
    def copy(a):
        return a.copy()

    @staticmethod
    def astype(a, dtype):
        return a.astype(dtype)

    @staticmethod
    def max(a):
        return a.max()

    @staticmethod
    def min(a):
        return a.min()


NUMPY = _NumpyNamespace()


def xp_of(*arrays):
    """The namespace for these arrays: TORCH if any is a tensor, else NUMPY."""
    if torch is not None:
        for a in arrays:
            if isinstance(a, torch.Tensor):
                return TORCH
    return NUMPY


def to_numpy(a):
    """A NumPy copy of an array or tensor (scalars pass through)."""
    if torch is not None and isinstance(a, torch.Tensor):
        return a.detach().cpu().numpy().copy()
    return a


def set_threads(n):
    """Torch intra-op threads. The measured sweet spot on the Xeon is ~8."""
    if torch is None:
        raise RuntimeError("PyTorch is not installed; use the numpy backend")
    torch.set_num_threads(int(n))
    return torch.get_num_threads()

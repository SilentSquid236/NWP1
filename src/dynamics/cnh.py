"""
Build and load the compiled kernels of the non-hydrostatic core (CAM stage S4a).

nh3d_kernels.c is compiled with the machine's own gcc the first time it is
needed and cached by the hash of its source and flags; nothing is installed.
Python calls the kernels through ctypes. The NumPy code in nh3d.py remains the
reference; test_nh3d_c.py checks every kernel against it.

    import cnh
    lib = cnh.load()          # raises cnh.Unavailable if there is no compiler
    lib.set_threads(26)

Where the library is cached: $NWP_CBUILD if set, otherwise
<tempdir>/nwp1_cbuild. Delete the directory to force a rebuild.
"""
import ctypes
import hashlib
import os
import subprocess
import sys
import tempfile
from pathlib import Path

import numpy as np

SRC = Path(__file__).with_name("nh3d_kernels.c")
FLAGS = ["-O3", "-march=native", "-fopenmp", "-ffp-contract=off", "-std=gnu99"]
# -ffp-contract=off: no fused multiply-add, so a kernel rounds exactly as the
# NumPy expression it transcribes and the tests can compare to round-off.

_LIB = None


class Unavailable(RuntimeError):
    """No working C compiler (or the build failed); the NumPy core still works."""


def _build(real):
    suffix = ".dll" if sys.platform == "win32" else ".so"
    extra = [] if real == "f64" else ["-DREAL=float", "-DSFX=_f32"]
    flags = FLAGS + ([] if sys.platform == "win32" else ["-fPIC"]) + extra
    key = hashlib.sha1(SRC.read_bytes() + " ".join(flags).encode()).hexdigest()[:12]
    cache = Path(os.environ.get("NWP_CBUILD", Path(tempfile.gettempdir()) / "nwp1_cbuild"))
    cache.mkdir(parents=True, exist_ok=True)
    out = cache / f"nh3d_kernels_{real}_{key}{suffix}"
    if not out.exists():
        tmp = out.with_name(out.name + f".{os.getpid()}.tmp")
        cmd = [os.environ.get("NWP_CC", "gcc"), *flags, "-shared", str(SRC), "-o", str(tmp)]
        try:
            r = subprocess.run(cmd, capture_output=True, text=True)
        except OSError as e:          # not found, or cannot be started on this host
            raise Unavailable(f"no usable C compiler ({cmd[0]}): {e}") from None
        if r.returncode:
            raise Unavailable(f"build failed:\n{' '.join(cmd)}\n{r.stderr[-2000:]}")
        os.replace(tmp, out)          # atomic: concurrent runs never load a half file
    return ctypes.CDLL(str(out))


_I = ctypes.POINTER(ctypes.c_int)


def _ptr(a, ctype):
    return a.ctypes.data_as(ctypes.POINTER(ctype))


class Kernels:
    """The kernel set for one precision (f64 or f32)."""

    def __init__(self, real="f64"):
        self.real = real
        self.dtype = np.float64 if real == "f64" else np.float32
        self.ct = ctypes.c_double if real == "f64" else ctypes.c_float
        self.lib = _build(real)
        s = "_" + real
        self._threads = getattr(self.lib, "nwp_set_threads" + s)
        self._threads.restype = ctypes.c_int
        self._pd = getattr(self.lib, "ac_pd" + s)
        self._uv = getattr(self.lib, "ac_uv" + s)
        self._col = getattr(self.lib, "ac_col" + s)
        self._td_col = getattr(self.lib, "td_col" + s)
        self._td_uv = getattr(self.lib, "td_uv" + s)
        self._td_tw = getattr(self.lib, "td_tw" + s)
        self._as_a = getattr(self.lib, "as_a" + s)
        self._as_b = getattr(self.lib, "as_b" + s)
        self._add = getattr(self.lib, "add_into" + s)

    def set_threads(self, n):
        """Set the OpenMP thread count (n <= 0 leaves it); returns the count in use."""
        return int(self._threads(ctypes.c_int(int(n))))

    # --- argument checking: a wrong dtype or layout must fail loudly --------
    def _a(self, a):
        if not (isinstance(a, np.ndarray) and a.dtype == self.dtype and a.flags.c_contiguous):
            raise TypeError(f"kernel argument must be C-contiguous {np.dtype(self.dtype).name}, "
                            f"got {getattr(a, 'dtype', type(a))} "
                            f"contiguous={getattr(getattr(a, 'flags', None), 'c_contiguous', None)}")
        return _ptr(a, self.ct)

    def _r(self, x):
        return self.ct(float(x))

    @staticmethod
    def _idx(a):
        if not (a.dtype == np.int32 and a.flags.c_contiguous):
            raise TypeError("index tables must be C-contiguous int32")
        return a.ctypes.data_as(_I)

    # --- kernels -----------------------------------------------------------
    def ac_pd(self, Qc, T2, E, f2, p_old, damp, pd):
        nz, ny, nx = T2.shape
        self._pd(nz, ny, nx, self._a(Qc), self._a(T2), self._a(E), self._a(f2),
                 self._a(p_old), self._r(damp), self._a(pd))

    def ac_uv(self, xm1, ym1, dtau, dx, dy, FU, FV, cu, mu_u, cv, mu_v, pd, f2, U2, V2):
        nz, ny, nx = U2.shape
        self._uv(nz, ny, nx, self._idx(xm1), self._idx(ym1), self._r(dtau), self._r(dx),
                 self._r(dy), self._a(FU), self._a(FV), self._a(cu), self._a(mu_u),
                 self._a(cv), self._a(mu_v), self._a(pd), self._a(f2), self._a(U2), self._a(V2))

    def ac_col(self, xp1, yp1, dtau, dx, dy, G, am, ds, dsw, Gk, dmuF, FTh, Fphi, FW,
               thw, thu, thv, Qc, E, dphidx, dphidy, dphi_ds, mus, Bv, lower, diag, upper,
               pd, U2, V2, m2, T2, W2, f2):
        nz, ny, nx = T2.shape
        a = self._a
        self._col(nz, ny, nx, self._idx(xp1), self._idx(yp1), self._r(dtau), self._r(dx),
                  self._r(dy), self._r(G), self._r(am), a(ds), a(dsw), a(Gk),
                  a(dmuF), a(FTh), a(Fphi), a(FW), a(thw), a(thu), a(thv), a(Qc), a(E),
                  a(dphidx), a(dphidy), a(dphi_ds), a(mus), a(Bv), a(lower), a(diag), a(upper),
                  a(pd), a(U2), a(V2), a(m2), a(T2), a(W2), a(f2))

    # --- tendencies (nh3d.tendencies) ----------------------------------------
    def td_col(self, idx, dx, dy, P0, RD, GAMMA, ds, sf, pb, alb, mub, phib,
               mu, U, V, W, Th, phi, out):
        nz, ny, nx = U.shape
        a, i = self._a, self._idx
        o = out
        self._td_col(nz, ny, nx, i(idx["xm1"]), i(idx["xp1"]), i(idx["ym1"]), i(idx["yp1"]),
                     self._r(dx), self._r(dy), self._r(P0), self._r(RD), self._r(GAMMA),
                     a(ds), a(sf), a(pb), a(alb), a(mub), a(phib),
                     a(mu), a(U), a(V), a(W), a(Th), a(phi),
                     a(o["pp"]), a(o["dpp"]), a(o["phim"]), a(o["phipm"]), a(o["A1"]), a(o["A2"]),
                     a(o["u"]), a(o["v"]), a(o["th"]), a(o["sdm"]),
                     a(o["Om"]), a(o["sd"]), a(o["w"]), a(o["dmu"]), a(o["mup"]))

    def td_uv(self, idx, dx, dy, sf, f_u, f_v, mu, mub, pb, w_, U, V, FU, FV):
        nz, ny, nx = U.shape
        a, i = self._a, self._idx
        self._td_uv(nz, ny, nx, i(idx["xm1"]), i(idx["xp1"]), i(idx["xm2"]), i(idx["xp2"]),
                    i(idx["ym1"]), i(idx["yp1"]), i(idx["ym2"]), i(idx["yp2"]),
                    self._r(dx), self._r(dy), a(sf), a(f_u), a(f_v), a(mu), a(mub), a(pb),
                    a(w_["pp"]), a(w_["dpp"]), a(w_["phim"]), a(w_["phipm"]), a(w_["A1"]),
                    a(w_["A2"]), a(U), a(V), a(w_["u"]), a(w_["v"]), a(w_["sdm"]), a(FU), a(FV))

    def td_tw(self, idx, dx, dy, G, ds, dsw, sh, mu, U, V, w_, phi, FTh, FW, Fphi):
        nz, ny, nx = U.shape
        a, i = self._a, self._idx
        self._td_tw(nz, ny, nx, i(idx["xm1"]), i(idx["xp1"]), i(idx["xm2"]), i(idx["xp2"]),
                    i(idx["ym1"]), i(idx["yp1"]), i(idx["ym2"]), i(idx["yp2"]),
                    self._r(dx), self._r(dy), self._r(G), a(ds), a(dsw), a(sh),
                    a(mu), a(U), a(V), a(w_["th"]), a(w_["Om"]), a(w_["pp"]), a(w_["mup"]),
                    a(w_["u"]), a(w_["v"]), a(w_["w"]), a(w_["sd"]), a(phi),
                    a(FTh), a(FW), a(Fphi))

    # --- acoustic set-up ------------------------------------------------------
    def as_a(self, P0, RD, GAMMA, ds, mus, Ths, phis, al, Qc, E, ths):
        nz, ny, nx = Ths.shape
        a = self._a
        self._as_a(nz, ny, nx, self._r(P0), self._r(RD), self._r(GAMMA), a(ds), a(mus),
                   a(Ths), a(phis), a(al), a(Qc), a(E), a(ths))

    def as_b(self, idx, dx, dy, bvc, sh, Gk, mus, al, ths, phis, Xs, X0, Qc, E, o):
        nz, ny, nx = ths.shape
        a, i = self._a, self._idx
        self._as_b(nz, ny, nx, i(idx["xm1"]), i(idx["xp1"]), i(idx["ym1"]), i(idx["yp1"]),
                   self._r(dx), self._r(dy), self._r(bvc), a(sh), a(Gk),
                   a(mus), a(al), a(ths), a(phis), a(Xs[1]), a(Xs[2]), a(Xs[3]), a(Xs[4]),
                   a(X0[0]), a(X0[1]), a(X0[2]), a(X0[3]), a(X0[4]), a(X0[5]), a(Qc), a(E),
                   a(o["cu"]), a(o["mu_u"]), a(o["cv"]), a(o["mu_v"]), a(o["thu"]), a(o["thv"]),
                   a(o["thw"]), a(o["dphidx"]), a(o["dphidy"]), a(o["dphi_ds"]), a(o["Bv"]),
                   a(o["lower"]), a(o["diag"]), a(o["upper"]),
                   a(o["m2"]), a(o["U2"]), a(o["V2"]), a(o["W2"]), a(o["T2"]), a(o["f2"]),
                   a(o["p_old"]))

    def add_into(self, x, y, out):
        if not (x.shape == y.shape == out.shape):
            raise ValueError("add_into: shapes differ")
        self._add(ctypes.c_size_t(x.size), self._a(x), self._a(y), self._a(out))


_CACHE = {}


def load(real="f64"):
    """The kernel set for 'f64' or 'f32' (built on first use)."""
    if real not in _CACHE:
        _CACHE[real] = Kernels(real)
    return _CACHE[real]


def available():
    try:
        load("f64")
        return True
    except Unavailable:
        return False


def shift_index(n_points, n, edge_mode):
    """int32 table t with a[t] == grid.shift(a, n, axis) along one axis."""
    i = np.arange(n_points) + n
    if edge_mode == "periodic":
        i = i % n_points
    else:
        i = np.clip(i, 0, n_points - 1)
    return np.ascontiguousarray(i, dtype=np.int32)

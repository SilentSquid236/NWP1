"""Bind a process's memory to the NUMA node of its CPUs (Linux only).

Why (RESEARCH_LOG 2026-10-11, operations): on the two-socket server the page
cache fills each node's memory, and with vm.zone_reclaim_mode = 0 the kernel
places new allocations on the *other* node rather than dropping local cache.
A 3 km run pinned to node 1 had 12 of its 19 GB on node 0 and ran 1.4-1.7x
slower than the same run on node 0 (IPC 0.50 vs 0.58).

With MPOL_BIND to the node(s) the process is pinned to, the kernel reclaims
local page cache instead. Nothing numerical changes. The policy applies to this
thread and every thread created after it, so call bind_to_cpu_nodes() before
OpenMP or torch start their pools.

Opt out with NWP_NUMA_BIND=0. It is a no-op when the affinity spans every
node, off Linux, or when the syscall fails.

WHOSE affinity (corrected 2026-10-12): the launching shell's, not this
process's. With OMP_PROC_BIND=close, the OpenMP runtime that torch loads at
import pins the Python main thread to its first place (CPU 0) before main()
runs, so this process's own affinity said "node 0" even for an unpinned cron
run whose 26 threads span both sockets. Binding all memory to node 0 then made
that run 22 % slower (A/B, 2026-10-12). `taskset -c MASK bash cam_cycle.sh`
puts the mask on the shell, which is the parent of every forecast.py it starts.
"""
import ctypes
import os

_SYS_SET_MEMPOLICY = {"x86_64": 238, "aarch64": 237}
_MPOL_BIND = 2


def cpu_nodes(cpus=None, sysfs="/sys/devices/system/node"):
    """Map the given CPUs (default: this process's affinity) to their NUMA nodes."""
    if cpus is None:
        cpus = os.sched_getaffinity(0)
    nodes, all_nodes = set(), set()
    for d in os.listdir(sysfs):
        if not (d.startswith("node") and d[4:].isdigit()):
            continue
        n = int(d[4:])
        all_nodes.add(n)
        try:
            with open(os.path.join(sysfs, d, "cpulist")) as f:
                txt = f.read().strip()
        except OSError:
            continue
        for part in filter(None, txt.split(",")):
            lo, _, hi = part.partition("-")
            rng = range(int(lo), int(hi or lo) + 1)
            if any(c in cpus for c in rng):
                nodes.add(n)
    return nodes, all_nodes


def launcher_cpus():
    """CPUs the launching shell may use (the intended placement), else our own."""
    try:
        return os.sched_getaffinity(os.getppid())
    except (OSError, AttributeError):
        return os.sched_getaffinity(0)


def bind_to_cpu_nodes(verbose=True, cpus=None):
    """Return a one-line description of what was done (for the run log)."""
    if os.environ.get("NWP_NUMA_BIND", "1") == "0":
        return "numa bind     : off (NWP_NUMA_BIND=0)"
    if not hasattr(os, "sched_getaffinity") or not os.path.isdir("/sys/devices/system/node"):
        return "numa bind     : not applicable on this system"
    try:
        nodes, all_nodes = cpu_nodes(launcher_cpus() if cpus is None else cpus)
    except OSError as e:
        return f"numa bind     : skipped ({e})"
    if len(all_nodes) < 2 or not nodes or nodes == all_nodes:
        return f"numa bind     : not needed (CPUs on nodes {sorted(nodes)} of {sorted(all_nodes)})"
    nr = _SYS_SET_MEMPOLICY.get(os.uname().machine)
    if nr is None:
        return "numa bind     : skipped (unknown architecture)"
    mask = 0
    for n in nodes:
        mask |= 1 << n
    buf = (ctypes.c_ulong * 2)(mask, 0)
    libc = ctypes.CDLL(None, use_errno=True)
    rc = libc.syscall(nr, _MPOL_BIND, buf, ctypes.c_ulong(128))
    if rc != 0:
        return f"numa bind     : set_mempolicy failed (errno {ctypes.get_errno()})"
    return f"numa bind     : memory bound to node(s) {sorted(nodes)} (launcher pinned there)"

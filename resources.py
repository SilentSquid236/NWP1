"""
Shared-server resource governor, load-aware.

This box is shared with ~30 researchers. Policy:

  * Never exceed a hard ceiling (default 50% of cores), even on an idle box.
  * Below that ceiling, scale to what other users are actually leaving free,
    re-checked between epochs so a long run yields as the box fills up.
  * Never drop below a floor, so a run always makes progress.

IMPORTANT: import this module and call apply() BEFORE torch. OpenMP and MKL
read their thread limits at import time; setting them afterwards is ignored.

    import resources
    PLAN = resources.apply()   # first
    import torch               # second

The env cap is set to the CEILING, and the governor then moves torch's
thread count up and down within that range at runtime via set_num_threads().

Environment overrides:
    NWP_RESOURCE_FRACTION=0.75   raise the ceiling (default 0.5)
    NWP_ADAPTIVE=0               pin to the ceiling, ignore other users
    NWP_MIN_CORES=8              floor (default 4)
"""

import os
import time

DEFAULT_FRACTION = 0.5
DEFAULT_MIN_CORES = 4
RETUNE_THRESHOLD = 4      # ignore adjustments smaller than this many threads
LOADER_WORKERS = 2


def cgroup_cpu_limit():
    """
    CPU quota imposed by cgroups, in cores, or None if unlimited.

    Affinity alone does NOT reveal this: an admin can leave every core
    visible while capping actual CPU time. Checks cgroup v2 then v1.
    """
    # cgroup v2: "<quota> <period>", or "max <period>" when unlimited
    try:
        with open("/sys/fs/cgroup/cpu.max") as f:
            quota, period = f.read().split()
            if quota != "max":
                return int(quota) / int(period)
    except (OSError, ValueError):
        pass

    # cgroup v1
    try:
        with open("/sys/fs/cgroup/cpu/cpu.cfs_quota_us") as f:
            quota = int(f.read().strip())
        with open("/sys/fs/cgroup/cpu/cpu.cfs_period_us") as f:
            period = int(f.read().strip())
        if quota > 0 and period > 0:
            return quota / period
    except (OSError, ValueError):
        pass

    return None


def cgroup_memory_limit_gb():
    """Memory ceiling from cgroups, in GB, or None if unlimited."""
    for path in ("/sys/fs/cgroup/memory.max",
                 "/sys/fs/cgroup/memory/memory.limit_in_bytes"):
        try:
            with open(path) as f:
                v = f.read().strip()
            if v == "max":
                return None
            v = int(v)
            # v1 reports an absurd sentinel when unlimited
            if v < (1 << 62):
                return v / 1024**3
        except (OSError, ValueError):
            continue
    return None


def total_cores() -> int:
    """
    Cores actually usable: the smaller of CPU affinity and any cgroup quota.
    """
    try:
        affinity = len(os.sched_getaffinity(0))
    except AttributeError:
        affinity = os.cpu_count() or 1

    quota = cgroup_cpu_limit()
    if quota is not None:
        return max(1, min(affinity, int(quota)))
    return affinity


def visible_cores() -> int:
    """Cores on the machine, ignoring any quota -- used for load context."""
    return os.cpu_count() or 1


def total_memory_gb() -> float:
    cg = cgroup_memory_limit_gb()
    if cg is not None:
        return cg
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return float("nan")


def load_average() -> float:
    """1-minute load average. Includes our own threads."""
    try:
        return os.getloadavg()[0]
    except (OSError, AttributeError):
        return 0.0


def resolve_fraction(explicit=None) -> float:
    frac = float(explicit) if explicit is not None else \
        float(os.environ.get("NWP_RESOURCE_FRACTION", DEFAULT_FRACTION))
    if not 0 < frac <= 1.0:
        raise ValueError(f"resource fraction must be in (0, 1], got {frac}")
    return frac


def adaptive_enabled() -> bool:
    return os.environ.get("NWP_ADAPTIVE", "1") not in ("0", "false", "False")


def min_cores() -> int:
    return max(1, int(os.environ.get("NWP_MIN_CORES", DEFAULT_MIN_CORES)))


def plan(fraction=None) -> dict:
    """Static plan: the ceiling, and what we start at."""
    frac = resolve_fraction(fraction)
    cores = total_cores()
    ceiling = max(1, int(cores * frac))
    workers = min(LOADER_WORKERS, max(0, ceiling - 1))

    return {
        "fraction": frac,
        "total_cores": cores,
        "ceiling_cores": ceiling,
        "torch_threads": max(1, ceiling - workers),
        "dataloader_workers": workers,
        "adaptive": adaptive_enabled(),
        "min_cores": min_cores(),
        "total_memory_gb": total_memory_gb(),
        "memory_budget_gb": total_memory_gb() * frac,
        "cpu_quota": cgroup_cpu_limit(),
        "mem_quota_gb": cgroup_memory_limit_gb(),
        "visible_cores": visible_cores(),
    }


def apply(fraction=None) -> dict:
    """Set env thread caps to the CEILING. Must run before torch import."""
    p = plan(fraction)
    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = str(p["torch_threads"])
    return p


class LoadGovernor:
    """
    Re-tunes thread count from observed load on the rest of the machine.

    others_load = load1 - (our own running threads)
    free        = cores - others_load
    target      = clamp(free * fraction, floor, ceiling)
    """

    def __init__(self, plan_dict):
        self.p = plan_dict
        self.current = plan_dict["torch_threads"]
        self.ceiling = plan_dict["torch_threads"]
        self.floor = min(plan_dict["min_cores"], self.ceiling)
        self.history = []

    def observe(self) -> dict:
        cores = self.p["total_cores"]
        load1 = load_average()
        others = max(0.0, load1 - self.current)
        free = max(0.0, cores - others)
        target = int(free * self.p["fraction"])
        target = max(self.floor, min(self.ceiling, target))
        return {"load1": load1, "others": others, "free": free, "target": target}

    def update(self, torch_module=None) -> dict:
        """Sample load and adjust threads if the change is worth making."""
        obs = self.observe()
        obs["previous"] = self.current
        obs["changed"] = False

        if not self.p["adaptive"]:
            obs["target"] = self.current
            return obs

        if abs(obs["target"] - self.current) >= RETUNE_THRESHOLD:
            self.current = obs["target"]
            obs["changed"] = True
            if torch_module is not None:
                torch_module.set_num_threads(self.current)

        obs["current"] = self.current
        self.history.append((time.time(), obs["load1"], self.current))
        return obs

    def format(self, obs) -> str:
        arrow = f"{obs['previous']} -> {obs['current']}" if obs.get("changed") \
            else f"{self.current}"
        return (f"load {obs['load1']:.1f} (others ~{obs['others']:.0f} cores) "
                f"| threads {arrow}")


def describe(p=None) -> str:
    p = p or plan()
    mode = "adaptive" if p["adaptive"] else "fixed"
    cpu_q = "none" if p["cpu_quota"] is None else f"{p['cpu_quota']:.1f} cores"
    mem_q = "none" if p["mem_quota_gb"] is None else f"{p['mem_quota_gb']:.0f} GB"
    return (
        f"  resource cap   : {p['fraction']:.0%} ceiling of {p['total_cores']} cores "
        f"-> {p['ceiling_cores']} core budget ({mode})\n"
        f"  torch threads  : {p['torch_threads']} max "
        f"(+{p['dataloader_workers']} loader workers), "
        f"floor {min(p['min_cores'], p['torch_threads'])}\n"
        f"  memory         : {p['total_memory_gb']:.0f} GB total, "
        f"stay under ~{p['memory_budget_gb']:.0f} GB\n"
        f"  quotas         : cpu {cpu_q}, memory {mem_q}"
    )


if __name__ == "__main__":
    p = apply()
    print("Resource plan:")
    print(describe(p))
    g = LoadGovernor(p)
    obs = g.observe()
    print(f"\nCurrent conditions:")
    print(f"  1-min load     : {obs['load1']:.2f}")
    print(f"  others using   : ~{obs['others']:.0f} cores")
    print(f"  free           : ~{obs['free']:.0f} cores")
    print(f"  would run with : {obs['target']} threads")

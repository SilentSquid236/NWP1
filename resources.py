"""
Shared-server resource governor.

This box is shared with ~30 researchers, so the project takes HALF of the
visible resources by default and never more unless explicitly told to.

IMPORTANT: import this module BEFORE torch. The OpenMP/MKL thread limits are
read by those libraries at import time; setting them afterwards is ignored.

    import resources; resources.apply()   # first
    import torch                          # second

Override the fraction when you know the box is quiet:

    export NWP_RESOURCE_FRACTION=0.75     # take 75%
    python src/train_autoregressive.py --resource-fraction 0.9
"""

import os

DEFAULT_FRACTION = 0.5
_applied = None


def total_cores() -> int:
    """Cores actually available to this process, honouring cpuset/affinity."""
    try:
        return len(os.sched_getaffinity(0))       # Linux, respects cgroups
    except AttributeError:
        return os.cpu_count() or 1


def total_memory_gb() -> float:
    try:
        with open("/proc/meminfo") as f:
            for line in f:
                if line.startswith("MemTotal:"):
                    return int(line.split()[1]) / 1024 / 1024
    except OSError:
        pass
    return float("nan")


def resolve_fraction(explicit=None) -> float:
    """CLI flag wins, then env var, then the 50% default."""
    if explicit is not None:
        frac = float(explicit)
    else:
        frac = float(os.environ.get("NWP_RESOURCE_FRACTION", DEFAULT_FRACTION))
    if not 0 < frac <= 1.0:
        raise ValueError(f"resource fraction must be in (0, 1], got {frac}")
    return frac


def plan(fraction=None) -> dict:
    """Work out the thread/worker budget without applying it."""
    frac = resolve_fraction(fraction)
    cores = total_cores()

    budget = max(1, int(cores * frac))

    # DataLoader workers are separate processes and count against the budget.
    # Keep them modest so most of the allowance goes to the math threads.
    workers = min(2, max(0, budget - 1))
    threads = max(1, budget - workers)

    return {
        "fraction": frac,
        "total_cores": cores,
        "budget_cores": budget,
        "torch_threads": threads,
        "dataloader_workers": workers,
        "total_memory_gb": total_memory_gb(),
        "memory_budget_gb": total_memory_gb() * frac,
    }


def apply(fraction=None) -> dict:
    """Set the environment thread caps. Must run before torch is imported."""
    global _applied
    p = plan(fraction)

    for var in ("OMP_NUM_THREADS", "MKL_NUM_THREADS", "OPENBLAS_NUM_THREADS",
                "NUMEXPR_NUM_THREADS", "VECLIB_MAXIMUM_THREADS"):
        os.environ[var] = str(p["torch_threads"])

    _applied = p
    return p


def describe(p=None) -> str:
    p = p or _applied or plan()
    return (
        f"  resource cap   : {p['fraction']:.0%} of {p['total_cores']} cores "
        f"-> {p['budget_cores']} core budget\n"
        f"  torch threads  : {p['torch_threads']}  "
        f"(+{p['dataloader_workers']} loader workers)\n"
        f"  memory         : {p['total_memory_gb']:.0f} GB total, "
        f"stay under ~{p['memory_budget_gb']:.0f} GB"
    )


if __name__ == "__main__":
    p = apply()
    print("Resource plan:")
    print(describe(p))

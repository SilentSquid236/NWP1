"""Tests for numa_bind (runs anywhere; the bind itself only on Linux)."""
import os, sys, tempfile
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import numa_bind

res = []
def check(name, ok):
    res.append(ok); print(("  [PASS] " if ok else "  [FAIL] ") + name)

# fake sysfs: node0 = CPUs 0,2,4  node1 = 1,3,5 (written as a list and a range)
d = tempfile.mkdtemp()
for n, txt in [(0, "0,2,4"), (1, "1,3-5")]:
    os.makedirs(os.path.join(d, f"node{n}"))
    open(os.path.join(d, f"node{n}", "cpulist"), "w").write(txt + "\n")
os.makedirs(os.path.join(d, "power"))            # non-node entries are ignored
check("CPUs {1,3} map to node 1 only", numa_bind.cpu_nodes({1, 3}, d) == ({1}, {0, 1}))
check("CPUs {0,1} span both nodes", numa_bind.cpu_nodes({0, 1}, d)[0] == {0, 1})
check("a range in cpulist is expanded (CPU 4 is in '3-5' on node 1 and '4' on node 0)",
      numa_bind.cpu_nodes({4}, d)[0] == {0, 1})
os.environ["NWP_NUMA_BIND"] = "0"
check("NWP_NUMA_BIND=0 turns it off", "off" in numa_bind.bind_to_cpu_nodes())
del os.environ["NWP_NUMA_BIND"]
msg = numa_bind.bind_to_cpu_nodes()
check("bind_to_cpu_nodes never raises and reports: " + msg, msg.startswith("numa bind"))
print(f"{sum(res)}/{len(res)} passed")
sys.exit(0 if all(res) else 1)

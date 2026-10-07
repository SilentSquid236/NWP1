"""Pool verification archives over cycles: per-lead bias/RMSE, A against B."""
import json, sys, glob, os, collections
import numpy as np

def load_tree(root):
    out = {}
    for f in glob.glob(os.path.join(root, "*", "matches.jsonl")):
        out[os.path.basename(os.path.dirname(f))] = [json.loads(l) for l in open(f)]
    return out

def errs(rows, var):
    d = collections.defaultdict(list)
    for m in rows:
        if m["variable"] == var: d[int(round(m["lead_hours"]))].append(m["error"])
    return d

def main(a_spec, b_spec, grid=None):
    la, ra = a_spec.split(":", 1); lb, rb = b_spec.split(":", 1)
    A, B = load_tree(ra), load_tree(rb)
    stamps = sorted(set(A) & set(B))
    print(f"A = {la} ({ra}), B = {lb} ({rb}); {len(stamps)} cycles in both: {stamps[0] if stamps else '-'} .. {stamps[-1] if stamps else '-'}")
    for var in ("TMP", "UGRD", "VGRD"):
        pa, pb = collections.defaultdict(list), collections.defaultdict(list)
        per_cycle = []
        for s in stamps:
            ea, eb = errs(A[s], var), errs(B[s], var)
            for L in ea: pa[L] += ea[L]
            for L in eb: pb[L] += eb[L]
            ca = [np.sqrt(np.mean(np.square(ea[L]))) for L in ea if L in eb]
            cb = [np.sqrt(np.mean(np.square(eb[L]))) for L in eb if L in ea]
            if ca: per_cycle.append((s, float(np.mean(ca)), float(np.mean(cb))))
        print(f"\n{var} pooled over cycles\n  lead      n   bias A   rmse A   bias B   rmse B    B - A")
        diffs = []
        for L in sorted(set(pa) & set(pb)):
            a, b = np.array(pa[L]), np.array(pb[L])
            d = np.sqrt(np.mean(b**2)) - np.sqrt(np.mean(a**2)); diffs.append(d)
            print(f"  {L:4d}  {len(b):6d}  {a.mean():+6.2f}  {np.sqrt(np.mean(a**2)):6.2f}  {b.mean():+6.2f}  {np.sqrt(np.mean(b**2)):6.2f}  {d:+6.2f}")
        if diffs:
            print(f"  mean over leads B - A {np.mean(diffs):+.3f}; B better at {sum(x < 0 for x in diffs)} of {len(diffs)} leads; "
                  f"leads 6-24 mean {np.mean(diffs[5:]):+.3f}")
            better = sum(cb < ca for _, ca, cb in per_cycle)
            print(f"  per cycle (mean over leads): B better in {better} of {len(per_cycle)}; " +
                  " ".join(f"{s[4:]}:{cb-ca:+.2f}" for s, ca, cb in per_cycle))

if __name__ == "__main__":
    main(sys.argv[1], sys.argv[2])

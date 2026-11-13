"""The per-cell species matrix (rck.ecdna.species_matrix) of every bench case at one λ, offline and
read-only (no Gurobi; gurobipy is blocked at import).

    python diag/species_matrix_bench.py [--lambda 0.2] [--aa /n/fs/yubi-iw/runs/bench_nested_aa]
        [--out /n/fs/yubi-iw/runs/bench_nested_v2/species_matrix] [--plot CASE ...]

Reads the AA exports that ``diag/aa_regression.py`` wrote (``AA/<case>/lambda_<λ>/<name>_cycles.txt``
and ``<name>_graph.txt``, ``<name>`` = ``union`` on the bench: rank 1 at λ, low_weight species flagged ``Status=low_weight``), and the
truth as that script does (``run_case.cases()``, ``--truth-sample``). The scope regions come from
the graph, as there. Writes:
- ``OUT/lambda_<λ>/<case>/aa_species_matrix.tsv``, and the heatmaps of the ``--plot`` cases;
- ``OUT/lambda_<λ>/aa_species_matrix.tsv``, all cases;
- ``OUT/lambda_<λ>/checks.txt``: every cell with ``exact_cycle`` = 1 should be (1, 1) blind
  (J_seg >= 0.999, because exact_cycle allows 100 bp at each end); cells that are not are listed.
"""
import sys

sys.modules["gurobipy"] = None  # no licence use possible

import argparse  # noqa: E402
import glob  # noqa: E402
import os  # noqa: E402

import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_case  # noqa: E402  (puts the RCK root on sys.path)
from rck.ecdna import aa_eval  # noqa: E402
from rck.ecdna import evaluate as ev  # noqa: E402
from rck.ecdna.species_matrix import plot_species_matrix  # noqa: E402


def truth_sample(case):
    opts = case.get("options", [])
    return opts[opts.index("--truth-sample") + 1] if "--truth-sample" in opts else None


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--lambda", dest="lam", type=float, default=0.2)
    ap.add_argument("--aa", default="/n/fs/yubi-iw/runs/bench_nested_aa")
    ap.add_argument("--out", default="/n/fs/yubi-iw/runs/bench_nested_v2/species_matrix")
    ap.add_argument("--plot", nargs="*", default=[])
    a = ap.parse_args()
    out = os.path.join(a.out, "lambda_%g" % a.lam)
    frames = []
    for case in run_case.cases():
        folder = os.path.join(a.aa, case["name"], "lambda_%g" % a.lam)
        found = [f for f in glob.glob(os.path.join(folder, "*_cycles.txt")) if "_alt" not in os.path.basename(f)]
        if len(found) != 1:
            print("no single export for", case["name"], found)
            continue
        cycles = found[0]
        graph = cycles[:-len("_cycles.txt")] + "_graph.txt"
        truth = ev.load_truth(case["truth"], None, truth_sample(case))
        result = aa_eval.score(cycles, truth, graph if os.path.exists(graph) else None, case=case["name"])
        m = result.matrix
        d = os.path.join(out, case["name"])
        os.makedirs(d, exist_ok=True)
        m.to_csv(os.path.join(d, "aa_species_matrix.tsv"), sep="\t", index=False)
        if case["name"] in a.plot:
            for hap in m["haplotype"].unique():
                plot_species_matrix(m, os.path.join(d, "aa_species_matrix_%s.png" % hap), hap, "%s, λ = %g" % (case["name"], a.lam))
        frames.append(m)
    allm = pd.concat(frames, ignore_index=True)
    allm.to_csv(os.path.join(out, "aa_species_matrix.tsv"), sep="\t", index=False)
    blind = allm[allm.haplotype == "blind"]
    exact = blind[blind.exact_cycle == 1]
    bad = exact[(exact.j_seg < 0.999) | (exact.j_adj != 1.0)]
    lines = ["%d cases, %d cells (blind), %d with exact_cycle = 1; not (1, 1): %d" % (allm.case.nunique(), len(blind), len(exact), len(bad))]
    lines += ["  %s %s x %s: j_seg %.6f j_adj %s" % (r.case, r.truth, r.inferred, r.j_seg, r.j_adj) for r in bad.itertuples()]
    near = exact[(exact.j_seg < 1.0) & (exact.j_seg >= 0.999)]
    lines.append("exact cells with 0.999 <= j_seg < 1 (breakpoints within 100 bp): %d" % len(near))
    with open(os.path.join(out, "checks.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")
    print("\n".join(lines))


if __name__ == "__main__":
    main()

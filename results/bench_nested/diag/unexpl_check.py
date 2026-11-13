"""Step 3 check, offline: re-rank saved pools with RCK's own --unexplained-junction-penalty code
(rck.ecdna.decompose.rank_pool, which prunes with the full structure cost), no Gurobi.

    python diag/unexpl_check.py RCK_RUNS [--lambda 0.2] [--mu 1] [--suffix _b1] [--compare nested_prior.tsv] [-o out.tsv]

Pools written before the field existed have no ``unexplained`` count, so it is filled here per
entry with rck.ecdna.decompose.unexplained_junctions: each species of the entry (pool
decompositions.tsv: integer weight, novel junctions) against the run's ecDNA flow on novel
junctions (flow_adjacencies.tsv, summed over phasings). Then rank_pool at (λ, μ) and at (λ, 0).

Prints, as diag/nested_prior.py does: mean, failures, fixed / broken / down > 0.005 against μ = 0,
and checks that μ = 0 reproduces the run's ranking (rank_pool without counts). With --compare (a
nested_prior.py -o table at the same λ) it also lists the cases whose rank 1 differs from that
script's ``unexpl μ`` (pruned on steps) and ``unexpl μ np`` (not pruned) picks.
"""
import sys

sys.modules["gurobipy"] = None  # no licence use possible

import argparse  # noqa: E402
import os  # noqa: E402

import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from junction_rules import FAIL_COMPOSITE, Case  # noqa: E402  (also puts the RCK root on sys.path)
from rck.ecdna.decompose import Species, rank_pool, unexplained_junctions  # noqa: E402
from rck.ecdna.phased import ADJ  # noqa: E402
from rck.ecdna.rerank import pool_params  # noqa: E402


def fill_unexplained(case):
    """Entries with ``unexplained`` filled from the saved species and flow."""
    adj = pd.read_csv(os.path.join(case.pool_dir, "flow_adjacencies.tsv"), sep="\t")
    jflow = {str(k): float(v) for k, v in adj[adj.kind == "novel"].groupby("aid")["ecdna"].sum().items()}
    out = []
    for e in case.entries:
        g = case.sp[case.sp["decomposition_rank"] == int(e["index"]) + 1]
        species = [Species(weight=int(r.weight_int),
                           counts={(ADJ, aid, "-"): 1 for aid in (r.junctions.split(",") if isinstance(r.junctions, str)
                                                                    and r.junctions != "." else [])})
                   for r in g.itertuples()]
        out.append(dict(e, unexplained=unexplained_junctions(species, jflow, case.params)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rck_runs")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.2)
    ap.add_argument("--mu", type=float, default=1.0)
    ap.add_argument("--suffix", default="")
    ap.add_argument("--compare")
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    rows, mismatched_default = [], []
    for folder in sorted(os.listdir(a.rck_runs)):
        if a.suffix and not folder.endswith(a.suffix):
            continue
        if not os.path.exists(os.path.join(a.rck_runs, folder, "pool", "pool", "pool.json")):
            continue
        name = folder[:-len(a.suffix)] if a.suffix else folder
        name = name[:-3] if name.endswith("_m4") else name
        case = Case(os.path.join(a.rck_runs, folder))
        off, on = pool_params(case.pool, a.lam, 0.0), pool_params(case.pool, a.lam, a.mu)
        case.params = off
        filled = fill_unexplained(case)
        run = rank_pool(case.entries, off)[0]
        if rank_pool(filled, off)[0] != run:
            mismatched_default.append(name)
        bench = not name.startswith("demo")
        ux = {int(e["index"]): e["unexplained"] for e in filled}
        for rule, i in (("default", run), ("rck mu=%g" % a.mu, rank_pool(filled, on)[0])):
            comp, merged, nrep = case.composite(i), case.merged(i), case.n_reported(i)
            rows.append({"case": name, "bench": bench, "rule": rule, "pool_index": i, "composite": round(comp, 4),
                         "merged": merged, "n_reported": nrep, "n_truth": case.n_truth, "unexplained": ux[i],
                         "fail": bench and (not comp >= FAIL_COMPOSITE or merged or nrep != case.n_truth)})
    t = pd.DataFrame(rows)
    if a.out:
        t.to_csv(a.out, sep="\t", index=False)
    base = t[t.rule == "default"].set_index("case")
    print("%s, λ = %g, μ = %g: %d cases (%d bench)" % (a.rck_runs, a.lam, a.mu, len(base), int(base.bench.sum())))
    print("μ = 0 with the filled counts reproduces the run's ranking: %s\n"
          % ("yes" if not mismatched_default else "NO: " + ", ".join(mismatched_default)))
    for rule, g in t.groupby("rule", sort=False):
        g = g.set_index("case")
        d = g.composite - base.composite
        fixed = [c for c in g.index if base.at[c, "fail"] and not g.at[c, "fail"]]
        broken = [c for c in g.index if g.at[c, "fail"] and not base.at[c, "fail"]]
        down = list(d[d < -0.005].index)
        print("%-12s %8.4f %6d  + %s | - %s | down: %s" % (rule, g.composite.mean(), int(g.fail.sum()), ", ".join(fixed) or "none",
                                                     ", ".join(broken) or "none", ", ".join(down) or "none"))
    if a.compare:
        ref = pd.read_csv(a.compare, sep="\t")
        mine = t[t.rule != "default"].set_index("case").pool_index
        for rule in ("unexpl %g" % a.mu, "unexpl %g np" % a.mu):
            other = ref[ref.rule == rule].set_index("case").pool_index.reindex(mine.index)
            diff = list(mine.index[mine != other])
            print("rank 1 differs from nested_prior.py '%s' in %d case(s)%s" % (rule, len(diff), (": " + ", ".join(diff)) if diff else ""))
        ux = ref[ref.rule == "default"].set_index("case").unexplained
        mine_ux = t[t.rule == "default"].set_index("case").unexplained.reindex(ux.index)
        diff = list(ux.index[ux != mine_ux])
        print("unexplained count of the default rank 1 differs from nested_prior.py in %d case(s)%s"
              % (len(diff), (": " + ", ".join(diff)) if diff else ""))


if __name__ == "__main__":
    main()

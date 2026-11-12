"""Step 4, offline: ranking priors for the nested-species (prior-class) failures, on saved pools
(no Gurobi; gurobipy is blocked at import).

    python diag/nested_prior.py RCK_RUNS [--suffix _b1] [--lambda 0.2] [-o out.tsv]

RCK_RUNS holds <case>/pool/pool/pool.json; with --suffix only folders ending in it are read, and the
suffix is dropped from the case name (e.g. the depth-bias validation runs <case>_b1, <case>_m4_b1).

The prior class: on noise-free flows (sigma 0, or --depth-bias-window) rank 1 is "outer species +
1-3 copy low_weight leftovers", while the pool holds "outer + the nested species" at the same fit.
The leftovers then carry the nested species' own novel junctions, at a flow the size of a species,
and the ranking does not see it because it counts significant species only.

Rules (added to the structure cost of rank_pool's ranking_key; ``np`` = without the λ pool pruning):

* ``default``      the run's ranking (must reproduce the run / collect.py at λ).
* ``unexpl mu``    + mu per *unexplained junction*: a novel junction whose flow (flow_adjacencies.tsv,
                   ecDNA layer) is at least the decomposition's own low_weight threshold
                   max(1, 0.1 x heaviest integer weight) but that no significant species uses.
* ``mass``         leftover mass, as diag/grow_rules.py: floor(sum of low_weight copies / threshold) more
                   significant species.

Scoring as collect.py (composite < 0.95, merged, or a wrong reported species count fails a bench case),
via diag/junction_rules.py's Case.
"""
import sys

sys.modules["gurobipy"] = None  # no licence use possible

import argparse  # noqa: E402
import math  # noqa: E402
import os  # noqa: E402

import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from junction_rules import FAIL_COMPOSITE, Case, pool_keeps, pool_params, ranking_key  # noqa: E402

MUS = [0.25, 0.5, 1.0, 2.0]


def unexplained(case):
    """{pool index: number of unexplained junctions} (see the module docstring)."""
    adj = pd.read_csv(os.path.join(case.pool_dir, "flow_adjacencies.tsv"), sep="\t")
    flow = adj[adj.kind == "novel"].groupby("aid")["ecdna"].sum()
    out = {}
    for e in case.entries:
        i = int(e["index"])
        g = case.sp[case.sp["decomposition_rank"] == i + 1]
        if not len(g):
            out[i] = 0
            continue
        thr = max(1.0, 0.1 * max(e["weights_int"]))
        # the run's integer low_weight flag of each species row, by its integer weight (the flag
        # depends on the weight only, so equal weights share it)
        flags = dict(zip(e["weights_int"], e["low_weight"]))
        used = set()
        for r in g.itertuples():
            if not flags.get(int(r.weight_int), r.status == "low_weight") and isinstance(r.junctions, str) and r.junctions != ".":
                used.update(r.junctions.split(","))
        out[i] = int(sum(1 for aid, f in flow.items() if f >= thr and aid not in used))
    return out


def rank1(case, params, extra=None, prune=True, entries=None):
    entries = entries or case.entries
    best = entries[0]
    kept = [e for e in entries if not prune or len(entries) <= 1
            or pool_keeps(e["objective"], e["steps"], best["objective"], best["steps"], params)]

    def key(e):
        k = ranking_key(e["residual"], e["n_significant"], e["n_species"], e["steps"], e["objective"], params)
        return (k[0], round(k[1] + (extra(int(e["index"])) if extra else 0.0), 6)) + k[2:]

    return int(sorted(kept, key=key)[0]["index"])


def mass_entries(case):
    out = []
    for e in case.entries:
        e2 = dict(e)
        thr = max(1.0, 0.1 * max(e["weights_int"]))
        low = sum(w for w, f in zip(e["weights_int"], e["low_weight"]) if f)
        e2["n_significant"] = e["n_significant"] + int(math.floor(low / thr + 1e-9))
        out.append(e2)
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rck_runs")
    ap.add_argument("--suffix", default="")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.2)
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    rows = []
    for folder in sorted(os.listdir(a.rck_runs)):
        if a.suffix and not folder.endswith(a.suffix):
            continue
        if not os.path.exists(os.path.join(a.rck_runs, folder, "pool", "pool", "pool.json")):
            continue
        name = folder[:-len(a.suffix)] if a.suffix else folder
        name = name[:-3] if name.endswith("_m4") else name
        case = Case(os.path.join(a.rck_runs, folder))
        params = pool_params(case.pool, a.lam)
        ux = unexplained(case)
        bench = not name.startswith("demo")
        picks = {"default": rank1(case, params)}
        for mu in MUS:
            f = lambda i, mu=mu: mu * ux.get(i, 0)  # noqa: E731
            picks["unexpl %g" % mu] = rank1(case, params, f)
            picks["unexpl %g np" % mu] = rank1(case, params, f, prune=False)
        picks["mass"] = rank1(case, params, entries=mass_entries(case))
        for rule, i in picks.items():
            comp, merged, nrep = case.composite(i), case.merged(i), case.n_reported(i)
            rows.append({"case": name, "bench": bench, "rule": rule, "pool_index": i, "composite": round(comp, 4),
                         "merged": merged, "n_reported": nrep, "n_truth": case.n_truth, "unexplained": ux.get(i, 0),
                         "fail": bench and (not comp >= FAIL_COMPOSITE or merged or nrep != case.n_truth)})
    t = pd.DataFrame(rows)
    if a.out:
        t.to_csv(a.out, sep="\t", index=False)
    base = t[t.rule == "default"].set_index("case")
    print("%s%s, λ = %g: %d cases (%d bench)\n" % (a.rck_runs, " (" + a.suffix + ")" if a.suffix else "", a.lam, len(base),
                                                int(base.bench.sum())))
    print("%-16s %8s %6s  %s" % ("rule", "mean", "fails", "fixed (+) / broken (-) / down > 0.005"))
    for rule, g in t.groupby("rule", sort=False):
        g = g.set_index("case")
        d = g.composite - base.composite
        fixed = [c for c in g.index if base.at[c, "fail"] and not g.at[c, "fail"]]
        broken = [c for c in g.index if g.at[c, "fail"] and not base.at[c, "fail"]]
        down = list(d[d < -0.005].index)
        print("%-16s %8.4f %6d  + %s | - %s | down: %s" % (rule, g.composite.mean(), int(g.fail.sum()), ", ".join(fixed) or "none",
                                                         ", ".join(broken) or "none", ", ".join(down) or "none"))
    print("\nunexplained junctions of the default rank 1 (failing cases): " + ", ".join(
        "%s %d" % (c, base.at[c, "unexplained"]) for c in base.index if base.at[c, "fail"]))
    print("gurobipy loaded: %s" % (sys.modules.get("gurobipy") is not None))


if __name__ == "__main__":
    main()

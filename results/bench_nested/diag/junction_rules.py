"""Offline re-ranking of every saved bench pool under two new candidate priors (no Gurobi, no new runs).

    python diag/junction_rules.py RCK_RUNS [--lambda 0.2] [-o rules.tsv]

Idea 1, junction incidence (``jinc mu``). The ranking prior counts significant species only, so a
light nested species split into low_weight pieces is free, and a merged species that carries two
truth junction sets is not charged for the walks that carry them again. Both put one novel
junction into several walks. The rule adds

    mu x sum over novel junctions j of (number of species, low_weight included, whose walk uses j - 1)

to the structure cost (``species_penalty`` x significant + ``step_penalty`` x steps). It counts
species, not traversals, so a true duplication inside one species is not charged.
``jinc mu np`` is the same rule without the λ pool pruning (``pool_keeps``), i.e. every pool
solution may win; that is what a run would rank if the rule also replaced the pruning.

Idea 2, long-read paths. ``paths first``: rank by residual, then the share of multi-junction
long-read support found in order in some species (``reads_explained``, higher first), then as
before. ``paths sum nu``: structure cost + nu x (1 - reads_explained).

Every rule is scored as collect.py scores a bench case: rank-1 mean composite (pool
evaluation.tsv), the merged flag (collect.py::merged_flag, reimplemented on cached tables), and the
reported species count (refitted flags, decompositions.tsv) against the truth count. A bench case
fails when composite < 0.95, merged, or a wrong count. ``default`` must reproduce collect.py.

gurobipy is blocked at import (sys.modules), so this script cannot open a licence session.
"""
import sys

sys.modules["gurobipy"] = None  # any `import gurobipy` now raises ImportError: no licence use possible

import argparse  # noqa: E402
import os  # noqa: E402

import pandas as pd  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_case  # noqa: E402,F401  (puts the RCK root on sys.path)
from rck.ecdna.decompose import pool_keeps, ranking_key  # noqa: E402
from rck.ecdna.evaluate import junction_scores, read_inferred, read_truth_table  # noqa: E402
from rck.ecdna.rerank import load_pool, pool_params  # noqa: E402

FAIL_COMPOSITE = 0.95  # as collect.py
JUNCTION_MATCH = 0.5  # as collect.py
MUS = [0.1, 0.25, 0.5, 1.0, 2.0]
NUS = [0.5, 1.0, 2.0]


class Case:
    """The tables of one case's kept pool, read once."""

    def __init__(self, folder):
        self.pool_dir = os.path.join(folder, "pool")
        pdir = os.path.join(self.pool_dir, "pool")
        self.pool = load_pool(self.pool_dir)
        self.entries = self.pool["entries"]
        self.ev = pd.read_csv(os.path.join(pdir, "evaluation.tsv"), sep="\t")
        self.n_truth = len(self.ev["truth"].unique())
        dec = pd.read_csv(os.path.join(pdir, "decompositions.tsv"), sep="\t")
        self.sp = dec[dec["ordering_rank"] == 1].drop_duplicates(["decomposition_rank", "species"])
        self.inferred = read_inferred(os.path.join(pdir, "decompositions.tsv"))
        self.truth = read_truth_table(os.path.join(pdir, "truth_cycles.tsv"))
        self.jinc, self.explained = {}, {}
        for k, g in self.sp.groupby("decomposition_rank"):
            seen = []
            for j in g["junctions"].fillna("."):
                if j != ".":
                    seen.extend(set(j.split(",")))
            self.jinc[int(k) - 1] = len(seen) - len(set(seen))
            r = g["reads_explained"].iloc[0]
            self.explained[int(k) - 1] = float(r) if pd.notna(r) else 1.0  # no multi-junction reads: no signal

    def composite(self, i):
        m = self.ev[self.ev["decomposition_rank"] == i + 1]["composite"]
        return float(m.mean()) if len(m) else float("nan")

    def n_reported(self, i):
        g = self.sp[self.sp["decomposition_rank"] == i + 1]
        return int((g["status"] != "low_weight").sum())

    def merged(self, i):
        """collect.py::merged_flag for pool decomposition rank i + 1."""
        rank = i + 1
        mine = self.ev[self.ev["decomposition_rank"] == rank]
        scope = set(mine["truth"])
        truth = [t for t in self.truth if t.name in scope]
        cycles = self.inferred.get(rank, {"cycles": []})["cycles"]
        carries = {c.name: [t.name for t in truth
                            if junction_scores(t.steps, c.steps, 10_000, False, t.circular, c.circular)[0] >= JUNCTION_MATCH]
                   for c in cycles}
        if any(len(ts) >= 2 for ts in carries.values()):
            return True
        missed = list(mine.loc[mine["inferred"].isna(), "truth"])
        return any(t in ts for t in missed for ts in carries.values())


def ranked(case, params, extra=None, prune=True, paths_first=False):
    """Pool index of rank 1: pool_keeps (unless prune is False), then ranking_key with ``extra(i)``
    added to the structure cost, or with -reads_explained before it (paths_first)."""
    entries = case.entries
    best = entries[0]
    kept = [e for e in entries if not prune or len(entries) <= 1
            or pool_keeps(e["objective"], e["steps"], best["objective"], best["steps"], params)]

    def key(e):
        k = ranking_key(e["residual"], e["n_significant"], e["n_species"], e["steps"], e["objective"], params)
        i = int(e["index"])
        if paths_first:
            return (k[0], -round(case.explained.get(i, 1.0), 6)) + k[1:]
        if extra is not None:
            return (k[0], round(k[1] + extra(i), 6)) + k[2:]
        return k

    return int(sorted(kept, key=key)[0]["index"])


def rules(case, lam):
    params = pool_params(case.pool, lam)
    out = {"default": ranked(case, params)}
    for mu in MUS:
        f = lambda i, mu=mu: mu * case.jinc.get(i, 0)  # noqa: E731
        out["jinc %g" % mu] = ranked(case, params, f)
        out["jinc %g np" % mu] = ranked(case, params, f, prune=False)
    out["paths first"] = ranked(case, params, paths_first=True)
    for nu in NUS:
        out["paths sum %g" % nu] = ranked(case, params, lambda i, nu=nu: nu * (1.0 - case.explained.get(i, 1.0)))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rck_runs")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.2)
    ap.add_argument("-o", "--out", help="per case x rule TSV")
    a = ap.parse_args()
    rows = []
    for name in sorted(os.listdir(a.rck_runs)):
        folder = os.path.join(a.rck_runs, name)
        if not os.path.exists(os.path.join(folder, "pool", "pool", "pool.json")):
            continue
        case = Case(folder)
        bench = not name.startswith("demo")
        for rule, i in rules(case, a.lam).items():
            comp, merged, nrep = case.composite(i), case.merged(i), case.n_reported(i)
            fail = bench and (not comp >= FAIL_COMPOSITE or merged or nrep != case.n_truth)
            rows.append({"case": name, "bench": bench, "rule": rule, "pool_index": i, "composite": round(comp, 4),
                         "merged": merged, "n_reported": nrep, "n_truth": case.n_truth, "fail": fail,
                         "jinc": case.jinc.get(i, 0), "reads_explained": case.explained.get(i, 1.0)})
    t = pd.DataFrame(rows)
    if a.out:
        t.to_csv(a.out, sep="\t", index=False)
    base = t[t["rule"] == "default"].set_index("case")
    print("λ = %g, %d cases (%d bench)\n" % (a.lam, len(base), int(base["bench"].sum())))
    print("%-16s %8s %6s  %-4s %-4s  %s" % ("rule", "mean", "fails", "up", "down", "fixed (+) / broken (-), cases moving > 0.005"))
    for rule, g in t.groupby("rule", sort=False):
        g = g.set_index("case")
        d = g["composite"] - base["composite"]
        up, down = d[d > 0.005].index.tolist(), d[d < -0.005].index.tolist()
        fixed = [c for c in g.index if base.at[c, "fail"] and not g.at[c, "fail"]]
        broken = [c for c in g.index if g.at[c, "fail"] and not base.at[c, "fail"]]
        print("%-16s %8.4f %6d  %-4d %-4d  + %s | - %s | down: %s" % (
            rule, g["composite"].mean(), int(g["fail"].sum()), len(up), len(down),
            ", ".join(fixed) or "none", ", ".join(broken) or "none", ", ".join(down) or "none"))
    print("\ndefault failures: " + ", ".join(base.index[base["fail"]]))
    print("gurobipy loaded: %s" % (sys.modules.get("gurobipy") is not None))


if __name__ == "__main__":
    main()

"""Offline re-ranking of every saved bench pool under candidate rules (no Gurobi, no new runs).

    python diag/grow_rules.py RCK_RUNS [--lambda 0.2] [-o rules.tsv]

For each case with a kept pool (RCK_RUNS/<case>/pool/pool/pool.json), the pool is re-ranked with
RCK's own rank_pool at λ after changing what the ranking reads, and the rank-1 mean composite
(pool/evaluation.tsv, bench scoring: low_weight species left out) is reported per rule:

* ``default``   the run's own flags (must reproduce collect.py's λ column).
* ``refit``     significant species, all species and significant steps as the output reports them
                (pool/decompositions.tsv: refitted weights re-flagged, refit drops removed), instead
                of the integer flags the run's ranking reads. pool_keeps still sees these steps.
* ``global f``  low weight below max(1, f x heaviest) (RCK's --min-species-weight-frac f).
* ``mass``      the default flags, but the low-weight species' copies count as
                floor(sum of low-weight copies / threshold) extra significant species, where
                threshold = max(1, 0.1 x heaviest): leftovers that together weigh a species are
                not free. Steps are unchanged.
* ``set k``     not a ranking: the best composite over rank 1 and its first k recorded same-fit
                alternatives (pool.json ``same_fit``), each matched to its pool entry by objective,
                steps and significant weights. What step a ("ambiguity reporting") would score.
"""
import argparse
import json
import math
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_case  # noqa: E402,F401
from rck.ecdna.decompose import rank_pool  # noqa: E402
from rck.ecdna.rerank import load_pool, pool_params  # noqa: E402


def flags(weights, frac):
    thr = max(1.0, frac * max(weights))
    return [len(weights) > 1 and w < thr for w in weights], thr


def per_species_steps(dec):
    """{pool index: [steps of species in weight order]} from pool/decompositions.tsv (ordering 1)."""
    t = dec[dec["ordering_rank"] == 1].drop_duplicates(["decomposition_rank", "species"])
    out = {}
    for k, g in t.groupby("decomposition_rank"):
        g = g.sort_values(["weight_int", "n_steps"], ascending=[False, False])
        out[int(k) - 1] = list(zip(g["weight_int"], g["n_steps"]))
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rck_runs")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.2)
    ap.add_argument("--fracs", type=float, nargs="*", default=[0.05, 0.075, 0.15])
    ap.add_argument("-o", "--out")
    a = ap.parse_args()
    rows = []
    for case in sorted(os.listdir(a.rck_runs)):
        pdir = os.path.join(a.rck_runs, case, "pool")
        if not os.path.exists(os.path.join(pdir, "pool", "pool.json")):
            continue
        pool = load_pool(pdir)
        entries = pool["entries"]
        params = pool_params(pool, a.lam)
        ev = pd.read_csv(os.path.join(pdir, "pool", "evaluation.tsv"), sep="\t")
        mean = ev.groupby("decomposition_rank")["composite"].mean()
        comp = lambda i: float(mean.get(i + 1, float("nan")))  # noqa: E731
        dec = pd.read_csv(os.path.join(pdir, "pool", "decompositions.tsv"), sep="\t")
        sps = per_species_steps(dec)
        row = {"case": case}

        top = rank_pool(entries, params)[0]
        row["default"] = comp(top)

        # leftover mass
        mod = []
        for e in entries:
            e2 = dict(e)
            thr = max(1.0, 0.1 * max(e["weights_int"]))
            low = sum(w for w, f in zip(e["weights_int"], e["low_weight"]) if f)
            e2["n_significant"] = e["n_significant"] + int(math.floor(low / thr + 1e-9))
            mod.append(e2)
        row["mass"] = comp(rank_pool(mod, params)[0])

        # refit flags: rank on the significance the output reports (refitted weights, refit drops)
        t1 = dec[dec["ordering_rank"] == 1].drop_duplicates(["decomposition_rank", "species"])
        mod = []
        for e in entries:
            e2 = dict(e)
            g = t1[t1["decomposition_rank"] == e["index"] + 1]
            if len(g):
                sig = g[g["status"] != "low_weight"]
                e2["n_significant"] = len(sig)
                e2["n_species"] = len(g)
                e2["steps"] = int(sig["n_steps"].sum())
            mod.append(e2)
        top_refit = rank_pool(mod, params)[0]
        row["refit"] = comp(top_refit)
        nsig = lambda i: int((t1[t1["decomposition_rank"] == i + 1]["status"] != "low_weight").sum())  # noqa: E731
        row["n_truth"] = int((ev["decomposition_rank"] == top + 1).sum())
        row["default_nsig_reported"] = nsig(top)
        row["refit_nsig_reported"] = nsig(top_refit)

        # global thresholds (significant steps recomputed from per-species steps when available)
        for f in a.fracs:
            mod = []
            for e in entries:
                e2 = dict(e)
                low, _ = flags(e["weights_int"], f)
                e2["low_weight"] = low
                e2["n_significant"] = sum(1 for x in low if not x)
                sp = sps.get(e["index"])
                if sp is not None and len(sp) == len(e["weights_int"]):
                    l2, _ = flags([w for w, _ in sp], f)
                    e2["steps"] = sum(s for (w, s), x in zip(sp, l2) if not x)
                mod.append(e2)
            row["global %g" % f] = comp(rank_pool(mod, params)[0])

        # set metric: rank 1 and its first k same-fit alternatives
        alts = pool["same_fit"].get(str(top)) or []
        sig = lambda e: sorted(w for w, x in zip(e["weights_int"], e["low_weight"]) if not x)  # noqa: E731
        matched = []
        for alt in alts:
            hit = [e["index"] for e in entries if abs(e["objective"] - alt["objective"]) < 1e-3
                   and e["steps"] == alt["steps"] and sig(e) == sorted(alt["weights"])]
            matched.append(hit[0] if hit else None)
        for k in (1, 2, 5):
            vals = [row["default"]] + [comp(i) for i in matched[:k] if i is not None]
            row["set %d" % k] = max(vals)
        row["n_alts"] = len(alts)
        row["alts_matched"] = sum(1 for i in matched if i is not None)
        row["pool_best"] = float(mean.max())
        rows.append(row)
    t = pd.DataFrame(rows)
    if a.out:
        t.to_csv(a.out, sep="\t", index=False)
    cols = [c for c in t.columns if c not in ("case", "n_alts", "alts_matched", "n_truth")
            and not c.endswith("_reported")]
    with pd.option_context("display.width", 250, "display.max_rows", 100):
        print(t.round(3).to_string(index=False))
        bench = ~t["case"].str.startswith("demo")
        for r in ("default", "refit"):
            fail = bench & ((t[r] < 0.95) | (t[r + "_nsig_reported"] != t["n_truth"]))
            print("%s: %d of %d bench cases below 0.95 or with a wrong reported species count (merged flag not checked): %s"
                  % (r, fail.sum(), bench.sum(), ", ".join(t.loc[fail, "case"])))
        print("\nmean over %d cases:" % len(t))
        print(t[cols].mean().round(4).to_string())
        base = t["default"]
        for c in cols:
            if c in ("default", "pool_best"):
                continue
            up = t.loc[t[c] > base + 0.005, "case"].tolist()
            down = t.loc[t[c] < base - 0.005, "case"].tolist()
            print("%-12s up %2d down %2d | down: %s" % (c, len(up), len(down), ", ".join(down)))


if __name__ == "__main__":
    main()

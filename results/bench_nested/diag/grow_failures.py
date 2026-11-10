"""Where the grown bench's failures come from, from the saved pools alone (no Gurobi, no new runs).

    python diag/grow_failures.py RCK_RUNS CASE [CASE ...] [--lambda 0.2] [-o out.tsv] [--species out_species.tsv]

Per case, at λ (default 0.2), with the run's own ranking (rck.ecdna.decompose.rank_pool):
- rank 1: pool index, integer weights and low_weight flags, refit drops (summary.json ``pruned``),
  composites as the bench scores them (low_weight species left out) and with low_weight species
  included (``evaluate(include_low_weight=True)``);
- for each truth species, the inferred species of rank 1 (significant or low_weight) that overlap
  it best (haplotype-aware Jaccard), and the copies of every rank-1 species whose best truth match
  it is;
- the pool's best solution by the bench composite: index, weights, flags, structure cost at λ
  against rank 1's, whether the λ pruning keeps it, and whether it is a recorded same-fit
  alternative of rank 1;
- how many pool solutions have each truth species at composite ≥ 0.9 (bench scoring) and ≥ 0.9
  with low_weight species included.

``--species`` writes one row per rank-1 species (weight, status, best truth, Jaccard).
"""
import argparse
import json
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.dirname(HERE))
import run_case  # noqa: E402,F401  (puts the RCK root on sys.path)
from rck.ecdna.decompose import pool_keeps, rank_pool  # noqa: E402
from rck.ecdna.evaluate import evaluate, jaccard, read_inferred, read_truth_table  # noqa: E402
from rck.ecdna.rerank import load_pool, pool_params  # noqa: E402


def comps(frame, k):
    g = frame[frame["decomposition_rank"] == k]
    return {r.truth: float(r.composite) for r in g.itertuples()}


def fmt(d):
    return ",".join("%s=%.3f" % (t, c) for t, c in sorted(d.items()))


def case_row(run, lam, species_rows):
    name = os.path.basename(run.rstrip("/"))
    pdir = os.path.join(run, "pool")
    pool = load_pool(pdir)
    entries = pool["entries"]
    params = pool_params(pool, lam)
    order = rank_pool(entries, params)
    top = entries[order[0]]
    k1 = top["index"] + 1
    truth = read_truth_table(os.path.join(pdir, "truth_cycles.tsv"))
    dec = os.path.join(pdir, "pool", "decompositions.tsv")
    bench = pd.read_csv(os.path.join(pdir, "pool", "evaluation.tsv"), sep="\t")
    withlw, _ = evaluate(dec, truth, include_low_weight=True)
    inferred = read_inferred(dec, include_low_weight=True)
    table = pd.read_csv(dec, sep="\t")
    table = table[table["ordering_rank"] == 1]

    # rank-1 species and their best truth match
    mine = table[table["decomposition_rank"] == k1].drop_duplicates("species")
    cycles = {c.name: c for c in inferred[k1]["cycles"]}
    per_truth = {t.name: [] for t in truth}
    for r in mine.itertuples():
        c = cycles.get(r.species)
        js = {t.name: jaccard(t.steps, c.steps) for t in truth} if c is not None else {}
        best = max(js, key=js.get) if js else None
        if best is not None and js[best] > 0:
            per_truth[best].append("%s:%d%s" % (r.species, r.weight_int, "L" if r.status == "low_weight" else ""))
        species_rows.append({"case": name, "rank1_pool_index": top["index"], "species": r.species, "status": r.status,
                             "weight_int": r.weight_int, "weight_refit": r.weight, "n_steps": r.n_steps,
                             "best_truth": best, "jaccard_hap": round(js[best], 3) if best else None,
                             "jaccards": ",".join("%s=%.2f" % kv for kv in sorted(js.items()))})

    # refit drops of the rank-1 solution: summary.json holds the run's own rank list (pool λ)
    pruned = None
    for d in json.load(open(os.path.join(pdir, "summary.json")))["decompositions"]:
        if d.get("objective") is not None and abs(d["objective"] - top["objective"]) < 1e-4 and d.get("steps") == top["steps"]:
            pruned = d.get("pruned")
            break

    # pool's best by bench composite
    means = bench.groupby("decomposition_rank")["composite"].mean()
    kb = int(means.idxmax())
    best = entries[kb - 1]
    cost = lambda e: params.species_penalty * e["n_significant"] + lam * e["steps"]  # noqa: E731
    kept = pool_keeps(best["objective"], best["steps"], entries[0]["objective"], entries[0]["steps"], params)
    alts = pool["same_fit"].get(str(top["index"])) or []
    alt_ws = [sorted(a.get("weights", [])) for a in alts]  # significant weights only
    sig_ws = lambda e: sorted(w for w, f in zip(e["weights_int"], e["low_weight"]) if not f)  # noqa: E731
    rank_in = order.index(best["index"]) + 1 if best["index"] in order else None

    good = {t.name: int(sum(1 for k in means.index if comps(bench, k).get(t.name, 0) >= 0.9)) for t in truth}
    good_lw = {t.name: int(sum(1 for k in means.index if comps(withlw, k).get(t.name, 0) >= 0.9)) for t in truth}
    return {
        "case": name, "truth_weights": ",".join("%s=%.1f" % (t.name, t.weight) for t in truth),
        "pool_size": len(entries), "rank1_index": top["index"],
        "rank1_weights": ",".join(str(w) for w in top["weights_int"]),
        "rank1_low": ",".join("L" if f else "-" for f in top["low_weight"]),
        "rank1_sig": top["n_significant"], "rank1_steps": top["steps"], "rank1_cost": round(cost(top), 3),
        "rank1_refit_pruned": ";".join("w%d->%.2f" % (p["weight_int"], p["refit_weight"]) for p in pruned) if pruned else "",
        "rank1_composites": fmt(comps(bench, k1)), "rank1_mean": round(means[k1], 4),
        "rank1_composites_with_low": fmt(comps(withlw, k1)),
        "rank1_species_by_truth": " | ".join("%s<-%s" % (t, "+".join(v) or "none") for t, v in per_truth.items()),
        "best_index": best["index"], "best_mean": round(means[kb], 4), "best_composites": fmt(comps(bench, kb)),
        "best_weights": ",".join(str(w) for w in best["weights_int"]),
        "best_low": ",".join("L" if f else "-" for f in best["low_weight"]),
        "best_sig": best["n_significant"], "best_steps": best["steps"], "best_residual": best["residual"],
        "best_cost": round(cost(best), 3), "best_kept": kept, "best_rank_at_lambda": rank_in,
        "best_is_same_fit_alt": sig_ws(best) in alt_ws,
        "n_same_fit_alts": len(alts),
        "pool_ge0.9": fmt(good).replace(".000", ""), "pool_ge0.9_with_low": fmt(good_lw).replace(".000", ""),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("rck_runs")
    ap.add_argument("cases", nargs="+")
    ap.add_argument("--lambda", dest="lam", type=float, default=0.2)
    ap.add_argument("-o", "--out")
    ap.add_argument("--species")
    a = ap.parse_args()
    species = []
    rows = [case_row(os.path.join(a.rck_runs, c), a.lam, species) for c in a.cases]
    t = pd.DataFrame(rows)
    if a.out:
        t.to_csv(a.out, sep="\t", index=False)
    if a.species:
        pd.DataFrame(species).to_csv(a.species, sep="\t", index=False)
    with pd.option_context("display.width", 250, "display.max_colwidth", 120):
        for r in rows:
            print("==", r["case"])
            for k, v in r.items():
                if k != "case":
                    print("   %-26s %s" % (k, v))


if __name__ == "__main__":
    main()

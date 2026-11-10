"""Offline test of a noise-aware leftover prior (NEXT_STEPS.md step 2b) on kept decomposition pools.

    python leftover_recount.py RUN_DIR [RUN_DIR ...] [--lambda 0.2] [--sigmas 0.1 0.15 0.2]
                               [--global-fracs 0.15 0.2] [--min-overlap 0 0.5] > recount.tsv

Each RUN_DIR is an ``rck-ecdna --keep-pool`` output folder, as for ``python -m rck.ecdna.rerank``
(e.g. ``runs/<case>/pool``, with ``pool/{pool.json,decompositions.tsv,evaluation.tsv}``). Nothing is solved: every pool solution's species
are re-flagged as significant or ``low_weight`` under a rule, the pool's ``n_significant`` and
``steps`` are recomputed from the significant species, and the pool is re-ranked with RCK's own
:func:`rck.ecdna.decompose.rank_pool` at ``--lambda``. One row per run and rule: the rank-1 pool
solution, its integer weights (low-weight ones in brackets) and its mean composite, against the
pool's best.

Rules (integer weights, as the run's ranking reads them; a single species is always significant):

* ``global f``: low weight below max(1, f x the heaviest species), RCK's ``--min-species-weight-frac``
  (``global 0.1`` is the default and must reproduce the run's rank 1).
* ``overlap σ / hap|total / o``: low weight below max(1, σ x the heaviest *other* species that
  the species overlaps). Species s overlaps t when the segments they share cover more than a
  fraction o of s's own length (o = 0: any shared segment). ``hap`` compares phased segments
  (chromosome, interval, haplotype label); ``total`` ignores the label, so a species on haplotype B
  overlaps one on A over the same interval. A species that overlaps no other species is
  significant (weight ≥ 1). A species the weight refit dropped has no saved walk: ``--pruned all``
  (default) lets it overlap every species, so it is judged against the heaviest one, ``none`` lets
  it overlap none; a difference between the two shows where the missing walks matter.
* ``overlap+global``: low weight under either rule (the overlap rule on top of the default 0.1).
* ``novel σ``: the default 0.1 rule, and also low weight below σ x the heaviest species when the
  species uses no junction that a heavier species does not already use (it adds no structure, only
  copies on top of heavier walks). ``novel 1`` flags every lighter species that adds no junction.
"""
import argparse
import json
import os
import sys
from dataclasses import replace

import pandas as pd

sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.abspath(__file__)), "..", "..", "..")))
from rck.ecdna.decompose import rank_pool  # noqa: E402
from rck.ecdna.rerank import pool_params  # noqa: E402


def parse_walk(walk):
    """{(chrom, start, end, hap): bp} of the segments a walk traverses (each counted once)."""
    segments = {}
    for step in str(walk).split(","):
        parts = step.split(":")
        if len(parts) < 4:
            continue
        chrom, interval, hap = ":".join(parts[:-3]), parts[-3], parts[-2]
        start, end = (int(x) for x in interval.split("-"))
        segments[(chrom, start, end, hap)] = end - start + 1
    return segments


def load(run):
    """The pool, its solutions as lists of species {w, steps, seg}, and the mean composite per pool rank.

    The pool's decompositions.tsv holds the species that survive the weight refit; a species the
    refit dropped (``pruned``) is in pool.json's ``weights_int`` only, without its walk. It gets
    ``seg`` None, and its steps from the pool's recorded ``steps`` (the significant species' steps
    under the run's own flags) less those of the surviving significant species, shared equally
    among the dropped species the run counted significant (0 for one the run flagged low weight)."""
    pool = json.load(open(os.path.join(run, "pool", "pool.json")))
    dec = pd.read_csv(os.path.join(run, "pool", "decompositions.tsv"), sep="\t")
    dec = dec[dec["ordering_rank"] == 1]
    ev = pd.read_csv(os.path.join(run, "pool", "evaluation.tsv"), sep="\t")
    means = ev.groupby("decomposition_rank")["composite"].mean()
    frac = float(pool["params"].get("min_species_weight_frac", 0.10))
    solutions = []
    for entry in pool["entries"]:
        rows = dec[dec["decomposition_rank"] == entry["index"] + 1]
        species = [{"w": int(r.weight_int), "steps": int(r.n_steps), "seg": parse_walk(r.walk),
                    "junctions": set(str(r.junctions).split(",")) - {"nan", "."}} for r in rows.itertuples()]
        missing = [int(w) for w in entry["weights_int"]]
        for sp in species:
            if sp["w"] not in missing:
                raise ValueError("{r}: pool index {i}: species do not match pool.json".format(r=run, i=entry["index"]))
            missing.remove(sp["w"])
        if missing:
            threshold = max(1.0, frac * max(int(w) for w in entry["weights_int"])) if len(entry["weights_int"]) > 1 else 0
            known = sum(sp["steps"] for sp in species if sp["w"] >= threshold)
            sig_missing = [w for w in missing if w >= threshold]
            share = (entry["steps"] - known) / float(len(sig_missing)) if sig_missing else 0
            species += [{"w": w, "steps": share if w >= threshold else 0, "seg": None, "junctions": None} for w in missing]
        solutions.append(species)
    return pool, solutions, means


def overlaps(s, t, mode, min_overlap, pruned):
    """Whether s overlaps t (module docstring). A species without a walk (dropped by the refit)
    overlaps every species with ``pruned`` "all" and none with "none"."""
    if s["seg"] is None or t["seg"] is None:
        return pruned == "all"
    if mode == "hap":
        a, b = s["seg"], set(t["seg"])
    else:
        a = {}
        for (c, x, y, _h), bp in s["seg"].items():
            a[(c, x, y)] = bp
        b = {(c, x, y) for (c, x, y, _h) in t["seg"]}
    shared = sum(bp for k, bp in a.items() if k in b)
    total = sum(a.values())
    return shared > 0 and shared > min_overlap * total


def adds_junction(i, species, pruned):
    """Whether species i uses a junction that no heavier species uses (heavier: larger integer
    weight, or equal weight and earlier). A species without a saved walk adds none with
    ``pruned`` "all" and one with "none"; it does not count as a heavier species' junctions."""
    s = species[i]
    if s["junctions"] is None:
        return pruned == "none"
    heavier = set()
    for j, t in enumerate(species):
        if j != i and t["junctions"] is not None and (t["w"] > s["w"] or (t["w"] == s["w"] and j < i)):
            heavier |= t["junctions"]
    return bool(s["junctions"] - heavier)


def low_weight(species, rule):
    """Per species, whether ``rule`` flags it low weight."""
    if len(species) <= 1:
        return [False] * len(species)
    heaviest = max(s["w"] for s in species)
    flags = []
    for i, s in enumerate(species):
        threshold = 1.0
        if rule["kind"] in ("global", "overlap+global", "novel"):
            threshold = max(threshold, rule.get("frac", 0.10) * heaviest)
        if rule["kind"] in ("overlap", "overlap+global"):
            heavier = [t["w"] for j, t in enumerate(species) if j != i and overlaps(s, t, rule["mode"], rule["min_overlap"], rule["pruned"])]
            if heavier:
                threshold = max(threshold, rule["sigma"] * max(heavier))
        if rule["kind"] == "novel" and not adds_junction(i, species, rule["pruned"]):
            threshold = max(threshold, rule["sigma"] * heaviest)
        flags.append(s["w"] < threshold)
    return flags


def recount(pool, solutions, means, rule, lam):
    entries = []
    for entry, species in zip(pool["entries"], solutions):
        flags = low_weight(species, rule)
        sig = [s for s, f in zip(species, flags) if not f]
        entries.append(dict(entry, n_significant=len(sig), steps=sum(s["steps"] for s in sig),
                            _weights=",".join(("[{w}]" if f else "{w}").format(w=s["w"])
                                              for s, f in sorted(zip(species, flags), key=lambda x: -x[0]["w"]))))
    params = replace(pool_params(pool, lam), rank_residual=pool["params"].get("rank_residual", "first"))
    order = rank_pool(entries, params)
    top = entries[order[0]]
    return {"pool_index": top["index"], "n_significant": top["n_significant"], "weights": top["_weights"],
            "steps": top["steps"], "residual": top["residual"],
            "mean_composite": round(float(means.get(top["index"] + 1, float("nan"))), 4)}


def rules(args):
    out = [{"kind": "global", "frac": 0.10, "label": "global 0.1 (default)"}]
    out += [{"kind": "global", "frac": f, "label": "global {f:g}".format(f=f)} for f in args.global_fracs]
    for kind in ("overlap", "overlap+global"):
        for mode in ("hap", "total"):
            for o in args.min_overlap:
                for sigma in args.sigmas:
                    out.append({"kind": kind, "mode": mode, "min_overlap": o, "sigma": sigma,
                                "pruned": args.pruned,
                                "label": "{k} {s:g} / {m} / {o:g}".format(k=kind, s=sigma, m=mode, o=o)})
    for sigma in args.novel_sigmas:
        out.append({"kind": "novel", "sigma": sigma, "pruned": args.pruned, "label": "novel {s:g}".format(s=sigma)})
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("runs", nargs="+")
    parser.add_argument("--lambda", dest="lam", type=float, default=0.2)
    parser.add_argument("--sigmas", type=float, nargs="+", default=[0.1, 0.15, 0.2])
    parser.add_argument("--global-fracs", type=float, nargs="*", default=[0.15, 0.2])
    parser.add_argument("--min-overlap", type=float, nargs="+", default=[0.0, 0.5])
    parser.add_argument("--novel-sigmas", type=float, nargs="*", default=[0.2, 0.3, 0.5, 1.0])
    parser.add_argument("--pruned", choices=["all", "none"], default="all",
                        help="Species the refit dropped (no walk saved) overlap every species (all) or none")
    args = parser.parse_args()
    rows = []
    for run in args.runs:
        pool, solutions, means = load(run)
        name = os.path.basename(os.path.dirname(os.path.normpath(run)))
        best = means.idxmax()
        n_pruned = sum(any(sp["seg"] is None for sp in species) for species in solutions)
        for rule in rules(args):
            rows.append(dict({"run": name, "rule": rule["label"], "best_pool_composite": round(float(means.max()), 4),
                              "best_pool_index": int(best) - 1, "pool_size": len(solutions), "with_pruned": n_pruned}, **recount(pool, solutions, means, rule, args.lam)))
    pd.DataFrame(rows).to_csv(sys.stdout, sep="\t", index=False)


if __name__ == "__main__":
    main()

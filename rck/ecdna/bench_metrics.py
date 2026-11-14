"""Method-agnostic copy-number recall and ecDNA proportion scores of one answer (bench_IW).

Both read an AA cycles file (any method's answer, as ``aa_eval`` does) and the simulator truth.

Copy-number recall (ecDNA only, haplotype-blind, base-pair weighted). Over every reference base x:

    t(x) = sum_s PHI_s * (times truth species s traverses x)        (tumour sample, bulk copies)
    c(x) = sum_c w_c   * (times inferred cycle c traverses x)       (the answer's strict cycles)
    recall    = sum_x min(c, t) / sum_x t
    precision = sum_x min(c, t) / sum_x c
    within_tol = share of truth bases (t > 0) with |c - t| <= max(1, 0.1 t)

Chromosomal copies are left out: they are an input (clone profiles) or fitted background, and
not every method reports them. ``w_c`` must be bulk copies per cell, as ``PHI_s``: the bulk
route's and the graph-route methods' weights are; for the multi-clone route pass
``--multiclone-run`` so w_c = sum_k u_k psi_k (its own ``weight`` is sum_k psi_k).

ecDNA proportions. Truth p_s = PHI_s / sum PHI over the truth species; inferred p_c = w_c / sum w
over the cycles of one set. Truth and inferred are matched one-to-one, greedily by haplotype-blind
Jaccard (``evaluate._match``, as M1 and the strict scores). Per set:

    tv  = 1/2 (sum_matched |p_c - p_s| + sum_unmatched truth p_s + sum_unmatched inferred p_c)
    mae = mean over matched pairs of |p_c - p_s|

The sets are ``strict`` (the answer, ``aa_eval.strict_select``: no ``Status=low_weight`` for
rck-ecdna, cycles >= 0.1 x the heaviest for files without ``Status``) and ``all`` (low_weight
cycles included); the difference shows the mass the method itself discarded. Cycles flagged
``Status=chromosomal`` are never species (``aa_eval.species_cycles``).

    python -m rck.ecdna.bench_metrics --cycles ANSWER_cycles.txt --truth-dir SIM_RUN --outdir OUT
        [--multiclone-run RUN_DIR --clone-proportions DIAG/clone_proportions.tsv] [--case NAME]

Writes OUT/bench_metrics.json and OUT/bench_proportions.tsv (one row per truth or inferred species
and set).
"""
import argparse
import json
import logging
import os
from collections import defaultdict
from typing import Dict, List, Optional, Sequence

import pandas as pd

from rck.ecdna import aa_eval
from rck.ecdna import aa_format as aaf
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import Cycle

log = logging.getLogger(__name__)
TOL_ABS, TOL_REL = 1.0, 0.1


def copy_profile(cycles: Sequence[Cycle]) -> Dict[str, List[tuple]]:
    """Per chromosome, (position, change in copies) events of sum_c w_c * traversals(x); steps
    are 0-based inclusive (``evaluate.CStep``)."""
    events: Dict[str, List[tuple]] = defaultdict(list)
    for c in cycles:
        for s in c.steps:
            lo, hi = min(s.start, s.end), max(s.start, s.end)
            events[s.chrom].append((lo, float(c.weight)))
            events[s.chrom].append((hi + 1, -float(c.weight)))
    return events


def cn_scores(truth: Sequence[Cycle], inferred: Sequence[Cycle]) -> Dict[str, float]:
    """Recall, precision and within-tolerance share of the ecDNA copy profiles (module docstring)."""
    t_ev, c_ev = copy_profile(truth), copy_profile(inferred)
    mass_t = mass_c = overlap = within = support = 0.0
    for chrom in set(t_ev) | set(c_ev):
        points = sorted({p for p, _ in t_ev.get(chrom, [])} | {p for p, _ in c_ev.get(chrom, [])})
        dt, dc = defaultdict(float), defaultdict(float)
        for p, d in t_ev.get(chrom, []):
            dt[p] += d
        for p, d in c_ev.get(chrom, []):
            dc[p] += d
        t = c = 0.0
        for left, right in zip(points, points[1:]):
            t += dt[left]
            c += dc[left]
            t, c = (0.0 if abs(t) < 1e-9 else t), (0.0 if abs(c) < 1e-9 else c)
            width = right - left
            mass_t += t * width
            mass_c += c * width
            overlap += min(t, c) * width
            if t > 0:
                support += width
                if abs(c - t) <= max(TOL_ABS, TOL_REL * t):
                    within += width
    return {"cn_recall": overlap / mass_t if mass_t else None, "cn_precision": overlap / mass_c if mass_c else None,
            "cn_within_tol": within / support if support else None, "truth_mass": mass_t, "inferred_mass": mass_c,
            "overlap_mass": overlap, "truth_bp": support}


def proportion_scores(truth: Sequence[Cycle], inferred: Sequence[Cycle], label: str):
    """(summary, rows) of one inferred set against the truth species (module docstring)."""
    t_total = sum(t.weight for t in truth)
    c_total = sum(c.weight for c in inferred)
    p_t = {t.name: t.weight / t_total if t_total else 0.0 for t in truth}
    p_c = {c.name: c.weight / c_total if c_total else 0.0 for c in inferred}
    match = ev._match(truth, inferred, with_hap=False)
    rows, errors, tv = [], [], 0.0
    for t in truth:
        c = match.get(t.name)
        err = abs(p_c[c] - p_t[t.name]) if c else None
        if c:
            errors.append(err)
            tv += err
        else:
            tv += p_t[t.name]
        rows.append({"set": label, "truth": t.name, "inferred": c, "truth_proportion": p_t[t.name],
                     "inferred_proportion": p_c.get(c) if c else None, "abs_err": err,
                     "jaccard_blind": ev.jaccard(t.steps, next(x for x in inferred if x.name == c).steps, False) if c else None})
    matched = set(match.values())
    for c in inferred:
        if c.name not in matched:
            tv += p_c[c.name]
            rows.append({"set": label, "truth": None, "inferred": c.name, "truth_proportion": None,
                         "inferred_proportion": p_c[c.name], "abs_err": None, "jaccard_blind": None})
    summary = {"proportion_tv": tv / 2.0, "proportion_mae_matched": sum(errors) / len(errors) if errors else None,
               "n_matched": len(errors), "n_inferred": len(inferred), "n_truth": len(truth)}
    return summary, rows


def multiclone_bulk_weights(run_dir: str, proportions: Dict[str, float]) -> Dict[str, float]:
    """Species name -> sum_k u_k psi_k of a multi-clone run's rank-1 decomposition."""
    table = pd.read_csv(os.path.join(run_dir, "decompositions.tsv"), sep="\t")
    table = table[table["decomposition_rank"] == table["decomposition_rank"].min()]
    psi = [c for c in table.columns if c.startswith("psi_")]
    return {str(r["species"]): sum(proportions.get(c[len("psi_"):], 0.0) * float(r[c] or 0.0) for c in psi)
            for _, r in table.iterrows()}


def read_proportions(path: str) -> Dict[str, float]:
    """clone -> u from a DIAG clone_proportions.tsv (columns clone/id and proportion/u)."""
    frame = pd.read_csv(path, sep="\t")
    key = next(c for c in frame.columns if c.lower() in ("clone", "id", "clone_id"))
    value = next(c for c in frame.columns if c.lower() in ("proportion", "u", "fraction"))
    return {str(k): float(v) for k, v in zip(frame[key], frame[value])}


def reweight(cycles: Sequence[Cycle], weights: Dict[str, float]) -> List[Cycle]:
    out = []
    for c in cycles:
        if c.name not in weights:
            raise KeyError("no bulk weight for cycle {n}".format(n=c.name))
        out.append(Cycle(c.name, weights[c.name], c.orderings, c.circular, c.weight_int))
    return out


def score(cycles_path: str, truth: Sequence[Cycle], weights: Optional[Dict[str, float]] = None,
          strict_rule: str = "auto"):
    """(summary, proportion rows) of one answer; ``weights`` overrides cycle weights by name."""
    records = aaf.read_cycles(cycles_path)
    cycles = aa_eval.species_cycles(records)
    strict, rule = aa_eval.strict_select(records, cycles, strict_rule)
    if weights is not None:
        cycles, strict = reweight(cycles, weights), reweight(strict, weights)
    summary = {"strict_rule": rule, "n_truth": len(truth), "n_cycles": len(cycles), "n_strict": len(strict)}
    summary.update(cn_scores(truth, strict))
    summary.update({k + "_all": v for k, v in cn_scores(truth, cycles).items() if k.startswith("cn_")})
    rows = []
    for label, chosen in (("strict", strict), ("all", cycles)):
        s, r = proportion_scores(truth, chosen, label)
        summary.update({"{k}{x}".format(k=k, x="" if label == "strict" else "_all"): v for k, v in s.items()
                        if k.startswith("proportion_") or k == "n_matched"})
        rows += r
    return summary, rows


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rck.ecdna.bench_metrics", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cycles", required=True, help="the method's AA cycles file")
    parser.add_argument("--truth-dir", required=True, help="simulator run directory")
    parser.add_argument("--truth-sample", default=None)
    parser.add_argument("--multiclone-run", default=None, help="multi-clone run: weights become sum_k u_k psi_k")
    parser.add_argument("--clone-proportions", default=None, help="DIAG clone_proportions.tsv (with --multiclone-run)")
    parser.add_argument("--strict-rule", default="auto")
    parser.add_argument("--case", default=None)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    truth = ev.load_truth(args.truth_dir, None, args.truth_sample)
    weights = None
    if args.multiclone_run:
        if not args.clone_proportions:
            parser.error("--multiclone-run needs --clone-proportions")
        weights = multiclone_bulk_weights(args.multiclone_run, read_proportions(args.clone_proportions))
    summary, rows = score(args.cycles, truth, weights, args.strict_rule)
    summary["case"] = args.case
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "bench_metrics.json"), "w") as handle:
        json.dump(summary, handle, indent=2)
    pd.DataFrame(rows).assign(case=args.case).to_csv(os.path.join(args.outdir, "bench_proportions.tsv"), sep="\t", index=False)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

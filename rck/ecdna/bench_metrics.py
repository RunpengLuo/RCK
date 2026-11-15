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

Total copy-number recall (the two RCK routes only; Yubi, 2026-11-14). Background plus ecDNA,
haplotype-blind, base-pair weighted, over D = the bases of the union subgraph's segments (the DIAG's
``kind=flank`` segments left out):

    t(x) = CN_A + CN_B of truth/segments.tsv: the tumour sample's cell-averaged copy number, normal
           cells included (the frame the subgraph's depth-derived copy numbers are in)
    c(x) = rck_bulk: inferred_total of flow_segments.tsv (fitted background + ecDNA)
           rck_mc:   sum_k u_k (A_k + B_k) over the clones of its output rck.scnt.tsv, plus the normal
                     cells' (1 - sum_k u_k) x 2, with u from the DIAG's clone_proportions.tsv
    total_cn_recall = sum_D min(c, t) / sum_D t;  total_cn_precision = sum_D min(c, t) / sum_D c;
    total_cn_within_tol = share of D with |c - t| <= max(1, 0.1 t);  total_cn_mae = mean_D |c - t|
    total_cn_ecdna_*: the same over the bases of D that a truth species traverses

D is mostly background (tens of Mb around a few Mb of ecDNA), and t is cell-averaged: with purity p,
a background where the tumour has 1 copy reads 1 + (1 - p) or similar fractions, which one
whole-number bulk karyotype cannot match. total_cn_* over D is therefore dominated by background and
purity; total_cn_ecdna_* is the part where the ecDNA copies matter.

truth/segments.tsv is 0-based half-open ([START, END)); the methods' segments are in the walk frame,
which is this module's 0-based inclusive frame. CoRAL and CycleExtractor report no segment copy
numbers of their own, so they get no total score.

What the multi-clone total measures: the route is given each clone's allele copy numbers (the true
clone profiles with psi folded in, --clone-profiles truth) as RCK's segment copy numbers, which it
may move within their bounds (cnb) while it infers the haplotype labelling and adjacency copy
numbers. On bench_IW (array 32267060) it kept them on 99.98% of the (base, clone) pairs: 203 changed
cells were 1-bp slivers, and 28 longer ones were all in hard_S4_s16003, where about 322 kb per clone
of a 78-copy region (clone1) came out at 31. So its total recall is mostly a check of its input, the
cell-averaged frame and the normal share, but it does catch the cases where the model departs from
the profiles it was given.
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


def _intervals_total(rows) -> Dict[str, List[tuple]]:
    """[(chrom, lo, hi, value)] (0-based inclusive) -> per chromosome sorted, non-overlapping."""
    out: Dict[str, List[tuple]] = defaultdict(list)
    for chrom, lo, hi, value in rows:
        out[str(chrom)].append((int(lo), int(hi), float(value)))
    for chrom in out:
        out[chrom].sort()
        for (a_lo, a_hi, _), (b_lo, _b, _v) in zip(out[chrom], out[chrom][1:]):
            if b_lo <= a_hi:
                raise ValueError("overlapping segments on {c}: {a}-{b} and {d}".format(c=chrom, a=a_lo, b=a_hi, d=b_lo))
    return out


def truth_total_cn(truth_dir: str) -> List[tuple]:
    """truth/segments.tsv as [(chrom, lo, hi, CN_A + CN_B)], 0-based inclusive."""
    t = pd.read_csv(os.path.join(truth_dir, "truth", "segments.tsv"), sep="\t")
    return [(str(r["#CHR"]), int(r["START"]), int(r["END"]) - 1, float(r["CN_A"]) + float(r["CN_B"])) for r in t.to_dict("records")]


def bulk_total_cn(run_dir: str) -> List[tuple]:
    """rck_bulk's fitted total copy number per segment (flow_segments.tsv inferred_total)."""
    t = pd.read_csv(os.path.join(run_dir, "flow_segments.tsv"), sep="\t")
    return [(str(r["chrom"]), int(r["start"]), int(r["end"]), float(r["inferred_total"])) for r in t.to_dict("records")]


def multiclone_total_cn(run_dir: str, proportions: Dict[str, float], normal_cn: float = 2.0) -> List[tuple]:
    """sum_k u_k (A_k + B_k) + (1 - sum_k u_k) x normal_cn per segment of the multi-clone run's
    rck.scnt.tsv, its ``kind=flank`` segments left out (they are outside the subgraph)."""
    import ast
    import re

    rows = []
    unknown = set()
    normal = max(0.0, 1.0 - sum(proportions.values())) * normal_cn
    with open(os.path.join(run_dir, "rck.scnt.tsv")) as handle:
        next(handle)
        for line in handle:
            chrom, lo, hi, extra = line.rstrip("\n").split("\t")[:4]
            if re.search(r"(^|;)kind=flank(;|$)", extra):
                continue
            cn = ast.literal_eval(re.search(r"cn=(\{.*?\}\})", extra).group(1))
            unknown |= set(cn) - set(proportions)
            rows.append((chrom, int(lo), int(hi), normal + sum(proportions.get(k, 0.0) * (v["A"] + v["B"]) for k, v in cn.items())))
    if unknown:
        raise KeyError("clones {u} have no proportion".format(u=sorted(unknown)))
    return rows


def ecdna_bases(truth: Sequence[Cycle]) -> List[tuple]:
    """The reference bases any truth species traverses, as [(chrom, lo, hi, 1.0)] (0-based inclusive, merged)."""
    spans: Dict[str, List[tuple]] = defaultdict(list)
    for c in truth:
        for st in c.steps:
            spans[st.chrom].append((min(st.start, st.end), max(st.start, st.end)))
    out = []
    for chrom, items in spans.items():
        items.sort()
        lo, hi = items[0]
        for a, b in items[1:]:
            if a <= hi + 1:
                hi = max(hi, b)
            else:
                out.append((chrom, lo, hi, 1.0))
                lo, hi = a, b
        out.append((chrom, lo, hi, 1.0))
    return out


def _restrict(rows: Sequence[tuple], mask: Sequence[tuple]) -> List[tuple]:
    """rows cut to the bases of mask (both [(chrom, lo, hi, value)], 0-based inclusive)."""
    by = _intervals_total(mask)
    out = []
    for chrom, lo, hi, v in rows:
        for m_lo, m_hi, _ in by.get(str(chrom), []):
            a, b = max(int(lo), m_lo), min(int(hi), m_hi)
            if a <= b:
                out.append((chrom, a, b, v))
    return out


def total_cn_all(truth_rows: Sequence[tuple], inferred_rows: Sequence[tuple], truth: Sequence[Cycle]) -> Dict[str, float]:
    """total_cn_* over the subgraph's segments, and total_cn_*_ecdna over the bases a truth species
    traverses (where ecDNA errors are not diluted by megabases of background)."""
    out = total_cn_scores(truth_rows, inferred_rows)
    out.update({k.replace("total_cn_", "total_cn_ecdna_"): v
                for k, v in total_cn_scores(truth_rows, _restrict(inferred_rows, ecdna_bases(truth))).items()})
    return out


def total_cn_scores(truth_rows: Sequence[tuple], inferred_rows: Sequence[tuple]) -> Dict[str, float]:
    """Total copy-number recall / precision / within-tolerance / MAE over the inferred segments' bases
    (module docstring). A base the truth table does not cover is reported, not scored."""
    truth, inferred = _intervals_total(truth_rows), _intervals_total(inferred_rows)
    mass_t = mass_c = overlap = within = err = bases = uncovered = 0.0
    for chrom, segs in inferred.items():
        tsegs = truth.get(chrom, [])
        for lo, hi, c in segs:
            covered = 0
            for t_lo, t_hi, t in tsegs:
                a, b = max(lo, t_lo), min(hi, t_hi)
                if a > b:
                    continue
                n = b - a + 1
                covered += n
                mass_t += n * t
                mass_c += n * c
                overlap += n * min(c, t)
                err += n * abs(c - t)
                within += n * (abs(c - t) <= max(TOL_ABS, TOL_REL * t))
                bases += n
            uncovered += (hi - lo + 1) - covered
    return {"total_cn_recall": overlap / mass_t if mass_t else None,
            "total_cn_precision": overlap / mass_c if mass_c else None,
            "total_cn_within_tol": within / bases if bases else None,
            "total_cn_mae": err / bases if bases else None,
            "total_cn_bases": int(bases), "total_cn_bases_no_truth": int(uncovered)}


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rck.ecdna.bench_metrics", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cycles", required=True, help="the method's AA cycles file")
    parser.add_argument("--truth-dir", required=True, help="simulator run directory")
    parser.add_argument("--truth-sample", default=None)
    parser.add_argument("--multiclone-run", default=None, help="multi-clone run: weights become sum_k u_k psi_k")
    parser.add_argument("--clone-proportions", default=None, help="DIAG clone_proportions.tsv (with --multiclone-run)")
    parser.add_argument("--strict-rule", default="auto")
    parser.add_argument("--bulk-run", default=None, help="rck_bulk run: also total copy-number scores (flow_segments.tsv)")
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
    if args.bulk_run:
        summary.update(total_cn_all(truth_total_cn(args.truth_dir), bulk_total_cn(args.bulk_run), truth))
    elif args.multiclone_run:
        summary.update(total_cn_all(truth_total_cn(args.truth_dir),
                                    multiclone_total_cn(args.multiclone_run, read_proportions(args.clone_proportions)), truth))
    summary["case"] = args.case
    os.makedirs(args.outdir, exist_ok=True)
    with open(os.path.join(args.outdir, "bench_metrics.json"), "w") as handle:
        json.dump(summary, handle, indent=2)
    pd.DataFrame(rows).assign(case=args.case).to_csv(os.path.join(args.outdir, "bench_proportions.tsv"), sep="\t", index=False)
    print(json.dumps(summary, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

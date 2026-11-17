"""RCK-style evaluation of an ecDNA answer: the scores Aganezov & Raphael (2020) used for haplotype-specific
karyotypes, transferred to ecDNA cycles (results/bench_IW_allele/PLAN.md, Part A).

Reference: Aganezov S, Raphael BJ. Reconstruction of clone- and haplotype-specific cancer genome karyotypes
from bulk tumor samples. Genome Res 30:1274-1290 (2020), doi:10.1101/gr.256701.119; "Evaluation of RCK on
simulated data" (Fig. 2), "Heterogeneous tumor samples" (Fig. 4), Methods "Reconstructing karyotypes under
generalized infinite sites constraints" and Eq. 11.

A novel adjacency is a non-contiguous join of a walk (``evaluate.junctions``): (exit breakend, entry
breakend, haplotype of the exit step, haplotype of the entry step); a breakend is (chrom, position,
strand), "+" at a step's right end and "-" at its left end. Answers and truth are compared after the
answer's A / B is swapped once per phase block where that lowers the allele error (``_block``: RCK's
"up to haplotype symmetries"); the absolute scores are reported too. In that within-block frame, a
haplotype switch between reference-contiguous steps exactly at a phase-block boundary is dropped from the
truth's and the answer's adjacencies (``block_boundary_switch``): the phase across a block boundary is not
identifiable, and swapping one block but not its neighbour would otherwise turn a plain traversal into an
A-to-B join at one position (a homologous-reciprocal violation that no method made). The swap is applied
per breakend (``relabel_junctions``), since a bulk-graph segment may cross a block boundary.

1. Novel adjacencies (RCK Fig. 2A). The answer's distinct novel adjacencies are matched one to one to the
   truth's within ``tol`` bp (``aa_eval.adjacency_items``), labelled (``phased``: the haplotype pair must
   agree too) or unlabelled:

       adj_fnr = truth adjacencies left unmatched / truth adjacencies
       adj_fdr = answer adjacencies left unmatched / answer adjacencies

   RCK's FPR is over the input: spurious input adjacencies used / spurious input adjacencies. Here the
   input is the candidate graph (``--graph``, an AA graph): ``cand_fpr`` = candidates that match no truth
   adjacency and are used by the answer / candidates that match no truth adjacency (unlabelled: candidates
   carry no haplotype).

2. Copy-number distance (RCK Eq. 11), bulk level (one genome: the bulk ecDNA copies, as
   ``allele_metrics``). With fragments F (phase blocks) and t / c the true / inferred allele-specific copy
   profiles (``allele_metrics.allele_table``):

       cn_dist = sum_{f in F} min over swapping c_A and c_B in f of
                 sum_{x in f} (|c_A(x) - t_A(x)| + |c_B(x) - t_B(x)| + c_U(x))      [copies x bp]

   ``cn_dist_rel`` divides by sum_x t(x), so cases of different size compare; c_U (unphased copies) are all
   error. Without the minimum (absolute frame) the same sum is ``cn_dist_abs``. RCK sums over clones too;
   per-clone distances for the multi-clone route are PLAN.md item 5.

3. Generalized infinite-sites violations (RCK Fig. 4), over the union of the answer's labelled novel
   adjacencies (RCK imposes the constraints on the whole sample, so violations within and across species
   are both counted). Breakends within ``tol`` on one chromosome and strand are one extremity.

       is_ext       extremities in two or more distinct (unlabelled) novel adjacencies (extremity
                    exclusivity; ``aa_eval.is_violations`` on graph edges)
       is_hom_same  novel adjacencies used with two or more haplotype labellings (Fig. 4A / 4C: "a novel
                    adjacency a be assigned both a label A and a label B")
       is_hom       extremities used on A by one novel adjacency and on B by any (homologous-extremity
                    exclusivity; every is_hom_same adjacency gives one or two)
       is_recip     reference-adjacent extremity pairs (the "+" end at p, the "-" end at p + 1) used on
                    opposite haplotypes (homologous-reciprocal-extremity exclusivity, Fig. 4E)

   Unphased steps (``.``) count for ``is_ext`` only.

4. Per-clone copy numbers (multi-clone route; RCK Eq. 11 as defined, summed over clones, and Supplemental
   Fig. S2B,C). The inferred per-clone karyotype (the run's ``rck.scnt.tsv``) against the truth's clone
   karyotypes (``truth/rck/rck.scnt.tsv``), over the run's input fragments F:

       clone_cn_dist = sum_{f in F} min over swapping A and B in f (every clone alike) of
                       sum_i sum_{x in f} (|a_i(x) - t^A_i(x)| + |b_i(x) - t^B_i(x)|)      [copies x bp]

   over the whole input region (``_all``) and over the truth's ecDNA bases (``_ecdna``), each also relative
   to sum_i sum_x (t^A_i + t^B_i). ``clone_state_agree``: share of (clone, allele, bp) whose state (0 loss,
   1 neutral, >= 2 gain) matches the truth's, after the same per-fragment swaps.

5. Input use (RCK Fig. 3B) and linear walks. ``input_used``: input novel adjacencies the answer uses / all
   input novel adjacencies (multi-clone: the run's ``novel_used``; CoRAL / CE: their answer's adjacencies
   matched to the graph they decompose). ``input_fpr``: input adjacencies matching no truth adjacency
   (``tol``) that are used / all such (RCK's FPR; defined only when the input holds spurious ones).
   ``linear_walks``: the answer's strict walks that are not cycles (open paths), the ecDNA analogue of RCK's
   novel telomeres (Supplemental Figs. S6, S9).

    python -m rck.ecdna.rck_eval --cycles ANSWER_cycles.txt --truth-dir SIM_RUN --outdir OUT
        [--graph GRAPH.txt] [--phase-blocks rck_input/rck.fragments.tsv]
        [--weights-from-multiclone RUN_DIR --clone-proportions DIAG/clone_proportions.tsv]
    python -m rck.ecdna.rck_eval --truth-only --truth-dir SIM_RUN --outdir OUT   # violations of the truth

Writes OUT/rck_eval.json.
"""
from __future__ import annotations

import argparse
import json
import os
from collections import defaultdict
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from rck.ecdna import aa_eval
from rck.ecdna import aa_format as aaf
from rck.ecdna import allele_metrics as am
from rck.ecdna import bench_metrics as bm
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import Cycle

ADJ_TOL = aa_eval.GRAPH_TOL  # bp, as aa_eval's graph-level and adjacency scores
EXT_TOL = 2  # bp: breakends this close on one chromosome and strand are one extremity (evaluate._contiguous)
HAPS = ("A", "B")

Junction = aa_eval.Junction


# ---------------------------------------------------------------------------------------------- #
# novel adjacencies
# ---------------------------------------------------------------------------------------------- #

def novel_adjacencies(cycles: Sequence[Cycle], phased: bool = True) -> List[Junction]:
    """Distinct novel adjacencies over every cycle (labelled when ``phased``)."""
    return aa_eval._distinct([j for c in cycles for j in aa_eval.walk_adjacencies(c, phased)], phased)


def block_boundary_switch(j: Junction, blocks: Sequence[ev.Block], tol: int = EXT_TOL) -> bool:
    """True for a join of reference-contiguous steps (the "+" end at p joined to the "-" end at p + 1, either
    direction) that is not scorable in the within-block frame: one on a phase-block boundary, or one whose two
    ends carry the same haplotype after the block swaps (a plain reference traversal again)."""
    a, b, h1, h2 = j
    if a[0] != b[0] or {a[2], b[2]} != {"+", "-"}:
        return False
    left, right = (a, b) if a[2] == "+" else (b, a)
    if abs(right[1] - (left[1] + 1)) > tol:
        return False
    if h1 == h2:
        return True
    return any(chrom == left[0] and (abs(end - left[1]) <= tol or abs(start - right[1]) <= tol)
               for chrom, start, end, _ in blocks)


def _keep(js: Sequence[Junction], blocks: Optional[Sequence[ev.Block]]) -> List[Junction]:
    return [j for j in js if not (blocks and block_boundary_switch(j, blocks))]


def relabel_junctions(js: Sequence[Junction], blocks: Sequence[ev.Block], swapped: set) -> List[Junction]:
    """Each junction end's haplotype swapped when the phase block holding that breakend is in ``swapped``.
    Done per breakend, not per step (``evaluate.relabel`` uses a step's midpoint), because a bulk-graph
    segment may cross a block boundary and its two ends then lie in different blocks."""
    other = {"A": "B", "B": "A"}

    def block_of(be):
        for chrom, start, end, name in blocks:
            if chrom == be[0] and start <= be[1] <= end:
                return name
        return None

    out = []
    for a, b, h1, h2 in js:
        h1 = other.get(h1, h1) if block_of(a) in swapped else h1
        h2 = other.get(h2, h2) if block_of(b) in swapped else h2
        out.append((a, b, h1, h2))
    return aa_eval._distinct(out, True)


def adjacency_rates(truth: Sequence[Cycle], answer: Sequence[Cycle], phased: bool, tol: int = ADJ_TOL) -> Dict[str, object]:
    """FNR / FDR of the answer's novel adjacencies against the truth's (module docstring, 1)."""
    return junction_rates(novel_adjacencies(truth, phased), novel_adjacencies(answer, phased), phased, tol)


def junction_rates(t_js: Sequence[Junction], a_js: Sequence[Junction], phased: bool, tol: int = ADJ_TOL) -> Dict[str, object]:
    """FNR / FDR of distinct adjacencies ``a_js`` against ``t_js``, matched one to one, nearest first,
    within ``tol`` (as ``aa_eval.adjacency_items``)."""
    pairs = sorted((d, i, k) for i, j in enumerate(a_js) for k, u in enumerate(t_js)
                   if (d := aa_eval.junction_distance(j, u, phased)) is not None and d <= tol)
    used_a, used_t = set(), set()
    for _, i, k in pairs:
        if i not in used_a and k not in used_t:
            used_a.add(i)
            used_t.add(k)
    tp, n_t, n_a = len(used_t), len(t_js), len(a_js)
    return {"n_truth": n_t, "n_answer": n_a, "tp": tp,
            "fnr": (n_t - tp) / n_t if n_t else None,
            "fdr": (n_a - tp) / n_a if n_a else None}


def candidate_fpr(truth: Sequence[Cycle], answer: Sequence[Cycle], graph: aaf.AAGraph, tol: int = ADJ_TOL) -> Dict[str, object]:
    """Spurious candidates (graph edges matching no truth adjacency) used by the answer / all of them."""
    t_js = novel_adjacencies(truth, phased=False)
    a_js = novel_adjacencies(answer, phased=False)

    def hit(j, pool):
        return any((d := aa_eval.junction_distance(j, k, False)) is not None and d <= tol for k in pool)

    spurious = [j for j in aa_eval._distinct(aa_eval.graph_adjacencies(graph), False) if not hit(j, t_js)]
    used = sum(1 for j in spurious if hit(j, a_js))
    return {"n_candidates": len(aa_eval._distinct(aa_eval.graph_adjacencies(graph), False)), "n_spurious": len(spurious),
            "n_spurious_used": used, "cand_fpr": used / len(spurious) if spurious else None}


# ---------------------------------------------------------------------------------------------- #
# copy-number distance (Eq. 11)
# ---------------------------------------------------------------------------------------------- #

def cn_distance(table: pd.DataFrame) -> Tuple[float, float]:
    """(sum_x (|c_A - t_A| + |c_B - t_B| + c_U) x bp, sum_x t x bp) of an allele table."""
    if table.empty:
        return 0.0, 0.0
    w = (table["end"] - table["start"] + 1).astype(float)
    dist = ((table["c_A"] - table["t_A"]).abs() + (table["c_B"] - table["t_B"]).abs() + table["c_U"]) * w
    return float(dist.sum()), float(((table["t_A"] + table["t_B"]) * w).sum())


# ---------------------------------------------------------------------------------------------- #
# generalized infinite-sites violations
# ---------------------------------------------------------------------------------------------- #

def _extremities(junctions: Sequence[Junction], tol: int) -> Dict[Tuple[str, int, str], int]:
    """breakend -> extremity id; breakends of one (chrom, strand) within ``tol`` (single linkage) share one."""
    by: Dict[Tuple[str, str], List[int]] = defaultdict(list)
    for a, b, _, _ in junctions:
        for be in (a, b):
            by[(be[0], be[2])].append(int(be[1]))
    ids: Dict[Tuple[str, int, str], int] = {}
    next_id = 0
    for (chrom, strand), positions in sorted(by.items()):
        last = None
        for pos in sorted(set(positions)):
            if last is None or pos - last > tol:
                next_id += 1
            ids[(chrom, pos, strand)] = next_id
            last = pos
    return ids


def is_violations(junctions: Sequence[Junction], tol: int = EXT_TOL) -> Dict[str, int]:
    """The module docstring's is_* counts over a set of (labelled) novel adjacencies."""
    ids = _extremities(junctions, tol)
    adjs_of: Dict[int, set] = defaultdict(set)      # extremity -> unlabelled adjacencies
    labels_of: Dict[frozenset, set] = defaultdict(set)  # unlabelled adjacency -> labellings
    haps_of: Dict[int, set] = defaultdict(set)      # extremity -> haplotypes it is used on
    for a, b, h1, h2 in junctions:
        x, y = ids[a], ids[b]
        key = frozenset((x, y))
        adjs_of[x].add(key)
        adjs_of[y].add(key)
        if h1 in HAPS and h2 in HAPS:
            labels_of[key].add(frozenset({(x, h1), (y, h2)}))
        for e, h in ((x, h1), (y, h2)):
            if h in HAPS:
                haps_of[e].add(h)
    # reciprocal pairs: the "+" extremity ending at p and the "-" extremity starting at p + 1
    plus: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    minus: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for (chrom, pos, strand), e in ids.items():
        (plus if strand == "+" else minus)[chrom].append((pos, e))
    recip = set()
    for chrom, ends in plus.items():
        for p, x in ends:
            for q, y in minus.get(chrom, []):
                if abs(q - (p + 1)) <= tol and any(hx != hy for hx in haps_of[x] for hy in haps_of[y]):
                    recip.add((x, y))
    return {"n_adjacencies": len({k for s in adjs_of.values() for k in s}),
            "n_extremities": len(set(ids.values())),
            "is_ext": sum(1 for s in adjs_of.values() if len(s) > 1),
            "is_hom_same": sum(1 for s in labels_of.values() if len(s) > 1),
            "is_hom": sum(1 for s in haps_of.values() if len(s) > 1),
            "is_recip": len(recip)}


# ---------------------------------------------------------------------------------------------- #
# per-clone copy numbers (Eq. 11 over clones; copy-number states)
# ---------------------------------------------------------------------------------------------- #

def _intervals(path: str):
    """[(chrom, start, end, {clone: {A, B}})] of an RCK segments file (1-based inclusive, cn= in extra)."""
    from rck.ecdna.rck_inputs import read_segments
    return [(s.chrom, s.start, s.end, s.extra.get("cn", {})) for s in read_segments(path)]


def _state(x: float) -> int:
    return 0 if x < 0.5 else 1 if x < 1.5 else 2


def ecdna_regions(truth: Sequence[Cycle]) -> List[Tuple[str, int, int]]:
    """The truth's ecDNA bases as merged 1-based inclusive intervals (walk positions are 0-based)."""
    by: Dict[str, List[Tuple[int, int]]] = defaultdict(list)
    for c in truth:
        for x in c.steps:
            by[x.chrom].append((min(x.start, x.end) + 1, max(x.start, x.end) + 1))
    out = []
    for chrom, ivs in by.items():
        ivs.sort()
        cur = list(ivs[0])
        for a, b in ivs[1:]:
            if a <= cur[1] + 1:
                cur[1] = max(cur[1], b)
            else:
                out.append((chrom, cur[0], cur[1]))
                cur = [a, b]
        out.append((chrom, cur[0], cur[1]))
    return out


def clone_cn_scores(inferred_scnt: str, truth_scnt: str, fragments_path: str,
                    regions: Sequence[Tuple[str, int, int]] = ()) -> Dict[str, object]:
    """The module docstring's per-clone Eq. 11 and copy-number state agreement (4)."""
    inferred, truth = _intervals(inferred_scnt), _intervals(truth_scnt)
    frags = [(f[0], f[1], f[2]) for f in _intervals(fragments_path)]
    clones = sorted({c for _, _, _, cn in inferred for c in cn})
    by_chrom: Dict[str, list] = defaultdict(list)
    for t in truth:
        by_chrom[t[0]].append(t)

    def frag_of(chrom, start, end):
        mid = (start + end) // 2
        return next((k for k, (c, a, b) in enumerate(frags) if c == chrom and a <= mid <= b), None)

    def in_regions(chrom, a, b):
        return sum(max(0, min(b, y) - max(a, x) + 1) for c, x, y in regions if c == chrom)

    # pieces: (fragment, bp, bp on ecDNA, inferred cn, truth cn)
    pieces = []
    for chrom, a, b, cn in inferred:
        for tc, ta, tb, tcn in by_chrom.get(chrom, []):
            lo, hi = max(a, ta), min(b, tb)
            if lo <= hi:
                pieces.append((frag_of(chrom, a, b), hi - lo + 1, in_regions(chrom, lo, hi), cn, tcn))

    def cost(cn, tcn, swap):
        total = 0.0
        for c in clones:
            x, y = (cn[c]["B"], cn[c]["A"]) if swap else (cn[c]["A"], cn[c]["B"])
            total += abs(x - tcn.get(c, {}).get("A", 0)) + abs(y - tcn.get(c, {}).get("B", 0))
        return total

    out: Dict[str, object] = {"clones": clones}
    for scope, w_index in (("all", 1), ("ecdna", 2)):
        keep: Dict[object, float] = defaultdict(float)
        swap: Dict[object, float] = defaultdict(float)
        mass = 0.0
        for p in pieces:
            w = p[w_index]
            if not w:
                continue
            keep[p[0]] += cost(p[3], p[4], False) * w
            swap[p[0]] += cost(p[3], p[4], True) * w
            mass += sum(p[4].get(c, {}).get("A", 0) + p[4].get(c, {}).get("B", 0) for c in clones) * w
        swapped = {f for f in keep if swap[f] < keep[f] - 1e-9}
        dist = sum(min(keep[f], swap[f]) for f in keep)
        agree = total = 0.0
        for f, bp, bp_ec, cn, tcn in pieces:
            w = bp if scope == "all" else bp_ec
            if not w:
                continue
            for c in clones:
                x, y = (cn[c]["B"], cn[c]["A"]) if f in swapped else (cn[c]["A"], cn[c]["B"])
                for got, want in ((x, tcn.get(c, {}).get("A", 0)), (y, tcn.get(c, {}).get("B", 0))):
                    agree += w * (_state(got) == _state(want))
                    total += w
        out.update({"dist_" + scope: dist, "mass_" + scope: mass, "dist_rel_" + scope: dist / mass if mass else None,
                    "state_agree_" + scope: agree / total if total else None})
    return out


# ---------------------------------------------------------------------------------------------- #
# input use and linear walks
# ---------------------------------------------------------------------------------------------- #

def rck_adjacencies(path: str) -> Dict[str, Junction]:
    """aid -> unlabelled junction (walk frame: RCK coordinate - 1, as AA graph positions) of an rck.adj.tsv."""
    from rck.ecdna.rck_inputs import read_adjacencies
    return {a.aid: ((a.chrom1, a.coord1 - 1, a.strand1), (a.chrom2, a.coord2 - 1, a.strand2), ".", ".")
            for a in read_adjacencies(path)}


def input_use(inputs: Sequence[Junction], used: Sequence[bool], truth: Sequence[Cycle], tol: int = ADJ_TOL) -> Dict[str, object]:
    """Share of the input novel adjacencies used, and RCK's FPR over the spurious ones (module docstring, 5)."""
    t_js = novel_adjacencies(truth, phased=False)
    spurious = [not any((d := aa_eval.junction_distance(j, k, False)) is not None and d <= tol for k in t_js) for j in inputs]
    n_sp = sum(spurious)
    n_sp_used = sum(1 for s, u in zip(spurious, used) if s and u)
    return {"n_input": len(inputs), "n_used": int(sum(used)), "input_used": sum(used) / len(inputs) if inputs else None,
            "n_spurious": n_sp, "n_spurious_used": n_sp_used, "input_fpr": n_sp_used / n_sp if n_sp else None}


def multiclone_input_use(run_dir: str, rck_input: str, truth: Sequence[Cycle], tol: int = ADJ_TOL) -> Dict[str, object]:
    with open(os.path.join(run_dir, "summary.json")) as handle:
        used_aids = {u["aid"] for u in json.load(handle).get("novel_used", []) if not u.get("imputed")}
    adj = rck_adjacencies(os.path.join(rck_input, "rck.adj.tsv"))
    aids = sorted(adj)
    return input_use([adj[a] for a in aids], [a in used_aids for a in aids], truth, tol)


def graph_input_use(graph: aaf.AAGraph, strict: Sequence[Cycle], truth: Sequence[Cycle], tol: int = ADJ_TOL) -> Dict[str, object]:
    inputs = aa_eval._distinct(aa_eval.graph_adjacencies(graph), False)
    a_js = novel_adjacencies(strict, phased=False)
    used = [any((d := aa_eval.junction_distance(j, k, False)) is not None and d <= tol for k in a_js) for j in inputs]
    return input_use(inputs, used, truth, tol)


def linear_walks(strict: Sequence[Cycle]) -> int:
    return sum(1 for c in strict if not c.circular)


# ---------------------------------------------------------------------------------------------- #
# one answer
# ---------------------------------------------------------------------------------------------- #

def score_truth(truth: Sequence[Cycle]) -> Dict[str, object]:
    """Violations of the truth itself (a property of the simulation)."""
    return {"is": is_violations(novel_adjacencies(truth, phased=True))}


def score(strict: Sequence[Cycle], truth: Sequence[Cycle], blocks: Optional[Sequence[ev.Block]] = None,
          graph: Optional[aaf.AAGraph] = None, tol: int = ADJ_TOL) -> Dict[str, object]:
    """Every score of the module docstring for one answer (its strict cycles, weights already set)."""
    table = am.allele_table(truth, strict)
    swapped = am.swapped_blocks(table, blocks) if blocks else set()
    t_js = novel_adjacencies(truth, phased=True)
    a_js = novel_adjacencies(strict, phased=True)
    a_rel = relabel_junctions(a_js, blocks, swapped) if blocks else a_js
    within, _ = am.within_block(table, blocks) if blocks else (table, 0)
    dist_abs, mass = cn_distance(table)
    dist_block, _ = cn_distance(within)
    out: Dict[str, object] = {
        "n_strict": len(strict), "n_blocks_swapped": len(swapped),
        "cn_dist_abs": dist_abs, "cn_dist": dist_block, "truth_mass_bp": mass,
        "cn_dist_abs_rel": dist_abs / mass if mass else None, "cn_dist_rel": dist_block / mass if mass else None,
        "adj_blind": adjacency_rates(truth, strict, False, tol),
        "adj_phased_abs": junction_rates(t_js, a_js, True, tol),
        "adj_phased": junction_rates(_keep(t_js, blocks), _keep(a_rel, blocks), True, tol),
        "is_abs": is_violations(a_js),
        "is": is_violations(_keep(a_rel, blocks)),
    }
    if graph is not None:
        out["candidates"] = candidate_fpr(truth, strict, graph, tol)
    out["linear_walks"] = linear_walks(strict)
    return out


def flatten(result: Dict[str, object]) -> Dict[str, object]:
    """One-level dict: nested blocks become prefix_key columns."""
    row: Dict[str, object] = {}
    for key, value in result.items():
        if isinstance(value, dict):
            row.update({"{p}_{k}".format(p=key, k=k): v for k, v in value.items()})
        else:
            row[key] = value
    return row


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cycles", default=None)
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--truth-only", action="store_true", help="only the truth's generalized-IS violations")
    p.add_argument("--graph", default=None, help="candidate graph (AA format) for cand_fpr")
    p.add_argument("--phase-blocks", default=None, help="rck_input/rck.fragments.tsv: scores up to one A/B swap per block")
    p.add_argument("--weights-from-multiclone", default=None, help="multi-clone run: w_c = sum_k u_k psi_k")
    p.add_argument("--clone-proportions", default=None)
    a = p.parse_args(argv)
    truth = ev.load_truth(a.truth_dir)
    if a.truth_only:
        result = score_truth(truth)
    else:
        if not a.cycles:
            p.error("--cycles is required unless --truth-only")
        weights = None
        if a.weights_from_multiclone:
            weights = bm.multiclone_bulk_weights(a.weights_from_multiclone, bm.read_proportions(a.clone_proportions))
        strict, _ = am.strict_cycles(a.cycles, weights)
        blocks = ev.read_phase_blocks(a.phase_blocks) if a.phase_blocks else None
        graph = aaf.read_graph(a.graph) if a.graph else None
        result = score(strict, truth, blocks, graph)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "rck_eval.json"), "w") as handle:
        json.dump(result, handle, indent=2)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

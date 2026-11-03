"""Method-agnostic ecDNA scorer: AA-format output of any method against simulator truth.

    python -m rck.ecdna.aa_eval --cycles METHOD_cycles.txt --graph METHOD_graph.txt \\
        --truth-dir SIM_RUN [--alt-cycles ALT_cycles.txt ...] [--segment-cn cn.tsv] \\
        [--paths SUBGRAPH/paths.tsv] [--strict-rule auto] --outdir OUT

Inputs are AmpliconArchitect files only (:mod:`rck.ecdna.aa_format`), so AmpliconArchitect,
CoRAL, the original RCK (once converted) and ``rck-ecdna`` all go through the same code:

* ``--cycles``: the method's cycles file. Its ``Interval`` lines (or, without them, the graph's
  segments) define the region the method looked at, and so which truth species are scored.
  ``Haplotypes=`` makes the cycles phased; without it they are scored haplotype-blind;
* ``--alt-cycles``: other answers of the method, used by the best-of mode only (for
  ``rck-ecdna``, rank 1's first k same-fit alternatives);
* ``--graph``: the method's breakpoint graph, for the graph-level metrics;
* ``--segment-cn``: optional allele copy numbers per segment (``chrom start end cn_a cn_b``, AA
  coordinates); without it the copy-number error is on total copy number;
* ``--truth-dir``: a simulator run. Species truth is read as ``rck-ecdna-eval`` reads it
  (:func:`rck.ecdna.evaluate.truth_from_simulation`); the graph truth is the mixture graph of
  amplicon-allele-analysis ``simulations/truth_mixture.py`` (``truth/aa/mixture_amplicon<N>_*``,
  or ``--mixture-dir``).

Metrics (definitions and sources in docs/ecDNA.md, section 6, "Method-agnostic scorer"):

graph level, per truth amplicon, at ``--graph-tol`` (100 bp, the AA and CoRAL convention)
    ``error_free`` (AA) and ``error_free_merged`` (the same after merging benign copy-number
    splits), ``bp_accuracy`` (CoRAL), ``novel_adjacency_fnr`` / ``_fdr`` (RCK),
    ``segment_cn_mae`` and ``segment_cn_rel_err`` (allele-specific when ``--segment-cn`` is
    given), ``is_violations`` (extremities used by more than one novel adjacency);
species level, after one-to-one matching
    the ``rck-ecdna-eval`` scores at ``--species-tol`` (10 kb), unchanged; CoRAL's cycle metrics,
    headline as CoRAL's code computes them (``cycle_overlap_coral_code``, ``cyclic_lcs_coral_code``),
    the supplementary-methods versions as ``cycle_overlap_supp`` and ``cyclic_lcs_supp``, and
    ``length_error``; ``exact_cycle``; ``composite_h0``;
sample level
    species-count error, missed and false species, ``cn_explained_top1`` / ``_top3``
    (CoRAL), ``subwalks_satisfied``.

Two modes, both always reported: ``strict`` (the method's own answer, after ``--strict-rule``;
counts false species) and ``best_of`` (per truth species, the best of every cycle the method
returned, ``--cycles`` and ``--alt-cycles``, each metric maximised on its own as CoRAL does).
"""

from __future__ import annotations

import argparse
import glob
import json
import logging
import math
import os
from dataclasses import asdict, dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from rck.ecdna import aa_format as aaf
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import CStep, Cycle

log = logging.getLogger(__name__)

GRAPH_TOL = 100
SPECIES_TOL = 10_000
STRICT_CN_FRAC = 0.1  # rck-ecdna's own --min-species-weight-frac
CN_EXPLAINED_MIN_CN = 5.0  # CoRAL: sequence edges below this copy number are left out of the total

Breakend = aaf.Breakend


# ---------------------------------------------------------------------------------------------- #
# graph level
# ---------------------------------------------------------------------------------------------- #

def _be_close(x: Breakend, y: Breakend, tol: int) -> bool:
    return x[0] == y[0] and x[2] == y[2] and abs(x[1] - y[1]) <= tol


def edges_match(e: aaf.AAEdge, f: aaf.AAEdge, tol: int) -> bool:
    """Same orientations and both breakpoints within ``tol``, in either direction (CoRAL's
    ``is_breakpoint_match``, with <= instead of <)."""
    return (_be_close(e.a, f.a, tol) and _be_close(e.b, f.b, tol)) or (_be_close(e.a, f.b, tol) and _be_close(e.b, f.a, tol))


def _in_intervals(intervals, chrom, pos, pad=0) -> bool:
    return any(c == chrom and s - pad <= pos <= e + pad for c, s, e in intervals)


def _overlap_len(intervals, chrom, s, e) -> int:
    return sum(max(0, min(e, ie) - max(s, ib) + 1) for c, ib, ie in intervals if c == chrom)


def scoped_edges(graph: aaf.AAGraph, intervals, tol: int) -> List[aaf.AAEdge]:
    """Discordant edges with at least one breakend in ``intervals`` (padded by ``tol``)."""
    return [e for e in graph.discordant if _in_intervals(intervals, e.a[0], e.a[1], tol) or _in_intervals(intervals, e.b[0], e.b[1], tol)]


def scoped_sequences(graph: aaf.AAGraph, intervals, tol: int) -> List[aaf.AASequence]:
    """Sequence edges overlapping ``intervals`` by more than ``tol`` bp (a flank segment that
    only touches an amplicon boundary is not part of it)."""
    return [s for s in graph.sequences if _overlap_len(intervals, s.chrom, s.start, s.end) > tol]


def is_violations(edges: Sequence[aaf.AAEdge], tol: int) -> int:
    """Extremities used by more than one distinct novel adjacency: breakends of the same
    chromosome and strand within ``tol`` of each other form one extremity (single linkage)."""
    ends: Dict[Tuple[str, str], List[Tuple[int, int]]] = {}
    distinct: List[aaf.AAEdge] = []
    for e in edges:
        if not any((e.a, e.b) in ((f.a, f.b), (f.b, f.a)) for f in distinct):
            distinct.append(e)
    for k, e in enumerate(distinct):
        for b in {e.a, e.b}:
            ends.setdefault((b[0], b[2]), []).append((b[1], k))
    violations = 0
    for points in ends.values():
        points.sort()
        cluster: set = set()
        last = None
        for pos, k in points:
            if last is not None and pos - last > tol:
                violations += len(cluster) > 1
                cluster = set()
            cluster.add(k)
            last = pos
        violations += len(cluster) > 1
    return violations


def _piecewise(table: Sequence[Tuple[str, int, int, Tuple[float, ...]]], chrom, s, e, width):
    """[(length, values)] of ``table`` over [s, e]; bases no row covers get zeros."""
    rows = sorted((rs, re_, v) for c, rs, re_, v in table if c == chrom and rs <= e and re_ >= s)
    out, pos = [], s
    for rs, re_, v in rows:
        lo, hi = max(rs, pos), min(re_, e)
        if lo > pos:
            out.append((lo - pos, (0.0,) * width))
        if hi >= lo:
            out.append((hi - lo + 1, v))
            pos = hi + 1
    if pos <= e:
        out.append((e - pos + 1, (0.0,) * width))
    return out


def segment_cn_error(truth_segments: pd.DataFrame, inferred_graph: aaf.AAGraph, inferred_cn: Optional[pd.DataFrame]) -> Dict[str, object]:
    """Length-weighted copy-number error over the truth amplicon's segments.

    With allele copy numbers on both sides the error per base is |ΔA| + |ΔB| (``mode`` allele);
    otherwise |Δ total|. ``segment_cn_mae`` is the mean over bases, ``segment_cn_rel_err`` the
    summed error over the summed truth copy number (both length-weighted); bases the method
    gives no segment count as copy number 0."""
    allele = inferred_cn is not None and {"cn_a", "cn_b"} <= set(inferred_cn.columns) and truth_segments["cn_a"].notna().all()
    if inferred_cn is not None:
        total_rows = [(str(r["chrom"]), int(r["start"]), int(r["end"]),
                       (float(r["cn"]) if "cn" in r and pd.notna(r.get("cn")) else float(r["cn_a"]) + float(r["cn_b"]),))
                      for r in inferred_cn.to_dict("records")]
    else:
        total_rows = [(s.chrom, s.start, s.end, (s.cn,)) for s in inferred_graph.sequences]
    allele_rows = [(str(r["chrom"]), int(r["start"]), int(r["end"]), (float(r["cn_a"]), float(r["cn_b"])))
                   for r in inferred_cn.to_dict("records")] if allele else []
    length = err_total = err_allele = truth_total = 0.0
    for r in truth_segments.to_dict("records"):
        c, s, e = str(r["chrom"]), int(r["start"]), int(r["end"])
        length += e - s + 1
        truth_total += float(r["cn"]) * (e - s + 1)
        err_total += sum(n * abs(v[0] - float(r["cn"])) for n, v in _piecewise(total_rows, c, s, e, 1))
        if allele:
            err_allele += sum(n * (abs(v[0] - float(r["cn_a"])) + abs(v[1] - float(r["cn_b"])))
                              for n, v in _piecewise(allele_rows, c, s, e, 2))
    err = err_allele if allele else err_total
    return {"cn_error_mode": "allele" if allele else "total",
            "segment_cn_mae": err / length if length else None,
            "segment_cn_rel_err": err / truth_total if truth_total else None,
            "segment_cn_mae_total": err_total / length if length else None}


def merge_benign_splits(graph: aaf.AAGraph, tol: int = 1) -> List[aaf.AASequence]:
    """The graph's sequence edges with every benign split removed: neighbouring segments
    (``b.start == a.end + 1``) joined by a concordant edge ``a.end+ -> b.start-``, with no
    discordant or source edge at either side of that boundary, become one segment. Such a
    boundary is a copy-number change (or an imputed grid point) without a breakpoint, which AA
    and CoRAL graphs also write. Breakends match within ``tol`` bp. The merged copy number is the
    length-weighted mean (unused by ``error_free_merged``)."""
    concordant = [e for e in graph.edges if e.kind == "concordant"]
    breaks = [b for e in graph.edges if e.kind in ("discordant", "source") for b in (e.a, e.b)]

    def joined(a: aaf.AASequence, b: aaf.AASequence) -> bool:
        if a.chrom != b.chrom or b.start != a.end + 1:
            return False
        left, right = (a.chrom, a.end, "+"), (b.chrom, b.start, "-")
        if not any(edges_match(e, aaf.AAEdge("concordant", left, right, 0.0), tol) for e in concordant):
            return False
        return not any(_be_close(x, left, tol) or _be_close(x, right, tol) for x in breaks)

    out: List[aaf.AASequence] = []
    for s in sorted(graph.sequences, key=lambda x: (x.chrom, x.start, x.end)):
        if out and joined(out[-1], s):
            a = out[-1]
            la, ls = a.end - a.start + 1, s.end - s.start + 1
            out[-1] = aaf.AASequence(a.chrom, a.start, s.end, (a.cn * la + s.cn * ls) / (la + ls))
        else:
            out.append(aaf.AASequence(s.chrom, s.start, s.end, s.cn))
    return out


def _segment_match(true_segs, inf_segs, tol: int) -> Tuple[List[bool], list]:
    seg_match = lambda a, b: a.chrom == b.chrom and abs(a.start - b.start) <= tol and abs(a.end - b.end) <= tol  # noqa: E731
    found = [any(seg_match(t, i) for i in inf_segs) for t in true_segs]
    false = [i for i in inf_segs if not any(seg_match(i, t) for t in true_segs)]
    return found, false


@dataclass
class MixtureAmplicon:
    index: int
    species: List[str]
    intervals: List[Tuple[str, int, int]]  # AA, 1-based inclusive
    graph: aaf.AAGraph
    segments: pd.DataFrame


def load_mixture(mixture_dir: str) -> List[MixtureAmplicon]:
    """The truth mixture amplicons written by amplicon-allele-analysis ``simulations/truth_mixture.py``."""
    index_path = os.path.join(mixture_dir, "mixture_amplicons.tsv")
    if not os.path.exists(index_path):
        raise FileNotFoundError("{p} not found: build it with amplicon-allele-analysis "
                                "simulations/truth_mixture.py".format(p=index_path))
    out = []
    for r in pd.read_csv(index_path, sep="\t").to_dict("records"):
        cycles = aaf.read_cycles(os.path.join(mixture_dir, r["cycles"]))
        out.append(MixtureAmplicon(int(r["amplicon"]), str(r["species"]).split(","), cycles.intervals,
                                   aaf.read_graph(os.path.join(mixture_dir, r["graph"])),
                                   pd.read_csv(os.path.join(mixture_dir, r["segments"]), sep="\t")))
    return out


def graph_metrics(amp: MixtureAmplicon, inferred: aaf.AAGraph, inferred_cn: Optional[pd.DataFrame], tol: int = GRAPH_TOL) -> Dict[str, object]:
    """Graph-level metrics of one truth amplicon against the method's graph (module docstring)."""
    true_edges = amp.graph.discordant
    inf_edges = scoped_edges(inferred, amp.intervals, tol)
    true_segs = amp.graph.sequences
    inf_segs = scoped_sequences(inferred, amp.intervals, tol)
    found = [any(edges_match(t, i, tol) for i in inf_edges) for t in true_edges]
    false_edges = [i for i in inf_edges if not any(edges_match(i, t, tol) for t in true_edges)]
    segs_found, false_segs = _segment_match(true_segs, inf_segs, tol)
    # error_free_merged: the same rule after removing benign splits (merge_benign_splits) in both graphs
    true_merged = [x for x in merge_benign_splits(amp.graph) if _overlap_len(amp.intervals, x.chrom, x.start, x.end) > tol]
    inf_merged = [x for x in merge_benign_splits(inferred) if _overlap_len(amp.intervals, x.chrom, x.start, x.end) > tol]
    msegs_found, mfalse_segs = _segment_match(true_merged, inf_merged, tol)
    edges_ok = all(found) and not false_edges
    n_true, n_inf = len(true_edges), len(inf_edges)
    recall = sum(found) / n_true if n_true else None
    row = {"amplicon": amp.index, "species": ",".join(amp.species),
           "intervals": ",".join("{c}:{s}-{e}".format(c=c, s=s, e=e) for c, s, e in amp.intervals), "graph_tol_bp": tol,
           "n_true_segments": len(true_segs), "n_segments_found": sum(segs_found), "n_inferred_segments": len(inf_segs),
           "n_false_segments": len(false_segs), "n_true_edges": n_true, "n_edges_found": sum(found),
           "n_inferred_edges": n_inf, "n_false_edges": len(false_edges),
           "error_free": int(all(segs_found) and not false_segs and edges_ok),
           "n_true_segments_merged": len(true_merged), "n_segments_found_merged": sum(msegs_found),
           "n_inferred_segments_merged": len(inf_merged), "n_false_segments_merged": len(mfalse_segs),
           "error_free_merged": int(all(msegs_found) and not mfalse_segs and edges_ok),
           "bp_accuracy": recall,
           "novel_adjacency_fnr": None if recall is None else 1.0 - recall,
           "novel_adjacency_fdr": len(false_edges) / n_inf if n_inf else None,
           "is_violations": is_violations(inf_edges, tol), "is_violations_truth": is_violations(true_edges, tol)}
    row.update(segment_cn_error(amp.segments, inferred, inferred_cn))
    return row


# ---------------------------------------------------------------------------------------------- #
# species level: CoRAL's cycle metrics
# ---------------------------------------------------------------------------------------------- #

def strip_haplotypes(c: Cycle) -> Cycle:
    """The cycle with every haplotype set to ``.`` (haplotype-blind scoring)."""
    orderings = [[CStep(s.chrom, s.start, s.end, aaf.UNPHASED, s.strand) for s in o] for o in c.orderings]
    return Cycle(name=c.name, weight=c.weight, orderings=orderings, circular=c.circular, weight_int=c.weight_int)


def cycle_length(steps: Sequence[CStep]) -> int:
    return sum(s.end - s.start + 1 for s in steps)


def _merged(steps: Sequence[CStep]) -> Dict[str, List[Tuple[int, int]]]:
    """Per chromosome, the union of the steps' intervals, half-open, without multiplicity."""
    out: Dict[str, List[Tuple[int, int]]] = {}
    for chrom in sorted({s.chrom for s in steps}):
        ivs = sorted((s.start, s.end + 1) for s in steps if s.chrom == chrom)
        merged = [list(ivs[0])]
        for a, b in ivs[1:]:
            if a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        out[chrom] = [tuple(x) for x in merged]
    return out


def cycle_overlap(truth: Sequence[CStep], inferred: Sequence[CStep]) -> float:
    """CoRAL cycle-interval overlap: bases in both interval sets over bases in either, without
    multiplicity; intervals on different chromosomes never overlap."""
    mt, mi = _merged(truth), _merged(inferred)
    inter = 0
    for chrom in set(mt) & set(mi):
        for a, b in mt[chrom]:
            for c, d in mi[chrom]:
                inter += max(0, min(b, d) - max(a, c))
    total = sum(b - a for ivs in mt.values() for a, b in ivs) + sum(b - a for ivs in mi.values() for a, b in ivs)
    union = total - inter
    return inter / union if union else 0.0


def _bins(truth: Sequence[CStep], inferred: Sequence[CStep]) -> Dict[str, List[int]]:
    points: Dict[str, set] = {}
    for s in list(truth) + list(inferred):
        points.setdefault(s.chrom, set()).update((s.start, s.end + 1))
    return {c: sorted(p) for c, p in points.items()}


def _counts(steps: Sequence[CStep], bins: Dict[str, List[int]]) -> Dict[Tuple[str, int], int]:
    """{(chrom, bin index): how many steps cover the bin}."""
    out: Dict[Tuple[str, int], int] = {}
    for s in steps:
        edges = bins[s.chrom]
        lo, hi = edges.index(s.start), edges.index(s.end + 1)
        for k in range(lo, hi):
            out[(s.chrom, k)] = out.get((s.chrom, k), 0) + 1
    return out


def cycle_overlap_coral_code(truth: Sequence[CStep], inferred: Sequence[CStep]) -> float:
    """CoRAL's code (``scoring_utils.get_fragments_similarity_unweighted``): 1 - sum over bins of
    |e1 - e2| x length over sum of max(e1, e2) x length, e1 and e2 the true and inferred
    multiplicities. Equal to :func:`cycle_overlap` when no interval is used twice."""
    bins = _bins(truth, inferred)
    ct, ci = _counts(truth, bins), _counts(inferred, bins)
    diff = top = 0
    for key in set(ct) | set(ci):
        n = bins[key[0]][key[1] + 1] - bins[key[0]][key[1]]
        a, b = ct.get(key, 0), ci.get(key, 0)
        diff += abs(a - b) * n
        top += max(a, b) * n
    return 1.0 - diff / top if top else 0.0


def _labels(steps: Sequence[CStep], bins, keep) -> List[Tuple[Tuple[str, int], str]]:
    out = []
    for s in steps:
        edges = bins[s.chrom]
        ks = range(edges.index(s.start), edges.index(s.end + 1))
        for k in (reversed(ks) if s.strand == "-" else ks):
            if (s.chrom, k) in keep:
                out.append(((s.chrom, k), s.strand))
    return out


def _weighted_lcs(a, b, weight) -> int:
    prev = [0] * (len(b) + 1)
    for x in a:
        cur = [0] * (len(b) + 1)
        wx = weight[x[0]]
        for j, y in enumerate(b, start=1):
            cur[j] = prev[j - 1] + wx if x == y else (prev[j] if prev[j] > cur[j - 1] else cur[j - 1])
        prev = cur
    return prev[-1]


def cyclic_lcs(truth: Sequence[CStep], inferred: Sequence[CStep], inferred_circular: bool = True) -> Tuple[float, float, int]:
    """CoRAL's cyclic longest common subsequence: (normalised by the true cycle's length,
    normalised as CoRAL's code does, raw bp).

    Both cycles are cut into the bins of their joint breakpoints and written as sequences of
    (bin, orientation), keeping only bins both cover. The LCS, each match weighted by the bin's
    length, is maximised over every rotation of the inferred cycle and of its reverse (one
    rotation for a linear walk). ``cyclic_lcs`` divides it by the true cycle's length (with
    multiplicity, CoRAL supplementary methods); CoRAL's code
    (``compute_normalized_longest_cycle_subsequence``) divides by the sum over bins of
    max(e1, e2) x length instead."""
    bins = _bins(truth, inferred)
    ct, ci = _counts(truth, bins), _counts(inferred, bins)
    keep = set(ct) & set(ci)
    weight = {key: bins[key[0]][key[1] + 1] - bins[key[0]][key[1]] for key in keep}
    top = sum(max(ct.get(k, 0), ci.get(k, 0)) * (bins[k[0]][k[1] + 1] - bins[k[0]][k[1]]) for k in set(ct) | set(ci))
    if not keep:
        return 0.0, 0.0, 0
    t = _labels(truth, bins, keep)
    fw = _labels(inferred, bins, keep)
    rev = [(k, "-" if o == "+" else "+") for k, o in reversed(fw)]
    best = 0
    for seq in (fw, rev):
        for r in (range(len(seq)) if inferred_circular else (0,)):
            best = max(best, _weighted_lcs(t, seq[r:] + seq[:r], weight))
    return best / cycle_length(truth), (best / top if top else 0.0), best


def length_error(truth: Sequence[CStep], inferred: Sequence[CStep]) -> float:
    """CoRAL's length error log2(1 + R), R = (L_inferred - L_true) / L_true; that is
    log2(L_inferred / L_true), 0 for the right length."""
    return math.log2(cycle_length(inferred) / cycle_length(truth))


# ---------------------------------------------------------------------------------------------- #
# species level: pairs and modes
# ---------------------------------------------------------------------------------------------- #

PHASED_FIELDS = ["haplotype_jaccard", "phased_junction_recall", "phased_junction_precision", "phased_junction_f1",
                 "phasing_accuracy", "phasing_accuracy_swap_invariant", "composite", "exact_cycle_phased"]
SPECIES_COLUMNS = ["mode", "truth", "truth_weight", "truth_length_bp", "inferred", "inferred_weight", "inferred_topology",
                   "segment_recall", "segment_precision", "segment_jaccard", "haplotype_jaccard",
                   "junction_recall", "junction_precision", "junction_f1",
                   "phased_junction_recall", "phased_junction_precision", "phased_junction_f1",
                   "phasing_accuracy", "phasing_accuracy_swap_invariant", "order_match", "topology_match", "cn_rel_err",
                   "composite", "exact_cycle_phased",
                   "inferred_h0", "segment_jaccard_h0", "junction_f1_h0", "order_match_h0", "cn_rel_err_h0", "composite_h0",
                   "cycle_overlap_coral_code", "cyclic_lcs_coral_code", "cycle_overlap_supp", "cyclic_lcs_supp", "length_error",
                   "exact_cycle"]
# best-of: each metric on its own over the candidates; these are minimised, the rest maximised
MINIMISED = {"cn_rel_err", "cn_rel_err_h0"}
ABS_MINIMISED = {"length_error"}


def composite_h0(s: ev.PairScore) -> float:
    """The composite without phasing accuracy, on haplotype-blind cycles: the mean of segment
    Jaccard, junction F1, order match and max(0, 1 - copy-number error)."""
    if s.inferred is None:
        return 0.0
    cn = max(0.0, 1.0 - s.cn_rel_err) if s.cn_rel_err is not None else 0.0
    return (s.jaccard_total + s.junction_f1 + float(s.order_match) + cn) / 4


def _exact(t: Cycle, c: Cycle, tol: int) -> int:
    return int(any(ev.order_match(t.steps, o, tol, t.circular, c.circular) for o in c.orderings))


def phased_part(s: Optional[ev.PairScore], t: Optional[Cycle] = None, c: Optional[Cycle] = None, graph_tol: int = GRAPH_TOL) -> Dict[str, object]:
    """The rck-ecdna-eval scores of a phased pair, under the recovery-matrix names."""
    if s is None:
        return {}
    row = {ev._RENAMED.get(k, k): v for k, v in asdict(s).items()}
    row["exact_cycle_phased"] = _exact(t, c, graph_tol) if (t is not None and c is not None and s.inferred) else 0
    return row


def blind_part(s: ev.PairScore, t: Cycle, c: Optional[Cycle], graph_tol: int = GRAPH_TOL) -> Dict[str, object]:
    """Haplotype-blind scores of a pair of haplotype-stripped cycles (``c`` None: missed)."""
    row = {"inferred_h0": s.inferred, "segment_jaccard_h0": s.jaccard_total, "junction_f1_h0": s.junction_f1,
           "order_match_h0": s.order_match, "cn_rel_err_h0": s.cn_rel_err, "composite_h0": composite_h0(s)}
    if c is None:
        row.update({"cycle_overlap_coral_code": 0.0, "cyclic_lcs_coral_code": 0.0, "cycle_overlap_supp": 0.0, "cyclic_lcs_supp": 0.0,
                    "length_error": None, "exact_cycle": 0})
        return row
    lcs, lcs_code, _ = cyclic_lcs(ev.runs(t.steps), ev.runs(c.steps, circular=c.circular), c.circular)
    row.update({"cycle_overlap_coral_code": cycle_overlap_coral_code(t.steps, c.steps), "cyclic_lcs_coral_code": lcs_code,
                "cycle_overlap_supp": cycle_overlap(t.steps, c.steps), "cyclic_lcs_supp": lcs, "length_error": length_error(t.steps, c.steps),
                "exact_cycle": _exact(t, c, graph_tol)})
    return row


def _null_phased() -> Dict[str, object]:
    return {k: None for k in PHASED_FIELDS}


def strict_rows(truth: Sequence[Cycle], strict: Sequence[Cycle], phased: bool, species_tol: int, graph_tol: int) -> Tuple[List[dict], Dict[str, object]]:
    """Strict mode: one-to-one matching as rck-ecdna-eval does (by haplotype Jaccard when phased),
    and separately on haplotype-stripped cycles for the haplotype-blind scores. Returns the
    rows and the counts (missed and false species, from the phased matching when phased)."""
    t0 = [strip_haplotypes(t) for t in truth]
    c0 = [strip_haplotypes(c) for c in strict]
    blind = ev.score_decomposition(t0, c0, species_tol)
    by0 = {c.name: c for c in c0}
    phased_scores = ev.score_decomposition(truth, strict, species_tol) if phased else None
    by = {c.name: c for c in strict}
    rows = []
    for k, t in enumerate(truth):
        row = {"mode": "strict", "truth": t.name, "truth_weight": t.weight, "truth_length_bp": cycle_length(t.steps)}
        b = blind[k]
        if phased:
            p = phased_scores[k]
            row.update(phased_part(p, t, by.get(p.inferred), graph_tol))
        else:
            row.update(phased_part(b, t0[k], by0.get(b.inferred), graph_tol))
            row.update(_null_phased())
        row.update(blind_part(b, t0[k], by0.get(b.inferred), graph_tol))
        rows.append(row)
    matched = {r["inferred"] if phased else r["inferred_h0"] for r in rows} - {None}
    counts = {"n_missed": sum(1 for r in rows if (r["inferred"] if phased else r["inferred_h0"]) is None),
              "n_false": sum(1 for c in strict if c.name not in matched)}
    return rows, counts


def best_of_rows(truth: Sequence[Cycle], candidates: Sequence[Cycle], phased: bool, species_tol: int, graph_tol: int) -> List[dict]:
    """Best-of mode: per truth species, every candidate cycle that overlaps it (haplotype Jaccard
    > 0 when phased, segment Jaccard otherwise, as the strict matching requires); each metric
    takes its best value over them on its own (CoRAL). ``inferred`` / ``inferred_h0`` name the
    candidate with the best composite / composite_h0."""
    rows = []
    c0 = [strip_haplotypes(c) for c in candidates]
    for t in truth:
        t0 = strip_haplotypes(t)
        scored = []
        for c, cb in zip(candidates, c0):
            part = {}
            if phased:
                if ev.jaccard(t.steps, c.steps, True) <= 0:
                    p = None
                else:
                    p = ev.score_pair(t, c, species_tol)
                    part.update(phased_part(p, t, c, graph_tol))
            if ev.jaccard(t0.steps, cb.steps, True) > 0:
                b = ev.score_pair(t0, cb, species_tol)
                part.update(blind_part(b, t0, cb, graph_tol))
                if not phased:
                    part.update(phased_part(b, t0, cb, graph_tol))
                    part.update(_null_phased())
            if part:
                scored.append(part)
        row = {"mode": "best_of", "truth": t.name, "truth_weight": t.weight, "truth_length_bp": cycle_length(t.steps)}
        if not scored:
            row.update({"composite": 0.0 if phased else None, "composite_h0": 0.0, "exact_cycle": 0, "cycle_overlap_coral_code": 0.0,
                        "cyclic_lcs_coral_code": 0.0, "cycle_overlap_supp": 0.0, "cyclic_lcs_supp": 0.0})
            if phased:
                row["exact_cycle_phased"] = 0
            rows.append(row)
            continue
        for col in SPECIES_COLUMNS[4:]:
            if col in ("inferred", "inferred_h0", "inferred_weight", "inferred_topology"):
                continue
            vals = [s[col] for s in scored if s.get(col) is not None and not (isinstance(s[col], float) and math.isnan(s[col]))]
            if not vals:
                continue
            if col in MINIMISED:
                row[col] = min(vals)
            elif col in ABS_MINIMISED:
                row[col] = min(vals, key=abs)
            else:
                row[col] = max(vals)
        if phased:
            top = max((s for s in scored if s.get("composite") is not None), key=lambda s: s["composite"], default=None)
            if top is not None:
                row.update({"inferred": top["inferred"], "inferred_weight": top["inferred_weight"],
                            "inferred_topology": top["inferred_topology"]})
        top0 = max((s for s in scored if "composite_h0" in s), key=lambda s: s["composite_h0"], default=None)
        if top0 is not None:
            row["inferred_h0"] = top0["inferred_h0"]
            if not phased:
                row.update({"inferred": top0["inferred_h0"], "inferred_weight": top0.get("inferred_weight"),
                            "inferred_topology": top0.get("inferred_topology")})
        rows.append(row)
    return rows


# ---------------------------------------------------------------------------------------------- #
# sample level
# ---------------------------------------------------------------------------------------------- #

def cn_explained(graph: Optional[aaf.AAGraph], cycles: Sequence[Cycle], ks=(1, 3), min_cn: float = CN_EXPLAINED_MIN_CN) -> Dict[str, object]:
    """CoRAL's copy-number ratio: the share of the graph's length-weighted copy number
    (sequence edges with copy number >= ``min_cn``) explained by the k heaviest cycles (copy
    count x length, with multiplicity; paths left out, as CoRAL scores cyclic walks only).
    Unlike CoRAL's code, which returns 1 when there are fewer than k cycles, the share is
    computed whatever the number of cycles."""
    if graph is None:
        return {"cn_explained_top%d" % k: None for k in ks}
    total = sum(s.cn * (s.end - s.start + 1) for s in graph.sequences if s.cn >= min_cn)
    weights = sorted((c.weight * cycle_length(c.steps) for c in cycles if c.circular), reverse=True)
    return {"cn_explained_top%d" % k: (sum(weights[:k]) / total if total else None) for k in ks}


@dataclass
class ReadPath:
    pid: str
    junctions: List[Tuple[Breakend, Breakend]]  # (exit, entry) in the CStep frame
    support: float


def _eid(token: str) -> Breakend:
    eid = token.split("@")[0]
    chrom, rest = eid.rsplit(":", 1)
    return chrom, int(rest[:-1]), rest[-1]


def _reference_join(x: Breakend, y: Breakend, tol: int = 2) -> bool:
    if x[0] != y[0] or x[2] == y[2]:
        return False
    plus, minus = (x, y) if x[2] == "+" else (y, x)
    return abs(minus[1] - (plus[1] + 1)) <= tol


def read_paths(path: str) -> List[ReadPath]:
    """Multi-junction long-read paths (``paths.tsv``: ``PID PATH NSTEP #NOV``; ``PATH`` is
    ``EID@ORIENT`` steps, steps 2i+1 and 2i+2 one adjacency): those with at least two novel
    junctions, as ``rck-ecdna`` uses them."""
    out = []
    table = pd.read_csv(path, sep="\t")
    for r in table.to_dict("records"):
        steps = [_eid(t) for t in str(r["PATH"]).split(",") if t]
        joins = [(steps[i], steps[i + 1]) for i in range(1, len(steps) - 1, 2)]
        novel = [(a, b) for a, b in joins if not _reference_join(a, b)]
        if len(novel) >= 2:
            out.append(ReadPath(str(r["PID"]), novel, float(r.get("#NOV", 1) or 1)))
    return out


def _cycle_junctions(c: Cycle) -> List[Tuple[Breakend, Breakend]]:
    """The cycle's novel joins in traversal order (exit, entry), haplotype-blind, so a crossing
    between reference neighbours is a reference join."""
    steps = strip_haplotypes(c).steps
    return [(a, b) for a, b, _, __ in ev.junctions(steps, circular=c.circular)]


def path_in_cycle(p: ReadPath, junctions: Sequence[Tuple[Breakend, Breakend]], circular: bool, tol: int) -> bool:
    """The read's junctions appear consecutively, in order, in the cycle (either direction)."""
    n, m = len(junctions), len(p.junctions)
    if m > n and not circular or n == 0:
        return False
    for seq in (p.junctions, [(b, a) for a, b in reversed(p.junctions)]):
        for i in range(n if circular else n - m + 1):
            if all(_be_close(seq[k][0], junctions[(i + k) % n][0], tol) and _be_close(seq[k][1], junctions[(i + k) % n][1], tol)
                   for k in range(m)):
                return True
    return False


def subwalks_satisfied(paths: Sequence[ReadPath], cycles: Sequence[Cycle], tol: int) -> Dict[str, object]:
    """Share of multi-junction long-read paths found, in order, in some cycle (by count and
    weighted by read support)."""
    if not paths:
        return {"n_subwalks": 0, "subwalks_satisfied": None, "subwalks_satisfied_support": None}
    js = [(_cycle_junctions(c), c.circular) for c in cycles]
    ok = [any(path_in_cycle(p, j, circ, tol) for j, circ in js) for p in paths]
    support = sum(p.support for p in paths)
    return {"n_subwalks": len(paths), "subwalks_satisfied": sum(ok) / len(paths),
            "subwalks_satisfied_support": sum(p.support for p, k in zip(paths, ok) if k) / support if support else None}


# ---------------------------------------------------------------------------------------------- #
# strict rule
# ---------------------------------------------------------------------------------------------- #

def strict_select(records: aaf.AACycles, cycles: Sequence[Cycle], rule: str = "auto") -> Tuple[List[Cycle], str]:
    """The method's own answer among its cycles, and the rule applied.

    * ``status``: leave out cycles flagged ``Status=low_weight`` (rck-ecdna's own answer);
    * ``cn-frac:F``: cycles with copy count >= F x the heaviest copy count;
    * ``top-cn:X``: the fewest heaviest cycles (copy count x length) explaining >= X of the
      cycles' total;
    * ``all``: every cycle;
    * ``auto`` (default): ``status`` when the file carries ``Status=``, else ``cn-frac:0.1``.
    """
    fields = {rec.fields.get("Species") or "cycle{k}".format(k=rec.cycle_id): rec.fields for rec in records.cycles}
    if rule == "auto":
        rule = "status" if any("Status" in f for f in fields.values()) else "cn-frac:%g" % STRICT_CN_FRAC
    if rule == "all":
        return list(cycles), rule
    if rule == "status":
        return [c for c in cycles if fields.get(c.name, {}).get("Status") != "low_weight"], rule
    kind, _, value = rule.partition(":")
    x = float(value)
    if kind == "cn-frac":
        top = max((c.weight for c in cycles), default=0.0)
        return [c for c in cycles if c.weight >= x * top], rule
    if kind == "top-cn":
        mass = sorted(cycles, key=lambda c: c.weight * cycle_length(c.steps), reverse=True)
        total = sum(c.weight * cycle_length(c.steps) for c in mass)
        out, acc = [], 0.0
        for c in mass:
            if total and acc >= x * total:
                break
            out.append(c)
            acc += c.weight * cycle_length(c.steps)
        return out, rule
    raise ValueError("unknown --strict-rule {r}".format(r=rule))


# ---------------------------------------------------------------------------------------------- #
# driver
# ---------------------------------------------------------------------------------------------- #

def regions_of(records: aaf.AACycles, graph: Optional[aaf.AAGraph]) -> List[ev.Region]:
    """The method's regions in the CStep frame (0-based inclusive): its ``Interval`` lines, else
    its graph's segments merged across gaps up to 1 Mb (as rck-ecdna's regions)."""
    if records.intervals:
        return [(c, s - 1, e - 1) for c, s, e in records.intervals]
    out: List[ev.Region] = []
    if graph is None:
        return out
    for s in sorted(graph.sequences, key=lambda x: (x.chrom, x.start)):
        if out and out[-1][0] == s.chrom and s.start - 1 - out[-1][2] <= 1_000_000:
            out[-1] = (s.chrom, out[-1][1], max(out[-1][2], s.end - 1))
        else:
            out.append((s.chrom, s.start - 1, s.end - 1))
    return out


@dataclass
class ScoreResult:
    graph: pd.DataFrame
    species: pd.DataFrame
    sample: pd.DataFrame
    truth: List[Cycle] = field(default_factory=list)  # in-scope truth species
    strict: List[Cycle] = field(default_factory=list)
    candidates: List[Cycle] = field(default_factory=list)
    info: Dict[str, object] = field(default_factory=dict)


def score(cycles_path: str, truth: Sequence[Cycle], graph_path: Optional[str] = None, alt_cycles: Sequence[str] = (),
          mixture: Optional[Sequence[MixtureAmplicon]] = None, segment_cn: Optional[str] = None, paths: Optional[str] = None,
          strict_rule: str = "auto", species_tol: int = SPECIES_TOL, graph_tol: int = GRAPH_TOL,
          min_truth_fraction: float = 0.5, all_truth_species: bool = False, cn_min: float = CN_EXPLAINED_MIN_CN) -> ScoreResult:
    """Score one method answer. ``truth`` is the species truth (evaluate.load_truth)."""
    records = aaf.read_cycles(cycles_path)
    cycles = aaf.to_cycles(records)
    graph = aaf.read_graph(graph_path) if graph_path else None
    phased = records.phased
    strict, rule = strict_select(records, cycles, strict_rule)
    candidates = list(cycles)
    for k, path in enumerate(alt_cycles, start=1):
        alt = aaf.read_cycles(path)
        phased = phased and alt.phased
        candidates += aaf.to_cycles(alt, prefix="alt{k}.".format(k=k))
    regions = None if all_truth_species else regions_of(records, graph)
    scored = []
    for t in truth:
        frac = ev.fraction_in_regions(t.steps, regions, pad=species_tol) if regions else 1.0
        if frac >= min_truth_fraction:
            scored.append(t)
    rows, counts = strict_rows(scored, strict, phased, species_tol, graph_tol)
    rows += best_of_rows(scored, candidates, phased, species_tol, graph_tol)
    species = pd.DataFrame(rows, columns=SPECIES_COLUMNS)

    graph_rows, mixture_graphs = [], []
    if mixture is not None and graph is not None:
        cn = pd.read_csv(segment_cn, sep="\t") if segment_cn else None
        names = {t.name.split(".")[0] for t in scored}
        for amp in mixture:
            if names & set(amp.species):
                graph_rows.append(graph_metrics(amp, graph, cn, graph_tol))
                mixture_graphs.append(amp.graph.path)
    graph_frame = pd.DataFrame(graph_rows)

    read_paths_ = read_paths(paths) if paths else []
    explained = cn_explained(graph, cycles, min_cn=cn_min)
    sample_rows = []
    for mode, chosen in (("strict", strict), ("best_of", candidates)):
        mine = species[species["mode"] == mode]
        row = {"mode": mode, "strict_rule": rule, "phased": phased, "n_truth": len(scored), "n_inferred": len(chosen),
               "species_count_error": len(chosen) - len(scored) if mode == "strict" else None,
               "n_missed": counts["n_missed"] if mode == "strict" else int((mine["composite_h0"].fillna(0) <= 0).sum()),
               "n_false": counts["n_false"] if mode == "strict" else None,
               "mean_composite": mine["composite"].astype(float).mean() if phased and len(mine) else None,
               "mean_composite_h0": mine["composite_h0"].astype(float).mean() if len(mine) else None,
               "mean_cycle_overlap_coral_code": mine["cycle_overlap_coral_code"].astype(float).mean() if len(mine) else None,
               "mean_cyclic_lcs_coral_code": mine["cyclic_lcs_coral_code"].astype(float).mean() if len(mine) else None,
               "mean_cycle_overlap_supp": mine["cycle_overlap_supp"].astype(float).mean() if len(mine) else None,
               "mean_cyclic_lcs_supp": mine["cyclic_lcs_supp"].astype(float).mean() if len(mine) else None,
               "n_exact_cycles": int(mine["exact_cycle"].fillna(0).sum()) if len(mine) else 0}
        row.update(explained)
        row.update(subwalks_satisfied(read_paths_, chosen, graph_tol))
        if len(graph_frame):
            row.update({"error_free": int(graph_frame["error_free"].all()),
                        "error_free_merged": int(graph_frame["error_free_merged"].all()), "bp_accuracy": graph_frame["bp_accuracy"].mean(),
                        "novel_adjacency_fnr": graph_frame["novel_adjacency_fnr"].mean(),
                        "novel_adjacency_fdr": graph_frame["novel_adjacency_fdr"].mean(),
                        "segment_cn_mae": graph_frame["segment_cn_mae"].mean(), "is_violations": int(graph_frame["is_violations"].sum())})
        sample_rows.append(row)
    info = {"cycles": cycles_path, "graph": graph_path, "alt_cycles": list(alt_cycles), "strict_rule": rule, "phased": phased,
            "species_tol_bp": species_tol, "graph_tol_bp": graph_tol, "n_truth_total": len(truth),
            "truth_in_scope": [t.name for t in scored], "regions": regions,
            "mixture_graphs": mixture_graphs}
    return ScoreResult(graph_frame, species, pd.DataFrame(sample_rows), list(scored), list(strict), candidates, info)


def write(outdir: str, result: ScoreResult) -> Dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    paths = {"graph": os.path.join(outdir, "aa_graph_metrics.tsv"), "species": os.path.join(outdir, "aa_species_metrics.tsv"),
             "sample": os.path.join(outdir, "aa_sample_metrics.tsv"), "summary": os.path.join(outdir, "aa_eval.json")}
    result.graph.to_csv(paths["graph"], sep="\t", index=False)
    result.species.to_csv(paths["species"], sep="\t", index=False)
    result.sample.to_csv(paths["sample"], sep="\t", index=False)
    with open(paths["summary"], "w") as handle:
        json.dump({"info": result.info, "sample": ev._json_records(result.sample), "graph": ev._json_records(result.graph)},
                  handle, indent=2, default=str)
    return paths


def default_mixture_dir(truth_dir: str) -> Optional[str]:
    d = os.path.join(truth_dir, "truth", "aa")
    return d if glob.glob(os.path.join(d, "mixture_amplicons.tsv")) else None


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rck.ecdna.aa_eval", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--cycles", required=True, help="the method's AA cycles file (its answer)")
    parser.add_argument("--alt-cycles", nargs="*", default=[], help="other answers of the method, for best-of only")
    parser.add_argument("--graph", default=None, help="the method's AA graph file (graph-level metrics, cn_explained)")
    parser.add_argument("--segment-cn", default=None, help="allele copy numbers per segment: chrom start end [cn] cn_a cn_b (AA coordinates)")
    parser.add_argument("--paths", default=None, help="long-read paths (subgraph paths.tsv) for subwalks_satisfied")
    parser.add_argument("--truth-dir", default=None, help="simulator run directory")
    parser.add_argument("--truth-cycles", default=None, help="explicit truth table, as rck-ecdna-eval")
    parser.add_argument("--truth-sample", default=None)
    parser.add_argument("--mixture-dir", default=None, help="truth mixture files (default TRUTH_DIR/truth/aa when present)")
    parser.add_argument("--strict-rule", default="auto", help="auto, status, all, cn-frac:F or top-cn:X (default auto)")
    parser.add_argument("--species-tol", type=int, default=SPECIES_TOL, help="species-level tolerance, bp (default 10000)")
    parser.add_argument("--graph-tol", type=int, default=GRAPH_TOL, help="graph-level and exact-cycle tolerance, bp (default 100)")
    parser.add_argument("--min-truth-fraction", type=float, default=0.5)
    parser.add_argument("--all-truth-species", action="store_true")
    parser.add_argument("--cn-explained-min-cn", type=float, default=CN_EXPLAINED_MIN_CN)
    parser.add_argument("--outdir", required=True)
    args = parser.parse_args(argv)
    logging.basicConfig(level=logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    truth = ev.load_truth(args.truth_dir, args.truth_cycles, args.truth_sample)
    mixture_dir = args.mixture_dir or (default_mixture_dir(args.truth_dir) if args.truth_dir else None)
    mixture = load_mixture(mixture_dir) if mixture_dir else None
    if mixture is None:
        log.warning("no truth mixture graph: graph-level metrics skipped (build it with simulations/truth_mixture.py)")
    result = score(args.cycles, truth, args.graph, args.alt_cycles, mixture, args.segment_cn, args.paths, args.strict_rule,
                   args.species_tol, args.graph_tol, args.min_truth_fraction, args.all_truth_species, args.cn_explained_min_cn)
    write(args.outdir, result)
    with pd.option_context("display.width", 200, "display.max_columns", 40):
        print(result.sample.to_string(index=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

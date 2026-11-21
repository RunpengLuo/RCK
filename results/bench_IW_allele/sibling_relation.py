#!/usr/bin/env python3
"""PLAN.md item 1: the sibling-vs-outer haplotype relation, per base, truth vs each finished rck_mc answer.

The within-phase-block scores (recovery_matrix.tsv *_blockwise) pick the A/B swap per phase block for each
truth species separately, so they could forgive a sibling put on its outer species' allele. This checks the
relation itself: on every base the two species share, is the sibling on the opposite or the same haplotype as
the outer species, in the truth (truth_cycles.tsv) and in the answer (the inferred species matched to each
truth species in recovery_matrix.tsv, rank-1 decomposition)? A relation is unaffected by phase-block
orientation, so no swap is needed. Local, no Gurobi, seconds.

    python sibling_relation.py            # -> sibling_relation.tsv here, and a printed table

Columns: case, scenario, stratum, truth_opp_bp / truth_same_bp (shared bases), inferred_opp_bp / inferred_same_bp,
agree_bp / disagree_bp (shared bases where the inferred relation equals / differs from the truth's),
matched (both truth species matched to two different inferred species), note.
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allelic_bench as ab  # noqa: E402

bi = ab.bi
OUTER, SIBLING = "ecdna1", "ecdna2"  # make_scenarios iw_allelic: allelic_sibling [null, 0]


def _steps_from_walk(walk):
    out = []
    for tok in str(walk).split(","):
        chrom, span, hap, _strand = tok.rsplit(":", 3)
        a, b = span.split("-")
        out.append((chrom, int(a), int(b) + 1, hap))  # walk ends are inclusive
    return out


def _merge(intervals):
    """{(chrom, hap): merged half-open intervals} of (chrom, start, end, hap) steps (a base counted once)."""
    by = {}
    for chrom, a, b, hap in intervals:
        by.setdefault((chrom, hap), []).append((a, b))
    out = {}
    for key, spans in by.items():
        spans.sort()
        merged = []
        for a, b in spans:
            if merged and a <= merged[-1][1]:
                merged[-1][1] = max(merged[-1][1], b)
            else:
                merged.append([a, b])
        out[key] = merged
    return out


def _overlap(a, b):
    """Total length of the intersection of two merged interval lists."""
    i = j = total = 0
    while i < len(a) and j < len(b):
        lo, hi = max(a[i][0], b[j][0]), min(a[i][1], b[j][1])
        total += max(0, hi - lo)
        if a[i][1] < b[j][1]:
            i += 1
        else:
            j += 1
    return total


def relation_intervals(sib, outer):
    """{(chrom, 'opp' | 'same'): merged intervals} of the bases the two species share."""
    s, o = _merge(sib), _merge(outer)
    pieces = []
    for (chrom, hs), si in s.items():
        for (chrom2, ho), oi in o.items():
            if chrom != chrom2:
                continue
            i = j = 0
            while i < len(si) and j < len(oi):
                lo, hi = max(si[i][0], oi[j][0]), min(si[i][1], oi[j][1])
                if hi > lo:
                    pieces.append((chrom, lo, hi, "same" if hs == ho else "opp"))
                if si[i][1] < oi[j][1]:
                    i += 1
                else:
                    j += 1
    return _merge(pieces)


def bp(rel, kind):
    return sum(b - a for (chrom, k), spans in rel.items() if k == kind for a, b in spans)


def score_case(case):
    run = bi.run_dir(case, "rck_mc")
    row = {"case": case["name"], "scenario": case["scenario"], "stratum": case["stratum"]}
    if not bi.done(case, "rck_mc"):
        return dict(row, note="rck_mc not finished")
    truth = pd.read_csv(os.path.join(run, "truth_cycles.tsv"), sep="\t")
    t = {sp: [(r.chrom, int(r.start), int(r.end) + 1, r.hap) for r in g.itertuples()] for sp, g in truth.groupby("species")}
    dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t")
    dec = dec[(dec["decomposition_rank"] == 1) & (dec["kind"] == "cycle")]
    inferred = {r.species: _steps_from_walk(r.walk) for r in dec.itertuples()}
    rec = pd.read_csv(os.path.join(run, "recovery_matrix.tsv"), sep="\t")
    match = {r.truth: r.inferred for r in rec.itertuples() if isinstance(r.inferred, str)}
    rt = relation_intervals(t[SIBLING], t[OUTER])
    row.update(truth_opp_bp=bp(rt, "opp"), truth_same_bp=bp(rt, "same"))
    a, b = match.get(SIBLING), match.get(OUTER)
    row["matched"] = bool(a and b and a != b and a in inferred and b in inferred)
    if not row["matched"]:
        return dict(row, note="truth species not matched to two inferred species ({m})".format(m=match))
    ri = relation_intervals(inferred[a], inferred[b])
    row.update(inferred_opp_bp=bp(ri, "opp"), inferred_same_bp=bp(ri, "same"))
    agree = sum(_overlap(rt.get((c, k), []), ri.get((c, k), [])) for (c, k) in rt)
    disagree = sum(_overlap(rt.get((c, k), []), ri.get((c, other), [])) for (c, k) in rt
                   for other in ("opp", "same") if other != k)
    return dict(row, agree_bp=agree, disagree_bp=disagree, note="")


def main():
    rows = [score_case(c) for c in ab.cases()]
    table = pd.DataFrame(rows)
    table.to_csv(os.path.join(HERE, "sibling_relation.tsv"), sep="\t", index=False)
    with pd.option_context("display.width", 200, "display.max_columns", 20):
        print(table.to_string(index=False))
    return 0


if __name__ == "__main__":
    sys.exit(main())

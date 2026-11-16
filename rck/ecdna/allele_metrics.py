"""Allele-specific ecDNA copy numbers of one answer, and post-hoc phasing of unphased answers.

Allele-specific copy number (ecDNA only, base-pair weighted). For haplotype h in {A, B}, over every
reference base x:

    t_h(x) = sum_s PHI_s * (times truth species s traverses x on haplotype h)
    c_h(x) = sum_c w_c   * (times inferred cycle c traverses x on haplotype h)
    c_U(x) = the same over unphased steps (haplotype ``.``: CoRAL / CE without post-hoc phasing)

and t = t_A + t_B, c = c_A + c_B + c_U. Over the bases with t > 0 or c > 0:

    as_cn_recall      = sum_x [min(c_A, t_A) + min(c_B, t_B)] / sum_x t
    as_cn_precision   = sum_x [min(c_A, t_A) + min(c_B, t_B)] / sum_x c   (unphased copies count as wrong)
    as_cn_mae         = mean over truth bases (t > 0) of |c_A - t_A| + |c_B - t_B|
    as_cn_within_tol  = share of truth bases where both alleles are within max(1, 0.1 t_h)
    major_allele_acc  = share of truth bases with t_A != t_B where sign(c_A - c_B) = sign(t_A - t_B)
    baf_mae           = mean over truth bases with c_A + c_B > 0 of |c_B / (c_A + c_B) - t_B / t|
    mixed_bp_frac     = share of truth bases with t_A > 0 and t_B > 0 (a property of the case: on these
                        bases no rule that gives a segment one haplotype can be right for every copy)

``as_cn_recall`` is the phased version of ``bench_metrics.cn_recall``: the difference between the two
is the ecDNA copy mass the answer places on the wrong allele.

Within phase blocks (``within_block``): the observed A / B of a phase block is arbitrary, so the
answer's A and B may be swapped once per block (shared by every cycle), keeping the swap that gives
the lower allele error, as ``evaluate.block_swaps`` does for the species scores.

Post-hoc phasing (``rephase``): give every step of an unphased answer the haplotype a rule picks for
its reference bases, so a haplotype-blind method's answer can be scored with the phased metrics.
The rules give one haplotype per reference segment, shared by every cycle and every traversal:

    ``baf_major``    the major allele of the observed allele-specific copy numbers (the subgraph's
                     ``segment_blocks.tsv``: CN_A / CN_B from depth and phased SNP counts, fitted
                     before any decomposition; what a CoRAL user with phased SNPs would have)
    ``truth_major``  the major allele of the true ecDNA copies (t_A vs t_B): the best any
                     one-haplotype-per-segment rule can do, an upper bound for post-hoc phasing

A segment the rule cannot call (no SNPs, or no ecDNA copies in the truth) takes the haplotype of the
nearest called step of the same cycle (``rephase(fill=True)``, the default), as a user would carry a
phase across a SNP-less segment; the number of filled steps is reported.

    python -m rck.ecdna.allele_metrics --cycles ANSWER_cycles.txt --truth-dir SIM_RUN --outdir OUT
        [--weights-from-multiclone RUN_DIR --clone-proportions DIAG/clone_proportions.tsv]
        [--phase-blocks rck_input/rck.fragments.tsv]

Writes OUT/allele_metrics.json and OUT/allele_segments.tsv (constant intervals: t_A, t_B, c_A, c_B, c_U).
"""
import argparse
import json
import logging
import os
from collections import defaultdict
from typing import Callable, Dict, List, Optional, Sequence, Tuple

import pandas as pd

from rck.ecdna import aa_eval
from rck.ecdna import aa_format as aaf
from rck.ecdna import bench_metrics as bm
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import Cycle

log = logging.getLogger(__name__)

HAPS = ("A", "B")
COLS = ("t_A", "t_B", "c_A", "c_B", "c_U")


# ---------------------------------------------------------------------------------------- #
# allele-specific copy profiles
# ---------------------------------------------------------------------------------------- #

def _events(cycles: Sequence[Cycle], key: Callable[[str], Optional[str]]) -> Dict[str, List[tuple]]:
    """chrom -> [(position, column, change)] of sum_c w_c * traversals, steps 0-based inclusive."""
    out: Dict[str, List[tuple]] = defaultdict(list)
    for c in cycles:
        for s in c.steps:
            col = key(s.hap)
            if col is None:
                continue
            lo, hi = min(s.start, s.end), max(s.start, s.end)
            out[s.chrom].append((lo, col, float(c.weight)))
            out[s.chrom].append((hi + 1, col, -float(c.weight)))
    return out


def _truth_col(hap: str) -> Optional[str]:
    return "t_" + hap if hap in HAPS else None


def _inferred_col(hap: str) -> str:
    return "c_" + hap if hap in HAPS else "c_U"


def allele_table(truth: Sequence[Cycle], inferred: Sequence[Cycle]) -> pd.DataFrame:
    """Constant intervals (chrom, start, end: 0-based inclusive) of t_A, t_B, c_A, c_B, c_U where any is > 0."""
    ev_t, ev_c = _events(truth, _truth_col), _events(inferred, _inferred_col)
    rows = []
    for chrom in sorted(set(ev_t) | set(ev_c)):
        changes = ev_t.get(chrom, []) + ev_c.get(chrom, [])
        points = sorted({p for p, _, _ in changes})
        delta: Dict[int, Dict[str, float]] = defaultdict(lambda: defaultdict(float))
        for p, col, d in changes:
            delta[p][col] += d
        level = dict.fromkeys(COLS, 0.0)
        for left, right in zip(points, points[1:]):
            for col, d in delta[left].items():
                level[col] += d
                if abs(level[col]) < 1e-9:
                    level[col] = 0.0
            if any(level.values()):
                rows.append(dict(chrom=chrom, start=left, end=right - 1, **level))
    frame = pd.DataFrame(rows, columns=["chrom", "start", "end", *COLS])
    # merge neighbours with equal levels (cosmetic: fewer rows, same scores)
    if len(frame):
        same = (frame["chrom"] == frame["chrom"].shift()) & (frame["start"] == frame["end"].shift() + 1)
        for col in COLS:
            same &= frame[col] == frame[col].shift()
        frame = frame.assign(grp=(~same).cumsum()).groupby("grp").agg(
            chrom=("chrom", "first"), start=("start", "first"), end=("end", "last"),
            **{c: (c, "first") for c in COLS}).reset_index(drop=True)
    return frame


def allele_scores(table: pd.DataFrame) -> Dict[str, Optional[float]]:
    """The module docstring's scores of one allele_table."""
    if table.empty:
        return {}
    w = (table["end"] - table["start"] + 1).astype(float)
    t_a, t_b, c_a, c_b, c_u = (table[c] for c in COLS)
    t, c = t_a + t_b, c_a + c_b + c_u
    overlap = (pd.concat([c_a, t_a], axis=1).min(axis=1) + pd.concat([c_b, t_b], axis=1).min(axis=1)) * w
    blind_overlap = pd.concat([c, t], axis=1).min(axis=1) * w
    truth = t > 0
    tol_a = (c_a - t_a).abs() <= (bm.TOL_REL * t_a).clip(lower=bm.TOL_ABS)
    tol_b = (c_b - t_b).abs() <= (bm.TOL_REL * t_b).clip(lower=bm.TOL_ABS)
    unequal = truth & (t_a != t_b)
    major_ok = unequal & ((c_a - c_b) * (t_a - t_b) > 0)
    called = truth & ((c_a + c_b) > 0)
    baf_err = ((c_b / (c_a + c_b)) - (t_b / t)).abs()
    mixed = truth & (t_a > 0) & (t_b > 0)

    def wmean(values, mask):
        return float((values[mask] * w[mask]).sum() / w[mask].sum()) if w[mask].sum() else None

    mass_t, mass_c = float((t * w).sum()), float((c * w).sum())
    return {
        "as_cn_recall": float(overlap.sum()) / mass_t if mass_t else None,
        "as_cn_precision": float(overlap.sum()) / mass_c if mass_c else None,
        "blind_cn_recall": float(blind_overlap.sum()) / mass_t if mass_t else None,
        "as_cn_mae": wmean((c_a - t_a).abs() + (c_b - t_b).abs(), truth),
        "as_cn_within_tol": wmean((tol_a & tol_b).astype(float), truth),
        "major_allele_acc": float(w[major_ok].sum() / w[unequal].sum()) if w[unequal].sum() else None,
        "baf_mae": wmean(baf_err, called),
        "unphased_mass_frac": float((c_u * w).sum()) / mass_c if mass_c else None,
        "mixed_bp_frac": float(w[mixed].sum() / w[truth].sum()) if w[truth].sum() else None,
        "truth_bp": int(w[truth].sum()),
        "truth_mass": mass_t, "inferred_mass": mass_c,
    }


def _block_of(chrom: str, start: int, end: int, blocks: Sequence[ev.Block]) -> Optional[str]:
    mid = (start + end) // 2
    for b_chrom, b_start, b_end, name in blocks:
        if b_chrom == chrom and b_start <= mid <= b_end:
            return name
    return None


def within_block(table: pd.DataFrame, blocks: Sequence[ev.Block]) -> Tuple[pd.DataFrame, int]:
    """(table with c_A / c_B swapped in the blocks where that lowers the allele error, n swapped)."""
    if table.empty or not blocks:
        return table, 0
    t = table.copy()
    t["block"] = [_block_of(r.chrom, r.start, r.end, blocks) for r in t.itertuples()]
    w = t["end"] - t["start"] + 1
    keep = ((t["c_A"] - t["t_A"]).abs() + (t["c_B"] - t["t_B"]).abs()) * w
    swap = ((t["c_B"] - t["t_A"]).abs() + (t["c_A"] - t["t_B"]).abs()) * w
    gain = (keep - swap).groupby(t["block"]).sum()
    swapped = set(gain[gain > 1e-9].index)
    rows = t["block"].isin(swapped).to_numpy()
    c_a, c_b = t["c_A"].to_numpy(copy=True), t["c_B"].to_numpy(copy=True)
    # NB: not t.loc[rows, [A, B]] = t.loc[rows, [B, A]].values: on a table read from a TSV (pandas
    # copy-on-write) the right side can be a view that the assignment overwrites, zeroing both
    t["c_A"], t["c_B"] = [x.copy() for x in (c_a, c_b)]
    t.loc[rows, "c_A"] = c_b[rows]
    t.loc[rows, "c_B"] = c_a[rows]
    return t.drop(columns="block"), len(swapped)


# ---------------------------------------------------------------------------------------- #
# post-hoc phasing
# ---------------------------------------------------------------------------------------- #

Rule = Callable[[str, int, int], str]  # (chrom, start, end: 0-based inclusive) -> haplotype


def _vote_rule(rows: Sequence[tuple]) -> Rule:
    """rows [(chrom, lo, hi, score)] (0-based inclusive): score > 0 votes A, < 0 votes B, weighted by
    overlapping bases; no informative overlap gives ``.``."""
    by: Dict[str, List[tuple]] = defaultdict(list)
    for chrom, lo, hi, score in rows:
        by[str(chrom)].append((int(lo), int(hi), float(score)))

    def rule(chrom: str, start: int, end: int) -> str:
        total = 0.0
        for lo, hi, score in by.get(chrom, []):
            n = min(hi, end) - max(lo, start) + 1
            if n > 0 and score == score:  # not NaN
                total += n * (1.0 if score > 0 else -1.0 if score < 0 else 0.0)
        return "A" if total > 0 else "B" if total < 0 else aaf.UNPHASED
    return rule


def baf_major_rule(segment_blocks: str) -> Rule:
    """Major allele of a subgraph's segment_blocks.tsv (piece_start / piece_end 1-based, CN_A / CN_B)."""
    t = pd.read_csv(segment_blocks, sep="\t")
    t = t[t["CN_A"].notna() & t["CN_B"].notna()]
    return _vote_rule([(r["CHROM_1"], int(r["piece_start"]) - 1, int(r["piece_end"]) - 1, float(r["CN_A"]) - float(r["CN_B"]))
                       for r in t.to_dict("records")])


def truth_major_rule(truth: Sequence[Cycle]) -> Rule:
    """Major allele of the true ecDNA copies: t_A - t_B per base."""
    table = allele_table(truth, [])
    return _vote_rule([(r.chrom, r.start, r.end, r.t_A - r.t_B) for r in table.itertuples()])


def _fill(labels: List[str], circular: bool) -> Tuple[List[str], int]:
    """``.`` labels replaced by the nearest labelled step of the same walk (the previous one on a tie,
    cyclically on a cycle): a segment without informative SNPs takes its neighbours' haplotype."""
    known = [i for i, h in enumerate(labels) if h != aaf.UNPHASED]
    if not known or len(known) == len(labels):
        return labels, 0
    n, out, filled = len(labels), list(labels), 0
    for i, h in enumerate(labels):
        if h != aaf.UNPHASED:
            continue
        dist = (lambda j: min(abs(i - j), n - abs(i - j))) if circular else (lambda j: abs(i - j))
        best = min(known, key=lambda j: (dist(j), (i - j) % n if circular else j > i))
        out[i] = labels[best]
        filled += 1
    return out, filled


def rephase(cycles_path: str, rule: Rule, out_path: str, fill: bool = True) -> Dict[str, int]:
    """Write cycles_path with every Haplotypes= (and Alt_haplotypes=) set by ``rule``, and with ``fill``
    the steps the rule leaves at ``.`` given their nearest labelled neighbour's haplotype (``_fill``).
    Returns the count of steps per final label, plus ``filled``."""
    records = aaf.read_cycles(cycles_path)
    counts: Dict[str, int] = defaultdict(int)

    def haps(segments_text: str) -> str:
        tokens = [t.strip() for t in segments_text.split(",") if t.strip()]
        steps = [i for i, t in enumerate(tokens) if t[:-1] != "0"]
        labels = []
        for i in steps:
            chrom, start, end = records.segments[tokens[i][:-1]]
            labels.append(rule(chrom, start - 1, end - 1))
        if fill:
            labels, n = _fill(labels, circular=len(steps) == len(tokens))
            counts["filled"] += n
        out = [aaf.UNPHASED] * len(tokens)
        for i, h in zip(steps, labels):
            out[i] = h
            counts[h] += 1
        return ",".join(out)

    os.makedirs(os.path.dirname(os.path.abspath(out_path)), exist_ok=True)
    with open(cycles_path) as src, open(out_path, "w") as dst:
        for line in src:
            if not line.startswith("Cycle="):
                dst.write(line)
                continue
            parts = [p for p in line.rstrip("\n").split(";") if p]
            fields = dict(p.split("=", 1) for p in parts if "=" in p)
            fields["Haplotypes"] = haps(fields.get("Segments", ""))
            if fields.get("Alt_segments"):
                fields["Alt_haplotypes"] = "|".join(haps(s) for s in fields["Alt_segments"].split("|"))
            order = [p.split("=", 1)[0] for p in parts if "=" in p]
            order += [k for k in ("Haplotypes", "Alt_haplotypes") if k in fields and k not in order]
            dst.write(";".join("{k}={v}".format(k=k, v=fields[k]) for k in order) + "\n")
    return dict(counts)


# ---------------------------------------------------------------------------------------- #
# one answer
# ---------------------------------------------------------------------------------------- #

def strict_cycles(cycles_path: str, weights: Optional[Dict[str, float]] = None, strict_rule: str = "auto"):
    """(strict cycles, rule) as bench_metrics.score selects them, weights overridden by name."""
    records = aaf.read_cycles(cycles_path)
    cycles = aa_eval.species_cycles(records)
    strict, rule = aa_eval.strict_select(records, cycles, strict_rule)
    if weights is not None:
        strict = bm.reweight(strict, weights)
    return strict, rule


def score(cycles_path: str, truth: Sequence[Cycle], weights: Optional[Dict[str, float]] = None,
          blocks: Optional[Sequence[ev.Block]] = None) -> Tuple[Dict[str, object], pd.DataFrame]:
    """(summary, allele table) of one answer: absolute scores, and _block ones when blocks are given."""
    strict, rule = strict_cycles(cycles_path, weights)
    table = allele_table(truth, strict)
    summary: Dict[str, object] = {"strict_rule": rule, "n_strict": len(strict)}
    summary.update(allele_scores(table))
    if blocks:
        swapped, n = within_block(table, blocks)
        summary.update({k + "_block": v for k, v in allele_scores(swapped).items()
                        if k not in ("mixed_bp_frac", "truth_bp", "truth_mass", "inferred_mass", "blind_cn_recall")})
        summary["n_blocks_swapped"] = n
    return summary, table


def main(argv=None) -> int:
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--cycles", required=True)
    p.add_argument("--truth-dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--weights-from-multiclone", default=None, help="multi-clone run: w_c = sum_k u_k psi_k")
    p.add_argument("--clone-proportions", default=None)
    p.add_argument("--phase-blocks", default=None, help="rck_input/rck.fragments.tsv: also score within phase blocks")
    a = p.parse_args(argv)
    weights = None
    if a.weights_from_multiclone:
        weights = bm.multiclone_bulk_weights(a.weights_from_multiclone, bm.read_proportions(a.clone_proportions))
    blocks = ev.read_phase_blocks(a.phase_blocks) if a.phase_blocks else None
    summary, table = score(a.cycles, ev.load_truth(a.truth_dir), weights, blocks)
    os.makedirs(a.outdir, exist_ok=True)
    with open(os.path.join(a.outdir, "allele_metrics.json"), "w") as handle:
        json.dump(summary, handle, indent=2)
    table.to_csv(os.path.join(a.outdir, "allele_segments.tsv"), sep="\t", index=False)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

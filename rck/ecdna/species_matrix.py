"""Per-cell species matrix: (J_seg, J_adj) for every (truth species, inferred species) pair.

M1 and M2 (docs/pipeline/evaluation.md §4.7) score a whole answer. This matrix shows *where* an
answer goes wrong. Rows are the in-scope truth species, ordered by copy number, high to low.
Columns are the strict rank-1 species, ordered by refitted copy number, high to low, followed by
the non-strict (``low_weight``) species, also high to low. Each cell is a tuple:

* ``j_seg``, the M1 side: the §3.1 base-pair Jaccard with multiplicity
  (``evaluate.jaccard``: haplotype Jaccard when phased, segment Jaccard when blind);
* ``j_adj``, the M2 side: |ι(T) ∩ ι(C)| / |ι(T) ∪ ι(C)| over M2's novel-adjacency items, NaN when
  the union is empty.

Items are M2's universe (``aa_eval.adjacency_items`` on the strict species). The adjacencies of
non-strict columns are added afterwards (:func:`column_items`), so that the strict cells use exactly
M2's item sets. Optional columns:
- ``j_adj_pairs``: M2's pair form C(k,2) / (C(a,2) + C(b,2) − C(k,2));
- ``j_seg_cn``: copy-weighted Σ min(w_T m_T, w_C m_C) / Σ max(w_T m_T, w_C m_C).

The merge and split flags read strict columns only.

    python -m rck.ecdna.species_matrix AA_EVAL_OUTDIR [--haplotype blind] [-o matrix.png]
"""

from __future__ import annotations

import argparse
import math
import os
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

from rck.ecdna import aa_eval as ae
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import CStep, Cycle

MATRIX_COLUMNS = ["case", "haplotype", "tol_bp", "truth", "truth_rank", "truth_weight", "column_rank", "inferred",
                  "copy_number", "low_weight", "j_seg", "j_adj", "exact_cycle", "j_adj_pairs", "j_seg_cn",
                  "n_items_truth", "n_items_inferred", "n_items_shared", "merge_column", "split_row"]


def column_items(truth: Sequence[Cycle], strict: Sequence[Cycle], extra: Sequence[Cycle], tol: int,
                 phased: bool) -> Tuple[List, List[set], List[set], List[set]]:
    """M2's items for ``truth`` and ``strict`` (``aa_eval.adjacency_items``, unchanged), extended
    with the adjacencies of ``extra`` (the non-strict columns). Each adjacency of ``extra`` takes:
    1. an existing item it equals exactly; else
    2. a truth item that no strict adjacency took, one to one and nearest first within ``tol``; else
    3. a new false item.
    This pass changes no strict set. Returns (items, truth_sets, strict_sets, extra_sets)."""
    items, truth_sets, strict_sets, n_truth = ae.adjacency_items(truth, strict, tol, phased)
    taken = {i for s in strict_sets for i in s if i < n_truth}
    e_js = [ae.walk_adjacencies(c, phased) for c in extra]
    distinct = ae._distinct([j for js in e_js for j in js], phased)
    assign: Dict[int, int] = {}
    for a, j in enumerate(distinct):
        same = next((i for i, k in enumerate(items) if ae.junction_distance(j, k, phased) == 0), None)
        if same is not None:
            assign[a] = same
    candidates = []
    for a, j in enumerate(distinct):
        if a in assign:
            continue
        for b in range(n_truth):
            if b in taken:
                continue
            d = ae.junction_distance(j, items[b], phased)
            if d is not None and d <= tol:
                candidates.append((d, a, b))
    for d, a, b in sorted(candidates):
        if a in assign or b in taken:
            continue
        assign[a] = b
        taken.add(b)
    for a, j in enumerate(distinct):
        if a not in assign:
            assign[a] = len(items)
            items.append(j)

    def index(j):
        return next(i for i, k in enumerate(distinct) if ae.junction_distance(j, k, phased) == 0)

    extra_sets = [{assign[index(j)] for j in js} for js in e_js]
    return items, truth_sets, strict_sets, extra_sets


def item_jaccard(t: set, c: set) -> float:
    """|t ∩ c| / |t ∪ c|; NaN for an empty union (M2's convention: undefined, never 1)."""
    union = len(t | c)
    return len(t & c) / union if union else float("nan")


def pair_jaccard(t: set, c: set) -> float:
    """M2's pair form on one cell: C(k,2) / (C(a,2) + C(b,2) − C(k,2)), NaN when no pair exists."""
    k, a, b = len(t & c), len(t), len(c)
    den = math.comb(a, 2) + math.comb(b, 2) - math.comb(k, 2)
    return math.comb(k, 2) / den if den else float("nan")


def weighted_jaccard(t: Sequence[CStep], wt: float, c: Sequence[CStep], wc: float, with_hap: bool) -> float:
    """Σ_x min(w_T m_T, w_C m_C) / Σ_x max(w_T m_T, w_C m_C), on ``evaluate._coverage``; 0 for an empty union."""
    ca, cb = ev._coverage(t, with_hap), ev._coverage(c, with_hap)
    num = den = 0.0
    for key in set(ca) | set(cb):
        da: Dict[int, int] = {}
        db: Dict[int, int] = {}
        for p, d in ca.get(key, []):
            da[p] = da.get(p, 0) + d
        for p, d in cb.get(key, []):
            db[p] = db.get(p, 0) + d
        points = sorted(set(da) | set(db))
        xa = xb = 0
        for p, q in zip(points, points[1:]):
            xa += da.get(p, 0)
            xb += db.get(p, 0)
            num += min(wt * xa, wc * xb) * (q - p)
            den += max(wt * xa, wc * xb) * (q - p)
    return num / den if den > 0 else 0.0


def _order(cycles: Sequence[Cycle]) -> List[Cycle]:
    return sorted(cycles, key=lambda c: (-float(c.weight), c.name))


def species_matrix(truth: Sequence[Cycle], strict: Sequence[Cycle], others: Sequence[Cycle], phased: bool,
                   species_tol: int = ae.SPECIES_TOL, graph_tol: int = ae.GRAPH_TOL, case: str = "") -> pd.DataFrame:
    """The long-format matrix of one answer: blind always, phased when the answer is phased.
    ``strict`` are the method's strict species; ``others`` its remaining (``low_weight``) cycles."""
    rows_t = _order(truth)
    cols = [(c, False) for c in _order(strict)] + [(c, True) for c in _order(others)]
    s_cols = [c for c, low in cols if not low]
    o_cols = [c for c, low in cols if low]
    out = []
    for hap in (("blind", "phased") if phased else ("blind",)):
        ph = hap == "phased"
        tw = rows_t if ph else [ae.strip_haplotypes(t) for t in rows_t]
        sw = s_cols if ph else [ae.strip_haplotypes(c) for c in s_cols]
        ow = o_cols if ph else [ae.strip_haplotypes(c) for c in o_cols]
        _, t_sets, s_sets, o_sets = column_items(tw, sw, ow, species_tol, ph)
        c_walks, c_sets = sw + ow, s_sets + o_sets
        block = []
        for i, t in enumerate(tw):
            for k, (c, low) in enumerate(cols):
                cw = c_walks[k]
                ts, cs = t_sets[i], c_sets[k]
                block.append({
                    "case": case, "haplotype": hap, "tol_bp": species_tol, "truth": t.name, "truth_rank": i + 1,
                    "truth_weight": t.weight, "column_rank": k + 1, "inferred": c.name, "copy_number": c.weight,
                    "low_weight": low, "j_seg": ev.jaccard(t.steps, cw.steps, with_hap=ph),
                    "j_adj": item_jaccard(ts, cs), "exact_cycle": ae._exact(t, cw, graph_tol),
                    "j_adj_pairs": pair_jaccard(ts, cs),
                    "j_seg_cn": weighted_jaccard(t.steps, float(t.weight), cw.steps, float(c.weight), ph),
                    "n_items_truth": len(ts), "n_items_inferred": len(cs), "n_items_shared": len(ts & cs)})
        frame = pd.DataFrame(block, columns=[c for c in MATRIX_COLUMNS if c not in ("merge_column", "split_row")])
        nz = frame[(~frame["low_weight"]) & (frame["j_adj"].fillna(0) > 0)]
        merges = nz.groupby("column_rank")["truth"].nunique()
        splits = nz.groupby("truth_rank")["column_rank"].nunique()
        frame["merge_column"] = frame["column_rank"].map(lambda r: bool(merges.get(r, 0) >= 2) and r <= len(s_cols))
        frame["split_row"] = frame["truth_rank"].map(lambda r: bool(splits.get(r, 0) >= 2))
        out.append(frame)
    return pd.concat(out, ignore_index=True) if out else pd.DataFrame(columns=MATRIX_COLUMNS)


def plot_species_matrix(frame: pd.DataFrame, path: str, haplotype: str = "blind", title: Optional[str] = None) -> Optional[str]:
    """Heatmap of one haplotype mode: two panels (J_seg, J_adj), x = column rank (``low_weight``
    columns after a dashed line, labelled ``(lw)``), y = truth species; each cell is annotated. NaN
    (undefined J_adj) is drawn grey. Returns the path, or None without matplotlib."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return None
    f = frame[frame["haplotype"] == haplotype]
    if not len(f):
        return None
    rows = f.drop_duplicates("truth_rank").sort_values("truth_rank")
    cols = f.drop_duplicates("column_rank").sort_values("column_rank")
    n_strict = int((~cols["low_weight"]).sum())
    xlabels = ["%d %s\n%.3g%s" % (r.column_rank, r.inferred, r.copy_number, " (lw)" if r.low_weight else "")
               for r in cols.itertuples()]
    ylabels = ["%s (%.3g)" % (r.truth, r.truth_weight) for r in rows.itertuples()]
    fig, axes = plt.subplots(1, 2, figsize=(max(6.0, 1.1 * len(cols) * 2 + 3), 1.0 + 0.7 * len(rows)), squeeze=False)
    cmap = plt.get_cmap("viridis").copy()
    cmap.set_bad("#cccccc")
    for ax, key, name in ((axes[0][0], "j_seg", "J_seg"), (axes[0][1], "j_adj", "J_adj")):
        grid = f.pivot(index="truth_rank", columns="column_rank", values=key).reindex(index=rows["truth_rank"], columns=cols["column_rank"])
        ax.imshow(grid.values.astype(float), vmin=0, vmax=1, cmap=cmap, aspect="auto")
        for i in range(grid.shape[0]):
            for j in range(grid.shape[1]):
                v = grid.values[i, j]
                ax.text(j, i, "–" if pd.isna(v) else "%.2f" % v, ha="center", va="center", fontsize=8,
                        color="black" if (pd.isna(v) or v > 0.5) else "white")
        if 0 < n_strict < len(cols):
            ax.axvline(n_strict - 0.5, color="red", linestyle="--", linewidth=1)
        ax.set_xticks(range(len(cols)))
        ax.set_xticklabels(xlabels, fontsize=7)
        ax.set_yticks(range(len(rows)))
        ax.set_yticklabels(ylabels, fontsize=8)
        ax.set_title("%s (%s)" % (name, haplotype), fontsize=9)
    if title:
        fig.suptitle(title, fontsize=10)
    fig.tight_layout()
    fig.savefig(path, dpi=120)
    plt.close(fig)
    return path


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(prog="python -m rck.ecdna.species_matrix", description=__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("outdir", help="an aa_eval output folder (with aa_species_matrix.tsv)")
    parser.add_argument("--haplotype", default="blind", choices=("blind", "phased"))
    parser.add_argument("-o", "--out", default=None, help="heatmap path (default OUTDIR/aa_species_matrix_<haplotype>.png)")
    args = parser.parse_args(argv)
    frame = pd.read_csv(os.path.join(args.outdir, "aa_species_matrix.tsv"), sep="\t")
    path = args.out or os.path.join(args.outdir, "aa_species_matrix_%s.png" % args.haplotype)
    print(plot_species_matrix(frame, path, args.haplotype) or "matplotlib not available")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

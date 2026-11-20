#!/usr/bin/env python3
"""Species-matrix and coverage figures of the multi-clone route on one bench_IW_allele case.

Writes into topology_plots/8_hard_S4_s16003/ (PNG and PDF):

  rck_mc_species_matrix   annotated truth x inferred grid of the per-cell species matrix
                          (aa_eval/aa_species_matrix.tsv, docs/pipeline/evaluation.md §4.7): each cell
                          "J_seg / J_adj", shaded by J_seg (0-1). Rows and columns are ordered so the
                          matched pairs lie on the diagonal (copy-number order would cross ecdna1/ecdna2,
                          both 8.6, against ec3/ec4); matched cells are outlined; the unmatched ec5
                          column is set off by a dashed line. Only the strict (non-low_weight) columns.
  rck_mc_species_coverage one horizontal bar per truth species, laid out in the truth walk's step order
                          (a '-' step runs from its end to its start). Each base is coloured by the inferred
                          species covering it on the same haplotype (any haplotype with --blind), with
                          multiplicity as in §3.1's overlap: the k-th visit of a base by the truth walk is
                          covered by C when C visits that base at least k times, so C's covered length is
                          exactly the §3.1 overlap min(m_T, m_C). A stretch covered by several species is
                          split in height equally among them; the shares written beside the bar are these
                          areas (they sum to 100%). Uncovered bases are light grey.

The inferred walks are the cycles file scored by allele_bench.run_one (done.json "cycles"), read
with rck.ecdna.aa_format / aa_eval.species_cycles; the truth walks come from evaluate.load_truth.
With --block the scores come from aa_eval_block/ and the inferred walks are relabelled per phase
block exactly as aa_eval.score does (evaluate.block_swaps on the blind strict matches).

    python make_species_matrix_figure.py [--blind] [--block] [--outdir DIR]
"""
import argparse
import bisect
import json
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allele_bench as ab  # noqa: E402  (configures bench_iw for iw_intpsi)

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.colors import Normalize  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from rck.ecdna import aa_eval as ae  # noqa: E402
from rck.ecdna import aa_format as aaf  # noqa: E402
from rck.ecdna import evaluate as ev  # noqa: E402

bi = ab.bi
CASE = "hard_S4_s16003"
PLOT_DIR = os.path.join(HERE, "topology_plots", "8_hard_S4_s16003")
METHOD, PHASING = "rck_mc", "native"

# matched pairs (rck_mc_recovery_matrix.png, phased, native labels), on the diagonal
ROWS = ["ecdna3", "ecdna4", "ecdna1", "ecdna2"]
COLS = ["ec1", "ec2", "ec4", "ec3", "ec5"]
MATCHES = {"ecdna3": "ec1", "ecdna4": "ec2", "ecdna1": "ec4", "ecdna2": "ec3"}
UNMATCHED = ["ec5"]

# Paul Tol's "muted" scheme (colourblind-safe), avoiding the haplotype blue/orange of rck.ecdna.plot
SPECIES_COLORS = {"ec1": "#882255", "ec2": "#44AA99", "ec3": "#DDCC77", "ec4": "#332288", "ec5": "#999933"}
UNCOVERED = "#DDDDDD"
INK, INK_MUTED = "#0b0b0b", "#52514e"


def case_info():
    case = next(c for c in bi.cases() if c["name"] == CASE)
    return case, os.path.join(ab.task_dir(case, METHOD, None, PHASING))


def load(block, haplotype):
    case, folder = case_info()
    sub = "aa_eval_block" if block else "aa_eval"
    matrix = pd.read_csv(os.path.join(folder, sub, "aa_species_matrix.tsv"), sep="\t")
    matrix = matrix[matrix["haplotype"] == haplotype]
    with open(os.path.join(folder, "done.json")) as fh:
        cycles_path = json.load(fh)["cycles"]
    truth = {t.name: t for t in ev.load_truth(case["truth"])}
    records = aaf.read_cycles(cycles_path)
    cycles = ae.species_cycles(records)
    if block:
        strict, _rule = ae.strict_select(records, cycles)
        blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
        by = {c.name: c for c in strict}
        pairs = [(truth[t], by[c]) for t, c in ev._match(list(truth.values()), strict, False).items()]
        swaps = ev.block_swaps(pairs, blocks)
        cycles = [ev.relabel(c, blocks, swaps) for c in cycles]
    inferred = {c.name: c for c in cycles}
    missing = [n for n in ROWS if n not in truth] + [n for n in COLS if n not in inferred]
    missing += [n for n in ROWS if n not in set(matrix["truth"])] + [n for n in COLS if n not in set(matrix["inferred"])]
    if missing:
        raise SystemExit("species not found: %s" % sorted(set(missing)))
    return matrix, truth, inferred


def labels(matrix, truth, inferred):
    tw = matrix.drop_duplicates("truth").set_index("truth")["truth_weight"]
    cw = matrix.drop_duplicates("inferred").set_index("inferred")["copy_number"]
    rl = {n: "%s\n%.1f copies · %d kb" % (n, tw[n], round(ae.cycle_length(truth[n].steps) / 1000)) for n in ROWS}
    cl = {n: "%s\n%.1f copies\n%d kb" % (n, cw[n], round(ae.cycle_length(inferred[n].steps) / 1000)) for n in COLS}
    return rl, cl


def _fmt(v):
    return "–" if pd.isna(v) else "%.2f" % v


def figure_matrix(matrix, rl, cl, title, path):
    seg = matrix.pivot(index="truth", columns="inferred", values="j_seg").reindex(index=ROWS, columns=COLS)
    adj = matrix.pivot(index="truth", columns="inferred", values="j_adj").reindex(index=ROWS, columns=COLS)
    fig, ax = plt.subplots(figsize=(8.4, 5.6))
    cmap = plt.get_cmap("Greys")
    norm = Normalize(0, 1)
    im = ax.imshow(seg.values.astype(float), cmap=cmap, norm=norm, aspect="auto")
    for i, r in enumerate(ROWS):
        for j, c in enumerate(COLS):
            v = seg.loc[r, c]
            ax.text(j, i, "%s / %s" % (_fmt(v), _fmt(adj.loc[r, c])), ha="center", va="center", fontsize=10,
                    color="white" if (not pd.isna(v) and v > 0.55) else INK)
            if MATCHES.get(r) == c:
                ax.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False, edgecolor="#d1006b",
                                       linewidth=2.2, zorder=3))
    # gridlines between cells
    for k in range(len(COLS) + 1):
        ax.axvline(k - 0.5, color="white", linewidth=2)
    for k in range(len(ROWS) + 1):
        ax.axhline(k - 0.5, color="white", linewidth=2)
    j0 = COLS.index(UNMATCHED[0])
    ax.axvline(j0 - 0.5, color=INK, linestyle="--", linewidth=1.3, zorder=4)
    ax.text(j0, len(ROWS) - 0.5 + 0.12, "unmatched\n(false species)", ha="center", va="top", fontsize=9,
            color=INK_MUTED)
    # column swatches (species colour, shared with the coverage figure)
    for j, c in enumerate(COLS):
        ax.add_patch(Rectangle((j - 0.45, -0.5 - 0.16), 0.9, 0.12, facecolor=SPECIES_COLORS[c], edgecolor="none",
                               clip_on=False))
    ax.set_xticks(range(len(COLS)))
    ax.set_xticklabels([cl[c] for c in COLS], fontsize=9)
    ax.xaxis.tick_top()
    ax.tick_params(axis="x", length=0, pad=14)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([rl[r] for r in ROWS], fontsize=9)
    ax.tick_params(axis="y", length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    ax.set_xlim(-0.5, len(COLS) - 0.5)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_ylabel("truth species", fontsize=10, color=INK_MUTED)
    cb = fig.colorbar(im, ax=ax, fraction=0.04, pad=0.03)
    cb.set_label("J_seg (cell text: J_seg / J_adj)", fontsize=9)
    cb.outline.set_visible(False)
    fig.suptitle(title, fontsize=11)
    fig.text(0.01, 0.01, "outlined: matched pair · copies are cell averages", fontsize=8, color=INK_MUTED)
    fig.tight_layout(rect=(0, 0.03, 1, 0.97))
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (path, ext), dpi=200)
    plt.close(fig)


def coverage_pieces(t, inferred, with_hap):
    """[(offset, length, covering names)] along the truth walk t, merged where the covering set repeats."""
    key = (lambda s: (s.chrom, s.hap)) if with_hap else (lambda s: (s.chrom,))
    names = [n for n in COLS]
    bps = {}
    for s in list(t.steps) + [s for n in names for s in inferred[n].steps]:
        bps.setdefault(key(s), set()).update((s.start, s.end + 1))
    bps = {k: sorted(v) for k, v in bps.items()}

    def elems(s):
        b = bps[key(s)]
        i, j = bisect.bisect_left(b, s.start), bisect.bisect_left(b, s.end + 1)
        return [(b[x], b[x + 1]) for x in range(i, j)]

    mult = {n: {} for n in names}
    for n in names:
        for s in inferred[n].steps:
            for e in elems(s):
                mult[n][(key(s), e)] = mult[n].get((key(s), e), 0) + 1
    occ, pieces, x = {}, [], 0
    for s in t.steps:
        es = elems(s)
        if s.strand == "-":
            es = es[::-1]
        for e in es:
            k = (key(s), e)
            occ[k] = occ.get(k, 0) + 1
            cov = tuple(n for n in names if mult[n].get(k, 0) >= occ[k])
            ln = e[1] - e[0]
            if pieces and pieces[-1][2] == cov and pieces[-1][0] + pieces[-1][1] == x:
                pieces[-1] = (pieces[-1][0], pieces[-1][1] + ln, cov)
            else:
                pieces.append((x, ln, cov))
            x += ln
    return pieces


def shares(pieces):
    area = {}
    for _x, ln, cov in pieces:
        if not cov:
            area["uncovered"] = area.get("uncovered", 0) + ln
        for n in cov:
            area[n] = area.get(n, 0) + ln / len(cov)
    total = sum(ln for _x, ln, _c in pieces)
    order = [n for n in COLS if n in area] + (["uncovered"] if "uncovered" in area else [])
    return [(n, area[n] / total) for n in order]


def figure_coverage(truth, inferred, rl, with_hap, title, path):
    data = {r: coverage_pieces(truth[r], inferred, with_hap) for r in ROWS}
    xmax = max(sum(ln for _x, ln, _c in p) for p in data.values())
    fig, ax = plt.subplots(figsize=(10.5, 4.6))
    h = 0.62
    texts = {}
    for i, r in enumerate(ROWS):
        for x, ln, cov in data[r]:
            if not cov:
                ax.add_patch(Rectangle((x / 1e3, i - h / 2), ln / 1e3, h, facecolor=UNCOVERED, edgecolor="none"))
                continue
            sub = h / len(cov)
            for k, n in enumerate(cov):
                ax.add_patch(Rectangle((x / 1e3, i - h / 2 + k * sub), ln / 1e3, sub, facecolor=SPECIES_COLORS[n],
                                       edgecolor="none"))
        # step boundaries of the truth walk
        pos = 0
        for s in truth[r].steps[:-1]:
            pos += s.end - s.start + 1
            ax.plot([pos / 1e3] * 2, [i - h / 2, i + h / 2], color="white", linewidth=1.2)
        texts[r] = " · ".join("%s %s" % (n, "<1%" if f < 0.005 else "%d%%" % round(100 * f)) for n, f in shares(data[r]))
        ax.text(xmax / 1e3 * 1.015, i, texts[r], va="center", ha="left", fontsize=9, color=INK)
    ax.set_xlim(0, xmax / 1e3)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([rl[r] for r in ROWS], fontsize=9)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("position along the truth walk (kb, step order; white ticks: step boundaries)", fontsize=9)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    handles = [Patch(facecolor=SPECIES_COLORS[n], label=n) for n in COLS] + [Patch(facecolor=UNCOVERED, label="uncovered")]
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.75, -0.2), ncol=len(handles), frameon=False,
              fontsize=9)
    fig.suptitle(title, fontsize=11)
    fig.subplots_adjust(left=0.13, right=0.64, top=0.88, bottom=0.26)
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (path, ext), dpi=200)
    plt.close(fig)
    return texts


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--blind", action="store_true", help="haplotype-blind rows/coverage (default: phased)")
    p.add_argument("--block", action="store_true", help="scores within phase blocks (aa_eval_block/)")
    p.add_argument("--outdir", default=PLOT_DIR)
    args = p.parse_args()
    hap = "blind" if args.blind else "phased"
    labelling = "within-block labels" if args.block else "native labels"
    suffix = ("_blind" if args.blind else "") + ("_block" if args.block else "")
    matrix, truth, inferred = load(args.block, hap)
    rl, cl = labels(matrix, truth, inferred)
    os.makedirs(args.outdir, exist_ok=True)
    title = "%s — rck_mc, %s, %s" % (CASE, hap, labelling)
    m_path = os.path.join(args.outdir, "rck_mc_species_matrix" + suffix)
    c_path = os.path.join(args.outdir, "rck_mc_species_coverage" + suffix)
    figure_matrix(matrix, rl, cl, title, m_path)
    texts = figure_coverage(truth, inferred, rl, not args.blind, title, c_path)
    for r in ROWS:
        print("%s: %s" % (r, texts[r]))
    print("wrote %s.{png,pdf}\nwrote %s.{png,pdf}" % (m_path, c_path))


if __name__ == "__main__":
    main()

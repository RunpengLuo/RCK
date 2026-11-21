#!/usr/bin/env python3
"""Species-matrix and coverage figures of the multi-clone route on one bench_IW_allele case.

Writes into topology_plots/<n>_<case>/ (PNG and PDF; --plot-dir overrides):

  rck_mc_species_matrix   annotated truth x inferred grid of the per-cell species matrix
                          (aa_eval/aa_species_matrix.tsv, docs/pipeline/evaluation.md §4.7): each cell
                          "J_seg / J_adj", shaded by J_seg (0-1). Rows are the truth species by copy
                          number (descending, ties by name); columns are ordered so each row's match lies
                          on the diagonal, and matched cells are outlined. A truth species with no match
                          gets an empty placeholder column whose diagonal cell reads "missed". Unmatched
                          inferred species (any number, by column rank) come last, after a dashed line.
                          Only the strict (non-low_weight) columns, plus any matched low_weight one.
                          The matches are the run's own one-to-one matching: the "inferred" column of the
                          rck_mc run's recovery_matrix.tsv (rck-ecdna-eval; the pairs drawn in
                          rck_mc_recovery_matrix.png). If that file is missing, the best one-to-one
                          assignment on phased J_seg (pairs with J_seg > 0) is used and a warning printed.
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

Species colours come from a fixed colourblind-safe list (Paul Tol's "muted"), assigned in column order.

    python make_species_matrix_figure.py [--case hard_S4_s16003] [--blind] [--block] [--plot-dir DIR]

Added 2026-11-21 for the RCK MC vs CoRAL topology comparison (make_compare_plots.py):
  --set iw_allelic     the item 1 cases (allele_bench.use_set); default iw_intpsi
  --method coral|ce    a baseline's answer (linear walks kept, --variant), drawn from allele_bench's baf_major
                       task folder; use with --blind (CoRAL / CE have no alleles). Default rck_mc. Files are named
                       <method>_species_matrix* / <method>_species_coverage*.
  --match score        pair species as the scorer does (evaluate.score_decomposition: with --block, blind pairing
                       chooses one A/B swap per phase block, then pairing with haplotypes; with --blind, blind
                       pairing) on the strict species, instead of the run's recovery_matrix.tsv (absolute frame,
                       rck_mc only). Default run (unchanged behaviour).
"""
import argparse
import bisect
import glob
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
DEFAULT_CASE = "hard_S4_s16003"
METHOD, VARIANT, PHASING, MATCH = "rck_mc", None, "native", "run"
LABEL = {"rck_mc": "RCK MC", "coral": "CoRAL", "ce": "CE"}
MISSED = "(missed)"  # placeholder column name for a truth species with no match

# Paul Tol's "muted" scheme (colourblind-safe), avoiding the haplotype blue/orange of rck.ecdna.plot;
# assigned to the inferred species in column order (cycled if there are more than nine)
PALETTE = ["#882255", "#44AA99", "#332288", "#DDCC77", "#999933", "#117733", "#88CCEE", "#CC6677", "#AA4499"]
UNCOVERED = "#DDDDDD"
INK, INK_MUTED = "#0b0b0b", "#52514e"


def case_info(name):
    case = next((c for c in bi.cases() if c["name"] == name), None)
    if case is None:
        raise SystemExit("unknown case %r" % name)
    return case, os.path.join(ab.task_dir(case, METHOD, VARIANT, PHASING))


def default_plot_dir(name):
    found = sorted(glob.glob(os.path.join(HERE, "topology_plots", "*_" + name)))
    if len(found) != 1:
        raise SystemExit("expected one topology_plots/<n>_%s folder, found %d; pass --plot-dir" % (name, len(found)))
    return found[0]


def run_matches(case, matrix_phased):
    """truth -> inferred (None: missed) from the run's own one-to-one matching (recovery_matrix.tsv)."""
    path = os.path.join(bi.run_dir(case, METHOD), "recovery_matrix.tsv")
    if os.path.exists(path):
        rec = pd.read_csv(path, sep="\t")
        return {r: (None if pd.isna(c) else c) for r, c in zip(rec["truth"], rec["inferred"])}, path
    from scipy.optimize import linear_sum_assignment
    seg = matrix_phased.pivot(index="truth", columns="inferred", values="j_seg").fillna(0.0)
    rows, cols = linear_sum_assignment(-seg.values)
    out = {r: None for r in seg.index}
    for i, j in zip(rows, cols):
        if seg.values[i, j] > 0:
            out[seg.index[i]] = seg.columns[j]
    print("WARNING: %s not found; matches are the best one-to-one assignment on phased J_seg" % path)
    return out, None


def layout(matrix, matches):
    """ROWS (truth by copy number), COLS (diagonal matches, then unmatched), UNMATCHED, colours."""
    tw = matrix.drop_duplicates("truth").set_index("truth")["truth_weight"]
    rows = sorted(tw.index, key=lambda n: (-tw[n], n))
    cols = [matches.get(r) or "%s %s" % (MISSED, r) for r in rows]
    inf = matrix.drop_duplicates("inferred").sort_values("column_rank")
    unmatched = [c for c, lw in zip(inf["inferred"], inf["low_weight"]) if not lw and c not in cols]
    cols += unmatched
    species = [c for c in cols if not c.startswith(MISSED)]
    colors = {c: PALETTE[k % len(PALETTE)] for k, c in enumerate(species)}
    return rows, cols, unmatched, colors


def load(name, block, haplotype):
    case, folder = case_info(name)
    sub = "aa_eval_block" if block else "aa_eval"
    full = pd.read_csv(os.path.join(folder, sub, "aa_species_matrix.tsv"), sep="\t")
    matrix = full[full["haplotype"] == haplotype]
    if MATCH == "run":
        matches, source = run_matches(case, full[full["haplotype"] == "phased"])
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
    if MATCH == "score":
        strict_names = {c.name for c in ae.strict_select(records, ae.species_cycles(records))[0]}
        strict = [c for c in cycles if c.name in strict_names]  # relabelled within blocks if --block
        chosen = ev._match(list(truth.values()), strict, haplotype == "phased")
        matches = {t: chosen.get(t) for t in truth}
        source = "scorer pairing (%s)" % ("within phase blocks" if block else haplotype)
    rows, cols, unmatched, colors = layout(matrix, matches)
    species = [c for c in cols if c in colors]
    missing = [n for n in rows if n not in truth] + [n for n in species if n not in inferred]
    missing += [n for n in species if n not in set(matrix["inferred"])]
    missing += [n for n in matches if n not in set(matrix["truth"])]
    if missing:
        raise SystemExit("species not found: %s" % sorted(set(missing)))
    lay = dict(case=case, rows=rows, cols=cols, species=species, unmatched=unmatched, colors=colors,
               matches=matches, match_source=source, cycles=cycles_path, folder=folder)
    return matrix, truth, inferred, lay


def labels(matrix, truth, inferred, lay):
    tw = matrix.drop_duplicates("truth").set_index("truth")["truth_weight"]
    cw = matrix.drop_duplicates("inferred").set_index("inferred")["copy_number"]
    rl = {n: "%s\n%.1f copies · %d kb" % (n, tw[n], round(ae.cycle_length(truth[n].steps) / 1000)) for n in lay["rows"]}
    cl = {n: "%s\n%.1f copies\n%d kb" % (n, cw[n], round(ae.cycle_length(inferred[n].steps) / 1000))
          for n in lay["species"]}
    cl.update({n: "no match" for n in lay["cols"] if n not in cl})
    return rl, cl


def _fmt(v):
    return "–" if pd.isna(v) else "%.2f" % v


def figure_matrix(matrix, rl, cl, lay, title, path):
    ROWS, COLS, UNMATCHED, MATCHES = lay["rows"], lay["cols"], lay["unmatched"], lay["matches"]
    seg = matrix.pivot(index="truth", columns="inferred", values="j_seg").reindex(index=ROWS, columns=COLS)
    adj = matrix.pivot(index="truth", columns="inferred", values="j_adj").reindex(index=ROWS, columns=COLS)
    # 8.4 x 5.6 in for case 8's 5 columns x 4 rows
    fig, ax = plt.subplots(figsize=(max(5.6, 2.6 + 1.16 * len(COLS)), max(3.4, 1.6 + 1.0 * len(ROWS))))
    cmap = plt.get_cmap("Greys")
    norm = Normalize(0, 1)
    im = ax.imshow(seg.values.astype(float), cmap=cmap, norm=norm, aspect="auto")
    for i, r in enumerate(ROWS):
        for j, c in enumerate(COLS):
            if c not in lay["colors"]:  # placeholder column of a missed truth species
                if j == i:
                    ax.text(j, i, "missed", ha="center", va="center", fontsize=10, color=INK_MUTED,
                            style="italic")
                    ax.add_patch(Rectangle((j - 0.47, i - 0.47), 0.94, 0.94, fill=False, edgecolor="#d1006b",
                                           linewidth=1.4, linestyle=":", zorder=3))
                continue
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
    if UNMATCHED:
        j0 = COLS.index(UNMATCHED[0])
        ax.axvline(j0 - 0.5, color=INK, linestyle="--", linewidth=1.3, zorder=4)
        ax.text(j0 + (len(UNMATCHED) - 1) / 2, len(ROWS) - 0.5 + 0.12, "unmatched\n(false species)", ha="center",
                va="top", fontsize=9, color=INK_MUTED)
    # column swatches (species colour, shared with the coverage figure)
    for j, c in enumerate(COLS):
        if c in lay["colors"]:
            ax.add_patch(Rectangle((j - 0.45, -0.5 - 0.16), 0.9, 0.12, facecolor=lay["colors"][c],
                                   edgecolor="none", clip_on=False))
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


def coverage_pieces(t, inferred, names, with_hap):
    """[(offset, length, covering names)] along the truth walk t, merged where the covering set repeats."""
    key = (lambda s: (s.chrom, s.hap)) if with_hap else (lambda s: (s.chrom,))
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


def shares(pieces, names):
    area = {}
    for _x, ln, cov in pieces:
        if not cov:
            area["uncovered"] = area.get("uncovered", 0) + ln
        for n in cov:
            area[n] = area.get(n, 0) + ln / len(cov)
    total = sum(ln for _x, ln, _c in pieces)
    order = [n for n in names if n in area] + (["uncovered"] if "uncovered" in area else [])
    return [(n, area[n] / total) for n in order]


def figure_coverage(truth, inferred, rl, lay, with_hap, title, path):
    ROWS, names, colors = lay["rows"], lay["species"], lay["colors"]
    data = {r: coverage_pieces(truth[r], inferred, names, with_hap) for r in ROWS}
    xmax = max(sum(ln for _x, ln, _c in p) for p in data.values())
    # 10.5 x 4.6 in for case 8's 4 rows; the margins and the legend's offset are kept in inches
    height = max(2.6, 1.8 + 0.7 * len(ROWS))
    top_in, bottom_in, axes_in = 0.12 * 4.6, 0.26 * 4.6, (0.88 - 0.26) * 4.6
    fig, ax = plt.subplots(figsize=(10.5, height))
    h = 0.62
    texts = {}
    for i, r in enumerate(ROWS):
        for x, ln, cov in data[r]:
            if not cov:
                ax.add_patch(Rectangle((x / 1e3, i - h / 2), ln / 1e3, h, facecolor=UNCOVERED, edgecolor="none"))
                continue
            sub = h / len(cov)
            for k, n in enumerate(cov):
                ax.add_patch(Rectangle((x / 1e3, i - h / 2 + k * sub), ln / 1e3, sub, facecolor=colors[n],
                                       edgecolor="none"))
        # step boundaries of the truth walk
        pos = 0
        for s in truth[r].steps[:-1]:
            pos += s.end - s.start + 1
            ax.plot([pos / 1e3] * 2, [i - h / 2, i + h / 2], color="white", linewidth=1.2)
        texts[r] = " · ".join("%s %s" % (n, "<1%" if f < 0.005 else "%d%%" % round(100 * f)) for n, f in shares(data[r], names))
        ax.text(xmax / 1e3 * 1.015, i, texts[r], va="center", ha="left", fontsize=9, color=INK)
    ax.set_xlim(0, xmax / 1e3)
    ax.set_ylim(len(ROWS) - 0.5, -0.5)
    ax.set_yticks(range(len(ROWS)))
    ax.set_yticklabels([rl[r] for r in ROWS], fontsize=9)
    ax.tick_params(axis="y", length=0)
    ax.set_xlabel("position along the truth walk (kb, step order; white ticks: step boundaries)", fontsize=9)
    for side in ("top", "right", "left"):
        ax.spines[side].set_visible(False)
    handles = [Patch(facecolor=colors[n], label=n) for n in names] + [Patch(facecolor=UNCOVERED, label="uncovered")]
    top, bottom = 1 - top_in / height, bottom_in / height
    legend_y = -0.2 if height == 4.6 else -0.2 * axes_in / ((top - bottom) * height)
    ax.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.75, legend_y), ncol=len(handles),
              frameon=False, fontsize=9)
    fig.suptitle(title, fontsize=11)
    fig.subplots_adjust(left=0.13, right=0.64, top=top, bottom=bottom)
    # widen the left margin only if the row labels would be clipped (they fit for case 8)
    fig.canvas.draw()
    x0 = min(t.get_window_extent().x0 for t in ax.get_yticklabels()) / fig.dpi / fig.get_figwidth()
    if x0 < 0:
        fig.subplots_adjust(left=0.13 + 0.01 - x0)
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (path, ext), dpi=200)
    plt.close(fig)
    return texts


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--blind", action="store_true", help="haplotype-blind rows/coverage (default: phased)")
    p.add_argument("--block", action="store_true", help="scores within phase blocks (aa_eval_block/)")
    p.add_argument("--case", default=DEFAULT_CASE, help="bench_IW run name (default: %(default)s)")
    p.add_argument("--plot-dir", "--outdir", dest="plot_dir", default=None,
                   help="output folder (default: topology_plots/<n>_<case>)")
    p.add_argument("--set", choices=["iw_intpsi", "iw_allelic"], default="iw_intpsi")
    p.add_argument("--method", choices=["rck_mc", "coral", "ce"], default="rck_mc")
    p.add_argument("--variant", choices=["walks_kept", "walks_dropped"], default="walks_kept",
                   help="CoRAL / CE scoring variant (default walks_kept)")
    p.add_argument("--match", choices=["run", "score"], default="run")
    args = p.parse_args()
    global METHOD, VARIANT, PHASING, MATCH
    ab.use_set(args.set)
    METHOD, MATCH = args.method, args.match
    if METHOD != "rck_mc":
        VARIANT, PHASING = args.variant, "baf_major"
        if MATCH == "run":
            raise SystemExit("--method %s needs --match score (only rck_mc runs have recovery_matrix.tsv)" % METHOD)
    plot_dir = args.plot_dir or default_plot_dir(args.case)
    hap = "blind" if args.blind else "phased"
    labelling = "within-block labels" if args.block else ("haplotype-blind" if args.blind else "native labels")
    suffix = ("_blind" if args.blind else "") + ("_block" if args.block else "")
    matrix, truth, inferred, lay = load(args.case, args.block, hap)
    rl, cl = labels(matrix, truth, inferred, lay)
    os.makedirs(plot_dir, exist_ok=True)
    if METHOD == "rck_mc" and MATCH == "run":  # the original files keep their original title
        title = "%s — rck_mc, %s, %s" % (args.case, hap, "within-block labels" if args.block else "native labels")
    else:
        title = "%s — %s, %s" % (args.case, LABEL[METHOD], labelling if args.blind else "%s, %s" % (hap, labelling))
    m_path = os.path.join(plot_dir, METHOD + "_species_matrix" + suffix)
    c_path = os.path.join(plot_dir, METHOD + "_species_coverage" + suffix)
    print("matches (%s): %s" % (lay["match_source"] or "fallback: best phased J_seg",
                                ", ".join("%s->%s" % (r, lay["matches"].get(r) or "missed") for r in lay["rows"])))
    print("columns: %s; unmatched: %s" % (", ".join(lay["cols"]), ", ".join(lay["unmatched"]) or "none"))
    figure_matrix(matrix, rl, cl, lay, title, m_path)
    texts = figure_coverage(truth, inferred, rl, lay, not args.blind, title, c_path)
    for r in lay["rows"]:
        print("%s: %s" % (r, texts[r]))
    print("wrote %s.{png,pdf}\nwrote %s.{png,pdf}" % (m_path, c_path))


if __name__ == "__main__":
    main()

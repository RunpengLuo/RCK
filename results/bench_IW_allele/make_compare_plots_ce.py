#!/usr/bin/env python3
"""Truth vs RCK MC vs CoRAL vs CE topology comparison on chosen cases (Yubi, 2026-11-26).

make_compare_plots.topology_compare() with a fourth block: CE's (CycleExtractor 1.0.0) strict cycles, drawn exactly
as CoRAL's are (grey, no alleles; labels = bulk copies, cycle / linear walk, and the haplotype-blind pairing from
CE's own score/aa_eval/aa_species_metrics.tsv). CE as benchmarked: linear walks kept (`ce:walks_kept`, the score/
folder, not score_nolinear/), as CoRAL is drawn (`coral:walks_kept`). Panel order: truth, RCK MC, CoRAL, CE.

topology_compare() is copied here (the CoRAL block generalised to a list of haplotype-blind methods) rather than
refactored in make_compare_plots.py, which is unchanged; its helpers (strict_cycles, scorer_pairs, thin_ticks,
LONG_STEP, CORAL_VARIANT) are imported from it. Steps over LONG_STEP (background in linear walks) do not set the
axis and are drawn clipped to it, per method.

Writes topology_plots/compare_ce/<case>/: topology_compare_ce.png (150 dpi), topology_compare_ce.pdf, README.txt
(inputs used, the benchmark table rows). Run folders are only read; no Gurobi; ~1-2 min per case.

    python make_compare_plots_ce.py                       # the cases in CASES
    python make_compare_plots_ce.py --case SET CASE       # any case with rck_mc, CoRAL and CE runs
"""
import argparse
import os
import shutil
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_compare_plots as mcp  # noqa: E402
import make_topology_plots as mtp  # noqa: E402
import allele_bench as alb  # noqa: E402

OUT = os.path.join(HERE, "topology_plots", "compare_ce")
CASES = [("iw_intpsi", "hard_S2_s15001"), ("iw_allelic", "sib_opp_shared_s18004")]
CE_VARIANT = "walks_kept"
BLIND = [("CoRAL", "coral", mcp.CORAL_VARIANT), ("CE", "ce", CE_VARIANT)]  # haplotype-blind methods, drawn in this order
TITLES = {"truth": "truth (simulator)", "RCK MC": "RCK MC (strict species; haplotypes within phase blocks)",
          "CoRAL": "CoRAL (strict cycles; no alleles)", "CE": "CE (strict cycles; no alleles)"}


def inputs(case):
    """[(what, path)] of every file the figure reads."""
    out = [("truth", case["truth"]),
           ("RCK MC answer", mtp.answer_files(case, "rck_mc")[0]),
           ("RCK MC pairing", os.path.join(alb.bi.score_folder(case, "rck_mc"), "aa_eval", "aa_species_metrics.tsv")),
           ("RCK MC psi", os.path.join(alb.bi.run_dir(case, "rck_mc"), "decompositions.tsv")),
           ("phase blocks", os.path.join(case["rck_input"], "rck.fragments.tsv")),
           ("clone proportions", os.path.join(case["rck_input"], "clone_proportions.tsv"))]
    for label, method, variant in BLIND:
        out.append(("{l} answer ({m}:{v})".format(l=label, m=method, v=variant), alb.answer_files(case, method, variant)[0]))
        out.append(("{l} pairing".format(l=label), os.path.join(alb.bi.score_folder(case, method, variant), "aa_eval",
                                                              "aa_species_metrics.tsv")))
    return out


def topology_compare_ce(case, folder, work):
    """make_compare_plots.topology_compare with CoRAL and CE (BLIND) as haplotype-blind blocks; returns
    {label: {inferred: truth}} as drawn (the scorers' pairings)."""
    import pandas as pd
    from rck.ecdna import bench_metrics as bm, evaluate as ev
    from rck.ecdna import plot as rp

    plt, _, cnplot = rp.import_plotting(mtp.CNPLOT)
    bi = alb.bi
    truth = ev.load_truth(case["truth"])
    psi_t, bulk_t = mtp.clone_psi(case["truth"])
    run = bi.run_dir(case, "rck_mc")
    weights = bm.multiclone_bulk_weights(run, bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv")))
    rck = mcp.strict_cycles(mtp.answer_files(case, "rck_mc")[0], weights)
    blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
    by_t, by_c = {t.name: t for t in truth}, {c.name: c for c in rck}
    swaps = ev.block_swaps([(by_t[a], by_c[b]) for a, b in ev._match(truth, rck, False).items()], blocks)
    rck = [ev.relabel(c, blocks, swaps) for c in rck]
    pairs = {"RCK MC": mcp.scorer_pairs(case, "rck_mc")}  # the table scorer's pairing (within phase blocks)
    check = {v: k for k, v in ev._match(truth, rck, True).items()}
    if check != pairs["RCK MC"]:
        print("{n}: RCK MC pairing: scorer {s} vs re-matched {m}; drawing the scorer's".format(n=case["name"], s=pairs["RCK MC"], m=check))
    dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t")
    dec = dec[dec["decomposition_rank"] == dec["decomposition_rank"].min()].drop_duplicates("species").set_index("species")
    clones = [c[len("psi_"):] for c in dec.columns if c.startswith("psi_")]
    blind = {}
    for label, method, variant in BLIND:
        blind[label] = mcp.strict_cycles(alb.answer_files(case, method, variant)[0])  # the scored (normalised) answer
        pairs[label] = mcp.scorer_pairs(case, method, variant)  # the table scorer's pairing (haplotype-blind)
        check = {v: k for k, v in ev._match(truth, blind[label], False).items()}
        if check != pairs[label]:
            print("{n}: {l} pairing: scorer {s} vs re-matched {m}; drawing the scorer's".format(
                n=case["name"], l=label, s=pairs[label], m=check))

    def steps(c):
        """The walk for drawing: consecutive steps contiguous on the reference (same chromosome, allele and strand)
        merged into one, so 1-bp breakpoint segments do not add rows; the walk itself is unchanged."""
        out = []
        for x in c.steps:
            if out:
                p = out[-1]
                same = p.chrom == x.chrom and p.hap == x.hap and p.strand == x.strand
                if same and x.strand == "+" and x.start == p.end + 1:
                    out[-1] = rp.PStep(p.chrom, p.start, x.end, p.hap, p.strand)
                    continue
                if same and x.strand == "-" and x.end == p.start - 1:
                    out[-1] = rp.PStep(p.chrom, x.start, p.end, p.hap, p.strand)
                    continue
            out.append(rp.PStep(x.chrom, x.start, x.end, x.hap, x.strand))
        return out

    rows = []
    for t in truth:
        base = t.name.split(".")[0]
        rows.append(("truth", t, "{n}\n{b:.1f} bulk copies\npsi {p}".format(n=t.name, b=bulk_t.get(base, t.weight), p=", ".join(
            "{c} {v:g}".format(c=c, v=v) for c, v in sorted(psi_t.get(base, {}).items())))))
    for c in rck:
        psi = ", ".join("{k} {v:g}".format(k=k, v=float(dec.loc[c.name, "psi_" + k])) for k in clones) if c.name in dec.index else ""
        rows.append(("RCK MC", c, "{n} → {p}\n{b:.1f} bulk copies\npsi {s}".format(
            n=c.name, p=pairs["RCK MC"].get(c.name, "no pair"), b=c.weight, s=psi)))
    for label, _, _ in BLIND:
        for c in blind[label]:
            rows.append((label, c, "{n} → {p}\n{b:.1f} copies\n{k}".format(
                n=c.name, p=pairs[label].get(c.name, "no pair"), b=c.weight, k="cycle" if c.circular else "linear walk")))
    # The axis covers every step up to LONG_STEP; longer steps (background in linear walks) are drawn clipped to the
    # axis on their chromosome, and that method's heading says so.
    every = [st for _, c, _ in rows for st in steps(c)]
    short = [st for st in every if st.end - st.start + 1 <= mcp.LONG_STEP]
    covered = {st.chrom for st in short}
    regions = rp.axis_regions(short + [st for st in every if st.chrom not in covered])
    axis = rp.build_axis(cnplot, regions, work)
    span = {}
    for chrom, a, b in regions:
        lo, hi = span.get(chrom, (a, b))
        span[chrom] = (min(lo, a), max(hi, b))
    clipped = {}  # kind -> [names]

    def drawn(kind, c):
        out = []
        for st in steps(c):
            lo, hi = span[st.chrom]
            if st.start < lo or st.end > hi - 1:
                clipped.setdefault(kind, []).append(c.name)
                st = rp.PStep(st.chrom, max(st.start, lo), min(st.end, hi - 1), st.hap, st.strand)
            out.append(st)
        return out

    heights = [max(2.1, 0.2 * len(steps(c))) for _, c, _ in rows]  # 2.1 in: room for the three-line rotated label
    gap = 0.5
    height = min(sum(heights) + 2.5 + 3 * gap, 80)
    fig = plt.figure(figsize=(15, height))
    ratios, kinds = [], []
    for i, (kind, _, _) in enumerate(rows):
        if i and kind != rows[i - 1][0]:
            ratios.append(gap)
            kinds.append(None)
        ratios.append(heights[i])
        kinds.append(i)
    grid = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0.25, top=1.0 - 0.9 / height,
                            bottom=1.8 / height)  # 1.8 in: tick labels and chromosome names clear the legend
    axes, first_of = [], {}
    for g, i in enumerate(kinds):
        if i is None:
            continue
        kind, c, label = rows[i]
        ax = fig.add_subplot(grid[g, 0], sharex=axes[0] if axes else None)
        rp.draw_topology(ax, axis, drawn(kind, c), c.circular, label)
        ax.set_xlim(axis.ch_coords[0], axis.chr_end)
        for boundary in axis.chr_boundaries:
            ax.axvline(boundary, color=rp.INK_MUTED, linewidth=0.8, zorder=0)
        ax.tick_params(axis="x", labelbottom=False, bottom=False)
        if kind not in first_of:
            first_of[kind] = ax
        axes.append(ax)
    for kind, ax in first_of.items():  # headings after drawing, so the clipped note covers the whole block
        title = TITLES[kind]
        if clipped.get(kind):
            title += "; steps over {m} Mb (background in linear walks) clipped to the axis".format(m=mcp.LONG_STEP // 1_000_000)
        ax.set_title(title, loc="left", fontsize=9, color=rp.INK)
    axes[-1].tick_params(axis="x", labelbottom=True, bottom=True)
    rp.position_ticks(axes[-1], axis)
    mcp.thin_ticks(axes[-1], axis)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor=rp.HAP_COLORS["A"], label="haplotype A"),
                        Patch(facecolor=rp.HAP_COLORS["B"], label="haplotype B"),
                        Patch(facecolor="#9a9a9a", label="no allele (CoRAL, CE)")],
               loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8, frameon=False)
    fig.suptitle("{n}: truth vs RCK MC vs CoRAL vs CE, one row per step in walk order (reference-contiguous steps merged; "
                 "→ = the truth species each is paired with)".format(n=case["name"]), fontsize=10, color=rp.INK, x=0.02,
                 y=1.0 - 0.1 / height, ha="left", va="top")
    fig.savefig(os.path.join(folder, "topology_compare_ce.png"), dpi=150)
    fig.savefig(os.path.join(folder, "topology_compare_ce.pdf"))
    plt.close(fig)
    return pairs, {k: sorted(set(v)) for k, v in clipped.items()}


def score_lines(set_name, case):
    """The case's RCK MC, CoRAL and CE rows of the benchmark table (make_compare_plots.score_lines + CE)."""
    import pandas as pd

    table = os.path.join(HERE, "bench_iw_allelic_table.tsv") if set_name == "iw_allelic" else \
        os.path.join(HERE, "..", "bench_IW", "bench_iw_intpsi_table.tsv")
    t = pd.read_csv(table, sep="\t").set_index(["name", "method_variant"])
    cols = ["n_truth", "n_strict", "n_missed", "n_false", "blind_exact_cycles", "phased_exact_cycles", "blind_jseg_best_mean",
            "blind_jadj_best_mean", "cn_recall", "cn_precision", "proportion_tv"]
    lines = ["scores (" + os.path.normpath(os.path.relpath(table, HERE)) + "; phased = within phase blocks, RCK MC only):"]
    lines.append("  " + "  ".join("{:>22}".format(c) if i == 0 else c for i, c in enumerate(["method"] + cols)))
    for mv, label in [("rck_mc", "RCK MC")] + [("{m}:{v}".format(m=m, v=v), l) for l, m, v in BLIND]:
        r = t.loc[(case["name"], mv)]
        lines.append("  {:>22}  ".format(label) + "  ".join("{:.2f}".format(r[c]) if pd.notna(r[c]) else "-" for c in cols))
    return lines


def readme(set_name, case, folder, pairs, clipped):
    lines = ["{n} ({s}): truth vs RCK MC vs CoRAL vs CE topology comparison (make_compare_plots_ce.py, 2026-11-26)".format(
                 n=case["name"], s=set_name), "",
             "topology_compare_ce.png (150 dpi) / .pdf: one genome axis, one row per step in walk order; panels truth, RCK MC "
             "(strict species, haplotypes after the scorer's one A/B swap per phase block), CoRAL (strict cycles, "
             "coral:{cv}, decomposition-only, grey = no alleles), CE (CycleExtractor 1.0.0, strict cycles, ce:{ev}, "
             "decomposition-only, grey). -> = the table scorer's pairing (strict mode of each run's "
             "aa_species_metrics.tsv). Same drawing as make_compare_plots.topology_compare.".format(cv=mcp.CORAL_VARIANT, ev=CE_VARIANT),
             "", "pairings drawn:"]
    for label in ["RCK MC"] + [b[0] for b in BLIND]:
        lines.append("  {l}: {p}".format(l=label, p=", ".join("{a} -> {b}".format(a=a, b=b) for a, b in sorted(pairs[label].items()))
                                          or "none"))
    if clipped:
        lines.append("steps over {m} Mb clipped to the axis: {c}".format(m=mcp.LONG_STEP // 1_000_000, c="; ".join(
            "{k} {v}".format(k=k, v=", ".join(v)) for k, v in clipped.items())))
    lines += ["", "inputs (read only):"] + ["  {w}: {p}".format(w=w, p=p) for w, p in inputs(case)]
    lines += [""] + score_lines(set_name, case)
    with open(os.path.join(folder, "README.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")


def run(set_name, name):
    alb.use_set(set_name)
    case = next(c for c in alb.bi.cases() if c["name"] == name)
    folder = os.path.join(OUT, name)
    os.makedirs(folder, exist_ok=True)
    work = tempfile.mkdtemp(prefix="compare_ce_")
    try:
        pairs, clipped = topology_compare_ce(case, folder, work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    readme(set_name, case, folder, pairs, clipped)
    print("{n}: pairs {p}; clipped {c} -> {f}".format(n=name, p=pairs, c=clipped, f=folder))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--case", nargs=2, metavar=("SET", "CASE"), help="one case (default: CASES)")
    a = p.parse_args()
    import logging
    logging.disable(logging.INFO)
    for set_name, name in ([tuple(a.case)] if a.case else CASES):
        run(set_name, name)
    return 0


if __name__ == "__main__":
    sys.exit(main())

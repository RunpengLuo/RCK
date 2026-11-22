#!/usr/bin/env python3
"""RCK MC vs CoRAL topology comparison on chosen nested, both-allele cases (Yubi, 2026-11-21).

For each case in CASES, writes topology_plots/compare/<k>_<case>/:

  topology_compare.png         one genome axis, one row per species: the truth species, then RCK MC's strict
                               species, then CoRAL's strict cycles. Each row is the walk in step order (steps
                               contiguous on the reference merged for drawing; block per step at its reference position, arrow = strand, connectors in traversal order, dotted
                               = the closing join). Colour = haplotype (A blue, B orange); CoRAL has no alleles, so
                               its steps are grey. RCK MC's haplotypes are shown after the one-A/B-swap-per-phase-
                               block relabelling the scorer uses (A/B are identified only within a phase block), so
                               they are comparable with the truth's. Labels: bulk copies; for truth / RCK MC also psi
                               per clone; the scorer's pairing to a truth species.
  truth_topology.png,          the truth alone (make_topology_plots.plot_truth), and CoRAL's figure of each truth
  coral_truth_<sp>_cycles.png  species
  rck_mc_rank1_ec<k>.png       rck-ecdna-plot pages of the RCK MC run (walk by haplotype, junction arcs, allele copy
                               numbers, ring) WITHIN PHASE BLOCKS (block_copy / block_pages, 2026-11-21): drawn from a
                               copy of the run whose recovery_matrix.tsv / evaluation.tsv are the within-block ones
                               (species_matrix/<set>/<case>/rck_mc_recovery_matrix_block.tsv) and whose walks and
                               inferred karyotypes have the scorer's per-block A/B swaps applied; the header says
                               "within phase blocks" and names the scorer's pair and composite. rck/ecdna/plot.py is
                               unchanged (the header is set by a wrapper here).
  rck_mc_rank1_ec<k>_          the run's own pages (absolute-frame pairing from its evaluation.tsv), kept for comparison
    absolute.png
  rck_mc_recovery_matrix_      the run's recovery matrix page (absolute frame), renamed *_absolute
    absolute.png
  rck_mc_recovery_matrix_      the recovery matrix within phase blocks, paired as the benchmark table's scorer pairs
    block.{png,tsv}            (make_block_figures.py; its mean composite equals the table's on every case)
  coral_rck_mc_cycles.png      CoRAL's own figure style drawing RCK MC's species (haplotype-blind)
  coral_cycles.png             CoRAL's own figure of CoRAL's own answer (its cycles file, all walks)
  coral_graph.png              CoRAL's figure of the graph CoRAL decomposed, trimmed to the region its walks use
                               (plotting copy without the path constraints, which are not drawn)
  rck_mc_species_matrix_block, J_seg / J_adj of every truth x inferred species pair and the coverage of each truth
  rck_mc_species_coverage_block species by the inferred ones (make_species_matrix_figure.py): RCK MC within phase
  rck_mc_species_*_blind,      blocks (alleles must match) and haplotype-blind; CoRAL haplotype-blind. Species are
  coral_species_*_blind        paired as the scorer pairs them (--match score). The *_matrix_* files are two panels,
                               J_seg | J_adj (--layout pair); the earlier one-panel versions are *_onepanel.
  README.txt                   the case, why it was chosen, and its scores (RCK MC vs CoRAL)

CoRAL here is CoRAL 3.1.1 as benchmarked: decomposition-only, linear walks kept (`coral:walks_kept`); on
sib_opp_mirror_s17006 it ran on the zero-read-patched graph (allelic_bench.py), not one of these cases. Run
folders are only read. ~1-2 min per case, so the set runs as a Slurm array (compare_plots.sbatch), one case per
task; no Gurobi.

    python make_compare_plots.py --list
    python make_compare_plots.py --index K          # case K of CASES (the array task id)
    python make_compare_plots.py --readme           # topology_plots/compare/README.md only
    python make_compare_plots.py --index K --pages  # only the RCK MC pages (within phase blocks), checked vs the tables
    python make_compare_plots.py --topology-case SET CASE FOLDER   # topology_compare.png (+ README.txt section) of any
                                                    # case with rck_mc and CoRAL runs (numbered topology_plots/ folders)
    python make_compare_plots.py --pages-case SET CASE FOLDER   # --pages for any case (numbered folders)
    python make_compare_plots.py --topology-index   # that section of topology_plots/README.md

topology_compare labels pairs from the scorer's own aa_species_metrics.tsv (strict; prints a note if re-matching
differs: it did not on any of the 18 cases drawn). Steps over LONG_STEP (CoRAL linear walks' background) do not set
the axis and are drawn clipped to it.
"""
import argparse
import glob
import os
import shutil
import subprocess
import sys
import tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_topology_plots as mtp  # noqa: E402  (imports allele_bench: bench_iw on iw_intpsi)
import allele_bench as alb  # noqa: E402

OUT = os.path.join(HERE, "topology_plots", "compare")
# (set, case, why): the six chosen cases, then the two extras (PLAN.md "Topology comparison")
CASES = [
    ("iw_allelic", "sib_opp_shared_s18004", "outer + nested sibling on the opposite allele, equal clone weights (no "
     "clone-wise allele signal in the input), clean. RCK MC: both species exact. CoRAL: one cycle over both, sibling missed."),
    ("iw_allelic", "sib_opp_mirror_s17004", "outer + nested sibling on the opposite allele, mirrored clones (outer in "
     "clone 1, sibling in clone 2), clean; the pilot case. RCK MC: both exact, psi per clone exact. CoRAL: 1 of 2 exact, "
     "a false species, CN precision 0.35."),
    ("iw_allelic", "sib_mixed_mirror_s19012", "sibling on the opposite allele on some shared intervals and the same "
     "allele on others, mirrored clones, clean. RCK MC: both exact. CoRAL: sibling missed."),
    ("iw_allelic", "sib_opp_shared_s18006", "opposite-allele sibling, equal clone weights, reuse stratum (the outer walk "
     "reuses an extremity). RCK MC: both exact. CoRAL: none exact, CN precision 0.15."),
    ("iw_intpsi", "hard_S4_s16005", "4 nested species, both alleles amplified on 56% of ecDNA bases. RCK MC: 3 of 4 "
     "exact, J_seg 1.00. CoRAL: one cycle, 3 species missed."),
    ("iw_intpsi", "hard_S4_s16002", "4 nested species, both alleles amplified on 39%. RCK MC: 3 of 4 exact, no false "
     "species. CoRAL: 1 missed, none exact, CN precision 0.30."),
    ("iw_intpsi", "hard_S4_s16006", "(extra) 4 nested species, both alleles on 39%. RCK MC: 2 of 4 exact, J_seg 1.00. "
     "CoRAL: 2 missed, none exact."),
    ("iw_intpsi", "hard_S2_s15003", "(extra) 2 nested species, both alleles on 62%. RCK MC: 1 of 2 exact, J_seg 0.82. "
     "CoRAL: J_seg 0.21, none exact."),
]
CORAL_VARIANT = "walks_kept"
LONG_STEP = 5_000_000  # topology_compare: steps longer than this do not set the axis (drawn clipped)


def folder_of(k, name):
    return os.path.join(OUT, "{k}_{n}".format(k=k + 1, n=name))


def strict_cycles(path, weights=None):
    from rck.ecdna import aa_eval as ae, aa_format as aaf, bench_metrics as bm

    records = aaf.read_cycles(path)
    strict = ae.strict_select(records, ae.species_cycles(records))[0]
    return bm.reweight(strict, weights) if weights is not None else strict


def scorer_pairs(case, method, variant=None):
    """{inferred: truth} as the benchmark table's scorer pairs them: score/aa_eval/aa_species_metrics.tsv, strict mode
    (rck_mc within phase blocks, CoRAL haplotype-blind)."""
    import pandas as pd

    sm = pd.read_csv(os.path.join(alb.bi.score_folder(case, method, variant), "aa_eval", "aa_species_metrics.tsv"), sep="\t")
    sm = sm[(sm["mode"] == "strict") & sm["inferred"].notna()]
    return {str(r.inferred): str(r.truth) for r in sm.itertuples()}


def thin_ticks(ax, axis, min_gap=0.025):
    """Drop position tick labels closer than min_gap of the axis width to the previous kept one (narrow collapsed
    regions otherwise stack their labels)."""
    width = float(axis.chr_end - axis.ch_coords[0])
    kept, labels, last = [], [], None
    for x, label in zip(ax.get_xticks(), [t.get_text() for t in ax.get_xticklabels()]):
        if last is None or x - last >= min_gap * width:
            kept.append(x)
            labels.append(label)
            last = x
    ax.set_xticks(kept)
    ax.set_xticklabels(labels)


def topology_compare(case, folder, work):
    """topology_compare.png (module docstring), for any case of the current set (alb.use_set) with rck_mc and CoRAL
    runs; folder is where it goes."""
    import pandas as pd
    from rck.ecdna import bench_metrics as bm, evaluate as ev
    from rck.ecdna import plot as rp

    plt, _, cnplot = rp.import_plotting(mtp.CNPLOT)
    bi = alb.bi
    truth = ev.load_truth(case["truth"])
    psi_t, bulk_t = mtp.clone_psi(case["truth"])
    run = bi.run_dir(case, "rck_mc")
    weights = bm.multiclone_bulk_weights(run, bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv")))
    rck = strict_cycles(mtp.answer_files(case, "rck_mc")[0], weights)
    blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
    by_t, by_c = {t.name: t for t in truth}, {c.name: c for c in rck}
    swaps = ev.block_swaps([(by_t[a], by_c[b]) for a, b in ev._match(truth, rck, False).items()], blocks)
    rck = [ev.relabel(c, blocks, swaps) for c in rck]
    rck_pair = scorer_pairs(case, "rck_mc")  # the table scorer's pairing (within phase blocks)
    check = {v: k for k, v in ev._match(truth, rck, True).items()}
    if check != rck_pair:
        print("{n}: RCK MC pairing: scorer {s} vs re-matched {m}; drawing the scorer's".format(n=case["name"], s=rck_pair, m=check))
    dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t")
    dec = dec[dec["decomposition_rank"] == dec["decomposition_rank"].min()].drop_duplicates("species").set_index("species")
    clones = [c[len("psi_"):] for c in dec.columns if c.startswith("psi_")]
    coral = strict_cycles(alb.answer_files(case, "coral", CORAL_VARIANT)[0])  # the scored (normalised) answer
    coral_pair = scorer_pairs(case, "coral", CORAL_VARIANT)  # the table scorer's pairing (haplotype-blind)
    check = {v: k for k, v in ev._match(truth, coral, False).items()}
    if check != coral_pair:
        print("{n}: CoRAL pairing: scorer {s} vs re-matched {m}; drawing the scorer's".format(n=case["name"], s=coral_pair, m=check))

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
            n=c.name, p=rck_pair.get(c.name, "no pair"), b=c.weight, s=psi)))
    for c in coral:
        rows.append(("CoRAL", c, "{n} → {p}\n{b:.1f} copies\n{k}".format(
            n=c.name, p=coral_pair.get(c.name, "no pair"), b=c.weight, k="cycle" if c.circular else "linear walk")))
    # The axis covers every step up to LONG_STEP. CoRAL's linear walks can hold background steps of tens of Mb, which
    # would stretch the axis until every block is a sliver; those are drawn clipped to the axis on their chromosome
    # (they still cover every drawn region they span), and the CoRAL heading says so.
    every = [st for _, c, _ in rows for st in steps(c)]
    short = [st for st in every if st.end - st.start + 1 <= LONG_STEP]
    covered = {st.chrom for st in short}
    regions = rp.axis_regions(short + [st for st in every if st.chrom not in covered])
    axis = rp.build_axis(cnplot, regions, work)
    span = {}
    for chrom, a, b in regions:
        lo, hi = span.get(chrom, (a, b))
        span[chrom] = (min(lo, a), max(hi, b))
    clipped = []

    def drawn(c):
        out = []
        for st in steps(c):
            lo, hi = span[st.chrom]
            if st.start < lo or st.end > hi - 1:
                clipped.append(c.name)
                st = rp.PStep(st.chrom, max(st.start, lo), min(st.end, hi - 1), st.hap, st.strand)
            out.append(st)
        return out

    heights = [max(2.1, 0.2 * len(steps(c))) for _, c, _ in rows]  # 2.1 in: room for the three-line rotated label
    gap = 0.5
    height = min(sum(heights) + 2.5 + 2 * gap, 60)
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
        rp.draw_topology(ax, axis, drawn(c), c.circular, label)
        ax.set_xlim(axis.ch_coords[0], axis.chr_end)
        for boundary in axis.chr_boundaries:
            ax.axvline(boundary, color=rp.INK_MUTED, linewidth=0.8, zorder=0)
        ax.tick_params(axis="x", labelbottom=False, bottom=False)
        if kind not in first_of:
            first_of[kind] = ax
            ax.set_title({"truth": "truth (simulator)", "RCK MC": "RCK MC (strict species; haplotypes within phase blocks)",
                          "CoRAL": "CoRAL (strict cycles; no alleles)"}[kind], loc="left", fontsize=9, color=rp.INK)
        if kind == "CoRAL" and clipped and "clipped" not in first_of["CoRAL"].get_title(loc="left"):
            first_of["CoRAL"].set_title(first_of["CoRAL"].get_title(loc="left") + "; steps over {m} Mb (background in linear "
                                        "walks) clipped to the axis".format(m=LONG_STEP // 1_000_000), loc="left", fontsize=9, color=rp.INK)
        axes.append(ax)
    axes[-1].tick_params(axis="x", labelbottom=True, bottom=True)
    rp.position_ticks(axes[-1], axis)
    thin_ticks(axes[-1], axis)
    from matplotlib.patches import Patch
    fig.legend(handles=[Patch(facecolor=rp.HAP_COLORS["A"], label="haplotype A"),
                        Patch(facecolor=rp.HAP_COLORS["B"], label="haplotype B"),
                        Patch(facecolor="#9a9a9a", label="no allele (CoRAL)")],
               loc="lower center", bbox_to_anchor=(0.5, 0.0), ncol=3, fontsize=8, frameon=False)
    fig.suptitle("{n}: truth vs RCK MC vs CoRAL, one row per step in walk order (reference-contiguous steps merged; → = the "
                 "truth species each is paired with)".format(n=case["name"]), fontsize=10, color=rp.INK, x=0.02, y=1.0 - 0.1 / height,
                 ha="left", va="top")
    fig.savefig(os.path.join(folder, "topology_compare.png"), dpi=150)
    plt.close(fig)


def _absolute_pages(case, work):
    """rck-ecdna-plot of a copy of the RCK MC run as it is (the run's own, absolute-frame pairing): [png paths]."""
    run = os.path.join(work, "rck_mc")
    shutil.copytree(alb.bi.run_dir(case, "rck_mc"), run, ignore=shutil.ignore_patterns("score", "plots"))
    subprocess.run([sys.executable, "-m", "rck.ecdna.plot", "--outdir", run, "--cnplot-path", mtp.CNPLOT,
                    "--rck-input", case["rck_input"]], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   cwd=mtp.RCK_ROOT)
    return sorted(glob.glob(os.path.join(run, "plots", "*_rank*_*.png")))


def block_copy(set_name, case, work):
    """A copy of the RCK MC run in the benchmark table's frame (within phase blocks), for rck-ecdna-plot. The run folder
    is only read. In the copy:
      recovery_matrix.tsv  <- species_matrix/<set>/<case>/rck_mc_recovery_matrix_block.tsv (make_block_figures.py:
                              the table's pairing and composites)
      evaluation.tsv       <- one row per pair of that TSV (rck-ecdna-plot takes each page's "matched truth" and its
                              scores from evaluation.tsv)
      decompositions.tsv,  <- haplotypes swapped A <-> B on the steps / segments of every block the scorer swaps
      rck.scnt.tsv            (make_block_figures.rank_swaps = evaluate.block_swaps), so walk colours and the A / B
                              copy-number tracks compare with the truth's
    Returns (copy, matrix, blocks, swaps)."""
    import ast
    import pandas as pd
    import make_block_figures as mbf
    from rck.ecdna import evaluate as ev

    run = os.path.join(work, "rck_mc_block")
    shutil.copytree(alb.bi.run_dir(case, "rck_mc"), run, ignore=shutil.ignore_patterns("score", "plots"))
    blocks, swaps, rank = mbf.rank_swaps(case)
    matrix = pd.read_csv(os.path.join(mbf.OUT, set_name, case["name"], "rck_mc_recovery_matrix_block.tsv"), sep="\t")
    if int(matrix["n_blocks_swapped"].iloc[0]) != sum(swaps.values()):
        raise ValueError("{n}: {k} blocks swapped here vs {m} in the block TSV".format(
            n=case["name"], k=sum(swaps.values()), m=int(matrix["n_blocks_swapped"].iloc[0])))
    other = {"A": "B", "B": "A"}

    def swapped(chrom, start, end):
        return bool(swaps.get(ev._block_name(ev.CStep(chrom, start, end, "A", "+"), blocks)))

    def token(t):
        x = ev._parse_step(t)
        return "{c}:{a}-{b}:{h}:{s}".format(c=x.chrom, a=x.start, b=x.end, s=x.strand,
                                             h=other.get(x.hap, x.hap) if swapped(x.chrom, x.start, x.end) else x.hap)

    dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t", dtype=str, keep_default_na=False)
    dec["walk"] = [",".join(token(t) for t in w.split(",")) if w else w for w in dec["walk"]]
    dec.to_csv(os.path.join(run, "decompositions.tsv"), sep="\t", index=False)
    lines = open(os.path.join(run, "rck.scnt.tsv")).read().splitlines()
    out = [lines[0]]
    for line in lines[1:]:
        f = line.split("\t")
        if len(f) > 3 and swapped(f[0], int(f[1]), int(f[2])):
            parts = []
            for part in f[3].split(";"):
                if part.startswith("cn="):
                    cn = ast.literal_eval(part[3:])
                    part = "cn=" + repr({k: {other.get(h, h): v for h, v in a.items()} for k, a in cn.items()})
                parts.append(part)
            f[3] = ";".join(parts)
        out.append("\t".join(f))
    with open(os.path.join(run, "rck.scnt.tsv"), "w") as handle:
        handle.write("\n".join(out) + "\n")
    matrix.to_csv(os.path.join(run, "recovery_matrix.tsv"), sep="\t", index=False)
    paired = matrix[matrix["inferred"].notna()]
    pd.DataFrame({"truth": paired["truth"], "inferred": paired["inferred"], "decomposition_rank": rank,
                  "jaccard_hap": paired["haplotype_jaccard"], "jaccard_total": paired["segment_jaccard"],
                  "junction_f1": paired["junction_f1"], "phasing_accuracy": paired["phasing_accuracy"],
                  "order_match": paired["order_match"], "truth_weight": paired["truth_weight"],
                  "inferred_weight": paired["inferred_weight"], "cn_rel_err": paired["cn_rel_err"],
                  "composite": paired["composite"]}).to_csv(os.path.join(run, "evaluation.tsv"), sep="\t", index=False)
    return run, matrix, blocks, swaps


def block_pages(set_name, case, work, table_mean):
    """rck-ecdna-plot (rck.ecdna.plot.plot_run) of block_copy's copy, each species page retitled: "within phase blocks",
    the paired truth species and the scorer's composite. The page layout is rck-ecdna-plot's, unchanged. Returns
    ([png paths], {inferred: (truth, composite)} drawn)."""
    from rck.ecdna import plot as rp

    run, matrix, blocks, swaps = block_copy(set_name, case, work)
    drawn = {}
    k = sum(swaps.values())
    note = ("within phase blocks: {k} of {n} phase blocks A/B-swapped, haplotypes relabelled so they compare with the truth's"
            if k else "within phase blocks: 0 of {n} phase blocks A/B-swapped (haplotypes as in the run)").format(k=k, n=len(blocks))
    original = rp.species_page

    def titled(plt, cnplot, species, segments, name, scratch, truth=None, score=None, pieces=None):
        fig = original(plt, cnplot, species, segments, name, scratch, truth=truth, score=score, pieces=pieces)
        head = fig._suptitle.get_text().split("\n")[0]
        if truth and score:
            drawn[species.name] = (truth[0], float(score["composite"]))
            lines = [head, "{w}; paired as the benchmark table's scorer pairs: truth {t} ({c:.2f} copies)".format(
                w=note, t=truth[0], c=truth[1]),
                "scorer's composite {c:.3f} (case mean {m:.3f} = table), segment Jaccard {j:.3f}, junction F1 {f:.3f}, "
                "phasing accuracy {p:.3f}, order {o}".format(c=score["composite"], m=table_mean, j=score["jaccard_total"],
                                                             f=score["junction_f1"], p=score["phasing_accuracy"],
                                                             o="matched" if str(score["order_match"]) == "True" else "not matched")]
        else:
            drawn[species.name] = (None, None)
            lines = [head, note + "; not paired with a truth species by the benchmark table's scorer"]
        fig._suptitle.set_text("\n".join(lines))
        return fig

    rp.species_page = titled
    try:
        written = rp.plot_run(run, cnplot_path=mtp.CNPLOT, rck_input=case["rck_input"])
    finally:
        rp.species_page = original
    return [p for p in written.get("png", []) if "_rank" in os.path.basename(p)], drawn


def page_name(png):
    return "rck_mc_" + os.path.basename(png).replace("union_", "")


def rck_mc_pages(set_name, case, folder, work):
    """rck_mc_rank1_ec<k>.png within phase blocks (block_pages); the run's own absolute-frame pages are kept as
    rck_mc_rank1_ec<k>_absolute.png (renamed if present, else drawn). Returns {inferred: (truth, composite)}."""
    import pandas as pd
    import make_block_figures as mbf

    table = pd.read_csv(mbf.TABLES[set_name], sep="\t")
    mean = float(table[(table["name"] == case["name"]) & (table["method_variant"] == "rck_mc")]["mean_composite"].iloc[0])
    existing = sorted(glob.glob(os.path.join(folder, "rck_mc_rank*_*.png")))
    existing = [p for p in existing if not p.endswith("_absolute.png")]
    if existing and not glob.glob(os.path.join(folder, "rck_mc_rank*_absolute.png")):
        for p in existing:
            os.replace(p, p[:-len(".png")] + "_absolute.png")
    elif not existing:
        for png in _absolute_pages(case, work):
            shutil.copy(png, os.path.join(folder, page_name(png)[:-len(".png")] + "_absolute.png"))
    pngs, drawn = block_pages(set_name, case, work, mean)
    for png in pngs:
        shutil.copy(png, os.path.join(folder, page_name(png)))
    return drawn


def coral_figures(case, folder, work):
    bi = alb.bi
    species = os.path.join(work, "rck_mc_species_cycles.txt")
    mtp.species_only(mtp.answer_files(case, "rck_mc")[0], species)
    mtp.coral_plot(os.path.join(work, "coral_rck_mc"), cycles=species)
    shutil.copy(os.path.join(work, "coral_rck_mc_cycles.png"), folder)
    cycles = bi.cycles_file(case, "coral")  # CoRAL's own output, every walk
    mtp.coral_plot(os.path.join(work, "coral"), cycles=cycles)
    shutil.copy(os.path.join(work, "coral_cycles.png"), folder)
    small = os.path.join(work, "coral_graph_trimmed.txt")
    mtp.trimmed_graph(bi.inferred_graph(case), cycles, small)
    # the long-read path constraints index sequence edges by position (e34 ...), which no longer exist after
    # trimming (coral plot then fails with IndexError); they are not drawn, so the plotting copy drops them
    with open(small) as handle:
        kept = [l for l in handle if not l.startswith(("path_constraint", "PathConstraint"))]
    with open(small, "w") as handle:
        handle.writelines(kept)
    mtp.coral_plot(os.path.join(work, "coral"), graph=small)
    shutil.copy(os.path.join(work, "coral_graph.png"), folder)


def matrices(set_name, case, folder):
    """Species matrices and coverage. One-panel matrices go to *_onepanel (kept for comparison); the matrices under
    the plain names are the two-panel J_seg | J_adj figures (--layout pair, 2026-11-21)."""
    base = [sys.executable, os.path.join(HERE, "make_species_matrix_figure.py"), "--set", set_name, "--case", case["name"],
            "--match", "score"]
    runs = (["--method", "rck_mc", "--block"], ["--method", "rck_mc", "--blind"],
            ["--method", "coral", "--blind", "--variant", CORAL_VARIANT])
    tmp = tempfile.mkdtemp(prefix="matrices_")
    try:
        with open(os.path.join(folder, "species_matrix.log"), "a") as log:
            for extra in runs:
                subprocess.run(base + extra + ["--plot-dir", tmp], check=True, stdout=log, stderr=subprocess.STDOUT, cwd=HERE)
                subprocess.run(base + extra + ["--plot-dir", folder, "--layout", "pair"], check=True, stdout=log,
                               stderr=subprocess.STDOUT, cwd=HERE)
        for f in os.listdir(tmp):
            stem, ext = os.path.splitext(f)
            shutil.copy(os.path.join(tmp, f), os.path.join(folder, stem + "_onepanel" + ext if "_matrix_" in f else f))
    finally:
        shutil.rmtree(tmp, ignore_errors=True)


def recovery_block(set_name, case, folder):
    """rck_mc_recovery_matrix_block.{tsv,png}: the recovery matrix within phase blocks (make_block_figures.py); the
    run's own (absolute-frame) page is kept as rck_mc_recovery_matrix_absolute.png."""
    import pandas as pd
    import make_block_figures as mbf

    absolute = os.path.join(folder, "rck_mc_recovery_matrix.png")
    if os.path.exists(absolute):
        os.replace(absolute, os.path.join(folder, "rck_mc_recovery_matrix_absolute.png"))
    table = pd.read_csv(mbf.TABLES[set_name], sep="\t")
    mean = float(table[(table["name"] == case["name"]) & (table["method_variant"] == "rck_mc")]["mean_composite"].iloc[0])
    matrix = mbf.recovery_block(case)
    matrix.to_csv(os.path.join(folder, "rck_mc_recovery_matrix_block.tsv"), sep="\t", index=False)
    mbf.draw_recovery(matrix, case, mean, os.path.join(folder, "rck_mc_recovery_matrix_block.png"))


def score_lines(set_name, case):
    """The case's RCK MC and CoRAL rows of the benchmark table, as README.txt lines."""
    import pandas as pd

    table = os.path.join(HERE, "bench_iw_allelic_table.tsv") if set_name == "iw_allelic" else \
        os.path.join(HERE, "..", "bench_IW", "bench_iw_intpsi_table.tsv")
    t = pd.read_csv(table, sep="\t").set_index(["name", "method_variant"])
    cols = ["n_truth", "n_strict", "n_missed", "n_false", "blind_exact_cycles", "phased_exact_cycles", "blind_jseg_best_mean",
            "blind_jadj_best_mean", "cn_recall", "cn_precision", "proportion_tv"]
    lines = ["scores (" + os.path.relpath(table, HERE) + "; phased = within phase blocks, RCK MC only):"]
    lines.append("  " + "  ".join("{:>22}".format(c) if i == 0 else c for i, c in enumerate(["method"] + cols)))
    for mv, label in (("rck_mc", "RCK MC"), ("coral:" + CORAL_VARIANT, "CoRAL")):
        r = t.loc[(case["name"], mv)]
        lines.append("  {:>22}  ".format(label) + "  ".join("{:.2f}".format(r[c]) if pd.notna(r[c]) else "-" for c in cols))
    return lines


TOPO_MARK = "--- topology_compare.png (make_compare_plots.py --topology-case) ---"


def topology_readme(set_name, case, folder):
    """Add (or replace) the topology_compare.png section of folder/README.txt: what it is and the case's RCK MC vs CoRAL
    scores from the benchmark table."""
    path = os.path.join(folder, "README.txt")
    text = open(path).read() if os.path.exists(path) else ""
    text = text.split("\n" + TOPO_MARK)[0].rstrip("\n")
    lines = ["", TOPO_MARK, "topology_compare.png: one genome axis; truth species, RCK MC strict species (haplotypes after the "
             "scorer's one A/B swap per phase block), CoRAL strict cycles (coral:" + CORAL_VARIANT + ", decomposition-only); "
             "-> = the table scorer's pairing.", ""] + score_lines(set_name, case)
    if not text:
        text = "{n} ({s}); truth run {t}".format(n=case["name"], s=set_name, t=case["truth"])
    with open(path, "w") as handle:
        handle.write(text + "\n" + "\n".join(lines) + "\n")


def readme(set_name, case, why, folder):
    import pandas as pd

    lines = ["{n} ({s}): {w}".format(n=case["name"], s=set_name, w=why), "truth run: " + case["truth"], ""]
    lines += score_lines(set_name, case)
    if set_name == "iw_allelic":
        rel = pd.read_csv(os.path.join(HERE, "sibling_relation.tsv"), sep="\t").set_index("case").loc[case["name"]]
        lines += ["", "RCK MC sibling-vs-outer relation (sibling_relation.tsv): opposite {o} bp, same {s} bp in the truth; "
                  "{d} bp disagreeing in the answer".format(o=rel.truth_opp_bp, s=rel.truth_same_bp, d=rel.disagree_bp)]
    with open(os.path.join(folder, "README.txt"), "w") as handle:
        handle.write("\n".join(lines) + "\n")


def write_index():
    os.makedirs(OUT, exist_ok=True)
    text = ["# RCK MC vs CoRAL topology comparison (make_compare_plots.py)", "",
            "Cases chosen 2026-11-21 (Yubi: the six plus the two extras) to show RCK MC on nested ecDNA with both "
            "alleles amplified. make_compare_plots.py's docstring says what each file is; each folder's README.txt has "
            "the case's scores. Haplotype A blue, B orange; CoRAL has no alleles (grey / CoRAL's own style).", "",
            "**RCK MC pages are within phase blocks (2026-11-21).** `rck_mc_rank1_ec<k>.png` pairs each inferred species with "
            "the truth species the benchmark table's scorer pairs it with, haplotypes relabelled with the scorer's one A/B swap "
            "per phase block (walk colours and A/B copy-number tracks); the header says \"within phase blocks\" and gives the "
            "scorer's composite and the case mean (= the table's `mean_composite`). Checked on all 8 cases against "
            "`rck_mc_recovery_matrix_block.tsv`, the scorer's per-species rows and the table. The run's own absolute-frame "
            "pages are kept as `rck_mc_rank1_ec<k>_absolute.png`. Pairs that changed: s18004 ec1 ecdna1 → ecdna2, ec2 ecdna2 "
            "→ ecdna1; s18006 the same; hard_S4_s16005 ec3 ecdna1 → ecdna4, ec4 ecdna4 → ecdna1; hard_S4_s16002 ec4 (none → "
            "ecdna4); hard_S4_s16006 ec1 ecdna4 → ecdna3, ec2 (none → ecdna1), ec3 (none → ecdna2), ec4 ecdna3 → ecdna4 "
            "(`../rck_mc_pages_frame_check.tsv`). Redo: `python make_compare_plots.py --index K --pages`.", "",
            "| folder | set | why |", "|---|---|---|"]
    for k, (set_name, name, why) in enumerate(CASES):
        text.append("| `{f}` | {s} | {w} |".format(f=os.path.basename(folder_of(k, name)), s=set_name, w=why))
    tmp = os.path.join(OUT, ".README.md.{p}".format(p=os.getpid()))
    with open(tmp, "w") as handle:
        handle.write("\n".join(text) + "\n")
    os.replace(tmp, os.path.join(OUT, "README.md"))


def run_case(k):
    set_name, name, why = CASES[k]
    alb.use_set(set_name)
    case = next(c for c in alb.bi.cases() if c["name"] == name)
    folder = folder_of(k, name)
    os.makedirs(folder, exist_ok=True)
    work = tempfile.mkdtemp(prefix="compare_")
    steps = [("readme", lambda: readme(set_name, case, why, folder)),
             ("truth", lambda: mtp.plot_truth(case["truth"], folder)),
             ("topology_compare", lambda: topology_compare(case, folder, work)),
             ("rck_mc pages", lambda: rck_mc_pages(set_name, case, folder, work)),
             ("coral figures", lambda: coral_figures(case, folder, work)),
             ("species matrices", lambda: matrices(set_name, case, folder)),
             ("recovery matrix within blocks", lambda: recovery_block(set_name, case, folder))]
    failed = []
    try:
        for label, step in steps:
            try:
                step()
                print("{n}: {l} ok".format(n=name, l=label), flush=True)
            except Exception as error:  # noqa: BLE001 - the other figures still get drawn; reported at the end
                import traceback
                traceback.print_exc()
                failed.append(label)
                print("{n}: {l} FAILED: {e!r}".format(n=name, l=label, e=error), flush=True)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    write_index()
    print("{n}: done -> {f}{x}".format(n=name, f=folder, x="; FAILED: " + ", ".join(failed) if failed else ""), flush=True)
    return 1 if failed else 0


INDEX_MARK = "## topology_compare.png in the numbered folders (make_compare_plots.py --topology-index)"
EXTRA_TOPOLOGY = [("iw_allelic", "sib_opp_shared_s18002", "item1_dryrun/sib_opp_shared_s18002"),
                  ("iw_allelic", "sib_mixed_mirror_s19006", "item1_dryrun/sib_mixed_mirror_s19006")]


def topology_index():
    """Add (or replace) the topology_compare.png section of topology_plots/README.md: one line per figure with the case's
    RCK MC vs CoRAL scores from the benchmark tables, and the RCK MC page frame check (audit_page_pairing.py)."""
    import pandas as pd

    tables = {"iw_intpsi": os.path.join(HERE, "..", "bench_IW", "bench_iw_intpsi_table.tsv"),
              "iw_allelic": os.path.join(HERE, "bench_iw_allelic_table.tsv")}
    jobs = [("iw_intpsi", n, "{k}_{n}".format(k=k, n=n)) for k, (n, _d) in enumerate(mtp.CASES, start=1)] + EXTRA_TOPOLOGY
    lines = [INDEX_MARK, "",
             "Added 2026-11-21. `topology_compare.png` as in `compare/`: one genome axis; the truth species, RCK MC's strict "
             "species (haplotypes after the scorer's one A/B swap per phase block) and CoRAL's strict cycles (`coral:walks_kept`, "
             "decomposition-only); \"→ ecdnaN\" is the benchmark scorer's pairing (`score/aa_eval/aa_species_metrics.tsv`, "
             "strict). CoRAL steps over 5 Mb (background in its linear walks) do not set the axis and are drawn clipped to it "
             "(said in the CoRAL heading when it happens: s15001). Each folder's README.txt has the full score rows. RCK MC "
             "exact = phased exact cycles (within phase blocks); CoRAL exact = blind. Scores: `bench_IW/bench_iw_intpsi_table.tsv`, "
             "`bench_iw_allelic_table.tsv`.", "",
             "| folder | RCK MC | CoRAL |", "|---|---|---|"]
    for set_name, name, folder in jobs:
        t = pd.read_csv(tables[set_name], sep="\t").set_index(["name", "method_variant"])
        r, c = t.loc[(name, "rck_mc")], t.loc[(name, "coral:" + CORAL_VARIANT)]
        lines.append("| `{f}` | {re:.0f} / {n:.0f} exact, J_seg {rs:.2f}, J_adj {ra:.2f}, {rm:.0f} missed, {rf:.0f} false, mean composite "
                     "{rc:.3f} | {ce:.0f} / {n:.0f} exact, J_seg {cs:.2f}, J_adj {ca:.2f}, {cm:.0f} missed, {cf:.0f} false, CN precision "
                     "{cp:.2f} |".format(f=folder, n=r.n_truth, re=r.phased_exact_cycles, rs=r.blind_jseg_best_mean,
                                        ra=r.blind_jadj_best_mean, rm=r.n_missed, rf=r.n_false, rc=r.mean_composite,
                                        ce=c.blind_exact_cycles, cs=c.blind_jseg_best_mean, ca=c.blind_jadj_best_mean,
                                        cm=c.n_missed, cf=c.n_false, cp=c.cn_precision))
    lines += ["", "`item1_dryrun/`: the three dry-run candidates' truth figures; all three are finished iw_allelic cases (the "
              "same simulations), so s18002 and s19006 got a `topology_compare.png` there; s17004's is `compare/2_sib_opp_mirror_s17004/`. "
              "`species_matrix/`: matrix figures only (no topology figure needed).", "",
              "**RCK MC pages here are within phase blocks too (2026-11-21, `make_compare_plots.py --pages-case iw_intpsi "
              "<case> <folder>`):** `rck_mc_rank1_ec<k>.png` pairs and colours as the benchmark table (checked on all 8 against "
              "the block TSV, the scorer's rows and `mean_composite`); the run's own pages are `*_absolute.png`. "
              "`rck_mc_pages_frame_check.tsv` (audit_page_pairing.py) lists old vs new pair per page. Pair changed on "
              "`4_simple_S4_s13001` ec3 (none → ecdna3, 0.733) and `7_hard_S4_s16005` ec3 (ecdna1 → ecdna4) and ec4 "
              "(ecdna4 → ecdna1). Same pair, higher composite: s12003 ec2, s13001 ec1/ec2/ec4, s16005 ec1/ec2. Unchanged: "
              "s11001, s14003, s15001, s15003, s16003. `rck_mc_recovery_matrix.png` here is still the run's (absolute); the "
              "within-block matrix is in `species_matrix/iw_intpsi/<case>/`."]
    path = os.path.join(mtp.OUT, "README.md")
    text = open(path).read().split("\n" + INDEX_MARK)[0].rstrip("\n")
    with open(path, "w") as handle:
        handle.write(text + "\n\n" + "\n".join(lines) + "\n")


def check_pages(set_name, case, drawn):
    """[problems]: each drawn page's pair and composite vs rck_mc_recovery_matrix_block.tsv, the table scorer's
    per-species rows and the table's mean_composite."""
    import pandas as pd
    import make_block_figures as mbf

    m = pd.read_csv(os.path.join(mbf.OUT, set_name, case["name"], "rck_mc_recovery_matrix_block.tsv"), sep="\t")
    want = {r.inferred: (r.truth, float(r.composite), r.table_inferred, float(r.table_composite))
            for r in m.itertuples() if isinstance(r.inferred, str)}
    table = pd.read_csv(mbf.TABLES[set_name], sep="\t")
    mean = float(table[(table["name"] == case["name"]) & (table["method_variant"] == "rck_mc")]["mean_composite"].iloc[0])
    problems = []
    for sp, (t, c) in drawn.items():
        if sp not in want:
            if t is not None:
                problems.append("{s}: drawn paired with {t}, not paired in the block TSV".format(s=sp, t=t))
            continue
        wt, wc, tt, tc = want[sp]
        if t != wt or sp != tt or abs(c - wc) > 1e-9 or abs(c - tc) > 1e-6:
            problems.append("{s}: drawn {t} {c} vs block TSV {wt} {wc}; table scorer pairs {wt} with {tt} ({tc})".format(
                s=sp, t=t, c=c, wt=wt, wc=wc, tt=tt, tc=tc))
    for sp in set(want) - set(drawn):
        problems.append("{s}: paired with {t} in the block TSV but has no page".format(s=sp, t=want[sp][0]))
    if abs(float(m["composite"].mean()) - mean) > 1e-6:
        problems.append("mean composite {a} vs table {b}".format(a=float(m["composite"].mean()), b=mean))
    return problems


def pages_only(k):
    set_name, name, _ = CASES[k]
    return pages_case(set_name, name, folder_of(k, name))


def pages_case(set_name, name, folder):
    """rck_mc_pages for any case into folder (compare or numbered topology_plots/ folders), checked vs the tables."""
    alb.use_set(set_name)
    case = next(c for c in alb.bi.cases() if c["name"] == name)
    work = tempfile.mkdtemp(prefix="compare_")
    try:
        drawn = rck_mc_pages(set_name, case, folder, work)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    problems = check_pages(set_name, case, drawn)
    for sp, (t, c) in sorted(drawn.items()):
        print("{n}\t{s}\t{t}\t{c}".format(n=name, s=sp, t=t or "-", c="-" if c is None else "{:.3f}".format(c)))
    print("{n}: pages {r}".format(n=name, r="agree with the block TSV and the table" if not problems else
                                  "DISAGREE: " + "; ".join(problems)), flush=True)
    return 1 if problems else 0


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--index", type=int, help="case index in CASES (0-based; the Slurm array task id)")
    p.add_argument("--list", action="store_true")
    p.add_argument("--readme", action="store_true")
    p.add_argument("--coral-only", action="store_true", help="with --index: redo only the CoRAL figures")
    p.add_argument("--pages", action="store_true", help="with --index: redo only the RCK MC pages (within phase blocks; "
                   "the run's own pages kept as *_absolute.png) and check their pairing against the tables")
    p.add_argument("--topology-case", nargs=3, metavar=("SET", "CASE", "FOLDER"),
                   help="draw only topology_compare.png of any case with rck_mc and CoRAL runs into FOLDER, and its "
                   "section in FOLDER/README.txt (scores from the benchmark table)")
    p.add_argument("--pages-case", nargs=3, metavar=("SET", "CASE", "FOLDER"),
                   help="--pages for any case into FOLDER (e.g. the numbered topology_plots/ folders)")
    p.add_argument("--topology-index", action="store_true", help="write the topology_compare.png section of "
                   "topology_plots/README.md (numbered folders; scores from the tables)")
    a = p.parse_args()
    import logging
    logging.disable(logging.INFO)
    if a.topology_index:
        topology_index()
        return 0
    if a.topology_case:
        set_name, name, folder = a.topology_case
        alb.use_set(set_name)
        case = next(c for c in alb.bi.cases() if c["name"] == name)
        os.makedirs(folder, exist_ok=True)
        work = tempfile.mkdtemp(prefix="compare_")
        try:
            topology_compare(case, folder, work)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        topology_readme(set_name, case, folder)
        print("{n}: topology_compare.png and README.txt section -> {f}".format(n=name, f=folder))
        return 0
    if a.pages_case:
        return pages_case(*a.pages_case)
    if a.pages:
        if a.index is None:
            p.error("--pages needs --index")
        return pages_only(a.index)
    if a.list:
        for k, (s, n, _) in enumerate(CASES):
            print(k, s, n, sep="\t")
        return 0
    if a.readme:
        write_index()
        return 0
    if a.index is None:
        p.error("--index, --list or --readme")
    if a.coral_only:
        set_name, name, _ = CASES[a.index]
        alb.use_set(set_name)
        case = next(c for c in alb.bi.cases() if c["name"] == name)
        work = tempfile.mkdtemp(prefix="compare_")
        try:
            coral_figures(case, folder_of(a.index, name), work)
        finally:
            shutil.rmtree(work, ignore_errors=True)
        print("{n}: coral figures ok".format(n=name))
        return 0
    return run_case(a.index)


if __name__ == "__main__":
    sys.exit(main())

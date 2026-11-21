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
  rck_mc_rank1_ec<k>.png,      rck-ecdna-plot of the RCK MC run (walk by haplotype, junction arcs, allele copy
  rck_mc_recovery_matrix.png   numbers, ring; make_topology_plots.plot_case's RCK MC part)
  coral_rck_mc_cycles.png      CoRAL's own figure style drawing RCK MC's species (haplotype-blind)
  coral_cycles.png             CoRAL's own figure of CoRAL's own answer (its cycles file, all walks)
  coral_graph.png              CoRAL's figure of the graph CoRAL decomposed, trimmed to the region its walks use
                               (plotting copy without the path constraints, which are not drawn)
  rck_mc_species_matrix_block, J_seg / J_adj of every truth x inferred species pair and the coverage of each truth
  rck_mc_species_coverage_block species by the inferred ones (make_species_matrix_figure.py): RCK MC within phase
  rck_mc_species_*_blind,      blocks (alleles must match) and haplotype-blind; CoRAL haplotype-blind. Species are
  coral_species_*_blind        paired as the scorer pairs them (--match score).
  README.txt                   the case, why it was chosen, and its scores (RCK MC vs CoRAL)

CoRAL here is CoRAL 3.1.1 as benchmarked: decomposition-only, linear walks kept (`coral:walks_kept`); on
sib_opp_mirror_s17006 it ran on the zero-read-patched graph (allelic_bench.py), not one of these cases. Run
folders are only read. ~1-2 min per case, so the set runs as a Slurm array (compare_plots.sbatch), one case per
task; no Gurobi.

    python make_compare_plots.py --list
    python make_compare_plots.py --index K          # case K of CASES (the array task id)
    python make_compare_plots.py --readme           # topology_plots/compare/README.md only
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


def folder_of(k, name):
    return os.path.join(OUT, "{k}_{n}".format(k=k + 1, n=name))


def strict_cycles(path, weights=None):
    from rck.ecdna import aa_eval as ae, aa_format as aaf, bench_metrics as bm

    records = aaf.read_cycles(path)
    strict = ae.strict_select(records, ae.species_cycles(records))[0]
    return bm.reweight(strict, weights) if weights is not None else strict


def topology_compare(case, folder, work):
    """topology_compare.png (module docstring)."""
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
    rck_pair = {v: k for k, v in ev._match(truth, rck, True).items()}
    dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t")
    dec = dec[dec["decomposition_rank"] == dec["decomposition_rank"].min()].drop_duplicates("species").set_index("species")
    clones = [c[len("psi_"):] for c in dec.columns if c.startswith("psi_")]
    coral = strict_cycles(alb.answer_files(case, "coral", CORAL_VARIANT)[0])  # the scored (normalised) answer
    coral_pair = {v: k for k, v in ev._match(truth, coral, False).items()}

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
    regions = rp.axis_regions([st for _, c, _ in rows for st in steps(c)])
    axis = rp.build_axis(cnplot, regions, work)
    heights = [max(1.3, 0.2 * len(steps(c))) for _, c, _ in rows]
    gap = 0.5
    height = min(sum(heights) + 1.8 + 2 * gap, 60)
    fig = plt.figure(figsize=(15, height))
    ratios, kinds = [], []
    for i, (kind, _, _) in enumerate(rows):
        if i and kind != rows[i - 1][0]:
            ratios.append(gap)
            kinds.append(None)
        ratios.append(heights[i])
        kinds.append(i)
    grid = fig.add_gridspec(len(ratios), 1, height_ratios=ratios, hspace=0.25, top=1.0 - 0.9 / height,
                            bottom=1.1 / height)
    axes, first_of = [], {}
    for g, i in enumerate(kinds):
        if i is None:
            continue
        kind, c, label = rows[i]
        ax = fig.add_subplot(grid[g, 0], sharex=axes[0] if axes else None)
        rp.draw_topology(ax, axis, steps(c), c.circular, label)
        ax.set_xlim(axis.ch_coords[0], axis.chr_end)
        for boundary in axis.chr_boundaries:
            ax.axvline(boundary, color=rp.INK_MUTED, linewidth=0.8, zorder=0)
        ax.tick_params(axis="x", labelbottom=False, bottom=False)
        if kind not in first_of:
            first_of[kind] = ax
            ax.set_title({"truth": "truth (simulator)", "RCK MC": "RCK MC (strict species; haplotypes within phase blocks)",
                          "CoRAL": "CoRAL (strict cycles; no alleles)"}[kind], loc="left", fontsize=9, color=rp.INK)
        axes.append(ax)
    axes[-1].tick_params(axis="x", labelbottom=True, bottom=True)
    rp.position_ticks(axes[-1], axis)
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


def rck_mc_pages(case, folder, work):
    run = os.path.join(work, "rck_mc")
    shutil.copytree(alb.bi.run_dir(case, "rck_mc"), run, ignore=shutil.ignore_patterns("score"))
    subprocess.run([sys.executable, "-m", "rck.ecdna.plot", "--outdir", run, "--cnplot-path", mtp.CNPLOT,
                    "--rck-input", case["rck_input"]], check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL,
                   cwd=mtp.RCK_ROOT)
    for png in sorted(glob.glob(os.path.join(run, "plots", "*.png"))):
        shutil.copy(png, os.path.join(folder, "rck_mc_" + os.path.basename(png).replace("union_", "")))


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
    base = [sys.executable, os.path.join(HERE, "make_species_matrix_figure.py"), "--set", set_name, "--case", case["name"],
            "--match", "score", "--plot-dir", folder]
    for extra in (["--method", "rck_mc", "--block"], ["--method", "rck_mc", "--blind"],
                  ["--method", "coral", "--blind", "--variant", CORAL_VARIANT]):
        with open(os.path.join(folder, "species_matrix.log"), "a") as log:
            subprocess.run(base + extra, check=True, stdout=log, stderr=subprocess.STDOUT, cwd=HERE)


def readme(set_name, case, why, folder):
    import pandas as pd

    table = os.path.join(HERE, "bench_iw_allelic_table.tsv") if set_name == "iw_allelic" else \
        os.path.join(HERE, "..", "bench_IW", "bench_iw_intpsi_table.tsv")
    t = pd.read_csv(table, sep="\t").set_index(["name", "method_variant"])
    cols = ["n_truth", "n_strict", "n_missed", "n_false", "blind_exact_cycles", "phased_exact_cycles", "blind_jseg_best_mean",
            "blind_jadj_best_mean", "cn_recall", "cn_precision", "proportion_tv"]
    lines = ["{n} ({s}): {w}".format(n=case["name"], s=set_name, w=why), "truth run: " + case["truth"], "",
             "scores (" + os.path.relpath(table, HERE) + "; phased = within phase blocks, RCK MC only):"]
    lines.append("  " + "  ".join("{:>22}".format(c) if i == 0 else c for i, c in enumerate(["method"] + cols)))
    for mv, label in (("rck_mc", "RCK MC"), ("coral:" + CORAL_VARIANT, "CoRAL")):
        r = t.loc[(case["name"], mv)]
        lines.append("  {:>22}  ".format(label) + "  ".join("{:.2f}".format(r[c]) if pd.notna(r[c]) else "-" for c in cols))
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
             ("rck_mc pages", lambda: rck_mc_pages(case, folder, work)),
             ("coral figures", lambda: coral_figures(case, folder, work)),
             ("species matrices", lambda: matrices(set_name, case, folder))]
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


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--index", type=int, help="case index in CASES (0-based; the Slurm array task id)")
    p.add_argument("--list", action="store_true")
    p.add_argument("--readme", action="store_true")
    p.add_argument("--coral-only", action="store_true", help="with --index: redo only the CoRAL figures")
    a = p.parse_args()
    if a.list:
        for k, (s, n, _) in enumerate(CASES):
            print(k, s, n, sep="\t")
        return 0
    if a.readme:
        write_index()
        return 0
    if a.index is None:
        p.error("--index, --list or --readme")
    import logging
    logging.disable(logging.INFO)
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

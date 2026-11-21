#!/usr/bin/env python3
"""Topology figures of solved bench_IW_intpsi cases, for choosing presentation figures (PLAN.md item 7).

For each case in CASES (simple → hard; one species and one allele up to four nested species with both
alleles amplified), writes into topology_plots/<k>_<case>/:

  rck_bulk_<page>.png        rck-ecdna-plot (rck.ecdna.plot) of the RCK bulk run: per species, the walk
                             coloured by haplotype (A blue, B orange) with the matched truth species, the
                             junction arcs, allele-specific copy numbers (A up, B mirrored) and a ring view;
                             plus the recovery matrix.
  rck_mc_<page>.png          rck-ecdna-plot of the multi-clone run (the method under test), same layout:
                             walks coloured by haplotype, copy-number tracks in bulk terms (sum_k u_k x clone
                             profile + normal), psi per clone in the title; chromosomal background paths left
                             out (rck.ecdna.plot docstring).
  coral_rck_bulk_cycles.png  CoRAL 3.1.1 `coral plot` of the same answer (its AA cycles file): CoRAL's own
  coral_rck_mc_cycles.png    figure style, haplotype-blind (CoRAL ignores Haplotypes=); and of the multi-clone
                             route's answer (the method under test)
  coral_rck_bulk_graph.png   CoRAL's graph figure (copy numbers and junction arcs) of RCK bulk's graph

CoRAL is used unchanged. Its inputs are plotting copies, adapted for legibility only: the cycles files keep
only the walks marked Status=species (the answer every score uses; low_weight walks and the multi-clone
route's chromosomal background paths, which run through tens of Mb of background, are dropped); the
cycles figures are drawn without --graph (CoRAL then takes its x-axis intervals from the cycle segments, merged within 10 kb, instead
of from the graph's sequence edges, which span the background between amplicon pieces: tens of Mb), and the
graph figure gets a plotting copy of the graph with only the sequence edges at least half covered by
segments the species use, padded by 10 kb (and the edges between kept sequence edges). --hide-genes: the simulated amplicons sit on real hg38
loci whose gene tracks would cover the figure.

Truth only (`--truth RUN_DIR ...`, PLAN.md item 1 dry run): for each simulator run, into
topology_plots/<--out>/<run>/: truth_topology.png (every truth species as a topology row of rck-ecdna-plot,
coloured by haplotype, on one genome axis; bulk copies and psi per clone in the labels) and
coral_truth_<species>_cycles.png (CoRAL's cycles figure of each species' truth/aa cycles, haplotype-blind).

    python make_topology_plots.py --truth RUN_DIR [RUN_DIR ...] --out item1_dryrun

Run folders are only read: the RCK bulk run is copied (without pool/ and score/) to a temporary folder for
rck-ecdna-plot, which writes next to its input. ~1 min per case, so the full set runs as a Slurm job:

    python make_topology_plots.py [CASE ...]        # default: all of CASES
    sbatch --time=00:30:00 --mem=4G --cpus-per-task=2 --wrap "python make_topology_plots.py"
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
import allele_bench as ab  # noqa: E402  (configures bench_iw for iw_intpsi)

bi = ab.bi
OUT = os.path.join(HERE, "topology_plots")
RCK_ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
CORAL = "/n/fs/yubi-iw/conda/envs/coral-3.1.1/bin/coral"
CNPLOT = os.path.expanduser("~/amplicon-allele-analysis-cnp/third_party/cnplot")
PAD = 10_000
# (case, one-line description), simple → hard; mixed = share of ecDNA bases with both alleles amplified
CASES = [
    ("simple_S1_s11001", "1 species, one allele"),
    ("hard_S1_s14003", "1 species with internal duplications / fold-back, one allele"),
    ("simple_S2_s12003", "2 distinct species (no shared bases), one allele each"),
    ("simple_S4_s13001", "4 distinct species, one allele each"),
    ("hard_S2_s15001", "2 nested species, one allele"),
    ("hard_S2_s15003", "2 nested species, both alleles amplified (62% of ecDNA bases)"),
    ("hard_S4_s16005", "4 nested species, both alleles amplified (56%)"),
    ("hard_S4_s16003", "4 nested species, both alleles amplified (70%); the hardest case, multi-clone's weakest"),
]


def answer_files(case, method):
    folder = os.path.join(bi.score_folder(case, method), "answer")
    stem = os.path.join(folder, "{c}_{m}".format(c=case["name"], m=method))
    return stem + "_cycles.txt", stem + "_graph.txt"


def cycle_segments(cycles_path):
    """[(chrom, start, end)] of the Segment lines a cycles file lists."""
    out = []
    with open(cycles_path) as handle:
        for line in handle:
            f = line.split()
            if f and f[0] == "Segment":
                out.append((f[2], int(f[3]), int(f[4])))
    return out


def species_only(cycles_path, out_path):
    """A plotting copy of an AA cycles file with only its Status=species walks (renumbered 1, 2, ...)."""
    k = 0
    with open(cycles_path) as src, open(out_path, "w") as dst:
        for line in src:
            if line.startswith("Cycle="):
                fields = dict(p.split("=", 1) for p in line.rstrip("\n").split(";") if "=" in p)
                if fields.get("Status", "species") != "species":
                    continue
                k += 1
                line = "Cycle={k};".format(k=k) + line.split(";", 1)[1]
            dst.write(line)
    return k


def used_segments(cycles_path):
    """Segment ids the walks of a cycles file use."""
    used = set()
    with open(cycles_path) as handle:
        for line in handle:
            if line.startswith("Cycle="):
                fields = dict(p.split("=", 1) for p in line.rstrip("\n").split(";") if "=" in p)
                used.update(t[:-1] for t in fields.get("Segments", "").split(",") if t[:-1] not in ("", "0"))
    return used


def trimmed_graph(graph_path, cycles_path, out_path):
    """A plotting copy of an AA graph: sequence edges within PAD of a cycle segment, and the concordant /
    discordant / source edges whose breakends both lie in kept sequence edges."""
    used = used_segments(cycles_path)
    segs = []
    with open(cycles_path) as handle:
        for line in handle:
            f = line.split()
            if f and f[0] == "Segment" and f[1] in used:
                segs.append((f[2], int(f[3]), int(f[4])))

    def near(chrom, start, end):
        """at least half of the edge lies on (padded) segments the species use"""
        cover = sum(max(0, min(end, e + PAD) - max(start, s - PAD) + 1) for c, s, e in segs if c == chrom)
        return cover >= 0.5 * (end - start + 1)

    def pos(token):  # chr16:13455002- -> (chr16, 13455002)
        chrom, rest = token.rsplit(":", 1)
        return chrom, int(rest[:-1])

    with open(graph_path) as handle:
        lines = handle.readlines()
    kept = []
    for line in lines:
        f = line.split("\t")
        if f[0] == "sequence":
            (c1, a), (_, b) = pos(f[1]), pos(f[2])
            if near(c1, a, b):
                kept.append((c1, a, b))
    def inside(token):
        c, p = pos(token)
        return any(c == k[0] and k[1] <= p <= k[2] for k in kept)

    with open(out_path, "w") as out:
        for line in lines:
            f = line.split("\t")
            if f[0] == "sequence":
                (c1, a), (_, b) = pos(f[1]), pos(f[2])
                if (c1, a, b) in kept:
                    out.write(line)
            elif f[0] in ("concordant", "discordant") and "->" in f[1]:
                left, right = f[1].split("->")
                if inside(left) and inside(right):
                    out.write(line)
            elif f[0] == "source":
                continue
            else:
                out.write(line)


def coral_plot(prefix, cycles=None, graph=None):
    cmd = [CORAL, "plot", "--ref", "hg38", "--output-prefix", prefix, "--hide-genes"]
    cmd += ["--cycles", cycles] if cycles else []
    cmd += ["--graph", graph] if graph else []
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)


def plot_case(rank, case, description):
    name = case["name"]
    folder = os.path.join(OUT, "{k}_{n}".format(k=rank, n=name))
    os.makedirs(folder, exist_ok=True)
    with open(os.path.join(folder, "README.txt"), "w") as handle:
        handle.write("{n}: {d}\nclass {c}, S = {s}; truth run {t}\n".format(n=name, d=description, c=case["class"],
                                                                             s=case["S"], t=case["truth"]))
    work = tempfile.mkdtemp(prefix="topo_")
    try:
        # rck-ecdna-plot on a copy of the RCK bulk run
        run = os.path.join(work, "rck_bulk")
        shutil.copytree(bi.run_dir(case, "rck_bulk"), run, ignore=shutil.ignore_patterns("pool", "score"))
        subprocess.run([sys.executable, "-m", "rck.ecdna.plot", "--outdir", run, "--cnplot-path", CNPLOT],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=RCK_ROOT)
        for png in sorted(glob.glob(os.path.join(run, "plots", "*.png"))):
            shutil.copy(png, os.path.join(folder, "rck_bulk_" + os.path.basename(png).replace("union_", "")))
        # rck-ecdna-plot on a copy of the multi-clone run (reads its input DIAG from summary.json, read only)
        run = os.path.join(work, "rck_mc")
        shutil.copytree(bi.run_dir(case, "rck_mc"), run, ignore=shutil.ignore_patterns("score"))
        subprocess.run([sys.executable, "-m", "rck.ecdna.plot", "--outdir", run, "--cnplot-path", CNPLOT,
                        "--rck-input", case["rck_input"]],
                       check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL, cwd=RCK_ROOT)
        for png in sorted(glob.glob(os.path.join(run, "plots", "*.png"))):
            shutil.copy(png, os.path.join(folder, "rck_mc_" + os.path.basename(png).replace("union_", "")))
        # CoRAL: cycles of both RCK routes (no --graph: compact axis), graph of RCK bulk (trimmed copy)
        for method in ("rck_bulk", "rck_mc"):
            cycles, _ = answer_files(case, method)
            species = os.path.join(work, "{m}_species_cycles.txt".format(m=method))
            species_only(cycles, species)
            coral_plot(os.path.join(work, "coral_" + method), cycles=species)
            shutil.copy(os.path.join(work, "coral_{m}_cycles.png".format(m=method)), folder)
        _, graph = answer_files(case, "rck_bulk")
        small = os.path.join(work, "graph_trimmed.txt")
        trimmed_graph(graph, os.path.join(work, "rck_bulk_species_cycles.txt"), small)
        coral_plot(os.path.join(work, "coral_rck_bulk"), graph=small)
        shutil.copy(os.path.join(work, "coral_rck_bulk_graph.png"), folder)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return folder


def clone_psi(run_dir):
    """{species: {clone: psi}} and {species: bulk copies of tumor1} from truth/proportions.tsv."""
    import pandas as pd
    p = pd.read_csv(os.path.join(run_dir, "truth", "proportions.tsv"), sep="\t")
    psi, bulk = {}, {}
    for r in p.itertuples():
        if r.KIND == "ecdna_psi":
            psi.setdefault(str(r.ID), {})[str(r.CLONE)] = float(r.PSI)
        elif r.KIND == "ecdna_burden" and r.SAMPLE == "tumor1":
            bulk[str(r.ID)] = float(r.PHI)
    return psi, bulk


def plot_truth(run_dir, folder):
    """truth_topology.png and CoRAL's truth cycles figures of one simulator run (module docstring)."""
    from rck.ecdna import evaluate as ev
    from rck.ecdna import plot as rp

    name = os.path.basename(run_dir.rstrip("/"))
    os.makedirs(folder, exist_ok=True)
    plt, _, cnplot = rp.import_plotting(CNPLOT)
    truth = ev.load_truth(run_dir)
    psi, bulk = clone_psi(run_dir)
    species = [(c.name, [rp.PStep(s.chrom, s.start, s.end, s.hap, s.strand) for s in c.steps]) for c in truth]
    work = tempfile.mkdtemp(prefix="topo_truth_")
    try:
        regions = rp.axis_regions([st for _, steps in species for st in steps])
        axis = rp.build_axis(cnplot, regions, work)
        heights = [max(1.4, 0.22 * len(steps)) for _, steps in species]
        height = min(sum(heights) + 1.6, 30)
        fig = plt.figure(figsize=(15, height))
        grid = fig.add_gridspec(len(heights), 1, height_ratios=heights, hspace=0.25, top=1.0 - 0.8 / height, bottom=0.7 / height)
        axes = []
        for row, (sp, steps) in enumerate(species):
            ax = fig.add_subplot(grid[row, 0], sharex=axes[0] if axes else None)
            base = sp.split(".")[0]
            label = "{s}\n{b:.1f} bulk copies\npsi {p}".format(s=sp, b=bulk.get(base, float("nan")), p=", ".join(
                "{c} {v:g}".format(c=c, v=v) for c, v in sorted(psi.get(base, {}).items())))
            rp.draw_topology(ax, axis, steps, True, label)
            ax.set_xlim(axis.ch_coords[0], axis.chr_end)
            for boundary in axis.chr_boundaries:
                ax.axvline(boundary, color=rp.INK_MUTED, linewidth=0.8, zorder=0)
            ax.tick_params(axis="x", labelbottom=False, bottom=False)
            axes.append(ax)
        axes[-1].tick_params(axis="x", labelbottom=True, bottom=True)
        rp.position_ticks(axes[-1], axis)
        from matplotlib.patches import Patch
        fig.legend(handles=[Patch(facecolor=rp.HAP_COLORS["A"], label="haplotype A"),
                            Patch(facecolor=rp.HAP_COLORS["B"], label="haplotype B")],
                   loc="upper right", bbox_to_anchor=(0.98, 1.0), ncol=2, fontsize=8, frameon=False)
        fig.suptitle("{n}: truth (simulator), one row per step in walk order".format(n=name), fontsize=10, color=rp.INK,
                     x=0.02, y=1.0 - 0.1 / height, ha="left", va="top")
        fig.savefig(os.path.join(folder, "truth_topology.png"), dpi=150)
        plt.close(fig)
        for path in sorted(glob.glob(os.path.join(run_dir, "truth", "aa", "ecdna*_cycles.txt"))):
            sp = os.path.basename(path)[:-len("_cycles.txt")]
            coral_plot(os.path.join(work, "coral_truth_" + sp), cycles=path)
            shutil.copy(os.path.join(work, "coral_truth_{s}_cycles.png".format(s=sp)), folder)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return folder


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("cases", nargs="*", help="case names (default: all of CASES)")
    p.add_argument("--truth", nargs="+", default=None, metavar="RUN_DIR", help="draw these simulator runs' truth only")
    p.add_argument("--out", default="truth", help="with --truth: folder under topology_plots/ (default truth)")
    a = p.parse_args()
    if a.truth:
        for run_dir in a.truth:
            print(plot_truth(run_dir, os.path.join(OUT, a.out, os.path.basename(run_dir.rstrip("/")))), flush=True)
        return
    by_name = {c["name"]: c for c in bi.cases()}
    for rank, (name, description) in enumerate(CASES, start=1):
        if a.cases and name not in a.cases:
            continue
        print(plot_case(rank, by_name[name], description), flush=True)
    with open(os.path.join(OUT, "README.md"), "w") as handle:
        handle.write("# Topology figures (make_topology_plots.py)\n\nSimple → hard; the script's docstring says what "
                     "each file is. rck_bulk_* / rck_mc_* pages (rck-ecdna-plot, RCK bulk and multi-clone routes): "
                     "haplotype A blue, B orange. coral_* figures: CoRAL 3.1.1's own style, which does not show alleles.\n\n| folder | case |\n|---|---|\n")
        for rank, (name, description) in enumerate(CASES, start=1):
            handle.write("| `{k}_{n}` | {d} |\n".format(k=rank, n=name, d=description))


if __name__ == "__main__":
    main()

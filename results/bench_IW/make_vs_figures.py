#!/usr/bin/env python
"""Absolute ecDNA species copies per cell (v_s): inferred vs truth, RCK MC / CoRAL / CE, 36 cases.

Cases: the 24 bench_IW_intpsi cases (set iw_intpsi) and the 12 allelic-sibling cases of bench_IW_allele
(set iw_allelic). Baselines: CoRAL / CE with linear walks kept (coral:walks_kept, ce:walks_kept).

Where the numbers come from (read only; no method is rerun, no Gurobi). bench_metrics.proportion_scores
writes truth_proportion = PHI_s / sum PHI and inferred_proportion = w_c / sum w_c, where PHI_s is the truth
species' weight from evaluate.load_truth (bulk copies per cell, sum_k u_k psi*_{s,k}) and w_c the weight of
the matched cycle in the strict answer (bulk copies per cell; for rck_mc reweighted to sum_k u_k psi_k by
bench_metrics.multiclone_bulk_weights, with u from the DIAG's clone_proportions.tsv). The totals that divide
them are not stored (truth_mass / inferred_mass in the tables are base-pair-weighted, bp x copies), so this
script reads the same weights again, exactly as bench_iw.score_one does:

    PHI_s  evaluate.load_truth(<sim run>)                                 per truth species
    w_c    <run>/score/answer/*cycles*.txt, strict set (aa_eval.strict_select); rck_mc reweighted as above

and takes the pairs (truth species, matched cycle) from the collected proportions tables
(bench_iw_intpsi_proportions.tsv, bench_iw_allelic_proportions.tsv; set = strict, the scorer's one-to-one
haplotype-blind Jaccard matching). As a check, PHI_s / sum PHI and w_c / sum w_c are recomputed and must
equal the stored proportions (to 1e-9), or the script stops. rck_mc's evaluation_psi.tsv is not used: it
holds per-clone psi with its own (absolute-frame) pairing and exists for rck_mc only.

Per pair: rel_err = |w_c - PHI_s| / PHI_s. A truth species with no matched cycle is "missed": it has no
inferred v_s, appears in the per-case figure as an x at 0, and is counted (not plotted) in the summary.

Writes to figures_vs/ (or --out):
    vs_summary.{pdf,png}     inferred vs truth v_s (log-log, y = x), one panel per method, marker = set;
                             and per-species relative error per method (box + points)
    vs_per_case.{pdf,png}    one panel per case: inferred vs truth v_s, each method, species named
    figure_data/vs_pairs.tsv every (case, method, truth species) row the figures plot

Usage (from the repository root RCK-cnp, gurobi-env python; Gurobi itself is not used):
    python results/bench_IW/make_vs_figures.py
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, ROOT)
sys.path.insert(0, HERE)

from rck.ecdna import aa_eval, bench_metrics as bm, evaluate as ev  # noqa: E402
from rck.ecdna import aa_format as aaf  # noqa: E402

# Style as make_figures.py (same palette slots, labels, rcParams).
METHODS = ["rck_mc", "coral", "ce"]
VARIANT = {"rck_mc": "rck_mc", "coral": "coral:walks_kept", "ce": "ce:walks_kept"}
COLOR = {"rck_mc": "#eb6834", "coral": "#1baf7a", "ce": "#eda100"}
LABEL = {"rck_mc": "RCK MC", "coral": "CoRAL", "ce": "CE"}
SETS = ["iw_intpsi", "iw_allelic"]
SET_LABEL = {"iw_intpsi": "bench_IW_intpsi (24 cases)", "iw_allelic": "bench_IW_allele (12 cases)"}
MARKER = {"iw_intpsi": "o", "iw_allelic": "s"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"

plt.rcParams.update({
    "font.family": "sans-serif", "font.size": 8, "axes.titlesize": 9, "axes.labelsize": 8,
    "xtick.labelsize": 7, "ytick.labelsize": 7, "legend.fontsize": 7,
    "axes.edgecolor": MUTED, "axes.labelcolor": INK2, "xtick.color": INK2, "ytick.color": INK2,
    "text.color": INK, "axes.spines.top": False, "axes.spines.right": False,
    "axes.grid": True, "grid.color": GRID, "grid.linewidth": 0.6, "axes.axisbelow": True,
    "lines.linewidth": 2, "savefig.dpi": 300, "pdf.fonttype": 42, "ps.fonttype": 42,
    "figure.facecolor": "white", "axes.facecolor": "white",
})


# -- data -------------------------------------------------------------------------------------

def strict_weights(bi, case, method):
    """name -> w_c (bulk copies per cell) of the strict answer the scorer used."""
    folder = bi.score_folder(case, method)
    found = glob.glob(os.path.join(folder, "answer", "*cycles*.txt"))
    if len(found) != 1:
        raise FileNotFoundError("{f}: expected one answer cycles file, found {n}".format(f=folder, n=len(found)))
    records = aaf.read_cycles(found[0])
    strict, _ = aa_eval.strict_select(records, aa_eval.species_cycles(records), "auto")
    if method == "rck_mc":
        props = bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv"))
        strict = bm.reweight(strict, bm.multiclone_bulk_weights(bi.run_dir(case, method), props))
    return {c.name: float(c.weight) for c in strict}


def collect_set(bi, set_name, proportions_path):
    props = pd.read_csv(proportions_path, sep="\t")
    props = props[(props["set"] == "strict") & props["truth"].notna()]
    rows = []
    for case in bi.cases():
        truth = {t.name: float(t.weight) for t in ev.load_truth(case["truth"])}
        t_total = sum(truth.values())
        for method in METHODS:
            pairs = props[(props["case"] == case["name"]) & (props["method_variant"] == VARIANT[method])]
            if pairs.empty:
                print("no scored pairs: {c} {m}".format(c=case["name"], m=VARIANT[method]))
                continue
            w = strict_weights(bi, case, method)
            c_total = sum(w.values())
            for p in pairs.itertuples():
                phi = truth[p.truth]
                if abs(phi / t_total - p.truth_proportion) > 1e-9:
                    raise ValueError("truth proportion mismatch {c} {m} {t}".format(c=case["name"], m=method, t=p.truth))
                matched = isinstance(p.inferred, str) and p.inferred != ""
                wc = w[p.inferred] if matched else np.nan
                if matched and abs(wc / c_total - p.inferred_proportion) > 1e-9:
                    raise ValueError("inferred proportion mismatch {c} {m} {i}".format(c=case["name"], m=method, i=p.inferred))
                rows.append({"set": set_name, "case": case["name"], "scenario": case["scenario"], "method": method,
                             "method_variant": VARIANT[method], "truth": p.truth, "inferred": p.inferred if matched else "",
                             "truth_v": phi, "inferred_v": wc, "rel_err": abs(wc - phi) / phi if matched else np.nan,
                             "truth_proportion": p.truth_proportion,
                             "inferred_proportion": p.inferred_proportion if matched else np.nan,
                             "jaccard_blind": p.jaccard_blind if matched else np.nan})
    return rows


def load():
    import bench_iw as bi

    bi.configure("iw_intpsi")
    rows = collect_set(bi, "iw_intpsi", os.path.join(HERE, "bench_iw_intpsi_proportions.tsv"))
    sys.path.insert(0, os.path.join(HERE, "..", "bench_IW_allele"))
    import allelic_bench  # noqa: F401 - points bench_iw at iw_allelic (paths, cases, patched CoRAL / CE graph)

    rows += collect_set(bi, "iw_allelic", os.path.join(HERE, "..", "bench_IW_allele", "bench_iw_allelic_proportions.tsv"))
    return pd.DataFrame(rows)


# -- figures ----------------------------------------------------------------------------------

def save(fig, out, name):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, name + "." + ext), bbox_inches="tight")
    plt.close(fig)


def summary_figure(d, out):
    m = d[d["inferred_v"].notna()]
    lo = min(m["truth_v"].min(), m["inferred_v"].min()) / 1.3
    hi = max(m["truth_v"].max(), m["inferred_v"].max()) * 1.3
    fig = plt.figure(figsize=(7.2, 4.6))
    grid = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.75], hspace=0.45, wspace=0.3)
    for i, method in enumerate(METHODS):
        ax = fig.add_subplot(grid[0, i])
        ax.plot([lo, hi], [lo, hi], color=MUTED, lw=1, ls="--", zorder=1)
        sub = d[d["method"] == method]
        for s in SETS:
            pts = sub[(sub["set"] == s) & sub["inferred_v"].notna()]
            ax.scatter(pts["truth_v"], pts["inferred_v"], s=22, marker=MARKER[s], color=COLOR[method],
                       edgecolor="white", linewidth=0.6, alpha=0.9, zorder=3)
        n, missed = len(sub), int(sub["inferred_v"].isna().sum())
        ax.set_xscale("log")
        ax.set_yscale("log")
        ax.set_xlim(lo, hi)
        ax.set_ylim(lo, hi)
        ax.set_aspect("equal")
        ax.set_title(LABEL[method])
        ax.text(0.04, 0.96, "{k} / {n} matched".format(k=n - missed, n=n), transform=ax.transAxes,
                ha="left", va="top", fontsize=7, color=INK2)
        ax.set_xlabel("truth $v_s$ (copies / cell)")
        if i == 0:
            ax.set_ylabel("inferred $v_s$ (copies / cell)")
    handles = [Line2D([], [], ls="", marker=MARKER[s], color=INK2, label=SET_LABEL[s]) for s in SETS]
    handles.append(Line2D([], [], color=MUTED, lw=1, ls="--", label="y = x"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=3, frameon=False)

    ax = fig.add_subplot(grid[1, :])
    rng = np.random.default_rng(0)
    data = [m.loc[m["method"] == k, "rel_err"].values for k in METHODS]
    ax.boxplot(data, positions=range(len(METHODS)), widths=0.5, orientation="horizontal", showfliers=False,
               medianprops={"color": INK, "lw": 1.5}, boxprops={"color": MUTED}, whiskerprops={"color": MUTED},
               capprops={"color": MUTED})
    for i, method in enumerate(METHODS):
        for s in SETS:
            e = m.loc[(m["method"] == method) & (m["set"] == s), "rel_err"].values
            ax.scatter(e, i + rng.uniform(-0.18, 0.18, len(e)), s=12, marker=MARKER[s], color=COLOR[method],
                       edgecolor="white", linewidth=0.4, alpha=0.85, zorder=3)
        med = np.median(data[i]) if len(data[i]) else np.nan
        label = "median < 0.001" if med < 0.001 else "median {v:.3f}".format(v=med)
        ax.text(1.01, i, label, transform=ax.get_yaxis_transform(), ha="left",
                va="center", fontsize=7, color=INK2)
    ax.set_yticks(range(len(METHODS)))
    ax.set_yticklabels([LABEL[k] for k in METHODS])
    ax.invert_yaxis()
    ax.set_xscale("symlog", linthresh=0.01)
    ax.set_xlim(0, m["rel_err"].max() * 1.4)
    ax.set_xlabel("per-species relative error $|w_c - \\Phi_s| / \\Phi_s$ (matched pairs; symlog, linear below 0.01)")
    ax.grid(axis="y", visible=False)
    save(fig, out, "vs_summary")


def per_case_figure(d, out):
    cases = list(dict.fromkeys(d["case"]))
    ncol = 6
    nrow = int(np.ceil(len(cases) / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(ncol * 1.9, nrow * 1.9))
    offset = {"rck_mc": -1, "coral": 0, "ce": 1}
    for ax, case in zip(axes.flat, cases):
        sub = d[d["case"] == case]
        top = np.nanmax(np.r_[sub["truth_v"].values, sub["inferred_v"].values]) * 1.12
        ax.plot([0, top], [0, top], color=MUTED, lw=0.8, ls="--", zorder=1)
        species = sub.drop_duplicates("truth").sort_values("truth_v")
        for k, sp in enumerate(species.itertuples()):  # species named once, at its truth v_s
            ax.axvline(sp.truth_v, color=GRID, lw=0.8, zorder=0)
            ax.text(sp.truth_v, top * (1.02 + 0.07 * (k % 3)), sp.truth.replace("ecdna", "ec"), ha="center",
                    va="bottom", fontsize=5, color=INK2, clip_on=False)
        dx = top * 0.018
        for method in METHODS:
            pts = sub[sub["method"] == method]
            hit, miss = pts[pts["inferred_v"].notna()], pts[pts["inferred_v"].isna()]
            ax.scatter(hit["truth_v"] + offset[method] * dx, hit["inferred_v"], s=16, color=COLOR[method],
                       edgecolor="white", linewidth=0.5, zorder=3)
            ax.scatter(miss["truth_v"] + offset[method] * dx, np.zeros(len(miss)), s=18, marker="x",
                       color=COLOR[method], linewidth=1.1, zorder=3, clip_on=False)
        ax.set_xlim(0, top)
        ax.set_ylim(0, top)
        ax.set_aspect("equal")
        ax.set_title(case, fontsize=6.5, pad=24)
        ax.tick_params(labelsize=5.5)
    for ax in axes.flat[len(cases):]:
        ax.axis("off")
    for ax in axes[:, 0]:
        ax.set_ylabel("inferred $v_s$", fontsize=7)
    for ax in axes[-1, :]:
        ax.set_xlabel("truth $v_s$", fontsize=7)
    handles = [Line2D([], [], ls="", marker="o", color=COLOR[k], label=LABEL[k]) for k in METHODS]
    handles += [Line2D([], [], ls="", marker="x", color=INK2, label="missed (no matched cycle)"),
                Line2D([], [], color=MUTED, lw=0.8, ls="--", label="y = x")]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.0), ncol=5, frameon=False, fontsize=8)
    fig.suptitle("ecDNA copies per cell per species: rows 1-4 bench_IW_intpsi, rows 5-6 bench_IW_allele; "
                 "species named above their truth $v_s$; methods offset slightly in x", y=1.025, fontsize=8, color=INK2)
    fig.tight_layout(rect=(0, 0, 1, 0.985))
    save(fig, out, "vs_per_case")


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--out", default=os.path.join(HERE, "figures_vs"))
    args = parser.parse_args()
    os.makedirs(os.path.join(args.out, "figure_data"), exist_ok=True)
    d = load()
    d.to_csv(os.path.join(args.out, "figure_data", "vs_pairs.tsv"), sep="\t", index=False)
    summary_figure(d, args.out)
    per_case_figure(d, args.out)
    g = d.groupby(["method", "set"]).agg(n=("truth", "size"), matched=("inferred_v", "count"),
                                         median_rel_err=("rel_err", "median"), mean_rel_err=("rel_err", "mean"))
    print(g.to_string())


if __name__ == "__main__":
    main()

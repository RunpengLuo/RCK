#!/usr/bin/env python
"""Manuscript figures 2-6 and Table 1 for bench_IW_intpsi.

Reads only the collected tables in this folder (no runs are touched, no Gurobi):
    bench_iw_intpsi_table.tsv        one row per (case, method variant)
    bench_iw_intpsi_proportions.tsv  truth / inferred species proportion pairs

Writes to figures_intpsi/ (or --out):
    fig2_cn_recall_precision.{pdf,png}   ecDNA copy-number recall vs precision, one point per case
    fig3_structure_vs_S.{pdf,png}        blind J_seg / J_adj per truth species vs number of species
    fig4_mass_ratio.{pdf,png}            inferred / true ecDNA copy-number mass (log scale)
    fig5_proportions.{pdf,png}           species proportions: inferred vs truth, and TV per case
    fig6_allele_specific.{pdf,png}       RCK only: blind vs phased J_seg; within-block vs absolute counts
    table1.tsv / table1.md               mean +- SD per method variant, overall and per class
    table1_paired_tests.tsv              paired Wilcoxon tests, RCK routes vs each baseline variant
    figure_data/*.tsv                    the exact numbers each figure plots

Usage (from results/bench_IW):
    python make_figures.py                                  # baselines = CoRAL/CE with linear walks kept
    python make_figures.py --baseline-variant walks_dropped --out figures_intpsi_walks_dropped

Conventions (see RESULTS.md):
- Main figures show four methods: rck_bulk, rck_mc, CoRAL, CE. CoRAL / CE are shown as one walk
  variant (--baseline-variant); Table 1 always lists all six method variants.
- A run with no cycle at all (CoRAL with walks dropped returns only linear walks in 3 cases) has no
  best-match Jaccard or exact cycles in the collected table (NaN). Here they count as 0: the method
  returned no ecDNA, so nothing matches. The pooled TSV from `bench_iw.py collect` skips them instead.
- Exact-cycle rate = exact cycles summed over cases / truth species summed over cases.
"""
import argparse
import os

import numpy as np
import pandas as pd

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402
from scipy.stats import wilcoxon  # noqa: E402

HERE = os.path.dirname(os.path.abspath(__file__))

# Categorical slots 1-4 of the dataviz reference palette, in their validated order (adjacent
# pairs pass CVD >= 8 and normal-vision >= 15). Colour follows the method everywhere. Four colours
# fail the all-pairs check, so no panel relies on colour alone: every plot also separates methods
# by x position or by panel, and carries a legend or axis labels.
METHODS = ["rck_bulk", "rck_mc", "coral", "ce"]
COLOR = {"rck_bulk": "#2a78d6", "rck_mc": "#eb6834", "coral": "#1baf7a", "ce": "#eda100"}
LABEL = {"rck_bulk": "RCK-ecDNA (bulk)", "rck_mc": "RCK-ecDNA (multi-clone*)",
         "coral": "CoRAL†", "ce": "CycleExtractor†"}
SHORT = {"rck_bulk": "RCK bulk", "rck_mc": "RCK multi-clone*", "coral": "CoRAL†", "ce": "CE†"}
VARIANT_ORDER = ["rck_bulk", "rck_mc", "coral:walks_kept", "coral:walks_dropped",
                 "ce:walks_kept", "ce:walks_dropped"]
CLASS_ORDER = ["simple", "hard"]
CLASS_LABEL = {"simple": "simple", "hard": "nested"}
MARKER = {"simple": "o", "hard": "^"}
INK, INK2, MUTED, GRID = "#0b0b0b", "#52514e", "#8a8984", "#e6e5e0"
FOOT = ("* given the true clone profiles (upper bound for the multi-clone route).  "
        "† decomposition only, on RCK bulk's graph.")

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

def load(table_path, baseline_variant):
    t = pd.read_csv(table_path, sep="\t")
    no_cycles = t["n_cycles"].fillna(0) == 0
    for col in ["blind_jseg_best_mean", "blind_jadj_best_mean", "blind_exact_cycles"]:
        t.loc[no_cycles & t[col].isna(), col] = 0.0
    t["mass_ratio"] = t["inferred_mass"] / t["truth_mass"]
    t["errors_within"] = t["n_missed"] + t["n_false"]
    t["errors_absolute"] = t["n_missed_absolute"] + t["n_false_absolute"]
    t["variant"] = t["method_variant"]
    keep = {"rck_bulk": "rck_bulk", "rck_mc": "rck_mc",
            "coral:" + baseline_variant: "coral", "ce:" + baseline_variant: "ce"}
    main = t[t["variant"].isin(keep)].copy()
    main["m"] = main["variant"].map(keep)
    return t, main, int(no_cycles.sum())


def jitter(n, width=0.18, seed=0):
    return np.random.default_rng(seed).uniform(-width, width, n)


def save(fig, out, name, data=None):
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(out, "{}.{}".format(name, ext)), bbox_inches="tight")
    plt.close(fig)
    if data is not None:
        data.to_csv(os.path.join(out, "figure_data", name + ".tsv"), sep="\t", index=False)


def footnote(fig, y=-0.02):
    fig.text(0.01, y, FOOT, fontsize=6.5, color=INK2, ha="left", va="top")


def class_legend(ax, loc="lower left"):
    handles = [Line2D([], [], marker=MARKER[c], linestyle="", markersize=6, markerfacecolor="white",
                      markeredgecolor=INK2, label=CLASS_LABEL[c]) for c in CLASS_ORDER]
    ax.legend(handles=handles, loc=loc, frameon=False, handletextpad=0.3)


def dots(ax, x0, values, classes, color, seed):
    """One marker per case: filled for simple, triangles for nested; white ring for overlaps."""
    xs = x0 + jitter(len(values), seed=seed)
    for c in CLASS_ORDER:
        sel = np.asarray(classes) == c
        ax.scatter(xs[sel], np.asarray(values)[sel], s=26, marker=MARKER[c], color=color,
                   edgecolors="white", linewidths=0.8, alpha=0.9, zorder=3)


def mean_bar(ax, x0, values, half=0.3):
    v = np.nanmean(values)
    ax.plot([x0 - half, x0 + half], [v, v], color=INK, linewidth=1.5, zorder=4,
            solid_capstyle="round")
    return v


# -- figures ----------------------------------------------------------------------------------

def fig2(main, out):
    fig, axes = plt.subplots(1, 4, figsize=(7.2, 2.3), sharex=True, sharey=True)
    for ax, m in zip(axes, METHODS):
        d = main[main.m == m]
        for c in CLASS_ORDER:
            s = d[d["class"] == c]
            ax.scatter(s.cn_recall, s.cn_precision, s=28, marker=MARKER[c], color=COLOR[m],
                       edgecolors="white", linewidths=0.8, alpha=0.9, zorder=3)
        mr, mp = d.cn_recall.mean(), d.cn_precision.mean()
        ax.scatter([mr], [mp], marker="+", s=90, color=INK, linewidths=1.5, zorder=4)
        ax.set_title("{}\nmean {:.2f} / {:.2f}".format(SHORT[m], mr, mp), color=INK)
        ax.set_xlim(0, 1.04)
        ax.set_ylim(0, 1.04)
        ax.set_aspect("equal")
    axes[0].set_ylabel("ecDNA copy-number precision")
    fig.supxlabel("ecDNA copy-number recall", fontsize=8, color=INK2, y=0.02)
    handles = [Line2D([], [], marker=MARKER[c], linestyle="", markersize=6, markerfacecolor="white",
                      markeredgecolor=INK2, label=CLASS_LABEL[c]) for c in CLASS_ORDER]
    handles.append(Line2D([], [], marker="+", linestyle="", markersize=8, color=INK,
                          markeredgewidth=1.5, label="mean (recall / precision)"))
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 1.1), ncol=3,
               frameon=False)
    footnote(fig, y=-0.04)
    save(fig, out, "fig2_cn_recall_precision",
         main[["name", "class", "S", "m", "variant", "cn_recall", "cn_precision"]])


def fig3(main, out):
    metrics = [("blind_jseg_best_mean", "best-match J_seg (blind)"),
               ("blind_jadj_best_mean", "best-match J_adj (blind)")]
    Ss = sorted(main.S.unique())
    off = {m: (i - 1.5) * 0.12 for i, m in enumerate(METHODS)}
    fig, axes = plt.subplots(2, 2, figsize=(6.4, 4.6), sharex=True, sharey=True)
    rows = []
    for r, (col, ylab) in enumerate(metrics):
        for k, c in enumerate(CLASS_ORDER):
            ax = axes[r, k]
            for m in METHODS:
                d = main[(main.m == m) & (main["class"] == c)]
                g = d.groupby("S")[col].agg(["mean", "std"]).reindex(Ss)
                xs = np.arange(len(Ss)) + off[m]
                ax.scatter(np.concatenate([[x] * int((d.S == s).sum()) for x, s in zip(xs, Ss)]),
                           np.concatenate([d[d.S == s][col].values for s in Ss]),
                           s=10, color=COLOR[m], alpha=0.45, linewidths=0, zorder=2)
                ax.plot(xs, g["mean"], color=COLOR[m], marker="o", markersize=5,
                        markeredgecolor="white", markeredgewidth=0.8, zorder=3, label=LABEL[m])
                for s, row in g.iterrows():
                    rows.append({"metric": col, "class": c, "S": s, "m": m,
                                 "mean": row["mean"], "sd": row["std"]})
            ax.set_xticks(range(len(Ss)))
            ax.set_xticklabels([str(int(s)) for s in Ss])
            ax.set_ylim(-0.03, 1.05)
            if r == 0:
                ax.set_title("{} species".format(CLASS_LABEL[c]), color=INK)
            if r == 1:
                ax.set_xlabel("ecDNA species per case (S)")
            if k == 0:
                ax.set_ylabel(ylab)
    axes[0, 0].legend(loc="lower left", frameon=False, ncol=1)
    footnote(fig, y=0.0)
    save(fig, out, "fig3_structure_vs_S", pd.DataFrame(rows))


def fig4(main, out):
    fig, ax = plt.subplots(figsize=(4.2, 2.8))
    rows = []
    for i, m in enumerate(METHODS):
        d = main[main.m == m]
        dots(ax, i, d.mass_ratio.values, d["class"].values, COLOR[m], seed=i)
        med = np.nanmedian(d.mass_ratio)
        ax.plot([i - 0.3, i + 0.3], [med, med], color=INK, linewidth=1.5, zorder=4)
        rows.append({"m": m, "median": med, "n_over_1.5": int((d.mass_ratio > 1.5).sum()),
                     "n_over_10": int((d.mass_ratio > 10).sum()), "n": len(d)})
    ax.axhline(1, color=INK2, linewidth=0.8, linestyle="--", zorder=1)
    ax.set_yscale("log")
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels(["{}\nmedian {:.2f}".format(SHORT[r["m"]], r["median"]) for r in rows])
    ax.set_ylabel("inferred / true ecDNA copy-number mass")
    ax.set_xlim(-0.6, len(METHODS) - 0.4)
    class_legend(ax, loc="upper left")
    footnote(fig)
    save(fig, out, "fig4_mass_ratio",
         main[["name", "class", "S", "m", "variant", "inferred_mass", "truth_mass", "mass_ratio"]])
    pd.DataFrame(rows).to_csv(os.path.join(out, "figure_data", "fig4_mass_ratio_summary.tsv"),
                              sep="\t", index=False)


def fig5(main, props, baseline_variant, out):
    keep = {"rck_bulk": "rck_bulk", "rck_mc": "rck_mc",
            "coral:" + baseline_variant: "coral", "ce:" + baseline_variant: "ce"}
    p = props[(props.set == "strict") & props.method_variant.isin(keep)].copy()
    p["m"] = p.method_variant.map(keep)
    # unmatched truth species (missed) plot at inferred 0; unmatched inferred (false) at truth 0
    p["truth_proportion"] = p.truth_proportion.fillna(0.0)
    p["inferred_proportion"] = p.inferred_proportion.fillna(0.0)
    p["kind"] = np.where(p.truth.isna(), "false", np.where(p.inferred.isna(), "missed", "matched"))

    fig = plt.figure(figsize=(7.2, 4.4))
    gs = fig.add_gridspec(2, 4, height_ratios=[1, 0.9], hspace=0.55, wspace=0.35)
    for k, m in enumerate(METHODS):
        ax = fig.add_subplot(gs[0, k])
        d = p[p.m == m]
        ax.plot([0, 1], [0, 1], color=INK2, linewidth=0.8, linestyle="--", zorder=1)
        mt = d[d.kind == "matched"]
        un = d[d.kind != "matched"]
        ax.scatter(mt.truth_proportion, mt.inferred_proportion, s=16, color=COLOR[m],
                   edgecolors="white", linewidths=0.6, alpha=0.85, zorder=3)
        ax.scatter(un.truth_proportion, un.inferred_proportion, s=18, marker="x", color=COLOR[m],
                   linewidths=1.0, zorder=3)
        ax.set_xlim(-0.03, 1.03)
        ax.set_ylim(-0.03, 1.03)
        ax.set_aspect("equal")
        ax.set_title("{}\n{} missed, {} false".format(SHORT[m], (d.kind == "missed").sum(),
                                                      (d.kind == "false").sum()), color=INK)
        ax.set_xlabel("true proportion")
        if k == 0:
            ax.set_ylabel("inferred proportion")
    ax = fig.add_subplot(gs[1, :])
    for i, m in enumerate(METHODS):
        d = main[main.m == m]
        dots(ax, i, d.proportion_tv.values, d["class"].values, COLOR[m], seed=10 + i)
        v = mean_bar(ax, i, d.proportion_tv.values)
        ax.text(i + 0.33, v, "{:.2f}".format(v), fontsize=6.5, color=INK2, va="center")
    ax.set_xticks(range(len(METHODS)))
    ax.set_xticklabels([LABEL[m] for m in METHODS])
    ax.set_ylabel("proportion TV distance\n(strict cycles; lower is better)")
    ax.set_ylim(-0.02, max(0.6, main.proportion_tv.max() + 0.05))
    class_legend(ax, loc="upper left")
    footnote(fig, y=0.03)
    fig.text(0.01, 0.0, "Top: strict cycles pooled over all cases; × = missed (true species, "
             "inferred 0) or false (inferred species, truth 0).", fontsize=6.5, color=INK2, va="top")
    save(fig, out, "fig5_proportions",
         p[["case", "m", "method_variant", "truth", "inferred", "truth_proportion",
            "inferred_proportion", "kind"]])


def fig6(main, out):
    rck = ["rck_bulk", "rck_mc"]
    d = main[main.m.isin(rck)]
    fig, axes = plt.subplots(1, 3, figsize=(7.2, 2.5), gridspec_kw={"width_ratios": [1, 1.25, 1]})

    ax = axes[0]
    ax.plot([0, 1], [0, 1], color=INK2, linewidth=0.8, linestyle="--", zorder=1)
    for m in rck:
        s = d[d.m == m]
        for c in CLASS_ORDER:
            q = s[s["class"] == c]
            ax.scatter(q.blind_jseg_best_mean, q.phased_jseg_best_mean, s=24, marker=MARKER[c],
                       color=COLOR[m], edgecolors="white", linewidths=0.8, alpha=0.9, zorder=3)
    ax.set_xlim(0.3, 1.03)
    ax.set_ylim(0.3, 1.03)
    ax.set_aspect("equal")
    ax.set_xlabel("J_seg, haplotype-blind")
    ax.set_ylabel("J_seg, phased")
    ax.set_title("a  allele-specific structure", loc="left", color=INK)
    n11 = int(((d.blind_jseg_best_mean > 0.999) & (d.phased_jseg_best_mean > 0.999)).sum())
    ax.text(0.98, 0.03, "{} of {} runs at (1, 1)".format(n11, len(d)), transform=ax.transAxes,
            fontsize=6, color=INK2, ha="right")
    ax.legend(handles=[Line2D([], [], marker="s", linestyle="", color=COLOR[m], label=SHORT[m])
                       for m in rck], loc="upper left", frameon=False, handletextpad=0.2, fontsize=6)

    ax = axes[1]
    s = d[d.m == "rck_mc"].sort_values(["class", "S", "name"], ascending=[False, True, True])
    x = np.arange(len(s))
    for xi, a, w in zip(x, s.errors_absolute, s.errors_within):
        if a != w:
            ax.plot([xi, xi], [w, a], color=MUTED, linewidth=1, zorder=1)
    ax.scatter(x, s.errors_absolute, s=20, facecolors="white", edgecolors=COLOR["rck_mc"],
               linewidths=1.2, zorder=3, label="absolute A/B")
    ax.scatter(x, s.errors_within, s=20, color=COLOR["rck_mc"], edgecolors="white",
               linewidths=0.6, zorder=4, label="within phase blocks")
    groups = list(zip(s["class"], s.S))
    centers, names = [], []
    for g in dict.fromkeys(groups):
        idx = [i for i, h in enumerate(groups) if h == g]
        centers.append(np.mean(idx))
        names.append("{}\nS={}".format(CLASS_LABEL[g[0]], int(g[1])))
        if idx[0] > 0:
            ax.axvline(idx[0] - 0.5, color=GRID, linewidth=1.2, zorder=0)
    ax.grid(axis="x", visible=False)
    ax.set_xticks(centers)
    ax.set_xticklabels(names, fontsize=6)
    ax.tick_params(axis="x", length=0)
    ax.set_xlabel("RCK multi-clone* cases")
    ax.set_ylabel("missed + false species")
    ax.set_title("b  phase-block identifiability", loc="left", color=INK, pad=18)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0.99), ncol=2, frameon=False,
              borderaxespad=0, handletextpad=0.2, columnspacing=0.8)

    ax = axes[2]
    rows = []
    width = 0.36
    for j, m in enumerate(rck):
        for k, c in enumerate(CLASS_ORDER):
            q = d[(d.m == m) & (d["class"] == c)]
            for kind, col, face in [("blind", "blind_exact_cycles", COLOR[m]),
                                    ("phased", "phased_exact_cycles", "white")]:
                rate = q[col].sum() / q.n_truth.sum()
                rows.append({"m": m, "class": c, "kind": kind, "rate": rate})
            b = [r for r in rows if r["m"] == m and r["class"] == c]
            xpos = k + (j - 0.5) * width
            ax.bar(xpos, b[0]["rate"], width=width - 0.04, color=COLOR[m], zorder=2)
            ax.scatter([xpos], [b[1]["rate"]], marker="_", s=120, color=INK, linewidths=1.5,
                       zorder=3)
    ax.set_xticks(range(len(CLASS_ORDER)))
    ax.set_xticklabels([CLASS_LABEL[c] for c in CLASS_ORDER])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("exact cycles / true species")
    ax.set_title("c  exact cycles", loc="left", color=INK, pad=30)
    ax.legend(handles=[Line2D([], [], marker="s", linestyle="", color=COLOR[m], label=SHORT[m] +
                              " (blind)") for m in rck] +
              [Line2D([], [], marker="_", linestyle="", color=INK, markersize=10, label="phased")],
              loc="lower left", bbox_to_anchor=(0, 0.99), ncol=2, frameon=False, fontsize=6,
              borderaxespad=0, handletextpad=0.2, columnspacing=0.6)
    fig.tight_layout(w_pad=1.2)
    fig.text(0.01, -0.02, "* given the true clone profiles. CoRAL and CycleExtractor report no "
             "alleles, so they have no phased scores.", fontsize=6.5, color=INK2, va="top")
    save(fig, out, "fig6_allele_specific",
         d[["name", "class", "S", "m", "blind_jseg_best_mean", "phased_jseg_best_mean",
            "n_missed", "n_false", "n_missed_absolute", "n_false_absolute", "blind_exact_cycles",
            "phased_exact_cycles", "n_truth"]])
    pd.DataFrame(rows).to_csv(os.path.join(out, "figure_data", "fig6c_exact_rates.tsv"),
                              sep="\t", index=False)


# -- table 1 ----------------------------------------------------------------------------------

T1_COLS = [("cn_recall", "ecDNA CN recall", "up"),
           ("cn_precision", "ecDNA CN precision", "up"),
           ("proportion_tv", "proportion TV", "down"),
           ("species_detection_strict_blind_f1", "species F1 (blind)", "up"),
           ("blind_jseg_best_mean", "J_seg (blind)", "up"),
           ("blind_jadj_best_mean", "J_adj (blind)", "up")]
VARIANT_LABEL = {"rck_bulk": "RCK-ecDNA bulk", "rck_mc": "RCK-ecDNA multi-clone*",
                 "coral:walks_kept": "CoRAL† (linear walks kept)",
                 "coral:walks_dropped": "CoRAL† (linear walks dropped)",
                 "ce:walks_kept": "CycleExtractor† (linear walks kept)",
                 "ce:walks_dropped": "CycleExtractor† (linear walks dropped)"}


def table1(t, out):
    rows = []
    for scope, sub in [("all", t)] + [(CLASS_LABEL[c], t[t["class"] == c]) for c in CLASS_ORDER]:
        for v in VARIANT_ORDER:
            d = sub[sub.variant == v]
            r = {"cases": scope, "method": VARIANT_LABEL[v], "n": len(d)}
            for col, name, _ in T1_COLS:
                r[name + " mean"] = d[col].mean()
                r[name + " sd"] = d[col].std()
            r["exact-cycle rate"] = d.blind_exact_cycles.sum() / d.n_truth.sum()
            r["strict rule"] = "/".join(sorted(d.strict_rule.dropna().unique()))
            rows.append(r)
    tab = pd.DataFrame(rows)
    tab.to_csv(os.path.join(out, "table1.tsv"), sep="\t", index=False, float_format="%.4f")

    heads = ["method"] + [n for _, n, _ in T1_COLS] + ["exact-cycle rate"]
    lines = []
    for scope in tab.cases.unique():
        s = tab[tab.cases == scope]
        lines += ["", "**{} cases** (n = {} per method)".format(scope, int(s.n.iloc[0])), "",
                  "| " + " | ".join(heads) + " |", "|" + "---|" * len(heads)]
        best = {}
        for _, n, direction in T1_COLS:
            col = s[n + " mean"]
            best[n] = col.max() if direction == "up" else col.min()
        for _, r in s.iterrows():
            cells = [r["method"]]
            for _, n, _ in T1_COLS:
                cell = "{:.2f} ± {:.2f}".format(r[n + " mean"], r[n + " sd"])
                if np.isclose(r[n + " mean"], best[n]):
                    cell = "**" + cell + "**"
                cells.append(cell)
            cells.append("{:.2f}".format(r["exact-cycle rate"]))
            lines.append("| " + " | ".join(cells) + " |")
    lines += ["", "Mean ± SD over cases; best mean per column in bold. " + FOOT,
              "Strict cycles: RCK `status`, CoRAL / CE `cn-frac:0.1`. Exact-cycle rate = exact "
              "cycles / true species, summed over cases."]
    with open(os.path.join(out, "table1.md"), "w") as fh:
        fh.write("# Table 1: bench_IW_intpsi, per method variant\n" + "\n".join(lines) + "\n")

    tests = []
    cols = T1_COLS + [("blind_exact_cycles", "exact cycles", "up")]
    for ours in ["rck_bulk", "rck_mc"]:
        for base in VARIANT_ORDER[2:]:
            a = t[t.variant == ours].set_index("name")
            b = t[t.variant == base].set_index("name")
            for col, name, direction in cols:
                x, y = a[col].align(b[col], join="inner")
                ok = x.notna() & y.notna()
                diff = (x - y)[ok]
                try:
                    p = wilcoxon(x[ok], y[ok], zero_method="wilcox").pvalue if (diff != 0).any() else 1.0
                except ValueError:
                    p = np.nan
                better = (diff > 0) if direction == "up" else (diff < 0)
                tests.append({"ours": ours, "baseline": base, "metric": name, "n_pairs": int(ok.sum()),
                              "median_diff": diff.median(), "ours_better": int(better.sum()),
                              "ties": int((diff == 0).sum()), "wilcoxon_p": p})
    tt = pd.DataFrame(tests)
    # Holm correction within each (ours, baseline) family
    tt["p_holm"] = np.nan
    for _, idx in tt.groupby(["ours", "baseline"]).groups.items():
        ps = tt.loc[idx, "wilcoxon_p"]
        order = ps.sort_values().index
        m = ps.notna().sum()
        running = 0.0
        for rank, i in enumerate(order):
            if pd.isna(ps[i]):
                continue
            running = max(running, min(1.0, (m - rank) * ps[i]))
            tt.loc[i, "p_holm"] = running
    tt.to_csv(os.path.join(out, "table1_paired_tests.tsv"), sep="\t", index=False, float_format="%.4g")
    return tab, tt


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--table", default=os.path.join(HERE, "bench_iw_intpsi_table.tsv"))
    ap.add_argument("--proportions", default=os.path.join(HERE, "bench_iw_intpsi_proportions.tsv"))
    ap.add_argument("--baseline-variant", choices=["walks_kept", "walks_dropped"], default="walks_kept",
                    help="which CoRAL / CE scoring the main figures show (Table 1 lists both)")
    ap.add_argument("--out", default=os.path.join(HERE, "figures_intpsi"))
    args = ap.parse_args()
    os.makedirs(os.path.join(args.out, "figure_data"), exist_ok=True)

    t, main_df, n_empty = load(args.table, args.baseline_variant)
    props = pd.read_csv(args.proportions, sep="\t")
    assert len(t) == 144 and t.run_status.eq("ok").all(), "expected 144 scored rows"
    fig2(main_df, args.out)
    fig3(main_df, args.out)
    fig4(main_df, args.out)
    fig5(main_df, props, args.baseline_variant, args.out)
    fig6(main_df, args.out)
    tab, tt = table1(t, args.out)
    print("wrote figures 2-6 and Table 1 to {} (baselines: {}; {} runs with no cycles counted as 0 "
          "for J / exact cycles)".format(args.out, args.baseline_variant, n_empty))


if __name__ == "__main__":
    main()

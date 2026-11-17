#!/usr/bin/env python
"""Figure 10 and Table 4: the multi-clone route on RCK's degraded inputs (PLAN.md Part B; RCK 2020 Fig. 2,
Supplemental Figs. S2-S3).

Reads partB_table.tsv (score_partB.py collect); writes to figures/ here:
    fig10_partB.{pdf,png}      per input variant (labelled reference, chil, chil_avg, chil_err), single-allele
                               cases vs cases where both alleles are amplified: labelled novel-adjacency FNR
                               and FDR (RCK Fig. 2A), Eq. 11 summed over clones on the ecDNA bases (Fig. 2B),
                               allele-specific ecDNA copy-number recall (items 2-3)
    table4.tsv / table4.md     means per variant and case group (Eq. 11: median), with input use, RCK's FPR
                               over spurious inputs, copy-number state agreement, IS violations and run times
    figure_data/fig10_partB.tsv

Scores within phase blocks (score_partB.py docstring: the common frame); Eq. 11 also in the variant's own
fragments in Table 4. Style from bench_IW/make_figures.py (imported, not edited): one method, so one colour;
filled markers = input with junction labels, open = CHIL (no labels).

    python make_partB_figures.py
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_IW"))
import make_figures as mf  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

VARIANTS = ["labelled", "chil", "chil_avg", "chil_err"]
VLABEL = {"labelled": "labelled\n(reference)", "chil": "CHIL", "chil_avg": "CHIL +\n50 kb avg.", "chil_err": "CHIL +\nerrors"}
GROUPS = [("single", "single-allele cases"), ("mixed", "both alleles amplified")]
COLOR = mf.COLOR["rck_mc"]
METRICS = [("adj_phased_fnr", "labelled adjacency\nFNR", "linear"),
           ("adj_phased_fdr", "labelled adjacency\nFDR", "linear"),
           ("clone_cn_pb_dist_rel_ecdna", "Eq. 11 over clones,\necDNA bases (rel.)", "log"),
           ("allele_as_cn_recall_block", "allele-specific ecDNA\ncopy-number recall", "linear")]
FOOT10 = ("RCK-ecDNA multi-clone*, P = 0.9. labelled: the subgraph with the true clone profiles and long-read junction\n"
          "haplotype labels; CHIL: A/B flipped per fragment, junction labels removed (Aganezov & Raphael 2020); avg.: copy\n"
          "numbers averaged over 50 kb; errors: + ±1 on 5% of fragment values, ±50 bp on half the breakends, 10% spurious\n"
          "junctions. Within phase blocks. Eq. 11 log scale, 0 drawn at 0.0001.\n"
          "* given the true clone profiles (upper bound for the multi-clone route).")


def load():
    t = pd.read_csv(os.path.join(HERE, "partB_table.tsv"), sep="\t")
    return t[t["status"] == "ok"].copy()


def fig10(t, out):
    fig, axes = plt.subplots(len(METRICS), 2, figsize=(6.6, 8.0), sharey="row", gridspec_kw={"wspace": 0.06, "hspace": 0.45})
    rows = []
    for j, (g, title) in enumerate(GROUPS):
        n = t[t.group == g]["name"].nunique()
        for i, (metric, ylabel, scale) in enumerate(METRICS):
            ax = axes[i, j]
            for x, v in enumerate(VARIANTS):
                d = t[(t.group == g) & (t.variant == v)]
                values = d[metric].values.astype(float)
                if scale == "log":
                    values = np.clip(values, 1e-4, None)
                xs = x + mf.jitter(len(d), seed=x)
                for c in mf.CLASS_ORDER:
                    sel = d["class"].values == c
                    face = COLOR if v == "labelled" else "white"
                    ax.scatter(xs[sel], values[sel], s=24, marker=mf.MARKER[c], facecolors=face, edgecolors=COLOR,
                               linewidths=1.1, zorder=3)
                if scale == "log":
                    val = float(np.median(values)) if len(values) else np.nan
                    ax.plot([x - 0.3, x + 0.3], [val, val], color=mf.INK, linewidth=1.6, zorder=4)
                else:
                    val = mf.mean_bar(ax, x, values) if len(values) else np.nan
                rows.append({"group": g, "metric": metric, "variant": v, "n": len(d),
                             "summary": "median" if scale == "log" else "mean", "value": val})
            if scale == "log":
                ax.set_yscale("log")
                ax.set_ylim(6e-5, 3)
            else:
                ax.set_ylim(-0.03, 1.05)
            ax.set_xticks(range(len(VARIANTS)))
            ax.set_xticklabels([VLABEL[v] for v in VARIANTS], fontsize=6.2)
            ax.set_xlim(-0.6, len(VARIANTS) - 0.4)
            ax.grid(axis="x", visible=False)
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 0:
                ax.set_title("{t} (n = {n})".format(t=title, n=n))
    handles = [Line2D([], [], marker="s", linestyle="", markersize=6, color=COLOR, label="with junction labels"),
               Line2D([], [], marker="s", linestyle="", markersize=6, markerfacecolor="white", markeredgecolor=COLOR,
                      label="CHIL (no junction labels)")]
    handles += [Line2D([], [], marker=mf.MARKER[c], linestyle="", markersize=6, markerfacecolor="white",
                       markeredgecolor=mf.INK2, label=mf.CLASS_LABEL[c]) for c in mf.CLASS_ORDER]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.06), ncol=4, frameon=False, fontsize=6.5,
               handletextpad=0.3)
    fig.text(0.06, 0.03, FOOT10, fontsize=6.0, color=mf.INK2, ha="left", va="top")
    mf.save(fig, out, "fig10_partB", pd.DataFrame(rows))


TABLE4 = [("adj_blind_fnr", "adj. FNR (unlabelled)", "{:.2f}", "mean"),
          ("adj_phased_fnr", "labelled adj. FNR", "{:.2f}", "mean"),
          ("adj_phased_fdr", "labelled adj. FDR", "{:.2f}", "mean"),
          ("input_input_used", "input adj. used", "{:.2f}", "mean"),
          ("clone_cn_pb_dist_rel_ecdna", "Eq. 11 clones, ecDNA (pb, median)", "{:.3f}", "median"),
          ("clone_cn_frag_dist_rel_ecdna", "Eq. 11 clones, ecDNA (fragments, median)", "{:.3f}", "median"),
          ("clone_cn_pb_state_agree_ecdna", "CN state agreement, ecDNA", "{:.3f}", "mean"),
          ("allele_as_cn_recall_block", "allele CN recall", "{:.2f}", "mean"),
          ("is_is_hom", "homologous", "{:.2f}", "mean"),
          ("is_is_recip", "homologous-reciprocal", "{:.2f}", "mean"),
          ("run_n_cycles", "cycles", "{:.1f}", "mean"),
          ("run_runtime_rck_s", "solve s (median)", "{:.1f}", "median")]


def table4(t, out):
    md = ["# Table 4: the multi-clone route on RCK's degraded inputs (bench_IW_intpsi, within phase blocks)", "",
          "Input variants as Aganezov & Raphael (2020) degraded RCK's simulated inputs (PLAN.md Part B). Means over "
          "cases unless marked. P = 0.9.", ""]
    rows = []
    for g, title in GROUPS + [("all", "all cases")]:
        sub = t if g == "all" else t[t.group == g]
        md += ["**{t}** (n = {n} cases)".format(t=title, n=sub["name"].nunique()), "",
               "| input | runs (optimal) | " + " | ".join(lbl for _, lbl, _, _ in TABLE4) + " |",
               "|---|---|" + "---|" * len(TABLE4)]
        for v in VARIANTS:
            d = sub[sub.variant == v]
            row = {"group": g, "variant": v, "n": len(d), "optimal": int((d["run_status"] == 2).sum())}
            cells = []
            for col, _, fmt, how in TABLE4:
                val = d[col].median() if how == "median" else d[col].mean()
                row[col] = val
                cells.append(fmt.format(val) if pd.notna(val) else "–")
            rows.append(row)
            md.append("| {v} | {n} ({o}) | {c} |".format(v=VLABEL[v].replace("\n", " "), n=len(d), o=row["optimal"],
                                                         c=" | ".join(cells)))
        md.append("")
    err = t[t.variant == "chil_err"]
    used, spurious = int(err["input_n_spurious_used"].sum()), int(err["input_n_spurious"].sum())
    md += ["RCK's FPR over the input (spurious input adjacencies used / spurious input adjacencies), CHIL + errors, "
           "pooled over cases: {u} / {s}{r}.".format(u=used, s=spurious, r=" = {:.2f}".format(used / spurious) if spurious else ""),
           "", "pb: phase blocks (the common frame); fragments: the variant's own input fragments (RCK's Eq. 11 as "
           "written). n = 6 mixed cases: descriptive only.", "",
           "* given the true clone profiles (upper bound for the multi-clone route)."]
    pd.DataFrame(rows).to_csv(os.path.join(out, "table4.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "table4.md"), "w") as handle:
        handle.write("\n".join(md) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--out", default=os.path.join(HERE, "figures"))
    a = p.parse_args()
    os.makedirs(os.path.join(a.out, "figure_data"), exist_ok=True)
    t = load()
    fig10(t, a.out)
    table4(t, a.out)
    print("fig10 and table4 ->", a.out)


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Figure 9 and Table 3: RCK-style evaluation of bench_IW_intpsi (PLAN.md Part A).

Reads rck_style_table.tsv (rck_style.py collect); writes to figures/ here:
    fig9_rck_style.{pdf,png}   the scores Aganezov & Raphael (2020) used for RCK, per method and phasing,
                               single-allele cases vs cases where both alleles are amplified:
                               haplotype-specific novel-adjacency FNR and FDR (their Fig. 2A),
                               the Eq. 11 allele-specific copy-number distance (Fig. 2B), and
                               homologous generalized-infinite-sites violations (Fig. 4)
    table3.tsv / table3.md     mean (Eq. 11: median) per method and phasing, by case group, with the truth's
                               own violations
    figure_data/fig9_rck_style.tsv

Methods: RCK MC (rck_mc), CoRAL and CE only (Yubi, 2026-11-21); rck_bulk is not shown.

Every score is within phase blocks (RCK's "up to haplotype symmetries"; one A/B swap per block). Style,
colours and markers are bench_IW/make_figures.py's (imported, not edited); post-hoc phasing is drawn with
open markers as in Figure 7.

    python make_rck_style_figures.py [--baseline-variant walks_kept]
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_IW"))
import make_figures as mf  # noqa: E402  (style, colours, helpers)
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

PHASE_LABEL = {"native": "own", "baf_major": "BAF", "truth_major": "oracle"}
GROUPS = [("single", "single-allele cases"), ("mixed", "both alleles amplified")]
FOOT9 = ("Within phase blocks. own: the method's haplotypes; BAF / oracle: post-hoc phasing (Fig. 7). Adjacencies match within\n"
         "100 bp and must carry the true haplotype pair. Eq. 11: sum over phase blocks of the A/B-swap minimum of\n"
         "(|c_A - t_A| + |c_B - t_B| + unphased) x bp, over the true copies x bp; log scale, 0 drawn at 0.001.\n"
         "Violations: extremities used on both homologs, plus reference-neighbour extremity pairs used on opposite\n"
         "homologs (RCK's homologous and homologous-reciprocal exclusivity); dashed line: the truth's own mean.")


def load(baseline_variant):
    t = pd.read_csv(os.path.join(HERE, "rck_style_table.tsv"), sep="\t")
    t = t[t["status"] == "ok"].copy()
    t["group"] = np.where(t["mixed_bp_frac"] > 0, "mixed", "single")
    t["is_homologous"] = t["is_is_hom"] + t["is_is_recip"]
    truth = t[t.method == "truth"].copy()
    keep = {"rck_mc": "rck_mc", "coral:" + baseline_variant: "coral", "ce:" + baseline_variant: "ce"}
    t = t[t.method_variant.isin(keep)].copy()
    t["m"] = t.method_variant.map(keep)
    return t, truth


def columns():
    """(method, phasing) in plot order."""
    out = []
    for m in mf.METHODS:
        for p in (["native"] if m.startswith("rck") else []) + ["baf_major", "truth_major"]:
            out.append((m, p))
    return out


METRICS = [("adj_phased_fnr", "labelled adjacency\nFNR", "linear"),
           ("adj_phased_fdr", "labelled adjacency\nFDR", "linear"),
           ("cn_dist_rel", "Eq. 11 CN distance\n(relative)", "log"),
           ("is_homologous", "homologous IS\nviolations per case", "linear")]


def fig9(t, truth, out):
    cols = columns()
    fig, axes = plt.subplots(len(METRICS), 2, figsize=(7.4, 8.2), sharey="row",
                             gridspec_kw={"wspace": 0.06, "hspace": 0.4})
    rows = []
    for j, (g, title) in enumerate(GROUPS):
        n = t[t.group == g]["name"].nunique()
        for i, (metric, ylabel, scale) in enumerate(METRICS):
            ax = axes[i, j]
            for x, (m, p) in enumerate(cols):
                d = t[(t.group == g) & (t.m == m) & (t.phasing == p)]
                values = d[metric].values.astype(float)
                if scale == "log":
                    values = np.clip(values, 1e-3, None)  # exact answers (0) drawn at the axis floor
                if p == "native":
                    mf.dots(ax, x, values, d["class"].values, mf.COLOR[m], seed=x)
                else:
                    xs = x + mf.jitter(len(d), seed=x)
                    for c in mf.CLASS_ORDER:
                        sel = d["class"].values == c
                        ax.scatter(xs[sel], values[sel], s=24, marker=mf.MARKER[c], facecolors="white",
                                   edgecolors=mf.COLOR[m], linewidths=1.1, zorder=3)
                if scale == "log":
                    v = float(np.median(values)) if len(values) else np.nan
                    ax.plot([x - 0.3, x + 0.3], [v, v], color=mf.INK, linewidth=1.6, zorder=4)
                else:
                    v = mf.mean_bar(ax, x, values)
                rows.append({"group": g, "metric": metric, "method": m, "phasing": p, "n": len(d),
                             "summary": "median" if scale == "log" else "mean", "value": v})
            if metric == "is_homologous":
                ref = truth[truth.group == g][metric].mean()
                ax.axhline(ref, color=mf.MUTED, linewidth=0.9, linestyle="--", zorder=1)
                if j == 1:  # the left column's label would run into the right panel
                    ax.text(1.0, ref, " truth", transform=ax.get_yaxis_transform(), fontsize=6, color=mf.INK2, ha="left",
                            va="center")
                rows.append({"group": g, "metric": metric, "method": "truth", "phasing": "truth",
                             "n": int((truth.group == g).sum()), "summary": "mean", "value": ref})
            if scale == "log":
                ax.set_yscale("log")
                ax.set_ylim(7e-4, 100)
            elif metric != "is_homologous":
                ax.set_ylim(-0.03, 1.05)
            ax.set_xticks(range(len(cols)))
            ax.set_xticklabels([PHASE_LABEL[p] for _, p in cols], fontsize=6, rotation=40, ha="right",
                               rotation_mode="anchor")
            ax.grid(axis="x", visible=False)
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 0:
                ax.set_title("{t} (n = {n})".format(t=title, n=n))
            starts = {}
            for x, (m, _) in enumerate(cols):
                starts.setdefault(m, []).append(x)
            if i == len(METRICS) - 1:
                for m, xs in starts.items():
                    ax.text(np.mean(xs), -0.3, mf.SHORT[m], transform=ax.get_xaxis_transform(), ha="center", va="top",
                            fontsize=7, color=mf.INK)
            for m, xs in list(starts.items())[1:]:
                ax.axvline(xs[0] - 0.5, color=mf.GRID, linewidth=0.8, zorder=1)
    handles = [Line2D([], [], marker="s", linestyle="", markersize=6, color=mf.INK2, label="filled: own haplotypes"),
               Line2D([], [], marker="s", linestyle="", markersize=6, markerfacecolor="white", markeredgecolor=mf.INK2,
                      label="open: post-hoc phasing")]
    handles += [Line2D([], [], marker=mf.MARKER[c], linestyle="", markersize=6, markerfacecolor="white",
                       markeredgecolor=mf.INK2, label=mf.CLASS_LABEL[c]) for c in mf.CLASS_ORDER]
    fig.legend(handles=handles, loc="upper center", bbox_to_anchor=(0.5, 0.045), ncol=4, frameon=False, fontsize=6.5,
               handletextpad=0.3)
    fig.text(0.06, 0.012, FOOT9, fontsize=6.0, color=mf.INK2, ha="left", va="top")
    mf.save(fig, out, "fig9_rck_style", pd.DataFrame(rows))


TABLE3 = [("adj_blind_fnr", "adjacency FNR (unlabelled)", "{:.2f}", "mean"),
          ("adj_phased_fnr", "labelled adjacency FNR", "{:.2f}", "mean"),
          ("adj_phased_fdr", "labelled adjacency FDR", "{:.2f}", "mean"),
          ("cn_dist_rel", "Eq. 11 CN distance (rel., median)", "{:.3f}", "median"),
          ("is_is_ext", "extremity reuse", "{:.2f}", "mean"),
          ("is_is_hom", "homologous", "{:.2f}", "mean"),
          ("is_is_recip", "homologous-reciprocal", "{:.2f}", "mean")]


def table3(t, truth, out):
    md = ["# Table 3: RCK-style evaluation (bench_IW_intpsi, within phase blocks)", "",
          "Scores of Aganezov & Raphael (2020): haplotype-specific novel-adjacency rates (their Fig. 2A), the Eq. 11 "
          "allele-specific copy-number distance (Fig. 2B), and generalized infinite-sites violations per case (Fig. 4; "
          "the last three columns).", ""]
    rows = []
    for g, title in GROUPS:
        sub = t[t.group == g]
        md += ["**{t}** (n = {n} cases)".format(t=title, n=sub["name"].nunique()), "",
               "| method | phasing | " + " | ".join(lbl for _, lbl, _, _ in TABLE3) + " |",
               "|---|---|" + "---|" * len(TABLE3)]
        tr = truth[truth.group == g]
        cells = ["–"] * 4 + [fmt.format(tr[col].mean()) for col, _, fmt, _ in TABLE3[4:]]
        md.append("| truth | – | " + " | ".join(cells) + " |")
        rows.append({"group": g, "method": "truth", "phasing": "truth", "n": len(tr),
                     **{col: tr[col].mean() for col, _, _, _ in TABLE3[4:]}})
        for m, p in columns():
            d = sub[(sub.m == m) & (sub.phasing == p)]
            row = {"group": g, "method": m, "phasing": p, "n": len(d)}
            cells = []
            for col, _, fmt, how in TABLE3:
                v = d[col].median() if how == "median" else d[col].mean()
                row[col] = v
                cells.append(fmt.format(v))
            rows.append(row)
            md.append("| {m} | {p} | {c} |".format(m=mf.SHORT[m], p=PHASE_LABEL[p], c=" | ".join(cells)))
        md.append("")
    md += ["Means over cases (Eq. 11: median, as CoRAL / CE's over-counted copies make the mean uninformative). own: the "
           "method's haplotypes; BAF / oracle: post-hoc phasing (Fig. 7). The input subgraphs contain no spurious "
           "candidate junctions (0 of 9 median per case), so RCK's input-level FPR is undefined here (PLAN.md Part B "
           "adds spurious candidates). n = 6 mixed cases: descriptive only."]
    pd.DataFrame(rows).to_csv(os.path.join(out, "table3.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "table3.md"), "w") as handle:
        handle.write("\n".join(md) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--baseline-variant", default="walks_kept", choices=["walks_kept", "walks_dropped"])
    p.add_argument("--out", default=os.path.join(HERE, "figures"))
    a = p.parse_args()
    os.makedirs(os.path.join(a.out, "figure_data"), exist_ok=True)
    t, truth = load(a.baseline_variant)
    fig9(t, truth, a.out)
    table3(t, truth, a.out)
    print("fig9 and table3 ->", a.out)


if __name__ == "__main__":
    main()

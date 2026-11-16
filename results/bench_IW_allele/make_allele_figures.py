#!/usr/bin/env python
"""Figures 7-8 and Table 2: allele-specific results on bench_IW_intpsi (PLAN.md items 2 and 3).

Reads allele_table.tsv (allele_bench.py collect) and the per-run allele tables under
$BENCH_IW_ALLELE_OUT; writes to figures/ here:
    fig7_posthoc_phasing.{pdf,png}   phasing inside the decomposition (RCK) vs after it (post-hoc rules),
                                     allele-specific ecDNA copy-number recall and phased J_adj,
                                     single-allele cases vs cases where both alleles are amplified
    fig8_allele_cn.{pdf,png}         true vs inferred ecDNA copy number per allele and segment, pooled
    table2.tsv / table2.md           mean +- SD per method and phasing, by case group
    figure_data/*.tsv                the numbers each figure plots

Every score here is within phase blocks (one A/B swap per block, the frame rck_mc is scored in;
docs/pipeline/evaluation.md 3.10). Absolute scores are in allele_table.tsv (columns without _block).
Style, colours and markers are bench_IW/make_figures.py's (imported, not edited).

    python make_allele_figures.py [--baseline-variant walks_kept]
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
FOOT7 = ("own: the method's haplotypes. BAF / oracle: post-hoc, each segment given the major allele of the observed "
         "CN_A / CN_B (BAF) or of the true ecDNA copies (oracle, the best any one-allele-per-segment rule can do).\n"
         + mf.FOOT)


def load(baseline_variant):
    a = pd.read_csv(os.path.join(HERE, "allele_table.tsv"), sep="\t")
    a = a[a["status"] == "ok"].copy()
    mixed = a[(a.method_variant == "rck_bulk") & (a.phasing == "native")].set_index("name")["mixed_bp_frac"]
    a["group"] = np.where(a["name"].map(mixed) > 0, "mixed", "single")
    keep = {"rck_bulk": "rck_bulk", "rck_mc": "rck_mc", "coral:" + baseline_variant: "coral", "ce:" + baseline_variant: "ce"}
    a = a[a.method_variant.isin(keep)].copy()
    a["m"] = a.method_variant.map(keep)
    for col in ("phased_jseg_best_mean_block", "phased_jadj_best_mean_block", "phased_exact_cycles_block"):
        a[col] = a[col].fillna(0.0)  # no cycle at all: nothing matches (make_figures.load)
    return a


def columns():
    """(method, phasing) in plot order."""
    out = []
    for m in mf.METHODS:
        for p in (["native"] if m.startswith("rck") else []) + ["baf_major", "truth_major"]:
            out.append((m, p))
    return out


def fig7(a, out):
    metrics = [("as_cn_recall_block", "allele-specific ecDNA\ncopy-number recall"),
               ("phased_jadj_best_mean_block", "phased adjacency\nJaccard (J_adj)")]
    cols = columns()
    fig, axes = plt.subplots(2, 2, figsize=(8.4, 4.8), sharey="row", gridspec_kw={"wspace": 0.06, "hspace": 0.35})
    rows = []
    for j, (g, title) in enumerate(GROUPS):
        n = a[a.group == g]["name"].nunique()
        for i, (metric, ylabel) in enumerate(metrics):
            ax = axes[i, j]
            for x, (m, p) in enumerate(cols):
                d = a[(a.group == g) & (a.m == m) & (a.phasing == p)]
                if p == "native":
                    mf.dots(ax, x, d[metric].values, d["class"].values, mf.COLOR[m], seed=x)
                else:  # post-hoc: open markers in the method's colour
                    xs = x + mf.jitter(len(d), seed=x)
                    for c in mf.CLASS_ORDER:
                        sel = d["class"].values == c
                        ax.scatter(xs[sel], d[metric].values[sel], s=24, marker=mf.MARKER[c], facecolors="white",
                                   edgecolors=mf.COLOR[m], linewidths=1.1, zorder=3)
                v = mf.mean_bar(ax, x, d[metric].values)
                rows.append({"group": g, "metric": metric, "method": m, "phasing": p, "n": len(d), "mean": v})
            ax.set_xticks(range(len(cols)))
            ax.set_xticklabels([PHASE_LABEL[p] for _, p in cols], fontsize=6, rotation=40, ha="right", rotation_mode="anchor")
            ax.set_ylim(-0.03, 1.05)
            ax.grid(axis="x", visible=False)
            if j == 0:
                ax.set_ylabel(ylabel)
            if i == 0:
                ax.set_title("{t} (n = {n})".format(t=title, n=n))
            # method names under each group of columns
            starts = {}
            for x, (m, _) in enumerate(cols):
                starts.setdefault(m, []).append(x)
            if i == 1:
                for m, xs in starts.items():
                    ax.text(np.mean(xs), -0.32, mf.SHORT[m], transform=ax.get_xaxis_transform(), ha="center", va="top",
                            fontsize=7, color=mf.INK)
            for m, xs in list(starts.items())[1:]:
                ax.axvline(xs[0] - 0.5, color=mf.GRID, linewidth=0.8, zorder=1)
    handles = [Line2D([], [], marker="s", linestyle="", markersize=6, color=mf.INK2, label="filled: own haplotypes"),
               Line2D([], [], marker="s", linestyle="", markersize=6, markerfacecolor="white", markeredgecolor=mf.INK2,
                      label="open: post-hoc phasing")]
    handles += [Line2D([], [], marker=mf.MARKER[c], linestyle="", markersize=6, markerfacecolor="white",
                       markeredgecolor=mf.INK2, label=mf.CLASS_LABEL[c]) for c in mf.CLASS_ORDER]
    axes[0, 0].legend(handles=handles, loc="lower left", frameon=False, fontsize=6.5, handletextpad=0.3)
    fig.text(0.01, -0.1, FOOT7, fontsize=6.2, color=mf.INK2, ha="left", va="top")
    mf.save(fig, out, "fig7_posthoc_phasing", pd.DataFrame(rows))


def within_block_tables(a, phasing_of):
    """Per (case, method): its allele table swapped within phase blocks (as allele_metrics.score does)."""
    import allele_bench as ab
    from rck.ecdna import allele_metrics as am, evaluate as ev

    cases = {c["name"]: c for c in ab.bi.cases()}
    frames = []
    for r in a.drop_duplicates(["name", "m"]).itertuples():
        p = phasing_of[r.m]
        row = a[(a.name == r.name) & (a.m == r.m) & (a.phasing == p)].iloc[0]
        case = cases[r.name]
        method, _, variant = row.method_variant.partition(":")
        t = pd.read_csv(os.path.join(ab.task_dir(case, method, variant, p), "allele", "allele_segments.tsv"), sep="\t")
        t, _ = am.within_block(t, ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv")))
        frames.append(t.assign(name=r.name, m=r.m, phasing=p, group=row.group))
    return pd.concat(frames, ignore_index=True)


def fig8(a, out):
    phasing_of = {"rck_bulk": "native", "rck_mc": "native", "coral": "baf_major", "ce": "baf_major"}
    t = within_block_tables(a, phasing_of)
    t = t[(t["end"] - t["start"] + 1) >= 1000]  # 1-bp slivers from breakpoint offsets add nothing
    mixed_base = (t.t_A > 0) & (t.t_B > 0)
    long = pd.concat([t.assign(allele=h, truth=t["t_" + h], inferred=t["c_" + h], mixed_base=mixed_base) for h in ("A", "B")])
    long = long[(long.truth > 0) | (long.inferred > 0)]
    fig, axes = plt.subplots(1, 4, figsize=(7.4, 2.5), sharex=True, sharey=True, gridspec_kw={"wspace": 0.08})
    lim = max(long.truth.max(), long.inferred.max()) * 1.15
    for ax, m in zip(axes, mf.METHODS):
        d = long[long.m == m]
        for mixed, face in ((False, mf.COLOR[m]), (True, "white")):
            s = d[d.mixed_base == mixed]
            ax.scatter(s.truth, s.inferred, s=12, facecolors=face, edgecolors=mf.COLOR[m], linewidths=0.8,
                       alpha=0.85, zorder=3 if mixed else 2)
        ax.plot([0, lim], [0, lim], color=mf.MUTED, linewidth=0.8, linestyle="--", zorder=1)
        ax.set_xscale("symlog", linthresh=2)
        ax.set_yscale("symlog", linthresh=2)
        ax.set_xlim(-0.2, lim)
        ax.set_ylim(-0.2, lim)
        rec = a[(a.m == m) & (a.phasing == phasing_of[m])]
        mix = rec[rec.group == "mixed"]
        ax.set_title("{s} ({p})\n\n".format(s=mf.SHORT[m], p=PHASE_LABEL[phasing_of[m]]), fontsize=8)
        ax.text(0.5, 1.02, "both-allele cases:\nrecall {r:.2f}, BAF error {b:.2f}".format(
            b=mix.baf_mae_block.mean(), r=mix.as_cn_recall_block.mean()), transform=ax.transAxes, ha="center", va="bottom",
            fontsize=6, color=mf.INK2)
        ax.set_xticks([0, 1, 10, 100])
        ax.set_xticklabels(["0", "1", "10", "100"])
        ax.set_yticks([0, 1, 10, 100])
        ax.set_yticklabels(["0", "1", "10", "100"])
        ax.set_xlabel("true copies (allele)")
    axes[0].set_ylabel("inferred copies (allele)")
    handles = [Line2D([], [], marker="o", linestyle="", markersize=5, color=mf.INK2, label="one allele amplified"),
               Line2D([], [], marker="o", linestyle="", markersize=5, markerfacecolor="white", markeredgecolor=mf.INK2,
                      label="both alleles amplified")]
    fig.legend(handles=handles, loc="lower center", bbox_to_anchor=(0.5, -0.2), ncol=2, frameon=False, fontsize=6.5, handletextpad=0.2)
    fig.text(0.01, -0.27, "One point per ecDNA segment (>= 1 kb) and allele, all 24 cases; within phase blocks. "
             "CoRAL / CE phased post-hoc from the observed CN_A / CN_B.\n" + mf.FOOT, fontsize=6.2, color=mf.INK2,
             ha="left", va="top")
    mf.save(fig, out, "fig8_allele_cn", long[["name", "m", "phasing", "group", "chrom", "start", "end", "allele", "truth",
                                               "inferred", "mixed_base"]])


TABLE2 = [("as_cn_recall_block", "allele CN recall", "{:.2f}"), ("as_cn_precision_block", "allele CN precision", "{:.2f}"),
          ("as_cn_mae_block", "allele CN error (copies)", "{:.1f}"), ("baf_mae_block", "ecDNA BAF error", "{:.3f}"),
          ("major_allele_acc_block", "major allele right", "{:.2f}"), ("phased_jseg_best_mean_block", "phased J_seg", "{:.2f}"),
          ("phased_jadj_best_mean_block", "phased J_adj", "{:.2f}")]


def table2(a, out):
    rows, md = [], ["# Table 2: allele-specific ecDNA reconstruction (bench_IW_intpsi, within phase blocks)", ""]
    truth_species = {}
    import allele_bench as ab
    for c in ab.bi.cases():
        truth_species[c["name"]] = len(pd.read_csv(os.path.join(c["truth"], "truth", "ecdna.tsv"), sep="\t"))
    for g, title in GROUPS:
        sub = a[a.group == g]
        md += ["**{t}** (n = {n} cases)".format(t=title, n=sub["name"].nunique()), "",
               "| method | phasing | " + " | ".join(lbl for _, lbl, _ in TABLE2) + " | phased exact-cycle rate |",
               "|---|---|" + "---|" * (len(TABLE2) + 1)]
        for m, p in columns():
            d = sub[(sub.m == m) & (sub.phasing == p)]
            row = {"group": g, "method": m, "phasing": p, "n": len(d)}
            cells = []
            for col, lbl, fmt in TABLE2:
                row[col + "_mean"], row[col + "_sd"] = d[col].mean(), d[col].std()
                cells.append((fmt + " ± " + fmt).format(d[col].mean(), d[col].std()))
            rate = d.phased_exact_cycles_block.sum() / sum(truth_species[n] for n in d["name"])
            row["phased_exact_rate"] = rate
            rows.append(row)
            md.append("| {m} | {p} | {c} | {r:.2f} |".format(m=mf.SHORT[m], p=PHASE_LABEL[p], c=" | ".join(cells), r=rate))
        md.append("")
    md += ["Mean ± SD over cases. own: the method's haplotypes; BAF / oracle: post-hoc phasing (fig. 7 caption). "
           "Exact-cycle rate = phased exact cycles / truth species, summed over cases. n = 6 mixed cases is too few "
           "for tests; report it descriptively.", "", mf.FOOT]
    pd.DataFrame(rows).to_csv(os.path.join(out, "table2.tsv"), sep="\t", index=False)
    with open(os.path.join(out, "table2.md"), "w") as handle:
        handle.write("\n".join(md) + "\n")


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("--baseline-variant", default="walks_kept", choices=["walks_kept", "walks_dropped"])
    p.add_argument("--out", default=os.path.join(HERE, "figures"))
    a = p.parse_args()
    os.makedirs(os.path.join(a.out, "figure_data"), exist_ok=True)
    table = load(a.baseline_variant)
    fig7(table, a.out)
    fig8(table, a.out)
    table2(table, a.out)
    print("figures and table2 ->", a.out)


if __name__ == "__main__":
    main()

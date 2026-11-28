#!/usr/bin/env python
"""Manuscript figures and tables for the combined 36-case benchmark (bench_IW_intpsi + bench_IW_allelic).

The two simulation sets are reported as one benchmark with three case classes of 12 cases each:
    discrete  (bench_IW_intpsi class "simple",  S = 1, 2, 4; 4 cases each)
    nested    (bench_IW_intpsi class "hard",    S = 1, 2, 4; 4 cases each)
    sibling   (bench_IW_allelic, all 12 cases;  S = 2; 4 cases per scenario)
The allelic cases are nested by construction (class "hard" in their own table); here they are their
own class so that the per-class results stay balanced.

Reads only collected tables (no runs, no Gurobi). Nothing in bench_IW/ or bench_IW_allele/ is
written or edited: their figure functions are imported and given the combined tables.
    ../bench_IW/bench_iw_intpsi_table.tsv, ../bench_IW_allele/bench_iw_allelic_table.tsv
    ../bench_IW/bench_iw_intpsi_proportions.tsv, ../bench_IW_allele/bench_iw_allelic_proportions.tsv
    ../bench_IW_allele/allele_table{,_allelic}.tsv, allele_segments{,_allelic}.tsv
    ../bench_IW_allele/rck_style_table{,_allelic}.tsv
    phase blocks: $AAA_SUBGRAPHS/<set>/<case>/union/rck_input/rck.fragments.tsv (for Fig. 8 only)

Writes here:
    combined_*.tsv                     the merged input tables (column "class" holds the three classes)
    figures/fig2_..fig6_*, table1.*    as bench_IW/make_figures.py (Fig. 3: one column per class)
    figures/fig7_*, fig8_*, table2.*   as bench_IW_allele/make_allele_figures.py
    figures/fig9_*, table3.*           as bench_IW_allele/make_rck_style_figures.py
    figures/figure_data/*.tsv

    python make_combined_figures.py
"""
import os
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
IW = os.path.join(HERE, "..", "bench_IW")
AL = os.path.join(HERE, "..", "bench_IW_allele")
sys.path.insert(0, IW)
sys.path.insert(0, AL)
sys.path.insert(0, os.path.join(HERE, "..", ".."))  # the rck package
SUBGRAPHS = os.environ.get("AAA_SUBGRAPHS", os.path.join(
    HERE, "..", "..", "..", "amplicon-allele-analysis", "simulations", "subgraphs"))
SET_DIR = {"iw_intpsi": "bench_IW_intpsi_cnp", "iw_allelic": "bench_IW_allelic_cnp"}

import make_figures as mf  # noqa: E402
import make_allele_figures as maf  # noqa: E402
import make_rck_style_figures as mrs  # noqa: E402
import matplotlib.pyplot as plt  # noqa: E402
from matplotlib.lines import Line2D  # noqa: E402

# three classes; ordered categorical so that sorts (Fig. 6b) run discrete, nested, sibling
CLASSES = ["simple", "hard", "sibling"]
mf.CLASS_ORDER = CLASSES
mf.CLASS_LABEL = {"simple": "discrete", "hard": "nested", "sibling": "sibling"}
mf.MARKER = {"simple": "o", "hard": "^", "sibling": "D"}  # not "s": legends use squares for fill
SORT_CATS = ["sibling", "hard", "simple"]  # Fig. 6b sorts class descending

OUT = os.path.join(HERE, "figures")


def as_classes(df):
    df["class"] = pd.Categorical(df["class"], categories=SORT_CATS, ordered=True)
    return df


def combine_tables():
    t1 = pd.read_csv(os.path.join(IW, "bench_iw_intpsi_table.tsv"), sep="\t").assign(bench_set="iw_intpsi")
    t2 = pd.read_csv(os.path.join(AL, "bench_iw_allelic_table.tsv"), sep="\t").assign(bench_set="iw_allelic")
    t2["class"] = "sibling"
    t = pd.concat([t1, t2], ignore_index=True)
    assert not set(t1.name) & set(t2.name), "case names collide"
    assert len(t) == 216 and t.run_status.eq("ok").all(), "expected 36 cases x 6 variants, all ok"
    t.to_csv(os.path.join(HERE, "combined_table.tsv"), sep="\t", index=False)

    p1 = pd.read_csv(os.path.join(IW, "bench_iw_intpsi_proportions.tsv"), sep="\t")
    p2 = pd.read_csv(os.path.join(AL, "bench_iw_allelic_proportions.tsv"), sep="\t")
    p = pd.concat([p1, p2], ignore_index=True)
    p.to_csv(os.path.join(HERE, "combined_proportions.tsv"), sep="\t", index=False)

    a1 = pd.read_csv(os.path.join(AL, "allele_table.tsv"), sep="\t").assign(bench_set="iw_intpsi")
    a2 = pd.read_csv(os.path.join(AL, "allele_table_allelic.tsv"), sep="\t").assign(bench_set="iw_allelic")
    a2["class"] = "sibling"
    a = pd.concat([a1, a2], ignore_index=True)
    a.to_csv(os.path.join(HERE, "combined_allele_table.tsv"), sep="\t", index=False)

    s1 = pd.read_csv(os.path.join(AL, "allele_segments.tsv"), sep="\t", dtype={"chrom": str}).assign(bench_set="iw_intpsi")
    s2 = pd.read_csv(os.path.join(AL, "allele_segments_allelic.tsv"), sep="\t", dtype={"chrom": str}).assign(bench_set="iw_allelic")
    s2["class"] = "sibling"
    s = pd.concat([s1, s2], ignore_index=True)
    s.to_csv(os.path.join(HERE, "combined_allele_segments.tsv"), sep="\t", index=False)

    r1 = pd.read_csv(os.path.join(AL, "rck_style_table.tsv"), sep="\t").assign(bench_set="iw_intpsi")
    r2 = pd.read_csv(os.path.join(AL, "rck_style_table_allelic.tsv"), sep="\t").assign(bench_set="iw_allelic")
    r2["class"] = "sibling"
    r = pd.concat([r1, r2], ignore_index=True)
    r.to_csv(os.path.join(HERE, "combined_rck_style_table.tsv"), sep="\t", index=False)
    return t, p, a, s, r



def fig2(main, out):
    """make_figures.fig2 with the legend on one row (four entries with three classes)."""
    import inspect
    src = inspect.getsource(mf.fig2).replace("ncol=3", "ncol=4")
    ns = dict(vars(mf))
    exec(src, ns)
    ns["fig2"](main, out)


# -- Fig. 3 with one column per class (the sibling class has S = 2 only) ------------------------

def fig3(main, out):
    metrics = [("blind_jseg_best_mean", "best-match J_seg (blind)"),
               ("blind_jadj_best_mean", "best-match J_adj (blind)")]
    off = {m: (i - (len(mf.METHODS) - 1) / 2) * 0.12 for i, m in enumerate(mf.METHODS)}
    Ss = {c: sorted(main[main["class"] == c].S.unique()) for c in CLASSES}
    fig, axes = plt.subplots(2, 3, figsize=(7.2, 4.6), sharey=True,
                             gridspec_kw={"width_ratios": [len(Ss[c]) + 0.4 for c in CLASSES]})
    rows = []
    for r, (col, ylab) in enumerate(metrics):
        for k, c in enumerate(CLASSES):
            ax = axes[r, k]
            for m in mf.METHODS:
                d = main[(main.m == m) & (main["class"] == c)]
                g = d.groupby("S")[col].agg(["mean", "std"]).reindex(Ss[c])
                xs = np.arange(len(Ss[c])) + off[m]
                ax.scatter(np.concatenate([[x] * int((d.S == s).sum()) for x, s in zip(xs, Ss[c])]),
                           np.concatenate([d[d.S == s][col].values for s in Ss[c]]),
                           s=10, color=mf.COLOR[m], alpha=0.45, linewidths=0, zorder=2)
                ax.plot(xs, g["mean"], color=mf.COLOR[m], marker="o", markersize=5,
                        markeredgecolor="white", markeredgewidth=0.8, zorder=3, label=mf.LABEL[m])
                for s, row in g.iterrows():
                    rows.append({"metric": col, "class": c, "S": s, "m": m, "n": int((d.S == s).sum()),
                                 "mean": row["mean"], "sd": row["std"]})
            ax.set_xticks(range(len(Ss[c])))
            ax.set_xticklabels([str(int(s)) for s in Ss[c]])
            ax.set_xlim(-0.6, len(Ss[c]) - 0.4)
            ax.set_ylim(-0.03, 1.05)
            if r == 0:
                ax.set_title("{} ({} cases)".format(mf.CLASS_LABEL[c], main[main["class"] == c].name.nunique()),
                             color=mf.INK)
            if r == 1:
                ax.set_xlabel("species per case (S)")
            if k == 0:
                ax.set_ylabel(ylab)
    axes[0, 0].legend(loc="lower left", frameon=False, ncol=1)
    mf.save(fig, out, "fig3_structure_vs_S", pd.DataFrame(rows))


# -- Fig. 6, as make_figures.fig6 with a class legend in (a) and class names once per class in (b) --

def fig6(main, out):
    rck = ["rck_mc"]
    d = main[main.m.isin(rck)]
    fig, axes = plt.subplots(1, 3, figsize=(7.4, 2.6), gridspec_kw={"width_ratios": [1, 1.25, 1]})

    ax = axes[0]
    ax.plot([0, 1], [0, 1], color=mf.INK2, linewidth=0.8, linestyle="--", zorder=1)
    for m in rck:
        s = d[d.m == m]
        for c in mf.CLASS_ORDER:
            q = s[s["class"] == c]
            ax.scatter(q.blind_jseg_best_mean, q.phased_jseg_best_mean, s=24, marker=mf.MARKER[c],
                       color=mf.COLOR[m], edgecolors="white", linewidths=0.8, alpha=0.9, zorder=3)
    ax.set_xlim(0.3, 1.03)
    ax.set_ylim(0.3, 1.03)
    ax.set_aspect("equal")
    ax.set_xlabel("J_seg, haplotype-blind")
    ax.set_ylabel("J_seg, phased")
    ax.set_title("a  allele-specific structure", loc="left", color=mf.INK)
    n11 = int(((d.blind_jseg_best_mean > 0.999) & (d.phased_jseg_best_mean > 0.999)).sum())
    ax.text(0.98, 0.03, "{} of {} runs at (1, 1)".format(n11, len(d)), transform=ax.transAxes,
            fontsize=6, color=mf.INK2, ha="right")
    ax.legend(handles=[Line2D([], [], marker=mf.MARKER[c], linestyle="", markersize=5, markerfacecolor="white",
                              markeredgecolor=mf.INK2, label=mf.CLASS_LABEL[c]) for c in mf.CLASS_ORDER],
              loc="upper left", frameon=False, handletextpad=0.2, fontsize=6)

    ax = axes[1]
    s = d[d.m == "rck_mc"].sort_values(["class", "S", "name"], ascending=[False, True, True])
    x = np.arange(len(s))
    for xi, a, w in zip(x, s.errors_absolute, s.errors_within):
        if a != w:
            ax.plot([xi, xi], [w, a], color=mf.MUTED, linewidth=1, zorder=1)
    ax.scatter(x, s.errors_absolute, s=20, facecolors="white", edgecolors=mf.COLOR["rck_mc"],
               linewidths=1.2, zorder=3, label="absolute A/B")
    ax.scatter(x, s.errors_within, s=20, color=mf.COLOR["rck_mc"], edgecolors="white",
               linewidths=0.6, zorder=4, label="within phase blocks")
    groups = list(zip(s["class"], s.S))
    centers, names = [], []
    for g in dict.fromkeys(groups):
        idx = [i for i, h in enumerate(groups) if h == g]
        centers.append(np.mean(idx))
        names.append("S={}".format(int(g[1])))
        if idx[0] > 0:
            ax.axvline(idx[0] - 0.5, color=mf.GRID, linewidth=1.2, zorder=0)
    ax.grid(axis="x", visible=False)
    ax.set_xticks(centers)
    ax.set_xticklabels(names, fontsize=6)
    ax.tick_params(axis="x", length=0)
    for c in dict.fromkeys(s["class"]):
        idx = [i for i, h in enumerate(s["class"]) if h == c]
        ax.text(np.mean(idx), -0.16, mf.CLASS_LABEL[c], transform=ax.get_xaxis_transform(), ha="center", va="top",
                fontsize=6.5, color=mf.INK)
    ax.set_ylabel("missed + false species")
    ax.set_title("b  phase-block identifiability", loc="left", color=mf.INK, pad=18)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 0.99), ncol=2, frameon=False,
              borderaxespad=0, handletextpad=0.2, columnspacing=0.8)

    ax = axes[2]
    rows = []
    width = 0.36
    for j, m in enumerate(rck):
        for k, c in enumerate(mf.CLASS_ORDER):
            q = d[(d.m == m) & (d["class"] == c)]
            for kind, col, face in [("blind", "blind_exact_cycles", mf.COLOR[m]),
                                    ("phased", "phased_exact_cycles", "white")]:
                rate = q[col].sum() / q.n_truth.sum()
                rows.append({"m": m, "class": c, "kind": kind, "rate": rate})
            b = [r for r in rows if r["m"] == m and r["class"] == c]
            xpos = k + (j - (len(rck) - 1) / 2) * width
            ax.bar(xpos, b[0]["rate"], width=width - 0.04, color=mf.COLOR[m], zorder=2)
            ax.scatter([xpos], [b[1]["rate"]], marker="_", s=120, color=mf.INK, linewidths=1.5,
                       zorder=3)
    ax.set_xticks(range(len(mf.CLASS_ORDER)))
    ax.set_xticklabels([mf.CLASS_LABEL[c] for c in mf.CLASS_ORDER])
    ax.set_ylim(0, 1.05)
    ax.set_ylabel("exact cycles / true species")
    ax.set_title("c  exact cycles", loc="left", color=mf.INK, pad=30)
    ax.legend(handles=[Line2D([], [], marker="s", linestyle="", color=mf.COLOR[m], label=mf.SHORT[m] +
                              " (blind)") for m in rck] +
              [Line2D([], [], marker="_", linestyle="", color=mf.INK, markersize=10, label="phased")],
              loc="lower left", bbox_to_anchor=(0, 0.99), ncol=2, frameon=False, fontsize=6,
              borderaxespad=0, handletextpad=0.2, columnspacing=0.6)
    fig.tight_layout(w_pad=1.2)
    fig.text(0.01, -0.02, "CoRAL and CE report no alleles, so they have no phased scores.", fontsize=6.5,
             color=mf.INK2, va="top")
    mf.save(fig, out, "fig6_allele_specific",
         d[["name", "class", "S", "m", "blind_jseg_best_mean", "phased_jseg_best_mean",
            "n_missed", "n_false", "n_missed_absolute", "n_false_absolute", "blind_exact_cycles",
            "phased_exact_cycles", "n_truth"]])
    pd.DataFrame(rows).to_csv(os.path.join(out, "figure_data", "fig6c_exact_rates.tsv"),
                              sep="\t", index=False)



# -- Fig. 8, as make_allele_figures.fig8 with the case count from the data --

def fig8(a, out):
    phasing_of = {"rck_mc": "native", "coral": "baf_major", "ce": "baf_major"}
    t = maf.within_block_tables(a, phasing_of)
    t = t[(t["end"] - t["start"] + 1) >= 1000]  # 1-bp slivers from breakpoint offsets add nothing
    mixed_base = (t.t_A > 0) & (t.t_B > 0)
    long = pd.concat([t.assign(allele=h, truth=t["t_" + h], inferred=t["c_" + h], mixed_base=mixed_base) for h in ("A", "B")])
    long = long[(long.truth > 0) | (long.inferred > 0)]
    fig, axes = plt.subplots(1, len(mf.METHODS), figsize=(5.8, 2.5), sharex=True, sharey=True, gridspec_kw={"wspace": 0.08})
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
        ax.set_title("{s} ({p})\n\n".format(s=mf.SHORT[m], p=maf.PHASE_LABEL[phasing_of[m]]), fontsize=8)
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
    fig.text(0.01, -0.27, "One point per ecDNA segment (>= 1 kb) and allele, all {n} cases; within phase blocks. "
             "CoRAL / CE phased post-hoc from the observed CN_A / CN_B.".format(n=a["name"].nunique()), fontsize=6.2, color=mf.INK2,
             ha="left", va="top")
    mf.save(fig, out, "fig8_allele_cn", long[["name", "m", "phasing", "group", "chrom", "start", "end", "allele", "truth",
                                               "inferred", "mixed_base"]])



# -- allele figures: combined inputs, phase blocks read from the local subgraph folders ----------

def allele_load(a_all):
    def load(baseline_variant):
        a = a_all[a_all["status"] == "ok"].copy()
        mixed = a[(a.method_variant == "rck_mc") & (a.phasing == "native")].set_index("name")["mixed_bp_frac"]
        a["group"] = np.where(a["name"].map(mixed) > 0, "mixed", "single")
        keep = {"rck_mc": "rck_mc", "coral:" + baseline_variant: "coral", "ce:" + baseline_variant: "ce"}
        a = a[a.method_variant.isin(keep)].copy()
        a["m"] = a.method_variant.map(keep)
        for col in ("phased_jseg_best_mean_block", "phased_jadj_best_mean_block", "phased_exact_cycles_block"):
            a[col] = a[col].fillna(0.0)
        return as_classes(a)
    return load


def within_block_local(segments):
    from rck.ecdna import allele_metrics as am, evaluate as ev

    def tables(a, phasing_of):
        frames, n_swaps = [], {}
        for r in a.drop_duplicates(["name", "m"]).itertuples():
            p = phasing_of[r.m]
            row = a[(a.name == r.name) & (a.m == r.m) & (a.phasing == p)].iloc[0]
            t = segments[(segments.name == r.name) & (segments.method_variant == row.method_variant)
                         & (segments.phasing == p)].drop(columns=["name", "class", "S", "method_variant",
                                                                  "phasing", "bench_set"])
            frag = os.path.join(SUBGRAPHS, SET_DIR[row.bench_set], r.name, "union", "rck_input", "rck.fragments.tsv")
            t, k = am.within_block(t.reset_index(drop=True), ev.read_phase_blocks(frag))
            n_swaps[(r.name, r.m)] = (k, row.n_blocks_swapped)
            frames.append(t.assign(name=r.name, m=r.m, phasing=p, group=row.group))
        bad = {key: v for key, v in n_swaps.items() if v[0] != v[1]}
        assert not bad, "block swaps differ from the collected table: {}".format(bad)
        return pd.concat(frames, ignore_index=True)
    return tables


def table2_truth(t):
    n_truth = t.drop_duplicates("name").set_index("name")["n_truth"].to_dict()

    def table2(a, out):
        # make_allele_figures.table2 counts truth species from each case's truth folder (cluster);
        # the collected bench table already holds that count (n_truth).
        rows, md = [], ["# Table 2: allele-specific ecDNA reconstruction (combined benchmark, within phase blocks)", ""]
        for g, title in maf.GROUPS:
            sub = a[a.group == g]
            md += ["**{t}** (n = {n} cases)".format(t=title, n=sub["name"].nunique()), "",
                   "| method | phasing | " + " | ".join(lbl for _, lbl, _ in maf.TABLE2) + " | phased exact-cycle rate |",
                   "|---|---|" + "---|" * (len(maf.TABLE2) + 1)]
            for m, p in maf.columns():
                d = sub[(sub.m == m) & (sub.phasing == p)]
                row = {"group": g, "method": m, "phasing": p, "n": len(d)}
                cells = []
                for col, lbl, fmt in maf.TABLE2:
                    row[col + "_mean"], row[col + "_sd"] = d[col].mean(), d[col].std()
                    cells.append((fmt + " ± " + fmt).format(d[col].mean(), d[col].std()))
                rate = d.phased_exact_cycles_block.sum() / sum(n_truth[n] for n in d["name"])
                row["phased_exact_rate"] = rate
                rows.append(row)
                md.append("| {m} | {p} | {c} | {r:.2f} |".format(m=mf.SHORT[m], p=maf.PHASE_LABEL[p],
                                                                 c=" | ".join(cells), r=rate))
            md.append("")
        pd.DataFrame(rows).to_csv(os.path.join(out, "table2.tsv"), sep="\t", index=False)
        with open(os.path.join(out, "table2.md"), "w") as handle:
            handle.write("\n".join(md) + "\n")
    return table2


def rck_style_load(r_all):
    def load(baseline_variant):
        t = r_all[r_all["status"] == "ok"].copy()
        t["group"] = np.where(t["mixed_bp_frac"] > 0, "mixed", "single")
        t["is_homologous"] = t["is_is_hom"] + t["is_is_recip"]
        truth = t[t.method == "truth"].copy()
        keep = {"rck_mc": "rck_mc", "coral:" + baseline_variant: "coral", "ce:" + baseline_variant: "ce"}
        t = t[t.method_variant.isin(keep)].copy()
        t["m"] = t.method_variant.map(keep)
        return as_classes(t), as_classes(truth)
    return load


def fig9(t, truth, out):
    """make_rck_style_figures.fig9 with the legend on one row and "allele" for "homolog" in the footnote."""
    import inspect
    src = inspect.getsource(mrs.fig9).replace("ncol=4", "ncol=5")
    ns = dict(vars(mrs))
    ns["FOOT9"] = mrs.FOOT9.replace("both homologs", "both alleles").replace("opposite\nhomologs", "opposite\nalleles")
    exec(src, ns)
    ns["fig9"](t, truth, out)


def main():
    os.makedirs(os.path.join(OUT, "figure_data"), exist_ok=True)
    t_all, props, a_all, segs, r_all = combine_tables()

    # Figs. 2-6 and Table 1
    t_path = os.path.join(HERE, "combined_table.tsv")
    t, main_df, n_empty = mf.load(t_path, "walks_kept")
    t, main_df = as_classes(t), as_classes(main_df)
    fig2(main_df, OUT)
    fig3(main_df, OUT)
    mf.fig4(main_df, OUT)
    mf.fig5(main_df, props, "walks_kept", OUT)
    fig6(main_df, OUT)
    mf.table1(t, OUT)

    # Figs. 7-8 and Table 2
    a = allele_load(a_all)("walks_kept")
    maf.within_block_tables = within_block_local(segs)
    maf.fig7(a, OUT)
    fig8(a, OUT)
    table2_truth(t_all)(a, OUT)

    # Fig. 9 and Table 3
    r, truth = rck_style_load(r_all)("walks_kept")
    fig9(r, truth, OUT)
    mrs.table3(r, truth, OUT)
    # the imported writers title their markdown for bench_IW_intpsi; retitle for this benchmark
    for name, fixes in {"table1.md": [("bench_IW_intpsi", "combined 36-case benchmark")],
                        "table3.md": [("bench_IW_intpsi", "combined 36-case benchmark"),
                                      ("n = 6 mixed cases: descriptive only.",
                                       "n = 18 mixed cases (12 sibling, 6 nested): descriptive only.")]}.items():
        path = os.path.join(OUT, name)
        with open(path) as handle:
            text = handle.read()
        for a, b in fixes:
            assert a in text, (name, a)
            text = text.replace(a, b)
        with open(path, "w") as handle:
            handle.write(text)
    print("combined benchmark: {} cases ({}); {} runs with no cycles counted as 0 -> {}".format(
        t_all.name.nunique(), t_all.drop_duplicates("name")["class"].value_counts().to_dict(), n_empty, OUT))


if __name__ == "__main__":
    main()

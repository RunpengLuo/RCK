#!/usr/bin/env python3
"""bench_IW_allele item 1 (allelic siblings): species-coverage figure per scenario, rck_mc vs CoRAL vs CE.

One figure per scenario (figures/item1_species_<scenario>.{png,pdf}; numbers in
figures/figure_data/item1_species_<scenario>.tsv). Rows: the 4 cases, each with its two truth species (outer,
sibling) as bars laid out in the truth walk's step order (length normalised; true length in the row label).
Columns: rck_mc (haplotype-aware: a base counts as covered only by an inferred species on the same allele,
after the per-phase-block A/B swap the scorer chooses), CoRAL and CE (haplotype-blind, as they have no
alleles; linear walks kept). Each base is coloured by the inferred species covering it, with multiplicity
(make_species_matrix_figure.coverage_pieces); a base covered by several is split in height.

Colours follow the role of the inferred species in that panel, the same in every panel: the species paired
with the outer species, the one paired with the sibling, then any other strict species. Pairing is the
scorer's: within phase blocks with haplotypes for rck_mc (evaluate.score_decomposition), haplotype-blind for
CoRAL / CE. Strict species only (bench_metrics / aa_eval strict selection). Beside each bar: the paired
species' share and how many other species cover the bar, with J_seg / J_adj of the pair (the scored
aa_species_matrix, phased within blocks for rck_mc, blind for CoRAL / CE); for a truth species with no pair,
"missed" and the species that cover its bases. CoRAL / CE on sib_opp_mirror_s17006 ran on the zero-read-patched graph
(allelic_bench.py, INPUT PATCH). rck_bulk is not shown (graph builder only).

    python make_item1_figures.py     # local, seconds, no Gurobi
"""
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import make_species_matrix_figure as smf  # noqa: E402  (imports allele_bench)
import allele_bench as alb  # noqa: E402

alb.use_set("iw_allelic")
import item1_results as ir  # noqa: E402

import matplotlib  # noqa: E402

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import pandas as pd  # noqa: E402
from matplotlib.patches import Patch, Rectangle  # noqa: E402

from rck.ecdna import evaluate as ev  # noqa: E402

bi = alb.bi
FIG = os.path.join(HERE, "figures")
PANELS = [("rck_mc", "", "native", True, "rck_mc (allele-specific)"),
          ("coral", "walks_kept", "baf_major", False, "CoRAL (haplotype-blind)"),
          ("ce", "walks_kept", "baf_major", False, "CE (haplotype-blind)")]
ROLE_COLOURS = {"outer": smf.PALETTE[0], "sibling": smf.PALETTE[1]}
EXTRA_COLOURS = smf.PALETTE[2:]


def panel(case, method, variant, phasing, with_hap):
    """truth list, strict inferred cycles (relabelled within blocks for rck_mc), truth -> inferred pairing,
    and the scored matrix rows for the pairs."""
    truth = ev.load_truth(case["truth"])
    weights = None
    if method == "rck_mc":
        from rck.ecdna import bench_metrics as bm
        weights = bm.multiclone_bulk_weights(bi.run_dir(case, method), bm.read_proportions(
            os.path.join(case["rck_input"], "clone_proportions.tsv")))
    cycles = ir._answer_cycles(alb.answer_files(case, method, variant)[0], weights)
    if with_hap:
        blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
        cycles, pairing = ir._within_block(truth, cycles, blocks)
        sub, hap = "aa_eval_block", "phased"
    else:
        pairing = ev._match(truth, cycles, False)
        sub, hap = "aa_eval", "blind"
    m = pd.read_csv(os.path.join(alb.task_dir(case, method, variant, phasing), sub, "aa_species_matrix.tsv"), sep="\t")
    m = m[m["haplotype"] == hap].set_index(["truth", "inferred"])
    return truth, cycles, pairing, m


def draw_scenario(scenario):
    cases = [c for c in bi.cases() if c["scenario"] == scenario]
    fig, axes = plt.subplots(1, len(PANELS), figsize=(15, 1.15 + 1.25 * len(cases)), sharey=True)
    rows_out, used_extra = [], 0
    for ax, (method, variant, phasing, with_hap, title) in zip(axes, PANELS):
        ylabels = []
        for ci, case in enumerate(cases):
            truth, cycles, pairing, m = panel(case, method, variant, phasing, with_hap)
            inferred = {c.name: c for c in cycles}
            roles = {pairing[t.name]: role for t, role in zip(truth, ("outer", "sibling")) if t.name in pairing}
            extras = [c.name for c in cycles if c.name not in roles]
            colours = {n: ROLE_COLOURS[r] for n, r in roles.items()}
            colours.update({n: EXTRA_COLOURS[k % len(EXTRA_COLOURS)] for k, n in enumerate(extras)})
            used_extra = max(used_extra, len(extras))
            names = list(roles) + extras
            for ti, (t, role) in enumerate(zip(truth, ("outer", "sibling"))):
                y = ci * 2.6 + ti
                pieces = smf.coverage_pieces(t, inferred, names, with_hap)
                total = float(sum(ln for _x, ln, _c in pieces))
                h = 0.72
                for x, ln, cov in pieces:
                    if not cov:
                        ax.add_patch(Rectangle((x / total, y - h / 2), ln / total, h, facecolor=smf.UNCOVERED, edgecolor="none"))
                        continue
                    sub = h / len(cov)
                    for k, n in enumerate(cov):
                        ax.add_patch(Rectangle((x / total, y - h / 2 + k * sub), ln / total, sub,
                                               facecolor=colours[n], edgecolor="white", linewidth=0.4))
                sh = smf.shares(pieces, names)
                paired = pairing.get(t.name)
                js = m.loc[(t.name, paired)] if paired and (t.name, paired) in m.index else None
                score = "J %.2f / %.2f" % (js["j_seg"], js["j_adj"]) if js is not None else "missed"
                share = dict(sh)
                others = [n for n, f in sh if n not in (paired, "uncovered") and f >= 0.005]
                if paired:
                    text = "%s %d%%" % (paired, round(100 * share.get(paired, 0.0)))
                    text += " +%d other%s" % (len(others), "s" if len(others) > 1 else "") if others else ""
                else:
                    text, score = "missed", "covered by " + (", ".join(others) or "none")
                ax.text(1.03, y, "%s\n%s" % (text, score), va="center", ha="left", fontsize=7, color=smf.INK,
                        linespacing=1.3)
                ylabels.append((y, "%s %s\n%s, %d kb" % (case["name"].rsplit("_", 1)[1], case["stratum"], role,
                                                          round(total / 1e3))))
                rows_out.append({"scenario": scenario, "case": case["name"], "stratum": case["stratum"],
                                 "method": method, "variant": variant or "", "haplotype_aware": with_hap,
                                 "truth": t.name, "role": role, "truth_len_bp": int(total), "paired": paired,
                                 "j_seg": None if js is None else js["j_seg"], "j_adj": None if js is None else js["j_adj"],
                                 "shares": "; ".join("%s=%.4f" % (n, f) for n, f in sh)})
        ax.set_xlim(0, 1.42)
        ax.set_title(title, fontsize=10, loc="left", color=smf.INK)
        ax.set_xticks([0, 0.5, 1.0])
        ax.set_xticklabels(["0", "50%", "100%"], fontsize=7, color=smf.INK_MUTED)
        for side in ("top", "right", "left"):
            ax.spines[side].set_visible(False)
        ax.spines["bottom"].set_bounds(0, 1)
        ax.tick_params(axis="y", length=0)
    ys, labels = zip(*ylabels)
    axes[0].set_yticks(ys)
    axes[0].set_yticklabels(labels, fontsize=7)
    axes[0].set_ylim(max(ys) + 0.6, -0.6)
    handles = [Patch(facecolor=ROLE_COLOURS["outer"], label="species paired with the outer"),
               Patch(facecolor=ROLE_COLOURS["sibling"], label="species paired with the sibling")]
    handles += [Patch(facecolor=EXTRA_COLOURS[k], label="other species %d" % (k + 1)) for k in range(used_extra)]
    handles += [Patch(facecolor=smf.UNCOVERED, label="not covered")]
    fig.legend(handles=handles, loc="lower center", ncol=len(handles), frameon=False, fontsize=8)
    fig.suptitle("%s: truth species covered by each method's species (position along the truth walk)" % scenario,
                 fontsize=11, x=0.01, ha="left")
    fig.subplots_adjust(left=0.1, right=0.995, top=0.9, bottom=0.11, wspace=0.08)
    path = os.path.join(FIG, "item1_species_" + scenario)
    for ext in ("png", "pdf"):
        fig.savefig("%s.%s" % (path, ext), dpi=200)
    plt.close(fig)
    os.makedirs(os.path.join(FIG, "figure_data"), exist_ok=True)
    pd.DataFrame(rows_out).to_csv(os.path.join(FIG, "figure_data", "item1_species_%s.tsv" % scenario), sep="\t", index=False)
    print("wrote " + path + ".{png,pdf}")


def main():
    import logging
    logging.disable(logging.INFO)
    os.makedirs(FIG, exist_ok=True)
    for s in ir.SCENARIOS:
        draw_scenario(s)


if __name__ == "__main__":
    main()

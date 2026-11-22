#!/usr/bin/env python3
"""Within-phase-block recovery matrices and two-panel J_seg / J_adj species matrices (Yubi, 2026-11-21).

Why: the run's own recovery_matrix.tsv (rck-ecdna's evaluation, drawn in rck_mc_recovery_matrix.png) pairs truth
and inferred species, and scores haplotypes, in the absolute frame. A phase block's A/B orientation cannot be
identified from the data, so the benchmark tables score within phase blocks (aa_eval.score with --phase-blocks:
blind pairing of the strict species, one A/B swap per block, evaluate.block_swaps; relabel; pair again with
haplotypes). Siblings that differ only by allele are then mispaired by the absolute frame. This script redraws
everything in the table's frame. Run folders are only read.

Per rck_mc case (12 iw_allelic + 24 bench_IW_intpsi), into topology_plots/species_matrix/<set>/<case>/:
  rck_mc_recovery_matrix_block.tsv   evaluate.RECOVERY_COLUMNS, within phase blocks: for every rank of the run's
                                     decompositions.tsv (rck_mc has one), blind pairing -> block_swaps -> relabel
                                     -> evaluate.score_decomposition (pairing by haplotype Jaccard, every score on
                                     that pair, CN error included); best_* over the ranks. Plus table_composite (the
                                     table scorer's per-species composite, score/aa_eval/aa_species_metrics.tsv
                                     strict) and table_inferred.
  rck_mc_recovery_matrix_block.png   rck.ecdna.plot.matrix_page of that TSV, retitled "within phase blocks", with
                                     the block footnote.
  rck_mc_species_matrix_block.{png,pdf}  two panels, J_seg | J_adj (make_species_matrix_figure --layout pair),
                                     aa_eval_block/ phased cells, scorer pairing within phase blocks.
  rck_mc_species_matrix_blind.{png,pdf}  the same, haplotype-blind (blind cells and blind pairing).
Per CoRAL compare case (8): coral_species_matrix_blind.{png,pdf} (haplotype-blind; linear walks kept).

Tables, in topology_plots/species_matrix/:
  jmatrix_summary.tsv     one row per (set, case, frame, truth species): paired inferred species, paired J_seg /
                          J_adj, row maxima and whether they differ from the paired cell, exact_cycle of the paired
                          cell, the (1, 1) flag (both >= 0.995, i.e. 1.00 as printed), scenario / class.
  recovery_check.tsv      one row per case: mean of the new composite vs the table's mean_composite, and the
                          per-species composite / pairing vs the table scorer's.

    python make_block_figures.py run [--jobs 8] [--only CASE ...]
    python make_block_figures.py reasons    # jmatrix_not_1_1.tsv: each non-(1, 1) RCK MC cell and why
    python make_block_figures.py compare    # topology_plots/compare/: old files renamed, new ones copied in

`compare` renames, once, rck_mc_recovery_matrix.png -> rck_mc_recovery_matrix_absolute.png (it is the absolute
frame) and the one-panel species matrices <method>_species_matrix_<frame>.{png,pdf} ->
<method>_species_matrix_<frame>_onepanel.{png,pdf} (those were already scorer-paired, so not "absolute"), then
copies the new files in. No Gurobi.
"""
import argparse
import os
import shutil
import sys
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allele_bench as ab  # noqa: E402
import make_species_matrix_figure as msf  # noqa: E402

import pandas as pd  # noqa: E402

OUT = os.path.join(HERE, "topology_plots", "species_matrix")
COMPARE = os.path.join(HERE, "topology_plots", "compare")
TABLES = {"iw_allelic": os.path.join(HERE, "bench_iw_allelic_table.tsv"),
          "iw_intpsi": os.path.join(HERE, "..", "bench_IW", "bench_iw_intpsi_table.tsv")}
TOL = 10_000  # evaluate / aa_eval SPECIES_TOL
ONE = 0.995  # "1.00" as printed


def compare_cases():
    import make_compare_plots as mcp
    return [(s, n, mcp.folder_of(k, n)) for k, (s, n, _w) in enumerate(mcp.CASES)]


def case_of(name):
    return next(c for c in ab.bi.cases() if c["name"] == name)


# ------------------------------------------------------------------------------------------------ #
# task 1: recovery matrix within phase blocks
# ------------------------------------------------------------------------------------------------ #

def recovery_block(case):
    from rck.ecdna import evaluate as ev
    from dataclasses import asdict

    run = ab.bi.run_dir(case, "rck_mc")
    truth = ev.load_truth(case["truth"])
    blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
    inferred = ev.read_inferred(os.path.join(run, "decompositions.tsv"))
    tb = {t.name: t for t in truth}
    rows = []
    for rank in sorted(inferred):
        cycles = inferred[rank]["cycles"]
        by = {c.name: c for c in cycles}
        swaps = ev.block_swaps([(tb[a], by[b]) for a, b in ev._match(truth, cycles, False).items()], blocks)
        relabelled = [ev.relabel(c, blocks, swaps) for c in cycles]
        for s in ev.score_decomposition(truth, relabelled, TOL):
            r = asdict(s)
            r["decomposition_rank"] = rank
            r["n_blocks_swapped"] = int(sum(swaps.values()))
            rows.append(r)
    frame = pd.DataFrame(rows)
    info = [{"species": t.name, "weight": t.weight, "length_bp": sum(x.end - x.start + 1 for x in t.steps),
             "n_junctions": len(ev._unique(ev.junctions(t.steps, circular=t.circular), TOL)),
             "fraction_in_subgraph": 1.0, "in_scope": True} for t in truth]
    matrix = ev.recovery_matrix(frame, info, min(inferred))
    top = frame[frame["decomposition_rank"] == min(inferred)].set_index("truth")
    matrix["n_blocks_swapped"] = int(top["n_blocks_swapped"].iloc[0]) if len(top) else 0
    # the table scorer's per-species rows (aa_eval.score, within phase blocks)
    sm = pd.read_csv(os.path.join(ab.bi.score_folder(case, "rck_mc", None), "aa_eval", "aa_species_metrics.tsv"), sep="\t")
    sm = sm[sm["mode"] == "strict"].set_index("truth")
    matrix["table_inferred"] = [sm.loc[t, "inferred"] if t in sm.index else None for t in matrix["truth"]]
    matrix["table_composite"] = [float(sm.loc[t, "composite"]) if t in sm.index else None for t in matrix["truth"]]
    return matrix


def rank_swaps(case):
    """(blocks, swaps, rank) of the run's top decomposition: the per-block A/B swaps recovery_block uses (blind pairing of
    the decomposition's species, evaluate.block_swaps), for relabelling the run's own files (make_compare_plots pages)."""
    from rck.ecdna import evaluate as ev

    run = ab.bi.run_dir(case, "rck_mc")
    truth = ev.load_truth(case["truth"])
    blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
    inferred = ev.read_inferred(os.path.join(run, "decompositions.tsv"))
    rank = min(inferred)
    cycles = inferred[rank]["cycles"]
    tb, by = {t.name: t for t in truth}, {c.name: c for c in cycles}
    return blocks, ev.block_swaps([(tb[a], by[b]) for a, b in ev._match(truth, cycles, False).items()], blocks), rank


def draw_recovery(matrix, case, table_mean, path):
    from rck.ecdna import plot as rp

    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    fig = rp.matrix_page(plt, matrix, case["name"])
    ax = fig.axes[0]
    mean = pd.to_numeric(matrix["composite"], errors="coerce").mean()
    ax.set_title("{n}: RCK MC recovery of the truth species, within phase blocks (rank-1 decomposition; best = best over "
                 "all alternatives)\nmean composite {m:.3f} (benchmark table: {t:.3f}); species paired as the table's "
                 "scorer pairs them".format(n=case["name"], m=mean, t=table_mean), fontsize=10, color=rp.INK, loc="left")
    fig.text(0.01, 0.01, msf.BLOCK_NOTE, fontsize=8, color=rp.INK_MUTED, va="bottom")
    fig.set_size_inches(fig.get_figwidth(), fig.get_figheight() + 0.35)
    fig.tight_layout(rect=(0, 0.06, 1, 1))
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ------------------------------------------------------------------------------------------------ #
# task 2: two-panel J matrices
# ------------------------------------------------------------------------------------------------ #

def jmatrix(set_name, case, method, frame, outdir):
    """Draw <method>_species_matrix_<frame> (two panels) and return the summary rows."""
    msf.METHOD, msf.MATCH = method, "score"
    msf.VARIANT, msf.PHASING = (None, "native") if method == "rck_mc" else ("walks_kept", "baf_major")
    block = frame == "block"
    hap = "phased" if block else "blind"
    matrix, truth, inferred, lay = msf.load(case["name"], block, hap)
    rl, cl = msf.labels(matrix, truth, inferred, lay)
    rl = {k: v.replace(" copies", " bulk copies") for k, v in rl.items()}
    cl = {k: v.replace(" copies", " bulk copies") for k, v in cl.items()}
    path = os.path.join(outdir, "%s_species_matrix_%s" % (method, frame))
    maxima = msf.figure_matrix_pair(matrix, rl, cl, lay, msf.pair_title(case["name"], block, not block), path,
                                    msf.BLOCK_NOTE if block else None)
    rows = []
    for r in lay["rows"]:
        paired = lay["matches"].get(r)
        cell = matrix[(matrix["truth"] == r) & (matrix["inferred"] == paired)] if paired else matrix.iloc[0:0]
        js, ja = (float(cell["j_seg"].iloc[0]), float(cell["j_adj"].iloc[0])) if len(cell) else (float("nan"),) * 2
        ms, ma = maxima[r]["j_seg"], maxima[r]["j_adj"]
        rows.append({
            "set": set_name, "case": case["name"], "scenario": case.get("scenario"), "class": case.get("class"),
            "S": case.get("S"), "stratum": case.get("stratum", ""), "method": msf.LABEL[method],
            "frame": "within phase blocks" if block else "haplotype-blind", "truth": r,
            "truth_bulk_copies": float(matrix.loc[matrix["truth"] == r, "truth_weight"].iloc[0]),
            "paired": paired or "", "paired_j_seg": js, "paired_j_adj": ja,
            "rowmax_j_seg_species": ms[0] or "", "rowmax_j_seg": ms[1], "rowmax_j_seg_differs": ms[3],
            "rowmax_j_adj_species": ma[0] or "", "rowmax_j_adj": ma[1], "rowmax_j_adj_differs": ma[3],
            "exact_cycle": int(cell["exact_cycle"].iloc[0]) if len(cell) and pd.notna(cell["exact_cycle"].iloc[0]) else 0,
            "is_1_1": bool(len(cell) and js >= ONE and pd.notna(ja) and ja >= ONE),
            "n_inferred_drawn": len(lay["species"]),
            "n_unmatched": len(lay["unmatched"])})
    return rows


# ------------------------------------------------------------------------------------------------ #
# driver
# ------------------------------------------------------------------------------------------------ #

def init(set_name):
    import logging
    logging.disable(logging.INFO)
    ab.use_set(set_name)


def task(args):
    set_name, name, coral = args
    case = case_of(name)
    outdir = os.path.join(OUT, set_name, name)
    os.makedirs(outdir, exist_ok=True)
    summary, check = [], None
    table = pd.read_csv(TABLES[set_name], sep="\t")
    trow = table[(table["name"] == name) & (table["method_variant"] == "rck_mc")].iloc[0]
    matrix = recovery_block(case)
    matrix.to_csv(os.path.join(outdir, "rck_mc_recovery_matrix_block.tsv"), sep="\t", index=False)
    draw_recovery(matrix, case, float(trow["mean_composite"]), os.path.join(outdir, "rck_mc_recovery_matrix_block.png"))
    mean = float(matrix["composite"].mean())
    diff_c = (matrix["composite"] - matrix["table_composite"]).abs().max()
    same_pair = bool((matrix["inferred"].fillna("") == matrix["table_inferred"].fillna("")).all())
    check = {"set": set_name, "case": name, "n_truth": len(matrix), "table_n_truth": int(trow["n_truth"]),
             "mean_composite_block": mean, "table_mean_composite": float(trow["mean_composite"]),
             "abs_diff_mean": abs(mean - float(trow["mean_composite"])), "max_abs_diff_species": float(diff_c),
             "pairing_same_as_table": same_pair, "n_blocks_swapped": int(matrix["n_blocks_swapped"].iloc[0]),
             "agrees": bool(abs(mean - float(trow["mean_composite"])) < 1e-6 and diff_c < 1e-6 and same_pair)}
    for frame in ("block", "blind"):
        summary += jmatrix(set_name, case, "rck_mc", frame, outdir)
    if coral:
        summary += jmatrix(set_name, case, "coral", "blind", outdir)
    print("%s %s: done (%s)" % (set_name, name, "agrees" if check["agrees"] else "DISAGREES"), flush=True)
    return summary, check


def run(jobs, only):
    coral_names = {(s, n) for s, n, _f in compare_cases()}
    summary, checks = [], []
    for set_name in ("iw_allelic", "iw_intpsi"):
        init(set_name)
        names = [c["name"] for c in ab.bi.cases() if not only or c["name"] in only]
        tasks = [(set_name, n, (set_name, n) in coral_names) for n in names]
        with Pool(jobs, initializer=init, initargs=(set_name,)) as pool:
            for s, c in pool.imap(task, tasks):
                summary += s
                checks.append(c)
    os.makedirs(OUT, exist_ok=True)
    if not only:
        pd.DataFrame(summary).to_csv(os.path.join(OUT, "jmatrix_summary.tsv"), sep="\t", index=False, float_format="%.6g")
        pd.DataFrame(checks).to_csv(os.path.join(OUT, "recovery_check.tsv"), sep="\t", index=False, float_format="%.6g")
    bad = [c for c in checks if not c["agrees"]]
    print("recovery matrix vs table: %d of %d cases agree%s" % (len(checks) - len(bad), len(checks),
          "; DISAGREE: " + ", ".join(c["case"] for c in bad) if bad else ""))


def compare():
    for set_name, name, folder in compare_cases():
        src = os.path.join(OUT, set_name, name)
        renames = [("rck_mc_recovery_matrix.png", "rck_mc_recovery_matrix_absolute.png")]
        for stem in ("rck_mc_species_matrix_block", "rck_mc_species_matrix_blind", "coral_species_matrix_blind"):
            renames += [(stem + ext, stem + "_onepanel" + ext) for ext in (".png", ".pdf")]
        for old, new in renames:
            if os.path.exists(os.path.join(folder, old)) and not os.path.exists(os.path.join(folder, new)):
                os.rename(os.path.join(folder, old), os.path.join(folder, new))
        for f in ("rck_mc_recovery_matrix_block.png", "rck_mc_recovery_matrix_block.tsv",
                  "rck_mc_species_matrix_block.png", "rck_mc_species_matrix_block.pdf",
                  "rck_mc_species_matrix_blind.png", "rck_mc_species_matrix_blind.pdf",
                  "coral_species_matrix_blind.png", "coral_species_matrix_blind.pdf"):
            shutil.copy(os.path.join(src, f), os.path.join(folder, f))
        print("%s: updated" % folder)


# species-level structure of the non-balanced cases, from the saved results (decompositions.tsv, PLAN.md)
STRUCTURAL = {
    ("sib_opp_shared_s18007", "ecdna1"): "outer returned as two cycles meeting at a reused extremity (ec2 1410 kb psi "
                                         "(40, 20) and ec3 230 kb psi (40, 19); truth psi (40, 20))",
    ("simple_S2_s12004", "ecdna1"): "species returned as two cycles (ec1 J_seg 0.65, ec2 J_seg 0.35, both psi (17, 8))",
    ("simple_S4_s13002", "ecdna1"): "low-copy species (5.2 copies) spread over 4 inferred cycles (ec3, ec5, ec6, ec7); "
                                    "paired ec7 has J_seg 0.80, 3 of 9 adjacencies",
    ("hard_S2_s15003", "ecdna1"): "partial miss: ec2 covers 0.59 of the bases, 2 of 6 truth adjacencies plus 3 false "
                                  "ones (to chr16:13.4-13.6 Mb)",
    ("hard_S4_s16003", "ecdna3"): "known weak case (item 8): ec1 J_seg 0.40, one of 3 adjacencies, one false",
    ("hard_S4_s16003", "ecdna1"): "known weak case (item 8): spread over ec2 / ec4 / ec5 (J_seg 0.39 / 0.55 / 0.45) "
                                  "and a low_weight exact copy ec6",
}


def _reciprocal(js):
    """True when every breakend of the junctions ``js`` has a partner breakend 1 bp away on the opposite side
    among the others (a balanced inversion pair or reciprocal chain)."""
    ends = [(a[0], a[1], a[2]) for j in js for a in j[:2]]
    if len(js) < 2:
        return False
    for chrom, pos, strand in ends:
        if not any(c == chrom and abs(p - pos) == 1 and s != strand for c, p, s in ends):
            return False
    return True


def reasons():
    """jmatrix_not_1_1.tsv: every RCK MC within-block paired cell that is not (1, 1), with the truth novel adjacencies
    the paired walk misses / adds (unphased, 100 bp) and a reason."""
    import json
    from rck.ecdna import aa_eval as ae, aa_format as aaf, evaluate as ev

    s = pd.read_csv(os.path.join(OUT, "jmatrix_summary.tsv"), sep="\t")
    b = s[(s["method"] == "RCK MC") & (s["frame"] == "within phase blocks") & (~s["is_1_1"])]
    lw_rows = []
    fmt = lambda js: " ".join("%s:%d%s-%s:%d%s" % (a[0], a[1], a[2], c[0], c[1], c[2]) for a, c, *_ in js)  # noqa: E731
    for set_name in ("iw_allelic", "iw_intpsi"):
        init(set_name)
        for r in b[b["set"] == set_name].to_dict("records"):
            case = case_of(r["case"])
            folder = ab.task_dir(case, "rck_mc", None, "native")
            truth = {t.name: t for t in ev.load_truth(case["truth"])}
            with open(os.path.join(folder, "done.json")) as fh:
                records = aaf.read_cycles(json.load(fh)["cycles"])
            cycles = {c.name: c for c in ae.species_cycles(records)}
            t, c = truth[r["truth"]], cycles[r["paired"]]
            tj = ev._unique(ev.junctions(t.steps, circular=t.circular), 100)
            cj = ev._unique(ev.junctions(c.steps, circular=c.circular), 100)
            miss = [j for j in tj if not any(ev._junction_match(j, k, 100, False) for k in cj)]
            extra = [j for j in cj if not any(ev._junction_match(j, k, 100, False) for k in tj)]
            m = pd.read_csv(os.path.join(folder, "aa_eval_block", "aa_species_matrix.tsv"), sep="\t")
            m = m[(m["haplotype"] == "phased") & (m["truth"] == r["truth"]) & (m["inferred"] != r["paired"])]
            full = m[m["j_seg"] >= ONE].sort_values(["low_weight", "j_adj"], ascending=[True, False])  # same bases
            if (r["case"], r["truth"]) in STRUCTURAL:
                kind, why = "structural", STRUCTURAL[(r["case"], r["truth"])]
            elif r["paired_j_seg"] >= ONE and not extra and _reciprocal(miss):
                kind = "balanced rearrangement walked through the reference"
                if len(full) and not full["low_weight"].iloc[0]:
                    why = "rearranged form returned as a separate strict species (%s, J_adj %.2f)" % (
                        full["inferred"].iloc[0], full["j_adj"].iloc[0])
                elif len(full):
                    why = "rearranged form only as a low_weight species (%s, %.1f copies, J_adj %.2f)" % (
                        full["inferred"].iloc[0], full["copy_number"].iloc[0], full["j_adj"].iloc[0])
                else:
                    why = "rearranged form not returned"
                why = "all bases, no false adjacency; misses a reciprocal set of %d novel adjacencies; %s" % (len(miss), why)
            else:
                kind, why = "other", "see missing / extra adjacencies"
            lw_rows.append({k: r[k] for k in ("set", "case", "scenario", "class", "S", "stratum", "truth", "truth_bulk_copies",
                                              "paired", "paired_j_seg", "paired_j_adj", "exact_cycle")}
                           | {"n_truth_adjacencies": len(tj), "n_missing": len(miss), "n_extra": len(extra), "kind": kind,
                              "reason": why, "missing_adjacencies": fmt(miss), "extra_adjacencies": fmt(extra)})
    pd.DataFrame(lw_rows).to_csv(os.path.join(OUT, "jmatrix_not_1_1.tsv"), sep="\t", index=False, float_format="%.6g")
    print(pd.DataFrame(lw_rows)[["case", "truth", "paired", "kind", "reason"]].to_string())


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    p.add_argument("action", choices=["run", "compare", "reasons"])
    p.add_argument("--jobs", type=int, default=8)
    p.add_argument("--only", nargs="*", default=None, help="case names (tables are not written)")
    a = p.parse_args()
    if a.action == "run":
        run(a.jobs, a.only)
    elif a.action == "reasons":
        reasons()
    else:
        compare()


if __name__ == "__main__":
    main()

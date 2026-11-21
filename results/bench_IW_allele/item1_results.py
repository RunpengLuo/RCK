#!/usr/bin/env python3
"""bench_IW_allele item 1 (allelic siblings), NEXT_STEPS.md §5: the results tables.

Compares rck_mc (the method under test) with CoRAL and CE only (Yubi, 2026-11-21); rck_bulk, which only builds
the graph CoRAL / CE decompose, is left out everywhere. CoRAL / CE on sib_opp_mirror_s17006 ran on the
zero-read-patched graph (allelic_bench.py, INPUT PATCH). n = 4 per scenario, 2 per stratum: descriptive only.

Reads (all already written; run these first):
  bench_iw_allelic_table.tsv      allelic_bench.py score + collect
  sibling_relation.tsv            sibling_relation.py
  allele_table_allelic.tsv        allele_bench.py --set iw_allelic run + collect (post-hoc phasing)
  rck_style_table_allelic.tsv     rck_style.py --set iw_allelic run + collect (Part A metrics)
  the rck_mc runs (decompositions.tsv, evaluation_psi.tsv) and answers, the truth, rck_input/rck.fragments.tsv

Writes here, under item1/:
  table_by_scenario.tsv, table_by_stratum.tsv   §5.1 (mean, and n, per scenario [x stratum] x method variant)
  mirror_vs_shared.tsv                          §5.2 (sib_opp_mirror - sib_opp_shared, per method variant)
  psi.tsv, psi_summary.tsv                      §5.3 (rck_mc psi per truth species and clone)
  proportions_hap.tsv, proportions_hap_summary.tsv   §5.4 (proportions, species matched with haplotypes)
  allele_by_scenario.tsv, rck_style_by_scenario.tsv  §5.5 (post-hoc phasing comparison, Part A metrics)

    python item1_results.py     # local, seconds, no Gurobi
"""
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allelic_bench as ab  # noqa: E402  (configures bench_iw for iw_allelic)

bi = ab.bi
OUT = os.path.join(HERE, "item1")
SCENARIOS = ["sib_opp_mirror", "sib_opp_shared", "sib_mixed_mirror"]
VARIANTS = ["rck_mc", "coral:walks_kept", "coral:walks_dropped", "ce:walks_kept", "ce:walks_dropped"]
# bench_iw_allelic_table columns: haplotype-blind ones for every method; within-block phased ones are rck_mc's
# only (CoRAL / CE have no haplotypes)
METRICS = ["n_cycles", "n_strict", "n_missed", "n_false", "blind_exact_cycles", "blind_jseg_best_mean",
           "blind_jadj_best_mean", "species_detection_strict_blind_f1", "phased_exact_cycles", "phased_jseg_best_mean",
           "phased_jadj_best_mean", "species_detection_strict_phased_f1", "cn_recall", "cn_precision", "proportion_tv"]
NO_CYCLE_ZERO = ["blind_exact_cycles", "blind_jseg_best_mean", "blind_jadj_best_mean"]  # RESULTS.md convention
PSI_FAMILY_JACCARD = 0.9  # inferred species counted in a truth species' "family" (same molecule, split answers)


def case_meta():
    return pd.DataFrame([{"name": c["name"], "scenario": c["scenario"], "stratum": c["stratum"]} for c in bi.cases()])


def main_table():
    t = pd.read_csv(os.path.join(HERE, "bench_iw_allelic_table.tsv"), sep="\t")
    t = t[t["method_variant"].isin(VARIANTS)].copy()
    t = t.drop(columns=[c for c in ("stratum",) if c in t]).merge(case_meta()[["name", "stratum"]], on="name")
    ok = t["run_status"] == "ok"
    no_cycle = ok & (t["n_cycles"].fillna(0) == 0)
    for c in NO_CYCLE_ZERO:
        t.loc[no_cycle, c] = 0.0
    return t


def grouped(t, keys, metrics):
    g = t.groupby(keys, sort=False)
    mean = g[metrics].mean().round(3)
    mean.insert(0, "n", g.size())
    return mean.reset_index()


def order(frame, keys):
    frame = frame.copy()
    if "scenario" in keys:
        frame["_s"] = frame["scenario"].map({s: i for i, s in enumerate(SCENARIOS)})
    if "method_variant" in keys:
        frame["_v"] = frame["method_variant"].map({v: i for i, v in enumerate(VARIANTS)})
    cols = [c for c in ("_s", "stratum", "_v", "phasing") if c in frame]
    return frame.sort_values(cols).drop(columns=[c for c in ("_s", "_v") if c in frame])


def tables_51_52(t):
    rel = pd.read_csv(os.path.join(HERE, "sibling_relation.tsv"), sep="\t")
    t = t.merge(rel[["case", "disagree_bp", "matched"]].rename(columns={"case": "name"}), on="name", how="left")
    rck = t["method_variant"] == "rck_mc"
    t["sibling_matched"] = t["matched"].map({True: 1.0, False: 0.0, "True": 1.0, "False": 0.0}).where(rck)
    t["disagree_bp"] = t["disagree_bp"].astype(float).where(rck)
    metrics = METRICS + ["sibling_matched", "disagree_bp"]
    by_s = order(grouped(t, ["scenario", "method_variant"], metrics), ["scenario", "method_variant"])
    by_st = order(grouped(t, ["scenario", "stratum", "method_variant"], metrics), ["scenario", "method_variant"])
    by_s.to_csv(os.path.join(OUT, "table_by_scenario.tsv"), sep="\t", index=False)
    by_st.to_csv(os.path.join(OUT, "table_by_stratum.tsv"), sep="\t", index=False)
    m = by_s.set_index(["scenario", "method_variant"])
    rows = []
    for v in VARIANTS:
        a, b = m.loc[("sib_opp_mirror", v)], m.loc[("sib_opp_shared", v)]
        for k in METRICS:
            rows.append({"method_variant": v, "metric": k, "sib_opp_mirror": a[k], "sib_opp_shared": b[k],
                         "mirror_minus_shared": round(a[k] - b[k], 3) if pd.notna(a[k]) and pd.notna(b[k]) else None})
    pd.DataFrame(rows).to_csv(os.path.join(OUT, "mirror_vs_shared.tsv"), sep="\t", index=False)
    return t


def _answer_cycles(path, weights=None, strict=True):
    from rck.ecdna import aa_eval, aa_format as aaf, bench_metrics as bm

    records = aaf.read_cycles(path)
    cycles = aa_eval.species_cycles(records)
    chosen = aa_eval.strict_select(records, cycles, "auto")[0] if strict else cycles
    return bm.reweight(chosen, weights) if weights is not None else chosen


def _within_block(truth, inferred, blocks):
    """evaluate.score_decomposition's within-block pairing: blind pairing chooses one A/B swap per phase
    block, the inferred walks are relabelled, then paired again with haplotypes. -> (relabelled, pairing)."""
    from rck.ecdna import evaluate as ev

    by_t, by_c = {t.name: t for t in truth}, {c.name: c for c in inferred}
    provisional = ev._match(truth, inferred, False)
    swaps = ev.block_swaps([(by_t[tn], by_c[cn]) for tn, cn in provisional.items()], blocks)
    relabelled = [ev.relabel(c, blocks, swaps) for c in inferred]
    return relabelled, ev._match(truth, relabelled, True)


def table_53():
    """rck_mc psi per truth species and clone: the within-block one-to-one pair, and the 'family' (every
    inferred ecDNA species with haplotype Jaccard >= PSI_FAMILY_JACCARD to the truth species after the block
    relabelling: an answer that splits one molecule in two is counted as the sum)."""
    from rck.ecdna import bench_metrics as bm, evaluate as ev

    rows = []
    for case in bi.cases():
        run = bi.run_dir(case, "rck_mc")
        truth = ev.load_truth(case["truth"])
        blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
        weights = bm.multiclone_bulk_weights(run, bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv")))
        answer = os.path.join(bi.score_folder(case, "rck_mc"), "answer", "{c}_rck_mc_cycles.txt".format(c=case["name"]))
        inferred = _answer_cycles(answer, weights, strict=False)
        relabelled, pairing = _within_block(truth, inferred, blocks)
        dec = pd.read_csv(os.path.join(run, "decompositions.tsv"), sep="\t")
        dec = dec[dec["decomposition_rank"] == dec["decomposition_rank"].min()].drop_duplicates("species").set_index("species")
        clones = [c[len("psi_"):] for c in dec.columns if c.startswith("psi_")]
        tpsi = pd.read_csv(os.path.join(run, "evaluation_psi.tsv"), sep="\t").drop_duplicates(["truth", "clone"])
        tpsi = tpsi.set_index(["truth", "clone"])["truth_psi"]
        rel_by = {c.name: c for c in relabelled}
        for t in truth:
            paired = pairing.get(t.name)
            fam = [c.name for c in relabelled if c.name == paired or ev.jaccard(t.steps, c.steps, True) >= PSI_FAMILY_JACCARD]
            for k in clones:
                truth_psi = float(tpsi.get((t.name, k), float("nan")))
                inf = float(dec.loc[paired, "psi_" + k]) if paired in dec.index else 0.0
                fam_psi = float(sum(dec.loc[n, "psi_" + k] for n in fam if n in dec.index))
                rows.append({"name": case["name"], "scenario": case["scenario"], "stratum": case["stratum"],
                             "truth": t.name, "truth_role": "outer" if t.name == truth[0].name else "sibling",
                             "clone": k, "truth_psi": truth_psi, "inferred": paired,
                             "jaccard_hap": round(ev.jaccard(t.steps, rel_by[paired].steps, True), 4) if paired else None,
                             "inferred_psi": inf, "abs_err": abs(inf - truth_psi),
                             "family": ",".join(fam), "family_psi": fam_psi, "family_abs_err": abs(fam_psi - truth_psi)})
    psi = pd.DataFrame(rows)
    psi.to_csv(os.path.join(OUT, "psi.tsv"), sep="\t", index=False)
    g = psi.groupby(["scenario", "truth_role", "clone"], sort=False)
    s = g.agg(n=("abs_err", "size"), exact=("abs_err", lambda x: int((x == 0).sum())), mae=("abs_err", "mean"),
              max_err=("abs_err", "max"), family_exact=("family_abs_err", lambda x: int((x == 0).sum())),
              family_mae=("family_abs_err", "mean")).round(3).reset_index()
    order(s, ["scenario"]).to_csv(os.path.join(OUT, "psi_summary.tsv"), sep="\t", index=False)
    return psi


def _proportions(truth, inferred, blocks):
    """bench_metrics.proportion_scores with the species matched with haplotypes, within phase blocks."""
    from rck.ecdna import evaluate as ev

    relabelled, match = _within_block(truth, inferred, blocks)
    t_total, c_total = sum(t.weight for t in truth), sum(c.weight for c in inferred)
    p_t = {t.name: t.weight / t_total for t in truth}
    p_c = {c.name: c.weight / c_total if c_total else 0.0 for c in inferred}
    tv, errors = 0.0, []
    for t in truth:
        c = match.get(t.name)
        if c:
            errors.append(abs(p_c[c] - p_t[t.name]))
            tv += errors[-1]
        else:
            tv += p_t[t.name]
    tv += sum(p_c[c.name] for c in inferred if c.name not in set(match.values()))
    return {"proportion_tv_hap": tv / 2.0, "proportion_mae_matched_hap": sum(errors) / len(errors) if errors else None,
            "n_matched_hap": len(errors), "n_inferred": len(inferred), "n_truth": len(truth)}


def table_54():
    """Proportions matched with haplotypes: rck_mc's own phasing; CoRAL / CE phased post hoc (allele_bench.py
    baf_major and the truth_major upper bound), strict species as bench_metrics.score."""
    import allele_bench as alb
    from rck.ecdna import bench_metrics as bm, evaluate as ev

    alb.use_set("iw_allelic")
    blind = pd.read_csv(os.path.join(HERE, "bench_iw_allelic_table.tsv"), sep="\t").set_index(["name", "method_variant"])
    rows = []
    for case in bi.cases():
        truth = ev.load_truth(case["truth"])
        blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
        for method, variant in alb.VARIANTS:
            mv = bi.method_variant(method, variant)
            for phasing in (["native"] if method == "rck_mc" else ["baf_major", "truth_major"]):
                folder = alb.task_dir(case, method, variant, phasing)
                cycles = os.path.join(folder, "cycles.txt") if phasing != "native" else alb.answer_files(case, method, variant)[0]
                weights = None
                if method == "rck_mc":
                    weights = bm.multiclone_bulk_weights(bi.run_dir(case, method), bm.read_proportions(
                        os.path.join(case["rck_input"], "clone_proportions.tsv")))
                row = {"name": case["name"], "scenario": case["scenario"], "stratum": case["stratum"],
                       "method_variant": mv, "phasing": phasing}
                row.update(_proportions(truth, _answer_cycles(cycles, weights), blocks))
                row["proportion_tv_blind"] = blind.loc[(case["name"], mv), "proportion_tv"]
                rows.append(row)
    p = pd.DataFrame(rows)
    p.to_csv(os.path.join(OUT, "proportions_hap.tsv"), sep="\t", index=False)
    s = grouped(p, ["scenario", "method_variant", "phasing"], ["proportion_tv_blind", "proportion_tv_hap",
                                                                "proportion_mae_matched_hap", "n_matched_hap"])
    order(s, ["scenario", "method_variant"]).to_csv(os.path.join(OUT, "proportions_hap_summary.tsv"), sep="\t", index=False)
    return p


def table_55():
    meta = case_meta()
    a = pd.read_csv(os.path.join(HERE, "allele_table_allelic.tsv"), sep="\t").merge(meta, on="name")
    a = a[a["method_variant"].isin(VARIANTS)]
    am = ["as_cn_recall_block", "as_cn_precision_block", "as_cn_mae_block", "major_allele_acc_block", "baf_mae_block",
          "blind_cn_recall", "phased_exact_cycles_block", "phased_jseg_best_mean_block", "phased_jadj_best_mean_block",
          "mixed_bp_frac"]
    order(grouped(a, ["scenario", "method_variant", "phasing"], am), ["scenario", "method_variant"]).to_csv(
        os.path.join(OUT, "allele_by_scenario.tsv"), sep="\t", index=False)
    r = pd.read_csv(os.path.join(HERE, "rck_style_table_allelic.tsv"), sep="\t").merge(meta, on="name")
    r = r[r["method_variant"].isin(VARIANTS + ["truth"])]
    rm = ["adj_phased_fnr", "adj_phased_fdr", "adj_blind_fnr", "adj_blind_fdr", "cn_dist_rel", "is_is_ext",
          "is_is_hom_same", "is_is_hom", "is_is_recip", "input_input_fpr", "clone_cn_dist_rel_ecdna"]
    rr = grouped(r, ["scenario", "method_variant", "phasing"], rm)
    rr["_v"] = rr["method_variant"].map({v: i for i, v in enumerate(["truth"] + VARIANTS)})
    rr["_s"] = rr["scenario"].map({s: i for i, s in enumerate(SCENARIOS)})
    rr.sort_values(["_s", "_v", "phasing"]).drop(columns=["_s", "_v"]).to_csv(
        os.path.join(OUT, "rck_style_by_scenario.tsv"), sep="\t", index=False)


def main():
    import logging
    logging.disable(logging.INFO)
    os.makedirs(OUT, exist_ok=True)
    tables_51_52(main_table())
    table_53()
    table_54()
    table_55()
    print("-> " + OUT)


if __name__ == "__main__":
    main()

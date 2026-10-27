"""Collect bench_nested: re-rank every case over the λ grid, choose λ, write the summaries.

    python collect.py              # rerank_validation.tsv, bench_sweep.tsv, pipeline.tsv, README.md results
    python collect.py --validate   # only the re-ranking check against results/step_sweep/sweep.tsv

Reads runs/<case>/lambda_0/ (a λ = 0 run) and runs/<case>/pool/ (a --keep-pool run, re-ranked
with rck.ecdna.rerank for every λ > 0), and amplicon-allele-analysis
simulations/subgraphs/bench_nested/ for the subgraph-stage checks. Missing runs are skipped.
BENCH_NESTED_ROOT and BENCH_NESTED_OUT are read as in run_case.py; BENCH_NESTED_RESULTS writes the
tables and the README results section to another folder (README.md there holds only the results).
"""
import argparse
import json
import os
import re
import sys

import numpy as np
import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
RESULTS = os.environ.get("BENCH_NESTED_RESULTS", HERE)
import run_case  # noqa: E402  (puts the RCK root on sys.path)
from rck.ecdna.evaluate import junction_scores, read_inferred, read_truth_table  # noqa: E402
from rck.ecdna.rerank import rerank  # noqa: E402

sys.path.insert(0, run_case.SIM)
import bench_nested as aaa_bench  # noqa: E402  (amplicon-allele-analysis simulations/bench_nested.py)
import verify_subgraph  # noqa: E402

SWEEP = os.path.join(run_case.RCK, "results", "step_sweep", "sweep.tsv")
LAMBDAS = run_case.LAMBDAS
TOL_REGRESSION = 0.005  # a control or demo case may lose at most this much mean composite against λ = 0
TOL_PLATEAU = 0.001  # λ within this of the best mean composite over all cases are on the plateau
FAIL_COMPOSITE = 0.95  # a bench case below this mean composite (or merged, or a wrong species count) is a failure
JUNCTION_MATCH = 0.5  # an inferred species "carries" a truth species with at least this junction recall
COLUMNS = ["case", "kind", "scenario", "step_penalty", "source", "n_species", "n_truth", "weights", "steps", "residual",
           "composites", "mean_composite", "merged", "merged_note", "ambiguous", "decided_by", "gap", "n_alternatives",
           "runtime_s"]


# ---------------------------------------------------------------------------------------------- #
# per-decomposition scores
# ---------------------------------------------------------------------------------------------- #

def merged_flag(folder, rank):
    """(merged, note) for decomposition ``rank`` of the run (or pool) in ``folder``: an inferred
    species carries >= 2 truth species, or a truth species is missed while an inferred species
    carries it. "Carries" is junction recall >= JUNCTION_MATCH: nested species share their
    reference segments by construction, so segment coverage cannot tell them apart, while in
    these simulations they share no junction."""
    evaluation = pd.read_csv(os.path.join(folder, "evaluation.tsv"), sep="\t")
    mine = evaluation[evaluation["decomposition_rank"] == rank]
    scope = set(mine["truth"])
    truth = [t for t in read_truth_table(os.path.join(folder, "truth_cycles.tsv")) if t.name in scope]
    inferred = read_inferred(os.path.join(folder, "decompositions.tsv")).get(rank, {"cycles": []})["cycles"]
    carries = {c.name: [t.name for t in truth
                        if junction_scores(t.steps, c.steps, 10_000, False, t.circular, c.circular)[0] >= JUNCTION_MATCH]
               for c in inferred}
    notes = ["{c} carries {t}".format(c=c, t="+".join(ts)) for c, ts in carries.items() if len(ts) >= 2]
    missed = list(mine.loc[mine["inferred"].isna(), "truth"])
    for t in missed:
        holders = [c for c, ts in carries.items() if t in ts]
        if holders:
            notes.append("{t} missed, carried by {c}".format(t=t, c="+".join(holders)))
    return bool(notes), "; ".join(notes)


def run_row(case, folder):
    """The rank-1 row of an actual run (the λ = 0 runs)."""
    summary = json.load(open(os.path.join(folder, "summary.json")))
    top = summary["decompositions"][0]
    significant = [s for s in top["species"] if s["status"] != "low_weight"]
    matrix = pd.read_csv(os.path.join(folder, "recovery_matrix.tsv"), sep="\t")
    matrix = matrix[matrix["in_scope"].astype(bool)]
    amb = summary.get("ambiguity") or {}
    merged, note = merged_flag(folder, 1)
    return {"n_species": len(significant), "weights": ",".join("%g" % s["weight"] for s in significant),
            "steps": top["steps"], "residual": top["residual"],
            "composites": ",".join("%s=%.3f" % (r.truth, r.composite) for r in matrix.itertuples()),
            "mean_composite": round(float(matrix["composite"].mean()), 4) if len(matrix) else None,
            "n_truth": len(matrix), "merged": merged, "merged_note": note, "ambiguous": amb.get("ambiguous"),
            "decided_by": amb.get("decided_by"), "gap": amb.get("gap"), "n_alternatives": amb.get("n_alternatives"),
            "runtime_s": summary.get("runtime_s")}


def pool_rows(case, folder, lambdas):
    """One row per λ from the kept pool (rck.ecdna.rerank), with the merged flag."""
    pool = os.path.join(folder, "pool")
    runtime = json.load(open(os.path.join(folder, "summary.json"))).get("runtime_s")
    n_truth = len(pd.read_csv(os.path.join(pool, "evaluation.tsv"), sep="\t")["truth"].unique())
    rows = []
    for r in rerank(folder, lambdas).to_dict("records"):
        merged, note = merged_flag(pool, int(r["pool_index"]) + 1)
        r.update({"merged": merged, "merged_note": note, "n_truth": n_truth, "runtime_s": runtime})
        rows.append(r)
    return rows


def case_rows(case):
    base = os.path.join(run_case.OUT, case["name"])
    rows = []
    head = {"case": case["name"], "kind": case["kind"], "scenario": case["scenario"]}
    if run_case.done(os.path.join(base, "lambda_0")):
        rows.append(dict(head, step_penalty=0.0, source="run", **run_row(case, os.path.join(base, "lambda_0"))))
    if run_case.done(os.path.join(base, "pool"), pool=True):
        for r in pool_rows(case, os.path.join(base, "pool"), [x for x in LAMBDAS if x > 0]):
            rows.append(dict(head, **r))
    return rows


# ---------------------------------------------------------------------------------------------- #
# validation against the step sweep
# ---------------------------------------------------------------------------------------------- #

def validate(table):
    """Rank-1 composites and ambiguity fields, re-ranked (λ > 0, and λ = 0 emulated) or run
    (λ = 0), against results/step_sweep/sweep.tsv."""
    sweep = pd.read_csv(SWEEP, sep="\t")
    rows = []
    for case in sweep["case"].unique():
        base = os.path.join(run_case.OUT, case)
        ours = {}
        if run_case.done(os.path.join(base, "pool"), pool=True):
            for r in rerank(os.path.join(base, "pool"), sorted(sweep["step_penalty"].unique())).to_dict("records"):
                ours[(r["step_penalty"], "pool" if r["step_penalty"] > 0 else "pool (λ=0 emulated)")] = r
        mine = table[(table["case"] == case) & (table["source"] == "run")]
        for r in mine.to_dict("records"):
            ours[(0.0, "run")] = r
        for (lam, method), r in sorted(ours.items()):
            s = sweep[(sweep["case"] == case) & (np.isclose(sweep["step_penalty"], lam))]
            if not len(s):
                continue
            s = s.iloc[0]
            same_comp = str(r.get("composites")) == str(s["composites"])
            same_amb = (str(r.get("decided_by") or "") == str(s["decided_by"] if pd.notna(s["decided_by"]) else "")
                        and (pd.isna(s["gap"]) and r.get("gap") is None
                             or (r.get("gap") is not None and pd.notna(s["gap"]) and abs(float(r["gap"]) - float(s["gap"])) < 1e-3)))
            rows.append({"case": case, "step_penalty": lam, "method": method, "sweep_composites": s["composites"],
                         "composites": r.get("composites"), "composites_match": same_comp,
                         "sweep_decided_by": s["decided_by"], "decided_by": r.get("decided_by"), "sweep_gap": s["gap"],
                         "gap": r.get("gap"), "ambiguity_match": same_amb})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RESULTS, "rerank_validation.tsv"), sep="\t", index=False)
    return out


# ---------------------------------------------------------------------------------------------- #
# choice of λ
# ---------------------------------------------------------------------------------------------- #

def is_control(row):
    return row["kind"] == "regression" or str(row["scenario"]).startswith("ctrl_")


def choose(table):
    """The smallest λ on the plateau of best mean composite over all cases with no control or
    demo case below its λ = 0 score by more than TOL_REGRESSION; and the λ that give the same
    rank-1 answers as it."""
    complete = table.groupby("case")["step_penalty"].nunique()
    cases = list(complete[complete == len(LAMBDAS)].index)
    t = table[table["case"].isin(cases)]
    overall = t.groupby("step_penalty")["mean_composite"].mean()
    base = t[t["step_penalty"] == 0].set_index("case")["mean_composite"]
    regress = {}
    for lam in LAMBDAS:
        rows = t[t["step_penalty"] == lam]
        regress[lam] = [r["case"] for r in rows.to_dict("records")
                        if is_control(r) and r["mean_composite"] < base[r["case"]] - TOL_REGRESSION]
    best = overall.max()
    plateau = [lam for lam in LAMBDAS if overall[lam] >= best - TOL_PLATEAU]
    eligible = [lam for lam in plateau if not regress[lam]]
    chosen = min(eligible) if eligible else None
    same = []
    if chosen is not None:
        key = lambda lam: t[t["step_penalty"] == lam].set_index("case")[["n_species", "composites"]].astype(str).sort_index()
        same = [lam for lam in LAMBDAS if key(lam).equals(key(chosen))]
    return {"cases": cases, "overall": overall, "regress": regress, "plateau": plateau, "chosen": chosen, "same": same}


# ---------------------------------------------------------------------------------------------- #
# pipeline table (bench cases at the chosen λ)
# ---------------------------------------------------------------------------------------------- #

def flip_check(folder, truth_dir):
    """Block flips of the flow against the truth orientation, for phase blocks over ecDNA segments
    whose truth orientation is clear (purity >= 0.9) and that the flow resolved."""
    path = os.path.join(folder, "block_flips.tsv")
    if not os.path.exists(path):
        return {}
    flips = pd.read_csv(path, sep="\t")
    truth = verify_subgraph.truth_block_orientation(truth_dir, 1_000_000)
    iv = aaa_bench.merged_intervals(verify_subgraph.load_truth_segments(truth_dir, 0.0))
    flips = flips[[aaa_bench.block_overlaps(b, iv) for b in flips["phase_block"]]]
    flips = flips.merge(truth[["truth_orientation", "truth_purity"]], left_on="phase_block", right_index=True, how="left")
    eff = np.where(flips["orientation"] == 0, 1, flips["orientation"]) * np.where(flips["flip"].astype(bool), -1, 1)
    known = (flips["status"] != "unresolved") & (flips["truth_purity"] >= 0.9)
    return {"ecdna_blocks_flipped": int(flips["flip"].astype(bool).sum()),
            "ecdna_blocks_unresolved": int((flips["status"] == "unresolved").sum()),
            "ecdna_blocks_after_flow_correct": int(((eff == flips["truth_orientation"]) & known).sum()),
            "ecdna_blocks_after_flow_wrong": int(((eff != flips["truth_orientation"]) & known).sum())}


def failure_stage(row, sub):
    """Which stage explains a failing bench case: the subgraph when it misses a true junction or
    covers under 80% of a true segment, else RCK (with what went wrong)."""
    bad = (row["mean_composite"] is None or row["mean_composite"] < FAIL_COMPOSITE or row["merged"]
           or row["n_species"] != row["n_truth"])
    if not bad:
        return "", ""
    if sub.get("subgraph_fail"):
        return "subgraph", "missing junctions {j}; min segment coverage {c}%".format(
            j=sub.get("junctions_missing") or "none", c=sub.get("segment_coverage_min"))
    why = []
    if row["merged"]:
        why.append("merged: " + row["merged_note"])
    if row["n_species"] != row["n_truth"]:
        why.append("{n} species for {t} true".format(n=row["n_species"], t=row["n_truth"]))
    if sub.get("ecdna_blocks_after_flow_wrong"):
        why.append("{w} ecDNA phase block(s) wrong after the flow".format(w=sub["ecdna_blocks_after_flow_wrong"]))
    if sub.get("imputed_used"):
        why.append("imputed edges " + sub["imputed_used"])
    why.append("composites " + str(row["composites"]))
    return "RCK", "; ".join(why)


def pipeline(table, chosen):
    rows = []
    for case in [c for c in run_case.cases() if c["kind"] == "bench"]:
        r = table[(table["case"] == case["name"]) & np.isclose(table["step_penalty"], chosen)]
        if not len(r):
            continue
        r = r.iloc[0].to_dict()
        folder = os.path.join(run_case.OUT, case["name"], "pool")
        sub_path = os.path.join(run_case.SUBGRAPHS, case["name"], "subgraph_check.json")
        sub = json.load(open(sub_path)) if os.path.exists(sub_path) else {}
        summary = json.load(open(os.path.join(folder, "summary.json")))
        flow = {"flow_objective": round(summary["flow"]["objective"], 3), "imputed_used": ",".join(summary["flow"]["imputed_used"]),
                "n_imputed_used": len(summary["flow"]["imputed_used"])}
        flow.update(flip_check(folder, case["truth"]))
        stage, why = failure_stage(r, dict(sub, **flow))
        rows.append({"case": case["name"], "scenario": case["scenario"], "step_penalty": chosen,
                     "segments_covered": sub.get("segments_covered"), "segments_total": sub.get("segments_total"),
                     "segment_coverage_min": sub.get("segment_coverage_min"), "junctions_included": sub.get("junctions_included"),
                     "junctions_total": sub.get("junctions_total"), "subgraph_fail": sub.get("subgraph_fail"),
                     "ecdna_blocks_correct": sub.get("ecdna_blocks_correct"), "ecdna_blocks_wrong": sub.get("ecdna_blocks_wrong"),
                     "ecdna_blocks_unoriented": sub.get("ecdna_blocks_unoriented"), "theta_used": sub.get("theta_used"),
                     "theta_fitted": sub.get("theta_fitted"), **flow,
                     "n_species": r["n_species"], "n_truth": r["n_truth"], "composites": r["composites"],
                     "mean_composite": r["mean_composite"], "merged": r["merged"], "ambiguous": r["ambiguous"],
                     "decided_by": r["decided_by"], "gap": r["gap"], "failure_stage": stage, "failure": why})
    out = pd.DataFrame(rows)
    out.to_csv(os.path.join(RESULTS, "pipeline.tsv"), sep="\t", index=False)
    return out


# ---------------------------------------------------------------------------------------------- #
# README results section
# ---------------------------------------------------------------------------------------------- #

def fmt(x, digits=3):
    return "" if x is None or (isinstance(x, float) and np.isnan(x)) else ("{:.%df}" % digits).format(x)


def results_markdown(table, choice, validation, pipe):
    L = []
    if len(validation):
        pool = validation[validation["method"] == "pool"]
        run0 = validation[validation["method"] == "run"]
        emu = validation[validation["method"].str.contains("emulated")]
        L.append("**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): "
                 "re-ranked λ > 0 composites match at {a}/{b} (case, λ), ambiguity fields at {c}/{b}; "
                 "the λ = 0 runs match at {d}/{e}; λ = 0 emulated from the pool matches at {f}/{g}.".format(
                     a=int(pool["composites_match"].sum()), b=len(pool), c=int(pool["ambiguity_match"].sum()),
                     d=int(run0["composites_match"].sum()), e=len(run0), f=int(emu["composites_match"].sum()), g=len(emu)))
        L.append("")
    if choice["chosen"] is None:
        L.append("**No λ qualifies**: every λ on the plateau regresses a control or demo case (table below).")
    else:
        L.append("**Recommended `--step-penalty`: {c:g}**, the smallest λ on the plateau of best mean composite over "
                 "all {n} cases ({p}), with no control or demo case more than {t} below its λ = 0 score. "
                 "λ giving the same rank-1 answers on every case: {s}. The default stays 0.".format(
                     c=choice["chosen"], n=len(choice["cases"]), p=", ".join("%g" % x for x in choice["plateau"]),
                     t=TOL_REGRESSION, s=", ".join("%g" % x for x in choice["same"])))
    L.append("")
    L.append("| λ | " + " | ".join("%g" % x for x in LAMBDAS) + " |")
    L.append("|---|" + "---|" * len(LAMBDAS))
    L.append("| mean composite, all cases | " + " | ".join(fmt(choice["overall"].get(x), 4) for x in LAMBDAS) + " |")
    L.append("| control / demo cases below λ = 0 | " + " | ".join(", ".join(choice["regress"][x]) or "none" for x in LAMBDAS) + " |")
    L.append("")
    chosen = choice["chosen"] if choice["chosen"] is not None else max(LAMBDAS)
    t = table[table["case"].isin(choice["cases"])]
    by = t.groupby(["scenario", "step_penalty"])["mean_composite"].mean().unstack()
    merged = t.groupby(["scenario", "step_penalty"])["merged"].sum().unstack()
    L.append("Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):")
    L.append("")
    L.append("| scenario | cases | λ = 0 | λ = {c:g} | merged at 0 | merged at {c:g} |".format(c=chosen))
    L.append("|---|---|---|---|---|---|")
    for scenario in by.index:
        n = t[(t["scenario"] == scenario) & (t["step_penalty"] == 0)]["case"].nunique()
        L.append("| {s} | {n} | {a} | {b} | {m0:g} | {m1:g} |".format(s=scenario, n=n, a=fmt(by.loc[scenario, 0.0]),
                                                                    b=fmt(by.loc[scenario, chosen]), m0=merged.loc[scenario, 0.0],
                                                                    m1=merged.loc[scenario, chosen]))
    L.append("")
    if len(pipe):
        fails = pipe[pipe["failure_stage"] != ""]
        L.append("Remaining failures at λ = {c:g} (mean composite < {f}, merged, or a wrong species count; `pipeline.tsv`): "
                 "{n} of {m} bench cases.".format(c=chosen, f=FAIL_COMPOSITE, n=len(fails), m=len(pipe)))
        L.append("")
        if len(fails):
            L.append("| case | stage | why |")
            L.append("|---|---|---|")
            for r in fails.itertuples():
                L.append("| {c} | {s} | {w} |".format(c=r.case, s=r.failure_stage, w=str(r.failure).replace("|", "/")))
            L.append("")
        subfail = pipe[pipe["subgraph_fail"].astype(str) == "True"]
        L.append("Subgraph stage: {n} of {m} subgraphs miss a true junction or cover under 80% of a true segment{x}.".format(
            n=len(subfail), m=len(pipe), x=" (" + ", ".join(subfail["case"]) + ")" if len(subfail) else ""))
    runtime = t.drop_duplicates(["case", "source"])["runtime_s"].astype(float).sum()
    L.append("")
    L.append("Compute: {s:.0f} s of RCK wall time over {n} cases (one λ = 0 run and one pool run each, 4 threads).".format(
        s=runtime, n=t["case"].nunique()))
    return "\n".join(L)


def write_readme(text):
    path = os.path.join(RESULTS, "README.md")
    start, end = "<!-- results:start -->", "<!-- results:end -->"
    os.makedirs(RESULTS, exist_ok=True)
    readme = open(path).read() if os.path.exists(path) else start + "\n" + end + "\n"
    readme = re.sub(re.escape(start) + ".*?" + re.escape(end), start + "\n" + text + "\n" + end, readme, flags=re.S)
    with open(path, "w") as handle:
        handle.write(readme)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--validate", action="store_true", help="only check the re-ranking against the step sweep")
    args = parser.parse_args()
    rows = []
    for case in run_case.cases():
        if args.validate and case["name"] not in pd.read_csv(SWEEP, sep="\t")["case"].unique():
            continue
        rows.extend(case_rows(case))
    table = pd.DataFrame(rows, columns=COLUMNS)
    validation = validate(table)
    if len(validation):
        print(validation.to_string(index=False))
    if args.validate:
        return 0
    table.to_csv(os.path.join(RESULTS, "bench_sweep.tsv"), sep="\t", index=False)
    choice = choose(table)
    pipe = pipeline(table, choice["chosen"] if choice["chosen"] is not None else max(LAMBDAS))
    text = results_markdown(table, choice, validation, pipe)
    write_readme(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())

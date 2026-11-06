"""Score the graph-route head-to-head (priority 1 d) with rck.ecdna.aa_eval. Offline, no Gurobi.

    python score_h2h.py [--jobs 16]

Every run of graph_route.list (CoRAL, CycleExtractor, RCK-ecDNA on the truth graph_reads_paths and
inferred graphs) and the native arm (RCK-ecDNA's bench runs at lambda 0.2, re-scored from the AA
exports of results/bench_nested/diag/aa_regression.py) goes through the same call,
aa_eval.score, with the scorer's defaults: strict rule `auto` (CoRAL / CE: cn-frac:0.1; RCK-ecDNA:
`status`), species tolerance 10 kb, graph tolerance 100 bp. Cross-method comparisons use the
haplotype-blind columns (composite_h0, the blind one-to-one matching), since CoRAL and CE do not
phase and RCK-ecDNA's graph-route runs are total-only.

  * cycles:  CoRAL coral_amplicon<N>_cycles.txt (on the flanked graph); CE ce_amplicon<N>_cycles.txt;
             RCK-ecDNA the pool's rank 1 re-ranked at lambda 0.2 + 2 same-fit alternatives
             (aa_export.export, written to the score folder), as the bench's lambda 0.2 results;
  * graph:   the input graph every method was given (for CoRAL the unflanked one), for
             cn_explained only; graph-level metrics are computed once per input graph
             (h2h_graph_metrics.tsv), not per method;
  * paths:   not passed (subwalks_satisfied left out): aa_eval.read_paths needs chr:pos EIDs and the
             inferred graphs' paths.tsv carry the subgraph's imputed_ext_* ids.

Linear walks. CoRAL 3.1.1 and CycleExtractor write a linear walk as `Path=k;...;Segments=0+,...,0-`,
which aa_format.read_cycles skips (it reads `Cycle=` lines only), while rck-ecdna writes its linear
walks AA's way, as `Cycle=` lines with segment 0, which are scored. So that every method's linear
walks are treated alike, each answer is first copied to <score dir>/input/ with `Path=` read as
`Cycle=` (the default: linear walks scored, as the scorer documents for AA-style paths), or with
every linear walk removed (`--drop-paths`, a sensitivity run: outputs with suffix _cycles_only).

Per truth species (strict, blind matching): `missed` = no strict cycle matched; `merged` = a strict
cycle carries >= 0.5 of its junctions (recall, blind, 10 kb) and >= 0.5 of another truth species'
junctions (collect.py's merged rule, haplotype-blind). `n_false` counts unmatched strict cycles,
`n_false_linear` those of them that are paths (CoRAL / CE background walks through the flanks).

Writes, here: h2h_species.tsv, h2h_table.tsv (case x graph x method x mode), h2h_scenario.tsv,
h2h_graph_metrics.tsv; per run: $G/<case>/score/<graph>/<method>/amplicon<N>/ (aa_eval.write).
"""
import argparse
import json
import os
import sys
import traceback
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import graph_route as gr  # noqa: E402  (puts the RCK root and results/bench_nested on sys.path)
import run_case  # noqa: E402

from rck.ecdna import aa_eval, aa_export  # noqa: E402
from rck.ecdna import aa_format as aaf  # noqa: E402
from rck.ecdna import evaluate as ev  # noqa: E402
from rck.ecdna.rerank import rerank  # noqa: E402

NATIVE = "/n/fs/yubi-iw/runs/bench_nested_aa"
RUN_CASE = {"demo6": "demo6_union", "demo7": "demo7_totalonly"}  # graph-route case -> run_case case (truth, sample)
NATIVE_CASE = {"demo6": "demo6_union", "demo7": "demo7_blockflip"}  # as native_aa_metrics.tsv
JUNCTION_MATCH = 0.5
CASES = {c["name"]: c for c in run_case.cases()}


def truth_of(case):
    c = CASES[RUN_CASE.get(case, case)]
    opts = c["options"]
    sample = opts[opts.index("--truth-sample") + 1] if "--truth-sample" in opts else None
    return ev.load_truth(c["truth"], None, sample), c


def scenario_of(case):
    return "demo" if case.startswith("demo") else CASES[case]["scenario"]


def merged_flags(truth, strict):
    """{truth name: merged} on haplotype-stripped cycles (collect.py's rule, blind)."""
    t0 = [aa_eval.strip_haplotypes(t) for t in truth]
    carries = {}
    for c in strict:
        c0 = aa_eval.strip_haplotypes(c)
        carries[c.name] = [t.name for t in t0
                           if ev.junction_scores(t.steps, c0.steps, aa_eval.SPECIES_TOL, False, t.circular, c0.circular)[0] >= JUNCTION_MATCH]
    return {t.name: any(t.name in ts and len(ts) >= 2 for ts in carries.values()) for t in truth}


def normalise_cycles(src, folder, drop_paths):
    """Copy a cycles file to folder/: CoRAL's `Path=` lines read as `Cycle=` (they carry AA's
    0+ ... 0- segments, so the reader marks them linear), or every linear walk removed."""
    os.makedirs(folder, exist_ok=True)
    dst = os.path.join(folder, os.path.basename(src))
    out = []
    with open(src) as handle:
        for line in handle:
            if line.startswith(("Path=", "Cycle=")):
                fields = dict(f.split("=", 1) for f in line.strip().split(";") if "=" in f)
                linear = line.startswith("Path=") or fields.get("Segments", "").split(",")[0][:-1] == "0"
                if linear and drop_paths:
                    continue
                if line.startswith("Path="):
                    line = "Cycle=" + line[len("Path="):]
            out.append(line)
    with open(dst, "w") as handle:
        handle.writelines(out)
    return dst


def one(task):
    kind, case, arm, amp, method, cycles, alts, graph, outdir, drop_paths = task
    try:
        truth, _ = truth_of(case)
        if method == "rck_ecdna":
            index = int(rerank(cycles, [gr.SCORE_LAMBDA])["pool_index"].iloc[0])
            paths = aa_export.export(cycles, outdir, index, alternatives=2)
            cycles, alts = paths["cycles"], paths["alt_cycles"]
        folder = os.path.join(outdir, "input")
        cycles = normalise_cycles(cycles, folder, drop_paths)
        alts = [normalise_cycles(a, folder, drop_paths) for a in alts]
        result = aa_eval.score(cycles, truth, graph, alts)
        aa_eval.write(outdir, result)
        species = result.species.copy()
        flags = merged_flags(result.truth, result.strict)
        species["merged"] = species["truth"].map(flags).where(species["mode"] == "strict")
        matched = set(species.loc[species["mode"] == "strict", "inferred_h0"].dropna())
        unmatched = [c for c in result.strict if c.name not in matched]
        species = species.assign(case=case, graph=arm, method=method, amplicon=amp, scenario=scenario_of(case),
                                 n_strict=len(result.strict), n_strict_linear=sum(not c.circular for c in result.strict),
                                 n_false_h0=len(unmatched), n_false_linear=sum(not c.circular for c in unmatched),
                                 n_candidates=len(result.candidates), strict_rule=result.info["strict_rule"])
        sample = result.sample.set_index("mode")
        for col in ("cn_explained_top1", "cn_explained_top3"):
            species[col] = species["mode"].map(sample[col])
        return species, None
    except Exception:  # noqa: BLE001 - report every failing run, keep going
        return None, "{k} {c} {g} a{a} {m}\n{e}".format(k=kind, c=case, g=arm, a=amp, m=method, e=traceback.format_exc())


def tasks(drop_paths=False):
    out = []
    top = "score_cycles_only" if drop_paths else "score"
    for case, arm, amp, method, graph, _ri, rundir, _t in gr.read_list():
        score_dir = os.path.join(gr.G, case, top, arm, method, "amplicon" + amp)
        cycles = rundir if method == "rck_ecdna" else gr.cycles_file(method, rundir, int(amp))
        out.append(("h2h", case, arm, int(amp), method, cycles, [], graph, score_dir, drop_paths))
    for case, _sub, _truth in gr.roundtrip.cases():
        folder = os.path.join(NATIVE, NATIVE_CASE.get(case, case), "lambda_0.2")
        with open(os.path.join(folder, "aa_eval.json")) as handle:  # the files aa_regression.py scored
            info = json.load(handle)["info"]
        out.append(("native", case, "native", 1, "native", info["cycles"], info["alt_cycles"], info["graph"],
                    os.path.join(gr.G, case, top, "native"), drop_paths))
    return out


def graph_rows():
    """Graph-level metrics once per input graph: truth graph_reads_paths amplicons and the inferred
    export (against every truth mixture amplicon its in-scope species belong to), plus native."""
    rows = []
    seen = set()
    for case, arm, amp, _m, graph, _ri, _rd, _t in gr.read_list():
        if (case, arm, amp) in seen:
            continue
        seen.add((case, arm, amp))
        truth, c = truth_of(case)
        mixture = aa_eval.load_mixture(aa_eval.default_mixture_dir(c["truth"]))
        segment_cn = graph.replace("_graph.txt", "_segment_cn.tsv") if arm == "inferred" else None
        cn = pd.read_csv(segment_cn, sep="\t") if segment_cn and os.path.exists(segment_cn) else None
        g = aaf.read_graph(graph)
        for m in mixture:
            if arm != "inferred" and m.index != int(amp):
                continue
            r = aa_eval.graph_metrics(m, g, cn, aa_eval.GRAPH_TOL)
            rows.append(dict(case=case, scenario=scenario_of(case), graph=arm, input_graph=graph, **r))
    return pd.DataFrame(rows)


def summarise(species):
    keys = ["case", "scenario", "graph", "method", "mode"]
    sp = species.copy()
    sp["exact_cycle"] = sp["exact_cycle"].fillna(0).astype(int)
    sp["missed"] = sp["inferred_h0"].isna()  # strict: no blind one-to-one match
    best = sp["mode"] == "best_of"
    sp.loc[best, "missed"] = sp.loc[best, "composite_h0"].fillna(0) <= 0  # aa_eval's best-of rule
    sp["missed_or_merged"] = sp["missed"].astype(bool) | sp["merged"].fillna(False).astype(bool)
    per_run = sp.drop_duplicates(keys + ["amplicon"])  # run-level counts: once per amplicon run
    run_counts = per_run.groupby(keys)[["n_strict", "n_strict_linear", "n_false_h0", "n_false_linear"]].sum()
    table = sp.groupby(keys).agg(
        n_truth=("truth", "size"), mean_composite_h0=("composite_h0", "mean"),
        mean_cycle_overlap_coral_code=("cycle_overlap_coral_code", "mean"), mean_cyclic_lcs_coral_code=("cyclic_lcs_coral_code", "mean"),
        n_exact=("exact_cycle", "sum"), n_missed=("missed", "sum"), n_merged=("merged", lambda x: int(x.fillna(False).astype(bool).sum())),
        n_missed_or_merged=("missed_or_merged", "sum"),
        composites_h0=("composite_h0", lambda x: ",".join("%.3f" % v for v in x)))
    table = table.join(run_counts).reset_index()
    table.loc[table["mode"] == "best_of", ["n_merged", "n_missed_or_merged", "n_false_h0", "n_false_linear"]] = None
    return table


def scenario_summary(table):
    t = table.copy()
    agg = t.groupby(["scenario", "graph", "method", "mode"]).agg(
        n_cases=("case", "size"), n_truth=("n_truth", "sum"), mean_composite_h0=("mean_composite_h0", "mean"),
        n_exact=("n_exact", "sum"), n_missed=("n_missed", "sum"), n_merged=("n_merged", "sum"),
        n_missed_or_merged=("n_missed_or_merged", "sum"),
        n_false_circular=("n_false_h0", lambda x: x.sum()), n_false_linear=("n_false_linear", "sum")).reset_index()
    agg["n_false_circular"] = agg["n_false_circular"] - agg["n_false_linear"]
    agg.loc[agg["mode"] == "best_of", ["n_merged", "n_missed_or_merged", "n_false_circular", "n_false_linear"]] = None
    return agg


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--jobs", type=int, default=16)
    ap.add_argument("--drop-paths", action="store_true", help="sensitivity: remove every linear walk before scoring")
    a = ap.parse_args()
    sfx = "_cycles_only" if a.drop_paths else ""
    todo = tasks(a.drop_paths)
    with Pool(a.jobs) as pool:
        results = pool.map(one, todo)
    errors = [e for _, e in results if e]
    for e in errors:
        print(e, file=sys.stderr)
    species = pd.concat([s for s, _ in results if s is not None], ignore_index=True)
    dup = species.duplicated(["case", "graph", "method", "mode", "truth"])
    if dup.any():  # a truth species in scope of two amplicon runs: keep its best row, and say so
        print("truth species in scope of two amplicon runs:", species.loc[dup, ["case", "graph", "method", "truth"]].drop_duplicates().to_string(), file=sys.stderr)
    front = ["case", "scenario", "graph", "method", "amplicon", "mode", "truth"]
    species = species[front + [c for c in species.columns if c not in front]]
    species.to_csv(os.path.join(HERE, "h2h_species%s.tsv" % sfx), sep="\t", index=False)
    table = summarise(species)
    table.to_csv(os.path.join(HERE, "h2h_table%s.tsv" % sfx), sep="\t", index=False, float_format="%.4f")
    scen = scenario_summary(table)
    scen.to_csv(os.path.join(HERE, "h2h_scenario%s.tsv" % sfx), sep="\t", index=False, float_format="%.4f")
    if not a.drop_paths:
        graph_rows().to_csv(os.path.join(HERE, "h2h_graph_metrics.tsv"), sep="\t", index=False, float_format="%.4f")
    print("{n} runs scored, {e} errors".format(n=len(todo) - len(errors), e=len(errors)))
    return 1 if errors else 0


if __name__ == "__main__":
    sys.exit(main())

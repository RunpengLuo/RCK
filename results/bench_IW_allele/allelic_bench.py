#!/usr/bin/env python3
"""bench_IW_allele item 1: run every method on the 12 allelic-sibling cases (bench_IW_allelic), score, collect.

bench_iw.py (results/bench_IW, only read) is imported and pointed at a new set, iw_allelic: the cases are the
accepted seeds of amplicon-allele-analysis configs/bench_IW_allelic/seeds.tsv (make_scenarios.py --set
iw_allelic; PLAN.md "Item 1 in detail"), in IW_ALLELIC_SCENARIOS and slot order. Methods, settings, graph
source (iw_intpsi's: rck_bulk's graph + the long-read path constraints for CoRAL / CE), retries and scoring
are bench_iw's, unchanged. Outputs: /n/fs/yubi-iw/runs/bench_IW_allelic/methods/<case>/<method>/; the
collected tables (bench_iw_allelic_{table,pooled,proportions}.tsv) here, never in results/bench_IW/.

    python allelic_bench.py list [--method M]        # index, case, method (IDX = case with --method)
    python allelic_bench.py run --method M IDX       # one run (Gurobi; the arrays run these)
    python allelic_bench.py status
    python allelic_bench.py score [--jobs N]         # offline, after the runs
    python allelic_bench.py collect

Arrays: submit_allelic.sh (prints the commands unless --submit). Gurobi WLS licence: 2 sessions, %1 per array.
"""
import argparse
import os
import sys

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_IW"))
import bench_iw as bi  # noqa: E402

SET, FOLDER, PREFIX = "iw_allelic", "bench_IW_allelic", "bench_iw_allelic"
bi.SETS[SET] = (FOLDER, PREFIX)
bi.GRAPH_SOURCES[SET] = "bulk+paths"  # as iw_intpsi: CoRAL / CE decompose rck_bulk's graph (+ path constraints)


def cases():
    """One dict per accepted seed, with bench_iw's fields (class hard, S 2) plus stratum and sibling scenario."""
    if not os.path.exists(bi.SEEDS):
        return []
    seeds = pd.read_csv(bi.SEEDS, sep="\t")
    seeds = seeds[seeds["accepted"].astype(str).str.lower() == "true"]
    sys.path.insert(0, bi.SIM)
    from make_scenarios import IW_ALLELIC_SCENARIOS

    order = {s: i for i, s in enumerate(IW_ALLELIC_SCENARIOS)}
    seeds = seeds.assign(o=seeds["scenario"].map(order)).sort_values(["o", "slot"])
    return [{"name": r.run, "scenario": r.scenario, "class": "hard", "S": 2, "replicate": int(r.slot),
             "seed": int(r.seed), "stratum": r.ext_stratum, "truth": os.path.join(bi.SIM_RUNS, r.run),
             "subgraph": os.path.join(bi.SUBGRAPHS, r.run, "union"),
             "rck_input": os.path.join(bi.SUBGRAPHS_CNP, r.run, "union", "rck_input")} for r in seeds.itertuples()]



# INPUT PATCH for CoRAL / CE (Yubi, 2026-11-21). CoRAL 3.1.1 and CycleExtractor cannot take a discordant
# edge with 0 reads: infer_discordant_edge_multiplicities takes log2 of the read counts, every partition
# score becomes NaN and it fails with IndexError before any solve. rck_bulk's graph can carry such edges
# (imputed adjacencies it used, from copy-number change points with no split reads): 4 of them on
# sib_opp_mirror_s17006, none on the other 11 cases. For a case with any, CoRAL / CE get a copy of the graph
# with those read counts set to ZERO_READ_PSEUDOCOUNT (copy numbers and everything else unchanged), in
# <case>/graph_bulk_reads1/; scoring uses the same copy. rck_bulk's own graph (graph_bulk/) is not changed,
# and neither are CoRAL / CE. The patched edges are listed in the copy's .json ("zero_read_patch") and so in
# each run.json (input_graph_counts.zero_read_patch).
ZERO_READ_PSEUDOCOUNT = 1.0
_bulk_inferred_graph = bi.inferred_graph


def zero_read_discordant(graph):
    """The discordant lines of an AA graph whose read count (4th column) is 0."""
    with open(graph) as handle:
        return [line.rstrip("\n") for line in handle
                if line.startswith("discordant") and float(line.split("\t")[3]) == 0]


def patched_graph(case):
    """bench_iw's inferred graph, or for a case with zero-read discordant edges the patched copy (above)."""
    import json
    import shutil

    graph = _bulk_inferred_graph(case)
    zeros = zero_read_discordant(graph)
    if not zeros:
        return graph
    folder = os.path.join(bi.OUT, case["name"], "graph_bulk_reads1")
    out = os.path.join(folder, os.path.basename(graph))
    if os.path.exists(out):
        return out
    os.makedirs(folder, exist_ok=True)
    with open(graph) as src, open(out + ".tmp", "w") as dst:
        for line in src:
            fields = line.rstrip("\n").split("\t")
            if fields[0] == "discordant" and float(fields[3]) == 0:
                fields[3] = str(ZERO_READ_PSEUDOCOUNT)
            dst.write("\t".join(fields) + "\n")
    with open(os.path.splitext(graph)[0] + ".json") as handle:
        meta = json.load(handle)
    meta["counts"]["zero_read_patch"] = {
        "note": "INPUT PATCH: CoRAL / CE cannot take 0-read discordant edges (log2(0) in "
                "infer_discordant_edge_multiplicities); read count set to {p:g} on these edges, nothing else "
                "changed".format(p=ZERO_READ_PSEUDOCOUNT),
        "source_graph": graph, "pseudocount": ZERO_READ_PSEUDOCOUNT, "edges": zeros}
    with open(os.path.splitext(out)[0] + ".json", "w") as handle:
        json.dump(meta, handle, indent=2)
    shutil.move(out + ".tmp", out)
    return out


def configure():
    """Point bench_iw at iw_allelic: paths, tables, cases() and the zero-read-patched CoRAL / CE graph (bench_iw's
    runs(), status(), run(), answer(), score() and collect() look these up at call time). Called on import, and
    again by allele_bench.use_set, since importing allele_bench configures bench_iw for iw_intpsi."""
    bi.configure(SET)
    bi.TABLE, bi.POOLED, bi.PROPORTIONS = (os.path.join(HERE, PREFIX + suffix) for suffix in
                                           ("_table.tsv", "_pooled.tsv", "_proportions.tsv"))
    bi.cases = cases
    bi.inferred_graph = patched_graph


configure()


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("list")
    p.add_argument("--method", choices=bi.METHODS, default=None)
    p = sub.add_parser("run")
    p.add_argument("index", type=int)
    p.add_argument("--method", choices=bi.METHODS, default=None)
    sub.add_parser("status")
    p = sub.add_parser("score")
    p.add_argument("--jobs", type=int, default=8)
    sub.add_parser("collect")
    a = parser.parse_args()
    if a.cmd == "list":
        for i, (case, method) in enumerate(bi.runs_of(a.method)):
            print(i, case["name"], method, case["stratum"], sep="\t")
    elif a.cmd == "run":
        return bi.run(a.index, a.method)
    elif a.cmd == "status":
        bi.status()
    elif a.cmd == "score":
        bi.score(a.jobs)
    else:
        bi.collect()
    return 0


if __name__ == "__main__":
    sys.exit(main())

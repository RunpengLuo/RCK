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
bi.configure(SET)
bi.TABLE, bi.POOLED, bi.PROPORTIONS = (os.path.join(HERE, PREFIX + suffix) for suffix in
                                       ("_table.tsv", "_pooled.tsv", "_proportions.tsv"))


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


bi.cases = cases  # bench_iw's runs(), runs_of(), status(), score() and collect() look it up at call time


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

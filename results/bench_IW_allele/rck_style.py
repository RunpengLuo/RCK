#!/usr/bin/env python3
"""bench_IW_allele Part A: RCK-style evaluation (rck.ecdna.rck_eval) of the answers allele_bench.py scored.

For every (case, method variant, phasing) that ``allele_bench.py run`` finished, scores its cycles (the
method's own, or rephased post hoc) in the terms Aganezov & Raphael (2020) used for RCK: haplotype-specific
novel-adjacency FNR / FDR (and FPR over the input candidates: the union subgraph's junctions), the Eq. 11 allele-specific
copy-number distance, and generalized infinite-sites violations; absolute and within phase blocks. The
truth of each case is scored too (``truth`` rows: its violations, a property of the simulation).

    python rck_style.py run [--jobs 16]     # finished tasks are skipped
    python rck_style.py collect             # rck_style_table.tsv here

Outputs: /n/fs/yubi-iw/runs/bench_IW_intpsi_allele/<case>/<method_variant>/<phasing>/rck_style/rck_eval.json
and <case>/truth/rck_style/rck_eval.json. Local, no Gurobi; reads allele_bench.py's outputs only.
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
import allele_bench as ab  # noqa: E402  (configures bench_iw for iw_intpsi)

bi = ab.bi
TABLE = os.path.join(HERE, "rck_style_table.tsv")


def tasks():
    """allele_bench's tasks, plus one truth task per case (method None)."""
    out = [(case, None, None, "truth") for case in bi.cases()]
    return out + ab.tasks()


def candidate_graph(case):
    """The union subgraph as an AA graph, as bench_iw.py exported it (bench_iw.subgraph_graph's path; never
    written here, so the bench_IW_intpsi runs stay read-only)."""
    path = os.path.join(bi.OUT, case["name"], "graph", "{c}_amplicon1_graph.txt".format(c=case["name"]))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def decomposed_graph(case):
    """The graph CoRAL / CE decompose (bench_iw.bulk_graph's path: RCK bulk's graph + path constraints); read
    only, never built here."""
    path = os.path.join(bi.OUT, case["name"], "graph_bulk", "{c}_amplicon1_graph.txt".format(c=case["name"]))
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return path


def out_dir(task):
    case, method, variant, phasing = task
    if method is None:
        return os.path.join(ab.OUT, case["name"], "truth", "rck_style")
    return os.path.join(ab.task_dir(case, method, variant, phasing), "rck_style")


def run_one(task):
    case, method, variant, phasing = task
    from rck.ecdna import aa_format as aaf, allele_metrics as am, bench_metrics as bm, evaluate as ev, rck_eval

    folder = out_dir(task)
    label = "truth" if method is None else "{m}:{p}".format(m=bi.method_variant(method, variant), p=phasing)
    if os.path.exists(os.path.join(folder, "rck_eval.json")):
        return case["name"], label, "done already"
    try:
        truth = ev.load_truth(case["truth"])
        if method is None:
            result = rck_eval.score_truth(truth)
        else:
            done = os.path.join(ab.task_dir(case, method, variant, phasing), "done.json")
            if not os.path.exists(done):
                return case["name"], label, "no allele_bench answer"
            with open(done) as handle:
                cycles = json.load(handle)["cycles"]
            weights = None
            if method == "rck_mc":
                weights = bm.multiclone_bulk_weights(bi.run_dir(case, method),
                                                     bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv")))
            strict, _ = am.strict_cycles(cycles, weights)
            blocks = ev.read_phase_blocks(os.path.join(case["rck_input"], "rck.fragments.tsv"))
            # the input candidates: the union subgraph's export (every junction the subgraph kept), read
            # only; rck_mc's DIAG is built from the same observations with the same junction filters
            graph = aaf.read_graph(candidate_graph(case))
            result = rck_eval.score(strict, truth, blocks, graph)
            result["cycles"] = cycles
            # RCK's Fig. 3B (input adjacencies used, FPR over spurious ones) and, for the multi-clone route,
            # Eq. 11 summed over clones with copy-number states (rck_eval docstring, 4-5)
            if method == "rck_mc":
                run = bi.run_dir(case, method)
                result["input"] = rck_eval.multiclone_input_use(run, case["rck_input"], truth)
                result["clone_cn"] = rck_eval.clone_cn_scores(
                    os.path.join(run, "rck.scnt.tsv"), os.path.join(case["truth"], "truth", "rck", "rck.scnt.tsv"),
                    os.path.join(case["rck_input"], "rck.fragments.tsv"), rck_eval.ecdna_regions(truth))
            elif method in ("coral", "ce"):
                result["input"] = rck_eval.graph_input_use(aaf.read_graph(decomposed_graph(case)), strict, truth)
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "rck_eval.json"), "w") as handle:
            json.dump(result, handle, indent=2)
        return case["name"], label, "ok"
    except Exception:  # noqa: BLE001 - recorded, the other tasks go on
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "error.txt"), "w") as handle:
            handle.write(traceback.format_exc())
        return case["name"], label, "error"


def run(jobs):
    import logging
    logging.disable(logging.INFO)
    states = {}
    with Pool(jobs) as pool:
        for name, label, state in pool.imap_unordered(run_one, tasks()):
            states[state] = states.get(state, 0) + 1
            if state not in ("ok", "done already"):
                print(name, label, state, flush=True)
    print(states)


def collect():
    from rck.ecdna import rck_eval
    rows = []
    for task in tasks():
        case, method, variant, phasing = task
        row = {k: case[k] for k in ("name", "class", "S", "replicate")}
        row.update(method=method or "truth", method_variant="truth" if method is None else bi.method_variant(method, variant),
                   phasing=phasing)
        path = os.path.join(out_dir(task), "rck_eval.json")
        if not os.path.exists(path):
            row["status"] = "error" if os.path.exists(os.path.join(out_dir(task), "error.txt")) else "missing"
        else:
            with open(path) as handle:
                result = json.load(handle)
            result.pop("cycles", None)
            row["status"] = "ok"
            result.get("clone_cn", {}).pop("clones", None)
            row.update(rck_eval.flatten(result))
        rows.append(row)
    frame = pd.DataFrame(rows)
    # mixed-allele cases (both alleles amplified on some ecDNA base), from items 2-3's table
    allele = pd.read_csv(ab.TABLE, sep="\t")
    mixed = allele[allele["phasing"] == "native"].groupby("name")["mixed_bp_frac"].first()
    frame["mixed_bp_frac"] = frame["name"].map(mixed)
    frame.to_csv(TABLE, sep="\t", index=False)
    print("{n} rows ({ok} ok) -> {t}".format(n=len(frame), ok=int((frame["status"] == "ok").sum()), t=TABLE))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=16)
    sub.add_parser("collect")
    a = p.parse_args()
    if a.cmd == "run":
        run(a.jobs)
    else:
        collect()


if __name__ == "__main__":
    main()

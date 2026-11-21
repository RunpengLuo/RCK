#!/usr/bin/env python3
"""bench_IW_allele: allele-specific analyses of the finished bench_IW_intpsi runs (PLAN.md items 2 and 3).

Reads the answers bench_iw.py scored (``<run>/score[_nolinear]/answer/``), never writes next to them.
For each case, method variant and phasing:

  phasing   native       the method's own haplotypes (rck_bulk, rck_mc; CoRAL / CE have none)
            baf_major    post-hoc: every step gets the major allele of the subgraph's observed CN_A / CN_B
                         (segment_blocks.tsv; depth + phased SNPs, before any decomposition)
            truth_major  post-hoc: the major allele of the true ecDNA copies (upper bound of any
                         one-haplotype-per-segment rule)
  scoring   aa_eval absolute and within phase blocks (--phase-blocks rck_input/rck.fragments.tsv), and
            rck.ecdna.allele_metrics (allele-specific ecDNA copy numbers, absolute and _block)

    python allele_bench.py run [--jobs 16]     # phase + score every (case, variant, phasing); finished ones skipped
    python allele_bench.py collect             # allele_table.tsv, allele_segments.tsv here
    python allele_bench.py --set iw_allelic run|collect   # item 1 cases (use_set); rck_mc, CoRAL, CE only

Outputs: /n/fs/yubi-iw/runs/bench_IW_intpsi_allele/<case>/<method_variant>/<phasing>/
  cycles.txt (+ alt<k>_cycles.txt), aa_eval/, aa_eval_block/, allele/ (allele_metrics.json, allele_segments.tsv)
No Gurobi. bench_iw.py is imported for its case list and paths only (configure("iw_intpsi")).
"""
import argparse
import json
import os
import sys
import traceback
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "bench_IW"))
import bench_iw as bi  # noqa: E402

bi.configure("iw_intpsi")

OUT = os.environ.get("BENCH_IW_ALLELE_OUT", "/n/fs/yubi-iw/runs/bench_IW_intpsi_allele")
TABLE = os.path.join(HERE, "allele_table.tsv")
SEGMENTS = os.path.join(HERE, "allele_segments.tsv")
VARIANTS = [("rck_bulk", ""), ("rck_mc", ""), ("coral", "walks_kept"), ("coral", "walks_dropped"),
            ("ce", "walks_kept"), ("ce", "walks_dropped")]
PHASINGS = ["native", "baf_major", "truth_major"]
SET = "iw_intpsi"


def use_set(name):
    """Point this module at a case set: iw_intpsi (default, items 2-3 and Part A) or iw_allelic (item 1: the
    allelic-sibling cases through allelic_bench.py, its outputs in /n/fs/yubi-iw/runs/bench_IW_allelic_allele/
    and *_allelic.tsv here; rck_mc, CoRAL and CE only, as Yubi asked on 2026-11-21). Call before run/collect."""
    global OUT, TABLE, SEGMENTS, VARIANTS, SET
    if name == "iw_intpsi":
        return
    if name != "iw_allelic":
        raise ValueError(name)
    import allelic_bench

    allelic_bench.configure()  # bench_iw -> iw_allelic (cases, zero-read-patched CoRAL / CE graph)

    SET = name
    OUT = os.environ.get("BENCH_IW_ALLELE_OUT", "/n/fs/yubi-iw/runs/bench_IW_allelic_allele")
    TABLE = os.path.join(HERE, "allele_table_allelic.tsv")
    SEGMENTS = os.path.join(HERE, "allele_segments_allelic.tsv")
    VARIANTS = [v for v in VARIANTS if v[0] != "rck_bulk"]


def answer_files(case, method, variant):
    """(main cycles file, [alternative cycles files]) of a scored run."""
    folder = os.path.join(bi.score_folder(case, method, variant), "answer")
    if method in ("rck_bulk", "rck_mc"):
        stem = os.path.join(folder, "{c}_{m}".format(c=case["name"], m=method))
        alts = [stem + "_alt{k}_cycles.txt".format(k=k) for k in (1, 2)]
        return stem + "_cycles.txt", [a for a in alts if os.path.exists(a)]
    return os.path.join(folder, "{m}_amplicon1_cycles.txt".format(m=method)), []


def tasks():
    out = []
    for case in bi.cases():
        for method, variant in VARIANTS:
            for phasing in PHASINGS:
                if phasing == "native" and method not in ("rck_bulk", "rck_mc"):
                    continue
                out.append((case, method, variant, phasing))
    return out


def task_dir(case, method, variant, phasing):
    return os.path.join(OUT, case["name"], bi.method_variant(method, variant), phasing)


def run_one(task):
    case, method, variant, phasing = task
    from rck.ecdna import aa_eval, allele_metrics as am, bench_metrics as bm, evaluate as ev

    folder = task_dir(case, method, variant, phasing)
    label = "{m}:{p}".format(m=bi.method_variant(method, variant), p=phasing)
    if os.path.exists(os.path.join(folder, "done.json")):
        return case["name"], label, "done already"
    try:
        os.makedirs(folder, exist_ok=True)
        main, alts = answer_files(case, method, variant)
        truth = ev.load_truth(case["truth"])
        if phasing == "native":
            cycles, alt_cycles, counts = main, alts, None
        else:
            rule = am.baf_major_rule(os.path.join(case["subgraph"], "segment_blocks.tsv")) if phasing == "baf_major" \
                else am.truth_major_rule(truth)
            cycles = os.path.join(folder, "cycles.txt")
            counts = am.rephase(main, rule, cycles)
            alt_cycles = []
            for k, a in enumerate(alts, start=1):
                alt_cycles.append(os.path.join(folder, "alt{k}_cycles.txt".format(k=k)))
                am.rephase(a, rule, alt_cycles[-1])
        fragments = os.path.join(case["rck_input"], "rck.fragments.tsv")
        for sub, blocks in (("aa_eval", None), ("aa_eval_block", fragments)):
            argv = ["--cycles", cycles, "--truth-dir", case["truth"], "--outdir", os.path.join(folder, sub),
                    "--case", "{c}:{l}".format(c=case["name"], l=label)]
            argv += ["--alt-cycles"] + alt_cycles if alt_cycles else []
            argv += ["--graph", bi.inferred_graph(case)] if method != "rck_mc" else []
            argv += ["--phase-blocks", blocks] if blocks else []
            aa_eval.main(argv)
        weights = None
        if method == "rck_mc":
            weights = bm.multiclone_bulk_weights(bi.run_dir(case, method),
                                                 bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv")))
        summary, table = am.score(cycles, truth, weights, ev.read_phase_blocks(fragments))
        os.makedirs(os.path.join(folder, "allele"), exist_ok=True)
        summary.update(case=case["name"], method=method, variant=variant, phasing=phasing, rephase_counts=counts)
        with open(os.path.join(folder, "allele", "allele_metrics.json"), "w") as handle:
            json.dump(summary, handle, indent=2)
        table.to_csv(os.path.join(folder, "allele", "allele_segments.tsv"), sep="\t", index=False)
        with open(os.path.join(folder, "done.json"), "w") as handle:
            json.dump({"cycles": cycles, "source": main}, handle)
        return case["name"], label, "ok"
    except Exception:  # noqa: BLE001 - recorded, the other tasks go on
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
            if state == "error":
                print(name, label, "error", flush=True)
    print(states)


def _aa_fields(ev_dir, suffix):
    row = {}
    sample = bi._read(os.path.join(ev_dir, "aa_sample_metrics.tsv"))
    if sample is not None:
        s = sample[sample["mode"] == "strict"].iloc[0]
        row.update({k + suffix: s.get(k) for k in ("n_inferred", "n_missed", "n_false", "mean_composite", "mean_composite_h0")})
    part = bi._read(os.path.join(ev_dir, "aa_partition_metrics.tsv"))
    if part is not None:
        head = part[part["headline"].astype(str).str.lower() == "true"]
        for r in head.to_dict("records"):
            row["{m}_{o}_{h}_f1{s}".format(m=r["metric"], o=r["mode"], h=r["haplotype"], s=suffix)] = r["f1"]
    matrix = bi._read(os.path.join(ev_dir, "aa_species_matrix.tsv"))
    for hap in ("blind", "phased"):
        row.update({k + suffix: v for k, v in bi.matrix_summary(matrix, hap).items()})
    return row


def collect():
    rows, segs = [], []
    for case, method, variant, phasing in tasks():
        folder = task_dir(case, method, variant, phasing)
        row = {k: case[k] for k in ("name", "class", "S", "replicate")}
        row.update(method=method, method_variant=bi.method_variant(method, variant), phasing=phasing)
        metrics = os.path.join(folder, "allele", "allele_metrics.json")
        if not os.path.exists(metrics):
            row["status"] = "error" if os.path.exists(os.path.join(folder, "error.txt")) else "missing"
            rows.append(row)
            continue
        row["status"] = "ok"
        with open(metrics) as handle:
            m = json.load(handle)
        counts = m.pop("rephase_counts", None) or {}
        row.update({k: v for k, v in m.items() if k not in ("case", "method", "variant", "phasing")})
        row["rephase_unphased_steps"] = counts.get(".", 0) if counts else None
        row["rephase_filled_steps"] = counts.get("filled", 0) if counts else None
        row.update(_aa_fields(os.path.join(folder, "aa_eval"), ""))
        row.update(_aa_fields(os.path.join(folder, "aa_eval_block"), "_block"))
        rows.append(row)
        s = pd.read_csv(os.path.join(folder, "allele", "allele_segments.tsv"), sep="\t")
        segs.append(s.assign(name=case["name"], **{"class": case["class"]}, S=case["S"],
                             method_variant=row["method_variant"], phasing=phasing))
    pd.DataFrame(rows).to_csv(TABLE, sep="\t", index=False)
    if segs:
        pd.concat(segs, ignore_index=True).to_csv(SEGMENTS, sep="\t", index=False)
    print("{n} rows -> {t}".format(n=len(rows), t=TABLE))


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    p.add_argument("--set", choices=["iw_intpsi", "iw_allelic"], default="iw_intpsi")
    r = sub.add_parser("run")
    r.add_argument("--jobs", type=int, default=16)
    sub.add_parser("collect")
    a = p.parse_args()
    use_set(a.set)
    if a.cmd == "run":
        run(a.jobs)
    else:
        collect()


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Summarise an RCK ecDNA benchmark: one row per (scenario, sample, amplicon), plus run errors.

Run from the RCK root after ``docs/examples/ecdna_bench.sbatch``:

    python docs/examples/collect_bench.py                       # results/bench, prints a table
    python docs/examples/collect_bench.py --out bench_summary.tsv
    python docs/examples/collect_bench.py --logs 'logs/ecbench_31948920_*.out'   # also tally log errors

Columns: scenario, sample, amplicon, call (ecDNA / linear / mixed / none / FAILED), rank (first
decomposition with every in-scope truth species in order), then per scored truth species
``truth:jaccard_hap:order_match:refit_copies(integer)/true_copies``. ``FAILED`` means RCK wrote no
summary.json (it crashed; see the log), ``EVAL_ERROR`` that scoring failed (message in note).
"""

import argparse
import collections
import glob
import json
import math
import os
import re


def _fmt(x, digits=2):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return "NA"
    if isinstance(x, float):
        return str(int(round(x))) if digits == 0 else str(round(x, digits))
    return str(x)


def collect(root):
    rows = []
    for run_dir in sorted(glob.glob(os.path.join(root, "*", "*", "*"))):
        if not os.path.isdir(run_dir):
            continue
        scenario, sample, amplicon = run_dir.rstrip("/").split(os.sep)[-3:]
        row = {"scenario": scenario, "sample": sample, "amplicon": amplicon, "call": "FAILED", "rank": "NA",
               "species": "", "note": ""}
        summary_path = os.path.join(run_dir, "summary.json")
        if os.path.exists(summary_path):
            summary = json.load(open(summary_path))
            row["call"] = summary.get("amplicon_call", "?")
            decs = summary.get("decompositions") or []
            if decs:
                row["note"] = "inferred " + ",".join("{t}:{w}".format(t=s["topology"][0], w=_fmt(s["weight"]))
                                                     for s in decs[0]["species"] if s.get("status") != "low_weight")
        eval_path = os.path.join(run_dir, "evaluation.json")
        if os.path.exists(eval_path):
            e = json.load(open(eval_path))
            if "error" in e:
                row["call"] = row["call"] + "/EVAL_ERROR"
                row["note"] = e["error"].splitlines()[0][:120]
            else:
                row["rank"] = _fmt(e.get("first_rank_with_all_orders_matched"))
                row["species"] = " ".join(
                    "{t}:{j}:{o}:{w}({wi})/{tw}".format(t=t["truth"], j=_fmt(t.get("jaccard_hap")), o="Y" if t.get("order_match") else "N",
                                                        w=_fmt(t.get("inferred_weight")), wi=_fmt(t.get("inferred_weight_int"), 0),
                                                        tw=_fmt(t.get("truth_weight")))
                    for t in e.get("top_decomposition", []) if t.get("inferred"))
        rows.append(row)
    return rows


def tally_logs(pattern):
    """Count distinct final error lines (the line after each traceback) across Slurm logs."""
    errors = collections.Counter()
    failed_logs = 0
    for path in sorted(glob.glob(pattern)):
        text = open(path, errors="replace").read()
        # final traceback lines, including module-qualified ones (pandas.errors.EmptyDataError,
        # gurobipy._exception.GurobiError), Slurm kills, and the script's own "RCK failed" lines
        found = re.findall(r"^(?:[\w.]*(?:Error|Exception)\b.*|.*(?:oom-kill|Killed|DUE TO TIME LIMIT).*)$", text, flags=re.M)
        found = [re.sub(r"/\S*/bench/[^/\s]+/", "<scenario>/", line) for line in found]
        found += ["RCK failed (run count)"] * len(re.findall(r"^RCK failed:", text, flags=re.M))
        if found:
            failed_logs += 1
        for line in found:
            errors[line[:160]] += 1
    return failed_logs, errors


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--root", default="results/bench", help="results root: <root>/<scenario>/<sample>/<amplicon>/")
    parser.add_argument("--scenario", default=None, help="only scenarios whose name starts with this, e.g. baseline")
    parser.add_argument("--out", default=None, help="also write the table as TSV")
    parser.add_argument("--logs", default=None, help="glob of Slurm logs to scan for errors, e.g. 'logs/ecbench_<jobid>_*.out'")
    args = parser.parse_args()

    rows = collect(args.root)
    if args.scenario:
        rows = [r for r in rows if r["scenario"].startswith(args.scenario)]
    cols = ["scenario", "sample", "amplicon", "call", "rank", "species", "note"]
    print("\t".join(cols))
    for r in rows:
        print("\t".join(str(r[c]) for c in cols))
    if args.out:
        with open(args.out, "w") as handle:
            handle.write("\t".join(cols) + "\n")
            for r in rows:
                handle.write("\t".join(str(r[c]) for c in cols) + "\n")

    calls = collections.Counter(r["call"] for r in rows)
    print("\n{n} runs in {s} scenarios; calls: {c}".format(n=len(rows), s=len({r['scenario'] for r in rows}),
                                                          c=", ".join("{k} {v}".format(k=k, v=v) for k, v in calls.most_common())))
    if args.logs:
        failed, errors = tally_logs(args.logs)
        print("\nlogs with errors: {f}".format(f=failed))
        for line, count in errors.most_common(15):
            print("{c:5d}  {l}".format(c=count, l=line))


if __name__ == "__main__":
    main()

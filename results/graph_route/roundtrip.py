"""Round-trip test of the graph route: subgraph -> AA graph -> subgraph -> RCK-ecDNA total-only.

    python roundtrip.py prepare [--arms original rt_comment rt_none_profile rt_none]   # offline, no Gurobi
    python roundtrip.py list                                                            # the LIST, one run per line
    python roundtrip.py run IDX                                                         # one RCK-ecDNA run (Slurm array task)
    python roundtrip.py compare                                                         # offline, after the array

Cases: demo5-demo8 and the 32 bench cases of results/bench_nested/run_case.py. Per case, under
$GRAPH_ROUTE_OUT/roundtrip/<case>/ (default /n/fs/yubi-iw/runs/graph_route):

  aa/            python -m rck.ecdna.aa_export subgraph (imputed candidates as #rck_imputed lines)
  rt_comment/    aa_import --imputed comment --cn-profile <original> --background-cn <original>:
                 the exact inverse; RCK-ecDNA must give the original's answer
  rt_none_profile/ aa_import --imputed none, original cn_profile and background: only the
                 imputed candidates are gone
  rt_none/       aa_import --imputed none, nothing but the graph (flank background from the
                 graph: 2): what CoRAL and CycleExtractor see
  rt_regenerate/ (optional arm) aa_import --imputed regenerate, original cn_profile and background
  runs/<arm>/    rck-ecdna on the original subgraph (arm "original") or on the arm's folder

Every run is total-only on both sides: --min-snps-allele 100000000 (no allele observation) and
--soft-junction-hap --soft-junction-penalty 0 (every junction labelling allowed at no cost,
which is what an imported junction without a haplotype label gets), --block-flip off (the
default). Plus --step-penalty 0, --no-plots, --threads 4, --time-limit 1800 and --truth-dir.

prepare also compares the graphs offline (rck.ecdna.aa_import.compare_subgraphs) and writes
prepare.tsv; compare reads the runs and writes compare.tsv (both here and in the output folder).
"""
import argparse
import json
import os
import sys
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RCK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, RCK)
sys.path.insert(0, os.path.join(RCK, "results", "bench_nested"))
OUT = os.path.join(os.environ.get("GRAPH_ROUTE_OUT", "/n/fs/yubi-iw/runs/graph_route"), "roundtrip")
LIST = os.path.join(HERE, "roundtrip.list")
ARMS = ["original", "rt_comment", "rt_none_profile", "rt_none"]
TOTAL_ONLY = ["--min-snps-allele", "100000000", "--soft-junction-hap", "--soft-junction-penalty", "0", "--step-penalty", "0"]
COMMON = ["--no-plots", "--threads", os.environ.get("SLURM_CPUS_PER_TASK", "4"), "--time-limit", "1800"]
RETRIES = 5


def cases():
    """(name, subgraph, truth) for demo5-demo8 and the bench cases, from run_case.py."""
    import run_case

    names = {"demo5": "demo5", "demo6_union": "demo6", "demo7_totalonly": "demo7", "demo8": "demo8"}
    out = []
    for c in run_case.cases():
        if c["kind"] == "regression" and c["name"] in names:
            out.append((names[c["name"]], c["subgraph"], c["truth"]))
        elif c["kind"] == "bench":
            out.append((c["name"], c["subgraph"], c["truth"]))
    return out


def import_args(arm, subgraph):
    with open(os.path.join(subgraph, "scale.json")) as handle:
        scale = json.load(handle)
    profile = ["--cn-profile", os.path.join(subgraph, "cn_profile.tsv"), "--background-cn", repr(float(scale["background_cn"]))]
    sample = ["--sample", str(scale.get("sample"))]
    return {"rt_comment": ["--imputed", "comment"] + profile, "rt_none_profile": ["--imputed", "none"] + profile,
            "rt_none": ["--imputed", "none"], "rt_regenerate": ["--imputed", "regenerate"] + profile}[arm] + sample


def prepare(arms):
    from rck.ecdna.aa_export import export_subgraph
    from rck.ecdna.aa_import import compare_subgraphs, main as import_main

    rows, lines = [], []
    for name, subgraph, truth in cases():
        base = os.path.join(OUT, name)
        info = export_subgraph(subgraph, os.path.join(base, "aa"), name=name)["info"]
        graph = os.path.join(base, "aa", "{n}_amplicon1_graph.txt".format(n=name))
        row = {"case": name, "foldbacks": ",".join(info["foldbacks"]) or ".", "imputed": info["counts"]["imputed"],
               "path_constraints": info["counts"]["path_constraint"], "rck_paths": info["counts"]["rck_path"]}
        for arm in arms:
            if arm == "original":
                lines.append("\t".join([name, arm, subgraph, os.path.join(base, "runs", arm), truth]))
                continue
            folder = os.path.join(base, arm)
            import_main([graph, folder] + import_args(arm, subgraph))
            diffs = compare_subgraphs(subgraph, folder, ignore_imputed=arm != "rt_comment")
            row[arm] = "exact" if not diffs else "; ".join(diffs)
            lines.append("\t".join([name, arm, folder, os.path.join(base, "runs", arm), truth]))
        rows.append(row)
        print(name, {k: v for k, v in row.items() if k.startswith("rt_")}, flush=True)
    with open(LIST, "w") as handle:
        handle.write("\n".join(lines) + "\n")
    table = pd.DataFrame(rows)
    for path in (os.path.join(HERE, "roundtrip_prepare.tsv"), os.path.join(OUT, "prepare.tsv")):
        table.to_csv(path, sep="\t", index=False)
    print("{n} runs in {f}".format(n=len(lines), f=LIST))


def read_list():
    with open(LIST) as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip()]


def done(outdir):
    return all(os.path.exists(os.path.join(outdir, f)) for f in ("summary.json", "decompositions.tsv"))


def run(index):
    from rck.ecdna_run import build_parser, run as rck_run

    entries = read_list()
    if index >= len(entries):
        print("no run {i} ({n} in the LIST)".format(i=index, n=len(entries)))
        return 0
    name, arm, subgraph, outdir, truth = entries[index]
    if done(outdir):
        print("done already: " + outdir)
        return 0
    argv = ["--subgraph-dir", subgraph, "--inputs", "subgraph", "--name", name, "--truth-dir", truth, "--outdir", outdir] + COMMON + TOTAL_ONLY
    for attempt in range(1, RETRIES + 1):
        started = time.time()
        try:
            rck_run(build_parser().parse_args(argv))
            print("{c} {a}: {t:.1f} s".format(c=name, a=arm, t=time.time() - started), flush=True)
            break
        except Exception as error:  # noqa: BLE001 - licence hiccups: retry, then report
            print("attempt {k} failed: {c} {a}: {e!r}".format(k=attempt, c=name, a=arm, e=error), flush=True)
            if attempt == RETRIES:
                return 1
            time.sleep(120)
    try:
        import gurobipy

        gurobipy.disposeDefaultEnv()
    except Exception:  # noqa: BLE001
        pass
    return 0


def answer(outdir, haps=True):
    """Rank-1 decomposition as a sorted tuple of (orderings, weight, status, topology)."""
    t = pd.read_csv(os.path.join(outdir, "decompositions.tsv"), sep="\t")
    t = t[t["decomposition_rank"] == t["decomposition_rank"].min()]
    out = []
    for sp, g in t.groupby("species", sort=False):
        walks = []
        for w in g.sort_values("ordering_rank")["walk"]:
            toks = str(w).split(",")
            if not haps:
                toks = [":".join(x.rsplit(":", 2)[0:1] + x.rsplit(":", 2)[2:]) for x in toks]
            walks.append(",".join(toks))
        r = g.iloc[0]
        out.append((tuple(sorted(walks)), round(float(r["weight"]), 6), str(r["status"]), str(r["topology"])))
    return tuple(sorted(out)), float(t["objective"].iloc[0]) if len(t) else None


def compare():
    entries = read_list()
    by_case = {}
    for name, arm, _s, outdir, _t in entries:
        by_case.setdefault(name, {})[arm] = outdir
    rows = []
    for name, arms in by_case.items():
        orig = arms.get("original")
        row = {"case": name}
        if not orig or not done(orig):
            row["status"] = "original missing"
            rows.append(row)
            continue
        with open(os.path.join(orig, "summary.json")) as handle:
            s = json.load(handle)
        row["imputed_used"] = ",".join((s.get("flow") or {}).get("imputed_used") or []) or "."
        flow_obj = (s.get("flow") or {}).get("objective")
        a, a_obj = answer(orig)
        ab, _ = answer(orig, haps=False)
        row["species"] = len(a)
        for arm, outdir in arms.items():
            if arm == "original":
                continue
            if not done(outdir):
                row[arm] = "missing"
                continue
            b, b_obj = answer(outdir)
            bb, _ = answer(outdir, haps=False)
            with open(os.path.join(outdir, "summary.json")) as handle:
                sb = json.load(handle)
            row[arm] = "identical" if a == b else ("same up to haplotype labels" if ab == bb else "differs")
            row[arm + "_dflow"] = None if flow_obj is None else (sb.get("flow") or {}).get("objective", float("nan")) - flow_obj
            row[arm + "_ddecomp"] = None if a_obj is None or b_obj is None else b_obj - a_obj
            row[arm + "_imputed_used"] = ",".join((sb.get("flow") or {}).get("imputed_used") or []) or "."
        rows.append(row)
    table = pd.DataFrame(rows)
    for path in (os.path.join(HERE, "roundtrip_compare.tsv"), os.path.join(OUT, "compare.tsv")):
        table.to_csv(path, sep="\t", index=False)
    print(table.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--arms", nargs="+", default=ARMS, choices=ARMS + ["rt_regenerate"])
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("index", type=int)
    sub.add_parser("compare")
    args = parser.parse_args()
    if args.cmd == "prepare":
        prepare(args.arms)
    elif args.cmd == "list":
        for i, e in enumerate(read_list()):
            print(i, *e, sep="\t")
    elif args.cmd == "run":
        return run(args.index)
    else:
        compare()
    return 0


if __name__ == "__main__":
    sys.exit(main())

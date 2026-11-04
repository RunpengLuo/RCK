"""Graph route (Task 2) head-to-head runs: every decomposer on the same AA-format graph.

    python graph_route.py prepare [--graphs truth_rp inferred]   # offline, no Gurobi: graphs, imports, LIST
    python graph_route.py list                                    # the LIST, one run per line
    python graph_route.py run IDX                                 # one run (Slurm array task)
    python graph_route.py status                                  # offline: done / missing / failed per run

Cases: demo5-demo8 and the 32 bench cases (roundtrip.cases()). Graph arms:

  truth_rp   the truth mixture graph with synthetic reads and path constraints
             ($G/<case>/truth/graph_reads_paths/mixture_amplicon<N>_graph.txt, one run per
             amplicon): the oracle headline
  inferred   our inferred union subgraph as an AA graph ($G/<case>/inferred/graph/, written here
             by aa_export subgraph; RCK-ecDNA's imputed candidates are #rck_imputed comment
             lines, which CoRAL and CycleExtractor skip)
  truth_plain (secondary arm, priority 2; not in the default LIST)

Methods, each with what it is given:

  coral      CoRAL 3.1.1 `coral cycle --solver gurobi_direct --solver-time-limit 1800 --threads T`
             (other options at CoRAL's defaults: --cycle-decomp-mode max_weight, alpha 0.01)
  ce         CycleExtractor 1.0.0 `CE.py --solver gurobi` (its defaults; its time limit is
             hard-coded, 7200 s per MILP, no flag); CE_SOLVER=highs runs it without the licence.
             The ce runs are the last block of the LIST, so they can be run apart
  rck_ecdna  rck-ecdna on `aa_import GRAPH --imputed none` (the graph only, as CoRAL and CE see
             it: flank background 2, no cn_profile, no allele data, so --block-flip does not
             apply), total-only (--min-snps-allele 100000000 --soft-junction-hap
             --soft-junction-penalty 0), the bench's step-penalty procedure: --step-penalty 0.1
             --keep-pool --pool-lambdas 0 0.02 0.05 0.1 0.2 0.5, answer re-ranked at λ = 0.2 when
             scoring (aa_export --lambda 0.2), as the bench's λ = 0.2 results are; --time-limit
             1800 per MILP, --threads T. Single seed (Gurobi's default).

Outputs under $GRAPH_ROUTE_OUT (default /n/fs/yubi-iw/runs/graph_route):
  <case>/inferred/graph/                     the exported inferred graph
  <case>/runs/<graph>/rck_input/amplicon<N>/ aa_import --imputed none of the graph
  <case>/runs/<graph>/<method>/amplicon<N>/  the run, with run.json (method, version, solver,
                                             time limit, input graph, command, elapsed, attempts)
A run is skipped when its outputs exist, so the array can be resubmitted. Gurobi "too many
sessions" errors (the academic WLS licence allows 2, and an exited process's session stays
counted for a few minutes) are retried up to 5 times, 2 min apart.
"""
import argparse
import json
import os
import shutil
import socket
import subprocess
import sys
import tarfile
import time

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RCK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, RCK)
sys.path.insert(0, HERE)
import roundtrip  # noqa: E402  (also puts results/bench_nested on sys.path)

G = os.environ.get("GRAPH_ROUTE_OUT", "/n/fs/yubi-iw/runs/graph_route")
LIST = os.path.join(HERE, "graph_route.list")
CORAL_ENV = "/n/fs/yubi-iw/conda/envs/coral-3.1.1"
CE_DIR = "/n/fs/yubi-iw/tools/CycleExtractor"
CORAL_VERSION = "3.1.1 (commit 7f3bf317)"
CE_VERSION = "1.0.0 (commit 1b29b33)"
GRAPHS = {"truth_rp": "graph_reads_paths", "truth_plain": "graph_plain"}
DEFAULT_GRAPHS = ["truth_rp", "inferred"]
METHODS = ["coral", "rck_ecdna", "ce"]  # LIST order: per graph coral, rck_ecdna; then every ce run as the last block
CE_SOLVER = os.environ.get("CE_SOLVER", "gurobi")  # gurobi (CE's default) or highs (no licence)
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "4")
TIME_LIMIT = 1800
POOL_LAMBDA = "0.1"
SCORE_LAMBDA = 0.2
LAMBDAS = ["0", "0.02", "0.05", "0.1", "0.2", "0.5"]
TOTAL_ONLY = ["--min-snps-allele", "100000000", "--soft-junction-hap", "--soft-junction-penalty", "0"]
RETRIES = 5
LICENCE_ERRORS = ("Too many sessions", "10030", "Error 10009", "License")


def graphs_of(case, subgraph, arm):
    """[(amplicon, graph file)] of one case and graph arm (inferred: exported here if missing)."""
    if arm == "inferred":
        from rck.ecdna.aa_export import export_subgraph

        outdir = os.path.join(G, case, "inferred", "graph")
        graph = os.path.join(outdir, "{c}_amplicon1_graph.txt".format(c=case))
        if not os.path.exists(graph):
            export_subgraph(subgraph, outdir, name=case)
        return [(1, graph)]
    folder = os.path.join(G, case, "truth", GRAPHS[arm])
    index = pd.read_csv(os.path.join(folder, "mixture_amplicons.tsv"), sep="\t")
    return [(int(r.amplicon), os.path.join(folder, r.graph)) for r in index.itertuples()]


def prepare(arms):
    from rck.ecdna.aa_import import main as import_main

    lines, rows = [], []
    for case, subgraph, truth in roundtrip.cases():
        with open(os.path.join(subgraph, "scale.json")) as handle:
            sample = str(json.load(handle).get("sample"))
        for arm in arms:
            for amp, graph in graphs_of(case, subgraph, arm):
                base = os.path.join(G, case, "runs", arm)
                rck_input = os.path.join(base, "rck_input", "amplicon{a}".format(a=amp))
                if not os.path.exists(os.path.join(rck_input, "import.json")):
                    import_main([graph, rck_input, "--imputed", "none", "--sample", sample])
                with open(os.path.join(rck_input, "import.json")) as handle:
                    info = json.load(handle)
                rows.append({"case": case, "graph": arm, "amplicon": amp, "input": graph,
                             "import": json.dumps({k: info[k] for k in ("counts", "scale") if k in info})})
                for method in METHODS:
                    outdir = os.path.join(base, method, "amplicon{a}".format(a=amp))
                    lines.append("\t".join([case, arm, str(amp), method, graph,
                                            rck_input if method == "rck_ecdna" else ".", outdir, truth]))
        print(case, flush=True)
    with open(LIST, "w") as handle:
        handle.write("\n".join([l for l in lines if l.split("\t")[3] != "ce"] + [l for l in lines if l.split("\t")[3] == "ce"]) + "\n")
    pd.DataFrame(rows).to_csv(os.path.join(HERE, "graph_route_prepare.tsv"), sep="\t", index=False)
    print("{n} runs in {f}".format(n=len(lines), f=LIST))


def read_list():
    with open(LIST) as handle:
        return [line.rstrip("\n").split("\t") for line in handle if line.strip()]


def cycles_file(method, outdir, amp):
    return os.path.join(outdir, "{m}_amplicon{a}_cycles.txt".format(m=method, a=amp))


def done(method, outdir, amp):
    if method == "rck_ecdna":
        need = ["summary.json", "decompositions.tsv", "pool/pool.json"]
        return all(os.path.exists(os.path.join(outdir, f)) for f in need)
    return os.path.exists(cycles_file(method, outdir, amp)) and os.path.exists(os.path.join(outdir, "run.json"))


def git_commit():
    try:
        head = subprocess.run(["git", "-C", RCK, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
        dirty = subprocess.run(["git", "-C", RCK, "status", "--porcelain", "--untracked-files=no"],
                               capture_output=True, text=True).stdout.strip()
        return head + ("-dirty" if dirty else "")
    except OSError:
        return "unknown"


def subprocess_run(argv, outdir, log_name):
    """Run a solver command; return (returncode, licence_error) with stdout+stderr in outdir/log_name."""
    log = os.path.join(outdir, log_name)
    with open(log, "w") as handle:
        rc = subprocess.run(argv, stdout=handle, stderr=subprocess.STDOUT, cwd=outdir).returncode
    with open(log, errors="replace") as handle:
        text = handle.read()
    return rc, any(e in text for e in LICENCE_ERRORS[:2])


def run_coral(graph, outdir, amp):
    argv = [os.path.join(CORAL_ENV, "bin", "coral"), "cycle", "--bp-graph", graph, "--output-prefix", os.path.join(outdir, "coral"),
            "--solver", "gurobi_direct", "--solver-time-limit", str(TIME_LIMIT), "--threads", THREADS]
    meta = {"version": CORAL_VERSION, "solver": "gurobi_direct", "time_limit_s": TIME_LIMIT,
            "time_limit_note": "--solver-time-limit per MILP; --global-time-limit default 21600", "command": argv}
    return argv, meta


def run_ce(graph, outdir, amp):
    argv = [os.path.join(CORAL_ENV, "bin", "python"), "-I", os.path.join(CE_DIR, "CE.py"), "--graph", graph,
            "--output", cycles_file("ce", outdir, amp), "--solver", CE_SOLVER, "--log_file", os.path.join(outdir, "ce.log")]
    meta = {"version": CE_VERSION, "solver": CE_SOLVER, "time_limit_s": 7200,
            "time_limit_note": "hard-coded TIME_LIMIT_SECONDS per MILP (cycle and path models), no flag",
            "command": argv}
    return argv, meta


def rck_argv(case, arm, amp, rck_input, outdir, truth):
    return (["--subgraph-dir", rck_input, "--inputs", "subgraph", "--name", "{c}_{g}_a{a}".format(c=case, g=arm, a=amp),
             "--truth-dir", truth, "--outdir", outdir, "--no-plots", "--threads", THREADS, "--time-limit", str(TIME_LIMIT)]
            + TOTAL_ONLY + ["--step-penalty", POOL_LAMBDA, "--keep-pool", "--pool-lambdas"] + LAMBDAS)


def gurobi_version(python):
    try:
        out = subprocess.run([python, "-c", "import gurobipy; print('.'.join(map(str, gurobipy.gurobi.version())))"],
                             capture_output=True, text=True, timeout=120)
        return out.stdout.strip() or None
    except (OSError, subprocess.TimeoutExpired):
        return None


def run(index):
    entries = read_list()
    if index >= len(entries):
        print("no run {i} ({n} in the LIST)".format(i=index, n=len(entries)))
        return 0
    case, arm, amp, method, graph, rck_input, outdir, truth = entries[index]
    amp = int(amp)
    if done(method, outdir, amp):
        print("done already: " + outdir)
        return 0
    os.makedirs(outdir, exist_ok=True)
    record = {"case": case, "graph": arm, "amplicon": amp, "method": method, "input_graph": graph,
              "host": socket.gethostname(), "slurm_job": "{j}_{t}".format(j=os.environ.get("SLURM_ARRAY_JOB_ID", ""),
                                                                         t=os.environ.get("SLURM_ARRAY_TASK_ID", "")),
              "threads": int(THREADS), "rck_commit": git_commit()}
    ok, attempts, started = False, 0, time.time()
    if method == "rck_ecdna":
        from rck.ecdna_run import build_parser, run as rck_run

        argv = rck_argv(case, arm, amp, rck_input, outdir, truth)
        import gurobipy

        record.update({"version": "rck-ecdna " + record["rck_commit"], "solver": "gurobi",
                       "gurobi": ".".join(map(str, gurobipy.gurobi.version())), "time_limit_s": TIME_LIMIT,
                       "time_limit_note": "--time-limit per MILP", "rck_input": rck_input,
                       "command": ["python", "-m", "rck.ecdna_run"] + argv,
                       "score_lambda": SCORE_LAMBDA, "seed": "Gurobi default (single seed)"})
        for attempts in range(1, RETRIES + 1):
            try:
                rck_run(build_parser().parse_args(argv))
                ok = True
                break
            except Exception as error:  # noqa: BLE001 - licence hiccups: retry, then report
                print("attempt {k} failed: {e!r}".format(k=attempts, e=error), flush=True)
                if attempts < RETRIES:
                    time.sleep(120)
        try:
            gurobipy.disposeDefaultEnv()
        except Exception:  # noqa: BLE001
            pass
    else:
        argv, meta = (run_coral if method == "coral" else run_ce)(graph, outdir, amp)
        record.update(meta)
        record["gurobi"] = None if meta["solver"] == "highs" else gurobi_version(os.path.join(CORAL_ENV, "bin", "python"))
        record["seed"] = "Gurobi default (single seed)"
        for attempts in range(1, RETRIES + 1):
            rc, licence = subprocess_run(argv, outdir, "{m}.stdout".format(m=method))
            ok = rc == 0 and os.path.exists(cycles_file(method, outdir, amp))
            if ok:
                break
            print("attempt {k} failed: rc {r}{l}".format(k=attempts, r=rc, l=" (licence)" if licence else ""), flush=True)
            if not licence or attempts == RETRIES:
                break
            time.sleep(120)
        models = os.path.join(outdir, "models")
        if method == "coral" and os.path.isdir(models):  # CoRAL's .lp / .nl files: one archive, not many files
            with tarfile.open(models + ".tar.gz", "w:gz") as tar:
                tar.add(models, arcname="models")
            shutil.rmtree(models)
    record.update({"status": "ok" if ok else "failed", "attempts": attempts, "elapsed_s": round(time.time() - started, 1)})
    with open(os.path.join(outdir, "run.json"), "w") as handle:
        json.dump(record, handle, indent=2)
    print("{c} {g} a{a} {m}: {s} in {t:.1f} s".format(c=case, g=arm, a=amp, m=method, s=record["status"],
                                                      t=record["elapsed_s"]), flush=True)
    return 0 if ok else 1


def status():
    rows = []
    for i, (case, arm, amp, method, _g, _r, outdir, _t) in enumerate(read_list()):
        state = "done" if done(method, outdir, int(amp)) else "missing"
        run_json = os.path.join(outdir, "run.json")
        elapsed = None
        if os.path.exists(run_json):
            with open(run_json) as handle:
                r = json.load(handle)
            elapsed = r.get("elapsed_s")
            if r.get("status") == "failed":
                state = "failed"
        rows.append({"idx": i, "case": case, "graph": arm, "amplicon": amp, "method": method, "state": state, "elapsed_s": elapsed})
    table = pd.DataFrame(rows)
    print(table.groupby(["method", "graph", "state"]).size().to_string())
    print("licence-bound time: {t:.0f} s".format(t=table["elapsed_s"].fillna(0).sum()))
    bad = table[table["state"] != "done"]
    if len(bad):
        print(bad.to_string(index=False))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    p = sub.add_parser("prepare")
    p.add_argument("--graphs", nargs="+", default=DEFAULT_GRAPHS, choices=list(GRAPHS) + ["inferred"])
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("index", type=int)
    sub.add_parser("status")
    args = parser.parse_args()
    if args.cmd == "prepare":
        prepare(args.graphs)
    elif args.cmd == "list":
        for i, e in enumerate(read_list()):
            print(i, *e, sep="\t")
    elif args.cmd == "run":
        return run(args.index)
    else:
        status()
    return 0


if __name__ == "__main__":
    sys.exit(main())

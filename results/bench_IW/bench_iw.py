"""bench_IW: run every method on every accepted bench_IW case, score each answer, collect the table.

    python bench_iw.py [--set iw_intpsi] CMD ...   # --set before the command; default iw (bench_IW)

    python bench_iw.py list                  # index, case, method (one Slurm array task per line)
    python bench_iw.py run IDX               # one (case, method) run (SLURM_ARRAY_TASK_ID)
    python bench_iw.py status                # done / missing / failed per run
    python bench_iw.py score [--jobs N]      # offline: export + aa_eval + bench_metrics per finished run
    python bench_iw.py collect               # offline: bench_iw_table.tsv and bench_iw_pooled.tsv here

Sets (--set): iw is bench_IW; iw_intpsi is bench_IW_intpsi, bench_IW re-simulated with whole-number
psi per clone (l3_mixture.integer_psi), so the multi-clone input and the scoring truth agree; it is the
benchmark from 2026-11-14 evening (NEXT_STEPS.md). A set chooses the case list, the simulator runs,
both subgraph folders, the method output folder and the results TSV names (SETS below).

Cases are the accepted seeds of amplicon-allele-analysis configs/<set folder>/seeds.tsv ($AAA, default
~/amplicon-allele-analysis-cnp), in make_scenarios.IW_SCENARIOS and slot order. Methods, each on the
case's inputs from amplicon-allele-analysis bench_IW.py:

  rck_bulk   rck-ecdna, bulk route, on subgraphs/bench_IW/<run>/union: --block-flip unoriented
             --step-penalty 0.1 --keep-pool --pool-lambdas 0 0.02 0.05 0.1 0.2 0.5, answer re-ranked at
             lambda 0.2 (rck.ecdna.rerank) as bench_nested's; --time-limit 1800 per MILP
  rck_mc     rck-ecdna --model multiclone on subgraphs/bench_IW_cnp/<run>/union/rck_input (the DIAG built
             with --clone-profiles truth --diag-labeling-groups block); --time-limit 60
  coral      CoRAL 3.1.1 on the case's graph (below), with 10 kb flanks (results/graph_route/graph_route.py
             run_coral: --solver gurobi_direct, 1800 s per MILP)
  ce         CycleExtractor 1.0.0 on the same graph (graph_route.py run_ce: its defaults, gurobi; Gurobi's
             threads capped at --threads by a gurobi.env in the run folder, CE has no thread flag)

The graph CoRAL and CE decompose (GRAPH_SOURCE, by set):
  iw         the union subgraph exported directly as an AA graph (aa_export subgraph): observed copy
             numbers, long-read path constraints. As array 32267060 ran.
  iw_intpsi  decomposition-only baselines (Yubi, 2026-11-14): the case's rck_bulk graph (aa/*_graph.txt,
             its fitted copy numbers and every edge it used) plus the long-read path constraints of the
             direct export, renumbered onto it (aa_export.baseline_graph; the reasons are in its
             docstring). A coral / ce run refuses to start until the case's rck_bulk run is done; the
             Slurm arrays also chain them (bench_iw_method.sbatch, submit_intpsi.sh).

Every run has --threads $SLURM_CPUS_PER_TASK (default 4) and keeps a run.json (command, elapsed, status).
Gurobi licence errors are retried up to LICENCE_RETRIES times (15), 2 min apart: the academic WLS
licence allows 2 sessions, and in array 32267060 two CoRAL runs used up 5 tries. Other errors of the
rck runs are retried up to RETRIES times (5), as before. A finished run is skipped, so an array can be
resubmitted. `list --method M` and `run --method M IDX` index the cases of one method (IDX = case),
for the one-array-per-method layout of bench_iw_method.sbatch.

Scoring (per run, into <run>/score/): the answer as an AA cycles file (rck_bulk: aa_export at the
lambda 0.2 pool index; rck_mc: aa_export rank 1; coral / ce: their cycles with Path= read as Cycle=,
results/graph_route/score_h2h.normalise_cycles), then
  aa_eval      species matrix (aa_species_matrix.tsv, J_seg / J_adj, blind and phased rows), M1/M2
               (aa_partition_metrics.tsv), strict and best-of species scores; graph = the exported
               inferred graph (coral, ce, rck_bulk); no truth mixture graph, so no graph-level metrics
  bench_metrics  ecDNA copy-number recall / precision (bulk copies; rck_mc weights are sum_k u_k psi_k
               with u from the DIAG) and ecDNA proportions, strict and with low_weight cycles; for
               rck_bulk and rck_mc also total copy number (background + ecDNA, total_cn_*, and
               total_cn_ecdna_* over the truth ecDNA bases; bench_metrics docstring)

Score variants (Yubi, 2026-11-14):
  rck_mc         phased scores within phase blocks (aa_eval --phase-blocks <DIAG>/rck.fragments.tsv; one
                 A/B swap per block, as the original RCK's per-fragment allele flip leaves it free), so its
                 n_missed / n_false, M1 and phased matrix are within-block; n_missed_absolute /
                 n_false_absolute keep the absolute pairing. rck_bulk is scored as before (absolute).
  coral, ce      scored twice: linear walks kept (score/, method_variant coral:walks_kept) and dropped
                 (score_nolinear/, coral:walks_dropped), reported side by side.

collect: one row per (case, method, variant) with the case fields (segment_max_traversals, species sizes,
footprint), strict_rule (Status or cn-frac:0.1, from bench_metrics), baseline (decomposition-only for
coral / ce in iw_intpsi) and every score.

Outputs under $BENCH_IW_OUT (default /n/fs/yubi-iw/runs/<set folder>/methods):
  <case>/graph/<case>_amplicon1_graph.txt   the exported union subgraph (coral and ce for iw; path
                                            constraints for iw_intpsi)
  <case>/graph_bulk/<case>_amplicon1_graph.txt  iw_intpsi: rck_bulk's graph + path constraints (coral, ce),
                                            with <case>_amplicon1_graph.json (counts, dropped edges)
  <case>/<method>/                          rck runs; coral / ce runs are in <method>/amplicon1/
  <case>/<method>/score/                    answer/, aa_eval/, bench_metrics.json, bench_proportions.tsv
"""
import argparse
import json
import os
import socket
import sys
import time
import traceback
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RCK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, RCK)
sys.path.insert(0, os.path.join(RCK, "results", "graph_route"))
AAA = os.environ.get("AAA", os.path.expanduser("~/amplicon-allele-analysis-cnp"))
SIM = os.path.join(AAA, "simulations")
# --set: folder name (simulator side and /n/fs/yubi-iw/runs/<folder>) and results TSV prefix
SETS = {"iw": ("bench_IW", "bench_iw"), "iw_intpsi": ("bench_IW_intpsi", "bench_iw_intpsi")}
GRAPH_SOURCES = {"iw": "subgraph", "iw_intpsi": "bulk+paths"}  # the graph coral / ce decompose (module docstring)
SET = GRAPH_SOURCE = SEEDS = SUBGRAPHS = SUBGRAPHS_CNP = SIM_RUNS = OUT = TABLE = POOLED = PROPORTIONS = None
METHODS = ["rck_bulk", "rck_mc", "coral", "ce"]
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "4")
TIME_LIMIT = 1800
MC_TIME_LIMIT = 60
POOL_LAMBDA, SCORE_LAMBDA = "0.1", 0.2
LAMBDAS = ["0", "0.02", "0.05", "0.1", "0.2", "0.5"]
RETRIES = 5
LICENCE_RETRIES = 15


def configure(set_name):
    """Point every path at one set (SETS); bench_IW ("iw") keeps its original paths and file names."""
    global SET, GRAPH_SOURCE, SEEDS, SUBGRAPHS, SUBGRAPHS_CNP, SIM_RUNS, OUT, TABLE, POOLED, PROPORTIONS
    folder, prefix = SETS[set_name]
    SET, GRAPH_SOURCE = set_name, GRAPH_SOURCES[set_name]
    SEEDS = os.path.join(SIM, "configs", folder, "seeds.tsv")
    SUBGRAPHS = os.path.join(SIM, "subgraphs", folder)
    SUBGRAPHS_CNP = os.path.join(SIM, "subgraphs", folder + "_cnp")
    SIM_RUNS = os.path.join(SIM, "runs", folder)
    OUT = os.environ.get("BENCH_IW_OUT", "/n/fs/yubi-iw/runs/{f}/methods".format(f=folder))
    TABLE = os.path.join(HERE, prefix + "_table.tsv")
    POOLED = os.path.join(HERE, prefix + "_pooled.tsv")
    PROPORTIONS = os.path.join(HERE, prefix + "_proportions.tsv")


configure("iw")


# ---------------------------------------------------------------------------------------------- #
# cases and runs
# ---------------------------------------------------------------------------------------------- #

def cases():
    """One dict per accepted bench_IW seed: name, scenario, class, S, replicate, sim run, subgraphs."""
    if not os.path.exists(SEEDS):
        return []
    seeds = pd.read_csv(SEEDS, sep="\t")
    seeds = seeds[seeds["accepted"].astype(str).str.lower() == "true"]
    sys.path.insert(0, SIM)
    from make_scenarios import IW_SCENARIOS

    order = {s: i for i, s in enumerate(IW_SCENARIOS)}
    seeds = seeds.assign(o=seeds["scenario"].map(order)).sort_values(["o", "slot"])
    out = []
    for r in seeds.itertuples():
        cls, count = r.scenario.split("_S")
        out.append({"name": r.run, "scenario": r.scenario, "class": cls, "S": int(count), "replicate": int(r.slot),
                    "seed": int(r.seed), "truth": os.path.join(SIM_RUNS, r.run),
                    "subgraph": os.path.join(SUBGRAPHS, r.run, "union"),
                    "rck_input": os.path.join(SUBGRAPHS_CNP, r.run, "union", "rck_input")})
    return out


def runs():
    """(case, method) in array order: every case's methods together."""
    return [(c, m) for c in cases() for m in METHODS]


def run_dir(case, method):
    base = os.path.join(OUT, case["name"], method)
    return os.path.join(base, "amplicon1") if method in ("coral", "ce") else base


def cycles_file(case, method):
    import graph_route as gr

    return gr.cycles_file(method, run_dir(case, method), 1)


def done(case, method):
    outdir = run_dir(case, method)
    if method == "rck_bulk":
        need = ["summary.json", "decompositions.tsv", "pool/pool.json"]
    elif method == "rck_mc":
        need = ["summary.json", "decompositions.tsv"]
    else:
        return os.path.exists(cycles_file(case, method)) and os.path.exists(os.path.join(outdir, "run.json"))
    return all(os.path.exists(os.path.join(outdir, f)) for f in need)


def subgraph_graph(case):
    """The union subgraph as an AA graph (aa_export subgraph), written once per case."""
    from rck.ecdna.aa_export import export_subgraph

    folder = os.path.join(OUT, case["name"], "graph")
    graph = os.path.join(folder, "{c}_amplicon1_graph.txt".format(c=case["name"]))
    if not os.path.exists(graph):
        export_subgraph(case["subgraph"], folder, name=case["name"])
    return graph


class BulkNotDone(RuntimeError):
    """A decomposition-only baseline asked for before the case's rck_bulk run finished."""


def bulk_graph(case):
    """rck_bulk's graph + the subgraph export's path constraints (aa_export.baseline_graph), written once
    per case. Written into a private folder and moved into place, so coral and ce of one case can
    ask for it at the same time."""
    import shutil
    import tempfile

    from rck.ecdna.aa_export import baseline_graph

    folder = os.path.join(OUT, case["name"], "graph_bulk")
    graph = os.path.join(folder, "{c}_amplicon1_graph.txt".format(c=case["name"]))
    if os.path.exists(graph):
        return graph
    if not done(case, "rck_bulk"):
        raise BulkNotDone("{c}: rck_bulk is not done ({d}); the baselines decompose its graph".format(
            c=case["name"], d=run_dir(case, "rck_bulk")))
    os.makedirs(folder, exist_ok=True)
    work = tempfile.mkdtemp(dir=folder, prefix=".build_")
    try:
        tmp = os.path.join(work, os.path.basename(graph))
        baseline_graph(run_dir(case, "rck_bulk"), subgraph_graph(case), tmp)
        os.replace(os.path.splitext(tmp)[0] + ".json", os.path.splitext(graph)[0] + ".json")
        os.replace(tmp, graph)
    finally:
        shutil.rmtree(work, ignore_errors=True)
    return graph


def inferred_graph(case):
    """The graph coral / ce decompose and aa_eval scores the graph-route methods and rck_bulk with
    (GRAPH_SOURCE): the direct export (iw) or rck_bulk's graph + path constraints (iw_intpsi)."""
    return subgraph_graph(case) if GRAPH_SOURCE == "subgraph" else bulk_graph(case)


def rck_argv(case, method, outdir):
    common = ["--truth-dir", case["truth"], "--outdir", outdir, "--no-plots", "--threads", THREADS]
    if method == "rck_bulk":
        return (["--subgraph-dir", case["subgraph"]] + common + ["--time-limit", str(TIME_LIMIT), "--block-flip", "unoriented",
                                                                 "--step-penalty", POOL_LAMBDA, "--keep-pool", "--pool-lambdas"] + LAMBDAS)
    return ["--model", "multiclone", "--rck-input", case["rck_input"]] + common + ["--time-limit", str(MC_TIME_LIMIT)]


def runs_of(method=None):
    """runs(), or with a method the cases in order (index = case) for that method's array."""
    return runs() if method is None else [(c, method) for c in cases()]


def _licence(text):
    import graph_route as gr

    return any(e in text for e in gr.LICENCE_ERRORS[:2])


def run(index, method=None):
    import graph_route as gr

    table = runs_of(method)
    if index >= len(table):
        print("no run {i} ({n} runs)".format(i=index, n=len(table)))
        return 0
    case, method = table[index]
    outdir = run_dir(case, method)
    if done(case, method):
        print("done already: " + outdir)
        return 0
    os.makedirs(outdir, exist_ok=True)
    record = {"case": case["name"], "method": method, "host": socket.gethostname(), "threads": int(THREADS),
              "slurm_job": "{j}_{t}".format(j=os.environ.get("SLURM_ARRAY_JOB_ID", ""), t=os.environ.get("SLURM_ARRAY_TASK_ID", "")),
              "rck_commit": gr.git_commit()}
    ok, attempts, started = False, 0, time.time()
    if method in ("rck_bulk", "rck_mc"):
        import gurobipy
        from rck.ecdna_run import build_parser, run as rck_run

        argv = rck_argv(case, method, outdir)
        record.update({"command": ["python", "-m", "rck.ecdna_run"] + argv, "solver": "gurobi",
                       "gurobi": ".".join(map(str, gurobipy.gurobi.version()))})
        attempts, other = 0, 0
        while True:
            attempts += 1
            try:
                rck_run(build_parser().parse_args(argv))
                ok = True
                break
            except Exception as error:  # noqa: BLE001 - licence hiccups: retry, then report
                licence = _licence(repr(error))
                other += 0 if licence else 1
                print("attempt {k} failed{l}: {e!r}".format(k=attempts, l=" (licence)" if licence else "", e=error), flush=True)
                if attempts >= LICENCE_RETRIES or other >= RETRIES:
                    break
                time.sleep(120)
        try:
            gurobipy.disposeDefaultEnv()
        except Exception:  # noqa: BLE001
            pass
    else:
        try:
            graph = inferred_graph(case)
        except BulkNotDone as error:  # no run.json: the run never started
            print(error, flush=True)
            return 2
        argv, meta = (gr.run_coral if method == "coral" else gr.run_ce)(graph, outdir, 1)
        record.update(meta)
        record.update(input_graph=graph, graph_source=GRAPH_SOURCE)
        if GRAPH_SOURCE != "subgraph":
            record["baseline"] = "decomposition-only"
            with open(os.path.splitext(graph)[0] + ".json") as handle:
                record["input_graph_counts"] = json.load(handle)["counts"]
        if method == "ce":  # CE has no thread flag; Gurobi reads gurobi.env in the working folder (outdir)
            with open(os.path.join(outdir, "gurobi.env"), "w") as handle:
                handle.write("Threads {t}\n".format(t=THREADS))
            record["gurobi_env"] = "Threads {t}".format(t=THREADS)
        for attempts in range(1, LICENCE_RETRIES + 1):
            rc, licence = gr.subprocess_run(argv, outdir, "{m}.stdout".format(m=method))
            ok = rc == 0 and os.path.exists(cycles_file(case, method))
            if ok:
                break
            print("attempt {k} failed: rc {r}{l}".format(k=attempts, r=rc, l=" (licence)" if licence else ""), flush=True)
            if not licence or attempts == LICENCE_RETRIES:
                break
            time.sleep(120)
    record.update({"status": "ok" if ok else "failed", "attempts": attempts, "elapsed_s": round(time.time() - started, 1)})
    with open(os.path.join(outdir, "run.json"), "w") as handle:
        json.dump(record, handle, indent=2)
    print("{c} {m}: {s} in {t:.1f} s".format(c=case["name"], m=method, s=record["status"], t=record["elapsed_s"]), flush=True)
    return 0 if ok else 1


def elapsed(case, method):
    path = os.path.join(run_dir(case, method), "run.json")
    if not os.path.exists(path):
        return None, None
    with open(path) as handle:
        record = json.load(handle)
    return record.get("status"), record.get("elapsed_s")


def status():
    rows = []
    for i, (case, method) in enumerate(runs()):
        state, seconds = elapsed(case, method)
        rows.append({"idx": i, "case": case["name"], "method": method,
                     "state": "done" if done(case, method) else (state or "missing"), "elapsed_s": seconds})
    table = pd.DataFrame(rows)
    if table.empty:
        print("no cases (is {s} there?)".format(s=SEEDS))
        return
    print(table.groupby(["method", "state"]).size().to_string())
    bad = table[table["state"] != "done"]
    if len(bad):
        print(bad.to_string(index=False))


# ---------------------------------------------------------------------------------------------- #
# scoring
# ---------------------------------------------------------------------------------------------- #

# Score variants per method (Yubi, 2026-11-14): CoRAL / CE answers scored twice, side by side,
# with their linear walks (Path= / 0+ ... 0-) kept as candidate species, and with them dropped
# (score_h2h.normalise_cycles drop_paths). Neither removes background closed into cycles.
VARIANTS = {"rck_bulk": [""], "rck_mc": [""], "coral": ["walks_kept", "walks_dropped"], "ce": ["walks_kept", "walks_dropped"]}


def score_folder(case, method, variant=""):
    return os.path.join(run_dir(case, method), "score_nolinear" if variant == "walks_dropped" else "score")


def method_variant(method, variant=""):
    return method + (":" + variant if variant else "")


def answer(case, method, folder, drop_paths=False):
    """{cycles, alt_cycles, graph?} of one finished run, as AA files in folder."""
    from rck.ecdna import aa_export
    from rck.ecdna.rerank import rerank
    import score_h2h

    outdir = run_dir(case, method)
    if method == "rck_bulk":
        index = int(rerank(outdir, [SCORE_LAMBDA])["pool_index"].iloc[0])
        paths = aa_export.export(outdir, folder, index, alternatives=2, name=case["name"] + "_" + method)
        paths["graph"] = inferred_graph(case)
        return paths
    if method == "rck_mc":
        return aa_export.export(outdir, folder, None, alternatives=0, name=case["name"] + "_" + method)
    return {"cycles": score_h2h.normalise_cycles(cycles_file(case, method), folder, drop_paths=drop_paths), "alt_cycles": [],
            "graph": inferred_graph(case)}


def score_one(task):
    case, method, variant = task
    folder = score_folder(case, method, variant)
    label = method_variant(method, variant)
    try:
        from rck.ecdna import aa_eval, bench_metrics as bm
        from rck.ecdna import evaluate as ev

        paths = answer(case, method, os.path.join(folder, "answer"), drop_paths=variant == "walks_dropped")
        argv = ["--cycles", paths["cycles"], "--truth-dir", case["truth"], "--outdir", os.path.join(folder, "aa_eval"),
                "--case", "{c}:{m}".format(c=case["name"], m=label)]
        argv += ["--alt-cycles"] + paths["alt_cycles"] if paths.get("alt_cycles") else []
        argv += ["--graph", paths["graph"]] if paths.get("graph") else []
        fragments = os.path.join(case["rck_input"], "rck.fragments.tsv")
        if method == "rck_mc" and os.path.exists(fragments):  # phase within phase blocks (aa_eval.score docstring)
            argv += ["--phase-blocks", fragments]
        aa_eval.main(argv)
        weights = props = None
        if method == "rck_mc":
            props = bm.read_proportions(os.path.join(case["rck_input"], "clone_proportions.tsv"))
            weights = bm.multiclone_bulk_weights(run_dir(case, method), props)
        truth = ev.load_truth(case["truth"])
        summary, rows = bm.score(paths["cycles"], truth, weights)
        if method == "rck_bulk":  # total copy number (background + ecDNA): the two RCK routes only
            summary.update(bm.total_cn_all(bm.truth_total_cn(case["truth"]), bm.bulk_total_cn(run_dir(case, method)), truth))
        elif method == "rck_mc":
            summary.update(bm.total_cn_all(bm.truth_total_cn(case["truth"]),
                                           bm.multiclone_total_cn(run_dir(case, method), props), truth))
        summary.update(case=case["name"], method=method, variant=variant)
        with open(os.path.join(folder, "bench_metrics.json"), "w") as handle:
            json.dump(summary, handle, indent=2)
        pd.DataFrame(rows).assign(case=case["name"], method=method, method_variant=label).to_csv(
            os.path.join(folder, "bench_proportions.tsv"), sep="\t", index=False)
        return case["name"], label, "ok"
    except Exception:  # noqa: BLE001 - recorded, the other runs go on
        os.makedirs(folder, exist_ok=True)
        with open(os.path.join(folder, "score_error.txt"), "w") as handle:
            handle.write(traceback.format_exc())
        return case["name"], label, "error"


def score(jobs):
    tasks = [(c, m, v) for c, m in runs() if done(c, m) for v in VARIANTS[m]]
    with Pool(jobs) as pool:
        for name, method, state in pool.imap_unordered(score_one, tasks):
            print(name, method, state, flush=True)


# ---------------------------------------------------------------------------------------------- #
# table
# ---------------------------------------------------------------------------------------------- #

def _read(path):
    return pd.read_csv(path, sep="\t") if os.path.exists(path) else None


def matrix_summary(matrix, haplotype):
    """Mean over truth species of the best strict J_seg / J_adj, and exact cycles found."""
    if matrix is None or matrix.empty:
        return {}
    m = matrix[(matrix["haplotype"] == haplotype) & (~matrix["low_weight"].astype(bool))]
    if m.empty:
        return {}
    best = m.groupby("truth").agg(j_seg=("j_seg", "max"), j_adj=("j_adj", "max"), exact=("exact_cycle", "max"))
    return {"{h}_jseg_best_mean".format(h=haplotype): best["j_seg"].mean(),
            "{h}_jadj_best_mean".format(h=haplotype): best["j_adj"].mean(),
            "{h}_exact_cycles".format(h=haplotype): int(best["exact"].fillna(0).astype(bool).sum())}


def case_fields(case):
    """Per-case fields Yubi asked every row to carry: segment_max_traversals (seeds.tsv, the seed
    check's count over the truth walks) and the species sizes (truth/ecdna.tsv SIZE, bp)."""
    seeds = pd.read_csv(SEEDS, sep="\t")
    row = seeds[seeds["run"] == case["name"]].iloc[0]
    ecdna = pd.read_csv(os.path.join(case["truth"], "truth", "ecdna.tsv"), sep="\t")
    sizes = [int(x) for x in ecdna["SIZE"]]
    return {"segment_max_traversals": int(row["segment_max_traversals"]),
            "species_size_bp": ",".join(str(x) for x in sizes), "species_size_min_bp": min(sizes),
            "species_size_max_bp": max(sizes), "footprint_bp": int(row["footprint_bp"])}


def collect():
    rows, props = [], []
    fields = {c["name"]: case_fields(c) for c in cases()}
    for case, method, variant in [(c, m, v) for c, m in runs() for v in VARIANTS[m]]:
        state, seconds = elapsed(case, method)
        row = {k: case[k] for k in ("name", "scenario", "class", "S", "replicate", "seed")}
        row.update(fields[case["name"]])
        row["baseline"] = "decomposition-only" if method in ("coral", "ce") and GRAPH_SOURCE != "subgraph" else ""
        row.update(method=method, method_variant=method_variant(method, variant),
                   linear_walks={"walks_kept": "kept", "walks_dropped": "dropped"}.get(variant, ""),
                   run_status=state, elapsed_s=seconds)
        folder = score_folder(case, method, variant)
        metrics = os.path.join(folder, "bench_metrics.json")
        if os.path.exists(metrics):
            with open(metrics) as handle:
                m = json.load(handle)
            row.update({k: v for k, v in m.items() if k not in ("case", "method", "variant")})
            ev_dir = os.path.join(folder, "aa_eval")
            sample = _read(os.path.join(ev_dir, "aa_sample_metrics.tsv"))
            if sample is not None:
                strict = sample[sample["mode"] == "strict"].iloc[0]
                row.update({"n_inferred_strict": strict["n_inferred"], "n_missed": strict["n_missed"],
                            "n_false": strict["n_false"], "mean_composite_h0": strict["mean_composite_h0"],
                            "mean_composite": strict.get("mean_composite"), "phase_frame": strict.get("phase_frame"),
                            "n_missed_absolute": strict.get("n_missed_absolute"), "n_false_absolute": strict.get("n_false_absolute")})
            partition = _read(os.path.join(ev_dir, "aa_partition_metrics.tsv"))
            if partition is not None:  # headline rows: M1 species_detection, M2 adjacency_pairs
                head = partition[partition["headline"].astype(str).str.lower() == "true"]
                for r in head.to_dict("records"):
                    tag = "{m}_{o}_{h}".format(m=r["metric"], o=r["mode"], h=r["haplotype"])
                    for key in ("precision", "recall", "f1"):
                        row["{t}_{k}".format(t=tag, k=key)] = r[key]
            matrix = _read(os.path.join(ev_dir, "aa_species_matrix.tsv"))
            row.update(matrix_summary(matrix, "blind"))
            row.update(matrix_summary(matrix, "phased"))
            p = _read(os.path.join(folder, "bench_proportions.tsv"))
            if p is not None:
                props.append(p)
        elif os.path.exists(os.path.join(folder, "score_error.txt")):
            row["score_status"] = "error"
        row.setdefault("strict_rule", None)
        rows.append(row)
    table = pd.DataFrame(rows)
    table.to_csv(TABLE, sep="\t", index=False)
    if props:
        pd.concat(props, ignore_index=True).to_csv(PROPORTIONS, sep="\t", index=False)
    case_cols = ("replicate", "seed", "S", "segment_max_traversals", "species_size_min_bp", "species_size_max_bp", "footprint_bp")
    numeric = [c for c in table.columns if c not in case_cols and pd.api.types.is_numeric_dtype(table[c])]
    pooled = []
    for keys in (["method_variant", "class", "S"], ["method_variant", "class"], ["method_variant"]):
        g = table.groupby(keys)[numeric].agg(["mean", "median"])
        g.columns = ["{c}_{s}".format(c=c, s=s) for c, s in g.columns]
        scored = table["cn_recall"].notna() if "cn_recall" in table else pd.Series(False, index=table.index)
        g = g.reset_index().assign(n=table.groupby(keys).size().values, n_scored=scored.groupby([table[k] for k in keys]).sum().values,
                                   level="+".join(keys))
        pooled.append(g)
    pd.concat(pooled, ignore_index=True).to_csv(POOLED, sep="\t", index=False)
    print("{n} rows -> {t}; pooled -> {p}".format(n=len(table), t=TABLE, p=POOLED))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--set", choices=sorted(SETS), default="iw",
                        help="iw: bench_IW (default); iw_intpsi: bench_IW_intpsi (whole-number psi)")
    sub = parser.add_subparsers(dest="cmd", required=True)
    lst = sub.add_parser("list")
    lst.add_argument("--method", choices=METHODS, default=None, help="only this method; the index is then the case")
    r = sub.add_parser("run")
    r.add_argument("index", type=int)
    r.add_argument("--method", choices=METHODS, default=None, help="IDX is the case index of this method's array")
    sub.add_parser("status")
    s = sub.add_parser("score")
    s.add_argument("--jobs", type=int, default=8)
    sub.add_parser("collect")
    args = parser.parse_args()
    configure(args.set)
    if args.cmd == "list":
        for i, (case, method) in enumerate(runs_of(args.method)):
            print(i, case["name"], method, sep="\t")
    elif args.cmd == "run":
        return run(args.index, args.method)
    elif args.cmd == "status":
        status()
    elif args.cmd == "score":
        score(args.jobs)
    else:
        collect()
    return 0


if __name__ == "__main__":
    sys.exit(main())

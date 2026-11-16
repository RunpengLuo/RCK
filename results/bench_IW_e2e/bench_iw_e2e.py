"""bench_IW_e2e: CoRAL and CycleExtractor end to end (reads -> their own graph -> cycles) on bench_IW_intpsi.

    python bench_iw_e2e.py list [--method M]          # index, case (one Slurm array task per line)
    python bench_iw_e2e.py reads IDX                  # simulated long reads, BAMs, CNVkit for case IDX
    python bench_iw_e2e.py run --method M IDX         # coral_e2e or ce_e2e for case IDX
    python bench_iw_e2e.py status
    python bench_iw_e2e.py mixture [--jobs N]         # truth mixture graphs for graph-level metrics
    python bench_iw_e2e.py inputs IDX                 # RCK-ecDNA's inputs from the BAMs (obs/, subgraph, DIAG)
    python bench_iw_e2e.py score [--jobs N] [--case IDX]   # offline: aa_eval (+ mixture graph) + bench_metrics
    python bench_iw_e2e.py collect                    # bench_iw_e2e_table.tsv, _pooled.tsv, _proportions.tsv here

Why a separate folder: results/bench_IW/ holds the finished bench_IW_intpsi comparison (decomposition-only
baselines on RCK bulk's graph, RESULTS_intpsi.md, figures). Nothing here writes there or into the bench_IW_intpsi
run folders. This driver imports results/bench_IW/bench_iw.py for the case list, the scorer and the table code,
and results/graph_route/graph_route.py for the CoRAL / CE settings; it does not modify either.

Cases: the 24 accepted bench_IW_intpsi seeds, same truth (the simulator runs are only read).

Data (E2E, default /n/fs/yubi-iw/runs/bench_IW_intpsi_e2e, env BENCH_IW_E2E):
  reads/<case>/      amplicon-allele-analysis simulations/e2e/emit_reads.py (footprint.bed, groups.tsv,
                     manifest.json) and e2e/align_cnv.sh (tumor.bam, normal.bam, cnvkit/tumor.cns)
  mixture/<case>/    truth mixture graph (simulations/truth_mixture.py --outdir here)
  methods/<case>/coral_e2e/   coral seed (seed_CNV_SEEDS.bed) + coral reconstruct (coral_amplicon<N>_{graph,cycles}.txt)
  methods/<case>/ce_e2e/amplicon<N>/   CE on CoRAL's amplicon<N> graph
  methods/<case>/<method>/score[_nolinear]/   answer/ (amplicons pooled), aa_eval/, bench_metrics.json

Methods:
  coral_e2e  coral seed (defaults: --gain 6, --min-seed-size 100000, --max-seg-gap 300000; GRCh38 centromeres)
             then coral reconstruct (defaults, but --solver-time-limit 1800 and --solver-threads $THREADS as
             graph_route.run_coral; --global-time-limit default 21600)
  ce_e2e     CycleExtractor on every coral_e2e amplicon graph (graph_route.run_ce: defaults, gurobi; gurobi.env
             Threads cap); waits for coral_e2e of the case
Scoring: CoRAL / CE report one cycles file per amplicon; they are pooled into one answer per case
(pool_amplicons: segments renumbered, every walk kept), then scored exactly as bench_iw.score_one, with linear
walks kept and dropped, plus --mixture-dir for the graph-level metrics, against the pooled graph.
"""
import argparse
import glob
import json
import os
import re
import shutil
import socket
import subprocess
import sys
import time
import traceback
from multiprocessing import Pool

import pandas as pd

HERE = os.path.dirname(os.path.abspath(__file__))
RCK = os.path.abspath(os.path.join(HERE, "..", ".."))
sys.path.insert(0, os.path.join(RCK, "results", "bench_IW"))
sys.path.insert(0, os.path.join(RCK, "results", "graph_route"))
import bench_iw as bi  # noqa: E402
import graph_route as gr  # noqa: E402

bi.configure("iw_intpsi")
E2E = os.environ.get("BENCH_IW_E2E", "/n/fs/yubi-iw/runs/bench_IW_intpsi_e2e")
READS, MIXTURE, OUT, OBS_RUNS, SUBGRAPHS, SUBGRAPHS_CNP = (
    os.path.join(E2E, d) for d in ("reads", "mixture", "methods", "obs_runs", "subgraphs", "subgraphs_cnp"))
METHODS = ["coral_e2e", "ce_e2e", "rck_bulk_e2e", "rck_mc_e2e"]
RCK_METHODS = {"rck_bulk_e2e": "rck_bulk", "rck_mc_e2e": "rck_mc"}
ALIGN_PYTHON = os.environ.get("ALIGN_PYTHON", "/n/fs/yubi-iw/conda/envs/e2e-align/bin/python")
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "8")
CORAL = os.path.join(gr.CORAL_ENV, "bin", "coral")
CENTROMERES = os.path.join(bi.AAA, "third_party", "ecSimulator", "low_comp_regions", "GRCh38_centromere.bed")
REFERENCE = os.environ.get("E2E_REF", "/n/fs/yubi-iw/ref/GRCh38_no_alt.fa")
SIM_PYTHON = os.environ.get("SIM_PYTHON", "/n/fs/yubi-iw/conda/envs/gurobi-env/bin/python")
FLANK_BP = 5_000_000
TABLE, POOLED, PROPORTIONS = (os.path.join(HERE, "bench_iw_e2e_" + n + ".tsv") for n in ("table", "pooled", "proportions"))
VARIANTS = {"coral_e2e": ["walks_kept", "walks_dropped"], "ce_e2e": ["walks_kept", "walks_dropped"],
            "rck_bulk_e2e": [""], "rck_mc_e2e": [""]}


def cases():
    return bi.cases()


def e2e_case(case):
    """The case with RCK-ecDNA's inputs pointed at the BAM-derived subgraph and DIAG (bench_iw.rck_argv reads these)."""
    return dict(case, subgraph=os.path.join(SUBGRAPHS, case["name"], "union"),
                rck_input=os.path.join(SUBGRAPHS_CNP, case["name"], "union", "rck_input"))


def inputs_done(case):
    c = e2e_case(case)
    return all(os.path.exists(p) for p in (os.path.join(c["subgraph"], "scale.json"), os.path.join(c["rck_input"], "rck.scnt.tsv"),
                                           os.path.join(SUBGRAPHS, case["name"], "build_e2e.json")))


def reads_dir(case):
    return os.path.join(READS, case["name"])


def reads_done(case):
    return os.path.exists(os.path.join(reads_dir(case), "align_cnv.json"))


def run_dir(case, method):
    return os.path.join(OUT, case["name"], method)


def coral_amplicons(case):
    """{N: (graph, cycles)} of a finished coral_e2e run."""
    out = {}
    for graph in glob.glob(os.path.join(run_dir(case, "coral_e2e"), "coral_amplicon*_graph.txt")):
        n = int(re.search(r"amplicon(\d+)_graph\.txt$", graph).group(1))
        out[n] = (graph, graph[:-len("_graph.txt")] + "_cycles.txt")
    return dict(sorted(out.items()))


def answer_files(case, method):
    """Per-amplicon cycles files of a finished run (CE writes graph_route.cycles_file per amplicon)."""
    if method == "coral_e2e":
        return [c for g, c in coral_amplicons(case).values() if os.path.exists(c)]
    return [gr.cycles_file("ce", os.path.join(run_dir(case, method), "amplicon{n}".format(n=n)), n)
            for n in coral_amplicons(case)]


def done(case, method):
    if method in RCK_METHODS:
        need = ["summary.json", "decompositions.tsv"] + (["pool/pool.json"] if method == "rck_bulk_e2e" else [])
        return all(os.path.exists(os.path.join(run_dir(case, method), f)) for f in need)
    record = os.path.join(run_dir(case, method), "run.json")
    if not os.path.exists(record):
        return False
    with open(record) as handle:
        return json.load(handle).get("status") == "ok"


# ---------------------------------------------------------------------------------------------- #
# reads
# ---------------------------------------------------------------------------------------------- #

def make_reads(index):
    table = cases()
    if index >= len(table):
        print("no case {i}".format(i=index))
        return 0
    case = table[index]
    out = reads_dir(case)
    if reads_done(case):
        print("done already: " + out)
        return 0
    config = os.path.join(bi.SIM, "configs", "bench_IW_intpsi", case["name"] + ".yaml")
    if not os.path.exists(os.path.join(out, "manifest.json")):
        subprocess.run([SIM_PYTHON, os.path.join(bi.SIM, "e2e", "emit_reads.py"), "--run", case["truth"], "--config", config,
                        "--out", out, "--reference", REFERENCE, "--flank-bp", str(FLANK_BP), "--threads", THREADS], check=True)
    subprocess.run(["bash", os.path.join(bi.SIM, "e2e", "align_cnv.sh"), out, THREADS], check=True)
    return 0


def make_inputs(index):
    """RCK-ecDNA's inputs from the BAMs: obs_from_bam.py, then build_subgraphs_e2e.py (union subgraph and DIAG)."""
    case = cases()[index]
    if inputs_done(case):
        print("done already: " + case["name"])
        return 0
    if not reads_done(case):
        print("{c}: reads not ready".format(c=case["name"]))
        return 2
    obs_run = os.path.join(OBS_RUNS, case["name"])
    env = dict(os.environ, PATH=os.path.dirname(ALIGN_PYTHON) + os.pathsep + os.environ.get("PATH", ""))
    subprocess.run([ALIGN_PYTHON, os.path.join(bi.SIM, "e2e", "obs_from_bam.py"), "--reads", reads_dir(case), "--run", case["truth"],
                    "--out", obs_run, "--threads", THREADS], check=True, env=env)
    for out, extra in ((os.path.join(SUBGRAPHS, case["name"]), []), (os.path.join(SUBGRAPHS_CNP, case["name"]), ["--cnp"])):
        subprocess.run([SIM_PYTHON, os.path.join(bi.SIM, "e2e", "build_subgraphs_e2e.py"), "--obs-run", obs_run,
                        "--out", out] + extra, check=True)
    return 0


# ---------------------------------------------------------------------------------------------- #
# methods
# ---------------------------------------------------------------------------------------------- #

def _licence_retry(argv, outdir, log):
    rc, licence = 1, False
    for attempt in range(1, bi.LICENCE_RETRIES + 1):
        rc, licence = gr.subprocess_run(argv, outdir, log)
        if rc == 0 or not licence:
            return rc, attempt
        print("attempt {k}: licence; retry in 2 min".format(k=attempt), flush=True)
        time.sleep(120)
    return rc, bi.LICENCE_RETRIES


def run(index, method):
    table = cases()
    if index >= len(table):
        print("no case {i}".format(i=index))
        return 0
    case = table[index]
    outdir = run_dir(case, method)
    if done(case, method):
        print("done already: " + outdir)
        return 0
    if not reads_done(case):
        print("{c}: reads not ready ({r})".format(c=case["name"], r=reads_dir(case)))
        return 2
    if method in RCK_METHODS:
        return run_rck(case, method)
    if method == "ce_e2e" and not done(case, "coral_e2e"):
        print("{c}: coral_e2e not done; ce_e2e decomposes its graphs".format(c=case["name"]))
        return 2
    os.makedirs(outdir, exist_ok=True)
    rd = reads_dir(case)
    record = {"case": case["name"], "method": method, "host": socket.gethostname(), "threads": int(THREADS),
              "slurm_job": "{j}_{t}".format(j=os.environ.get("SLURM_ARRAY_JOB_ID", ""), t=os.environ.get("SLURM_ARRAY_TASK_ID", "")),
              "rck_commit": gr.git_commit(), "reads": rd, "end_to_end": True}
    started, ok, attempts = time.time(), False, 0
    if method == "coral_e2e":
        cns, bam = os.path.join(rd, "cnvkit", "tumor.cns"), os.path.join(rd, "tumor.bam")
        seed = [CORAL, "seed", "--cn-seg", cns, "--lr-bam", bam, "--centromere-file", CENTROMERES,
                "--output-prefix", os.path.join(outdir, "seed")]
        rc = gr.subprocess_run(seed, outdir, "coral_seed.stdout")[0]
        seeds = os.path.join(outdir, "seed_CNV_SEEDS.bed")
        n_seeds = sum(1 for _ in open(seeds)) if rc == 0 and os.path.exists(seeds) else 0
        recon = [CORAL, "reconstruct", "--lr-bam", bam, "--cnv-seed", seeds, "--cn-seg", cns,
                 "--output-prefix", os.path.join(outdir, "coral"), "--solver", "gurobi_direct",
                 "--solver-time-limit", str(gr.TIME_LIMIT), "--solver-threads", THREADS]
        record.update(version=gr.CORAL_VERSION, seed_command=seed, command=recon, n_seeds=n_seeds,
                      seed_defaults={"gain": 6.0, "min_seed_size": 100000, "max_seg_gap": 300000})
        if rc == 0 and n_seeds:
            rc, attempts = _licence_retry(recon, outdir, "coral_reconstruct.stdout")
        record["n_amplicons"] = len(coral_amplicons(case))
        ok = rc == 0  # no seed / no amplicon is a valid answer: nothing reconstructed
    else:
        meta_all, ok = [], True
        for n, (graph, _) in coral_amplicons(case).items():
            sub = os.path.join(outdir, "amplicon{n}".format(n=n))
            os.makedirs(sub, exist_ok=True)
            with open(os.path.join(sub, "gurobi.env"), "w") as handle:
                handle.write("Threads {t}\n".format(t=THREADS))
            argv, meta = gr.run_ce(graph, sub, n)
            rc, k = _licence_retry(argv, sub, "ce.stdout")
            attempts = max(attempts, k)
            meta.update(input_graph=graph, status="ok" if rc == 0 and os.path.exists(gr.cycles_file("ce", sub, n)) else "failed")
            meta_all.append(meta)
            ok = ok and meta["status"] == "ok"
        record.update(version=gr.CE_VERSION, amplicons=meta_all, gurobi_env="Threads {t}".format(t=THREADS))
    record.update(status="ok" if ok else "failed", attempts=attempts, elapsed_s=round(time.time() - started, 1))
    with open(os.path.join(outdir, "run.json"), "w") as handle:
        json.dump(record, handle, indent=2)
    print("{c} {m}: {s} in {t:.1f} s".format(c=case["name"], m=method, s=record["status"], t=record["elapsed_s"]), flush=True)
    return 0 if ok else 1


def run_rck(case, method):
    """rck-ecdna on the BAM-derived subgraph (bulk) or DIAG (multi-clone), with bench_iw.rck_argv's flags."""
    import gurobipy
    from rck.ecdna_run import build_parser, run as rck_run

    if not inputs_done(case):
        print("{c}: RCK inputs not built (bench_iw_e2e.py inputs)".format(c=case["name"]))
        return 2
    outdir = run_dir(case, method)
    os.makedirs(outdir, exist_ok=True)
    bi.THREADS = THREADS
    argv = bi.rck_argv(e2e_case(case), RCK_METHODS[method], outdir)
    record = {"case": case["name"], "method": method, "host": socket.gethostname(), "threads": int(THREADS),
              "slurm_job": "{j}_{t}".format(j=os.environ.get("SLURM_ARRAY_JOB_ID", ""), t=os.environ.get("SLURM_ARRAY_TASK_ID", "")),
              "rck_commit": gr.git_commit(), "end_to_end": True, "command": ["python", "-m", "rck.ecdna_run"] + argv,
              "gurobi": ".".join(map(str, gurobipy.gurobi.version()))}
    started, ok, attempts, other = time.time(), False, 0, 0
    while True:
        attempts += 1
        try:
            rck_run(build_parser().parse_args(argv))
            ok = True
            break
        except Exception as error:  # noqa: BLE001 - licence hiccups: retry, then report
            licence = bi._licence(repr(error))
            other += 0 if licence else 1
            print("attempt {k} failed{l}: {e!r}".format(k=attempts, l=" (licence)" if licence else "", e=error), flush=True)
            if attempts >= bi.LICENCE_RETRIES or other >= bi.RETRIES:
                break
            time.sleep(120)
    try:
        gurobipy.disposeDefaultEnv()
    except Exception:  # noqa: BLE001
        pass
    record.update(status="ok" if ok else "failed", attempts=attempts, elapsed_s=round(time.time() - started, 1))
    with open(os.path.join(outdir, "run.json"), "w") as handle:
        json.dump(record, handle, indent=2)
    print("{c} {m}: {s} in {t:.1f} s".format(c=case["name"], m=method, s=record["status"], t=record["elapsed_s"]), flush=True)
    return 0 if ok else 1


def status():
    rows = []
    for c in cases():
        row = {"case": c["name"], "reads": "done" if reads_done(c) else
               ("emitted" if os.path.exists(os.path.join(reads_dir(c), "manifest.json")) else "missing")}
        row["inputs"] = "done" if inputs_done(c) else "missing"
        for m in METHODS:
            path = os.path.join(run_dir(c, m), "run.json")
            row[m] = json.load(open(path)).get("status") if os.path.exists(path) else "missing"
        rows.append(row)
    print(pd.DataFrame(rows).to_string(index=False))


# ---------------------------------------------------------------------------------------------- #
# scoring
# ---------------------------------------------------------------------------------------------- #

def pool_amplicons(files, dst, drop_paths):
    """One AA cycles file from several amplicons' files: Interval and Segment lines concatenated with
    segment ids renumbered, every Cycle/Path kept (Path= read as Cycle=, as score_h2h.normalise_cycles;
    linear walks removed when drop_paths)."""
    intervals, segments, walks, offset = [], [], [], 0
    for path in files:
        local = {}
        for line in open(path):
            if line.startswith("Interval"):
                parts = line.rstrip("\n").split("\t")
                intervals.append(parts[1:])
            elif line.startswith("Segment"):
                parts = line.rstrip("\n").split("\t")
                local[parts[1]] = str(offset + len(local) + 1)
                segments.append([local[parts[1]]] + parts[2:])
            elif line.startswith(("Cycle=", "Path=")):
                fields = dict(f.split("=", 1) for f in line.strip().split(";") if "=" in f)
                key = "Cycle" if "Cycle" in fields else "Path"
                segs = fields.get("Segments", "").split(",")
                linear = key == "Path" or segs[0][:-1] == "0"
                if linear and drop_paths:
                    continue
                fields.pop(key)
                fields["Segments"] = ",".join(s if s[:-1] == "0" else local[s[:-1]] + s[-1] for s in segs if s)
                walks.append(fields)
        offset += len(local)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as out:
        for i, iv in enumerate(intervals, 1):
            out.write("\t".join(["Interval", str(i)] + iv[1:]) + "\n")
        out.write("List of cycle segments\n")
        for s in segments:
            out.write("\t".join(["Segment"] + s) + "\n")
        for i, f in enumerate(walks, 1):
            out.write(";".join(["Cycle={i}".format(i=i)] + ["{k}={v}".format(k=k, v=v) for k, v in f.items()]) + "\n")
    return dst


def pool_graphs(graphs, dst):
    """One AA graph from several amplicons' graphs: the edge lines are coordinate-based, so they concatenate."""
    seq, bp, header_seq, header_bp = [], [], None, None
    for path in graphs:
        section = None
        for line in open(path):
            if line.startswith("SequenceEdge:"):
                header_seq, section = line, "seq"
            elif line.startswith("BreakpointEdge:"):
                header_bp, section = line, "bp"
            elif line.startswith("PathConstraint:"):
                section = "path"  # path constraints refer to per-amplicon edges; not pooled
            elif line.strip() and section == "seq":
                seq.append(line)
            elif line.strip() and section == "bp":
                bp.append(line)
    os.makedirs(os.path.dirname(dst), exist_ok=True)
    with open(dst, "w") as out:
        out.writelines([header_seq or "SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped\n"] + seq)
        out.writelines([header_bp or "BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence\n"] + bp)
    return dst


def mixture_dir(case):
    return os.path.join(MIXTURE, case["name"])


def make_mixture(case):
    out = mixture_dir(case)
    if os.path.exists(os.path.join(out, "mixture_amplicons.tsv")):
        return case["name"], "done already"
    subprocess.run([SIM_PYTHON, os.path.join(bi.SIM, "truth_mixture.py"), case["truth"], "--outdir", out],
                   check=True, stdout=subprocess.DEVNULL)
    return case["name"], "ok"


def score_folder(case, method, variant):
    return os.path.join(run_dir(case, method), "score_nolinear" if variant == "walks_dropped" else "score")


def score_one(task):
    case, method, variant = task
    folder = score_folder(case, method, variant)
    label = bi.method_variant(method, variant)
    try:
        from rck.ecdna import aa_eval, bench_metrics as bm
        from rck.ecdna import evaluate as ev

        if os.path.exists(folder):
            shutil.rmtree(folder)
        if method in RCK_METHODS:
            return score_rck(case, method, folder, label)
        cycles = pool_amplicons(answer_files(case, method), os.path.join(folder, "answer", case["name"] + "_" + method + "_cycles.txt"),
                                drop_paths=variant == "walks_dropped")
        graphs = [g for g, _ in coral_amplicons(case).values()]
        argv = ["--cycles", cycles, "--truth-dir", case["truth"], "--outdir", os.path.join(folder, "aa_eval"),
                "--case", "{c}:{m}".format(c=case["name"], m=label)]
        if graphs:
            argv += ["--graph", pool_graphs(graphs, os.path.join(folder, "answer", case["name"] + "_" + method + "_graph.txt"))]
        if os.path.exists(os.path.join(mixture_dir(case), "mixture_amplicons.tsv")):
            argv += ["--mixture-dir", mixture_dir(case)]
        aa_eval.main(argv)
        truth = ev.load_truth(case["truth"])
        summary, rows = bm.score(cycles, truth, None)
        summary.update(case=case["name"], method=method, variant=variant, n_amplicons=len(graphs))
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


def score_rck(case, method, folder, label):
    """As bench_iw.score_one for rck_bulk / rck_mc, on the e2e run and the BAM-derived inputs; plus --mixture-dir."""
    from rck.ecdna import aa_eval, aa_export, bench_metrics as bm
    from rck.ecdna import evaluate as ev
    from rck.ecdna.rerank import rerank

    c = e2e_case(case)
    outdir = run_dir(case, method)
    answer = os.path.join(folder, "answer")
    if method == "rck_bulk_e2e":
        index = int(rerank(outdir, [bi.SCORE_LAMBDA])["pool_index"].iloc[0])
        paths = aa_export.export(outdir, answer, index, alternatives=2, name=case["name"] + "_" + method)
        sub = aa_export.export_subgraph(c["subgraph"], os.path.join(folder, "graph_export"), name=case["name"])
        graph = os.path.join(folder, "graph", case["name"] + "_amplicon1_graph.txt")
        os.makedirs(os.path.dirname(graph), exist_ok=True)
        aa_export.baseline_graph(outdir, sub["graph"] if isinstance(sub, dict) else sub, graph)
        paths["graph"] = graph
    else:
        paths = aa_export.export(outdir, answer, None, alternatives=0, name=case["name"] + "_" + method)
    argv = ["--cycles", paths["cycles"], "--truth-dir", case["truth"], "--outdir", os.path.join(folder, "aa_eval"),
            "--case", "{c}:{m}".format(c=case["name"], m=label)]
    argv += ["--alt-cycles"] + paths["alt_cycles"] if paths.get("alt_cycles") else []
    argv += ["--graph", paths["graph"]] if paths.get("graph") else []
    fragments = os.path.join(c["rck_input"], "rck.fragments.tsv")
    if method == "rck_mc_e2e" and os.path.exists(fragments):
        argv += ["--phase-blocks", fragments]
    if os.path.exists(os.path.join(mixture_dir(case), "mixture_amplicons.tsv")):
        argv += ["--mixture-dir", mixture_dir(case)]
    aa_eval.main(argv)
    weights = props = None
    if method == "rck_mc_e2e":
        props = bm.read_proportions(os.path.join(c["rck_input"], "clone_proportions.tsv"))
        weights = bm.multiclone_bulk_weights(outdir, props)
    truth = ev.load_truth(case["truth"])
    summary, rows = bm.score(paths["cycles"], truth, weights)
    if method == "rck_bulk_e2e":
        summary.update(bm.total_cn_all(bm.truth_total_cn(case["truth"]), bm.bulk_total_cn(outdir), truth))
    else:
        summary.update(bm.total_cn_all(bm.truth_total_cn(case["truth"]), bm.multiclone_total_cn(outdir, props), truth))
    summary.update(case=case["name"], method=method, variant="")
    with open(os.path.join(folder, "bench_metrics.json"), "w") as handle:
        json.dump(summary, handle, indent=2)
    pd.DataFrame(rows).assign(case=case["name"], method=method, method_variant=label).to_csv(
        os.path.join(folder, "bench_proportions.tsv"), sep="\t", index=False)
    return case["name"], label, "ok"


def score(jobs, index=None):
    """Score every finished run, or with index only that case's (so concurrent pilot tasks never share a folder)."""
    chosen = cases() if index is None else [cases()[index]]
    tasks = [(c, m, v) for c in chosen for m in METHODS if done(c, m) for v in VARIANTS[m]]
    with Pool(jobs) as pool:
        for name, method, state in pool.imap_unordered(score_one, tasks):
            print(name, method, state, flush=True)


def collect():
    """bench_iw.collect on this set's runs and folders (its globals pointed here, nothing written there)."""
    bi.METHODS, bi.VARIANTS = METHODS, VARIANTS
    bi.TABLE, bi.POOLED, bi.PROPORTIONS = TABLE, POOLED, PROPORTIONS
    bi.GRAPH_SOURCE = "end-to-end"
    bi.score_folder = score_folder
    bi.run_dir = run_dir
    bi.collect()
    table = pd.read_csv(TABLE, sep="\t")
    table["baseline"] = "end-to-end"
    extra = []
    for _, r in table.iterrows():
        case = next(c for c in cases() if c["name"] == r["name"])
        ev = os.path.join(score_folder(case, r["method"], r["method_variant"].split(":")[-1] if ":" in r["method_variant"] else ""), "aa_eval")
        gpath = os.path.join(ev, "aa_graph_metrics.tsv")
        g = bi._read(gpath) if os.path.exists(gpath) and os.path.getsize(gpath) > 1 else None  # rck_mc: no AA graph, empty file
        rec = json.load(open(os.path.join(run_dir(case, r["method"]), "run.json"))) if os.path.exists(os.path.join(run_dir(case, r["method"]), "run.json")) else {}
        row = {"n_seeds": rec.get("n_seeds"), "n_amplicons": rec.get("n_amplicons")}
        if g is not None and len(g):
            row.update({"graph_" + k: g[k].iloc[0] for k in g.columns if k not in ("case", "amplicon")})
        extra.append(row)
    pd.concat([table, pd.DataFrame(extra)], axis=1).to_csv(TABLE, sep="\t", index=False)


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = parser.add_subparsers(dest="cmd", required=True)
    lst = sub.add_parser("list")
    lst.add_argument("--method", choices=METHODS + ["reads", "inputs"], default="reads")
    rd = sub.add_parser("reads")
    rd.add_argument("index", type=int)
    ip = sub.add_parser("inputs")
    ip.add_argument("index", type=int)
    r = sub.add_parser("run")
    r.add_argument("index", type=int)
    r.add_argument("--method", choices=METHODS, required=True)
    sub.add_parser("status")
    for name in ("mixture", "score"):
        s = sub.add_parser(name)
        s.add_argument("--jobs", type=int, default=8)
        if name == "score":
            s.add_argument("--case", type=int, default=None, help="only this case index")
    sub.add_parser("collect")
    args = parser.parse_args()
    if args.cmd == "list":
        for i, case in enumerate(cases()):
            print(i, case["name"], args.method, sep="\t")
    elif args.cmd == "reads":
        return make_reads(args.index)
    elif args.cmd == "inputs":
        return make_inputs(args.index)
    elif args.cmd == "run":
        return run(args.index, args.method)
    elif args.cmd == "status":
        status()
    elif args.cmd == "mixture":
        with Pool(args.jobs) as pool:
            for name, state in pool.imap_unordered(make_mixture, cases()):
                print(name, state, flush=True)
    elif args.cmd == "score":
        score(args.jobs, args.case)
    else:
        collect()
    return 0


if __name__ == "__main__":
    sys.exit(main())

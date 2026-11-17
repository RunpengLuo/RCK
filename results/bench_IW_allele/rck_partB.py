#!/usr/bin/env python3
"""bench_IW_allele Part B, steps 3-4: the multi-clone route on RCK's degraded inputs (PLAN.md "Part B").

One task per (variant, case): ``rck-ecdna --model multiclone`` on
/n/fs/yubi-iw/runs/bench_IW_intpsi_rckstyle/<case>/<variant>/rck_input/ (make_rck_inputs.py), with the
bench_IW_intpsi rck_mc settings (bench_iw.rck_argv) except the time limit per MILP (PART_B_TIME_LIMIT; the
reference runs solved optimally in <= 5 s against 60 s, the label-free inputs may need longer; every run
records its solver status). P = 0.9 is the default --mc-nas-fp 0.1.

    python rck_partB.py list                 # task index, variant, case
    python rck_partB.py run IDX              # one task (Gurobi; the array runs these, submit_partB.sh)
    python rck_partB.py status               # done / failed / missing per variant
    python rck_partB.py export [--jobs 8]    # AA cycles of every finished run, for scoring (no Gurobi)

Outputs: <case>/<variant>/rck_mc/ (the run) and <case>/<variant>/rck_mc/score/answer/ (exported cycles).
Gurobi licence errors are retried as bench_iw.py does (LICENCE_RETRIES x 2 min).
"""
import argparse
import json
import os
import socket
import sys
import time
import traceback
from multiprocessing import Pool

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import allele_bench as ab  # noqa: E402  (configures bench_iw for iw_intpsi)
import make_rck_inputs as mi  # noqa: E402

bi = ab.bi
PART_B_TIME_LIMIT = 600  # s per MILP
THREADS = os.environ.get("SLURM_CPUS_PER_TASK", "4")


def tasks():
    """Variant-major order, so a partial array covers whole variants first."""
    from rck.ecdna.rck_inputs import VARIANTS
    return [(variant, case) for variant in VARIANTS for case in bi.cases()]


def run_dir(case, variant):
    return os.path.join(mi.OUT, case["name"], variant, "rck_mc")


def done(case, variant):
    d = run_dir(case, variant)
    return all(os.path.exists(os.path.join(d, f)) for f in ("summary.json", "decompositions.tsv"))


def argv_of(case, variant):
    return ["--model", "multiclone", "--rck-input", mi.input_dir(case, variant), "--truth-dir", case["truth"],
            "--outdir", run_dir(case, variant), "--no-plots", "--threads", THREADS, "--time-limit", str(PART_B_TIME_LIMIT)]


def run(index):
    table = tasks()
    if index >= len(table):
        print("no task {i} ({n} tasks)".format(i=index, n=len(table)))
        return 0
    variant, case = table[index]
    outdir = run_dir(case, variant)
    if done(case, variant):
        print("done already: " + outdir)
        return 0
    os.makedirs(outdir, exist_ok=True)
    import gurobipy
    from rck.ecdna_run import build_parser, run as rck_run

    argv = argv_of(case, variant)
    record = {"case": case["name"], "variant": variant, "method": "rck_mc", "host": socket.gethostname(),
              "threads": int(THREADS), "command": ["python", "-m", "rck.ecdna_run"] + argv,
              "slurm_job": "{j}_{t}".format(j=os.environ.get("SLURM_ARRAY_JOB_ID", ""), t=os.environ.get("SLURM_ARRAY_TASK_ID", "")),
              "gurobi": ".".join(map(str, gurobipy.gurobi.version()))}
    ok, attempts, other, started = False, 0, 0, time.time()
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
    print("{c} {v}: {s} in {t:.1f} s".format(c=case["name"], v=variant, s=record["status"], t=record["elapsed_s"]), flush=True)
    return 0 if ok else 1


def status():
    counts = {}
    for variant, case in tasks():
        d = run_dir(case, variant)
        state = "done" if done(case, variant) else "failed" if os.path.exists(os.path.join(d, "run.json")) else "missing"
        if state == "done":
            with open(os.path.join(d, "summary.json")) as handle:
                s = json.load(handle)
            state = "done (optimal)" if s.get("status") == 2 else "done (status {x})".format(x=s.get("status"))
        counts.setdefault(variant, {}).setdefault(state, 0)
        counts[variant][state] += 1
    for variant, c in counts.items():
        print(variant, c)


def export_one(task):
    variant, case = task
    if not done(case, variant):
        return case["name"], variant, "not done"
    try:
        from rck.ecdna import aa_export
        folder = os.path.join(run_dir(case, variant), "score", "answer")
        aa_export.export(run_dir(case, variant), folder, None, alternatives=0, name=case["name"] + "_rck_mc")
        return case["name"], variant, "ok"
    except Exception:  # noqa: BLE001
        with open(os.path.join(run_dir(case, variant), "export_error.txt"), "w") as handle:
            handle.write(traceback.format_exc())
        return case["name"], variant, "error"


def main():
    p = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    sub = p.add_subparsers(dest="cmd", required=True)
    sub.add_parser("list")
    r = sub.add_parser("run")
    r.add_argument("index", type=int)
    sub.add_parser("status")
    e = sub.add_parser("export")
    e.add_argument("--jobs", type=int, default=8)
    a = p.parse_args()
    if a.cmd == "list":
        for i, (variant, case) in enumerate(tasks()):
            print(i, variant, case["name"], sep="\t")
    elif a.cmd == "run":
        return run(a.index)
    elif a.cmd == "status":
        status()
    else:
        import logging
        logging.disable(logging.INFO)
        states = {}
        with Pool(a.jobs) as pool:
            for name, variant, state in pool.imap_unordered(export_one, tasks()):
                states[state] = states.get(state, 0) + 1
                if state == "error":
                    print(name, variant, state)
        print(states)
    return 0


if __name__ == "__main__":
    sys.exit(main())

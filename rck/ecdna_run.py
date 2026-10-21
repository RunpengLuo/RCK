"""``rck-ecdna``: reconstruct allele-specific ecDNA structures from one breakpoint subgraph.

Example (demo3 simulation)::

    rck-ecdna \\
        --subgraph-dir subgraphs/demo3/amplicon_0 \\
        --obs-dir simulations/runs/demo3_clone1_ecdna/obs \\
        --truth-dir simulations/runs/demo3_clone1_ecdna \\
        --outdir rck_ecdna_out/amplicon_0

See docs/ecDNA.md for the model and every option.
"""

from __future__ import annotations

import argparse
import logging
import os
import sys
import time
from typing import Dict, List, Optional

from rck.ecdna.decompose import DecompositionParams, decompose
from rck.ecdna.gurobi_env import default_threads, quiet_gurobi_logger
from rck.ecdna.ilp_circular import CircularParams, solve_flow
from rck.ecdna.output import write_all
from rck.ecdna.scale import ScaleParams, build_observations
from rck.ecdna.subgraph_io import read_subgraph
from rck.ecdna.walk import read_paths

log = logging.getLogger("rck.ecdna")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rck-ecdna", description="Allele-specific ecDNA reconstruction (circular RCK).",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    io = parser.add_argument_group("input / output")
    io.add_argument("--subgraph-dir", required=True, help="Directory with nodes.tsv / edges.tsv of one amplicon")
    io.add_argument("--obs-dir", default=None, help="Observation directory; sets --windows-dir OBS/windows, --extremities OBS/extremity.tsv "
                                                    "and --paths OBS/path.tsv when those exist")
    io.add_argument("--windows-dir", default=None, help="UGP windows directory (window.tsv.gz, window.depth.npz, sample_ids.tsv)")
    io.add_argument("--extremities", default=None, help="extremity.tsv with genome-wide #REF counts (for junction reads per copy)")
    io.add_argument("--paths", default=None, help="path.tsv with long-read paths (ranks alternative orderings)")
    io.add_argument("--hatchet-seg", default=None, help="HATCHet best.seg.ucn; fits depth per copy when there is no matched normal")
    io.add_argument("--sample", default=None, help="Tumor sample name in sample_ids.tsv / HATCHet SAMPLE column")
    io.add_argument("--outdir", required=True)
    io.add_argument("--name", default=None, help="Name used in output files (default: subgraph directory name)")

    sc = parser.add_argument_group("copy-number scale")
    sc.add_argument("--depth-per-copy", type=float, default=None, help="Override: tumor depth of one copy")
    sc.add_argument("--junction-reads-per-copy", type=float, default=None, help="Override: reads spanning one junction copy")
    sc.add_argument("--depth-mode", choices=["tumor", "ratio"], default="tumor", help="tumor depth / gamma, or tumor/normal ratio")
    sc.add_argument("--baseline-cn", type=float, default=2.0, help="Copy number of the dominant balanced background")
    sc.add_argument("--block-bp", type=int, default=1_000_000)
    sc.add_argument("--min-snps-allele", type=int, default=5, help="SNPs needed for an allele-specific segment observation")

    gr = parser.add_argument_group("graph")
    gr.add_argument("--phantom-max-bp", type=int, default=2, help="Contract segments up to this length between two reference joins (0 = off)")
    gr.add_argument("--soft-junction-hap", action="store_true", help="Junction haplotype labels become penalties instead of constraints")

    fl = parser.add_argument_group("flow MILP")
    fl.add_argument("--imputed-penalty", type=float, default=2.0)
    fl.add_argument("--novel-penalty", type=float, default=0.05)
    fl.add_argument("--alpha-novel", type=float, default=0.5)
    fl.add_argument("--alpha-ref", type=float, default=0.25)
    fl.add_argument("--bg-max", type=int, default=4, help="Per-haplotype background copy-number cap")
    fl.add_argument("--bg-prior-weight", type=float, default=0.1)
    fl.add_argument("--soft-junction-penalty", type=float, default=1.0)
    fl.add_argument("--enforce-is", action="store_true", help="Per-haplotype extremity exclusivity (RCK infinite sites)")
    fl.add_argument("--no-background-novel", action="store_true", help="Background layer may use reference adjacencies only")
    fl.add_argument("--linear-ends", choices=["boundaries", "all", "none"], default="boundaries",
                    help="Where amplified copies may start or stop: copy-number boundaries and flanks, any extremity, "
                         "or nowhere (strictly circular)")
    fl.add_argument("--linear-end-penalty", type=float, default=0.25, help="Cost per copy per linear end")
    fl.add_argument("--boundary-min-copies", type=float, default=2.0, help="Copy-number step that marks a boundary")
    fl.add_argument("--boundary-frac", type=float, default=0.25, help="... or this fraction of the larger copy number")

    dc = parser.add_argument_group("decomposition")
    dc.add_argument("--max-species", type=int, default=4)
    dc.add_argument("--max-multiplicity", type=int, default=3, help="Largest number of times one species may use an edge")
    dc.add_argument("--species-penalty", type=float, default=1.0)
    dc.add_argument("--min-species-weight", type=int, default=None,
                    help="Species lighter than this are flagged low_weight (default: --min-species-weight-frac of the heaviest)")
    dc.add_argument("--min-species-weight-frac", type=float, default=0.10)
    dc.add_argument("--min-ecdna-copies", type=int, default=2, help="Circular species lighter than this do not make an ecDNA call")
    dc.add_argument("--ecdna-call-fraction", type=float, default=0.8,
                    help="Circular share of copy-number x length needed for an ecDNA call (linear below 1 - this)")
    dc.add_argument("--num-alternatives", type=int, default=10)
    dc.add_argument("--pool-gap", type=float, default=0.10)
    dc.add_argument("--alternatives-mode", choices=["pool", "nogood"], default="pool")
    dc.add_argument("--max-orderings", type=int, default=20, help="Alternative orderings enumerated per species")

    ev = parser.add_argument_group("evaluation (optional)")
    ev.add_argument("--truth-dir", default=None, help="Simulator run directory with truth/ and obs/ (amplicon-allele-analysis)")
    ev.add_argument("--truth-cycles", default=None, help="Explicit truth table: species order chrom start end strand hap weight")
    ev.add_argument("--truth-sample", default=None, help="Sample in truth/proportions.tsv for the truth copy numbers (default: --sample, else the first tumor)")
    ev.add_argument("--eval-tol-bp", type=int, default=10_000, help="Coordinate tolerance for junction and order matching")
    ev.add_argument("--min-truth-fraction", type=float, default=0.5,
                    help="Score only truth species with at least this share of their bases inside this subgraph")
    ev.add_argument("--all-truth-species", action="store_true", help="Score every truth species, inside this subgraph or not")

    gu = parser.add_argument_group("Gurobi")
    gu.add_argument("--threads", type=int, default=None, help="Default: SLURM_CPUS_PER_TASK, else the CPUs this process may use")
    gu.add_argument("--time-limit", type=float, default=None, help="Seconds, per MILP")
    gu.add_argument("--mip-gap", type=float, default=1e-4)
    gu.add_argument("-v", "--verbose", action="store_true")
    return parser


def _resolve_inputs(args) -> Dict[str, Optional[str]]:
    windows, extremities, paths = args.windows_dir, args.extremities, args.paths
    if args.obs_dir:
        candidate = os.path.join(args.obs_dir, "windows")
        windows = windows or (candidate if os.path.isdir(candidate) else None)
        candidate = os.path.join(args.obs_dir, "extremity.tsv")
        extremities = extremities or (candidate if os.path.exists(candidate) else None)
        candidate = os.path.join(args.obs_dir, "path.tsv")
        paths = paths or (candidate if os.path.exists(candidate) else None)
    return {"windows": windows, "extremities": extremities, "paths": paths}


def run(args) -> Dict[str, object]:
    started = time.time()
    args.threads = default_threads(args.threads)
    quiet_gurobi_logger(getattr(args, "verbose", False))
    os.makedirs(args.outdir, exist_ok=True)
    inputs = _resolve_inputs(args)
    graph = read_subgraph(args.subgraph_dir, phantom_max_bp=args.phantom_max_bp, soft_junction_hap=args.soft_junction_hap, name=args.name)
    log.info("graph %s: %d segments, %d adjacencies; Gurobi threads %d", graph.name, len(graph.segments), len(graph.adjacencies), args.threads)

    scale_params = ScaleParams(depth_per_copy=args.depth_per_copy, junction_reads_per_copy=args.junction_reads_per_copy,
                               baseline_cn=args.baseline_cn, depth_mode=args.depth_mode, sample=args.sample,
                               hatchet_seg=args.hatchet_seg, block_bp=args.block_bp, min_snps_allele=args.min_snps_allele)
    obs = build_observations(graph, scale_params, windows_dir=inputs["windows"], extremities_path=inputs["extremities"])

    flow_params = CircularParams(imputed_penalty=args.imputed_penalty, novel_penalty=args.novel_penalty, alpha_novel=args.alpha_novel,
                                 alpha_ref=args.alpha_ref, bg_max=args.bg_max, bg_prior_weight=args.bg_prior_weight,
                                 soft_junction_penalty=args.soft_junction_penalty, enforce_is=args.enforce_is,
                                 allow_background_novel=not args.no_background_novel, linear_ends=args.linear_ends,
                                 linear_end_penalty=args.linear_end_penalty, boundary_min_copies=args.boundary_min_copies,
                                 boundary_frac=args.boundary_frac, time_limit=args.time_limit,
                                 threads=args.threads, mip_gap=args.mip_gap, log_file=os.path.join(args.outdir, "gurobi_flow.log"),
                                 log_to_console=args.verbose)
    sol = solve_flow(graph, obs, flow_params)
    log.info("flow: %s, objective %.4f, %d ecDNA edge copies, %d linear-end copies, imputed used: %s", sol.status, sol.objective,
             sum(sol.ecdna_edges().values()), sum(sol.linear_ends.values()), ",".join(sol.imputed_used) or "none")

    paths = read_paths(inputs["paths"], graph) if inputs["paths"] else None
    dec_params = DecompositionParams(max_species=args.max_species, max_multiplicity=args.max_multiplicity,
                                     species_penalty=args.species_penalty, min_species_weight=args.min_species_weight,
                                     min_species_weight_frac=args.min_species_weight_frac, min_ecdna_copies=args.min_ecdna_copies,
                                     ecdna_call_fraction=args.ecdna_call_fraction, num_alternatives=args.num_alternatives,
                                     pool_gap=args.pool_gap, mode=args.alternatives_mode, max_orderings=args.max_orderings,
                                     time_limit=args.time_limit, threads=args.threads, mip_gap=args.mip_gap,
                                     log_file=os.path.join(args.outdir, "gurobi_decompose.log"), log_to_console=args.verbose)
    decompositions = decompose(graph, sol.ecdna_edges(), dec_params, paths=paths)
    if decompositions and decompositions[0].residual > 0:
        log.warning("best decomposition leaves residual %.3f; consider raising --max-species or --max-multiplicity",
                    decompositions[0].residual)
    outputs = write_all(args.outdir, graph, obs, sol, decompositions,
                        extra={"inputs": inputs, "runtime_s": round(time.time() - started, 2),
                               "paths_loaded": 0 if paths is None else len(paths)})
    evaluation = None
    if args.truth_dir or args.truth_cycles:
        evaluation = _evaluate(args, graph, outputs)
    for d in decompositions[:3]:
        log.info("decomposition %d: %s, objective %.4f, %d species (%s)", d.rank, d.call, d.objective, len(d.species),
                 ", ".join("{t} w={w} len={l}{f}".format(t=s.topology, w=s.weight, l=s.length_bp, f=" low_weight" if s.low_weight else "")
                           for s in d.species))
    return {"graph": graph, "observations": obs, "flow": sol, "decompositions": decompositions, "outputs": outputs, "inputs": inputs,
            "evaluation": evaluation}


def _evaluate(args, graph, outputs) -> Dict[str, object]:
    """Score the run against truth. A problem with the truth files is logged and recorded in
    evaluation.json instead of failing a run whose reconstruction already succeeded."""
    import json
    import traceback

    from rck.ecdna.evaluate import evaluate, load_truth, log_summary, write_evaluation, write_truth_table

    try:
        truth = load_truth(args.truth_dir, args.truth_cycles, args.truth_sample or args.sample)
        write_truth_table(os.path.join(args.outdir, "truth_cycles.tsv"), truth)
        regions = None if args.all_truth_species else graph.regions()
        frame, evaluation = evaluate(outputs["decompositions"], truth, tol=args.eval_tol_bp, regions=regions,
                                     min_truth_fraction=args.min_truth_fraction)
        outputs.update({"evaluation_" + k: v for k, v in write_evaluation(args.outdir, frame, evaluation).items()})
        log_summary(evaluation)
        return evaluation
    except Exception as error:  # noqa: BLE001 - report, do not lose the reconstruction
        log.error("evaluation failed: %s", error)
        evaluation = {"error": str(error), "traceback": traceback.format_exc()}
        path = os.path.join(args.outdir, "evaluation.json")
        with open(path, "w") as handle:
            json.dump(evaluation, handle, indent=2)
        outputs["evaluation_summary"] = path
        return evaluation


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    quiet_gurobi_logger(args.verbose)
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())

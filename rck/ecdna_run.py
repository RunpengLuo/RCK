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
from rck.ecdna.refit import RefitParams, refit_decompositions
from rck.ecdna.scale import ScaleParams, build_observations, load_subgraph_scale
from rck.ecdna.subgraph_io import read_segment_blocks, read_subgraph
from rck.ecdna.walk import read_paths

log = logging.getLogger("rck.ecdna")


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(prog="rck-ecdna", description="Allele-specific ecDNA reconstruction (circular RCK).",
                                     formatter_class=argparse.ArgumentDefaultsHelpFormatter)
    io = parser.add_argument_group("input / output")
    io.add_argument("--subgraph-dir", required=True,
                    help="Directory with nodes.tsv / edges.tsv of one amplicon. When it also has scale.json (written by "
                         "subgraph_construct.py) it is all RCK needs: see --inputs")
    io.add_argument("--inputs", choices=["auto", "subgraph", "windows"], default="auto",
                    help="Where copy numbers come from. subgraph: scale.json, CN_* edge columns, cn_profile.tsv and paths.tsv "
                         "of the subgraph folder (no windows needed). windows: estimate from --windows-dir / --obs-dir. "
                         "auto: subgraph when scale.json exists, else windows")
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
    sc.add_argument("--baseline-cn", type=float, default=None,
                    help="Copy number of the dominant balanced background (default: scale.json's cell ploidy in subgraph mode, else 2)")
    sc.add_argument("--block-bp", type=int, default=1_000_000)
    sc.add_argument("--min-snps-allele", type=int, default=5, help="SNPs needed for an allele-specific segment observation")
    sc.add_argument("--block-flip", choices=["off", "unoriented"], default="off",
                    help="unoriented: allele observations per phase-block piece (segment_blocks.tsv), and the flow MILP chooses "
                         "an A/B flip for every phase block the subgraph left unoriented")
    sc.add_argument("--flip-oriented-penalty", type=float, default=None,
                    help="With --block-flip: let blocks the subgraph oriented flip too, at this cost (default: never)")
    sc.add_argument("--flip-tol", type=float, default=0.01,
                    help="With --block-flip: a block whose reversed flip changes the flow objective by at most this is unresolved")

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
    dc.add_argument("--step-penalty", type=float, default=0.0,
                    help="Cost per segment traversal of a species (walk complexity): favours separate simple walks over "
                         "walks that re-traverse shared sequence, e.g. merged nested cycles. 0 = tie-break only")
    dc.add_argument("--ambiguity-window", type=float, default=2.0,
                    help="Report alternatives to the top decomposition that fit the flow exactly as well and whose objective "
                         "is within this window (summary.json ambiguity); 0 = off")
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

    lr = parser.add_argument_group("Long reads in the decomposition (--paths)")
    lr.add_argument("--path-weight", type=float, default=0.5,
                    help="Cost of leaving all multi-junction read support unexplained, scaled by the unexplained share; 0 = reads only rank orderings")
    lr.add_argument("--path-fraction", type=float, default=0.0,
                    help="Required explained share of that support (CoRAL uses 0.9); dropped with a warning if infeasible")
    lr.add_argument("--min-path-support", type=float, default=1.0, help="Reads behind a path before it constrains the decomposition")
    lr.add_argument("--min-path-junctions", type=int, default=2, help="Junctions a path must span to constrain the decomposition")

    rf = parser.add_argument_group("Copy-number refit")
    rf.add_argument("--no-refit", action="store_true", help="Report the integer weights of the structure search")
    rf.add_argument("--refit-min-copies", type=float, default=0.5,
                    help="Species refitted below this many copies are dropped as unsupported")

    ev = parser.add_argument_group("evaluation (optional)")
    ev.add_argument("--truth-dir", default=None, help="Simulator run directory with truth/ and obs/ (amplicon-allele-analysis)")
    ev.add_argument("--truth-cycles", default=None, help="Explicit truth table: species order chrom start end strand hap weight")
    ev.add_argument("--truth-sample", default=None, help="Sample in truth/proportions.tsv for the truth copy numbers (default: --sample, else the first tumor)")
    ev.add_argument("--eval-tol-bp", type=int, default=10_000, help="Coordinate tolerance for junction and order matching")
    ev.add_argument("--min-truth-fraction", type=float, default=0.5,
                    help="Score only truth species with at least this share of their bases inside this subgraph")
    ev.add_argument("--all-truth-species", action="store_true", help="Score every truth species, inside this subgraph or not")

    pl = parser.add_argument_group("figures (need matplotlib and cnplot)")
    pl.add_argument("--no-plots", action="store_true", help="Do not draw plots/<name>_ecdna.pdf")
    pl.add_argument("--plot-ranks", type=int, default=1, help="Draw the species of this many top decompositions")
    pl.add_argument("--plot-low-weight", action="store_true", help="Also draw species flagged low_weight")
    pl.add_argument("--cnplot-path", default=None,
                    help="cnplot checkout (amplicon-allele-analysis third_party/cnplot) when cnplot is not installed; "
                         "default $CNPLOT_PATH")

    gu = parser.add_argument_group("Gurobi")
    gu.add_argument("--threads", type=int, default=None, help="Default: SLURM_CPUS_PER_TASK, else the CPUs this process may use")
    gu.add_argument("--time-limit", type=float, default=None, help="Seconds, per MILP")
    gu.add_argument("--mip-gap", type=float, default=1e-4)
    gu.add_argument("-v", "--verbose", action="store_true")
    return parser


def _resolve_inputs(args) -> Dict[str, Optional[str]]:
    windows, extremities, paths = args.windows_dir, args.extremities, args.paths
    mode = getattr(args, "inputs", "auto")
    has_scale = os.path.exists(os.path.join(args.subgraph_dir, "scale.json"))
    if mode == "subgraph" and not has_scale:
        raise SystemExit("--inputs subgraph: {d} has no scale.json; rebuild the subgraph with subgraph_construct.py "
                         "(without --no-rck-inputs) or use --inputs windows".format(d=args.subgraph_dir))
    if mode == "subgraph" or (mode == "auto" and has_scale):
        local = os.path.join(args.subgraph_dir, "paths.tsv")
        return {"mode": "subgraph", "windows": None, "extremities": None,
                "paths": paths or (local if os.path.exists(local) and os.path.getsize(local) > 0 else None)}
    if args.obs_dir:
        candidate = os.path.join(args.obs_dir, "windows")
        windows = windows or (candidate if os.path.isdir(candidate) else None)
        candidate = os.path.join(args.obs_dir, "extremity.tsv")
        extremities = extremities or (candidate if os.path.exists(candidate) else None)
        candidate = os.path.join(args.obs_dir, "path.tsv")
        paths = paths or (candidate if os.path.exists(candidate) else None)
    return {"mode": "windows", "windows": windows, "extremities": extremities, "paths": paths}


def run(args) -> Dict[str, object]:
    started = time.time()
    args.threads = default_threads(args.threads)
    quiet_gurobi_logger(getattr(args, "verbose", False))
    os.makedirs(args.outdir, exist_ok=True)
    inputs = _resolve_inputs(args)
    graph = read_subgraph(args.subgraph_dir, phantom_max_bp=args.phantom_max_bp, soft_junction_hap=args.soft_junction_hap, name=args.name)
    log.info("graph %s: %d segments, %d adjacencies; Gurobi threads %d", graph.name, len(graph.segments), len(graph.adjacencies), args.threads)

    subgraph_scale = load_subgraph_scale(args.subgraph_dir) if inputs["mode"] == "subgraph" else None
    baseline_cn = args.baseline_cn
    if baseline_cn is None:
        baseline_cn = float(subgraph_scale.get("background_cn", 2.0)) if subgraph_scale is not None else 2.0
    if subgraph_scale is not None:
        built_for = str(subgraph_scale.get("sample", "?"))
        log.info("inputs: subgraph folder only (scale.json for sample %s; %s)", built_for, subgraph_scale.get("method", ""))
        if not getattr(args, "truth_sample", None) and not args.sample:
            args.truth_sample = built_for  # score against the sample the copy numbers describe
        if args.sample and args.sample != built_for:
            log.warning("--sample %s ignored: the subgraph's copy numbers were computed for %s; rebuild it with "
                        "subgraph_construct.py --sample %s for that sample", args.sample, built_for, args.sample)
        if args.windows_dir or args.obs_dir or args.extremities or args.hatchet_seg:
            log.info("--windows-dir / --obs-dir / --extremities / --hatchet-seg are not used with subgraph inputs "
                     "(pass --inputs windows to use them)")
    scale_params = ScaleParams(depth_per_copy=args.depth_per_copy, junction_reads_per_copy=args.junction_reads_per_copy,
                               baseline_cn=baseline_cn, depth_mode=args.depth_mode, sample=args.sample,
                               hatchet_seg=args.hatchet_seg, block_bp=args.block_bp, min_snps_allele=args.min_snps_allele,
                               block_flip=getattr(args, "block_flip", "off"))
    segment_blocks = read_segment_blocks(args.subgraph_dir, graph) if scale_params.block_flip != "off" else None
    obs = build_observations(graph, scale_params, windows_dir=inputs["windows"], extremities_path=inputs["extremities"],
                             subgraph_scale=subgraph_scale, segment_blocks=segment_blocks)

    flow_params = CircularParams(imputed_penalty=args.imputed_penalty, novel_penalty=args.novel_penalty, alpha_novel=args.alpha_novel,
                                 alpha_ref=args.alpha_ref, bg_max=args.bg_max, bg_prior_weight=args.bg_prior_weight,
                                 soft_junction_penalty=args.soft_junction_penalty, enforce_is=args.enforce_is,
                                 allow_background_novel=not args.no_background_novel, linear_ends=args.linear_ends,
                                 linear_end_penalty=args.linear_end_penalty, boundary_min_copies=args.boundary_min_copies,
                                 boundary_frac=args.boundary_frac, flip_oriented_penalty=getattr(args, "flip_oriented_penalty", None),
                                 flip_tol=getattr(args, "flip_tol", 0.01), time_limit=args.time_limit,
                                 threads=args.threads, mip_gap=args.mip_gap, log_file=os.path.join(args.outdir, "gurobi_flow.log"),
                                 log_to_console=args.verbose)
    sol = solve_flow(graph, obs, flow_params)
    log.info("flow: %s, objective %.4f, %d ecDNA edge copies, %d linear-end copies, imputed used: %s", sol.status, sol.objective,
             sum(sol.ecdna_edges().values()), sum(sol.linear_ends.values()), ",".join(sol.imputed_used) or "none")
    if sol.blocks:
        log.info("phase blocks: %s", ", ".join("{b} {s}".format(b=b["phase_block"], s=b["status"]) for b in sol.blocks
                                                if b["status"] != "fixed" or b["flip"]))

    paths = read_paths(inputs["paths"], graph) if inputs["paths"] else None
    dec_params = DecompositionParams(max_species=args.max_species, max_multiplicity=args.max_multiplicity,
                                     species_penalty=args.species_penalty, step_penalty=getattr(args, "step_penalty", 0.0),
                                     ambiguity_window=getattr(args, "ambiguity_window", 2.0), min_species_weight=args.min_species_weight,
                                     min_species_weight_frac=args.min_species_weight_frac, min_ecdna_copies=args.min_ecdna_copies,
                                     ecdna_call_fraction=args.ecdna_call_fraction, num_alternatives=args.num_alternatives,
                                     pool_gap=args.pool_gap, mode=args.alternatives_mode, max_orderings=args.max_orderings,
                                     path_weight=args.path_weight, path_fraction=args.path_fraction,
                                     min_path_support=args.min_path_support, min_path_junctions=args.min_path_junctions,
                                     time_limit=args.time_limit, threads=args.threads, mip_gap=args.mip_gap,
                                     log_file=os.path.join(args.outdir, "gurobi_decompose.log"), log_to_console=args.verbose)
    decompositions = decompose(graph, sol.ecdna_edges(), dec_params, paths=paths)
    if decompositions and decompositions[0].residual > 0:
        log.warning("best decomposition leaves residual %.3f; consider raising --max-species or --max-multiplicity",
                    decompositions[0].residual)
    refit_decompositions(decompositions, graph, obs, sol, flow_params, dec_params,
                         RefitParams(enabled=not args.no_refit, min_copies=args.refit_min_copies))
    outputs = write_all(args.outdir, graph, obs, sol, decompositions,
                        extra={"inputs": inputs, "runtime_s": round(time.time() - started, 2),
                               "paths_loaded": 0 if paths is None else len(paths)})
    evaluation = None
    if args.truth_dir or args.truth_cycles:
        evaluation = _evaluate(args, graph, outputs)
    if not getattr(args, "no_plots", False):
        _plot(args, outputs)
    for d in decompositions[:3]:
        log.info("decomposition %d: %s, objective %.4f, %d species (%s)%s%s", d.rank, d.call, d.objective, len(d.species),
                 ", ".join("{t} w={w:g} (int {wi}) len={l}{f}".format(t=s.topology, w=s.weight, wi=s.weight_int, l=s.length_bp,
                                                                     f=" low_weight" if s.low_weight else "") for s in d.species),
                 "; pruned " + ", ".join("w={w:g}".format(w=p["refit_weight"]) for p in d.pruned) if d.pruned else "",
                 "; reads explained {r:.0%} of {t:.0f}".format(r=d.reads_explained, t=d.reads_total) if d.reads_total else "")
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


def _plot(args, outputs) -> None:
    """Draw the species pages and the recovery matrix. Plotting never fails a finished run."""
    from rck.ecdna.plot import plot_run

    try:
        written = plot_run(args.outdir, ranks=getattr(args, "plot_ranks", 1), include_low_weight=getattr(args, "plot_low_weight", False),
                           cnplot_path=getattr(args, "cnplot_path", None))
        if written:
            outputs["plots"] = written["pdf"]
    except Exception as error:  # noqa: BLE001 - a figure problem must not lose the reconstruction
        log.warning("plots failed: %s", error, exc_info=log.isEnabledFor(logging.DEBUG))


def main(argv: Optional[List[str]] = None) -> int:
    args = build_parser().parse_args(argv)
    logging.basicConfig(level=logging.DEBUG if args.verbose else logging.INFO, format="%(asctime)s %(name)s %(levelname)s %(message)s")
    quiet_gurobi_logger(args.verbose)
    run(args)
    return 0


if __name__ == "__main__":
    sys.exit(main())

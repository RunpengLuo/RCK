"""The multi-clone RCK-ecDNA route: ``rck-ecdna --model multiclone --rck-input <dir>``.

Reads the annotated DIAG folder that amplicon-allele-analysis ``subgraph_construct.py
--clone-profiles`` writes (``<subgraph>/rck_input/``, RCK's native formats), builds
:class:`rck.core.ilp_gurobi_ecdna.OptModelMultiCloneEcDNA` on it and writes RCK's native per-clone
``rck.scnt.tsv`` / ``rck.acnt.tsv`` exactly as ``rck_run.py`` writes them, plus ``summary.json``.

The inputs are read the way ``rck/rck_run.py --no-pre`` reads them, with RCK's own readers
(``read_segments_from_file``, ``extract_scnt_from_segments``, ``extract_scnb_from_segments``,
``read_scnb_from_file``, ``read_adjacencies_from_file``, ``read_positions_from_file``,
``read_adjacency_groups_from_file``, ``get_segments_for_fragments_ids_dict``,
``get_ref_adjacencies_from_segments``, ``projected_groups``) and the same alignment checks, so
the model sees exactly what plain RCK sees. The decomposition of every clone's karyotype into
cycles (ecDNA) and telomere-to-telomere paths (chromosomes), heaviest first (plan 4.3, step 4),
is not here yet.
"""

from __future__ import annotations

import json
import logging
import os
import time
from dataclasses import dataclass, field
from typing import Dict, List, Optional

from rck.core.io import read_segments_from_file, extract_scnt_from_segments, extract_scnb_from_segments, read_scnb_from_file, \
    read_adjacencies_from_file, read_positions_from_file, read_adjacency_groups_from_file, write_scnt_to_file, write_acnt_to_file, \
    remove_cn_data_from_segments, remove_cnb_data_from_segments, remove_cn_data_from_adjacencies, EXTERNAL_NA_ID
from rck.core.process import positions_aligned, adj_groups_concur
from rck.core.structures import get_segments_for_fragments_ids_dict, get_ref_telomeres_from_segments, get_ref_adjacencies_from_segments, \
    SegmentCopyNumberBoundaries, AdjacencyType, AdjacencyGroupType, Haplotype, Phasing
from rck.utils.adj.adjacency_group_process import projected_groups, refined_labeling_groups

log = logging.getLogger("rck.ecdna")

SCNT, ADJ, ADJ_IMPUTED, FRAGMENTS, TELOMERES, ADJGROUPS, SCNB, PROPORTIONS = (
    "rck.scnt.tsv", "rck.adj.tsv", "rck.adj.imputed.tsv", "rck.fragments.tsv", "rck.telomeres.tsv", "rck.adjgroups.tsv",
    "rck.scnb.tsv", "clone_proportions.tsv")


@dataclass
class MultiCloneParams:
    """Settings of the multi-clone route. RCK's run defaults (``rck_run.py``) unless noted."""
    nas_fp: float = 0.1  # --run-nas-fp
    group_m_fp: float = 0.1  # --run-group-m-default-fp
    group_n_fp: float = 0.1  # --run-group-n-default-fp
    segment_length_attr: str = "length_10"  # --run-segment-length-attr
    # ecDNA additions, off by default (docs/ecDNA.md, multi-clone deviation table)
    imputed_penalty: Optional[float] = None  # None: rck.adj.imputed.tsv is not read
    novel_penalty: float = 0.0
    junction_weight: float = 0.0
    junction_reads_per_copy: Optional[float] = None  # theta; default scale.json next to the DIAG
    linear_ends: str = "none"  # none | boundaries | all
    linear_end_penalty: float = 0.25
    boundary_min_copies: float = 2.0
    boundary_frac: float = 0.25
    # Gurobi
    time_limit: Optional[float] = None
    threads: int = 2
    mip_gap: float = 1e-4
    log_file: Optional[str] = None
    log_to_console: bool = False


@dataclass
class MultiCloneInput:
    rck_dir: str
    segments: list
    scnt: dict
    scnb: dict
    clone_ids: List[str]
    input_adjacencies: list  # observed novel adjacencies (+ imputed when read)
    adjacencies: list  # input + reference adjacencies derived from the segments
    imputed_aids: List[str]
    telomeres: list
    segments_to_fragments: dict
    adjacency_groups: list
    clone_proportions: Dict[str, float] = field(default_factory=dict)


def _path(rck_dir: str, name: str) -> Optional[str]:
    path = os.path.join(rck_dir, name)
    return path if os.path.exists(path) and os.path.getsize(path) > 0 else None


def read_clone_proportions(path: Optional[str]) -> Dict[str, float]:
    if path is None:
        return {}
    result = {}
    with open(path) as handle:
        header = handle.readline().rstrip("\n").split("\t")
        for line in handle:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            result[row["clone"]] = float(row["u"])
    return result


def load_rck_input(rck_dir: str, read_imputed: bool = False) -> MultiCloneInput:
    """Read a DIAG folder as ``rck_run.py --no-pre`` does (same readers, same checks)."""
    scnt_path = _path(rck_dir, SCNT)
    if scnt_path is None:
        raise SystemExit("--rck-input {d}: no {f}; build it with subgraph_construct.py --clone-profiles".format(d=rck_dir, f=SCNT))
    segments = read_segments_from_file(file_name=scnt_path)
    scnt = extract_scnt_from_segments(segments=segments)
    clone_ids = sorted(scnt.keys())
    try:
        scnb = extract_scnb_from_segments(segments=segments, clone_ids=clone_ids, allow_missing=True)
    except ValueError:
        scnb = {cid: SegmentCopyNumberBoundaries() for cid in clone_ids}
    scnb_path = _path(rck_dir, SCNB)
    if scnb_path is not None:
        _, scnb_from_file = read_scnb_from_file(file_name=scnb_path, clone_ids=clone_ids, allow_missing=True)
        for cid in clone_ids:
            scnb[cid].update(scnb_from_file[cid]._records)
    for cid in clone_ids:  # rck_run.py fills missing bounds with its defaults (uniform-min-max 0..10)
        scnb[cid].fill(segments=segments, scnp=scnt[cid], missing_only=True)

    input_adjacencies = read_adjacencies_from_file(file_name=_path(rck_dir, ADJ)) if _path(rck_dir, ADJ) else []
    imputed_aids = []
    if read_imputed and _path(rck_dir, ADJ_IMPUTED):
        imputed = read_adjacencies_from_file(file_name=_path(rck_dir, ADJ_IMPUTED))
        imputed_aids = [a.stable_id_non_phased for a in imputed]
        input_adjacencies = input_adjacencies + imputed
    telomere_positions = read_positions_from_file(file_name=_path(rck_dir, TELOMERES)) if _path(rck_dir, TELOMERES) else []
    adjacency_groups = read_adjacency_groups_from_file(file_name=_path(rck_dir, ADJGROUPS)) if _path(rck_dir, ADJGROUPS) else []
    fragments = read_segments_from_file(file_name=_path(rck_dir, FRAGMENTS)) if _path(rck_dir, FRAGMENTS) else None
    if fragments is None:
        fragments = list(segments)

    # rck_run.py's alignment checks
    segments_to_fragments = get_segments_for_fragments_ids_dict(segments=segments, fragments=fragments, allow_non_covered=False)
    segments_positions = [p for s in segments for p in (s.start_position, s.end_position)]
    if not positions_aligned(segments_positions=segments_positions,
                             other_positions=[p for a in input_adjacencies for p in (a.position1, a.position2)]):
        raise ValueError("novel adjacency positions are not aligned with segment extremities in {d}".format(d=rck_dir))
    telomeres = sorted(get_ref_telomeres_from_segments(segments=segments) + telomere_positions, key=lambda p: p.coordinate)
    if not positions_aligned(segments_positions=segments_positions, other_positions=telomeres):
        raise ValueError("telomere positions are not aligned with segment extremities in {d}".format(d=rck_dir))
    if adjacency_groups and not adj_groups_concur(adj_groups=adjacency_groups, adjacencies=input_adjacencies):
        raise ValueError("an adjacency group refers to an adjacency missing from {f}".format(f=ADJ))

    ref_adjacencies = get_ref_adjacencies_from_segments(segments=segments, assign_external_ids=True)
    adjacencies = input_adjacencies + ref_adjacencies
    by_type = {t: projected_groups(groups=[ag for ag in adjacency_groups if ag.group_type == t], adjacencies=adjacencies)
               for t in (AdjacencyGroupType.MOLECULE, AdjacencyGroupType.LABELING, AdjacencyGroupType.GENERAL)}
    if by_type[AdjacencyGroupType.LABELING]:
        from rck.core.graph import IntervalAdjacencyGraph
        iag = IntervalAdjacencyGraph(segments=segments, adjacencies=adjacencies)
        iag.build_graph()
        by_type[AdjacencyGroupType.LABELING] = refined_labeling_groups(adj_groups=by_type[AdjacencyGroupType.LABELING], iag=iag, adjacencies=adjacencies,
                                                                       gid_suffix="rck-run-L", retain_source_gids=True)
    adjacency_groups = by_type[AdjacencyGroupType.MOLECULE] + by_type[AdjacencyGroupType.LABELING] + by_type[AdjacencyGroupType.GENERAL]
    by_external = {a.extra.get(EXTERNAL_NA_ID, a.stable_id_non_phased): a for a in input_adjacencies}
    for ag in adjacency_groups:
        ag.populate_adjacencies_via_ids(source=input_adjacencies, source_by_ids=by_external)
    return MultiCloneInput(rck_dir=rck_dir, segments=segments, scnt=scnt, scnb=scnb, clone_ids=clone_ids,
                           input_adjacencies=input_adjacencies, adjacencies=adjacencies, imputed_aids=imputed_aids,
                           telomeres=telomeres, segments_to_fragments=segments_to_fragments, adjacency_groups=adjacency_groups,
                           clone_proportions=read_clone_proportions(_path(rck_dir, PROPORTIONS)))


def linear_end_candidates(inp: MultiCloneInput, mode: str, min_copies: float = 2.0, frac: float = 0.25) -> list:
    """Extra extremities where copies may end at a cost (``--mc-linear-ends``): ``all``
    segment extremities, or ``boundaries``, those whose total copy number (any clone) steps by at
    least max(min_copies, frac x the larger) to the reference neighbour. RCK's telomeres are
    already free ends and are skipped by the model."""
    mode = (mode or "none").lower()
    if mode == "none":
        return []
    ordered = sorted(inp.segments, key=lambda s: (s.chromosome, s.start_coordinate))
    if mode == "all":
        return [p for s in ordered for p in (s.start_position, s.end_position)]
    if mode != "boundaries":
        raise ValueError("--mc-linear-ends must be none, boundaries or all, not {m}".format(m=mode))

    def total(cid, s):
        return sum(inp.scnt[cid].get_cn(sid=s.stable_id_non_hap, haplotype=h, default=0) for h in (Haplotype.A, Haplotype.B))

    result = []
    for left, right in zip(ordered, ordered[1:]):
        if left.chromosome != right.chromosome or right.start_coordinate != left.end_coordinate + 1:
            continue
        for cid in inp.clone_ids:
            a, b = total(cid, left), total(cid, right)
            if abs(a - b) >= max(min_copies, frac * max(a, b)):
                result.extend([left.end_position, right.start_position])
                break
    return result


def junction_copies(inp: MultiCloneInput, theta: Optional[float]) -> Dict[str, float]:
    """#NOV / theta (bulk copies) for every observed novel adjacency with read support."""
    if not theta:
        return {}
    result = {}
    imputed = set(inp.imputed_aids)
    for a in inp.input_adjacencies:
        if a.adjacency_type != AdjacencyType.NOVEL or a.stable_id_non_phased in imputed or "#NOV" not in a.extra:
            continue
        try:
            result[a.stable_id_non_phased] = float(a.extra["#NOV"]) / float(theta)
        except (TypeError, ValueError):
            continue
    return result


def build_model(inp: MultiCloneInput, params: MultiCloneParams, theta: Optional[float] = None):
    from rck.core.ilp_gurobi import DEFAULT_GROUP_M_FP, DEFAULT_GROUP_N_FP, SEGMENT_LENGTH_ATTRIBUTE
    from rck.core.ilp_gurobi_ecdna import EcDNAOptions, OptModelMultiCloneEcDNA

    if params.junction_weight and not theta:
        raise SystemExit("--mc-junction-weight needs junction reads per copy: --junction-reads-per-copy or scale.json next to the DIAG")
    options = EcDNAOptions(imputed_aids=inp.imputed_aids, imputed_penalty=params.imputed_penalty or 0.0,
                           novel_penalty=params.novel_penalty, junction_weight=params.junction_weight,
                           junction_copies=junction_copies(inp, theta) if params.junction_weight else {},
                           clone_proportions=inp.clone_proportions,
                           linear_end_nodes=linear_end_candidates(inp, params.linear_ends, params.boundary_min_copies, params.boundary_frac),
                           linear_end_penalty=params.linear_end_penalty)
    if params.junction_weight and not inp.clone_proportions:
        raise SystemExit("--mc-junction-weight needs {f} in the DIAG folder".format(f=PROPORTIONS))
    extra = {DEFAULT_GROUP_M_FP: params.group_m_fp, DEFAULT_GROUP_N_FP: params.group_n_fp, SEGMENT_LENGTH_ATTRIBUTE: params.segment_length_attr}
    model = OptModelMultiCloneEcDNA(hapl_segments=inp.segments, hapl_adjacencies=inp.adjacencies, scnt=inp.scnt, hapl_telomeres=inp.telomeres,
                                    hapl_segments_to_fragments=inp.segments_to_fragments, hapl_adjacencies_groups=inp.adjacency_groups,
                                    scnb=inp.scnb, hapl_nov_adjacencies_fp=params.nas_fp, extra=extra, ecdna=options)
    model.build_gurobi_model()
    gm = model.gm
    gm.Params.LogToConsole = 1 if params.log_to_console else 0
    if params.log_file:
        gm.Params.LogFile = params.log_file
    gm.Params.Threads = params.threads
    gm.Params.MIPGap = params.mip_gap
    if params.time_limit:
        gm.Params.TimeLimit = params.time_limit
    return model


def solve(inp: MultiCloneInput, params: MultiCloneParams, theta: Optional[float] = None):
    import gurobipy as g

    model = build_model(inp, params, theta)
    model.solve_model()
    status = model.gm.Status
    if status == g.GRB.INFEASIBLE:
        raise RuntimeError("multi-clone model is infeasible (check the copy-number bounds and --mc-nas-fp)")
    if model.gm.SolCount == 0:
        raise RuntimeError("Gurobi found no solution (status {s})".format(s=status))
    return model


def write_outputs(outdir: str, inp: MultiCloneInput, model, extra_summary: Optional[dict] = None) -> Dict[str, str]:
    """RCK's native per-clone ``rck.scnt.tsv`` and ``rck.acnt.tsv``, written as ``rck_run.py``
    writes them, and ``summary.json`` with ``model: multiclone``."""
    os.makedirs(outdir, exist_ok=True)
    scnt = model.get_scnt_from_model()
    acnt = model.get_acnt_from_model()
    remove_cn_data_from_segments(segments=inp.segments)
    remove_cnb_data_from_segments(segments=inp.segments)
    scnt_path = os.path.join(outdir, "rck.scnt.tsv")
    write_scnt_to_file(file_name=scnt_path, segments=inp.segments, scnt=scnt)
    remove_cn_data_from_adjacencies(adjacencies=inp.adjacencies)
    acnt_path = os.path.join(outdir, "rck.acnt.tsv")
    write_acnt_to_file(file_name=acnt_path, acnt=acnt, adjacencies=inp.adjacencies, output_reference=True)

    gm = model.gm
    flips = {str(fid): model.variables["fragment_flipping"][fid].X > 0.5 for fid in model.fragments_ids}
    summary = {
        "model": "multiclone",
        "rck_input": os.path.abspath(inp.rck_dir),
        "clones": inp.clone_ids,
        "clone_proportions": inp.clone_proportions,
        "status": int(gm.Status),
        "objective": gm.ObjVal,
        "objective_bound": gm.ObjBound,
        "mip_gap": gm.MIPGap if gm.IsMIP else 0.0,
        "segments": len(inp.segments),
        "novel_adjacencies": sum(1 for a in inp.input_adjacencies if a.adjacency_type == AdjacencyType.NOVEL),
        "imputed_adjacencies": len(inp.imputed_aids),
        "telomeres": len(inp.telomeres),
        "fragments": len(model.fragments_ids),
        "fragments_flipped": sorted(fid for fid, kept in flips.items() if not kept),
        "novel_used": novel_usage(inp, acnt),
        "linear_ends": {"{c}:{v}:{h}".format(c=c, v=v, h=h.value): n for (c, v, h), n in model.get_linear_ends_from_model().items()},
    }
    if extra_summary:
        summary.update(extra_summary)
    summary_path = os.path.join(outdir, "summary.json")
    with open(summary_path, "w") as handle:
        json.dump(summary, handle, indent=2, sort_keys=True, default=str)
    return {"scnt": scnt_path, "acnt": acnt_path, "summary": summary_path}


def novel_usage(inp: MultiCloneInput, acnt) -> List[dict]:
    rows = []
    imputed = set(inp.imputed_aids)
    for a in inp.input_adjacencies:
        aid = a.stable_id_non_phased
        for cid in inp.clone_ids:
            for ph in (Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB):
                cn = acnt[cid].get_cn(aid=aid, phasing=ph, default=0)
                if cn:
                    rows.append({"aid": a.extra.get(EXTERNAL_NA_ID, aid), "clone": cid, "phasing": str(ph), "cn": cn, "imputed": aid in imputed})
    return rows


def params_from_args(args) -> MultiCloneParams:
    return MultiCloneParams(nas_fp=args.mc_nas_fp, group_m_fp=args.mc_group_m_fp, segment_length_attr=args.mc_segment_length_attr,
                            imputed_penalty=args.mc_imputed_penalty, novel_penalty=args.mc_novel_penalty,
                            junction_weight=args.mc_junction_weight, junction_reads_per_copy=args.junction_reads_per_copy,
                            linear_ends=args.mc_linear_ends, linear_end_penalty=args.mc_linear_end_penalty,
                            boundary_min_copies=args.boundary_min_copies, boundary_frac=args.boundary_frac,
                            time_limit=args.time_limit, threads=args.threads, mip_gap=args.mip_gap,
                            log_file=os.path.join(args.outdir, "gurobi_multiclone.log"), log_to_console=args.verbose)


def theta_for(rck_dir: str, override: Optional[float]) -> Optional[float]:
    """Junction reads per copy: the override, else scale.json of the subgraph holding the DIAG."""
    if override:
        return float(override)
    path = os.path.join(os.path.dirname(os.path.abspath(rck_dir)), "scale.json")
    if os.path.exists(path):
        with open(path) as handle:
            value = json.load(handle).get("junction_reads_per_copy")
        return float(value) if value else None
    return None


def run_multiclone(args) -> Dict[str, object]:
    from rck.ecdna.gurobi_env import close_model

    started = time.time()
    os.makedirs(args.outdir, exist_ok=True)
    params = params_from_args(args)
    inp = load_rck_input(args.rck_input, read_imputed=params.imputed_penalty is not None)
    theta = theta_for(args.rck_input, params.junction_reads_per_copy) if params.junction_weight else None
    log.info("multiclone: %d segments, %d novel adjacencies (%d imputed), %d telomeres, %d fragments, clones %s; Gurobi threads %d",
             len(inp.segments), len(inp.input_adjacencies), len(inp.imputed_aids), len(inp.telomeres),
             len(set(inp.segments_to_fragments.values())), ",".join(inp.clone_ids), params.threads)
    model = solve(inp, params, theta)
    try:
        outputs = write_outputs(args.outdir, inp, model,
                                extra_summary={"runtime_s": round(time.time() - started, 2), "junction_reads_per_copy": theta,
                                               "params": {k: v for k, v in vars(params).items() if k not in ("log_file", "log_to_console")}})
        with open(outputs["summary"]) as handle:
            summary = json.load(handle)
        log.info("multiclone: status %s, objective %.4f, %d fragments flipped", summary["status"], summary["objective"],
                 len(summary["fragments_flipped"]))
    finally:
        close_model(model.gm)
    return {"input": inp, "outputs": outputs, "summary": summary}

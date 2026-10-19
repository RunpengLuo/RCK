"""Write flow, decomposition and AA-format outputs."""

from __future__ import annotations

import json
import os
from typing import Dict, List, Optional

import pandas as pd

from rck.ecdna.decompose import Decomposition
from rck.ecdna.ilp_circular import FlowSolution
from rck.ecdna.scale import Observations
from rck.ecdna.subgraph_io import HAPLOTYPES, REFERENCE, ECGraph
from rck.ecdna.walk import describe_join, describe_step

DECOMPOSITION_COLUMNS = ["decomposition_rank", "objective", "residual", "n_species", "species", "status", "weight", "length_bp",
                         "ordering_rank", "n_orderings", "path_support", "paths_matched", "paths_applicable",
                         "n_steps", "uses_imputed", "junctions", "walk", "joins"]


def decomposition_table(graph: ECGraph, decompositions: List[Decomposition]) -> pd.DataFrame:
    """One row per (decomposition, species, ordering).

    ``status`` is ``species`` or ``low_weight`` (below the minimum species weight; usually flow
    residue rather than a molecule). ``walk`` lists the steps as ``chrom:start-end:hap:strand``; ``joins`` lists, for each step,
    the adjacency taken after it as ``kind:aid:phasing`` (the last one closes the circle).
    """
    rows = []
    for d in decompositions:
        for i, sp in enumerate(d.species, start=1):
            for j, o in enumerate(sp.orderings, start=1):
                rows.append({
                    "decomposition_rank": d.rank, "objective": d.objective, "residual": d.residual, "n_species": len(d.species),
                    "species": "ec{i}".format(i=i), "status": sp.status, "weight": sp.weight, "length_bp": sp.length_bp, "ordering_rank": j,
                    "n_orderings": len(sp.orderings), "path_support": o.path_score, "paths_matched": o.paths_matched,
                    "paths_applicable": o.paths_applicable, "n_steps": len(o.steps), "uses_imputed": ",".join(sp.imputed) or ".",
                    "junctions": ",".join(sp.junctions) or ".",
                    "walk": ",".join(describe_step(graph, s) for s in o.steps),
                    "joins": ",".join(describe_join(graph, s) for s in o.steps),
                })
    return pd.DataFrame(rows, columns=DECOMPOSITION_COLUMNS)


def _side_coordinate(graph: ECGraph, nid: str):
    node = graph.nodes[nid]
    seg = graph.segment_at(nid)
    if seg is None:
        return node.chrom, node.pos, node.side or "+"
    if nid == seg.tail:
        return seg.chrom, seg.start, "-"
    return seg.chrom, seg.end, "+"


def write_aa_graph(path: str, graph: ECGraph, obs: Observations, sol: FlowSolution) -> None:
    """AmpliconArchitect-style breakpoint graph (total copy number; haplotypes summed)."""
    lines = ["SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped"]
    for sid, s in sorted(graph.segments.items(), key=lambda kv: (kv[1].chrom, kv[1].start)):
        cn = sum(sol.seg_ec[(sid, h)] + sol.seg_bg[(sid, h)] for h in HAPLOTYPES)
        lines.append("sequence\t{c}:{s}-\t{c}:{e}+\t{cn:.6f}\t{cov:.6f}\t{size}\t0".format(
            c=s.chrom, s=s.start, e=s.end, cn=float(cn), cov=obs.segments[sid].depth, size=s.length))
    lines.append("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, "
                 "HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence")
    for aid, a in sorted(graph.adjacencies.items()):
        cn = sum(v for (x, _), v in sol.adj_ec.items() if x == aid) + sum(v for (x, _), v in sol.adj_bg.items() if x == aid)
        if cn <= 0 and a.kind != REFERENCE:
            continue
        c1, p1, s1 = _side_coordinate(graph, a.u)
        c2, p2, s2 = _side_coordinate(graph, a.v)
        kind = "concordant" if a.kind == REFERENCE else "discordant"
        lines.append("{k}\t{c1}:{p1}{s1}->{c2}:{p2}{s2}\t{cn:.6f}\t{r}\tNone\tNone".format(
            k=kind, c1=c1, p1=p1, s1=s1, c2=c2, p2=p2, s2=s2, cn=float(cn), r=int(a.reads or 0)))
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def write_aa_cycles(path: str, graph: ECGraph, decomposition: Decomposition) -> None:
    """AmpliconArchitect-style cycles file for one decomposition (top ordering per species)."""
    order = sorted(graph.segments, key=lambda sid: (graph.segments[sid].chrom, graph.segments[sid].start))
    number = {sid: i for i, sid in enumerate(order, start=1)}
    lines = []
    for i, (c, s, e) in enumerate(graph.regions(), start=1):
        lines.append("Interval\t{i}\t{c}\t{s}\t{e}".format(i=i, c=c, s=s, e=e))
    lines.append("List of cycle segments")
    for sid in order:
        seg = graph.segments[sid]
        lines.append("Segment\t{n}\t{c}\t{s}\t{e}".format(n=number[sid], c=seg.chrom, s=seg.start, e=seg.end))
    for i, sp in enumerate(decomposition.species, start=1):
        steps = sp.orderings[0].steps
        segs = ",".join("{n}{st}".format(n=number[s.sid], st=s.strand) for s in steps)
        haps = ",".join(s.hap for s in steps)
        lines.append("Cycle={i};Copy_count={w};Segments={segs};Haplotypes={h}".format(i=i, w=float(sp.weight), segs=segs, h=haps))
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def write_all(outdir: str, graph: ECGraph, obs: Observations, sol: FlowSolution, decompositions: List[Decomposition],
              extra: Optional[Dict[str, object]] = None) -> Dict[str, str]:
    os.makedirs(outdir, exist_ok=True)
    aa_dir = os.path.join(outdir, "aa")
    os.makedirs(aa_dir, exist_ok=True)
    paths = {
        "scale": os.path.join(outdir, "observations.json"),
        "flow_segments": os.path.join(outdir, "flow_segments.tsv"),
        "flow_adjacencies": os.path.join(outdir, "flow_adjacencies.tsv"),
        "decompositions": os.path.join(outdir, "decompositions.tsv"),
        "aa_graph": os.path.join(aa_dir, "{n}_graph.txt".format(n=graph.name)),
        "summary": os.path.join(outdir, "summary.json"),
    }
    obs.write_json(paths["scale"])
    sol.segment_table(graph, obs).to_csv(paths["flow_segments"], sep="\t", index=False)
    sol.adjacency_table(graph, obs).to_csv(paths["flow_adjacencies"], sep="\t", index=False)
    decomposition_table(graph, decompositions).to_csv(paths["decompositions"], sep="\t", index=False)
    write_aa_graph(paths["aa_graph"], graph, obs, sol)
    for d in decompositions:
        write_aa_cycles(os.path.join(aa_dir, "{n}_decomp{r}_cycles.txt".format(n=graph.name, r=d.rank)), graph, d)
    summary = {
        "graph": graph.summary(),
        "scale": obs.scale.to_dict(),
        "flow": {"status": sol.status, "objective": sol.objective, "mip_gap": sol.mip_gap, "imputed_used": sol.imputed_used,
                 "ecdna_copies_total": int(sum(sol.ecdna_edges().values()))},
        "decompositions": [
            {"rank": d.rank, "objective": d.objective, "residual": d.residual,
             "species": [{"weight": s.weight, "status": s.status, "length_bp": s.length_bp, "n_orderings": len(s.orderings),
                          "imputed": s.imputed, "junctions": s.junctions} for s in d.species]}
            for d in decompositions
        ],
    }
    if extra:
        summary.update(extra)
    with open(paths["summary"], "w") as handle:
        json.dump(summary, handle, indent=2, default=str)
    return paths

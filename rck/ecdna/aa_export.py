"""Export an ``rck-ecdna`` answer, or a subgraph, to AA format.

    python -m rck.ecdna.aa_export RUN_DIR OUT [--pool-index I | --lambda L] [--alternatives 2]
    python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name N] [--amplicon 1]
                                  [--imputed comment|drop] [--min-path-junctions 1]

The ``subgraph`` subcommand (:func:`export_subgraph`) writes the breakpoint graph RCK-ecDNA
solves on as an AA graph that CoRAL 3.1.1 ``cycle`` and CycleExtractor read, for the graph
route of the benchmark (results/graph_route/); see its docstring. The rest of this docstring is
about exporting an answer.

``rck-ecdna`` already writes ``aa/<name>_graph.txt`` and ``aa/<name>_decomp<k>_cycles.txt``,
but those files lack what the scorer needs to reproduce ``rck-ecdna-eval`` exactly: the
``low_weight`` status, a species' other orderings, and the solutions of a kept pool
(``--keep-pool``, which the λ > 0 results come from). This module writes, from
``decompositions.tsv`` (or ``pool/decompositions.tsv``):

* ``<name>_cycles.txt``: the answer, one cycle per species with ``Copy_count`` = the refitted
  weight, ``Haplotypes=``, ``Status=``, ``Species=`` and the other orderings
  (``Alt_segments=`` / ``Alt_haplotypes=``); ``Interval`` lines are the run's regions;
* ``<name>_alt<k>_cycles.txt``: the first k same-fit alternatives of the answer (the set
  metric of results/bench_nested/diag/grow_rules.py), for best-of scoring;
* ``<name>_graph.txt``: the run's AA graph;
* ``<name>_segment_cn.tsv``: allele copy numbers per segment (ecDNA + background, from
  ``flow_segments.tsv``).

Coordinates: the scorer reads AA segments as ``start - 1`` and ``end - 1`` (0-based inclusive,
how the simulator's 1-based AA truth is read). ``rck-ecdna``'s walks and ``aa/`` files carry the
subgraph's segment coordinates unchanged, so every position is written here as the walk's
coordinate + 1: the scorer then reads back exactly the steps ``rck-ecdna-eval`` reads from
``decompositions.tsv``, and its scores are identical.
"""

from __future__ import annotations

import argparse
import glob
import json
import os
import sys
from collections import Counter
from typing import Dict, List, Optional, Sequence, Tuple

import pandas as pd

SHIFT = 1  # walk coordinate -> AA coordinate (see the module docstring)


def _token(text: str):
    chrom, span, hap, strand = text.rsplit(":", 3)
    start, end = span.split("-")
    return chrom, int(start), int(end), hap, strand


def write_cycles(path: str, rows: pd.DataFrame, regions: Sequence[str]) -> None:
    """One decomposition's species (``rows``: its decompositions.tsv rows, every ordering)."""
    walks = {}
    for (species, _), group in rows.groupby(["species", "ordering_rank"], sort=False):
        walks.setdefault(species, []).append([_token(t) for t in str(group["walk"].iloc[0]).split(",")])
    segs = sorted({(c, s, e) for ws in walks.values() for w in ws for c, s, e, _, __ in w})
    number = {seg: i for i, seg in enumerate(segs, start=1)}
    lines = ["# rck-ecdna answer for the method-agnostic scorer (rck.ecdna.aa_export); coordinates = walk + {s}".format(s=SHIFT)]
    for i, text in enumerate(regions, start=1):
        chrom, span = text.rsplit(":", 1)
        a, b = span.split("-")
        lines.append("Interval\t{i}\t{c}\t{a}\t{b}".format(i=i, c=chrom, a=int(a) + SHIFT, b=int(b) + SHIFT))
    lines.append("List of cycle segments")
    for seg in segs:
        lines.append("Segment\t{n}\t{c}\t{s}\t{e}".format(n=number[seg], c=seg[0], s=seg[1] + SHIFT, e=seg[2] + SHIFT))
    first = rows.drop_duplicates("species")
    for k, r in enumerate(first.to_dict("records"), start=1):
        linear = str(r.get("topology", "circular")) == "linear"

        def encode(walk):
            segtext = ",".join("{n}{st}".format(n=number[(c, s, e)], st=st) for c, s, e, _, st in walk)
            haps = ",".join(h for _, __, ___, h, ____ in walk)
            return ("0+," + segtext + ",0-", ".," + haps + ",.") if linear else (segtext, haps)

        main_segs, main_haps = encode(walks[r["species"]][0])
        fields = ["Cycle={k}".format(k=k), "Copy_count={w!r}".format(w=float(r["weight"])), "Segments=" + main_segs,
                  "Haplotypes=" + main_haps, "Status=" + str(r.get("status", "species")), "Species=" + str(r["species"])]
        alts = [encode(w) for w in walks[r["species"]][1:]]
        if alts:
            fields += ["Alt_segments=" + "|".join(a for a, _ in alts), "Alt_haplotypes=" + "|".join(h for _, h in alts)]
        lines.append(";".join(fields))
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def shift_graph(src: str, dst: str) -> None:
    """Copy an ``rck-ecdna`` AA graph with every position + SHIFT."""
    from rck.ecdna.aa_format import parse_breakend

    def be(text):
        c, p, s = parse_breakend(text)
        return "{c}:{p}{s}".format(c=c, p=p + SHIFT, s=s)

    out = []
    with open(src) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] == "sequence":
                f[1], f[2] = be(f[1]), be(f[2])
            elif f[0] in ("concordant", "discordant", "source"):
                a, b = f[1].split("->")
                f[1] = be(a) + "->" + be(b)
            out.append("\t".join(f))
    with open(dst, "w") as handle:
        handle.write("\n".join(out) + "\n")


def write_segment_cn(src: str, dst: str) -> None:
    t = pd.read_csv(src, sep="\t")
    a = t["ecdna_A"] + t["bg_A"]
    b = t["ecdna_B"] + t["bg_B"]
    pd.DataFrame({"chrom": t["chrom"], "start": t["start"] + SHIFT, "end": t["end"] + SHIFT, "cn": a + b, "cn_a": a,
                  "cn_b": b}).to_csv(dst, sep="\t", index=False)


def run_alternatives(run_dir: str, k: int) -> List[int]:
    """Decomposition ranks of rank 1's first k same-fit alternatives in a run (summary.json
    ``ambiguity``); alternatives the run did not report (rank null) are skipped, as unmatched
    ones are in grow_rules.py's set metric."""
    with open(os.path.join(run_dir, "summary.json")) as handle:
        amb = json.load(handle).get("ambiguity") or {}
    return [a["rank"] for a in (amb.get("alternatives") or [])[:k] if a.get("rank") is not None]


def pool_alternatives(pool: Dict[str, object], top: int, k: int) -> List[int]:
    """Pool indices of pool entry ``top``'s first k same-fit alternatives (pool.json
    ``same_fit``), matched to entries by objective, steps and significant weights as
    results/bench_nested/diag/grow_rules.py's set metric does."""
    entries = pool["entries"]
    sig = lambda e: sorted(w for w, x in zip(e["weights_int"], e["low_weight"]) if not x)  # noqa: E731
    out = []
    for alt in (pool["same_fit"].get(str(top)) or [])[:k]:
        hit = [e["index"] for e in entries if abs(e["objective"] - alt["objective"]) < 1e-3
               and e["steps"] == alt["steps"] and sig(e) == sorted(alt["weights"])]
        if hit:
            out.append(hit[0])
    return out


def export(run_dir: str, outdir: str, pool_index: Optional[int] = None, alternatives: int = 2, name: Optional[str] = None) -> Dict[str, object]:
    """Write the answer of ``run_dir`` (rank 1; or pool entry ``pool_index`` of its kept pool)
    and its alternatives to ``outdir``. Returns the paths for :func:`rck.ecdna.aa_eval.score`."""
    os.makedirs(outdir, exist_ok=True)
    with open(os.path.join(run_dir, "summary.json")) as handle:
        summary = json.load(handle)
    regions = summary["graph"]["regions"]
    name = name or summary["graph"]["name"]
    if pool_index is None:
        table = pd.read_csv(os.path.join(run_dir, "decompositions.tsv"), sep="\t")
        ranks = [int(table["decomposition_rank"].min())] if len(table) else []
        ranks += run_alternatives(run_dir, alternatives)
    else:
        table = pd.read_csv(os.path.join(run_dir, "pool", "decompositions.tsv"), sep="\t")
        with open(os.path.join(run_dir, "pool", "pool.json")) as handle:
            pool = json.load(handle)
        ranks = [pool_index + 1] + [i + 1 for i in pool_alternatives(pool, pool_index, alternatives)]
    paths = {"cycles": os.path.join(outdir, name + "_cycles.txt"), "alt_cycles": [], "ranks": ranks}
    for k, rank in enumerate(ranks):
        path = paths["cycles"] if k == 0 else os.path.join(outdir, "{n}_alt{k}_cycles.txt".format(n=name, k=k))
        write_cycles(path, table[table["decomposition_rank"] == rank], regions)
        if k:
            paths["alt_cycles"].append(path)
    graphs = glob.glob(os.path.join(run_dir, "aa", "*_graph.txt"))
    if graphs:
        paths["graph"] = os.path.join(outdir, name + "_graph.txt")
        shift_graph(graphs[0], paths["graph"])
    if os.path.exists(os.path.join(run_dir, "flow_segments.tsv")):
        paths["segment_cn"] = os.path.join(outdir, name + "_segment_cn.tsv")
        write_segment_cn(os.path.join(run_dir, "flow_segments.tsv"), paths["segment_cn"])
    return paths


# ------------------------------------------------------------------------------------------ #
# subgraph -> AA graph (graph route)
# ------------------------------------------------------------------------------------------ #

GRAPH_FORMAT = "rck-ecdna-aa-graph/1"


def _na(x) -> str:
    return "NA" if x is None else repr(float(x))


def subgraph_observations(subgraph_dir: str, graph, depth_per_copy: Optional[float] = None, junction_reads_per_copy: Optional[float] = None,
                          baseline_cn: Optional[float] = None):
    """The observations ``rck-ecdna`` builds for ``graph`` (``--inputs auto``, no windows): from
    ``scale.json`` / ``cn_profile.tsv`` when the folder has them, else from the overrides."""
    from rck.ecdna.scale import ScaleParams, build_observations, load_subgraph_scale

    sub = load_subgraph_scale(subgraph_dir)
    if baseline_cn is None:
        baseline_cn = float(sub.get("background_cn", 2.0)) if sub is not None else 2.0
    params = ScaleParams(depth_per_copy=depth_per_copy, junction_reads_per_copy=junction_reads_per_copy, baseline_cn=baseline_cn)
    return build_observations(graph, params, subgraph_scale=sub), sub


def _runs(sequences) -> List[Tuple[str, int, int]]:
    """Runs of contiguous sequence edges (next start == previous end + 1), CoRAL 3.1.1's
    ``infer_amplicon_intervals``; ``sequences`` sorted by (chrom, start)."""
    out: List[Tuple[str, int, int]] = []
    for s in sequences:
        if out and out[-1][0] == s.chrom and s.start == out[-1][2] + 1:
            out[-1] = (s.chrom, out[-1][1], s.end)
        else:
            out.append((s.chrom, s.start, s.end))
    return out


def export_subgraph(subgraph_dir: str, outdir: str, name: Optional[str] = None, amplicon: int = 1, imputed: str = "comment",
                    min_path_junctions: int = 1, depth_per_copy: Optional[float] = None, junction_reads_per_copy: Optional[float] = None,
                    baseline_cn: Optional[float] = None, phantom_max_bp: int = 2) -> Dict[str, object]:
    """Write ``OUTDIR/<name>_amplicon<N>_graph.txt`` (+ ``_segment_cn.tsv``, ``_export.json``).

    The graph is the one ``rck-ecdna`` solves on: :func:`rck.ecdna.subgraph_io.read_subgraph`
    (phantom segments contracted, parallel imputed edges merged), with the copy numbers it uses
    (:func:`rck.ecdna.scale.build_observations`): segment total CN, and junction CN = reads /
    junction reads per copy. Positions are subgraph coordinates + ``SHIFT`` (as for answers).

    * ``sequence``: one per segment, sorted; AverageCoverage = the segment's depth (WEIGHT).
    * ``concordant``: every reference adjacency between two segments; CN = observed CN, or for
      an unobserved join (at an imputed grid point) the smaller segment CN (``export.json``).
    * ``discordant``: every observed junction between two segments; fold-backs (both ends on
      one node) are written as they are and listed in ``export.json``.
    * imputed candidates (RCK-ecDNA's own, no read evidence) are not CoRAL / CycleExtractor
      edges: ``imputed='comment'`` writes them as ``#rck_imputed`` lines, ``'drop'`` leaves them out.
    * adjacencies with an end on a node that has no segment (the reference join at a region's
      outer flank) cannot be AA edges (CoRAL adds its own s/t ends there): ``#rck_outside``.
    * ``path_constraint``: the read paths of ``paths.tsv`` that only use segments, reference
      joins and observed junctions and have at least ``min_path_junctions`` junctions, as CoRAL
      3.1.1 tokens (:func:`rck.ecdna.aa_format.coral_path_string`); support = ``#NOV`` rounded.
      Paths ``rck-ecdna`` reads that cannot be written so (they use an imputed edge, or a join
      the graph lacks) become ``#rck_path`` lines.
    * ``interval``: the runs of contiguous segments, which is what CoRAL infers; a candidate
      region with a gap between segments is an error.
    * ``#rck_node`` / ``#rck_seq`` / ``#rck_edge`` / ``#rck_pc`` lines keep the subgraph's node,
      segment, adjacency and path ids and orders, so :mod:`rck.ecdna.aa_import` restores the
      graph exactly. CoRAL and CycleExtractor skip every ``#`` line.
    """
    from rck.ecdna import aa_format as aaf
    from rck.ecdna.subgraph_io import IMPUTED, NOVEL, REFERENCE, read_subgraph

    if imputed not in ("comment", "drop"):
        raise ValueError("imputed must be comment or drop")
    os.makedirs(outdir, exist_ok=True)
    name = name or os.path.basename(os.path.normpath(subgraph_dir))
    graph = read_subgraph(subgraph_dir, phantom_max_bp=phantom_max_bp, name=name)
    obs, sub = subgraph_observations(subgraph_dir, graph, depth_per_copy, junction_reads_per_copy, baseline_cn)
    theta = obs.scale.junction_reads_per_copy

    def coord(nid):
        node = graph.nodes[nid]
        seg = graph.segment_at(nid)
        if seg is None:
            return node.chrom, node.pos + SHIFT, node.side or "+"
        return (seg.chrom, seg.start + SHIFT, "-") if nid == seg.tail else (seg.chrom, seg.end + SHIFT, "+")

    order = sorted(graph.segments, key=lambda k: (aaf.chrom_key(graph.segments[k].chrom), graph.segments[k].start))
    seq_index = {sid: i for i, sid in enumerate(order)}
    sequences = []
    for sid in order:
        s = graph.segments[sid]
        sequences.append(aaf.AASequence(s.chrom, s.start + SHIFT, s.end + SHIFT, obs.segments[sid].total, coverage=s.depth,
                                        size=s.length, reads=0))
    runs = _runs(sequences)
    regions = [(c, a + SHIFT, b + SHIFT) for c, a, b in graph.regions()]
    if sorted(runs) != sorted(regions):
        raise ValueError("{n}: candidate regions {r} are not runs of contiguous segments {u}; CoRAL 3.1.1 would infer other "
                         "amplicon intervals".format(n=name, r=regions, u=runs))

    adj_order = {aid: i for i, aid in enumerate(graph.adjacencies)}
    stats: Dict[str, object] = {"concordant_cn_from_segments": [], "foldbacks": [], "outside": [], "imputed": [],
                                "reference_between_noncontiguous": []}
    conc, disc, outside, imp = [], [], [], []
    for aid, a in graph.adjacencies.items():
        su, sv = graph.segment_at(a.u), graph.segment_at(a.v)
        bu, bv = coord(a.u), coord(a.v)
        if a.kind == IMPUTED:
            imp.append((aid, a, bu, bv))
            continue
        if su is None or sv is None:
            outside.append((aid, a, bu, bv))
            continue
        cn = obs.adjacencies.get(aid)
        if a.kind == REFERENCE:
            if cn is None:
                cn = min(obs.segments[su.sid].total, obs.segments[sv.sid].total)
                stats["concordant_cn_from_segments"].append(aid)
            lo, hi = aaf.lower_first(bu, bv)
            if not (lo[0] == hi[0] and lo[2] == "+" and hi[2] == "-" and hi[1] == lo[1] + 1):
                stats["reference_between_noncontiguous"].append(aid)
                raise ValueError("{n}: reference adjacency {a} joins {x} and {y}, which are not neighbouring segment ends".format(
                    n=name, a=aid, x=aaf.format_breakend(bu), y=aaf.format_breakend(bv)))
            conc.append((aid, a, bu, bv, cn))
        elif a.kind == NOVEL:
            if a.is_self_loop:
                stats["foldbacks"].append(aid)
            disc.append((aid, a, bu, bv, 0.0 if cn is None else cn))

    def oriented(bu, bv):
        return aaf.lower_first(bu, bv)

    conc.sort(key=lambda r: aaf.edge_key(*oriented(r[2], r[3])))  # CoRAL re-sorts concordant edges in this order
    disc.sort(key=lambda r: aaf.edge_key(*oriented(r[2], r[3])))
    edges = [aaf.AAEdge("concordant", *oriented(bu, bv), cn=cn, reads=a.reads) for _, a, bu, bv, cn in conc]
    edges += [aaf.AAEdge("discordant", *oriented(bu, bv), cn=cn, reads=a.reads) for _, a, bu, bv, cn in disc]
    edge_index = {aid: ("c", i) for i, (aid, *_r) in enumerate(conc)}
    edge_index.update({aid: ("d", i) for i, (aid, *_r) in enumerate(disc)})

    # ---- read paths -> path constraints ---- #
    constraints, rck_paths = [], []
    pstats = Counter()
    paths_tsv = os.path.join(subgraph_dir, "paths.tsv")
    if os.path.exists(paths_tsv) and os.path.getsize(paths_tsv) > 0:
        frame = pd.read_csv(paths_tsv, sep="\t")
        for k, row in enumerate(frame.to_dict("records")):
            tokens = [t.strip() for t in str(row["PATH"]).split(",") if t.strip()]
            nodes = [graph.resolve("{e}_{o}".format(e=t.rpartition("@")[0], o=t.rpartition("@")[2])) for t in tokens]
            # rck-ecdna's read_paths resolves only the junction steps: a read may start or end
            # outside the graph. Its outer run is then the one segment it leaves / enters by.
            for end, inner in ((0, 1), (len(nodes) - 1, len(nodes) - 2)):
                if len(nodes) >= 2 and nodes[end] is None and nodes[inner] is not None and graph.segment_at(nodes[inner]) is not None:
                    nodes[end] = graph.segment_at(nodes[inner]).other_end(nodes[inner])
                    pstats["outer_end_outside_graph"] += 1
            if any(n is None for n in nodes) or len(nodes) < 2:
                pstats["outside_graph"] += 1  # rck-ecdna's read_paths drops these too
                continue
            path, reason, n_novel, n_rck = _path_edges(graph, nodes, seq_index, edge_index)
            support = float(row.get("#NOV", 0) or 0)
            if reason is None and n_novel >= max(min_path_junctions, 1):
                flipped = int(path[0][1] > path[-1][1])  # coral_path_string writes it the other way round
                constraints.append((aaf.coral_path_string(path), int(round(support)), k, str(row.get("PID", "")), support, flipped))
                pstats["written"] += 1
            elif n_rck == 0:
                pstats["no_junction"] += 1  # only reference joins: rck-ecdna ignores it as well
            elif reason is None:
                pstats["too_few_junctions"] += 1
                rck_paths.append(("few_junctions", k, row))
            elif reason == "imputed" and imputed == "drop":
                pstats["imputed_dropped"] += 1
            else:
                pstats[reason] += 1
                rck_paths.append((reason, k, row))

    # ---- comments ---- #
    comments = ["#rck_node\t{n}\t{c}\t{p}\t{s}\t{src}".format(n=nid, c=coord(nid)[0], p=coord(nid)[1],
                                                                s=graph.nodes[nid].side or ".", src=graph.nodes[nid].source or ".")
                for nid in graph.nodes]
    comments += ["#rck_seq\t{k}\t{sid}\t{t}\t{h}".format(k=i + 1, sid=sid, t=graph.segments[sid].tail, h=graph.segments[sid].head)
                 for i, sid in enumerate(order)]
    for kind, rows in (("concordant", conc), ("discordant", disc)):
        for i, (aid, a, _bu, _bv, _cn) in enumerate(rows):
            comments.append("#rck_edge\t{k}\t{i}\t{aid}\t{u}\t{v}\t{r}\t{o}".format(k=kind, i=i + 1, aid=aid, u=a.u, v=a.v,
                                                                                    r=_na(a.reads), o=adj_order[aid]))
    for aid, a, bu, bv in outside:
        cn = obs.adjacencies.get(aid)
        comments.append("#rck_outside\t{bp}\t{cn}\t{r}\t{kind}\t{aid}\t{u}\t{v}\t{o}".format(
            bp=aaf.format_breakend(bu) + "->" + aaf.format_breakend(bv), cn=_na(cn), r=_na(a.reads), kind=a.kind, aid=aid, u=a.u, v=a.v,
            o=adj_order[aid]))
        stats["outside"].append(aid)
    if imputed == "comment":
        for aid, a, bu, bv in imp:
            comments.append("#rck_imputed\t{bp}\t{w}\t{aid}\t{u}\t{v}\t{o}".format(
                bp=aaf.format_breakend(bu) + "->" + aaf.format_breakend(bv), w=_na(a.imputed_support), aid=aid, u=a.u, v=a.v,
                o=adj_order[aid]))
    stats["imputed"] = [aid for aid, *_r in imp]
    for i, (_t, _s, k, pid, support, flipped) in enumerate(constraints):
        comments.append("#rck_pc\t{i}\t{k}\t{pid}\t{s}\t{f}".format(i=i + 1, k=k, pid=pid, s=repr(support), f=flipped))
    for reason, k, row in rck_paths:
        comments.append("#rck_path\t{r}\t{k}\t{pid}\t{p}\t{n}\t{s}".format(r=reason, k=k, pid=row.get("PID", ""), p=row["PATH"],
                                                                          n=row.get("NSTEP", ""), s=repr(float(row.get("#NOV", 0) or 0))))

    stem = "{n}_amplicon{a}".format(n=name, a=amplicon)
    paths = {"graph": os.path.join(outdir, stem + "_graph.txt"), "segment_cn": os.path.join(outdir, stem + "_segment_cn.tsv"),
             "export": os.path.join(outdir, stem + "_export.json")}
    header = "#rck-ecdna aa graph {f} source=subgraph name={n} shift={s} imputed={i}".format(f=GRAPH_FORMAT, n=name, s=SHIFT, i=imputed)
    aaf.write_graph(paths["graph"], sequences, edges, [(t, s) for t, s, *_r in constraints], runs, header=header, comments=comments)

    rows = []
    for sid in order:
        s, so = graph.segments[sid], obs.segments[sid]
        rows.append({"chrom": s.chrom, "start": s.start + SHIFT, "end": s.end + SHIFT, "cn": so.total, "cn_a": so.a, "cn_b": so.b})
    pd.DataFrame(rows, columns=["chrom", "start", "end", "cn", "cn_a", "cn_b"]).to_csv(paths["segment_cn"], sep="\t", index=False)

    info = {
        "format": GRAPH_FORMAT, "created_by": "rck.ecdna.aa_export subgraph", "subgraph_dir": os.path.abspath(subgraph_dir), "name": name,
        "amplicon": amplicon, "shift": SHIFT, "imputed": imputed, "min_path_junctions": min_path_junctions,
        "counts": {"sequence": len(sequences), "concordant": len(conc), "discordant": len(disc), "path_constraint": len(constraints),
                   "interval": len(runs), "imputed": len(imp), "outside": len(outside), "rck_path": len(rck_paths)},
        "intervals": ["{c}:{a}-{b}".format(c=c, a=a, b=b) for c, a, b in runs],
        "scale": {"depth_per_copy": obs.scale.depth_per_copy, "junction_reads_per_copy": theta, "baseline_cn": obs.scale.baseline_cn,
                  "depth_method": obs.scale.depth_method, "junction_method": obs.scale.junction_method,
                  "background_cn": float(sub.get("background_cn")) if sub is not None and sub.get("background_cn") is not None else None},
        "paths": dict(pstats),
        "foldbacks": stats["foldbacks"],
        "foldback_note": "fold-back junctions (both ends on one node) are written as discordant self-loops; see results/graph_route/README.md "
                         "for how each tool treats them" if stats["foldbacks"] else None,
        "concordant_cn_from_segments": stats["concordant_cn_from_segments"],
        "outside": stats["outside"], "imputed_edges": stats["imputed"],
        "normalization": list(graph.log_messages),
        "files": {k: os.path.basename(v) for k, v in paths.items()},
    }
    with open(paths["export"], "w") as handle:
        json.dump(info, handle, indent=2)
    return dict(paths, info=info)


def _path_edges(graph, nodes, seq_index, edge_index):
    """A read path's node steps as a CoRAL walk: (path, reason it cannot be one or None,
    observed junctions on it, junctions ``rck-ecdna`` reads from it).

    Steps (2i, 2i+1) run along segments and reference joins, steps (2i+1, 2i+2) are one
    adjacency (:func:`rck.ecdna.walk.read_paths`). The walk is CoRAL's layout,
    ``[e, exit, bp, entry, e, exit, bp, ...]``: one node strand between two edges, the exit
    node of a sequence edge before a breakpoint edge and the entry node after it."""
    from rck.ecdna.subgraph_io import IMPUTED, NOVEL, REFERENCE

    def between(x, y, kind):
        hits = [a for a in graph.adjacencies_at(x) if a.kind == kind and a.other_end(x) == y]
        return hits[0] if hits else None

    path: List[object] = []
    reason, n_novel, n_rck = None, 0, 0

    def bp(edge, at):  # a breakpoint edge leaving node ``at``
        path.append("+" if at == graph.segment_at(at).head else "-")
        path.append(edge)

    for i in range(0, len(nodes) - 1, 2):
        cur, end = nodes[i], nodes[i + 1]
        for _ in range(len(graph.segments) + 1):  # the run cur -> end
            seg = graph.segment_at(cur)
            if seg is None:
                return path, "unwalkable", n_novel, n_rck + 1
            if path:
                path.append("-" if cur == seg.tail else "+")
            path.append(("e", seq_index[seg.sid]))
            out = seg.other_end(cur)
            if out == end:
                break
            refs = [a for a in graph.adjacencies_at(out) if a.kind == REFERENCE and graph.segment_at(a.other_end(out)) is not None]
            if len(refs) != 1:
                return path, "unwalkable", n_novel, n_rck + 1
            bp(edge_index[refs[0].aid], out)
            cur = refs[0].other_end(out)
        else:
            return path, "unwalkable", n_novel, n_rck + 1
        if i + 2 >= len(nodes):
            break
        nxt = nodes[i + 2]  # the join end -> next run
        ref = between(end, nxt, REFERENCE)
        if ref is not None:  # a reference join: rck-ecdna does not count it as a junction
            if ref.aid not in edge_index:
                return path, "unwalkable", n_novel, n_rck
            bp(edge_index[ref.aid], end)
            continue
        n_rck += 1
        nov = between(end, nxt, NOVEL)
        if nov is not None and nov.aid in edge_index:
            bp(edge_index[nov.aid], end)
            n_novel += 1
            continue
        reason = reason or ("imputed" if between(end, nxt, IMPUTED) is not None else "no_edge")
        bp(("x", -1), end)
    return path, reason, n_novel, n_rck


def main_subgraph(argv) -> int:
    parser = argparse.ArgumentParser(prog="python -m rck.ecdna.aa_export subgraph", description=export_subgraph.__doc__,
                                     formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("subgraph_dir")
    parser.add_argument("outdir")
    parser.add_argument("--name", default=None, help="file name stem (default: the folder name)")
    parser.add_argument("--amplicon", type=int, default=1, help="N in <name>_amplicon<N>_graph.txt (CoRAL reads it from the name)")
    parser.add_argument("--imputed", choices=["comment", "drop"], default="comment",
                        help="imputed candidates as #rck_imputed comment lines (default) or left out")
    parser.add_argument("--min-path-junctions", type=int, default=1, help="observed junctions a read path needs to become a path_constraint")
    parser.add_argument("--depth-per-copy", type=float, default=None, help="as rck-ecdna (needed without scale.json)")
    parser.add_argument("--junction-reads-per-copy", type=float, default=None, help="as rck-ecdna (needed without scale.json)")
    parser.add_argument("--baseline-cn", type=float, default=None, help="as rck-ecdna")
    parser.add_argument("--phantom-max-bp", type=int, default=2, help="as rck-ecdna")
    args = parser.parse_args(argv)
    out = export_subgraph(args.subgraph_dir, args.outdir, name=args.name, amplicon=args.amplicon, imputed=args.imputed,
                          min_path_junctions=args.min_path_junctions, depth_per_copy=args.depth_per_copy,
                          junction_reads_per_copy=args.junction_reads_per_copy, baseline_cn=args.baseline_cn,
                          phantom_max_bp=args.phantom_max_bp)
    json.dump({k: v for k, v in out.items() if k != "info"} | {"counts": out["info"]["counts"]}, sys.stdout, indent=2)
    print()
    return 0


def main(argv=None) -> int:
    argv = sys.argv[1:] if argv is None else list(argv)
    if argv and argv[0] == "subgraph" and not os.path.exists(os.path.join("subgraph", "summary.json")):
        return main_subgraph(argv[1:])
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("run_dir")
    parser.add_argument("outdir")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--pool-index", type=int, default=None, help="export this entry of the kept pool (0-based)")
    group.add_argument("--lambda", dest="lam", type=float, default=None, help="export the pool's rank 1 at this --step-penalty")
    parser.add_argument("--alternatives", type=int, default=2, help="same-fit alternatives to export for best-of (default 2)")
    args = parser.parse_args(argv)
    index = args.pool_index
    if args.lam is not None:
        from rck.ecdna.rerank import rerank

        index = int(rerank(args.run_dir, [args.lam])["pool_index"].iloc[0])
    paths = export(args.run_dir, args.outdir, index, args.alternatives)
    json.dump(paths, sys.stdout, indent=2)
    print()
    return 0


if __name__ == "__main__":
    sys.exit(main())

"""Check AA graphs against CoRAL 3.1.1's own parser, sorter and writer (no solver, no licence).

Run with the coral-3.1.1 environment's python and RCK on the path:

    PYTHONPATH=~/RCK ~/.conda/envs/coral-3.1.1/bin/python check_coral_format.py GRAPH... [--paths-from-export]
    PYTHONPATH=~/RCK ~/.conda/envs/coral-3.1.1/bin/python check_coral_format.py --golden OUT.tsv

For every graph (NEXT_STEPS.md step 2c):

1. ``parse_breakpoint_graph`` reads it; its sequence / concordant / discordant / path-constraint
   counts equal ours (``rck.ecdna.aa_format.read_graph``) and ``*_export.json``'s;
2. ``infer_amplicon_intervals`` (what ``coral cycle`` calls) gives our ``interval`` lines and
   does not reorder the edges: it sorts sequence and concordant edges without renumbering the
   path constraints, so a graph out of CoRAL's order would silently corrupt them;
3. every ``path_constraint`` comes back from CoRAL's ``path_to_str`` unchanged;
4. with ``--paths-from-export`` (graphs from ``aa_export subgraph``): every path constraint is
   rebuilt independently, by walking the subgraph's read path (``paths.tsv`` row named by its
   ``#rck_pc`` line, nodes placed by the ``#rck_node`` lines) over CoRAL's own adjacency
   structures, and written with CoRAL's ``path_to_str``: the tokens must be ours. For a truth
   graph with ``*_path_constraints.tsv`` (truth_mixture.py ``--path-constraints``), the
   discordant edges CoRAL reads from each constraint must be the junctions its read crosses;
5. CoRAL's writer (``output_breakpoint_graph_lr``) writes it back and our reader reads that
   file: same sequence edges, concordant and discordant edges, CN and path constraints.

``--golden OUT.tsv`` writes the reference tokens for tests/ecdna/test_aa_import.py: the
amplicon_0 test fixture with the synthetic read paths of tests/ecdna/data/graph_route/amplicon0_paths.tsv, exported
and then encoded by CoRAL itself (step 4).

``--ce GRAPH...`` runs CycleExtractor (~/tools/CycleExtractor/CE.py, this python, HiGHS) on each
graph and reads its cycles back with ``aa_format.read_cycles``.
"""
import argparse
import json
import os
import shutil
import subprocess
import sys
import tempfile

import pandas as pd
from coral import core_utils
from coral.breakpoint.breakpoint_utilities import output_breakpoint_graph_lr
from coral.breakpoint.parse_graph import parse_breakpoint_graph
from coral.datatypes import EdgeId, EdgeType, Node, OutputPCOptions, Strand

from rck.ecdna import aa_format as aaf

HERE = os.path.dirname(os.path.abspath(__file__))
RCK = os.path.abspath(os.path.join(HERE, "..", ".."))
CE = os.path.expanduser("~/tools/CycleExtractor/CE.py")


def coral_walk(g, coords, concordant_pairs):
    """CoRAL Walk of a read path given as node coordinates (AA frame), built on CoRAL's own
    ``node_adjacencies``: steps (2i, 2i+1) run along sequence and concordant edges, steps
    (2i+1, 2i+2) are one concordant or discordant edge. None if the graph cannot carry it."""
    nodes = [Node(c, p, Strand(s)) for c, p, s in coords]
    walk = []

    def seq_at(n):
        idx = g.node_adjacencies[n].sequence
        return idx[0] if idx else None

    def other(i, n):
        e = g.sequence_edges[i]
        return e.end_node if n == e.start_node else e.start_node

    for i in range(0, len(nodes) - 1, 2):
        cur, end = nodes[i], nodes[i + 1]
        for _ in range(len(g.sequence_edges) + 1):
            k = seq_at(cur)
            if k is None:
                return None
            if walk:
                walk.append(cur)
            walk.append(EdgeId(EdgeType.SEQUENCE, k))
            out = other(k, cur)
            if out == end:
                break
            cs = [j for j in g.node_adjacencies[out].concordant]
            if len(cs) != 1:
                return None
            walk.append(out)
            walk.append(EdgeId(EdgeType.CONCORDANT, cs[0]))
            ce = g.concordant_edges[cs[0]]
            cur = ce.node2 if out == ce.node1 else ce.node1
        else:
            return None
        if i + 2 < len(nodes):
            y, z = nodes[i + 1], nodes[i + 2]
            hit = None
            for j in g.node_adjacencies[y].concordant:
                ce = g.concordant_edges[j]
                if {ce.node1, ce.node2} == {y, z}:
                    hit = EdgeId(EdgeType.CONCORDANT, j)
            if hit is None:
                for j in g.node_adjacencies[y].discordant:
                    de = g.discordant_edges[j]
                    if {de.node1, de.node2} == {y, z} and (de.node1 == de.node2) == (y == z):
                        hit = EdgeId(EdgeType.DISCORDANT, j)
            if hit is None:
                return None
            walk.append(y)
            walk.append(hit)
    return walk


def export_paths(graph_path, aa):
    """(constraint number, node coordinates of its read path) for a graph from aa_export subgraph."""
    with open(graph_path.replace("_graph.txt", "_export.json")) as handle:
        info = json.load(handle)
    rows = pd.read_csv(os.path.join(info["subgraph_dir"], "paths.tsv"), sep="\t")
    where = {f[0]: (f[1], int(f[2]), "+" if f[3] == "." else f[3]) for f in aa.rck_lines("rck_node")}
    ends = {}
    for f in aa.rck_lines("rck_seq"):
        ends[f[2]], ends[f[3]] = f[3], f[2]  # the other end of each node's segment
    out, skipped = [], 0
    for f in aa.rck_lines("rck_pc"):
        k, row = int(f[0]), rows.iloc[int(f[1])]
        steps = ["{e}_{o}".format(e=t.rpartition("@")[0], o=t.rpartition("@")[2]) for t in str(row["PATH"]).split(",")]
        # an outer end outside the graph: the segment the read leaves / enters by (as the exporter does)
        if steps[0] not in where and steps[1] in ends:
            steps[0] = ends[steps[1]]
        if steps[-1] not in where and steps[-2] in ends:
            steps[-1] = ends[steps[-2]]
        if any(s not in where for s in steps):
            skipped += 1  # a phantom-contracted node: the exporter followed the alias
            continue
        out.append((k, [where[s] for s in steps]))
    return out, skipped


def check(graph_path, paths_from_export):
    problems = []
    ours = aaf.read_graph(graph_path)
    with open(graph_path) as handle:
        g = parse_breakpoint_graph(handle)
    counts = {"sequence": len(g.sequence_edges), "concordant": len(g.concordant_edges), "discordant": len(g.discordant_edges),
              "path_constraint": len(g.path_constraints)}
    mine = {"sequence": len(ours.sequences), "concordant": len(ours.concordant), "discordant": len(ours.discordant),
            "path_constraint": len(ours.path_constraints)}
    if counts != mine:
        problems.append("counts: CoRAL {c} vs ours {m}".format(c=counts, m=mine))
    export = graph_path.replace("_graph.txt", "_export.json")
    if os.path.exists(export):
        with open(export) as handle:
            ec = json.load(handle)["counts"]
        if any(ec[k] != v for k, v in counts.items()):
            problems.append("counts: CoRAL {c} vs export.json {e}".format(c=counts, e=ec))
    before = ([(e.chr, e.start, e.end) for e in g.sequence_edges], [(e.node1, e.node2) for e in g.concordant_edges])
    g.infer_amplicon_intervals(0)
    after = ([(e.chr, e.start, e.end) for e in g.sequence_edges], [(e.node1, e.node2) for e in g.concordant_edges])
    if before != after:
        problems.append("CoRAL's sort_edges reordered the edges: path constraint indices would be wrong")
    intervals = [(i.chr, i.start, i.end) for i in g.amplicon_intervals]
    if ours.intervals and intervals != ours.intervals:
        problems.append("intervals: CoRAL {c} vs ours {m}".format(c=intervals, m=ours.intervals))
    for i, pc in enumerate(g.path_constraints):
        text = core_utils.path_to_str(pc.path, core_utils.path_to_edge_count(pc.path))
        if text != ours.path_constraints[i][0]:
            problems.append("path constraint {i}: CoRAL writes {c}, file has {m}".format(i=i + 1, c=text, m=ours.path_constraints[i][0]))
    sidecar = graph_path.replace("_graph.txt", "_path_constraints.tsv")
    n_sidecar = 0
    if os.path.exists(sidecar):  # truth_mixture.py --path-constraints: the junctions each constraint must cross
        table = pd.read_csv(sidecar, sep="\t", dtype=str)
        for i, (pc, want) in enumerate(zip(g.path_constraints, table["junctions"])):
            got = []
            for j in range(1, len(pc.path) - 1, 2):
                e = pc.path[j + 1] if isinstance(pc.path[j], Node) and j + 1 < len(pc.path) else None
                if isinstance(e, tuple) and getattr(e, "type", None) == EdgeType.DISCORDANT:
                    de = g.discordant_edges[e.idx]
                    got.append(aaf.lower_first((de.node1.chr, de.node1.pos, de.node1.strand.value), (de.node2.chr, de.node2.pos, de.node2.strand.value)))
            want_l = [aaf.lower_first(*(aaf.parse_breakend(x) for x in w.split("->"))) for w in str(want).split(",")]
            if got != want_l and got[::-1] != want_l:
                problems.append("path constraint {i}: CoRAL reads junctions {g}, the read crosses {w}".format(i=i + 1, g=got, w=want_l))
            n_sidecar += 1
    n_rebuilt = n_skipped = 0
    if paths_from_export:
        pairs = {(e.node1, e.node2) for e in g.concordant_edges}
        rebuilt, n_skipped = export_paths(graph_path, ours)
        for k, coords in rebuilt:
            walk = coral_walk(g, coords, pairs)
            if walk is None:
                problems.append("path constraint {k}: its read path does not follow CoRAL's edges".format(k=k))
                continue
            text = core_utils.path_to_str(walk, core_utils.path_to_edge_count(walk))
            if text != ours.path_constraints[k - 1][0]:
                problems.append("path constraint {k}: rebuilt by CoRAL as {c}, file has {m}".format(k=k, c=text, m=ours.path_constraints[k - 1][0]))
            n_rebuilt += 1
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "rewrite_amplicon1_graph.txt")
        output_breakpoint_graph_lr(g, out, OutputPCOptions.ALL)
        back = aaf.read_graph(out)
    seq = lambda gr: [(s.chrom, s.start, s.end, s.cn) for s in gr.sequences]  # noqa: E731
    edg = lambda gr, k: [(aaf.lower_first(e.a, e.b), e.cn) for e in gr.edges if e.kind == k]  # noqa: E731
    if seq(back) != seq(ours):
        problems.append("CoRAL rewrite: sequence edges differ")
    for kind in ("concordant", "discordant"):
        if sorted(edg(back, kind)) != sorted(edg(ours, kind)):
            problems.append("CoRAL rewrite: {k} edges differ".format(k=kind))
    if [p for p, _ in back.path_constraints] != [p for p, _ in ours.path_constraints]:
        problems.append("CoRAL rewrite: path constraints differ")
    return {"graph": graph_path, **counts, "intervals": len(intervals), "rebuilt": n_rebuilt + n_sidecar, "rebuild_skipped": n_skipped,
            "problems": problems}


def golden(out_tsv):
    """Reference tokens for the unit test: amplicon_0 + tests/ecdna/data/graph_route/amplicon0_paths.tsv, exported,
    and the paths encoded by CoRAL (coral_walk + path_to_str)."""
    from rck.ecdna.aa_export import export_subgraph

    with tempfile.TemporaryDirectory() as tmp:
        sub = os.path.join(tmp, "amplicon_0")
        shutil.copytree(os.path.join(RCK, "tests", "ecdna", "data", "amplicon_0"), sub)
        shutil.copyfile(os.path.join(RCK, "tests", "ecdna", "data", "graph_route", "amplicon0_paths.tsv"), os.path.join(sub, "paths.tsv"))
        res = export_subgraph(sub, os.path.join(tmp, "aa"), name="amplicon0", depth_per_copy=15.0, junction_reads_per_copy=6.0)
        result = check(res["graph"], True)
        aa = aaf.read_graph(res["graph"])
        with open(res["graph"]) as handle:
            g = parse_breakpoint_graph(handle)
        g.infer_amplicon_intervals(0)
        rows = pd.read_csv(os.path.join(sub, "paths.tsv"), sep="\t")
        rebuilt, _ = export_paths(res["graph"], aa)
        by_pid = {}
        for k, coords in rebuilt:
            walk = coral_walk(g, coords, None)
            pid = aa.rck_lines("rck_pc")[k - 1][2]
            by_pid[pid] = core_utils.path_to_str(walk, core_utils.path_to_edge_count(walk))
    table = rows.assign(coral_tokens=[by_pid.get(p, "") for p in rows["PID"]])
    table.to_csv(out_tsv, sep="\t", index=False)
    print(table.to_string(index=False))
    print("problems:", result["problems"] or "none")
    return 0 if not result["problems"] else 1


def run_ce(graphs, outdir):
    os.makedirs(outdir, exist_ok=True)
    for graph in graphs:
        stem = os.path.basename(graph).replace("_graph.txt", "")
        cycles = os.path.join(outdir, stem + "_ce_highs_cycles.txt")
        log = os.path.join(outdir, stem + "_ce_highs.log")
        with open(os.path.join(outdir, stem + "_ce_highs.stdout"), "w") as handle:
            rc = subprocess.call(["nice", "-n", "10", sys.executable, "-I", CE, "--graph", graph, "--output", cycles, "--log_file", log,
                                  "--solver", "highs"], stdout=handle, stderr=subprocess.STDOUT)
        n = len(aaf.read_cycles(cycles).cycles) if rc == 0 and os.path.exists(cycles) else None
        print("{s}: CE exit {rc}, {n} cycles/paths read back -> {c}".format(s=stem, rc=rc, n=n, c=cycles))


def main():
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("graphs", nargs="*")
    parser.add_argument("--paths-from-export", action="store_true", help="rebuild every path constraint from the subgraph's read paths")
    parser.add_argument("--golden", default=None, metavar="OUT.tsv")
    parser.add_argument("--ce", default=None, metavar="OUTDIR", help="run CycleExtractor (HiGHS) on the graphs, outputs in OUTDIR")
    parser.add_argument("--json", default=None, help="write the per-graph results here")
    args = parser.parse_args()
    if args.golden:
        return golden(args.golden)
    if args.ce:
        run_ce(args.graphs, args.ce)
        return 0
    results, bad = [], 0
    for path in args.graphs:
        r = check(path, args.paths_from_export)
        results.append(r)
        bad += bool(r["problems"])
        print("{g}: {s}".format(g=os.path.basename(path), s="OK" if not r["problems"] else "; ".join(r["problems"])),
              "(seq {sequence}, conc {concordant}, disc {discordant}, pc {path_constraint}, intervals {intervals}, rebuilt {rebuilt}, "
              "skipped {rebuild_skipped})".format(**r))
    if args.json:
        with open(args.json, "w") as handle:
            json.dump(results, handle, indent=2)
    print("{n} graphs, {b} with problems".format(n=len(results), b=bad))
    return 1 if bad else 0


if __name__ == "__main__":
    sys.exit(main())

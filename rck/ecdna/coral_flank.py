"""Add flank sequence edges to an AA graph so that CoRAL 3.1.1 ``cycle`` can use breakpoints at its ends.

    python -m rck.ecdna.coral_flank GRAPH OUT_GRAPH [--flank 10000]

CoRAL infers amplicon intervals from runs of contiguous sequence edges and makes the two outer
nodes of each run *end nodes* (``BreakpointGraph.infer_amplicon_intervals``). At an end node its
Eulerian constraint is ``x(s,v) + x(v,t) == x(sequence edge)`` (``models/constraints.py``,
``get_eulerian_node_constraints``): discordant and concordant edges there are left out of the
balance, so a cycle through a breakpoint at an interval end is infeasible, and CoRAL returns
zero-weight paths. CoRAL's own graphs extend each amplicon beyond its breakpoints; graphs cut at
the breakpoints (the truth mixture graphs, ``aa_export subgraph``) do not.

This adds, at both ends of every run, a flank sequence edge of ``--flank`` bp (shorter when the
gap to the next run on the chromosome is small; the flanks of two runs never touch) and a
concordant edge joining it to the run. Its copy number is the node-balance residual at the old end
node: the segment's CN minus the CN of the discordant edges there (a self-loop counted twice, as
CoRAL lists it twice), floored at 0; AverageCoverage, NumberReadsMapped and the concordant read
count are 0. Discordant edges and their order are unchanged; sequence and concordant edges are
written in CoRAL's sort order and the path constraints' ``e``/``c`` indices renumbered. ``#``
comment lines and ``interval`` lines are not copied (CoRAL ignores them; the ``#rck_*`` ids would
no longer match). ``source`` edges are kept.

RCK-ecDNA sees the same picture without this (its region-flank joins and flank background), and
CycleExtractor connects s/t at every node, so the graph route gives the flanked graph to CoRAL only.
"""
import argparse
import json
import re
from typing import Dict, List, Tuple

from rck.ecdna import aa_format as aaf

DEFAULT_FLANK = 10000
TOKEN = re.compile(r"^([ecds])(\d+)([+-]):(\d+)$")


def _runs(sequences: List[aaf.AASequence]) -> List[List[int]]:
    order = sorted(range(len(sequences)), key=lambda i: aaf.sequence_key(sequences[i]))
    runs: List[List[int]] = []
    for i in order:
        s = sequences[i]
        if runs:
            p = sequences[runs[-1][-1]]
            if p.chrom == s.chrom and s.start == p.end + 1:
                runs[-1].append(i)
                continue
        runs.append([i])
    return runs


def _renumber(tokens: str, maps: Dict[str, Dict[int, int]]) -> str:
    out = []
    for t in tokens.split(","):
        m = TOKEN.match(t.strip())
        if not m:
            raise ValueError("unreadable path-constraint token {t!r}".format(t=t))
        kind, idx, strand, count = m.group(1), int(m.group(2)), m.group(3), m.group(4)
        idx = maps.get(kind, {}).get(idx, idx)
        out.append("{k}{i}{s}:{c}".format(k=kind, i=idx, s=strand, c=count))
    return ",".join(out)


def add_flanks(graph: aaf.AAGraph, flank: int = DEFAULT_FLANK):
    """(sequences, edges, path_constraints, info) of the flanked graph, in CoRAL's order."""
    if flank < 1:
        raise ValueError("--flank must be at least 1 bp")
    seqs = list(graph.sequences)
    disc = [e for e in graph.edges if e.kind == "discordant"]
    conc = [e for e in graph.edges if e.kind == "concordant"]
    other = [e for e in graph.edges if e.kind not in ("discordant", "concordant")]

    def disc_cn(node) -> float:
        return sum(e.cn for e in disc for end in (e.a, e.b) if end == node)

    runs = _runs(seqs)
    bounds = [(seqs[r[0]].chrom, seqs[r[0]].start, seqs[r[-1]].end) for r in runs]
    new_seqs, new_conc, added = list(seqs), list(conc), []
    for k, (chrom, start, end) in enumerate(bounds):
        prev_end = max([e for c, _s, e in bounds[:k] if c == chrom] or [0])
        later = [s for c, s, _e in bounds[k + 1:] if c == chrom]
        next_start = min(later) if later else None
        # each flank takes at most (gap - 1) // 2 of a gap, so two flanks never touch (CoRAL would merge them)
        left = min(flank, (start - prev_end - 2) // 2 if prev_end else start - 1)
        right = min(flank, (next_start - end - 2) // 2) if next_start is not None else flank
        for side, length in (("left", left), ("right", right)):
            node = (chrom, start, "-") if side == "left" else (chrom, end, "+")
            run_seq = seqs[runs[k][0]] if side == "left" else seqs[runs[k][-1]]
            if length < 1:
                raise ValueError("no room for a {s} flank at {c}:{p} (the next run is too close)".format(s=side, c=chrom, p=node[1]))
            cn = max(0.0, run_seq.cn - disc_cn(node))
            if side == "left":
                fs, fe = start - length, start - 1
                join = aaf.AAEdge("concordant", (chrom, fe, "+"), (chrom, start, "-"), cn, 0.0)
            else:
                fs, fe = end + 1, end + length
                join = aaf.AAEdge("concordant", (chrom, end, "+"), (chrom, fs, "-"), cn, 0.0)
            new_seqs.append(aaf.AASequence(chrom, fs, fe, cn, coverage=0.0, size=length, reads=0.0))
            new_conc.append(join)
            added.append({"chrom": chrom, "start": fs, "end": fe, "cn": cn, "side": side,
                          "residual": run_seq.cn - disc_cn(node)})
    new_seqs.sort(key=aaf.sequence_key)
    new_conc.sort(key=lambda e: aaf.edge_key(*aaf.lower_first(e.a, e.b)))
    seq_map = {i + 1: next(j for j, t in enumerate(new_seqs) if t is s) + 1 for i, s in enumerate(seqs)}
    conc_map = {i + 1: next(j for j, t in enumerate(new_conc) if t is e) + 1 for i, e in enumerate(conc)}
    pcs = [(_renumber(tokens, {"e": seq_map, "c": conc_map}), support) for tokens, support in graph.path_constraints]
    info = {"flank": flank, "runs": len(runs), "flanks": added,
            "negative_residuals": [a for a in added if a["residual"] < 0]}
    return new_seqs, new_conc + disc + other, pcs, info


def flank_graph(src: str, dst: str, flank: int = DEFAULT_FLANK) -> dict:
    """Write the flanked graph to ``dst`` (its name must keep CoRAL's ``_amplicon<N>_graph.txt``)."""
    graph = aaf.read_graph(src)
    seqs, edges, pcs, info = add_flanks(graph, flank)
    aaf.write_graph(dst, seqs, edges, pcs, header="# coral_flank: {s} with {f} bp flanks (rck.ecdna.coral_flank)".format(s=src, f=flank))
    info["source"] = src
    return info


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("graph")
    parser.add_argument("out")
    parser.add_argument("--flank", type=int, default=DEFAULT_FLANK)
    args = parser.parse_args(argv)
    print(json.dumps(flank_graph(args.graph, args.out, args.flank), indent=2))


if __name__ == "__main__":
    main()

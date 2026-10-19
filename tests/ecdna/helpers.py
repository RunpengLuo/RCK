"""Builders for small synthetic inputs in the formats rck.ecdna reads.

``SubgraphBuilder`` writes ``nodes.tsv`` / ``edges.tsv`` in the layout of
amplicon-allele-analysis ``subgraph_construct.py``. Breakpoints use the real-data naming
(one extremity per breakpoint, nodes ``<chrom>:<pos>_+`` at ``pos`` and ``<chrom>:<pos>_-`` at
``pos + 1``), so there are no phantom segments unless a test adds them on purpose.
"""

import os
import unittest

import numpy as np
import pandas as pd

try:  # pragma: no cover - depends on the environment
    import gurobipy  # noqa: F401

    HAVE_GUROBI = True
except ImportError:  # pragma: no cover
    HAVE_GUROBI = False

requires_gurobi = unittest.skipUnless(HAVE_GUROBI, "gurobipy is not installed")

DATA_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

NODE_COLUMNS = ["NODE_ID", "EID", "#CHR", "POS", "ORIENT", "SOURCE"]
EDGE_COLUMNS = ["NODE_1", "NODE_2", "EDGE_TYPE", "WEIGHT", "CHROM_1", "POS_1", "CHROM_2", "POS_2", "AID",
                "ORIENT_1", "ORIENT_2", "n_snps", "BAF", "NOV_A", "NOV_B", "NOV_U", "junction_hap", "seg_hap"]


def head(chrom, pos):
    """Node id of the head (right end) of the segment ending at ``pos``."""
    return "{c}:{p}_+".format(c=chrom, p=pos)


def tail(chrom, pos):
    """Node id of the tail (left end) of the segment starting at ``pos + 1``."""
    return "{c}:{p}_-".format(c=chrom, p=pos)


class SubgraphBuilder:
    def __init__(self):
        self.nodes = {}
        self.edges = []

    def node(self, nid, chrom, pos, orient, source="extremity", eid=None):
        self.nodes[nid] = {"NODE_ID": nid, "EID": eid or nid.rsplit("_", 1)[0], "#CHR": chrom, "POS": pos,
                           "ORIENT": orient, "SOURCE": source}
        return nid

    def edge(self, n1, n2, etype, weight, **extra):
        a, b = self.nodes[n1], self.nodes[n2]
        row = {c: None for c in EDGE_COLUMNS}
        row.update({"NODE_1": n1, "NODE_2": n2, "EDGE_TYPE": etype, "WEIGHT": weight,
                    "CHROM_1": a["#CHR"], "POS_1": a["POS"], "CHROM_2": b["#CHR"], "POS_2": b["POS"],
                    "ORIENT_1": a["ORIENT"], "ORIENT_2": b["ORIENT"]})
        row.update(extra)
        self.edges.append(row)

    def chain(self, chrom, breakpoints, depths, bafs=None, ref_reads=10, n_snps=50):
        """Breakpoints b0 < b1 < ... ; segment i spans (b_i, b_{i+1}]. Returns segment (tail, head) node pairs."""
        assert len(depths) == len(breakpoints) - 1
        bafs = bafs or [0.05] * len(depths)
        if not isinstance(ref_reads, (list, tuple)):
            ref_reads = [ref_reads] * len(breakpoints)
        for bp, reads in zip(breakpoints, ref_reads):
            self.node(head(chrom, bp), chrom, bp, "+")
            self.node(tail(chrom, bp), chrom, bp + 1, "-")
            self.edge(head(chrom, bp), tail(chrom, bp), "extremity", reads)
        segs = []
        for i, depth in enumerate(depths):
            t, h = tail(chrom, breakpoints[i]), head(chrom, breakpoints[i + 1])
            self.edge(t, h, "reference", depth, n_snps=n_snps if bafs[i] is not None else 0, BAF=bafs[i],
                      seg_hap="A" if (bafs[i] or 0) < 0.5 else "B")
            segs.append((t, h))
        return segs

    def variant(self, n1, n2, reads, hap="A", aid=None):
        aid = aid or "v{i}".format(i=len(self.edges))
        nov_a = reads if hap == "A" else 0 if hap == "B" else reads / 2
        nov_b = reads if hap == "B" else 0 if hap == "A" else reads / 2
        self.edge(n1, n2, "variant", reads, AID=aid, NOV_A=nov_a, NOV_B=nov_b, NOV_U=0, junction_hap=hap)
        return aid

    def imputed(self, n1, n2, weight, aid):
        self.edge(n1, n2, "imputed_variant", weight, AID=aid, junction_hap="imputed")
        return aid

    def write(self, directory):
        os.makedirs(directory, exist_ok=True)
        pd.DataFrame(list(self.nodes.values()), columns=NODE_COLUMNS).to_csv(os.path.join(directory, "nodes.tsv"), sep="\t", index=False)
        pd.DataFrame(self.edges, columns=EDGE_COLUMNS).to_csv(os.path.join(directory, "edges.tsv"), sep="\t", index=False)
        return directory


def single_circle(directory, depth=180.0, reads=60.0):
    """chr1 with three segments; an ecDNA joins the head of the last back to the tail of the first.

    With depth per copy 15 and junction reads per copy 6 this is 10 ecDNA copies on A plus one
    background copy per haplotype: depth 180 = 12 copies at BAF 1/12, junction 60 reads = 10
    copies, interior reference joins 72 reads = 12 copies, flank joins 12 reads = 2 copies.
    """
    b = SubgraphBuilder()
    segs = b.chain("chr1", [1000, 101000, 201000, 301000], [depth] * 3, bafs=[1 / 12.0] * 3,
                   ref_reads=[12, 12 + reads, 12 + reads, 12])
    b.variant(segs[-1][1], segs[0][0], reads, hap="A", aid="circ")
    return b.write(directory), segs


def write_windows(directory, chroms, sample_rows=None):
    """Write a UGP-style windows directory.

    ``chroms`` maps chromosome -> list of (n_windows, tumor_depth, normal_depth, A_count, B_count)
    runs, laid out left to right in 1 kb windows starting at 0.
    """
    os.makedirs(directory, exist_ok=True)
    rows, depth, snp_rows, a_counts, b_counts = [], [], [], [], []
    for chrom, runs in chroms.items():
        pos = 0
        for n, t, nd, a, bb in runs:
            for _ in range(n):
                rows.append((chrom, pos, pos + 1000))
                depth.append((nd, t))
                snp_rows.append((chrom, pos + 500))
                a_counts.append((15, a))
                b_counts.append((15, bb))
                pos += 1000
    pd.DataFrame(rows, columns=["#CHR", "START", "END"]).to_csv(os.path.join(directory, "window.tsv.gz"), sep="\t", index=False)
    np.savez(os.path.join(directory, "window.depth.npz"), mat=np.array(depth, dtype=float))
    pd.DataFrame(snp_rows, columns=["#CHR", "POS"]).to_csv(os.path.join(directory, "snp.tsv.gz"), sep="\t", index=False)
    np.savez(os.path.join(directory, "snp.Aallele.npz"), mat=np.array(a_counts, dtype=float))
    np.savez(os.path.join(directory, "snp.Ballele.npz"), mat=np.array(b_counts, dtype=float))
    sample_rows = sample_rows or [("sim_normal", "sim", "normal", "normal"), ("sim_tumor1", "sim", "tumor1", "tumor")]
    pd.DataFrame(sample_rows, columns=["SAMPLE", "sample_id", "dataset_id", "sample_type"]).to_csv(
        os.path.join(directory, "sample_ids.tsv"), sep="\t", index=False)
    return directory


def loops_graph(directory, loops=("a", "b", "c")):
    """A hub segment ``x`` (chr1) and one loop segment per name, each on its own chromosome.

    Junctions ``x>NAME`` (x head -> NAME tail) and ``NAME>x`` (NAME head -> x tail) let a walk
    visit ``x`` once per loop. Returns {name: (tail, head)} including ``x``.
    """
    b = SubgraphBuilder()
    ends = {"x": b.chain("chr1", [0, 1000], [100.0])[0]}
    for i, name in enumerate(loops, start=2):
        ends[name] = b.chain("chr{i}".format(i=i), [0, 1000], [100.0])[0]
    for name in loops:
        b.variant(ends["x"][1], ends[name][0], 50, hap="A", aid="x>{n}".format(n=name))
        b.variant(ends[name][1], ends["x"][0], 50, hap="A", aid="{n}>x".format(n=name))
    b.write(directory)
    return ends


def loop_counts(graph, ends, uses, weight=1):
    """Phased (haplotype A) edge multiset for a walk that visits x once per loop in ``uses``."""
    from rck.ecdna.phased import adj_key, seg_key

    def sid(name):
        tail_node = ends[name][0]
        return graph.nodes[tail_node].segment

    counts = {seg_key(sid("x"), "A"): weight * len(uses)}
    for name in uses:
        counts[seg_key(sid(name), "A")] = counts.get(seg_key(sid(name), "A"), 0) + weight
        for aid in ("x>{n}".format(n=name), "{n}>x".format(n=name)):
            counts[adj_key(aid, "AA")] = counts.get(adj_key(aid, "AA"), 0) + weight
    return counts

"""rck.ecdna.aa_export.baseline_graph: the bulk route's graph plus the subgraph export's long-read path
constraints, the input of the decomposition-only baselines (CoRAL, CycleExtractor) in bench_IW_intpsi.
No solver, no CoRAL."""
import json
import os
import tempfile
import unittest

from rck.ecdna import aa_format as aaf
from rck.ecdna.aa_export import SHIFT, baseline_graph


def seq(c, s, e, cn):
    return aaf.AASequence(c, s, e, cn, coverage=cn * 10, size=e - s + 1, reads=0)


def edge(kind, a, b, cn, reads=0.0):
    return aaf.AAEdge(kind, a, b, cn, reads)


def direct_export():
    """A subgraph export (AA frame): demo5-like cycle, a direct-only discordant edge, 3 constraints."""
    seqs = [seq("chr17", 3906001, 5783000, 18.0), seq("chr17", 5783001, 5786000, 18.0), seq("chr17", 61941001, 62060000, 18.0)]
    edges = [edge("concordant", ("chr17", 5783000, "+"), ("chr17", 5783001, "-"), 2.0),
             edge("discordant", ("chr17", 3906001, "-"), ("chr17", 62060000, "+"), 16.0),
             edge("discordant", ("chr17", 5783000, "+"), ("chr17", 5786000, "+"), 16.0),
             edge("discordant", ("chr17", 5783001, "-"), ("chr17", 61941001, "-"), 16.0),
             edge("discordant", ("chr17", 5786000, "+"), ("chr17", 62060000, "+"), 1.0)]  # not in the bulk graph
    pcs = [("e1-:1,d1+:1,e3-:1", "101"), ("e1+:1,d2+:1,e2-:1,d3-:1,e3+:1", "80"), ("e2+:1,d4+:1,e3+:1", "7")]
    return seqs, edges, pcs


def bulk_graph():
    """The bulk route's aa/ graph (walk frame = AA - SHIFT), its edges in another order and direction,
    with an edge to a node outside the sequences and an imputed edge the export only commented."""
    def w(b):
        return b[0], b[1] - SHIFT, b[2]

    seqs = [seq("chr17", 61941001 - SHIFT, 62060000 - SHIFT, 18.0), seq("chr17", 3906001 - SHIFT, 5783000 - SHIFT, 18.0),
            seq("chr17", 5783001 - SHIFT, 5786000 - SHIFT, 0.0)]
    edges = [edge("discordant", w(("chr17", 61941001, "-")), w(("chr17", 5783001, "-")), 16.0),  # written high end first
             edge("discordant", w(("chr17", 3906001, "-")), w(("chr17", 62060000, "+")), 16.0),
             edge("concordant", w(("chr17", 3906000, "+")), w(("chr17", 3906001, "-")), 2.0),  # outside the sequences
             edge("discordant", w(("chr17", 5783000, "+")), w(("chr17", 5786000, "+")), 16.0),
             edge("concordant", w(("chr17", 5783000, "+")), w(("chr17", 5783001, "-")), 2.0),
             edge("discordant", w(("chr17", 5786000, "+")), w(("chr17", 61941001, "-")), 3.0)]  # imputed, used by RCK
    return seqs, edges


class TestBaselineGraph(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = self.tmp.name
        os.makedirs(os.path.join(root, "bulk", "aa"))
        seqs, edges = bulk_graph()
        aaf.write_graph(os.path.join(root, "bulk", "aa", "union_graph.txt"), seqs, edges)
        seqs, edges, pcs = direct_export()
        self.direct = os.path.join(root, "direct_graph.txt")
        aaf.write_graph(self.direct, seqs, edges, pcs)
        self.dst = os.path.join(root, "out", "case_amplicon1_graph.txt")
        os.makedirs(os.path.dirname(self.dst))
        self.info = baseline_graph(os.path.join(root, "bulk"), self.direct, self.dst)
        self.graph = aaf.read_graph(self.dst)

    def tearDown(self):
        self.tmp.cleanup()

    def test_shifted_sorted_bulk_copy_numbers(self):
        self.assertEqual([(s.start, s.end, s.cn) for s in self.graph.sequences],
                         [(3906001, 5783000, 18.0), (5783001, 5786000, 0.0), (61941001, 62060000, 18.0)])
        disc = [e for e in self.graph.edges if e.kind == "discordant"]
        self.assertEqual([(e.a, e.b) for e in disc], sorted([(e.a, e.b) for e in disc], key=lambda p: aaf.edge_key(*p)))
        for e in self.graph.edges:
            self.assertEqual((e.a, e.b), aaf.lower_first(e.a, e.b))
        self.assertIn((("chr17", 5786000, "+"), ("chr17", 61941001, "-")), [(e.a, e.b) for e in disc])  # imputed edge kept

    def test_outside_edge_left_out_and_recorded(self):
        self.assertEqual(self.info["counts"]["outside_dropped"], 1)
        self.assertEqual(self.info["outside_edges"][0]["edge"], "chr17:3906000+->chr17:3906001-")
        self.assertNotIn(("chr17", 3906000, "+"), [e.a for e in self.graph.edges])
        self.assertTrue(any(c.startswith("#rck_outside_bulk") for c in self.graph.comments))

    def test_constraints_renumbered_to_same_coordinates(self):
        direct = aaf.read_graph(self.direct)

        def coords(g, kind, idx):
            if kind == "e":
                s = g.sequences[idx]
                return s.chrom, s.start, s.end
            pool = [e for e in g.edges if e.kind == {"c": "concordant", "d": "discordant"}[kind]]
            return aaf.lower_first(pool[idx].a, pool[idx].b)

        self.assertEqual(len(self.graph.path_constraints), 2)
        for (t_old, s_old), (t_new, s_new) in zip(direct.path_constraints[:2], self.graph.path_constraints):
            self.assertEqual(s_old, s_new)
            for a, b in zip(aaf.parse_path_tokens(t_old), aaf.parse_path_tokens(t_new)):
                self.assertEqual((a[0], a[2], a[3]), (b[0], b[2], b[3]))
                self.assertEqual(coords(direct, a[0], a[1]), coords(self.graph, b[0], b[1]))

    def test_unmappable_constraint_dropped_and_counted(self):
        c = self.info["counts"]
        self.assertEqual((c["path_constraints_in"], c["path_constraints_written"], c["path_constraints_dropped"]), (3, 2, 1))
        self.assertEqual(self.info["path_constraints_dropped_by_reason"], {"no_discordant_edge": 1})

    def test_zero_copy_edges_counted(self):
        self.assertEqual(self.info["counts"]["path_constraints_through_zero_cn"], 1)  # the 80-read path uses e2 (0 copies)

    def test_intervals_and_info_file(self):
        self.assertEqual(self.graph.intervals, [("chr17", 3906001, 5786000), ("chr17", 61941001, 62060000)])
        with open(os.path.splitext(self.dst)[0] + ".json") as handle:
            self.assertEqual(json.load(handle)["role"], "decomposition-only baseline input")


if __name__ == "__main__":
    unittest.main()

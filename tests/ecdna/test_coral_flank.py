"""rck.ecdna.coral_flank: flank sequence edges at interval ends for CoRAL 3.1.1 (no solver, no CoRAL)."""
import os
import tempfile
import unittest

from rck.ecdna import aa_format as aaf
from rck.ecdna.coral_flank import add_flanks, flank_graph


def seq(c, s, e, cn):
    return aaf.AASequence(c, s, e, cn, coverage=cn * 10, size=e - s + 1, reads=0)


def demo5_like():
    """demo5's truth graph: one cycle 1+ d2 2- d3 3+ d1 at CN 16 over a background of 2."""
    g = aaf.AAGraph("x")
    g.sequences = [seq("chr17", 3906001, 5783000, 18.0), seq("chr17", 5783001, 5786000, 18.0), seq("chr17", 61941001, 62060000, 18.0)]
    g.edges = [aaf.AAEdge("concordant", ("chr17", 5783000, "+"), ("chr17", 5783001, "-"), 2.0, 13.0),
               aaf.AAEdge("discordant", ("chr17", 3906001, "-"), ("chr17", 62060000, "+"), 16.0, 104.0),
               aaf.AAEdge("discordant", ("chr17", 5783000, "+"), ("chr17", 5786000, "+"), 16.0, 104.0),
               aaf.AAEdge("discordant", ("chr17", 5783001, "-"), ("chr17", 61941001, "-"), 16.0, 104.0)]
    g.path_constraints = [("e1-:1,d1+:1,e3-:1", "101"), ("e1+:1,d2+:1,e2-:1,d3-:1,e3+:1", "80"), ("e1+:1,c1-:1,e2+:1", "5")]
    return g


class TestCoralFlank(unittest.TestCase):
    def test_flanks_at_every_run_end(self):
        seqs, edges, pcs, info = add_flanks(demo5_like(), 10000)
        self.assertEqual(info["runs"], 2)
        self.assertEqual([(s.start, s.end, s.cn) for s in seqs],
                         [(3896001, 3906000, 2.0), (3906001, 5783000, 18.0), (5783001, 5786000, 18.0), (5786001, 5796000, 2.0),
                          (61931001, 61941000, 2.0), (61941001, 62060000, 18.0), (62060001, 62070000, 2.0)])
        self.assertTrue(all(s.coverage == 0.0 and s.reads == 0.0 for s in seqs if s.size == 10000))
        conc = [e for e in edges if e.kind == "concordant"]
        self.assertEqual([(e.a[1], e.b[1], e.cn) for e in conc],
                         [(3906000, 3906001, 2.0), (5783000, 5783001, 2.0), (5786000, 5786001, 2.0),
                          (61941000, 61941001, 2.0), (62060000, 62060001, 2.0)])
        self.assertEqual(conc[1].reads, 13.0)  # the original join keeps its reads
        self.assertEqual([e for e in edges if e.kind == "discordant"], [e for e in demo5_like().edges if e.kind == "discordant"])
        self.assertEqual(info["negative_residuals"], [])

    def test_path_constraints_renumbered(self):
        _seqs, _edges, pcs, _info = add_flanks(demo5_like(), 10000)
        self.assertEqual(pcs, [("e2-:1,d1+:1,e6-:1", "101"), ("e2+:1,d2+:1,e3-:1,d3-:1,e6+:1", "80"), ("e2+:1,c2-:1,e3+:1", "5")])

    def test_close_runs_get_short_flanks_that_never_touch(self):
        g = aaf.AAGraph("x")
        g.sequences = [seq("chr1", 1001, 2000, 10.0), seq("chr1", 2101, 3000, 10.0)]  # 100 bp gap
        g.edges = [aaf.AAEdge("discordant", ("chr1", 2000, "+"), ("chr1", 2101, "-"), 8.0, 50.0)]
        seqs, _edges, _pcs, _info = add_flanks(g, 10000)
        spans = [(s.start, s.end) for s in seqs]
        self.assertEqual(spans, [(1, 1000), (1001, 2000), (2001, 2049), (2052, 2100), (2101, 3000), (3001, 13000)])
        g.sequences = [seq("chr1", 1001, 2000, 10.0), seq("chr1", 2003, 3000, 10.0)]  # 2 bp gap: no room
        with self.assertRaises(ValueError):
            add_flanks(g, 10000)

    def test_self_loop_counts_twice(self):
        g = aaf.AAGraph("x")
        g.sequences = [seq("chr2", 5001, 9000, 12.0)]
        g.edges = [aaf.AAEdge("discordant", ("chr2", 9000, "+"), ("chr2", 9000, "+"), 4.0, 26.0),
                   aaf.AAEdge("discordant", ("chr2", 5001, "-"), ("chr2", 9000, "+"), 2.0, 13.0)]
        seqs, _edges, _pcs, info = add_flanks(g, 1000)
        right = [s for s in seqs if s.start == 9001][0]
        left = [s for s in seqs if s.end == 5000][0]
        self.assertEqual(right.cn, 12.0 - 2 * 4.0 - 2.0)
        self.assertEqual(left.cn, 10.0)

    def test_negative_residual_floored_and_reported(self):
        g = aaf.AAGraph("x")
        g.sequences = [seq("chr3", 1001, 2000, 3.0)]
        g.edges = [aaf.AAEdge("discordant", ("chr3", 1001, "-"), ("chr3", 2000, "+"), 5.0, 30.0)]
        seqs, _edges, _pcs, info = add_flanks(g, 100)
        self.assertEqual([s.cn for s in seqs], [0.0, 3.0, 0.0])
        self.assertEqual(len(info["negative_residuals"]), 2)

    def test_written_graph_reads_back(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "x_amplicon1_graph.txt")
            g = demo5_like()
            aaf.write_graph(src, g.sequences, g.edges, g.path_constraints, comments=["#rck_node\tfoo"])
            dst = os.path.join(tmp, "coral_amplicon1_graph.txt")
            info = flank_graph(src, dst, 10000)
            back = aaf.read_graph(dst)
            self.assertEqual(len(back.sequences), 7)
            self.assertEqual(len(back.concordant), 5)
            self.assertEqual(len(back.discordant), 3)
            self.assertEqual([t for t, _s in back.path_constraints][0], "e2-:1,d1+:1,e6-:1")
            self.assertEqual(back.rck_lines("rck_node"), [])  # #rck_* lines are not copied
            self.assertEqual(info["source"], src)


if __name__ == "__main__":
    unittest.main()

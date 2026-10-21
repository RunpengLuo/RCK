import os
import unittest

from rck.core.structures import Phasing
from rck.ecdna.subgraph_io import IMPUTED, NOVEL, REFERENCE, read_subgraph
from tests.ecdna.helpers import DATA_DIR, head, single_circle, SubgraphBuilder, tail, temp_dir


class TestChain(unittest.TestCase):
    def test_single_circle_structure(self):
        with temp_dir() as tmp:
            directory, segs = single_circle(tmp)
            g = read_subgraph(directory)
        self.assertEqual(len(g.segments), 3)
        kinds = sorted(a.kind for a in g.adjacencies.values())
        self.assertEqual(kinds, [NOVEL] + [REFERENCE] * 4)
        self.assertEqual(g.validate(), [])
        # the two flank sides with no segment are the only terminals
        self.assertEqual(set(g.terminals()), {head("chr1", 1000), tail("chr1", 301000)})
        circ = g.adjacencies["circ"]
        self.assertEqual(circ.allowed, (Phasing.AA,))
        for nid in g.nodes:
            if nid not in g.terminals():
                self.assertIsNotNone(g.nodes[nid].segment)

    def test_segment_coordinates_and_baf(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g = read_subgraph(directory)
        spans = sorted((s.start, s.end) for s in g.segments.values())
        self.assertEqual(spans, [(1001, 101000), (101001, 201000), (201001, 301000)])
        for s in g.segments.values():
            self.assertAlmostEqual(s.baf, 1 / 12.0)
            self.assertEqual(s.n_snps, 50)


class TestPhantomContraction(unittest.TestCase):
    """Simulator-style naming: breakpoint x has extremities 'x+' and 'x+1-', giving a 1 bp segment."""

    def build(self, tmp):
        b = SubgraphBuilder()
        c = "chr2"
        # left segment ends at 5000 (EID chr2:5000+), right segment starts at 5001 (EID chr2:5001-)
        b.node("chr2:1000-_+", c, 1000, "+")
        b.node("chr2:1000-_-", c, 1001, "-")
        b.node("chr2:5000+_+", c, 5000, "+")
        b.node("chr2:5000+_-", c, 5001, "-")  # phantom tail
        b.node("chr2:5001-_+", c, 5001, "+")  # phantom head
        b.node("chr2:5001-_-", c, 5002, "-")
        b.node("chr2:9000+_+", c, 9000, "+")
        b.node("chr2:9000+_-", c, 9001, "-")
        b.edge("chr2:1000-_+", "chr2:1000-_-", "extremity", 10)
        b.edge("chr2:5000+_+", "chr2:5000+_-", "extremity", 8)
        b.edge("chr2:5001-_+", "chr2:5001-_-", "extremity", 12)
        b.edge("chr2:9000+_+", "chr2:9000+_-", "extremity", 10)
        b.edge("chr2:1000-_-", "chr2:5000+_+", "reference", 100, n_snps=10, BAF=0.1)
        b.edge("chr2:5000+_-", "chr2:5001-_+", "reference", 100, seg_hap="no_span")
        b.edge("chr2:5001-_-", "chr2:9000+_+", "reference", 100, n_snps=10, BAF=0.1)
        # an imputed edge on the phantom head (moves to the true head chr2:5000+_+)
        b.imputed("chr2:9000+_+", "chr2:5001-_+", 50, aid="imp_anchor")
        return b.write(tmp)

    def test_contraction(self):
        with temp_dir() as tmp:
            g = read_subgraph(self.build(tmp))
        self.assertEqual(len(g.segments), 2)
        self.assertNotIn("chr2:5000+_-", g.nodes)
        self.assertNotIn("chr2:5001-_+", g.nodes)
        merged = [a for a in g.adjacencies.values() if {a.u, a.v} == {"chr2:5000+_+", "chr2:5001-_-"}]
        self.assertEqual(len(merged), 1)
        self.assertEqual(merged[0].kind, REFERENCE)
        self.assertAlmostEqual(merged[0].reads, 10.0)  # mean of 8 and 12
        right = [s for s in g.segments.values() if s.tail == "chr2:5001-_-"][0]
        self.assertEqual(right.start, 5001)  # phantom base absorbed
        imp = g.adjacencies["imp_anchor"]
        self.assertEqual({imp.u, imp.v}, {"chr2:9000+_+", "chr2:5000+_+"})
        self.assertEqual(g.resolve("chr2:5001-_+"), "chr2:5000+_+")
        self.assertEqual(g.resolve("chr2:5000+_-"), "chr2:5001-_-")
        self.assertEqual(g.validate(), [])

    def test_contraction_disabled(self):
        with temp_dir() as tmp:
            g = read_subgraph(self.build(tmp), phantom_max_bp=0)
        self.assertEqual(len(g.segments), 3)


class TestParallelAndPhasing(unittest.TestCase):
    def test_parallel_imputed_merged_into_reference(self):
        with temp_dir() as tmp:
            b = SubgraphBuilder()
            b.chain("chr3", [0, 1000, 2000, 3000], [100, 100, 100])
            b.imputed(head("chr3", 1000), tail("chr3", 1000), 90, aid="inter")  # same pair as the reference join
            b.imputed(head("chr3", 3000), tail("chr3", 0), 90, aid="intra")  # a real candidate
            g = read_subgraph(b.write(tmp))
        self.assertNotIn("inter", g.adjacencies)
        ref = [a for a in g.adjacencies.values() if {a.u, a.v} == {head("chr3", 1000), tail("chr3", 1000)}]
        self.assertEqual(len(ref), 1)
        self.assertEqual(ref[0].kind, REFERENCE)
        self.assertEqual(ref[0].imputed_support, 90)
        intra = g.adjacencies["intra"]
        self.assertEqual(intra.kind, IMPUTED)
        self.assertEqual(intra.allowed, (Phasing.AA, Phasing.BB))
        self.assertIsNone(intra.reads)

    def test_junction_labels(self):
        with temp_dir() as tmp:
            b = SubgraphBuilder()
            segs = b.chain("chr5", [0, 1000, 2000, 3000, 4000], [100] * 4)
            b.variant(segs[0][1], segs[2][0], 50, hap="A", aid="ja")
            b.variant(segs[1][1], segs[3][0], 50, hap="B", aid="jb")
            b.variant(segs[2][1], segs[1][0], 50, hap="mixed", aid="jm")
            b.variant(segs[3][1], segs[3][1], 50, hap="mixed", aid="jself")
            directory = b.write(tmp)
            hard = read_subgraph(directory)
            soft = read_subgraph(directory, soft_junction_hap=True)
        self.assertEqual(hard.adjacencies["ja"].allowed, (Phasing.AA,))
        self.assertEqual(hard.adjacencies["jb"].allowed, (Phasing.BB,))
        self.assertEqual(set(hard.adjacencies["jm"].allowed), {Phasing.AB, Phasing.BA})
        self.assertEqual(set(hard.adjacencies["jself"].allowed), {Phasing.AA, Phasing.BB})
        self.assertEqual(set(soft.adjacencies["ja"].allowed), {Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB})
        self.assertEqual(soft.adjacencies["ja"].preferred, (Phasing.AA,))


class TestAmplicon0Fixture(unittest.TestCase):
    """demo3 amplicon_0 from amplicon-allele-analysis (chr4 + chr6 ecDNA)."""

    @classmethod
    def setUpClass(cls):
        cls.g = read_subgraph(os.path.join(DATA_DIR, "amplicon_0"))

    def test_counts(self):
        g = self.g
        self.assertEqual(len(g.segments), 18)  # 28 reference edges minus 10 phantom 1 bp segments
        kinds = {}
        for a in g.adjacencies.values():
            kinds[a.kind] = kinds.get(a.kind, 0) + 1
        self.assertEqual(kinds, {REFERENCE: 20, NOVEL: 14, IMPUTED: 12})
        self.assertEqual(len(g.nodes), 40)
        self.assertEqual(g.validate(), [])

    def test_terminals_are_region_flanks(self):
        self.assertEqual(set(self.g.terminals()), {"chr4:6019000-_+", "imputed_ext_chr4_4_-", "chr6:69527000-_+", "chr6:70856999+_-"})

    def test_haplotype_crossings(self):
        g = self.g
        for aid in ("sim00008", "sim00009", "sim00010", "sim00022"):
            self.assertEqual(set(g.adjacencies[aid].allowed), {Phasing.AB, Phasing.BA}, aid)
        for aid in ("sim00005", "sim00006", "sim00007", "sim00011", "sim00012"):
            self.assertEqual(g.adjacencies[aid].allowed, (Phasing.AA,), aid)
        # sim00008 joins the same two extremities as the reference adjacency it crosses
        refs = [a for a in g.adjacencies.values() if a.kind == REFERENCE and {a.u, a.v} == {g.adjacencies["sim00008"].u, g.adjacencies["sim00008"].v}]
        self.assertEqual(len(refs), 1)

    def test_parallel_imputed_edges_merged(self):
        for aid in ("imputed_adj_inter_chr4_2", "imputed_adj_anchor_chr4_3", "imputed_adj_inter_chr4_5", "imputed_adj_inter_chr4_8",
                    "imputed_adj_inter_chr6_24", "imputed_adj_anchor_chr6_25", "imputed_adj_inter_chr6_27"):
            self.assertNotIn(aid, self.g.adjacencies)
        for aid in ("imputed_adj_intra_chr4_1", "imputed_adj_anchor_chr4_6", "imputed_adj_intra_chr6_29"):
            self.assertEqual(self.g.adjacencies[aid].kind, IMPUTED)

    def test_anchor_on_phantom_moved_to_true_head(self):
        a = self.g.adjacencies["imputed_adj_anchor_chr4_6"]
        self.assertEqual({a.u, a.v}, {"imputed_ext_chr4_1_-", "chr4:6344999+_+"})


class TestRegions(unittest.TestCase):
    def test_distant_segments_on_one_chromosome_are_separate_regions(self):
        g = read_subgraph(os.path.join(DATA_DIR, "demo4", "amplicon_0"))
        chr16 = [r for r in g.regions() if r[0] == "chr16"]
        self.assertEqual(len(chr16), 2)  # 25.4 Mb and 46.7 Mb, not one 21 Mb span
        self.assertTrue(all(e - s < 1_000_000 for _, s, e in chr16))

    def test_demo3_regions(self):
        g = read_subgraph(os.path.join(DATA_DIR, "amplicon_0"))
        self.assertEqual([r[0] for r in g.regions()], ["chr4", "chr6"])


if __name__ == "__main__":
    unittest.main()

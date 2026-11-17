"""RCK-style evaluation (rck.ecdna.rck_eval) on hand-made cycles.

Run from the RCK root: python -m unittest tests.ecdna.test_rck_eval
"""
import unittest

from rck.ecdna import aa_format as aaf
from rck.ecdna import rck_eval as re_
from rck.ecdna.evaluate import CStep, Cycle


def cyc(name, weight, *steps, chrom="chr1"):
    """A circular cycle over [(start, end, hap)] 0-based inclusive on chrom, forward."""
    return Cycle(name, weight, [[CStep(chrom, s, e, h, "+") for s, e, h in steps]])


class Violations(unittest.TestCase):
    def test_twin_junction_on_both_homologs(self):
        """One novel adjacency labelled AA in one species and BB in another (an allelic twin): RCK Fig. 4A/4C."""
        js = re_.novel_adjacencies([cyc("a", 5, (0, 99, "A"), (300, 399, "A")), cyc("b", 2, (0, 99, "B"), (300, 399, "B"))])
        got = re_.is_violations(js)
        self.assertEqual(got["is_hom_same"], 2)  # both joins of the cycle: 99+ -> 300-, 399+ -> 0-
        self.assertEqual(got["is_hom"], 4)       # each of their 4 extremities used on A and on B
        self.assertEqual(got["is_ext"], 0)       # unlabelled, each extremity has one adjacency
        self.assertEqual(got["n_adjacencies"], 2)

    def test_sibling_with_own_breakpoints_is_clean(self):
        """Same locus, other allele, different breakpoints: no violation."""
        js = re_.novel_adjacencies([cyc("a", 5, (0, 99, "A"), (300, 399, "A")), cyc("b", 2, (20, 149, "B"), (320, 449, "B"))])
        self.assertEqual({k: v for k, v in re_.is_violations(js).items() if k.startswith("is_")},
                         {"is_ext": 0, "is_hom_same": 0, "is_hom": 0, "is_recip": 0})

    def test_reciprocal_extremities_on_opposite_homologs(self):
        """99+ used on A and 100- (its reference neighbour) used on B: RCK Fig. 4E."""
        js = re_.novel_adjacencies([cyc("a", 5, (0, 99, "A"), (300, 399, "A")), cyc("b", 2, (100, 199, "B"), (500, 599, "B"))])
        got = re_.is_violations(js)
        self.assertEqual(got["is_recip"], 1)
        self.assertEqual(got["is_hom"], 0)

    def test_reciprocal_on_one_homolog_is_allowed(self):
        js = re_.novel_adjacencies([cyc("a", 5, (0, 99, "A"), (300, 399, "A")), cyc("b", 2, (100, 199, "A"), (500, 599, "A"))])
        self.assertEqual(re_.is_violations(js)["is_recip"], 0)

    def test_extremity_in_two_adjacencies(self):
        js = re_.novel_adjacencies([cyc("a", 5, (0, 99, "A"), (300, 399, "A")), cyc("b", 2, (0, 99, "A"), (500, 599, "A"))])
        got = re_.is_violations(js)
        self.assertEqual(got["is_ext"], 2)  # 99+ (to 300- and to 500-) and 0- (from 399+ and from 599+)
        self.assertEqual(got["is_hom"], 0)

    def test_unphased_counts_for_extremity_exclusivity_only(self):
        js = re_.novel_adjacencies([cyc("a", 5, (0, 99, "."), (300, 399, ".")), cyc("b", 2, (0, 99, "."), (500, 599, "."))])
        got = re_.is_violations(js)
        self.assertEqual((got["is_ext"], got["is_hom"], got["is_hom_same"]), (2, 0, 0))


class Scores(unittest.TestCase):
    def test_adjacency_rates_up_to_symmetry(self):
        """The answer has the right junctions on the swapped allele: phased FNR 1 absolute, 0 within the block."""
        truth = [cyc("t", 5, (0, 99, "A"), (300, 399, "A"))]
        answer = [cyc("c", 5, (0, 99, "B"), (300, 399, "B"))]
        got = re_.score(answer, truth, blocks=[("chr1", 0, 1000, "b1")])
        self.assertEqual(got["adj_blind"]["fnr"], 0.0)
        self.assertEqual(got["adj_phased_abs"]["fnr"], 1.0)
        self.assertEqual(got["adj_phased_abs"]["fdr"], 1.0)
        self.assertEqual(got["adj_phased"]["fnr"], 0.0)
        self.assertEqual(got["n_blocks_swapped"], 1)
        self.assertEqual(got["cn_dist"], 0.0)
        self.assertGreater(got["cn_dist_abs"], 0.0)

    def test_block_boundary_switch_is_not_scored(self):
        """A traversal across a phase-block boundary, the answer's second block swapped: no adjacency, no violation
        within blocks (the phase across blocks is unidentifiable)."""
        truth = [cyc("t", 5, (0, 99, "A"), (100, 199, "A"), (300, 399, "A"))]
        answer = [cyc("c", 5, (0, 99, "A"), (100, 199, "B"), (300, 399, "B"))]
        blocks = [("chr1", 0, 99, "b1"), ("chr1", 100, 1000, "b2")]
        got = re_.score(answer, truth, blocks=blocks)
        self.assertEqual(got["n_blocks_swapped"], 1)
        self.assertEqual(got["is"]["is_recip"], 0)
        self.assertEqual(got["adj_phased"]["fnr"], 0.0)
        self.assertEqual(got["adj_phased"]["fdr"], 0.0)
        self.assertEqual(got["is_abs"]["is_recip"], 1)  # the absolute frame still sees the A-to-B join at 99|100

    def test_switch_inside_a_block_is_scored(self):
        """Post-hoc phasing that switches allele at a reference join inside one block: a real violation."""
        truth = [cyc("t1", 5, (0, 99, "A"), (100, 199, "A"), (300, 399, "A"))]
        answer = [cyc("c", 5, (0, 99, "A"), (100, 199, "B"), (300, 399, "A"))]
        got = re_.score(answer, truth, blocks=[("chr1", 0, 1000, "b1")])
        self.assertEqual(got["is"]["is_recip"], 1)
        self.assertGreater(got["adj_phased"]["fdr"], 0.0)

    def test_eq11_distance_one_allele_per_segment(self):
        """Truth 6 A + 2 B copies on 1000 bp, answer 8 on A: |8-6| + |0-2| = 4 per bp, relative 4/8."""
        truth = [cyc("t1", 6, (0, 999, "A")), cyc("t2", 2, (0, 999, "B"))]
        got = re_.score([cyc("m", 8, (0, 999, "A"))], truth, blocks=[("chr1", 0, 999, "b1")])
        self.assertAlmostEqual(got["cn_dist"], 4000.0)
        self.assertAlmostEqual(got["cn_dist_rel"], 0.5)

    def test_unphased_copies_are_all_distance(self):
        truth = [cyc("t1", 6, (0, 999, "A"))]
        got = re_.score([cyc("m", 6, (0, 999, "."))], truth)
        self.assertAlmostEqual(got["cn_dist_abs_rel"], 2.0)  # 6 missing on A + 6 unphased, over 6

    def test_candidate_fpr(self):
        """Graph: the 2 true joins and 2 spurious ones (AA positions are walk positions + 1); answer uses one spurious."""
        truth = [cyc("t", 5, (0, 99, "A"), (300, 399, "A"))]
        answer = [cyc("c", 5, (0, 99, "A"), (300, 399, "A")), cyc("x", 1, (600, 699, "A"), (800, 899, "A"))]
        edges = [aaf.AAEdge("discordant", ("chr1", 100, "+"), ("chr1", 301, "-"), 5),
                 aaf.AAEdge("discordant", ("chr1", 400, "+"), ("chr1", 1, "-"), 5),
                 aaf.AAEdge("discordant", ("chr1", 700, "+"), ("chr1", 801, "-"), 1),
                 aaf.AAEdge("discordant", ("chr1", 2000, "+"), ("chr1", 3001, "-"), 1)]
        got = re_.candidate_fpr(truth, answer, aaf.AAGraph("g", edges=edges))
        self.assertEqual((got["n_candidates"], got["n_spurious"], got["n_spurious_used"]), (4, 2, 1))
        self.assertAlmostEqual(got["cand_fpr"], 0.5)


if __name__ == "__main__":
    unittest.main()

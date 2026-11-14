"""The per-cell species matrix (J_seg, J_adj) of rck.ecdna.species_matrix, on small hand-built cases.
Definitions: docs/pipeline/evaluation.md, section 4.7, "Per-cell species matrix".

Run from the RCK root: python -m unittest tests.ecdna.test_species_matrix
"""
import math
import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna import aa_eval
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import CStep, Cycle
from rck.ecdna.species_matrix import (MATRIX_COLUMNS, column_items, item_jaccard, pair_jaccard, plot_species_matrix,
                                      species_matrix, weighted_jaccard)
from tests.ecdna.test_aa_partition import cycles_text

MB = 1_000_000


def walk(*segs, hap="A", chrom="chr1"):
    return [CStep(chrom, s, e, hap, "+") for s, e in segs]


def cyc(name, weight, steps, circular=True):
    return Cycle(name, weight, [steps], circular=circular)


# The worked example. Outer O (2 Mb): [0, 1.5 Mb) then [2, 2.5 Mb); novel adjacencies j1 = 1.5 Mb | 2 Mb
# and j2 = 2.5 Mb | 0 (wrap-around). Nested I (0.5 Mb), inside O in base space: [0.2, 0.45 Mb) then
# [2.1, 2.35 Mb); novel adjacencies j3 = 0.45 Mb | 2.1 Mb and j4 = 2.35 Mb | 0.2 Mb. The merged walk M is
# their Euler union (2.5 Mb): it uses j1..j4 and covers I's bases twice.
OUTER = walk((0, MB - 1), (MB, int(1.5 * MB) - 1), (2 * MB, int(2.5 * MB) - 1))
INNER = walk((200_000, 449_999), (2_100_000, 2_349_999))
MERGED = walk((0, 449_999), (2_100_000, 2_349_999), (200_000, int(1.5 * MB) - 1), (2 * MB, int(2.5 * MB) - 1))


def truth_pair():
    return [cyc("O", 8.0, OUTER), cyc("I", 4.0, INNER)]


def cell(frame, truth, inferred, haplotype="blind"):
    r = frame[(frame.truth == truth) & (frame.inferred == inferred) & (frame.haplotype == haplotype)]
    assert len(r) == 1, (truth, inferred, haplotype)
    return r.iloc[0]


def nan(x):
    return isinstance(x, float) and math.isnan(x)


class TestWorkedExample(unittest.TestCase):
    def test_lengths(self):
        self.assertEqual(aa_eval.cycle_length(OUTER), 2 * MB)
        self.assertEqual(aa_eval.cycle_length(INNER), MB // 2)
        self.assertEqual(aa_eval.cycle_length(MERGED), int(2.5 * MB))
        self.assertEqual(len(aa_eval.walk_adjacencies(cyc("m", 1, MERGED), False)), 4)

    def test_correct_answer(self):
        f = species_matrix(truth_pair(), truth_pair(), [], phased=False)
        self.assertEqual(list(f.columns), MATRIX_COLUMNS)
        self.assertEqual((cell(f, "O", "O").j_seg, cell(f, "O", "O").j_adj), (1.0, 1.0))
        self.assertEqual((cell(f, "I", "I").j_seg, cell(f, "I", "I").j_adj), (1.0, 1.0))
        self.assertAlmostEqual(cell(f, "O", "I").j_seg, 0.25)
        self.assertEqual(cell(f, "O", "I").j_adj, 0.0)
        self.assertAlmostEqual(cell(f, "I", "O").j_seg, 0.25)
        self.assertEqual(cell(f, "I", "O").j_adj, 0.0)
        self.assertEqual(cell(f, "O", "O").exact_cycle, 1)
        self.assertEqual(cell(f, "O", "I").exact_cycle, 0)
        # rows and columns by copy number, high to low: the correct answer is the diagonal
        diag = f[f.truth_rank == f.column_rank]
        self.assertTrue(((diag.j_seg == 1.0) & (diag.j_adj == 1.0)).all())
        self.assertFalse(f.merge_column.any() or f.split_row.any())

    def test_merge(self):
        f = species_matrix(truth_pair(), [cyc("M", 4.0, MERGED)], [], phased=False)
        self.assertAlmostEqual(cell(f, "O", "M").j_seg, 0.8)
        self.assertAlmostEqual(cell(f, "O", "M").j_adj, 0.5)
        self.assertAlmostEqual(cell(f, "I", "M").j_seg, 0.2)
        self.assertAlmostEqual(cell(f, "I", "M").j_adj, 0.5)
        self.assertTrue(f.merge_column.all())  # one column, two rows with J_adj > 0
        self.assertFalse(f.split_row.any())

    def test_phased_equals_blind_on_one_haplotype(self):
        f = species_matrix(truth_pair(), [cyc("M", 4.0, MERGED)], [], phased=True)
        self.assertEqual(set(f.haplotype), {"blind", "phased"})
        for t in ("O", "I"):
            self.assertAlmostEqual(cell(f, t, "M", "phased").j_seg, cell(f, t, "M").j_seg)
            self.assertAlmostEqual(cell(f, t, "M", "phased").j_adj, cell(f, t, "M").j_adj)


class TestPatterns(unittest.TestCase):
    def test_split(self):
        # the truth is the single figure-8 walk M; the answer reports it as O and I
        f = species_matrix([cyc("M", 4.0, MERGED)], truth_pair(), [], phased=False)
        self.assertAlmostEqual(cell(f, "M", "O").j_adj, 0.5)
        self.assertAlmostEqual(cell(f, "M", "I").j_adj, 0.5)
        self.assertTrue(f.split_row.all())
        self.assertFalse(f.merge_column.any())

    def test_single_junction_species(self):
        # a simple circle: its only novel adjacency is the closing one. The item form is defined
        # (1); M2's pair form is not (no pair), so it is NaN.
        circle = walk((0, 99_999))
        f = species_matrix([cyc("T", 5.0, circle)], [cyc("C", 5.0, circle)], [], phased=False)
        r = cell(f, "T", "C")
        self.assertEqual((r.n_items_truth, r.n_items_inferred, r.j_seg, r.j_adj), (1, 1, 1.0, 1.0))
        self.assertTrue(nan(r.j_adj_pairs))

    def test_empty_union_is_nan(self):
        # two linear one-segment walks have no novel adjacency: J_adj is undefined, never 1
        line = walk((0, 99_999))
        f = species_matrix([cyc("T", 5.0, line, circular=False)], [cyc("C", 5.0, line, circular=False)], [], phased=False)
        r = cell(f, "T", "C")
        self.assertEqual((r.n_items_truth, r.n_items_inferred), (0, 0))
        self.assertTrue(nan(r.j_adj))
        self.assertEqual(r.j_seg, 1.0)
        self.assertTrue(nan(item_jaccard(set(), set())))
        self.assertEqual(item_jaccard({0}, set()), 0.0)

    def test_pair_form(self):
        self.assertAlmostEqual(pair_jaccard({0, 1}, {0, 1, 2, 3}), 1 / 6)  # C(2,2) / (C(2,2) + C(4,2) - C(2,2))
        f = species_matrix(truth_pair(), [cyc("M", 4.0, MERGED)], [], phased=False)
        self.assertAlmostEqual(cell(f, "O", "M").j_adj_pairs, 1 / 6)

    def test_copy_weighted(self):
        # same walk, copy numbers 8 and 6: 6/8 on every base
        self.assertAlmostEqual(weighted_jaccard(OUTER, 8.0, OUTER, 6.0, False), 0.75)
        f = species_matrix(truth_pair(), [cyc("O", 6.0, OUTER), cyc("I", 4.0, INNER)], [], phased=False)
        self.assertAlmostEqual(cell(f, "O", "O").j_seg_cn, 0.75)
        self.assertEqual(cell(f, "O", "O").j_seg, 1.0)  # J_seg itself ignores copy number


class TestLowWeight(unittest.TestCase):
    """The 6:1 failure: the outer species alone is strict; the nested species sits in low_weight walks."""

    def setUp(self):
        self.truth = [cyc("O", 36.0, OUTER), cyc("I", 6.0, INNER)]
        self.strict = [cyc("ec1", 36.0, OUTER)]
        # two 3-copy leftovers (below 0.1 x 36): the nested circle, and the same circle read from another start
        self.low = [cyc("ec2", 3.0, INNER), cyc("ec3", 3.0, walk((2_100_000, 2_349_999), (200_000, 449_999)))]

    def test_columns_and_items(self):
        f = species_matrix(self.truth, self.strict, self.low, phased=False)
        cols = f.drop_duplicates("column_rank").sort_values("column_rank")
        self.assertEqual(list(cols.inferred), ["ec1", "ec2", "ec3"])
        self.assertEqual(list(cols.low_weight), [False, True, True])
        # the leftovers carry the nested species' adjacencies: (1, 1) on the I row, though low_weight
        self.assertEqual((cell(f, "I", "ec2").j_seg, cell(f, "I", "ec2").j_adj), (1.0, 1.0))
        # low_weight columns do not make a split
        self.assertFalse(f.split_row.any())

    def test_strict_sets_are_m2s(self):
        t0 = [aa_eval.strip_haplotypes(t) for t in self.truth]
        s0 = [aa_eval.strip_haplotypes(c) for c in self.strict]
        o0 = [aa_eval.strip_haplotypes(c) for c in self.low]
        items, ts, ss, _ = aa_eval.adjacency_items(t0, s0, 10_000, False)
        items2, ts2, ss2, os2 = column_items(t0, s0, o0, 10_000, False)
        self.assertEqual((ts, ss), (ts2, ss2))
        self.assertEqual(items2[:len(items)], items)
        self.assertEqual(os2[0], ts[1])  # the leftover took the nested species' unmatched truth items

    def test_low_weight_does_not_take_a_strict_item(self):
        # a leftover adjacency nearer to a truth item than the strict one: the strict set stays M2's
        near = walk((0, MB - 1), (MB, int(1.5 * MB) - 1), (2 * MB + 5_000, int(2.5 * MB) - 1))  # j1 shifted 5 kb
        exact = walk((0, MB - 1), (MB, int(1.5 * MB) - 1), (2 * MB + 1, int(2.5 * MB) - 1))  # j1 shifted 1 bp
        f = species_matrix(self.truth, [cyc("ec1", 36.0, near)], [cyc("lw", 1.0, exact)], phased=False)
        self.assertEqual(cell(f, "O", "ec1").j_adj, 1.0)
        self.assertLess(cell(f, "O", "lw").j_adj, 1.0)  # its j1 becomes a false item


class TestAgainstRecoveryMatrix(unittest.TestCase):
    def test_j_seg_is_the_matched_pairs_jaccard(self):
        truth = [cyc("O", 8.0, OUTER), cyc("I", 4.0, walk((200_000, 449_999), (2_100_000, 2_349_999), hap="B"))]
        answer = [cyc("c1", 7.5, walk((0, MB - 1), (MB, int(1.5 * MB) - 1), (2 * MB, int(2.4 * MB) - 1))),
                  cyc("c2", 4.2, walk((200_000, 449_999), hap="B") + walk((2_100_000, 2_349_999), hap="A"))]
        f = species_matrix(truth, answer, [], phased=True)
        for s in ev.score_decomposition(truth, answer, 10_000):
            self.assertIsNotNone(s.inferred)
            self.assertAlmostEqual(cell(f, s.truth, s.inferred, "phased").j_seg, s.jaccard_hap)
        t0 = [aa_eval.strip_haplotypes(t) for t in truth]
        a0 = [aa_eval.strip_haplotypes(c) for c in answer]
        for s in ev.score_decomposition(t0, a0, 10_000):
            self.assertAlmostEqual(cell(f, s.truth, s.inferred, "blind").j_seg, s.jaccard_total)


class TestEndToEnd(unittest.TestCase):
    """aa_eval.score writes aa_species_matrix.tsv; the low_weight cycle of the file is a column."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_score_writes_the_matrix(self):
        from tests.ecdna.test_aa_partition import S1, S2, S3, S4, S5, truth_pair as pair

        segs = [("chr1", s + 1, e + 1) for s, e in (S1, S2, S3, S4, S5)]
        outer = "Copy_count=8.0;Segments=1+,2+,3+,5+;Haplotypes=A,A,A,A;Status=species;Species=ec1"
        inner = "Copy_count=1.0;Segments=2+,4+;Haplotypes=A,A;Status=low_weight;Species=ec2"
        cycles = os.path.join(self.dir, "a_cycles.txt")
        with open(cycles, "w") as handle:
            handle.write(cycles_text([("chr1", 1, 6000)], segs, [outer, inner]))
        result = aa_eval.score(cycles, pair(), case="a")
        m = result.matrix
        self.assertEqual(list(m.columns), MATRIX_COLUMNS)
        self.assertEqual(set(m.haplotype), {"blind", "phased"})
        self.assertEqual(set(m.case), {"a"})
        self.assertEqual(len(m), 2 * 2 * 2)  # 2 haplotype modes x 2 truth x 2 columns
        self.assertTrue(m[m.inferred == "ec2"].low_weight.all())
        self.assertEqual((cell(m, "outer", "ec1").j_seg, cell(m, "outer", "ec1").j_adj), (1.0, 1.0))
        self.assertEqual(cell(m, "inner", "ec2").j_adj, 1.0)
        out = os.path.join(self.dir, "out")
        paths = aa_eval.write(out, result)
        back = pd.read_csv(paths["matrix"], sep="\t")
        self.assertEqual(len(back), len(m))
        self.assertNotIn("j_seg", pd.read_csv(paths["species"], sep="\t").columns)  # existing tables unchanged
        self.assertEqual(list(pd.read_csv(paths["partition"], sep="\t").columns), aa_eval.PARTITION_COLUMNS)
        png = plot_species_matrix(m, os.path.join(out, "m.png"), "blind", "a")
        if png is not None:
            self.assertTrue(os.path.getsize(png) > 0)

    def test_chromosomal_paths_are_not_species(self):
        """The multi-clone route's telomere-to-telomere paths (Status=chromosomal) are left out of
        M1, M2 and the matrix, as evaluate.read_inferred leaves them out of the composite."""
        from tests.ecdna.test_aa_partition import S1, S2, S3, S4, S5, truth_pair as pair

        segs = [("chr1", s + 1, e + 1) for s, e in (S1, S2, S3, S4, S5)]
        outer = "Copy_count=8.0;Segments=1+,2+,3+,5+;Haplotypes=A,A,A,A;Status=species;Species=ec1"
        inner = "Copy_count=4.0;Segments=2+,4+;Haplotypes=A,A;Status=species;Species=ec2"
        chrom = "Copy_count=2.0;Segments=0+,1+,2+,3+,4+,5+,0-;Haplotypes=.,B,B,B,B,B,.;Status=chromosomal;Species=path1"
        cycles = os.path.join(self.dir, "mc_cycles.txt")
        with open(cycles, "w") as handle:
            handle.write(cycles_text([("chr1", 1, 6000)], segs, [outer, inner, chrom]))
        result = aa_eval.score(cycles, pair(), case="mc")
        self.assertNotIn("path1", set(result.matrix.inferred))
        self.assertEqual({c.name for c in result.strict}, {"ec1", "ec2"})
        self.assertNotIn("path1", {c.name for c in result.candidates})
        m1 = result.partition[(result.partition["metric"] == "species_detection") & (result.partition["mode"] == "strict")]
        self.assertTrue(len(m1))
        self.assertTrue((m1["n_predicted"] == 2).all())  # ec1 and ec2; the chromosome is not a predicted species


if __name__ == "__main__":
    unittest.main()

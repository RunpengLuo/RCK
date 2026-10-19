import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna.evaluate import (CStep, Cycle, evaluate, jaccard, junction_scores, order_match, read_truth_table, runs,
                                truth_from_simulation, write_truth_table)
from tests.ecdna.helpers import DATA_DIR

TRUTH = [CStep("chr1", 0, 999, "A", "+"), CStep("chr1", 5000, 5999, "A", "-"), CStep("chr2", 100, 1099, "B", "+")]


def rotate(steps, k):
    return steps[k:] + steps[:k]


def reverse(steps):
    return [CStep(s.chrom, s.start, s.end, s.hap, "-" if s.strand == "+" else "+") for s in reversed(steps)]


class TestMetrics(unittest.TestCase):
    def test_identical(self):
        self.assertEqual(jaccard(TRUTH, TRUTH), 1.0)
        self.assertEqual(junction_scores(TRUTH, TRUTH, tol=0, phased=True), (1.0, 1.0))
        self.assertTrue(order_match(TRUTH, TRUTH, tol=0))

    def test_rotation_and_reversal(self):
        self.assertTrue(order_match(TRUTH, rotate(TRUTH, 1), tol=0))
        self.assertTrue(order_match(TRUTH, reverse(TRUTH), tol=0))
        self.assertEqual(junction_scores(TRUTH, reverse(TRUTH), tol=0, phased=True), (1.0, 1.0))

    def test_finer_segmentation_is_equivalent(self):
        split = [CStep("chr1", 0, 499, "A", "+"), CStep("chr1", 500, 999, "A", "+")] + TRUTH[1:]
        self.assertEqual(len(runs(split)), 3)
        self.assertEqual(jaccard(TRUTH, split), 1.0)
        self.assertTrue(order_match(TRUTH, split, tol=0))

    def test_haplotype_errors(self):
        wrong = [CStep(s.chrom, s.start, s.end, "A", s.strand) for s in TRUTH]
        self.assertLess(jaccard(TRUTH, wrong, with_hap=True), 1.0)
        self.assertEqual(jaccard(TRUTH, wrong, with_hap=False), 1.0)
        self.assertEqual(junction_scores(TRUTH, wrong, tol=0, phased=False), (1.0, 1.0))
        self.assertLess(junction_scores(TRUTH, wrong, tol=0, phased=True)[0], 1.0)
        self.assertFalse(order_match(TRUTH, wrong, tol=0))

    def test_tolerance(self):
        shifted = [CStep(s.chrom, s.start + 3, s.end + 3, s.hap, s.strand) for s in TRUTH]
        self.assertFalse(order_match(TRUTH, shifted, tol=0))
        self.assertTrue(order_match(TRUTH, shifted, tol=5))
        self.assertEqual(junction_scores(TRUTH, shifted, tol=5, phased=True), (1.0, 1.0))

    def test_multiplicity_counts_in_jaccard(self):
        doubled = TRUTH + [CStep("chr1", 0, 999, "A", "+")]
        self.assertAlmostEqual(jaccard(TRUTH, doubled), 3000 / 4000.0)


class TestTables(unittest.TestCase):
    def test_truth_table_round_trip(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "truth.tsv")
            write_truth_table(path, [Cycle("ec", 12.0, [TRUTH])])
            back = read_truth_table(path)
        self.assertEqual(back[0].steps, TRUTH)
        self.assertEqual(back[0].weight, 12.0)

    def test_evaluate_decompositions_table(self):
        def walk(steps):
            return ",".join("{c}:{s}-{e}:{h}:{st}".format(c=x.chrom, s=x.start, e=x.end, h=x.hap, st=x.strand) for x in steps)

        rows = [
            {"decomposition_rank": 1, "objective": 1.0, "species": "ec1", "weight": 11, "ordering_rank": 1, "walk": walk(rotate(TRUTH, 2))},
            {"decomposition_rank": 2, "objective": 2.0, "species": "ec1", "weight": 6, "ordering_rank": 1, "walk": walk(TRUTH[:1] + TRUTH[2:])},
            {"decomposition_rank": 2, "objective": 2.0, "species": "ec2", "weight": 6, "ordering_rank": 1, "walk": walk(TRUTH[1:2])},
        ]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "decompositions.tsv")
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            frame, summary = evaluate(path, [Cycle("truth1", 12.0, [TRUTH])], tol=0)
        self.assertEqual(summary["first_rank_with_all_orders_matched"], 1)
        self.assertEqual(summary["best_mean_jaccard_rank"], 1)
        top = frame[frame.decomposition_rank == 1].iloc[0]
        self.assertAlmostEqual(top.cn_rel_err, 1 / 12.0)
        self.assertFalse(bool(frame[frame.decomposition_rank == 2].iloc[0].order_match))


class TestSimulatorTruth(unittest.TestCase):
    """demo3_clone1_ecdna: one episome species, 14 segments, 16 steps, phi 16."""

    @classmethod
    def setUpClass(cls):
        cls.cycles = truth_from_simulation(os.path.join(DATA_DIR, "demo3_truth"))

    def test_species(self):
        self.assertEqual(len(self.cycles), 1)
        c = self.cycles[0]
        self.assertEqual(c.name, "ecdna1")
        self.assertEqual(c.weight, 16.0)
        self.assertEqual(len(c.steps), 16)
        self.assertEqual(sum(s.end - s.start + 1 for s in c.steps), 2_343_000)

    def test_haplotypes(self):
        steps = self.cycles[0].steps
        seg4 = [s for s in steps if (s.chrom, s.start) == ("chr4", 6344000)]
        self.assertEqual(sorted(s.hap for s in seg4), ["A", "B"])
        seg5 = [s for s in steps if (s.chrom, s.start) == ("chr4", 6345000)]
        self.assertEqual([s.hap for s in seg5], ["A", "A"])
        seg10 = [s for s in steps if (s.chrom, s.start) == ("chr6", 70090000)]
        self.assertEqual([s.hap for s in seg10], ["B"])
        others = [s for s in steps if s not in seg4 + seg10]
        self.assertTrue(all(s.hap == "A" for s in others))

    def test_truth_matches_itself(self):
        c = self.cycles[0]
        self.assertTrue(order_match(c.steps, rotate(c.steps, 5), tol=0))
        # seg1 / seg2 and seg12 / seg13 are joined by reference joins, so 16 steps form 14 runs
        self.assertEqual(len(runs(c.steps)), 14)


if __name__ == "__main__":
    unittest.main()

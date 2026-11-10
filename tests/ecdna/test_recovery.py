import json
import os
import unittest

import pandas as pd

from rck.ecdna.evaluate import RECOVERY_COLUMNS, CStep, Cycle, evaluate, score_pair, write_evaluation
from tests.ecdna.helpers import temp_dir

TRUTH = [CStep("chr1", 0, 999, "A", "+"), CStep("chr1", 5000, 5999, "A", "-"), CStep("chr2", 100, 1099, "B", "+")]
SMALL = [CStep("chr3", 0, 999, "A", "+")]
FAR = [CStep("chr9", 0, 999, "A", "+")]


def walk(steps):
    return ",".join("{c}:{s}-{e}:{h}:{st}".format(c=x.chrom, s=x.start, e=x.end, h=x.hap, st=x.strand) for x in steps)


def swapped(steps):
    return [CStep(s.chrom, s.start, s.end, {"A": "B", "B": "A"}[s.hap], s.strand) for s in steps]


class TestPairScores(unittest.TestCase):
    def test_exact_recovery_scores_one(self):
        s = score_pair(Cycle("t", 12.0, [TRUTH]), Cycle("ec1", 12.0, [TRUTH]), tol=0)
        self.assertEqual((s.segment_recall, s.segment_precision, s.junction_f1, s.phased_junction_f1), (1.0, 1.0, 1.0, 1.0))
        self.assertEqual(s.phasing_accuracy, 1.0)
        self.assertAlmostEqual(s.composite, 1.0)

    def test_copy_number_error_lowers_the_composite(self):
        s = score_pair(Cycle("t", 12.0, [TRUTH]), Cycle("ec1", 9.0, [TRUTH]), tol=0)
        self.assertAlmostEqual(s.composite, (4 + 0.75) / 5)

    def test_phasing_accuracy(self):
        one_hap = [CStep(s.chrom, s.start, s.end, "A", s.strand) for s in TRUTH]
        s = score_pair(Cycle("t", 12.0, [TRUTH]), Cycle("ec1", 12.0, [one_hap]), tol=0)
        self.assertAlmostEqual(s.phasing_accuracy, 2 / 3)
        self.assertAlmostEqual(s.phasing_accuracy_swap_invariant, 2 / 3)
        self.assertEqual(s.segment_recall, 1.0)
        # every haplotype swapped: wrong as given, right up to the A/B labels
        s = score_pair(Cycle("t", 12.0, [TRUTH]), Cycle("ec1", 12.0, [swapped(TRUTH)]), tol=0)
        self.assertEqual(s.phasing_accuracy, 0.0)
        self.assertEqual(s.phasing_accuracy_swap_invariant, 1.0)
        self.assertEqual(s.junction_f1, 1.0)
        self.assertEqual(s.phased_junction_f1, 0.0)

    def test_partial_species(self):
        s = score_pair(Cycle("t", 12.0, [TRUTH]), Cycle("ec1", 12.0, [TRUTH[:2]]), tol=0)
        self.assertAlmostEqual(s.segment_recall, 2 / 3)
        self.assertEqual(s.segment_precision, 1.0)
        self.assertAlmostEqual(s.junction_recall, 1 / 3)
        self.assertAlmostEqual(s.junction_precision, 1 / 2)
        self.assertAlmostEqual(s.junction_f1, 0.4)


class TestRecoveryMatrix(unittest.TestCase):
    def setUp(self):
        rows = [
            {"decomposition_rank": 1, "objective": 1.0, "species": "ec1", "weight": 9, "ordering_rank": 1, "walk": walk(TRUTH)},
            {"decomposition_rank": 2, "objective": 1.5, "species": "ec1", "weight": 12, "ordering_rank": 1, "walk": walk(TRUTH)},
            {"decomposition_rank": 2, "objective": 1.5, "species": "ec2", "weight": 5, "ordering_rank": 1, "walk": walk(SMALL)},
        ]
        self.truth = [Cycle("truth1", 12.0, [TRUTH]), Cycle("truth2", 5.0, [SMALL]), Cycle("truth3", 3.0, [FAR])]
        self.regions = [("chr1", 0, 10_000), ("chr2", 0, 2_000), ("chr3", 0, 2_000)]
        self.tmp = temp_dir()
        self.dir = self.tmp.name
        self.path = os.path.join(self.dir, "decompositions.tsv")
        pd.DataFrame(rows).to_csv(self.path, sep="\t", index=False)

    def tearDown(self):
        self.tmp.cleanup()

    def test_one_row_per_truth_species(self):
        _, summary = evaluate(self.path, self.truth, tol=0, regions=self.regions)
        matrix = {r["truth"]: r for r in summary["recovery_matrix"]}
        self.assertEqual(list(matrix), ["truth1", "truth2", "truth3"])
        t1, t2, t3 = matrix["truth1"], matrix["truth2"], matrix["truth3"]
        # rank 1 recovers truth1 with 9 of 12 copies and misses truth2
        self.assertEqual(t1["inferred"], "ec1")
        self.assertAlmostEqual(t1["composite"], 0.95)
        self.assertEqual((t1["best_rank"], t1["best_composite"]), (2, 1.0))
        self.assertEqual(t1["truth_n_junctions"], 3)
        self.assertIsNone(t2["inferred"])
        self.assertEqual(t2["composite"], 0.0)
        self.assertEqual((t2["best_rank"], t2["best_inferred"], t2["best_composite"]), (2, "ec2", 1.0))
        # out of scope: listed, not scored, not in the mean
        self.assertFalse(t3["in_scope"])
        self.assertIsNone(t3["composite"])
        self.assertAlmostEqual(summary["mean_composite"], 0.475)

    def test_written_as_tsv_and_strict_json(self):
        frame, summary = evaluate(self.path, self.truth, tol=0, regions=self.regions)
        paths = write_evaluation(self.dir, frame, summary)
        table = pd.read_csv(paths["recovery_matrix"], sep="\t")
        self.assertEqual(list(table.columns), RECOVERY_COLUMNS)
        self.assertEqual(len(table), 3)

        def reject(constant):
            raise ValueError("non-standard JSON constant " + constant)

        with open(paths["summary"]) as handle:
            loaded = json.loads(handle.read(), parse_constant=reject)
        self.assertEqual(len(loaded["recovery_matrix"]), 3)

    def test_no_decompositions(self):
        empty = os.path.join(self.dir, "empty.tsv")
        pd.DataFrame(columns=["decomposition_rank", "objective", "species", "weight", "ordering_rank", "walk"]).to_csv(
            empty, sep="\t", index=False)
        _, summary = evaluate(empty, self.truth[:1], tol=0)
        row = summary["recovery_matrix"][0]
        self.assertEqual(row["composite"], 0.0)
        self.assertIsNone(row["inferred"])


if __name__ == "__main__":
    unittest.main()

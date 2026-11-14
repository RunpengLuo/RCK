"""bench_IW metrics (rck.ecdna.bench_metrics): ecDNA copy-number recall / precision and the ecDNA
proportion scores, on hand-made cycles.

Run from the RCK root: python -m unittest tests.ecdna.test_bench_metrics
"""
import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna import bench_metrics as bm
from rck.ecdna.evaluate import CStep, Cycle


def cyc(name, weight, *spans, chrom="chr1"):
    """A cycle over 0-based inclusive spans [(start, end)] on chrom, haplotype A, forward."""
    return Cycle(name, weight, [[CStep(chrom, s, e, "A", "+") for s, e in spans]])


class CopyNumber(unittest.TestCase):
    def test_exact_answer(self):
        truth = [cyc("ecdna1", 8.0, (0, 999)), cyc("ecdna2", 4.0, (500, 1499))]
        got = bm.cn_scores(truth, [cyc("a", 8.0, (0, 999)), cyc("b", 4.0, (500, 1499))])
        self.assertAlmostEqual(got["cn_recall"], 1.0)
        self.assertAlmostEqual(got["cn_precision"], 1.0)
        self.assertAlmostEqual(got["cn_within_tol"], 1.0)
        self.assertAlmostEqual(got["truth_mass"], 8 * 1000 + 4 * 1000)

    def test_merged_species_keeps_recall(self):
        """One cycle over both truth species' union at the summed copy number: right where they
        overlap, too high elsewhere; the profile metric scores copies, not species."""
        truth = [cyc("ecdna1", 6.0, (0, 999)), cyc("ecdna2", 2.0, (500, 999))]
        got = bm.cn_scores(truth, [cyc("m", 8.0, (0, 999))])
        # truth: 6 over 0-499, 8 over 500-999 -> mass 7000; inferred 8000; overlap 6*500 + 8*500
        self.assertAlmostEqual(got["cn_recall"], 7000 / 7000)
        self.assertAlmostEqual(got["cn_precision"], 7000 / 8000)
        self.assertAlmostEqual(got["cn_within_tol"], 0.5)  # 0-499: |8-6| = 2 > max(1, 0.6)

    def test_traversals_count(self):
        truth = [cyc("ecdna1", 5.0, (0, 99), (0, 99))]  # traversed twice: 10 copies
        got = bm.cn_scores(truth, [cyc("a", 5.0, (0, 99))])
        self.assertAlmostEqual(got["cn_recall"], 0.5)
        self.assertAlmostEqual(got["cn_precision"], 1.0)

    def test_missed_and_other_chromosome(self):
        truth = [cyc("ecdna1", 4.0, (0, 99)), cyc("ecdna2", 4.0, (0, 99), chrom="chr2")]
        got = bm.cn_scores(truth, [cyc("a", 4.0, (0, 99)), cyc("b", 4.0, (0, 99), chrom="chr3")])
        self.assertAlmostEqual(got["cn_recall"], 0.5)
        self.assertAlmostEqual(got["cn_precision"], 0.5)

    def test_empty_answer(self):
        got = bm.cn_scores([cyc("ecdna1", 4.0, (0, 99))], [])
        self.assertEqual(got["cn_recall"], 0.0)
        self.assertIsNone(got["cn_precision"])


class Proportions(unittest.TestCase):
    def test_matched_and_unmatched(self):
        truth = [cyc("ecdna1", 6.0, (0, 999)), cyc("ecdna2", 2.0, (5000, 5999))]
        inferred = [cyc("a", 5.0, (0, 999)), cyc("b", 5.0, (9000, 9999))]  # b matches nothing
        summary, rows = bm.proportion_scores(truth, inferred, "strict")
        # ecdna1 0.75 vs a 0.5; ecdna2 0.25 unmatched; b 0.5 unmatched
        self.assertAlmostEqual(summary["proportion_tv"], (0.25 + 0.25 + 0.5) / 2)
        self.assertAlmostEqual(summary["proportion_mae_matched"], 0.25)
        self.assertEqual(summary["n_matched"], 1)
        self.assertEqual(len(rows), 3)

    def test_perfect(self):
        truth = [cyc("ecdna1", 6.0, (0, 999)), cyc("ecdna2", 2.0, (5000, 5999))]
        summary, _ = bm.proportion_scores(truth, [cyc("x", 3.0, (0, 999)), cyc("y", 1.0, (5000, 5999))], "strict")
        self.assertAlmostEqual(summary["proportion_tv"], 0.0)


class AnswerFile(unittest.TestCase):
    """score() on an AA cycles file: strict leaves out low_weight, chromosomal is never a species,
    and the 'all' set keeps the low_weight mass."""

    def test_strict_and_all(self):
        lines = ["Interval\t1\tchr1\t1\t100000",
                 "List of cycle segments",
                 "Segment\t1\tchr1\t1\t1000",
                 "Segment\t2\tchr1\t5001\t6000",
                 "Segment\t3\tchr1\t1\t100000",
                 "Cycle=1;Copy_count=6.0;Segments=1+;Species=ec1;Status=ok",
                 "Cycle=2;Copy_count=1.0;Segments=2+;Species=ec2;Status=low_weight",
                 "Cycle=3;Copy_count=2.0;Segments=0+,3+,0-;Species=path1;Status=chromosomal"]
        truth = [cyc("ecdna1", 6.0, (0, 999)), cyc("ecdna2", 1.0, (5000, 5999))]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "x_cycles.txt")
            with open(path, "w") as handle:
                handle.write("\n".join(lines) + "\n")
            summary, rows = bm.score(path, truth)
        self.assertEqual(summary["n_cycles"], 2)
        self.assertEqual(summary["n_strict"], 1)
        self.assertAlmostEqual(summary["cn_recall"], 6000 / 7000)
        self.assertAlmostEqual(summary["cn_recall_all"], 1.0)
        self.assertAlmostEqual(summary["proportion_tv"], 1 / 7)  # truth 6/7 vs 1, plus ecdna2 1/7 missed: (1/7 + 1/7) / 2
        self.assertAlmostEqual(summary["proportion_tv_all"], 0.0)
        self.assertEqual(set(pd.DataFrame(rows)["set"]), {"strict", "all"})

    def test_multiclone_weights(self):
        with tempfile.TemporaryDirectory() as tmp:
            pd.DataFrame({"decomposition_rank": [1, 1], "species": ["ec1", "path1"], "psi_clone1": [10.0, 1.0],
                          "psi_clone2": [4.0, 1.0]}).to_csv(os.path.join(tmp, "decompositions.tsv"), sep="\t", index=False)
            got = bm.multiclone_bulk_weights(tmp, {"normal": 0.2, "clone1": 0.5, "clone2": 0.3})
        self.assertAlmostEqual(got["ec1"], 0.5 * 10 + 0.3 * 4)


if __name__ == "__main__":
    unittest.main()

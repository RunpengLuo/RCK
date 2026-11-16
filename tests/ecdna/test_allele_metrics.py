"""Allele-specific ecDNA copy numbers and post-hoc phasing (rck.ecdna.allele_metrics), on hand-made cycles.

Run from the RCK root: python -m unittest tests.ecdna.test_allele_metrics
"""
import os
import tempfile
import unittest

from rck.ecdna import aa_format as aaf
from rck.ecdna import allele_metrics as am
from rck.ecdna.evaluate import CStep, Cycle


def cyc(name, weight, *steps, chrom="chr1"):
    """A cycle over [(start, end, hap)] 0-based inclusive on chrom, forward."""
    return Cycle(name, weight, [[CStep(chrom, s, e, h, "+") for s, e, h in steps]])


class Scores(unittest.TestCase):
    def test_exact_phased_answer(self):
        truth = [cyc("t1", 6.0, (0, 999, "A")), cyc("t2", 2.0, (0, 999, "B"))]
        got = am.allele_scores(am.allele_table(truth, [cyc("a", 6.0, (0, 999, "A")), cyc("b", 2.0, (0, 999, "B"))]))
        for key in ("as_cn_recall", "as_cn_precision", "blind_cn_recall", "as_cn_within_tol", "major_allele_acc"):
            self.assertAlmostEqual(got[key], 1.0, msg=key)
        self.assertAlmostEqual(got["as_cn_mae"], 0.0)
        self.assertAlmostEqual(got["baf_mae"], 0.0)
        self.assertAlmostEqual(got["mixed_bp_frac"], 1.0)  # both alleles amplified on every base

    def test_unphased_answer_keeps_blind_recall_only(self):
        truth = [cyc("t1", 6.0, (0, 999, "A"))]
        got = am.allele_scores(am.allele_table(truth, [cyc("a", 6.0, (0, 999, "."))]))
        self.assertAlmostEqual(got["blind_cn_recall"], 1.0)
        self.assertAlmostEqual(got["as_cn_recall"], 0.0)
        self.assertAlmostEqual(got["as_cn_precision"], 0.0)
        self.assertAlmostEqual(got["unphased_mass_frac"], 1.0)
        self.assertIsNone(got["baf_mae"])  # nothing phased to compare

    def test_one_label_per_segment_loses_the_minor_allele(self):
        """Truth 6 copies on A and 2 on B of the same bases; all 8 called on A: recall 6/8, major right."""
        truth = [cyc("t1", 6.0, (0, 999, "A")), cyc("t2", 2.0, (0, 999, "B"))]
        got = am.allele_scores(am.allele_table(truth, [cyc("m", 8.0, (0, 999, "A"))]))
        self.assertAlmostEqual(got["blind_cn_recall"], 1.0)
        self.assertAlmostEqual(got["as_cn_recall"], 6 / 8)
        self.assertAlmostEqual(got["as_cn_mae"], 4.0)  # |8-6| + |0-2|
        self.assertAlmostEqual(got["major_allele_acc"], 1.0)
        self.assertAlmostEqual(got["baf_mae"], 0.25)

    def test_within_block_swaps_once_per_block(self):
        truth = [cyc("t1", 5.0, (0, 99, "A"), (200, 299, "A"))]
        answer = [cyc("a", 5.0, (0, 99, "B"), (200, 299, "A"))]
        table = am.allele_table(truth, answer)
        self.assertAlmostEqual(am.allele_scores(table)["as_cn_recall"], 0.5)
        swapped, n = am.within_block(table, [("chr1", 0, 149, "b1"), ("chr1", 150, 400, "b2")])
        self.assertEqual(n, 1)
        self.assertAlmostEqual(am.allele_scores(swapped)["as_cn_recall"], 1.0)

    def test_within_block_on_a_table_read_back_from_tsv(self):
        """allele_segments.tsv read back (string dtype, copy-on-write) swaps the same way."""
        truth = [cyc("t1", 7.0, (0, 99, "B"), (200, 299, "A"))]
        answer = [cyc("a", 7.0, (0, 99, "A"), (200, 299, "B"))]
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "allele_segments.tsv")
            am.allele_table(truth, answer).to_csv(path, sep="\t", index=False)
            import pandas as pd
            table = pd.read_csv(path, sep="\t")
        swapped, n = am.within_block(table, [("chr1", 0, 149, "b1"), ("chr1", 150, 400, "b2")])
        self.assertEqual(n, 2)
        self.assertAlmostEqual(am.allele_scores(swapped)["as_cn_recall"], 1.0)
        self.assertEqual(float(swapped["c_A"].sum() + swapped["c_B"].sum()), 14.0)


class Rephase(unittest.TestCase):
    def test_rules_and_file_rewrite(self):
        truth = [cyc("t1", 6.0, (0, 999, "A")), cyc("t2", 2.0, (0, 999, "B")), cyc("t3", 3.0, (2000, 2999, "B"))]
        rule = am.truth_major_rule(truth)
        self.assertEqual(rule("chr1", 0, 999), "A")
        self.assertEqual(rule("chr1", 2000, 2999), "B")
        self.assertEqual(rule("chr2", 0, 10), aaf.UNPHASED)
        with tempfile.TemporaryDirectory() as tmp:
            src, dst = os.path.join(tmp, "in_cycles.txt"), os.path.join(tmp, "out_cycles.txt")
            with open(src, "w") as handle:
                handle.write("Interval\t1\tchr1\t1\t3000\n")
                handle.write("Segment\t1\tchr1\t1\t1000\nSegment\t2\tchr1\t2001\t3000\n")
                handle.write("Cycle=1;Copy_count=8.0;Segments=1+,2-\n")
                handle.write("Cycle=2;Copy_count=1.0;Segments=0+,2+,0-;Path_constraints_satisfied=\n")
            counts = am.rephase(src, rule, dst)
            self.assertEqual(counts, {"A": 1, "B": 2, "filled": 0})
            recs = aaf.read_cycles(dst)
            self.assertEqual(recs.cycles[0].haplotypes, ["A", "B"])
            self.assertEqual(recs.cycles[1].haplotypes, ["B"])
            self.assertTrue(recs.cycles[1].linear)

    def test_fill_takes_the_nearest_called_step(self):
        self.assertEqual(am._fill([".", "A", ".", ".", "B"], circular=False), (["A", "A", "A", "B", "B"], 3))
        # cyclically the first step's nearest called neighbour is the last one
        self.assertEqual(am._fill([".", "A", "A", "A", "B"], circular=True), (["B", "A", "A", "A", "B"], 1))
        self.assertEqual(am._fill([".", "."], circular=True), ([".", "."], 0))


if __name__ == "__main__":
    unittest.main()

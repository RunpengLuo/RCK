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


class TestTotalCopyNumber(unittest.TestCase):
    """Total copy number (background + ecDNA) of the two RCK routes: total_cn_* and total_cn_ecdna_*."""

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        root = self.tmp.name
        os.makedirs(os.path.join(root, "sim", "truth"))
        # 0-based half-open: [100, 200) total 3, [200, 300) total 10 (5 + 5), [300, 400) total 2
        pd.DataFrame({"SID": [1, 2, 3], "#CHR": ["chr1"] * 3, "START": [100, 200, 300], "END": [200, 300, 400],
                      "CN_A": [1.5, 5.0, 1.0], "CN_B": [1.5, 5.0, 1.0]}).to_csv(
            os.path.join(root, "sim", "truth", "segments.tsv"), sep="\t", index=False)
        os.makedirs(os.path.join(root, "bulk"))
        pd.DataFrame({"sid": ["s1", "s2"], "chrom": ["chr1", "chr1"], "start": [100, 200], "end": [199, 399],
                      "inferred_total": [3, 8]}).to_csv(os.path.join(root, "bulk", "flow_segments.tsv"), sep="\t", index=False)
        os.makedirs(os.path.join(root, "mc"))
        with open(os.path.join(root, "mc", "rck.scnt.tsv"), "w") as handle:
            handle.write("chr\tstart\tend\textra\n")
            handle.write("chr1\t0\t99\tkind=flank;cn={'c1': {'A': 9, 'B': 9}, 'c2': {'A': 9, 'B': 9}}\n")
            handle.write("chr1\t100\t199\tkind=subgraph;cn={'c1': {'A': 2, 'B': 1}, 'c2': {'A': 1, 'B': 1}}\n")
            handle.write("chr1\t200\t299\tkind=subgraph;cn={'c1': {'A': 1, 'B': 12}, 'c2': {'A': 1, 'B': 6}}\n")
        self.root = root
        self.truth = [Cycle("e1", 5.0, [[CStep("chr1", 200, 299, "A", "+")]])]

    def tearDown(self):
        self.tmp.cleanup()

    def test_truth_half_open_to_inclusive(self):
        self.assertEqual(bm.truth_total_cn(os.path.join(self.root, "sim"))[0], ("chr1", 100, 199, 3.0))

    def test_bulk_scores(self):
        s = bm.total_cn_scores(bm.truth_total_cn(os.path.join(self.root, "sim")), bm.bulk_total_cn(os.path.join(self.root, "bulk")))
        # bases 100-199: c 3 t 3; 200-299: c 8 t 10; 300-399: c 8 t 2
        self.assertAlmostEqual(s["total_cn_recall"], (300 + 800 + 200) / (300 + 1000 + 200))
        self.assertAlmostEqual(s["total_cn_precision"], (300 + 800 + 200) / (300 + 800 + 800))
        self.assertAlmostEqual(s["total_cn_within_tol"], 100 / 300)  # |8-10| = 2 > max(1, 1); |8-2| > 1
        self.assertAlmostEqual(s["total_cn_mae"], (0 + 200 + 600) / 300)
        self.assertEqual((s["total_cn_bases"], s["total_cn_bases_no_truth"]), (300, 0))

    def test_multiclone_weights_flank_and_normal(self):
        rows = bm.multiclone_total_cn(os.path.join(self.root, "mc"), {"c1": 0.5, "c2": 0.3})
        self.assertEqual(len(rows), 2)  # the flank row is left out
        self.assertAlmostEqual(rows[0][3], 0.5 * 3 + 0.3 * 2 + 0.2 * 2)
        self.assertAlmostEqual(rows[1][3], 0.5 * 13 + 0.3 * 7 + 0.2 * 2)
        with self.assertRaises(KeyError):
            bm.multiclone_total_cn(os.path.join(self.root, "mc"), {"c1": 0.5})

    def test_ecdna_restricted(self):
        out = bm.total_cn_all(bm.truth_total_cn(os.path.join(self.root, "sim")), bm.bulk_total_cn(os.path.join(self.root, "bulk")),
                              self.truth)
        self.assertEqual(out["total_cn_ecdna_bases"], 100)
        self.assertAlmostEqual(out["total_cn_ecdna_recall"], 0.8)
        self.assertAlmostEqual(out["total_cn_ecdna_precision"], 1.0)

    def test_ecdna_bases_merged(self):
        cyc = [Cycle("a", 1.0, [[CStep("chr1", 10, 19, "A", "+"), CStep("chr1", 20, 29, "A", "-")]]),
               Cycle("b", 1.0, [[CStep("chr1", 25, 40, "A", "+"), CStep("chr2", 5, 6, "A", "+")]])]
        self.assertEqual(sorted(bm.ecdna_bases(cyc)), [("chr1", 10, 40, 1.0), ("chr2", 5, 6, 1.0)])

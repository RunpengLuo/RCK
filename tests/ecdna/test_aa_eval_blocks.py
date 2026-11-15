"""rck.ecdna.aa_eval.score(phase_blocks=...): phased scores within phase blocks (docs/pipeline/evaluation.md
§3.10), as bench_IW_intpsi scores the multi-clone route (Yubi, 2026-11-14). No solver."""
import os
import tempfile
import unittest

from rck.ecdna import aa_eval
from rck.ecdna.evaluate import CStep, Cycle

# two segments in two phase blocks (walk frame = 0-based inclusive; AA file = +1)
SEG1, SEG2 = (0, 999), (2000, 2999)
BLOCKS = [("chr1", 0, 1500, "b1"), ("chr1", 1501, 4000, "b2")]


def truth_cycle(name, weight, h1, h2):
    return Cycle(name, weight, [[CStep("chr1", SEG1[0], SEG1[1], h1, "+"), CStep("chr1", SEG2[0], SEG2[1], h2, "+")]])


def answer_text():
    """Inferred: the truth's two species with block b2's A/B swapped (an unoriented block)."""
    return "\n".join([
        "Interval\t1\tchr1\t1\t3000",
        "List of cycle segments",
        "Segment\t1\tchr1\t{a}\t{b}".format(a=SEG1[0] + 1, b=SEG1[1] + 1),
        "Segment\t2\tchr1\t{a}\t{b}".format(a=SEG2[0] + 1, b=SEG2[1] + 1),
        "Cycle=1;Copy_count=5.0;Segments=1+,2+;Haplotypes=A,B",
        "Cycle=2;Copy_count=3.0;Segments=1+,2+;Haplotypes=B,A",
    ]) + "\n"


class TestWithinBlock(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.path = os.path.join(self.tmp.name, "x_cycles.txt")
        with open(self.path, "w") as handle:
            handle.write(answer_text())
        self.truth = [truth_cycle("X", 5.0, "A", "A"), truth_cycle("Y", 3.0, "B", "B")]

    def tearDown(self):
        self.tmp.cleanup()

    def sample(self, result):
        return result.sample[result.sample["mode"] == "strict"].iloc[0]

    def test_absolute_by_default(self):
        s = self.sample(aa_eval.score(self.path, self.truth))
        self.assertEqual(s["phase_frame"], "absolute")
        self.assertLess(s["mean_composite"], 0.9)

    def test_within_block_pairs_and_phases(self):
        result = aa_eval.score(self.path, self.truth, phase_blocks=BLOCKS)
        s = self.sample(result)
        self.assertEqual(s["phase_frame"], "within_block")
        self.assertEqual(s["n_blocks_swapped"], 1)
        self.assertEqual((s["n_missed"], s["n_false"]), (0, 0))
        self.assertAlmostEqual(s["mean_composite"], s["mean_composite_h0"])
        m = result.matrix[(result.matrix["haplotype"] == "phased") & (~result.matrix["low_weight"].astype(bool))]
        best = m.groupby("truth")["j_seg"].max()
        self.assertTrue((best > 0.999).all())

    def test_absolute_counts_kept(self):
        s = self.sample(aa_eval.score(self.path, self.truth, phase_blocks=BLOCKS))
        absolute = self.sample(aa_eval.score(self.path, self.truth))
        self.assertEqual((s["n_missed_absolute"], s["n_false_absolute"]), (absolute["n_missed"], absolute["n_false"]))


if __name__ == "__main__":
    unittest.main()

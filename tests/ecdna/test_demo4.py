"""End-to-end runs on the demo4 subgraphs (two ecDNA species, three clones, two tumor samples).

amplicon_0 holds ecdna1 (haplotype B, 8 steps) plus a fragment of ecdna2; amplicon_1 holds
ecdna2 (haplotype A, 4 steps). These tests check that --truth-dir evaluation works on this
layout: it never errors, it scores each subgraph only against the species it contains, and
it uses the requested sample's copy numbers. Reconstruction accuracy is printed, not asserted:
demo4 is the first multi-species benchmark and its expected accuracy is not established yet.

The default tests use the vendored subgraphs and truth tables with depth per copy 13.5 (the
tumor's ~27x diploid median). Set RCK_ECDNA_DEMO4_OBS to
amplicon-allele-analysis/simulations/runs/demo4_clone3_ecdna/obs to also run with scales
estimated from the windows.
"""

import json
import os
import unittest

from tests.ecdna.helpers import DATA_DIR, requires_gurobi, temp_dir

TRUTH = os.path.join(DATA_DIR, "demo4_truth")
OBS = os.environ.get("RCK_ECDNA_DEMO4_OBS")
EXPECTED = {"amplicon_0": "ecdna1", "amplicon_1": "ecdna2"}
WEIGHT = {"tumor1": 6.2821, "tumor2": 7.4376}


def run(amplicon, outdir, extra):
    from rck.ecdna_run import build_parser, run as run_pipeline

    argv = ["--subgraph-dir", os.path.join(DATA_DIR, "demo4", amplicon), "--outdir", outdir, "--truth-dir", TRUTH,
            "--min-species-weight", "2", "--num-alternatives", "5"] + extra
    return run_pipeline(build_parser().parse_args(argv))


@requires_gurobi
class TestDemo4(unittest.TestCase):
    def check(self, amplicon, result, outdir, sample):
        evaluation = result["evaluation"]
        self.assertNotIn("error", evaluation, evaluation.get("traceback"))
        self.assertEqual(evaluation["n_truth_species"], 1)
        self.assertEqual(evaluation["n_truth_species_total"], 2)
        scored = [sp for sp in evaluation["truth_species"] if sp["in_scope"]]
        self.assertEqual([sp["species"] for sp in scored], [EXPECTED[amplicon]])
        self.assertAlmostEqual(scored[0]["weight"], WEIGHT[sample], places=3)
        for name in ("evaluation.json", "evaluation.tsv", "truth_cycles.tsv", "decompositions.tsv"):
            self.assertTrue(os.path.exists(os.path.join(outdir, name)), name)
        top = evaluation["top_decomposition"][0] if evaluation["top_decomposition"] else {}
        print("\ndemo4 {a} ({s}): order match at rank {r}; jaccard_hap {j}; junction recall {jr}; weight {w} vs {tw:.2f}".format(
            a=amplicon, s=sample, r=evaluation["first_rank_with_all_orders_matched"], j=top.get("jaccard_hap"),
            jr=top.get("junction_recall"), w=top.get("inferred_weight"), tw=WEIGHT[sample]))
        return evaluation

    def test_fixed_depth_scale(self):
        for amplicon in ("amplicon_0", "amplicon_1"):
            with self.subTest(amplicon=amplicon), temp_dir() as tmp:
                result = run(amplicon, tmp, ["--depth-per-copy", "13.5", "--truth-sample", "tumor1"])
                self.check(amplicon, result, tmp, "tumor1")

    @unittest.skipUnless(OBS and os.path.isdir(OBS or ""), "set RCK_ECDNA_DEMO4_OBS to the demo4 obs directory")
    def test_estimated_scales_both_samples(self):
        for sample in ("tumor1", "tumor2"):
            for amplicon in ("amplicon_0", "amplicon_1"):
                with self.subTest(sample=sample, amplicon=amplicon), temp_dir() as tmp:
                    result = run(amplicon, tmp, ["--obs-dir", OBS, "--sample", "sim_" + sample])
                    self.check(amplicon, result, tmp, sample)
                    print("  scales:", json.dumps({k: round(v, 3) for k, v in result["observations"].scale.to_dict().items()
                                                  if k in ("depth_per_copy", "junction_reads_per_copy")}))


if __name__ == "__main__":
    unittest.main()

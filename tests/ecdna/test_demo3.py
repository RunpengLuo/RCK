"""End-to-end run on the demo3 amplicon_0 subgraph, scored against the simulator truth.

The default test uses only files shipped in tests/ecdna/data (the subgraph and the small
truth tables) with fixed scales: depth per copy 15 (the matched normal's 30x diploid depth)
and junction reads per copy 6. Segment copy numbers then come from the subgraph's own
depth weights.

Set RCK_ECDNA_DEMO3_OBS to the full observation directory
(amplicon-allele-analysis/simulations/runs/demo3_clone1_ecdna/obs) to also run with scales
and segment depths estimated from the windows, extremities and long-read paths.
"""

import json
import os
import unittest

import pandas as pd

from tests.ecdna.helpers import DATA_DIR, requires_gurobi, temp_dir

SUBGRAPH = os.path.join(DATA_DIR, "amplicon_0")
TRUTH = os.path.join(DATA_DIR, "demo3_truth")
OBS = os.environ.get("RCK_ECDNA_DEMO3_OBS")


def run(outdir, extra):
    from rck.ecdna_run import build_parser, run as run_pipeline

    argv = ["--subgraph-dir", SUBGRAPH, "--outdir", outdir, "--truth-dir", TRUTH, "--num-alternatives", "5"] + extra
    return run_pipeline(build_parser().parse_args(argv))


@requires_gurobi
class TestDemo3(unittest.TestCase):
    def check(self, result, outdir):
        flow = result["flow"]
        self.assertGreater(sum(flow.ecdna_edges().values()), 0)
        evaluation = result["evaluation"]
        frame = pd.read_csv(os.path.join(outdir, "evaluation.tsv"), sep="\t")
        # judge the top-ranked decomposition only: a correct answer at a lower rank is not enough
        top = frame[frame["decomposition_rank"] == frame["decomposition_rank"].min()].iloc[0]
        self.assertEqual(evaluation["first_rank_with_all_orders_matched"], 1, frame.to_string())
        self.assertGreaterEqual(top.jaccard_hap, 0.95, frame.to_string())
        self.assertGreaterEqual(top.junction_recall, 0.9, frame.to_string())
        self.assertTrue(bool(top.order_match), frame.to_string())
        self.assertLessEqual(top.cn_rel_err, 0.3, frame.to_string())
        with open(os.path.join(outdir, "evaluation.json")) as handle:
            self.assertEqual(json.load(handle)["n_truth_species"], 1)
        return evaluation

    def test_fixed_scales_from_subgraph(self):
        with temp_dir() as tmp:
            result = run(tmp, ["--depth-per-copy", "15", "--junction-reads-per-copy", "6"])
            evaluation = self.check(result, tmp)
        print("\ndemo3 (fixed scales):", json.dumps(evaluation, indent=1, default=str))

    @unittest.skipUnless(OBS and os.path.isdir(OBS or ""), "set RCK_ECDNA_DEMO3_OBS to the demo3 obs directory")
    def test_scales_from_observations(self):
        with temp_dir() as tmp:
            result = run(tmp, ["--obs-dir", OBS])
            scale = result["observations"].scale
            self.assertAlmostEqual(scale.depth_per_copy, 15.0, delta=2.0)
            evaluation = self.check(result, tmp)
        print("\ndemo3 (estimated scales):", json.dumps(evaluation, indent=1, default=str))


if __name__ == "__main__":
    unittest.main()

"""--block-flip unoriented on the demo7 union subgraph (two nested ecDNA species on chr4).

subgraph_construct.py leaves phase blocks chr4:185 and chr4:186 unoriented, and against the
truth their allele observations have A and B swapped. Without block flips the flow fits badly
(objective 39.3, two imputed edges); with total-only observations it reaches 8.8 with none.

Needs a full Gurobi licence and a demo7 union subgraph rebuilt with segment_blocks.tsv (the
committed one predates it). Point RCK_ECDNA_DEMO7_SUBGRAPH at that folder, e.g. the output of

    python simulations/subgraph_construct.py --sim-dir simulations/runs/demo7_nestedcycles_ecdna \\
        --out-dir <new folder> --union --observed-seeds-only --strong-sv-support 20 \\
        --jump-sv-support 20 --rck-junction-reads-per-copy 6.5 \\
        --truth-segments simulations/runs/demo7_nestedcycles_ecdna/subgraph_union/truth_segments.bed

(run from amplicon-allele-analysis), then ``<new folder>/union``. The flow objective is printed,
not asserted.
"""

import os
import unittest

import pandas as pd

from tests.ecdna.helpers import have_full_gurobi, temp_dir

SUBGRAPH = os.environ.get("RCK_ECDNA_DEMO7_SUBGRAPH")


@unittest.skipUnless(SUBGRAPH and os.path.exists(os.path.join(SUBGRAPH or "", "segment_blocks.tsv")),
                     "set RCK_ECDNA_DEMO7_SUBGRAPH to a demo7 union subgraph with segment_blocks.tsv")
@unittest.skipUnless(have_full_gurobi(), "needs a full Gurobi licence")
class TestDemo7BlockFlip(unittest.TestCase):
    def test_unoriented_blocks_flipped(self):
        from rck.ecdna_run import build_parser, run

        with temp_dir() as tmp:
            result = run(build_parser().parse_args(["--subgraph-dir", SUBGRAPH, "--outdir", tmp, "--block-flip", "unoriented",
                                                    "--no-plots", "--num-alternatives", "2"]))
            blocks = pd.read_csv(os.path.join(tmp, "block_flips.tsv"), sep="\t")
        status = dict(zip(blocks["phase_block"], blocks["status"]))
        flow = result["flow"]
        print("\ndemo7 --block-flip unoriented: flow objective {o:.3f} (total-only run 8.8, alleles without flips 39.3); "
              "imputed used: {i}\n{b}".format(o=flow.objective, i=",".join(flow.imputed_used) or "none",
                                              b=blocks.to_string(index=False)))
        self.assertEqual(status["chr4:185"], "flipped")
        self.assertEqual(status["chr4:186"], "flipped")
        for block in ("chr4:171", "chr4:172", "chr4:187"):
            self.assertEqual(status[block], "fixed", block)
        self.assertEqual(flow.imputed_used, [])


if __name__ == "__main__":
    unittest.main()

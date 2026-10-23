import os
import unittest
from unittest import mock

import pandas as pd

from rck.ecdna import plot
from rck.ecdna.evaluate import CStep, Cycle, evaluate, write_evaluation, write_truth_table
from tests.ecdna.helpers import temp_dir

try:  # pragma: no cover - depends on the environment
    import cnplot  # noqa: F401
    import matplotlib  # noqa: F401

    HAVE_PLOTTING = True
except ImportError:  # pragma: no cover
    HAVE_PLOTTING = False

requires_plotting = unittest.skipUnless(HAVE_PLOTTING, "matplotlib or cnplot is not installed")

SEGMENTS = [("s1", "chr1", 100_000, 199_999), ("s2", "chr1", 200_000, 299_999), ("s3", "chr1", 300_000, 349_999),
            ("s4", "chr2", 1_000_000, 1_099_999)]
CIRCLE = [CStep("chr1", 100_000, 199_999, "A", "+"), CStep("chr1", 200_000, 299_999, "A", "+"),
          CStep("chr2", 1_000_000, 1_099_999, "B", "-")]


def walk(steps):
    return ",".join("{c}:{s}-{e}:{h}:{st}".format(c=x.chrom, s=x.start, e=x.end, h=x.hap, st=x.strand) for x in steps)


def write_run(outdir, linear=False):
    """A minimal rck-ecdna output folder: one circular species (or a linear one) over SEGMENTS."""
    joins = "reference:r1:AA,novel:n1:AB,novel:n2:BA" if not linear else "reference:r1:AA,novel:n1:AB,end:x:B"
    rows = [{"decomposition_rank": 1, "objective": 1.0, "amplicon_call": "linear" if linear else "ecDNA", "species": "ec1",
             "topology": "linear" if linear else "circular", "status": "species", "weight": 9.6, "weight_int": 10,
             "length_bp": 300_000, "ordering_rank": 1, "walk": walk(CIRCLE), "joins": joins}]
    pd.DataFrame(rows).to_csv(os.path.join(outdir, "decompositions.tsv"), sep="\t", index=False)
    seg_rows = []
    for sid, chrom, start, end in SEGMENTS:
        on_a = 10 if sid in ("s1", "s2") else 0
        on_b = 10 if sid == "s4" else 0
        seg_rows.append({"sid": sid, "chrom": chrom, "start": start, "end": end, "length": end - start + 1,
                         "obs_total": on_a + on_b + 2.1, "obs_A": on_a + 1.05, "obs_B": None if sid == "s3" else on_b + 1.0,
                         "obs_source": "test", "ecdna_A": on_a, "ecdna_B": on_b, "bg_A": 1, "bg_B": 1, "inferred_total": on_a + on_b + 2})
    pd.DataFrame(seg_rows).to_csv(os.path.join(outdir, "flow_segments.tsv"), sep="\t", index=False)


class TestAxisRegions(unittest.TestCase):
    def test_padded_and_merged(self):
        regions = plot.axis_regions(CIRCLE)
        self.assertEqual([r[0] for r in regions], ["chr1", "chr2"])
        chrom, start, end = regions[0]
        self.assertLess(start, 100_000)
        self.assertGreater(end, 300_000)  # adjacent steps merge into one padded interval

    def test_missing_cnplot_skips_plots(self):
        with temp_dir() as tmp:
            write_run(tmp)
            with mock.patch("importlib.import_module", side_effect=ImportError("no cnplot")):
                self.assertEqual(plot.plot_run(tmp, cnplot_path=tmp), {})


@requires_plotting
class TestPlotRun(unittest.TestCase):
    def test_species_pages_and_recovery_matrix(self):
        with temp_dir() as tmp:
            write_run(tmp)
            truth = [Cycle("truth1", 10.0, [CIRCLE])]
            write_truth_table(os.path.join(tmp, "truth_cycles.tsv"), truth)
            frame, summary = evaluate(os.path.join(tmp, "decompositions.tsv"), truth, tol=0)
            write_evaluation(tmp, frame, summary)
            written = plot.plot_run(tmp)
            self.assertTrue(os.path.getsize(written["pdf"]) > 0)
            names = sorted(os.path.basename(p) for p in written["png"])
            self.assertEqual(names, ["{n}_rank1_ec1.png".format(n=os.path.basename(tmp)),
                                     "{n}_recovery_matrix.png".format(n=os.path.basename(tmp))])

    def test_linear_species_without_truth(self):
        with temp_dir() as tmp:
            write_run(tmp, linear=True)
            written = plot.plot_run(tmp, png=False)
            self.assertTrue(os.path.exists(written["pdf"]))
            self.assertEqual(written["png"], [])


if __name__ == "__main__":
    unittest.main()

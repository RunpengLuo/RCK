import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna.scale import ScaleParams, Windows, _log_mode, build_observations
from rck.ecdna.subgraph_io import read_subgraph
from tests.ecdna.helpers import single_circle, write_windows


def background_windows(tmp, amp_depth=180.0, shared_bias=1.0):
    """chr9: 40 balanced 100 kb blocks at depth 30 (CN 2), 10 LOH blocks at depth 15 and 10
    blocks at depth 45 (CN 3, BAF 1/3). chr1: the single_circle amplicon at 1-301 kb, flanked
    by depth-30 background."""
    chroms = {
        "chr9": [(4_000_000 // 1000, 30.0, 30.0, 15, 15), (1_000_000 // 1000, 15.0, 30.0, 15, 0), (1_000_000 // 1000, 45.0, 30.0, 20, 10)],
        "chr1": [(1, 30.0 * shared_bias, 30.0 * shared_bias, 15, 15), (300, amp_depth * shared_bias, 30.0 * shared_bias, 28, 2),
                 (1_500, 30.0, 30.0, 15, 15)],
    }
    return write_windows(os.path.join(tmp, "windows"), chroms)


class TestDepthScale(unittest.TestCase):
    def test_log_mode(self):
        values = [30.0] * 40 + [15.0] * 10 + [45.0] * 10
        self.assertAlmostEqual(_log_mode(values), 30.0, delta=1.0)

    def test_log_mode_identical_values(self):
        # regression: the smoothing kernel is longer than the histogram here; this used to
        # raise IndexError or return a value shifted by one kernel offset (31.58 for 30)
        self.assertAlmostEqual(_log_mode([30.0] * 50), 30.0, places=6)
        self.assertAlmostEqual(_log_mode([1.0] * 20), 1.0, places=6)

    def test_log_mode_narrow_spread(self):
        values = [30.0 * (1 + 0.002 * i) for i in range(-5, 6)] * 4 + [45.0] * 5
        self.assertAlmostEqual(_log_mode(values), 30.0, delta=0.6)

    def test_gamma_from_balanced_blocks(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "sg"))
            g = read_subgraph(directory)
            windows = background_windows(tmp)
            obs = build_observations(g, ScaleParams(block_bp=100_000, exclude_pad_bp=100_000, junction_reads_per_copy=6.0), windows_dir=windows)
        self.assertAlmostEqual(obs.scale.depth_per_copy, 15.0, delta=0.75)
        self.assertEqual(obs.scale.depth_method, "tumor_mode")
        for so in obs.segments.values():
            self.assertEqual(so.source, "windows")
            self.assertAlmostEqual(so.total, 12.0, delta=0.7)
            self.assertAlmostEqual(so.a, so.total * 11 / 12.0, delta=1e-6)
            self.assertAlmostEqual(so.b, so.total / 12.0, delta=1e-6)
        # flanking background is depth 30 = 2 copies
        for value in obs.background.values():
            self.assertAlmostEqual(value, 2.0, delta=0.2)

    def test_override_and_subgraph_weight_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "sg"))
            g = read_subgraph(directory)
            obs = build_observations(g, ScaleParams(depth_per_copy=15.0, junction_reads_per_copy=6.0))
        self.assertEqual(obs.scale.depth_method, "override")
        for so in obs.segments.values():
            self.assertEqual(so.source, "subgraph_weight")
            self.assertAlmostEqual(so.total, 12.0)
        self.assertAlmostEqual(obs.adjacencies["circ"], 10.0)
        self.assertEqual(obs.background[next(iter(obs.background))], 2.0)

    def test_ratio_mode_removes_shared_bias(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "sg"))
            g = read_subgraph(directory)
            windows = background_windows(tmp, shared_bias=0.8)
            tumor = build_observations(g, ScaleParams(block_bp=100_000, exclude_pad_bp=100_000, junction_reads_per_copy=6.0), windows_dir=windows)
            ratio = build_observations(g, ScaleParams(block_bp=100_000, exclude_pad_bp=100_000, junction_reads_per_copy=6.0, depth_mode="ratio"),
                                       windows_dir=windows)
        for sid, so in ratio.segments.items():
            self.assertAlmostEqual(so.total, 12.0, delta=0.7)
            self.assertLess(tumor.segments[sid].total, 10.5)  # biased low without the normal

    def test_junction_scale_from_extremities(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "sg"))
            g = read_subgraph(directory)
            windows = background_windows(tmp)
            ext = os.path.join(tmp, "extremity.tsv")
            rows = [("chr9:{p}+".format(p=p), "chr9", p, p + 1, p, 12, 0) for p in range(150_000, 3_900_000, 200_000)]
            pd.DataFrame(rows, columns=["EID", "#CHR", "START", "END", "POS0", "#REF", "#NOV"]).to_csv(ext, sep="\t", index=False)
            obs = build_observations(g, ScaleParams(block_bp=100_000, exclude_pad_bp=100_000), windows_dir=windows, extremities_path=ext)
        self.assertEqual(obs.scale.junction_method, "extremities")
        self.assertAlmostEqual(obs.scale.junction_reads_per_copy, 6.0, delta=0.4)

    def test_junction_scale_fallback(self):
        with tempfile.TemporaryDirectory() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "sg"))
            g = read_subgraph(directory)
            obs = build_observations(g, ScaleParams(depth_per_copy=15.0))
        self.assertEqual(obs.scale.junction_method, "subgraph_reference_median")
        self.assertGreater(obs.scale.junction_reads_per_copy, 0)

    def test_sample_selection(self):
        with tempfile.TemporaryDirectory() as tmp:
            windows = background_windows(tmp)
            w = Windows(windows)
            self.assertEqual((w.normal_col, w.tumor_col), (0, 1))
            self.assertEqual(Windows(windows, sample="sim_tumor1").tumor_col, 1)
            with self.assertRaises(ValueError):
                Windows(windows, sample="missing")


if __name__ == "__main__":
    unittest.main()

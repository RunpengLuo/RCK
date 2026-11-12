"""--depth-bias-window: one multiplicative depth bias per window of nearby segments in the flow MILP.

The window grouping, the defaults, the output columns and the refit targets need no Gurobi; the
solves are marked ``requires_gurobi``.
"""
import unittest
from types import SimpleNamespace

from rck.core.structures import Haplotype
from rck.ecdna.ilp_circular import CircularParams, FlowSolution, depth_bias_windows
from rck.ecdna.refit import observation_rows
from rck.ecdna.scale import ScaleParams, build_observations
from rck.ecdna.subgraph_io import HAPLOTYPES, read_subgraph
from tests.ecdna.helpers import SubgraphBuilder, requires_gurobi, single_circle, temp_dir

FIXED = ScaleParams(depth_per_copy=15.0, junction_reads_per_copy=6.0)
CPARAMS = SimpleNamespace(seg_length_cap=50_000, min_seg_weight=0.02, alpha_ref=0.25, alpha_novel=0.5)


def fake_graph(*segments):
    """segments: (sid, chrom, start, end)."""
    return SimpleNamespace(segments={sid: SimpleNamespace(chrom=c, start=s, end=e, length=e - s + 1) for sid, c, s, e in segments})


class TestWindows(unittest.TestCase):
    def test_off(self):
        g = fake_graph(("s1", "chr1", 1, 1000))
        self.assertEqual(depth_bias_windows(g, 0), ({}, [], []))
        self.assertEqual(CircularParams().depth_bias_window, 0)

    def test_grouping_and_neighbours(self):
        g = fake_graph(("a", "chr1", 1, 100_000), ("b", "chr1", 100_001, 180_000),  # one window (180 kb)
                       ("c", "chr1", 180_001, 260_000),  # would make it 260 kb: new window, neighbour
                       ("d", "chr1", 260_001, 900_000),  # longer than the window: its own
                       ("e", "chr1", 5_000_001, 5_050_000),  # 4.1 Mb gap: new window, not a neighbour
                       ("f", "chr2", 1, 50_000))  # other chromosome
        seg, windows, pairs = depth_bias_windows(g, 200_000)
        self.assertEqual(seg["a"], seg["b"])
        self.assertEqual(len({seg["a"], seg["c"], seg["d"], seg["e"], seg["f"]}), 5)
        self.assertEqual(windows, ["chr1:0", "chr1:1", "chr1:2", "chr1:3", "chr2:0"])
        self.assertEqual(pairs, [("chr1:0", "chr1:1"), ("chr1:1", "chr1:2")])

    def test_unsorted_input(self):
        g = fake_graph(("b", "chr1", 100_001, 150_000), ("a", "chr1", 1, 100_000))
        seg, windows, pairs = depth_bias_windows(g, 200_000)
        self.assertEqual(seg, {"a": "chr1:0", "b": "chr1:0"})
        self.assertEqual(pairs, [])


class TestOutputsAndRefit(unittest.TestCase):
    def setUp(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            self.g = read_subgraph(directory)
        self.obs = build_observations(self.g, FIXED)

    def solution(self, seg_bias=None):
        z = {(sid, h): 0 for sid in self.g.segments for h in HAPLOTYPES}
        return FlowSolution(status="optimal", objective=0.0, mip_gap=0.0, seg_ec=dict(z), seg_bg=dict(z), adj_ec={}, adj_bg={},
                            terminal_bg={}, depth_bias={"w": 0.1} if seg_bias else {}, seg_bias=seg_bias or {})

    def test_default_segment_table_unchanged(self):
        table = self.solution().segment_table(self.g, self.obs)
        self.assertNotIn("depth_bias", table.columns)
        self.assertEqual(list(table.columns)[-1], "inferred_total")

    def test_segment_table_reports_bias(self):
        bias = {sid: 0.1 for sid in self.g.segments}
        table = self.solution(bias).segment_table(self.g, self.obs)
        self.assertEqual(list(table["depth_bias"]), [0.1] * len(self.g.segments))
        for r in table.itertuples():
            self.assertAlmostEqual(r.obs_corrected, r.obs_total * 0.9, places=9)

    def test_refit_targets(self):
        flow = SimpleNamespace(seg_bg={(sid, h): 1 for sid in self.g.segments for h in HAPLOTYPES}, adj_bg={})
        plain = observation_rows(self.g, self.obs, flow, CPARAMS)
        flow.seg_bias = {sid: 0.2 for sid in self.g.segments}
        biased = observation_rows(self.g, self.obs, flow, CPARAMS)
        self.assertEqual(len(plain), len(biased))
        for p, b in zip(plain, biased):
            if p.keys[0][0] == "S":  # segment rows: (value x 0.8) - background, with background 1 per haplotype
                self.assertAlmostEqual(b.target, (p.target + 1) * 0.8 - 1, places=9)
            else:  # junction rows are not touched
                self.assertEqual(b.target, p.target)

    def test_refit_without_attribute(self):
        """Flows built elsewhere (no seg_bias) keep the plain targets."""
        flow = SimpleNamespace(seg_bg={}, adj_bg={})
        self.assertEqual([r.target for r in observation_rows(self.g, self.obs, flow, CPARAMS)],
                         [r.target for r in observation_rows(self.g, self.obs, SimpleNamespace(seg_bg={}, adj_bg={}, seg_bias={}), CPARAMS)])


class TestCommandLine(unittest.TestCase):
    def test_flags(self):
        from rck.ecdna_run import build_parser

        args = build_parser().parse_args(["--subgraph-dir", "x", "--outdir", "y"])
        self.assertEqual((args.depth_bias_window, args.depth_bias_cost, args.depth_bias_smooth, args.depth_bias_max), (0, 1.0, 1.0, 0.3))
        args = build_parser().parse_args(["--subgraph-dir", "x", "--outdir", "y", "--depth-bias-window", "200000", "--depth-bias-cost", "2"])
        self.assertEqual((args.depth_bias_window, args.depth_bias_cost), (200_000, 2.0))


def biased_circle(directory, middle_depth):
    """single_circle with the middle segment's depth changed (no junction there, so balance keeps
    the copy number flat across all three segments)."""
    b = SubgraphBuilder()
    segs = b.chain("chr1", [1000, 101000, 201000, 301000], [180.0, middle_depth, 180.0], bafs=[1 / 12.0] * 3,
                   ref_reads=[12, 72, 72, 12])
    b.variant(segs[-1][1], segs[0][0], 60.0, hap="A", aid="circ")
    return b.write(directory)


@requires_gurobi
class TestSolve(unittest.TestCase):
    def solve(self, directory, **params):
        from rck.ecdna.ilp_circular import solve_flow

        g = read_subgraph(directory)
        obs = build_observations(g, FIXED)
        return g, solve_flow(g, obs, CircularParams(**params))

    def test_unbiased_data_fit_no_bias(self):
        with temp_dir() as tmp:
            g, sol = self.solve(single_circle(tmp)[0], depth_bias_window=50_000)
        self.assertEqual({sol.seg_ec[sid, Haplotype.A] for sid in g.segments}, {10})
        self.assertTrue(all(abs(v) < 1e-6 for v in sol.depth_bias.values()))
        self.assertEqual(set(sol.seg_bias), set(g.segments))

    def test_off_is_unchanged(self):
        with temp_dir() as tmp:
            directory = single_circle(tmp)[0]
            g, off = self.solve(directory)
            _, zero = self.solve(directory, depth_bias_window=0)
        self.assertEqual(off.seg_ec, zero.seg_ec)
        self.assertEqual(off.depth_bias, {})
        self.assertEqual(off.seg_bias, {})

    def test_low_middle_segment_is_absorbed_by_its_bias(self):
        """Middle segment read 17 % low (150 for 180): its window's bias takes the dip."""
        with temp_dir() as tmp:
            directory = biased_circle(tmp, 150.0)
            g, sol = self.solve(directory, depth_bias_window=50_000, depth_bias_cost=0.5, depth_bias_smooth=0.5)
        self.assertEqual({sol.seg_ec[sid, Haplotype.A] for sid in g.segments}, {10})
        middle = sorted(g.segments, key=lambda s: g.segments[s].start)[1]
        self.assertLess(sol.seg_bias[middle], -0.1)  # obs x (1 - beta) = 12 copies from 10 observed: beta about -0.2


if __name__ == "__main__":
    unittest.main()

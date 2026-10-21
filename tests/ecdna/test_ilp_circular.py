import os
import unittest

from rck.core.structures import Haplotype, Phasing
from rck.ecdna.phased import is_balanced
from rck.ecdna.scale import ScaleParams, build_observations
from rck.ecdna.subgraph_io import read_subgraph
from tests.ecdna.helpers import head, requires_gurobi, single_circle, SubgraphBuilder, tail, temp_dir

FIXED = ScaleParams(depth_per_copy=15.0, junction_reads_per_copy=6.0)


def solve(directory, **params):
    from rck.ecdna.ilp_circular import CircularParams, solve_flow

    g = read_subgraph(directory)
    obs = build_observations(g, FIXED)
    return g, obs, solve_flow(g, obs, CircularParams(**params))


@requires_gurobi
class TestSingleCircle(unittest.TestCase):
    def test_recovers_ten_copies_on_A(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g, obs, sol = solve(directory)
        for sid in g.segments:
            self.assertEqual(sol.seg_ec[sid, Haplotype.A], 10, sid)
            self.assertEqual(sol.seg_ec[sid, Haplotype.B], 0, sid)
            self.assertEqual(sol.seg_bg[sid, Haplotype.A] + sol.seg_bg[sid, Haplotype.B], 2, sid)
        self.assertEqual(sol.adj_ec["circ", Phasing.AA], 10)
        self.assertEqual(is_balanced(g, sol.ecdna_edges()), [])
        self.assertEqual(sol.imputed_used, [])
        self.assertEqual(sum(sol.linear_ends.values()), 0)  # the observed junction closes the circle

    def test_background_ends_only_at_terminals(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g, obs, sol = solve(directory)
        terminals = set(g.terminals())
        for (nid, h), value in sol.terminal_bg.items():
            self.assertIn(nid, terminals)
        self.assertEqual(sum(sol.terminal_bg.values()), 4)  # one copy per haplotype enters and leaves

    def test_model_is_disposed_after_solve(self):
        """Regression: an open Gurobi log file blocked temp-dir cleanup on NFS."""
        from rck.ecdna.ilp_circular import CircularParams, OptModelCircular

        with temp_dir() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "amp"))
            g = read_subgraph(directory)
            obs = build_observations(g, FIXED)
            log_file = os.path.join(tmp, "flow.log")
            model = OptModelCircular(g, obs, CircularParams(log_file=log_file))
            model.build()
            model.solve()
            self.assertIsNone(model.model)
            os.remove(log_file)  # the handle is closed, so this leaves nothing behind
            self.assertEqual([f for f in os.listdir(tmp) if f.startswith(".nfs")], [])

    def test_enforce_is_keeps_solution(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g, obs, sol = solve(directory, enforce_is=True)
        self.assertEqual(sol.adj_ec["circ", Phasing.AA], 10)


@requires_gurobi
class TestHaplotypeCrossing(unittest.TestCase):
    """ecDNA: s1 (A) -> crossing at 100 kb -> s2 (B) -> crossing at 200 kb -> s3 (A) -> back to s1."""

    def build(self, tmp):
        b = SubgraphBuilder()
        segs = b.chain("chr1", [0, 100_000, 200_000, 300_000], [180.0] * 3, bafs=[1 / 12.0, 11 / 12.0, 1 / 12.0],
                       ref_reads=[12, 12, 12, 12])
        b.variant(head("chr1", 100_000), tail("chr1", 100_000), 60, hap="mixed", aid="x1")
        b.variant(head("chr1", 200_000), tail("chr1", 200_000), 60, hap="mixed", aid="x2")
        b.variant(segs[-1][1], segs[0][0], 60, hap="A", aid="circ")
        return b.write(tmp)

    def test_crossing_labelings(self):
        with temp_dir() as tmp:
            g, obs, sol = solve(self.build(tmp))
        sids = sorted(g.segments, key=lambda s: g.segments[s].start)
        self.assertEqual(sol.seg_ec[sids[0], Haplotype.A], 10)
        self.assertEqual(sol.seg_ec[sids[1], Haplotype.B], 10)
        self.assertEqual(sol.seg_ec[sids[2], Haplotype.A], 10)
        a1 = g.adjacencies["x1"]
        expected = Phasing.AB if a1.u == head("chr1", 100_000) else Phasing.BA
        self.assertEqual(sol.adj_ec["x1", expected], 10)
        self.assertEqual(is_balanced(g, sol.ecdna_edges()), [])


@requires_gurobi
class TestImputedGap(unittest.TestCase):
    """The closing junction was not observed; only an imputed candidate can close the cycle."""

    def build(self, tmp):
        b = SubgraphBuilder()
        segs = b.chain("chr1", [1000, 101000, 201000, 301000], [180.0] * 3, bafs=[1 / 12.0] * 3, ref_reads=[12, 72, 72, 12])
        b.imputed(segs[-1][1], segs[0][0], 180.0, aid="imp_close")
        return b.write(tmp)

    def test_cheap_imputation_is_used(self):
        with temp_dir() as tmp:
            g, obs, sol = solve(self.build(tmp), imputed_penalty=0.5)
        self.assertEqual(sol.imputed_used, ["imp_close"])
        self.assertEqual(sol.adj_ec["imp_close", Phasing.AA], 10)

    def test_expensive_imputation_is_not_used(self):
        """Without an affordable closing edge, the amplified copies are explained as linear."""
        with temp_dir() as tmp:
            g, obs, sol = solve(self.build(tmp), imputed_penalty=1000.0)
        self.assertEqual(sol.imputed_used, [])
        self.assertEqual(sol.adj_ec["imp_close", Phasing.AA], 0)
        self.assertEqual(sum(sol.linear_ends.values()), 20)  # 10 copies, two ends each
        for sid in g.segments:
            self.assertEqual(sol.seg_ec[sid, Haplotype.A], 10, sid)

    def test_strictly_circular_mode(self):
        with temp_dir() as tmp:
            g, obs, sol = solve(self.build(tmp), imputed_penalty=1000.0, linear_ends="none")
        self.assertEqual(sum(sol.ecdna_edges().values()), 0)
        self.assertEqual(sol.linear_ends, {})


@requires_gurobi
class TestLinearEnds(unittest.TestCase):
    """chr1: two amplified segments (12 copies) then one background segment (2 copies)."""

    def build(self, tmp):
        b = SubgraphBuilder()
        segs = b.chain("chr1", [0, 100_000, 200_000, 300_000], [180.0, 180.0, 30.0], bafs=[1 / 12.0, 1 / 12.0, 0.5],
                       ref_reads=[12, 72, 12, 12])
        return b.write(tmp), segs

    def test_ends_at_boundaries_and_flanks_only(self):
        from rck.ecdna.ilp_circular import OptModelCircular

        with temp_dir() as tmp:
            directory, segs = self.build(tmp)
            g = read_subgraph(directory)
            obs = build_observations(g, FIXED)
            ends = set(OptModelCircular(g, obs).linear_end_nodes())
        self.assertIn(segs[0][0], ends)  # left flank of the region
        self.assertIn(segs[1][1], ends)  # copy number drops from 12 to 2 here ...
        self.assertIn(segs[2][0], ends)  # ... on both sides of the boundary
        self.assertIn(segs[2][1], ends)  # right flank
        self.assertNotIn(segs[0][1], ends)  # 12 -> 12: not a boundary
        self.assertNotIn(segs[1][0], ends)

    def test_linear_amplicon_has_linear_flow(self):
        with temp_dir() as tmp:
            directory, segs = self.build(tmp)
            g, obs, sol = solve(directory)
        sids = sorted(g.segments, key=lambda x: g.segments[x].start)
        self.assertEqual(sol.seg_ec[sids[0], Haplotype.A], 10)
        self.assertEqual(sol.seg_ec[sids[1], Haplotype.A], 10)
        self.assertEqual(sol.seg_ec[sids[2], Haplotype.A], 0)
        self.assertEqual(sol.linear_ends[segs[0][0], Haplotype.A], 10)
        self.assertEqual(sol.linear_ends[segs[1][1], Haplotype.A], 10)
        self.assertEqual(is_balanced(g, sol.ecdna_edges()), [])


if __name__ == "__main__":
    unittest.main()

"""--unexplained-junction-penalty: the count, the ranking and pruning, the pool record and the CLI.

Fixture (as tests/ecdna/test_step_penalty.py): a hub segment ``x`` and two loops ``a`` and ``b``,
each a circle through ``x``; the circle (x b) is nested in the sense that it shares ``x`` with
(x a). Outer (x a) at 36 copies and nested (x b) at 6, as bench_nested's 6:1 cases: the low_weight
threshold is max(1, 0.1 x 36) = 3.6, so the nested circle split into two 3-copy leftovers is
"one significant species", and the ranking (species count) prefers it to the right answer unless
the leftovers' junctions x>b and b>x (6 copies of flow each) are charged as unexplained.
"""

import unittest

from rck.ecdna.phased import adj_key, seg_key
from tests.ecdna.helpers import requires_gurobi, temp_dir
from tests.ecdna.test_step_penalty import add, circle_flow, nested_graph


def circle_counts(g, ends, name):
    """Edge multiplicities (1 each) of the circle x -> loop ``name`` -> x, haplotype A."""
    return {k: 1 for k in circle_flow(g, ends, name, 1)}


class Fixture(unittest.TestCase):
    def setUp(self):
        self.tmp = temp_dir()
        self.g, self.ends = nested_graph(self.tmp.name)
        self.flow = add(circle_flow(self.g, self.ends, "a", 36), circle_flow(self.g, self.ends, "b", 6))

    def tearDown(self):
        self.tmp.cleanup()

    def species(self, name, weight):
        from rck.ecdna.decompose import Species

        return Species(weight=weight, counts=circle_counts(self.g, self.ends, name))


class TestCount(Fixture):
    """unexplained_junctions and novel_junction_flow (no Gurobi)."""

    def test_junction_flow(self):
        from rck.ecdna.decompose import novel_junction_flow

        jflow = novel_junction_flow(self.g, self.flow)
        self.assertEqual(jflow, {"x>a": 36.0, "a>x": 36.0, "x>b": 6.0, "b>x": 6.0})  # reference adjacencies left out

    def test_counts(self):
        from rck.ecdna.decompose import DecompositionParams, novel_junction_flow, unexplained_junctions

        jflow, p = novel_junction_flow(self.g, self.flow), DecompositionParams()
        right = [self.species("a", 36), self.species("b", 6)]
        leftovers = [self.species("a", 36), self.species("b", 3), self.species("b", 3)]
        self.assertEqual(unexplained_junctions(right, jflow, p), 0)
        self.assertEqual(unexplained_junctions(leftovers, jflow, p), 2)  # x>b and b>x, 6 >= 3.6, low_weight only
        self.assertFalse(any(s.low_weight for s in leftovers))  # decided on a copy: the inputs are not flagged
        # a flow below the threshold is a leftover's business, not a species'
        small = add(circle_flow(self.g, self.ends, "a", 36), circle_flow(self.g, self.ends, "b", 3))
        self.assertEqual(unexplained_junctions(leftovers[:2], novel_junction_flow(self.g, small), p), 0)
        # the threshold follows --min-species-weight when it is set
        self.assertEqual(unexplained_junctions(leftovers, jflow, DecompositionParams(min_species_weight=5)), 2)
        self.assertEqual(unexplained_junctions(right, jflow, DecompositionParams(min_species_weight=5)), 0)
        self.assertEqual(unexplained_junctions(right, jflow, DecompositionParams(min_species_weight=7)), 0)  # 6 < 7
        self.assertEqual(unexplained_junctions(leftovers, jflow, DecompositionParams(min_species_weight=2)), 0)  # 3 >= 2
        self.assertEqual(unexplained_junctions([], jflow, p), 0)

    def test_single_species_threshold(self):
        from rck.ecdna.decompose import DecompositionParams, novel_junction_flow, unexplained_junctions

        # one species of 36 leaves x>b and b>x (6 >= 3.6) unexplained
        jflow = novel_junction_flow(self.g, self.flow)
        self.assertEqual(unexplained_junctions([self.species("a", 36)], jflow, DecompositionParams()), 2)

    def test_structure_cost(self):
        from rck.ecdna.decompose import DecompositionParams, structure_cost

        right = [self.species("a", 36), self.species("b", 6)]
        p = DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)
        self.assertAlmostEqual(structure_cost(right, p), 2 + 0.2 * 8)
        self.assertAlmostEqual(structure_cost(right, p, unexplained=3), 2 + 0.2 * 8 + 3)
        self.assertAlmostEqual(structure_cost(right, DecompositionParams(step_penalty=0.2), unexplained=3), 2 + 0.2 * 8)


class TestRanking(unittest.TestCase):
    """ranking_key, pool_keeps and rank_pool on recorded entries (no Gurobi)."""

    @staticmethod
    def entry(index, objective, n_significant, steps, unexplained=None, n_species=None, residual=0.0):
        e = {"index": index, "objective": objective, "residual": residual, "n_significant": n_significant,
             "n_species": n_significant if n_species is None else n_species, "steps": steps}
        if unexplained is not None:
            e["unexplained"] = unexplained
        return e

    def test_defaults_unchanged(self):
        from rck.ecdna.decompose import DecompositionParams, ranking_key

        p = DecompositionParams(step_penalty=0.2)
        self.assertEqual(DecompositionParams().unexplained_junction_penalty, 0.0)
        self.assertEqual(ranking_key(0.0, 1, 3, 4, 1.0, p), ranking_key(0.0, 1, 3, 4, 1.0, p, unexplained=5))

    def test_ranking_key(self):
        from rck.ecdna.decompose import DecompositionParams, ranking_key

        p = DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)
        leftovers = ranking_key(0.0, 1, 3, 4, 1.0, p, unexplained=2)  # 1 + 0.8 + 2
        right = ranking_key(0.0, 2, 2, 8, 2.0, p, unexplained=0)  # 2 + 1.6
        self.assertLess(right, leftovers)
        self.assertEqual(leftovers[1], 3.8)

    def test_rank_pool_prunes_with_the_prior(self):
        """The right answer lies beyond pool_gap and saves no steps: λ alone prunes it, the full
        structure cost keeps it (bench_nested s804), and then ranks it first."""
        from rck.ecdna.decompose import DecompositionParams, rank_pool

        entries = [self.entry(0, 1.0, 1, 4, unexplained=2, n_species=3), self.entry(1, 2.0, 2, 8, unexplained=0)]
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=0.2)), [0])
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)), [1, 0])
        # μ = 0.5 saves 1.0, less than the objective (1.0) and steps (0.8) it loses: pruned
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=0.5)), [0])
        # λ = 0 with the penalty: the prior alone keeps it too
        self.assertEqual(rank_pool(entries, DecompositionParams(unexplained_junction_penalty=1.0)), [1, 0])

    def test_old_pools_count_zero(self):
        from rck.ecdna.decompose import DecompositionParams, rank_pool

        entries = [self.entry(0, 1.0, 1, 4, n_species=3), self.entry(1, 2.0, 2, 8)]
        params = DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)
        self.assertEqual(rank_pool(entries, params), rank_pool(entries, DecompositionParams(step_penalty=0.2)))

    def test_pool_keeps(self):
        from rck.ecdna.decompose import DecompositionParams, pool_keeps

        p = DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)
        self.assertTrue(pool_keeps(2.0, 8, 1.0, 4, p, 0, 2))  # 1.6 + 1.0 < 0.8 + 2
        self.assertFalse(pool_keeps(2.0, 8, 1.0, 4, p))  # without counts: steps only
        self.assertTrue(pool_keeps(1.05, 40, 1.0, 4, p))  # within pool_gap

    def test_prune_pool(self):
        """prune_pool computes the counts itself from the junction flow."""
        from rck.ecdna.decompose import DecompositionParams, novel_junction_flow, prune_pool

        fx = Fixture()
        fx.setUp()
        try:
            jflow = novel_junction_flow(fx.g, fx.flow)
            leftovers = [fx.species("a", 36), fx.species("b", 3), fx.species("b", 3)]
            right = [fx.species("a", 36), fx.species("b", 6)]
            pool = [(1.0, leftovers), (2.0, right)]
            self.assertEqual(len(prune_pool(pool, DecompositionParams(step_penalty=0.2), jflow)), 1)
            params = DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)
            self.assertEqual(len(prune_pool(pool, params, jflow)), 2)
            self.assertEqual(len(prune_pool(pool, params)), 1)  # no junction flow given: nothing is unexplained
        finally:
            fx.tearDown()


class TestAmbiguity(unittest.TestCase):
    """ambiguity_summary with the new term (no Gurobi)."""

    def test_decided_by(self):
        from rck.ecdna.decompose import DecompositionParams, ambiguity_summary

        alts = [{"objective": 1.0, "n_species": 2, "steps": 8, "unexplained": 2, "weights": [30, 6]}]
        on = DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0)
        report = ambiguity_summary(1.0, 8, 2, alts, on, warn=False, top_unexplained=0)
        self.assertEqual(report["decided_by"], "unexplained_junction_penalty")
        self.assertAlmostEqual(report["gap"], 2.0)
        self.assertEqual(report["unexplained"], 0)
        off = ambiguity_summary(1.0, 8, 2, alts, DecompositionParams(step_penalty=0.2), warn=False, top_unexplained=0)
        self.assertEqual(off["decided_by"], "reads or tie-break")
        self.assertNotIn("unexplained", off)
        # alternatives recorded before the field existed count 0
        old = [{k: v for k, v in alts[0].items() if k != "unexplained"}]
        self.assertEqual(ambiguity_summary(1.0, 8, 2, old, on, warn=False)["gap"], 0.0)


class TestCommandLine(unittest.TestCase):
    def test_flag_parses(self):
        from rck.ecdna_run import build_parser

        base = ["--subgraph-dir", "x", "--outdir", "y"]
        self.assertEqual(build_parser().parse_args(base).unexplained_junction_penalty, 0.0)
        self.assertEqual(build_parser().parse_args(base + ["--unexplained-junction-penalty", "1"]).unexplained_junction_penalty, 1.0)

    def test_rerank_flag_parses(self):
        import io
        import json
        import os
        from contextlib import redirect_stdout

        from rck.ecdna.rerank import main

        tmp = temp_dir()
        try:
            os.makedirs(os.path.join(tmp.name, "pool"))
            entries = [TestRanking.entry(0, 1.0, 1, 4, unexplained=2, n_species=3),
                       TestRanking.entry(1, 2.0, 2, 8, unexplained=0)]
            pool = {"params": {"species_penalty": 1.0, "step_penalty": 0.2, "pool_gap": 0.1, "ambiguity_window": 0.0},
                    "lambdas": [0.2], "rank1": {}, "entries": entries, "same_fit": {}}
            with open(os.path.join(tmp.name, "pool", "pool.json"), "w") as handle:
                json.dump(pool, handle)
            for mu, index in ((None, 0), ("1", 1)):
                out = io.StringIO()
                with redirect_stdout(out):
                    main([tmp.name, "--lambdas", "0.2"] + (["--unexplained-junction-penalty", mu] if mu else []))
                row = out.getvalue().splitlines()[1].split("\t")
                self.assertEqual(int(row[1]), index, mu)
        finally:
            tmp.cleanup()


@requires_gurobi
class TestSolve(Fixture):
    """With a solve: the default ranking takes "outer + leftovers"; the penalty restores the nested circle."""

    def test_six_to_one(self):
        from rck.ecdna.decompose import DecompositionParams, decompose

        default = decompose(self.g, self.flow, DecompositionParams(step_penalty=0.2))[0]
        on = decompose(self.g, self.flow, DecompositionParams(step_penalty=0.2, unexplained_junction_penalty=1.0))[0]
        significant = lambda d: sorted(s.weight for s in d.species if not s.low_weight)  # noqa: E731
        self.assertEqual(significant(on), [6, 36])
        self.assertEqual(on.unexplained, 0)
        # the default ranks "outer + leftover" first: one significant species (33, with a 3-copy
        # low_weight walk x a x b carrying the rest), and x>b, b>x left to the leftover
        self.assertEqual(len(significant(default)), 1)
        self.assertEqual(default.unexplained, 2)


if __name__ == "__main__":
    unittest.main()

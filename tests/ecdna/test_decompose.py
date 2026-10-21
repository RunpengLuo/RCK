import json
import os
import unittest

import pandas as pd

from rck.ecdna.subgraph_io import read_subgraph
from tests.ecdna.helpers import loop_counts, loops_graph, requires_gurobi, single_circle, SubgraphBuilder, temp_dir


def add(*multisets):
    total = {}
    for counts in multisets:
        for k, v in counts.items():
            total[k] = total.get(k, 0) + v
    return total


@requires_gurobi
class TestDecompose(unittest.TestCase):
    def run_decompose(self, g, flow, **params):
        from rck.ecdna.decompose import DecompositionParams, decompose

        return decompose(g, flow, DecompositionParams(**params))

    def test_single_species_with_repeated_segment(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        flow = loop_counts(g, ends, ["a", "b"], weight=16)  # x visited twice per copy
        best = self.run_decompose(g, flow)[0]
        self.assertEqual(len(best.species), 1)
        sp = best.species[0]
        self.assertEqual(sp.weight, 16)
        self.assertEqual(max(sp.counts.values()), 2)
        self.assertAlmostEqual(best.residual, 0.0)

    def test_two_species_sharing_a_segment(self):
        """(x a) at 10 copies plus (x b) at 7. With weights 10 and 7 no single species fits
        (it would need a multiplicity above 3), but several two-species answers do, e.g.
        {x a: 10, x b: 7} and {x a x b: 7, x a: 3}. Both explain the flow exactly with two species,
        so both must be reported; the tie-break on edge traversals ranks the simpler one
        (8 traversals against 12) first."""
        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        flow = add(loop_counts(g, ends, ["a"], weight=10), loop_counts(g, ends, ["b"], weight=7))
        decs = self.run_decompose(g, flow, num_alternatives=10)
        best = decs[0]
        self.assertEqual(len(best.species), 2)
        self.assertAlmostEqual(best.residual, 0.0)
        self.assertEqual(tuple(sorted(s.weight for s in best.species)), (7, 10))
        # copy-number-equivalent answers: same residual and same number of species
        equivalent = {tuple(sorted(s.weight for s in d.species)) for d in decs
                      if abs(d.residual - best.residual) < 1e-6 and len(d.species) == len(best.species)}
        self.assertIn((3, 7), equivalent)

    def test_disconnected_circles_are_separate_species(self):
        with temp_dir() as tmp:
            d1 = os.path.join(tmp, "one")
            b = SubgraphBuilder()
            s1 = b.chain("chr1", [0, 1000, 2000], [100.0, 100.0])
            b.variant(s1[1][1], s1[0][0], 50, aid="c1")
            s2 = b.chain("chr2", [0, 1000, 2000], [100.0, 100.0])
            b.variant(s2[1][1], s2[0][0], 50, aid="c2")
            g = read_subgraph(b.write(d1))
        from rck.ecdna.phased import adj_key, seg_key

        flow = {}
        for a in g.adjacencies.values():
            interior = a.kind == "reference" and g.nodes[a.u].segment is not None and g.nodes[a.v].segment is not None
            if a.aid in ("c1", "c2") or interior:
                flow[adj_key(a.aid, "AA")] = 7
        for sid in g.segments:
            flow[seg_key(sid, "A")] = 7
        best = self.run_decompose(g, flow)[0]
        self.assertEqual([s.weight for s in best.species], [7, 7])
        self.assertAlmostEqual(best.residual, 0.0)

    def test_alternatives_are_distinct(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        # 10 copies through a and 10 through b: one species (x a x b) at 10, or two species at 10 each
        flow = loop_counts(g, ends, ["a", "b"], weight=10)
        for mode in ("pool", "nogood"):
            decs = self.run_decompose(g, flow, num_alternatives=5, pool_gap=2.0, mode=mode)
            sigs = [d.signature for d in decs]
            self.assertEqual(len(sigs), len(set(sigs)), mode)
            self.assertEqual(len(decs[0].species), 1, mode)
            self.assertGreaterEqual(len(decs), 2, mode)
            self.assertEqual(len(decs[1].species), 2, mode)

    def test_repeated_circle_is_one_circle_at_double_weight(self):
        """Regression: going round one circle twice at weight w is that circle at weight 2w."""
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g = read_subgraph(directory)
        from rck.ecdna.phased import adj_key, seg_key

        flow = {seg_key(sid, "A"): 10 for sid in g.segments}
        for a in g.adjacencies.values():
            if a.aid == "circ" or (a.kind == "reference" and a.reads and a.reads > 50):
                flow[adj_key(a.aid, "AA")] = 10
        decs = self.run_decompose(g, flow, num_alternatives=10, pool_gap=1.0)
        self.assertEqual(len(decs[0].species), 1)
        self.assertEqual(decs[0].species[0].weight, 10)
        self.assertEqual(set(decs[0].species[0].counts.values()), {1})
        sigs = [d.signature for d in decs]
        self.assertEqual(len(sigs), len(set(sigs)))
        for d in decs:
            for sp in d.species:
                self.assertIn(1, sp.counts.values())  # no species is a multiple of a smaller walk

    def test_low_weight_species_are_flagged(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        flow = add(loop_counts(g, ends, ["a"], weight=20), loop_counts(g, ends, ["b"], weight=1))
        best = self.run_decompose(g, flow)[0]
        self.assertEqual([s.weight for s in best.species], [20, 1])
        self.assertEqual([s.status for s in best.species], ["species", "low_weight"])
        explicit = self.run_decompose(g, flow, min_species_weight=1)[0]
        self.assertEqual([s.low_weight for s in explicit.species], [False, False])

    def test_orderings_are_reported(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp)
            g = read_subgraph(tmp)
        flow = loop_counts(g, ends, ["a", "b", "c"], weight=8)
        best = self.run_decompose(g, flow)[0]
        self.assertEqual(len(best.species), 1)
        self.assertEqual(len(best.species[0].orderings), 2)


@requires_gurobi
class TestCommandLine(unittest.TestCase):
    def test_single_circle_end_to_end(self):
        from rck.ecdna_run import main

        with temp_dir() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "amp"))
            out = os.path.join(tmp, "out")
            self.assertEqual(main(["--subgraph-dir", directory, "--outdir", out, "--depth-per-copy", "15",
                                   "--junction-reads-per-copy", "6"]), 0)
            for name in ("flow_segments.tsv", "flow_adjacencies.tsv", "decompositions.tsv", "summary.json", "observations.json",
                         os.path.join("aa", "amp_graph.txt"), os.path.join("aa", "amp_decomp1_cycles.txt")):
                self.assertTrue(os.path.exists(os.path.join(out, name)), name)
            table = pd.read_csv(os.path.join(out, "decompositions.tsv"), sep="\t")
            top = table[(table.decomposition_rank == 1) & (table.ordering_rank == 1)]
            self.assertEqual(len(top), 1)
            self.assertEqual(int(top.weight.iloc[0]), 10)
            self.assertEqual(top.walk.iloc[0].count(":A:"), 3)
            with open(os.path.join(out, "summary.json")) as handle:
                summary = json.load(handle)
            self.assertEqual(summary["flow"]["imputed_used"], [])
            with open(os.path.join(out, "aa", "amp_decomp1_cycles.txt")) as handle:
                self.assertIn("Copy_count=10.0", handle.read())
            self.assertEqual(top.status.iloc[0], "species")
            # regression: gurobipy's logger must not flood the console at the default verbosity
            import logging

            self.assertGreaterEqual(logging.getLogger("gurobipy").level, logging.WARNING)
            with open(os.path.join(out, "gurobi_flow.log")) as handle:
                self.assertIn("Optimal solution found", handle.read())


@requires_gurobi
class TestEvaluationNeverFailsTheRun(unittest.TestCase):
    def test_bad_truth_dir_is_reported_not_raised(self):
        from rck.ecdna_run import main

        with temp_dir() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "amp"))
            out = os.path.join(tmp, "out")
            code = main(["--subgraph-dir", directory, "--outdir", out, "--depth-per-copy", "15", "--junction-reads-per-copy", "6",
                         "--truth-dir", os.path.join(tmp, "does_not_exist")])
            self.assertEqual(code, 0)
            self.assertTrue(os.path.exists(os.path.join(out, "decompositions.tsv")))
            with open(os.path.join(out, "evaluation.json")) as handle:
                self.assertIn("error", json.load(handle))


if __name__ == "__main__":
    unittest.main()

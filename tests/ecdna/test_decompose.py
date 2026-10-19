import json
import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna.subgraph_io import read_subgraph
from tests.ecdna.helpers import SubgraphBuilder, loop_counts, loops_graph, requires_gurobi, single_circle


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
        with tempfile.TemporaryDirectory() as tmp:
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
        {x a: 10, x b: 7} and {x a x b: 7, x a: 3}. All of them must be reported."""
        with tempfile.TemporaryDirectory() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        flow = add(loop_counts(g, ends, ["a"], weight=10), loop_counts(g, ends, ["b"], weight=7))
        decs = self.run_decompose(g, flow, num_alternatives=10)
        best = decs[0]
        self.assertEqual(len(best.species), 2)
        self.assertAlmostEqual(best.residual, 0.0)
        optimal = {tuple(sorted(s.weight for s in d.species)) for d in decs if abs(d.objective - best.objective) < 1e-6}
        self.assertIn((7, 10), optimal)
        self.assertIn((3, 7), optimal)

    def test_disconnected_circles_are_separate_species(self):
        with tempfile.TemporaryDirectory() as tmp:
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
        with tempfile.TemporaryDirectory() as tmp:
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

    def test_orderings_are_reported(self):
        with tempfile.TemporaryDirectory() as tmp:
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

        with tempfile.TemporaryDirectory() as tmp:
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


if __name__ == "__main__":
    unittest.main()

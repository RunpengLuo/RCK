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


def merged_read(ends, support=30.0, pid="r1"):
    """A read through a>x then x>b: only the merged circle (x a x b) contains both junctions."""
    from rck.ecdna.walk import ReadPath

    return ReadPath(pid=pid, support=support, junctions=[(ends["a"][1], ends["x"][0]), (ends["x"][1], ends["b"][0])])


class TestReadConstraints(unittest.TestCase):
    def test_mapping_merging_and_filters(self):
        from rck.ecdna.decompose import read_constraints
        from rck.ecdna.walk import ReadPath

        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b", "c"))
            g = read_subgraph(tmp)
        flow = add(loop_counts(g, ends, ["a"], weight=10), loop_counts(g, ends, ["b"], weight=7))
        forward = merged_read(ends, 30.0, "fwd")
        backward = ReadPath(pid="rev", support=5.0, junctions=[(b, a) for a, b in reversed(forward.junctions)])
        single = ReadPath(pid="one", support=50.0, junctions=[(ends["x"][1], ends["a"][0])])
        off_flow = ReadPath(pid="c", support=9.0, junctions=[(ends["c"][1], ends["x"][0]), (ends["x"][1], ends["a"][0])])
        weak = ReadPath(pid="weak", support=1.0, junctions=forward.junctions)
        reads, skipped = read_constraints(g, sorted(flow), [forward, backward, single, off_flow, weak], min_support=2.0)
        self.assertEqual(len(reads), 1)  # forward and its reverse are one junction sequence; single, weak dropped
        self.assertEqual(reads[0].support, 35.0)
        self.assertEqual(sorted(reads[0].pids), ["fwd", "rev"])
        self.assertEqual([[e[1] for e in grp] for grp in reads[0].groups], [["a>x"], ["x>b"]])
        self.assertEqual(skipped, 1)  # c>x carries no flow


@requires_gurobi
class TestReadsInDecomposition(unittest.TestCase):
    """(x a) 10 + (x b) 7 has two exact two-species answers; a read through a>x>b picks the merged one."""

    def setUp(self):
        self.tmp = temp_dir()
        directory = self.tmp.__enter__()
        self.ends = loops_graph(directory, loops=("a", "b"))
        self.g = read_subgraph(directory)
        self.flow = add(loop_counts(self.g, self.ends, ["a"], weight=10), loop_counts(self.g, self.ends, ["b"], weight=7))

    def tearDown(self):
        self.tmp.__exit__(None, None, None)

    def run_decompose(self, paths, **params):
        from rck.ecdna.decompose import DecompositionParams, decompose

        return decompose(self.g, self.flow, DecompositionParams(num_alternatives=10, **params), paths=paths)

    @staticmethod
    def weights(d):
        return tuple(sorted(s.weight for s in d.species))

    def test_without_reads_the_simpler_answer_wins(self):
        self.assertEqual(self.weights(self.run_decompose(None)[0]), (7, 10))

    def test_read_selects_the_merged_circle(self):
        best = self.run_decompose([merged_read(self.ends)])[0]
        self.assertEqual(self.weights(best), (3, 7))
        merged = [s for s in best.species if s.weight == 7][0]
        self.assertEqual(sum(merged.counts.values()), 8)  # x twice, a, b and four junctions
        self.assertAlmostEqual(best.residual, 0.0)
        self.assertEqual(best.reads_total, 30.0)
        self.assertEqual(best.reads_explained, 1.0)
        self.assertEqual(merged.orderings[0].paths_matched, 1)

    def test_required_fraction_without_reward(self):
        best = self.run_decompose([merged_read(self.ends)], path_weight=0.0, path_fraction=0.9)[0]
        self.assertEqual(self.weights(best), (3, 7))

    def test_reads_never_add_a_species(self):
        """The read reward is capped at path_weight (< species_penalty): a single circle stays single."""
        flow = loop_counts(self.g, self.ends, ["a", "b"], weight=10)
        from rck.ecdna.decompose import DecompositionParams, decompose
        from rck.ecdna.walk import ReadPath

        loop_a = ReadPath(pid="aa", support=40.0, junctions=[(self.ends["x"][1], self.ends["a"][0]), (self.ends["a"][1], self.ends["x"][0]),
                                                           (self.ends["x"][1], self.ends["a"][0])])
        best = decompose(self.g, flow, DecompositionParams(), paths=[loop_a])[0]
        self.assertEqual(len(best.species), 1)
        self.assertEqual(best.reads_explained, 0.0)  # x a x a is not a sub-walk of x a x b


@requires_gurobi
class TestLinearSpecies(unittest.TestCase):
    def run_decompose(self, g, flow, **params):
        from rck.ecdna.decompose import DecompositionParams, decompose

        return decompose(g, flow, DecompositionParams(**params))

    def test_linear_walk_is_one_linear_species(self):
        from tests.ecdna.test_walk import linear_chain_counts

        with temp_dir() as tmp:
            b = SubgraphBuilder()
            b.chain("chr1", [0, 1000, 2000, 3000], [100.0] * 3)
            g = read_subgraph(b.write(tmp))
        flow = {k: 8 * v for k, v in linear_chain_counts(g).items()}
        best = self.run_decompose(g, flow)[0]
        self.assertEqual(len(best.species), 1)
        self.assertEqual(best.species[0].topology, "linear")
        self.assertEqual(best.species[0].weight, 8)
        self.assertEqual(best.call, "linear")
        self.assertAlmostEqual(best.circular_fraction, 0.0)
        steps = best.species[0].orderings[0].steps
        self.assertEqual(len(steps), 3)
        self.assertEqual(steps[-1].join[0], "T")

    def test_circle_plus_linear_piece(self):
        """(x a) circle at 10 copies and a linear walk over b at 2: two species, called ecDNA."""
        from rck.ecdna.phased import end_key, seg_key

        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        flow = loop_counts(g, ends, ["a"], weight=10)
        b_sid = g.nodes[ends["b"][0]].segment
        flow[seg_key(b_sid, "A")] = 2
        flow[end_key(ends["b"][0], "A")] = 2
        flow[end_key(ends["b"][1], "A")] = 2
        best = self.run_decompose(g, flow)[0]
        self.assertEqual(sorted((s.topology, s.weight) for s in best.species), [("circular", 10), ("linear", 2)])
        self.assertEqual(best.call, "ecDNA")
        self.assertAlmostEqual(best.circular_fraction, 20 / 22.0, places=3)
        self.assertAlmostEqual(best.residual, 0.0)


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
            self.assertEqual(top.topology.iloc[0], "circular")
            self.assertEqual(top.amplicon_call.iloc[0], "ecDNA")
            self.assertEqual(summary["amplicon_call"], "ecDNA")
            # regression: gurobipy's logger must not flood the console at the default verbosity
            import logging

            self.assertGreaterEqual(logging.getLogger("gurobipy").level, logging.WARNING)
            with open(os.path.join(out, "gurobi_flow.log")) as handle:
                self.assertIn("Optimal solution found", handle.read())


@requires_gurobi
class TestLinearCommandLine(unittest.TestCase):
    def test_amplicon_without_closing_junction_is_called_linear(self):
        from rck.ecdna_run import main

        with temp_dir() as tmp:
            b = SubgraphBuilder()
            b.chain("chr1", [1000, 101000, 201000, 301000], [180.0] * 3, bafs=[1 / 12.0] * 3, ref_reads=[12, 72, 72, 12])
            directory = b.write(os.path.join(tmp, "lin"))
            out = os.path.join(tmp, "out")
            self.assertEqual(main(["--subgraph-dir", directory, "--outdir", out, "--depth-per-copy", "15",
                                   "--junction-reads-per-copy", "6"]), 0)
            table = pd.read_csv(os.path.join(out, "decompositions.tsv"), sep="\t")
            top = table[(table.decomposition_rank == 1) & (table.ordering_rank == 1)]
            self.assertEqual(list(top.topology), ["linear"])
            self.assertEqual(int(top.weight.iloc[0]), 10)
            self.assertEqual(top.amplicon_call.iloc[0], "linear")
            with open(os.path.join(out, "summary.json")) as handle:
                summary = json.load(handle)
            self.assertEqual(summary["amplicon_call"], "linear")
            self.assertEqual(summary["flow"]["linear_end_copies"], 20)
            with open(os.path.join(out, "aa", "lin_decomp1_cycles.txt")) as handle:
                self.assertIn("Segments=0+,", handle.read())
            # the strictly circular model finds no ecDNA here
            out2 = os.path.join(tmp, "out2")
            main(["--subgraph-dir", directory, "--outdir", out2, "--depth-per-copy", "15", "--junction-reads-per-copy", "6",
                  "--linear-ends", "none"])
            with open(os.path.join(out2, "summary.json")) as handle:
                self.assertEqual(json.load(handle)["amplicon_call"], "none")


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

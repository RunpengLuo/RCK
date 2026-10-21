import json
import os
import unittest
from types import SimpleNamespace

import pandas as pd

from rck.ecdna.refit import RefitParams, Row, misfit, observation_rows, refit_weights
from rck.ecdna.subgraph_io import HAPLOTYPES, read_subgraph
from tests.ecdna.helpers import requires_gurobi, single_circle, temp_dir

E1, E2, E3 = ("S", "s1", "A"), ("S", "s2", "A"), ("S", "s3", "A")
MAIN = {E1: 1, E2: 1, E3: 1}
LEFTOVER = {E1: 1}


def rows(*targets):
    return [Row(1.0, t, (k,)) for k, t in zip((E1, E2, E3), targets)]


class TestRefitWeights(unittest.TestCase):
    """The LP needs scipy only, no Gurobi."""

    def test_fractional_weight(self):
        result = refit_weights(rows(10.4, 10.4, 10.4), [MAIN], [10])
        self.assertAlmostEqual(result.weights[0], 10.4, places=6)
        self.assertLess(result.objective, result.objective_int)

    def test_unsupported_leftover_is_dropped(self):
        """Integer answer 9 + 1 on an uneven flow; the data say one molecule at about 10.2."""
        result = refit_weights(rows(10.2, 10.3, 10.1), [MAIN, LEFTOVER], [9, 1])
        self.assertEqual(result.kept, [0])
        self.assertEqual([i for i, _ in result.dropped], [1])
        self.assertAlmostEqual(result.weights[0], 10.2, places=6)  # the L1 fit is the median

    def test_supported_leftover_is_kept(self):
        result = refit_weights(rows(13.0, 10.3, 10.1), [MAIN, LEFTOVER], [10, 3])
        self.assertEqual(result.kept, [0, 1])
        self.assertAlmostEqual(sum(result.weights), 13.0, places=6)
        self.assertGreater(result.weights[1], 2.5)

    def test_min_copies_threshold(self):
        loose = refit_weights(rows(10.6, 10.2, 10.2), [MAIN, LEFTOVER], [10, 1], RefitParams(min_copies=0.1))
        self.assertEqual(loose.kept, [0, 1])
        strict = refit_weights(rows(10.6, 10.2, 10.2), [MAIN, LEFTOVER], [10, 1], RefitParams(min_copies=0.5))
        self.assertEqual(strict.kept, [0])

    def test_last_species_is_never_dropped(self):
        result = refit_weights(rows(0.2, 0.2, 0.2), [MAIN], [1])
        self.assertEqual(result.kept, [0])
        self.assertAlmostEqual(result.weights[0], 0.2, places=6)

    def test_observation_rows_use_the_flow_background(self):
        from rck.ecdna.scale import ScaleParams, build_observations

        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g = read_subgraph(directory)
        obs = build_observations(g, ScaleParams(depth_per_copy=15.0, junction_reads_per_copy=6.0))
        flow = SimpleNamespace(seg_bg={(sid, h): 1 for sid in g.segments for h in HAPLOTYPES}, adj_bg={})
        table = observation_rows(g, obs, flow, SimpleNamespace(seg_length_cap=50_000, min_seg_weight=0.02, alpha_ref=0.25,
                                                               alpha_novel=0.5))
        seg_a = [r for r in table if r.keys[0][0] == "S" and r.keys[0][2] == "A"]
        self.assertEqual(len(seg_a), 3)
        for r in seg_a:
            self.assertAlmostEqual(r.target, 10.0, places=6)  # 11 A copies minus 1 background
        circ = [r for r in table if r.keys[0][1] == "circ"]
        self.assertEqual(len(circ), 1)
        self.assertAlmostEqual(circ[0].target, 10.0, places=6)
        counts = {k: 1 for r in seg_a for k in r.keys}
        counts.update({k: 1 for k in circ[0].keys})
        self.assertGreater(misfit(table, [counts], [9]), misfit(table, [counts], [10]))


@requires_gurobi
class TestRefitCommandLine(unittest.TestCase):
    def run_cli(self, tmp, extra=()):
        from rck.ecdna_run import main

        # 12.4 copies at BAF 1/12 and 62.4 junction reads: about 10.4 ecDNA copies on A
        directory, _ = single_circle(os.path.join(tmp, "amp"), depth=186.0, reads=62.4)
        out = os.path.join(tmp, "out" + str(len(extra)))
        self.assertEqual(main(["--subgraph-dir", directory, "--outdir", out, "--depth-per-copy", "15",
                               "--junction-reads-per-copy", "6"] + list(extra)), 0)
        table = pd.read_csv(os.path.join(out, "decompositions.tsv"), sep="\t")
        with open(os.path.join(out, "summary.json")) as handle:
            summary = json.load(handle)
        return table[(table.decomposition_rank == 1) & (table.ordering_rank == 1)], summary

    def test_fractional_copy_number(self):
        with temp_dir() as tmp:
            top, summary = self.run_cli(tmp)
        self.assertEqual(len(top), 1)
        self.assertEqual(int(top.weight_int.iloc[0]), 10)
        self.assertAlmostEqual(float(top.weight.iloc[0]), 10.4, delta=0.1)
        d = summary["decompositions"][0]
        self.assertLessEqual(d["data_objective"], d["data_objective_int"])
        self.assertEqual(d["pruned"], [])

    def test_no_refit_keeps_integers(self):
        with temp_dir() as tmp:
            top, _ = self.run_cli(tmp, ["--no-refit"])
        self.assertEqual(float(top.weight.iloc[0]), 10.0)


if __name__ == "__main__":
    unittest.main()

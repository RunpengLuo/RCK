"""Running RCK ecDNA from the subgraph folder alone (scale.json, CN_* columns, cn_profile.tsv, paths.tsv)."""

import argparse
import json
import os
import unittest

import pandas as pd

from rck.ecdna.scale import ScaleParams, build_observations, load_subgraph_scale
from rck.ecdna.subgraph_io import read_subgraph
from tests.ecdna.helpers import requires_gurobi, single_circle, temp_dir


def add_rck_inputs(directory, gamma=15.0, theta=6.0, flank_cn=2.5, scale_gamma=None, sample="sim_tumor1"):
    """What subgraph_construct.py (rck_inputs.py) adds to a subgraph folder.

    ``scale_gamma`` lets scale.json disagree with the gamma used for the CN columns, to show
    which one RCK uses.
    """
    path = os.path.join(directory, "edges.tsv")
    edges = pd.read_csv(path, sep="\t")
    for c in ("CN_TOTAL", "CN_A", "CN_B", "CN"):
        edges[c] = float("nan")
    seg = edges["EDGE_TYPE"] == "reference"
    edges.loc[seg, "CN_TOTAL"] = edges.loc[seg, "WEIGHT"] / gamma
    edges.loc[seg, "CN_A"] = edges.loc[seg, "CN_TOTAL"] * (1 - edges.loc[seg, "BAF"])
    edges.loc[seg, "CN_B"] = edges.loc[seg, "CN_TOTAL"] * edges.loc[seg, "BAF"]
    junction = edges["EDGE_TYPE"].isin(["extremity", "variant"])
    edges.loc[junction, "CN"] = edges.loc[junction, "WEIGHT"] / theta
    edges.to_csv(path, sep="\t", index=False)
    pd.DataFrame([{"#CHR": "chr1", "START": s, "END": s + 100_000, "CN_TOTAL": flank_cn, "N_BINS": 100}
                  for s in range(0, 3_000_000, 100_000)]).to_csv(os.path.join(directory, "cn_profile.tsv"), sep="\t", index=False)
    with open(os.path.join(directory, "scale.json"), "w") as handle:
        json.dump({"format": "rck-ecdna-subgraph-scale/1", "sample": sample, "depth_per_copy": scale_gamma or gamma,
                   "junction_reads_per_copy": theta, "junction_sample": sample, "purity": 0.8, "tumor_ploidy": 2.0,
                   "cell_ploidy": flank_cn, "background_cn": flank_cn, "ploidy_class": 2, "method": "test",
                   "files": {"cn_profile": "cn_profile.tsv", "paths": None}}, handle)
    return directory


class TestSubgraphObservations(unittest.TestCase):
    def observations(self, **kwargs):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            add_rck_inputs(directory, **kwargs)
            g = read_subgraph(directory)
            scale = load_subgraph_scale(directory)
            obs = build_observations(g, ScaleParams(baseline_cn=2.5), subgraph_scale=scale)
        return g, obs

    def test_scales_and_copy_numbers_come_from_the_subgraph(self):
        g, obs = self.observations()
        self.assertEqual(obs.scale.depth_method, "subgraph_scale")
        self.assertEqual(obs.scale.junction_method, "subgraph_scale")
        self.assertAlmostEqual(obs.scale.depth_per_copy, 15.0)
        for so in obs.segments.values():
            self.assertEqual(so.source, "subgraph_cn")
            self.assertAlmostEqual(so.total, 12.0)
            self.assertAlmostEqual(so.a, 11.0)
            self.assertAlmostEqual(so.b, 1.0)
        self.assertAlmostEqual(obs.adjacencies["circ"], 10.0)

    def test_flank_background_from_profile(self):
        _, obs = self.observations(flank_cn=3.6)
        for value in obs.background.values():
            self.assertAlmostEqual(value, 3.6)

    def test_cn_columns_win_over_the_scale(self):
        """The CN columns are what subgraph_construct.py measured; scale.json's gamma is only reported."""
        _, obs = self.observations(scale_gamma=7.5)
        self.assertAlmostEqual(obs.scale.depth_per_copy, 7.5)
        for so in obs.segments.values():
            self.assertAlmostEqual(so.total, 12.0)

    def test_depth_override_recomputes_from_weight(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            add_rck_inputs(directory)
            g = read_subgraph(directory)
            obs = build_observations(g, ScaleParams(depth_per_copy=30.0), subgraph_scale=load_subgraph_scale(directory))
        self.assertEqual(obs.scale.depth_method, "override")
        for so in obs.segments.values():
            self.assertAlmostEqual(so.total, 6.0)  # 180 / 30
            self.assertEqual(so.source, "subgraph_weight")


class TestInputMode(unittest.TestCase):
    def resolve(self, directory, mode="auto", **extra):
        from rck.ecdna_run import _resolve_inputs

        args = argparse.Namespace(subgraph_dir=directory, inputs=mode, windows_dir=None, extremities=None, paths=None, obs_dir=None)
        vars(args).update(extra)
        return _resolve_inputs(args)

    def test_auto(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            self.assertEqual(self.resolve(directory)["mode"], "windows")
            add_rck_inputs(directory)
            pd.DataFrame([{"PID": "p1", "PATH": "x@+,y@-", "NSTEP": 2, "#NOV": 3}]).to_csv(os.path.join(directory, "paths.tsv"), sep="\t", index=False)
            inputs = self.resolve(directory)
            self.assertEqual(inputs["mode"], "subgraph")
            self.assertEqual(inputs["paths"], os.path.join(directory, "paths.tsv"))
            self.assertEqual(self.resolve(directory, mode="windows")["mode"], "windows")

    def test_subgraph_mode_needs_scale_json(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            with self.assertRaises(SystemExit):
                self.resolve(directory, mode="subgraph")


@requires_gurobi
class TestSubgraphOnlyRun(unittest.TestCase):
    def test_reconstruct_from_subgraph_folder_only(self):
        from rck.ecdna_run import main

        with temp_dir() as tmp:
            directory, _ = single_circle(os.path.join(tmp, "amp"))
            add_rck_inputs(directory, flank_cn=2.0)
            out = os.path.join(tmp, "out")
            self.assertEqual(main(["--subgraph-dir", directory, "--outdir", out]), 0)
            table = pd.read_csv(os.path.join(out, "decompositions.tsv"), sep="\t")
            top = table[(table.decomposition_rank == 1) & (table.ordering_rank == 1)]
            self.assertEqual(len(top), 1)
            self.assertAlmostEqual(float(top.weight.iloc[0]), 10.0, delta=0.2)
            with open(os.path.join(out, "summary.json")) as handle:
                summary = json.load(handle)
            self.assertEqual(summary["inputs"]["mode"], "subgraph")
            self.assertEqual(summary["amplicon_call"], "ecDNA")
            with open(os.path.join(out, "observations.json")) as handle:
                observed = json.load(handle)
            self.assertEqual(observed["scale"]["depth_method"], "subgraph_scale")


if __name__ == "__main__":
    unittest.main()

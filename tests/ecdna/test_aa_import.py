"""Graph route, offline part (no Gurobi): the AA graph writer, ``aa_export subgraph`` and
``aa_import`` (subgraph -> AA graph -> subgraph).

The round trip must give back the normalised graph ``rck-ecdna`` solves on, exactly: nodes,
segments, adjacencies (ids, ends, kinds, reads, orders), observations and read paths
(``rck.ecdna.aa_import.compare_subgraphs``). Path-constraint tokens are compared with the
ones CoRAL 3.1.1 writes itself (``tests/ecdna/data/graph_route/coral311_amplicon0_path_tokens.tsv``,
made by results/graph_route/check_coral_format.py --golden).

The demo5-demo8 and bench cases come from amplicon-allele-analysis ($AAA, default
~/amplicon-allele-analysis) and are skipped when it is not there.

Run from the RCK root: python -m unittest tests.ecdna.test_aa_import
"""
import json
import os
import shutil
import tempfile
import unittest

import pandas as pd

from rck.ecdna import aa_format as aaf
from rck.ecdna.aa_export import export_subgraph, main as export_main
from rck.ecdna.aa_import import compare_subgraphs, import_graph
from rck.ecdna.subgraph_io import IMPUTED, NOVEL, read_subgraph

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, "data")
AMP0 = os.path.join(DATA, "amplicon_0")
GOLDEN = os.path.join(DATA, "graph_route")
AMP0_SCALE = {"depth_per_copy": 15.0, "junction_reads_per_copy": 6.0}  # as results/step_sweep's demo3
AAA = os.environ.get("AAA", os.path.expanduser("~/amplicon-allele-analysis"))
SIM = os.path.join(AAA, "simulations")
DEMOS = {"demo5": "runs/demo5_cycle1_ecdna/subgraph/amplicon_0", "demo6": "runs/demo6_cycle2_ecdna/subgraph_union/union",
         "demo7": "runs/demo7_nestedcycles_ecdna/subgraph_union/union", "demo8": "runs/demo8_nestedcycles_ecdna/subgraph_union/union"}
BENCH = ["ctrl_dup_s710", "nest_samehap_6to1_s804", "ctrl_separate_s603", "nest_cross_1to1_s505"]  # imputed users, chr15+chr7, fold-back


def real_cases():
    out = {k: os.path.join(SIM, v) for k, v in DEMOS.items()}
    out.update({k: os.path.join(SIM, "subgraphs", "bench_nested", k, "union") for k in BENCH})
    return {k: v for k, v in out.items() if os.path.exists(os.path.join(v, "edges.tsv"))}


class TmpDir(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def path(self, *parts):
        return os.path.join(self.dir, *parts)

    def write(self, name, text):
        with open(self.path(name), "w") as handle:
            handle.write(text)
        return self.path(name)


class TestGraphWriter(TmpDir):
    def test_write_read_round_trip(self):
        seqs = [aaf.AASequence("chr2", 101, 200, 10.25, coverage=150.5, size=100, reads=0),
                aaf.AASequence("chr2", 201, 300, 4.0, coverage=60.0)]
        edges = [aaf.AAEdge("concordant", ("chr2", 200, "+"), ("chr2", 201, "-"), 1.0 / 3.0, 2.5),
                 aaf.AAEdge("discordant", ("chr2", 101, "-"), ("chr2", 300, "+"), 6.1, 40.0)]
        path = self.path("x_amplicon1_graph.txt")
        aaf.write_graph(path, seqs, edges, [("e1+:1,d1-:1,e2+:1", 7)], [("chr2", 101, 300)], header="#test", comments=["#rck_x\ta\tb"])
        g = aaf.read_graph(path)
        self.assertEqual([(s.chrom, s.start, s.end, s.cn, s.coverage, s.size) for s in g.sequences],
                         [("chr2", 101, 200, 10.25, 150.5, 100), ("chr2", 201, 300, 4.0, 60.0, 100)])
        self.assertEqual([(e.kind, e.a, e.b, e.cn, e.reads) for e in g.edges], [(e.kind, e.a, e.b, e.cn, e.reads) for e in edges])
        self.assertEqual(g.path_constraints, [("e1+:1,d1-:1,e2+:1", "7")])
        self.assertEqual(g.intervals, [("chr2", 101, 300)])
        self.assertEqual(g.rck_lines("rck_x"), [["a", "b"]])
        with open(path) as handle:
            lines = handle.read().split("\n")[:-1]
        self.assertTrue(all(line.strip() for line in lines))  # no blank line: CoRAL 3.1.1 crashes on one
        self.assertTrue(all(len(line.split("\t")) == 7 for line in lines if line.startswith("sequence")))
        self.assertTrue(all(len(line.split("\t")) == 6 for line in lines if line.startswith(("concordant", "discordant"))))
        self.assertTrue(all(len(line.split("\t")) == 3 for line in lines if line.startswith("path_constraint")))  # CycleExtractor

    def test_refuses_bad_comment(self):
        with self.assertRaises(ValueError):
            aaf.write_graph(self.path("y_graph.txt"), [], [], comments=["no hash"])

    def test_coral_path_string(self):
        # CoRAL reverses a path whose first edge index is larger than its last, and counts repeats
        self.assertEqual(aaf.coral_path_string([("e", 0), "+", ("d", 1), "-", ("e", 2)]), "e1+:1,d2-:1,e3+:1")
        self.assertEqual(aaf.coral_path_string([("e", 2), "-", ("d", 1), "+", ("e", 0)]), "e1+:1,d2-:1,e3+:1")
        self.assertEqual(aaf.coral_path_string([("e", 3), "-", ("d", 0), "+", ("e", 3)]), "e4-:2,d1+:1,e4-:2")
        self.assertEqual(aaf.parse_path_tokens("e4-:2,d1+:1"), [("e", 3, "-", 2), ("d", 0, "+", 1)])

    def test_chrom_key_is_corals(self):
        self.assertEqual(sorted(["chr15", "chr7", "chrX", "chr1", "chrM", "chrUn_x"], key=aaf.chrom_key),
                         ["chr1", "chr7", "chr15", "chrX", "chrUn_x", "chrM"])  # others before chrM, as CoRAL


class TestFixtureRoundTrip(TmpDir):
    """amplicon_0 (two chromosomes, imputed candidates, flanks) with synthetic read paths."""

    def setUp(self):
        super().setUp()
        self.sub = self.path("amplicon_0")
        shutil.copytree(AMP0, self.sub)
        shutil.copyfile(os.path.join(GOLDEN, "amplicon0_paths.tsv"), os.path.join(self.sub, "paths.tsv"))
        self.out = export_subgraph(self.sub, self.path("aa"), name="amplicon0", **AMP0_SCALE)
        self.graph = self.out["graph"]

    def test_tokens_match_coral(self):
        """Our path_constraint tokens are the ones CoRAL 3.1.1's path_to_str writes for the same reads."""
        golden = pd.read_csv(os.path.join(GOLDEN, "coral311_amplicon0_path_tokens.tsv"), sep="\t", dtype=str).fillna("")
        g = aaf.read_graph(self.graph)
        pids = [f[2] for f in g.rck_lines("rck_pc")]
        mine = dict(zip(pids, [t for t, _ in g.path_constraints]))
        for pid, tokens in zip(golden["PID"], golden["coral_tokens"]):
            self.assertEqual(mine.get(pid, ""), tokens, pid)
        self.assertEqual(self.out["info"]["paths"], {"written": 4, "outer_end_outside_graph": 2, "no_junction": 1})

    def test_layout(self):
        g = aaf.read_graph(self.graph)
        self.assertTrue(os.path.basename(self.graph).endswith("_amplicon1_graph.txt"))  # coral cycle reads N from the name
        self.assertEqual([s.chrom for s in g.sequences], sorted([s.chrom for s in g.sequences], key=aaf.chrom_key))
        self.assertEqual(g.intervals, [("chr4", 6019002, 6690001), ("chr6", 69527002, 70857000)])
        sub = read_subgraph(self.sub)
        self.assertEqual(len(g.rck_lines("rck_imputed")), sum(a.kind == IMPUTED for a in sub.adjacencies.values()))
        self.assertEqual(len(g.discordant), sum(a.kind == NOVEL for a in sub.adjacencies.values()))
        # junction CN is what rck-ecdna uses: reads / junction reads per copy
        self.assertTrue(all(abs(e.cn - e.reads / 6.0) < 1e-12 for e in g.discordant))
        info = json.load(open(self.out["export"]))
        self.assertEqual(info["counts"]["outside"], 4)  # the four region-flank joins

    def test_comment_is_exact(self):
        import_graph(self.graph, self.path("rt"), imputed="comment")
        self.assertEqual(compare_subgraphs(self.sub, self.path("rt"), original_scale=AMP0_SCALE), [])

    def test_none_drops_imputed_only(self):
        info = import_graph(self.graph, self.path("rt"), imputed="none")
        g = read_subgraph(self.path("rt"))
        self.assertFalse(any(a.kind == IMPUTED for a in g.adjacencies.values()))
        self.assertEqual(info["imputed_in_graph"], 12)
        self.assertEqual(compare_subgraphs(self.sub, self.path("rt"), ignore_imputed=True, original_scale=AMP0_SCALE), [])

    def test_drop_writes_no_imputed(self):
        out = export_subgraph(self.sub, self.path("aa2"), name="amplicon0", imputed="drop", **AMP0_SCALE)
        self.assertEqual(aaf.read_graph(out["graph"]).rck_lines("rck_imputed"), [])

    def test_regenerate(self):
        info = import_graph(self.graph, self.path("rt"), imputed="regenerate")
        g = read_subgraph(self.path("rt"))
        imputed = [a for a in g.adjacencies.values() if a.kind == IMPUTED]
        self.assertGreater(info["counts"]["regenerated_intra"], 0)
        for a in imputed:  # a closing candidate: from a segment end (+) back to a segment start (-)
            self.assertEqual({g.nodes[a.u].side, g.nodes[a.v].side}, {"+", "-"})
        self.assertEqual(g.validate(), [])

    def test_ids_ignored_graph_still_imports(self):
        """Without the #rck lines (any AA / CoRAL graph) the import names nodes by EID and keeps the structure."""
        with open(self.graph) as handle:
            plain = [line for line in handle if not line.startswith("#rck")]
        path = self.write("plain_amplicon1_graph.txt", "".join(plain))
        info = import_graph(path, self.path("plain"), imputed="none")
        self.assertFalse(info["rck_ids"])
        g, sub = read_subgraph(self.path("plain")), read_subgraph(self.sub)
        self.assertEqual(g.validate(), [])
        self.assertEqual(sorted((s.chrom, s.start, s.end) for s in g.segments.values()), sorted((s.chrom, s.start, s.end) for s in sub.segments.values()))
        self.assertEqual(sum(a.kind == NOVEL for a in g.adjacencies.values()), sum(a.kind == NOVEL for a in sub.adjacencies.values()))
        nodes = pd.read_csv(self.path("plain", "nodes.tsv"), sep="\t")
        self.assertTrue(all(n.endswith(("_+", "_-")) for n in nodes["NODE_ID"]))
        self.assertEqual(info["paths"]["path_constraint"], 4)


class TestHandBuilt(TmpDir):
    TEXT = ("SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped\n"
            "sequence\tchr3:1001-\tchr3:2000+\t20.0\t300.0\t1000\t0\n"
            "sequence\tchr3:2001-\tchr3:3000+\t10.0\t150.0\t1000\t0\n"
            "BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, H, S\n"
            "source\tchr3:-1-->chr3:1001-\t2.0\t3\tNone\tNone\n"
            "concordant\tchr3:2000+->chr3:2001-\t10.0\t{c}\tNone\tNone\n"
            "discordant\tchr3:1001-->chr3:2000+\t10.0\t{d}\tNone\tNone\n")

    def test_scale_and_source(self):
        info = import_graph(self.write("h_amplicon1_graph.txt", self.TEXT.format(c=50, d=60)), self.path("h"))
        self.assertEqual(info["scale"]["junction_reads_per_copy"], 6.0)
        self.assertEqual(info["scale"]["depth_per_copy"], 15.0)
        self.assertEqual(info["source_edges_dropped"], ["chr3:-1-->chr3:1001-"])
        edges = pd.read_csv(self.path("h", "edges.tsv"), sep="\t")
        self.assertEqual(sorted(edges["EDGE_TYPE"]), ["extremity", "reference", "reference", "variant"])
        g = read_subgraph(self.path("h"))
        self.assertEqual([a.u == a.v for a in g.adjacencies.values() if a.kind == NOVEL], [False])
        self.assertEqual(g.terminals(), sorted(n for n in g.nodes if n in (g.segments["s0001"].tail, g.segments["s0002"].head)))

    def test_no_reads_fallback(self):
        info = import_graph(self.write("h_amplicon1_graph.txt", self.TEXT.format(c=0, d=0)), self.path("h"))
        self.assertEqual(info["scale"]["reads_mode"], "cn")
        self.assertEqual(info["scale"]["junction_reads_per_copy"], 1.0)
        edges = pd.read_csv(self.path("h", "edges.tsv"), sep="\t")
        self.assertEqual(list(edges.loc[edges["EDGE_TYPE"] == "variant", "WEIGHT"]), [10.0])  # WEIGHT = CN

    def test_coral_rewrite_not_used_for_depth(self):
        info = import_graph(self.write("h_amplicon1_graph.txt", "#CoRAL 3.1.1\n" + self.TEXT.format(c=50, d=60)), self.path("h"))
        self.assertTrue(info["coral_rewritten"])
        self.assertIn("not fitted", info["scale"]["depth_method"])

    def test_noncontiguous_concordant_rejected(self):
        text = self.TEXT.format(c=5, d=6).replace("chr3:2000+->chr3:2001-", "chr3:2000+->chr3:2501-")
        with self.assertRaises(ValueError):
            import_graph(self.write("h_amplicon1_graph.txt", text), self.path("h"))

    def test_foldback_kept(self):
        text = self.TEXT.format(c=50, d=60) + "discordant\tchr3:3000+->chr3:3000+\t4.0\t24\tNone\tNone\n"
        import_graph(self.write("h_amplicon1_graph.txt", text), self.path("h"))
        g = read_subgraph(self.path("h"))
        self.assertEqual(sum(a.is_self_loop for a in g.adjacencies.values() if a.kind == NOVEL), 1)


@unittest.skipUnless(real_cases(), "amplicon-allele-analysis runs not present")
class TestRealRoundTrip(TmpDir):
    """demo5-demo8 and bench cases (incl. the imputed-edge users ctrl_dup_s710 and nest_samehap_6to1_s804)."""

    def round_trip(self, name, sub):
        out = export_subgraph(sub, self.path(name, "aa"), name=name)
        scale = json.load(open(os.path.join(sub, "scale.json")))
        common = dict(cn_profile=os.path.join(sub, "cn_profile.tsv"), background_cn=scale["background_cn"], sample=scale["sample"])
        import_graph(out["graph"], self.path(name, "comment"), imputed="comment", **common)
        import_graph(out["graph"], self.path(name, "none"), imputed="none", **common)
        return out

    def test_cases(self):
        for name, sub in real_cases().items():
            with self.subTest(case=name):
                out = self.round_trip(name, sub)
                self.assertEqual(compare_subgraphs(sub, self.path(name, "comment")), [])
                self.assertEqual(compare_subgraphs(sub, self.path(name, "none"), ignore_imputed=True), [])
                self.assertEqual(out["info"]["paths"].get("outside_graph", 0) + out["info"]["counts"]["rck_path"], 0)

    def test_imputed_users(self):
        used = {"ctrl_dup_s710": "imputed_adj_intra_chr2_10", "nest_samehap_6to1_s804": "imputed_adj_intra_chr3_10"}
        cases = real_cases()
        for name, aid in used.items():
            if name not in cases:
                continue
            with self.subTest(case=name):
                self.round_trip(name, cases[name])
                self.assertIn(aid, read_subgraph(self.path(name, "comment")).adjacencies)
                self.assertNotIn(aid, read_subgraph(self.path(name, "none")).adjacencies)

    def test_foldbacks_listed(self):
        cases = real_cases()
        if "ctrl_dup_s710" in cases:
            out = export_subgraph(cases["ctrl_dup_s710"], self.path("fb"), name="ctrl_dup_s710")
            self.assertEqual(out["info"]["foldbacks"], ["sim00005", "sim00008"])
            g = aaf.read_graph(out["graph"])
            self.assertEqual(sum(e.a == e.b for e in g.discordant), 2)

    def test_cli_dispatch(self):
        cases = real_cases()
        name, sub = next(iter(cases.items()))
        self.assertEqual(export_main(["subgraph", sub, self.path("cli"), "--name", name]), 0)
        self.assertTrue(os.path.exists(self.path("cli", name + "_amplicon1_graph.txt")))


if __name__ == "__main__":
    unittest.main()

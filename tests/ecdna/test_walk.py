import os
import unittest

import pandas as pd

from rck.ecdna.phased import adj_key, is_balanced, seg_key
from rck.ecdna.subgraph_io import read_subgraph
from rck.ecdna.walk import (ReadPath, alternating_circuit, canonical, enumerate_circuits, path_support, read_paths,
                            reverse_walk, walk_counts)
from tests.ecdna.helpers import head, loop_counts, loops_graph, single_circle, SubgraphBuilder, tail, temp_dir


def circle_counts(g):
    counts = {seg_key(sid, "A"): 1 for sid in g.segments}
    for a in g.adjacencies.values():
        if a.aid == "circ" or (a.kind == "reference" and a.reads and a.reads > 50):
            counts[adj_key(a.aid, "AA")] = 1
    return counts


class TestAlternatingCircuit(unittest.TestCase):
    def test_single_circle(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g = read_subgraph(directory)
        counts = circle_counts(g)
        self.assertEqual(is_balanced(g, counts), [])
        steps = alternating_circuit(g, counts)
        self.assertEqual(len(steps), 3)
        self.assertEqual(walk_counts(steps), counts)
        self.assertEqual(len({s.strand for s in steps}), 1)

    def test_canonical_is_rotation_and_reversal_invariant(self):
        with temp_dir() as tmp:
            directory, _ = single_circle(tmp)
            g = read_subgraph(directory)
        steps = alternating_circuit(g, circle_counts(g))
        key = canonical(steps)
        self.assertEqual(canonical(steps[1:] + steps[:1]), key)
        self.assertEqual(canonical(reverse_walk(steps)), key)
        self.assertEqual(walk_counts(reverse_walk(steps)), walk_counts(steps))

    def test_tandem_duplicated_segment(self):
        with temp_dir() as tmp:
            b = SubgraphBuilder()
            segs = b.chain("chr1", [0, 100_000, 200_000, 300_000], [100.0] * 3)
            b.variant(segs[1][1], segs[1][0], 50, aid="dup")
            b.variant(segs[2][1], segs[0][0], 50, aid="circ")
            g = read_subgraph(b.write(tmp))
        sids = sorted(g.segments, key=lambda s: g.segments[s].start)
        ref = {frozenset((a.u, a.v)): a.aid for a in g.adjacencies.values() if a.kind == "reference"}
        counts = {seg_key(sids[0], "A"): 1, seg_key(sids[1], "A"): 2, seg_key(sids[2], "A"): 1,
                  adj_key(ref[frozenset((head("chr1", 100_000), tail("chr1", 100_000)))], "AA"): 1,
                  adj_key(ref[frozenset((head("chr1", 200_000), tail("chr1", 200_000)))], "AA"): 1,
                  adj_key("dup", "AA"): 1, adj_key("circ", "AA"): 1}
        self.assertEqual(is_balanced(g, counts), [])
        steps = alternating_circuit(g, counts)
        self.assertEqual([s.sid for s in steps].count(sids[1]), 2)
        self.assertEqual(walk_counts(steps), counts)
        self.assertEqual(len(enumerate_circuits(g, counts)), 1)

    def test_three_loops_have_two_orderings(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp)
            g = read_subgraph(tmp)
        counts = loop_counts(g, ends, ["a", "b", "c"])
        self.assertEqual(is_balanced(g, counts), [])
        orderings = enumerate_circuits(g, counts)
        self.assertEqual(len(orderings), 2)
        for steps in orderings:
            self.assertEqual(walk_counts(steps), counts)
        self.assertEqual(len({canonical(s) for s in orderings}), 2)

    def test_balance_checks(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp, loops=("a", "b"))
            g = read_subgraph(tmp)
        counts = loop_counts(g, ends, ["a"])
        counts_b = loop_counts(g, ends, ["b"])
        merged = dict(counts)
        for k, v in counts_b.items():
            merged[k] = merged.get(k, 0) + v
        self.assertEqual(is_balanced(g, merged), [])
        self.assertEqual(len(alternating_circuit(g, merged)), 4)  # connected through x: x a x b
        unbalanced = dict(counts)
        unbalanced[seg_key(g.nodes[ends["a"][0]].segment, "A")] += 1
        self.assertNotEqual(is_balanced(g, unbalanced), [])
        with self.assertRaises(ValueError):
            alternating_circuit(g, unbalanced)


def linear_chain_counts(g):
    """A linear walk over every segment of a chain (haplotype A), ending at both flanks."""
    from rck.ecdna.phased import end_key

    sids = sorted(g.segments, key=lambda x: g.segments[x].start)
    counts = {seg_key(sid, "A"): 1 for sid in sids}
    for a in g.adjacencies.values():
        if a.kind == "reference" and g.nodes[a.u].segment is not None and g.nodes[a.v].segment is not None:
            counts[adj_key(a.aid, "AA")] = 1
    counts[end_key(g.segments[sids[0]].tail, "A")] = 1
    counts[end_key(g.segments[sids[-1]].head, "A")] = 1
    return counts


class TestLinearWalks(unittest.TestCase):
    def build(self, tmp):
        b = SubgraphBuilder()
        b.chain("chr1", [0, 1000, 2000, 3000], [100.0] * 3)
        return read_subgraph(b.write(tmp))

    def test_open_walk_from_end_to_end(self):
        from rck.ecdna.walk import is_linear

        with temp_dir() as tmp:
            g = self.build(tmp)
        counts = linear_chain_counts(g)
        self.assertEqual(is_balanced(g, counts), [])
        steps = alternating_circuit(g, counts)
        self.assertEqual(len(steps), 3)
        self.assertTrue(is_linear(steps))
        self.assertEqual(steps[0].entry[0], "T")
        self.assertEqual(steps[-1].join[0], "T")
        self.assertEqual(walk_counts(steps), counts)
        self.assertEqual(len({s.strand for s in steps}), 1)

    def test_linear_canonical_and_enumeration(self):
        with temp_dir() as tmp:
            g = self.build(tmp)
        counts = linear_chain_counts(g)
        steps = alternating_circuit(g, counts)
        self.assertEqual(canonical(reverse_walk(steps)), canonical(steps))
        self.assertEqual(walk_counts(reverse_walk(steps)), counts)
        self.assertEqual(len(enumerate_circuits(g, counts)), 1)

    def test_one_end_is_not_a_walk(self):
        with temp_dir() as tmp:
            g = self.build(tmp)
        counts = linear_chain_counts(g)
        broken = {k: v for k, v in counts.items() if not (k[0] == "T" and k[1] == max(x[1] for x in counts if x[0] == "T"))}
        with self.assertRaises(ValueError):
            alternating_circuit(g, broken)


class TestPathSupport(unittest.TestCase):
    def test_paths_rank_orderings(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp)
            g = read_subgraph(tmp)
            a_t, a_h = ends["a"]
            x_t, x_h = ends["x"]
            b_t, b_h = ends["b"]

            def tok(nid):
                eid, orient = nid.rsplit("_", 1)
                return "{e}@{o}".format(e=eid, o=orient)

            path_tsv = os.path.join(tmp, "path.tsv")
            pd.DataFrame([{"PID": "p1", "PATH": ",".join(tok(n) for n in (a_t, a_h, x_t, x_h, b_t, b_h)), "NSTEP": 6, "#NOV": 40},
                          {"PID": "p2", "PATH": ",".join(tok(n) for n in ("chr9:5_-", "chr9:7_+")), "NSTEP": 2, "#NOV": 3}]).to_csv(
                path_tsv, sep="\t", index=False)
            paths = read_paths(path_tsv, g)
        self.assertEqual([p.pid for p in paths], ["p1"])
        self.assertEqual(paths[0].junctions, [(a_h, x_t), (x_h, b_t)])
        counts = loop_counts(g, ends, ["a", "b", "c"])
        scores = sorted(path_support(g, steps, paths)[0] for steps in enumerate_circuits(g, counts))
        self.assertEqual(scores, [0.0, 40.0])

    def test_single_junction_paths_do_not_score(self):
        with temp_dir() as tmp:
            ends = loops_graph(tmp)
            g = read_subgraph(tmp)
        counts = loop_counts(g, ends, ["a", "b", "c"])
        steps = enumerate_circuits(g, counts)[0]
        p = ReadPath(pid="one", support=10, junctions=[(ends["x"][1], ends["a"][0])])
        score, matched, applicable = path_support(g, steps, [p])
        self.assertEqual((score, matched, applicable), (0.0, 1, 1))


if __name__ == "__main__":
    unittest.main()

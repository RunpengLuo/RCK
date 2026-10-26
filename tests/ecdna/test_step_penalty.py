"""--step-penalty (walk complexity) and the ambiguity report of the top decomposition.

Fixture: a hub segment ``x`` (chr1) and two loops ``a`` and ``b`` of three segments each, on
chr2 and chr3, joined to ``x`` by junctions ``x>a``, ``a>x``, ``x>b`` and ``b>x``. Each loop is a
circle through ``x``, so the two species share ``x`` as nested ecDNA share sequence. Steps
(segment traversals): (x a) 4, (x b) 4, (x a x a x b) 12, (x a x b) 8.

* demo7-like, 2:1 (``a`` at 16, ``b`` at 8): the merged walk (x a x a x b) at 8 is one species
  with 12 steps; the separate circles are two species with 8. Same flow, so the species penalty
  merges them unless step_penalty x 4 extra steps > 1.
* demo8-like, 3:1 (``a`` at 18, ``b`` at 6): separate {(x a) 18, (x b) 6} and merged
  {(x a) 12, (x a x b) 6} both have two species and the same flow; only steps (8 vs 12) differ.
  A one-species walk would need ``x`` four times, above --max-multiplicity 3.
"""

import json
import os
import unittest

import pandas as pd

from rck.ecdna.phased import adj_key, seg_key
from rck.ecdna.subgraph_io import read_subgraph
from tests.ecdna.helpers import requires_gurobi, SubgraphBuilder, temp_dir


def nested_graph(directory):
    b = SubgraphBuilder()
    ends = {"x": b.chain("chr1", [0, 100000], [300.0])}
    for chrom, name in (("chr2", "a"), ("chr3", "b")):
        ends[name] = b.chain(chrom, [0, 100000, 200000, 300000], [300.0] * 3)
    for name in ("a", "b"):
        b.variant(ends["x"][0][1], ends[name][0][0], 50, hap="A", aid="x>{n}".format(n=name))
        b.variant(ends[name][-1][1], ends["x"][0][0], 50, hap="A", aid="{n}>x".format(n=name))
    return read_subgraph(b.write(directory)), ends


def circle_flow(g, ends, name, weight):
    """Haplotype-A flow of the circle x -> loop ``name`` -> x at ``weight`` copies."""
    flow = {}
    sids = [g.nodes[t].segment for t, _ in ends["x"] + ends[name]]
    for sid in sids:
        flow[seg_key(sid, "A")] = weight
    for a in g.adjacencies.values():
        inside = a.kind == "reference" and g.nodes[a.u].segment in sids and g.nodes[a.v].segment in sids
        if inside or a.aid in ("x>{n}".format(n=name), "{n}>x".format(n=name)):
            flow[adj_key(a.aid, "AA")] = weight
    return flow


def add(*flows):
    total = {}
    for f in flows:
        for k, v in f.items():
            total[k] = total.get(k, 0) + v
    return total


@requires_gurobi
class TestStepPenalty(unittest.TestCase):
    def decompose(self, flow, g, **params):
        from rck.ecdna.decompose import DecompositionParams, decompose

        return decompose(g, flow, DecompositionParams(**params))

    def setUp(self):
        self.tmp = temp_dir()
        self.g, self.ends = nested_graph(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def flow(self, wa, wb):
        return add(circle_flow(self.g, self.ends, "a", wa), circle_flow(self.g, self.ends, "b", wb))

    def test_two_to_one_merges_without_step_penalty(self):
        best = self.decompose(self.flow(16, 8), self.g)[0]
        self.assertEqual([s.weight for s in best.species], [8])  # one walk: x a x a x b
        self.assertEqual(best.steps, 12)
        amb = best.ambiguity
        self.assertTrue(amb["ambiguous"])
        self.assertEqual(amb["decided_by"], "species_penalty")
        separate = [a for a in amb["alternatives"] if a["n_species"] == 2 and a["steps"] == 8]
        self.assertEqual(len(separate), 1)
        self.assertEqual(sorted(separate[0]["weights"]), [8, 16])
        self.assertAlmostEqual(separate[0]["gap"], 1.0, delta=0.01)  # one more species, minus the tie-break

    def test_two_to_one_separates_with_step_penalty(self):
        best = self.decompose(self.flow(16, 8), self.g, step_penalty=0.5)[0]
        self.assertEqual(sorted(s.weight for s in best.species), [8, 16])
        self.assertEqual(best.steps, 8)
        self.assertAlmostEqual(best.residual, 0.0)
        amb = best.ambiguity
        self.assertTrue(amb["ambiguous"])  # the merged walk still fits the flow exactly as well
        merged = [a for a in amb["alternatives"] if a["n_species"] == 1]
        self.assertEqual(merged[0]["steps"], 12)
        self.assertAlmostEqual(merged[0]["gap"], 0.5 * 4 - 1.0, delta=0.01)

    def test_three_to_one_is_decided_by_steps(self):
        best = self.decompose(self.flow(18, 6), self.g, num_alternatives=10)[0]
        self.assertEqual(sorted(s.weight for s in best.species), [6, 18])
        self.assertEqual(best.ambiguity["decided_by"], "tie-break")  # as on demo8
        merged = [a for a in best.ambiguity["alternatives"] if a["steps"] == 12]
        self.assertEqual(sorted(merged[0]["weights"]), [6, 12])
        self.assertLess(merged[0]["gap"], 0.01)
        strong = self.decompose(self.flow(18, 6), self.g, step_penalty=0.1)[0]
        self.assertEqual(sorted(s.weight for s in strong.species), [6, 18])
        self.assertEqual(strong.ambiguity["decided_by"], "step_penalty")
        gap = [a for a in strong.ambiguity["alternatives"] if a["steps"] == 12][0]["gap"]
        self.assertAlmostEqual(gap, 0.1 * 4, delta=0.01)

    def test_default_ranking_unchanged(self):
        """step_penalty 0 keeps today's ranking: fewest significant species first."""
        from rck.ecdna.decompose import DecompositionParams, structure_cost

        decs = self.decompose(self.flow(16, 8), self.g, num_alternatives=10, pool_gap=2.0)
        counts = [sum(1 for s in d.species if not s.low_weight) for d in decs]
        self.assertEqual(counts, sorted(counts))
        costs = [structure_cost(d.species, DecompositionParams()) for d in decs]
        self.assertEqual(costs, sorted(costs))

    def test_single_species_is_not_ambiguous(self):
        best = self.decompose(circle_flow(self.g, self.ends, "a", 10), self.g)[0]
        self.assertEqual([s.weight for s in best.species], [10])
        self.assertFalse(best.ambiguity["ambiguous"])
        self.assertEqual(best.ambiguity["alternatives"], [])

    def test_ambiguity_off(self):
        best = self.decompose(self.flow(18, 6), self.g, ambiguity_window=0)[0]
        self.assertIsNone(best.ambiguity)


class TestPoolRanking(unittest.TestCase):
    """rank_pool on recorded entries (no Gurobi)."""

    def entry(self, index, objective, n_significant, steps, residual=0.0, n_species=None):
        return {"index": index, "objective": objective, "residual": residual, "n_significant": n_significant,
                "n_species": n_significant if n_species is None else n_species, "steps": steps}

    def test_step_penalty_reorders_and_prunes(self):
        from rck.ecdna.decompose import DecompositionParams, rank_pool

        merged = self.entry(0, 1.001, 1, 12)  # one species, merged walk
        separate = self.entry(1, 2.000, 2, 8)  # two circles: one more species, 4 fewer steps
        far = self.entry(2, 3.5, 2, 8)  # 2.5 worse: beyond pool_gap, and no step saving pays for it
        entries = [merged, separate, far]
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=0.0)), [0])  # λ = 0: within pool_gap only
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=0.1))[0], 0)  # 0.4 < 1 species
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=0.5)), [1, 0])  # 2.0 > 1 species

    def test_residual_first(self):
        from rck.ecdna.decompose import DecompositionParams, rank_pool

        entries = [self.entry(0, 0.5, 1, 4, residual=0.5), self.entry(1, 0.52, 2, 30)]
        self.assertEqual(rank_pool(entries, DecompositionParams(step_penalty=1.0))[0], 1)


@requires_gurobi
class TestKeepPool(unittest.TestCase):
    """A pool kept at one λ re-ranks to what a run at another λ ranks first."""

    LAMBDAS = (0.02, 0.05, 0.1, 0.2, 0.5)

    def setUp(self):
        self.tmp = temp_dir()
        self.g, self.ends = nested_graph(self.tmp.name)

    def tearDown(self):
        self.tmp.cleanup()

    def check(self, flow):
        from rck.ecdna.decompose import (DecompositionParams, ambiguity_summary, decompose, decompose_with_pool,
                                         rank_pool)

        base = DecompositionParams(step_penalty=0.1)
        ranked, record = decompose_with_pool(self.g, flow, base, pool_lambdas=(0.0,) + self.LAMBDAS)
        self.assertEqual([d.signature for d in ranked], [d.signature for d in decompose(self.g, flow, base)])
        for lam in self.LAMBDAS:
            params = DecompositionParams(step_penalty=lam)
            direct = decompose(self.g, flow, params)[0]
            index = rank_pool(record.entries, params)[0]
            self.assertEqual(index, record.rank1[lam])
            self.assertEqual(record.decompositions[index].signature, direct.signature, lam)
            entry = record.entries[index]
            report = ambiguity_summary(entry["objective"], entry["steps"], entry["n_significant"], record.same_fit[index],
                                       params, warn=False)
            for key in ("ambiguous", "decided_by", "gap", "n_alternatives"):
                self.assertEqual(report.get(key), direct.ambiguity.get(key), (lam, key))
        return record

    def test_two_to_one(self):
        record = self.check(self.flow(16, 8))
        self.assertEqual(len({record.rank1[lam] for lam in self.LAMBDAS}), 2)  # merged at 0.02, separate from 0.5

    def test_three_to_one(self):
        self.check(self.flow(18, 6))

    def flow(self, wa, wb):
        return add(circle_flow(self.g, self.ends, "a", wa), circle_flow(self.g, self.ends, "b", wb))

    def test_command_line_and_rerank(self):
        from rck.ecdna.rerank import rerank
        from rck.ecdna_run import main

        out = os.path.join(self.tmp.name, "out")
        argv = ["--subgraph-dir", os.path.join(self.tmp.name), "--depth-per-copy", "15", "--junction-reads-per-copy", "6",
                "--no-plots", "--step-penalty", "0.1"]
        with self.assertRaises(SystemExit):
            main(argv[:-2] + ["--outdir", out, "--keep-pool"])  # needs a step penalty above 0
        self.assertEqual(main(argv + ["--outdir", out, "--keep-pool"]), 0)
        for name in ("pool.json", "pool.tsv", "decompositions.tsv"):
            self.assertTrue(os.path.exists(os.path.join(out, "pool", name)), name)
        with open(os.path.join(out, "pool", "pool.json")) as handle:
            pool = json.load(handle)
        self.assertEqual(len(pool["entries"]), len(pd.read_csv(os.path.join(out, "pool", "pool.tsv"), sep="\t")))
        table = rerank(out, [0.1])
        top = pd.read_csv(os.path.join(out, "decompositions.tsv"), sep="\t")
        top = top[(top.decomposition_rank == 1) & (top.ordering_rank == 1) & (top.status != "low_weight")]
        self.assertEqual(table.weights.iloc[0], ",".join("%g" % w for w in top.weight))
        with open(os.path.join(out, "summary.json")) as handle:
            ambiguity = json.load(handle)["ambiguity"]
        self.assertEqual(bool(table.ambiguous.iloc[0]), ambiguity["ambiguous"])


if __name__ == "__main__":
    unittest.main()

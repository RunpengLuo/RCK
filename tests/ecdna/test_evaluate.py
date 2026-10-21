import os
import unittest

import pandas as pd

from rck.ecdna.evaluate import (CStep, Cycle, evaluate, fraction_in_regions, jaccard, junction_scores, order_match,
                                read_aa_cycles, read_truth_table, runs,
                                truth_from_simulation, write_truth_table)
from tests.ecdna.helpers import DATA_DIR, temp_dir

TRUTH = [CStep("chr1", 0, 999, "A", "+"), CStep("chr1", 5000, 5999, "A", "-"), CStep("chr2", 100, 1099, "B", "+")]


def rotate(steps, k):
    return steps[k:] + steps[:k]


def reverse(steps):
    return [CStep(s.chrom, s.start, s.end, s.hap, "-" if s.strand == "+" else "+") for s in reversed(steps)]


class TestMetrics(unittest.TestCase):
    def test_identical(self):
        self.assertEqual(jaccard(TRUTH, TRUTH), 1.0)
        self.assertEqual(junction_scores(TRUTH, TRUTH, tol=0, phased=True), (1.0, 1.0))
        self.assertTrue(order_match(TRUTH, TRUTH, tol=0))

    def test_rotation_and_reversal(self):
        self.assertTrue(order_match(TRUTH, rotate(TRUTH, 1), tol=0))
        self.assertTrue(order_match(TRUTH, reverse(TRUTH), tol=0))
        self.assertEqual(junction_scores(TRUTH, reverse(TRUTH), tol=0, phased=True), (1.0, 1.0))

    def test_finer_segmentation_is_equivalent(self):
        split = [CStep("chr1", 0, 499, "A", "+"), CStep("chr1", 500, 999, "A", "+")] + TRUTH[1:]
        self.assertEqual(len(runs(split)), 3)
        self.assertEqual(jaccard(TRUTH, split), 1.0)
        self.assertTrue(order_match(TRUTH, split, tol=0))

    def test_haplotype_errors(self):
        wrong = [CStep(s.chrom, s.start, s.end, "A", s.strand) for s in TRUTH]
        self.assertLess(jaccard(TRUTH, wrong, with_hap=True), 1.0)
        self.assertEqual(jaccard(TRUTH, wrong, with_hap=False), 1.0)
        self.assertEqual(junction_scores(TRUTH, wrong, tol=0, phased=False), (1.0, 1.0))
        self.assertLess(junction_scores(TRUTH, wrong, tol=0, phased=True)[0], 1.0)
        self.assertFalse(order_match(TRUTH, wrong, tol=0))

    def test_tolerance(self):
        shifted = [CStep(s.chrom, s.start + 3, s.end + 3, s.hap, s.strand) for s in TRUTH]
        self.assertFalse(order_match(TRUTH, shifted, tol=0))
        self.assertTrue(order_match(TRUTH, shifted, tol=5))
        self.assertEqual(junction_scores(TRUTH, shifted, tol=5, phased=True), (1.0, 1.0))

    def test_multiplicity_counts_in_jaccard(self):
        doubled = TRUTH + [CStep("chr1", 0, 999, "A", "+")]
        self.assertAlmostEqual(jaccard(TRUTH, doubled), 3000 / 4000.0)


class TestTables(unittest.TestCase):
    def test_truth_table_round_trip(self):
        with temp_dir() as tmp:
            path = os.path.join(tmp, "truth.tsv")
            write_truth_table(path, [Cycle("ec", 12.0, [TRUTH])])
            back = read_truth_table(path)
        self.assertEqual(back[0].steps, TRUTH)
        self.assertEqual(back[0].weight, 12.0)

    def test_evaluate_decompositions_table(self):
        def walk(steps):
            return ",".join("{c}:{s}-{e}:{h}:{st}".format(c=x.chrom, s=x.start, e=x.end, h=x.hap, st=x.strand) for x in steps)

        rows = [
            {"decomposition_rank": 1, "objective": 1.0, "species": "ec1", "weight": 11, "ordering_rank": 1, "walk": walk(rotate(TRUTH, 2))},
            {"decomposition_rank": 2, "objective": 2.0, "species": "ec1", "weight": 6, "ordering_rank": 1, "walk": walk(TRUTH[:1] + TRUTH[2:])},
            {"decomposition_rank": 2, "objective": 2.0, "species": "ec2", "weight": 6, "ordering_rank": 1, "walk": walk(TRUTH[1:2])},
        ]
        with temp_dir() as tmp:
            path = os.path.join(tmp, "decompositions.tsv")
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            frame, summary = evaluate(path, [Cycle("truth1", 12.0, [TRUTH])], tol=0)
        self.assertEqual(summary["first_rank_with_all_orders_matched"], 1)
        self.assertEqual(summary["best_mean_jaccard_rank"], 1)
        top = frame[frame.decomposition_rank == 1].iloc[0]
        self.assertAlmostEqual(top.cn_rel_err, 1 / 12.0)
        self.assertFalse(bool(frame[frame.decomposition_rank == 2].iloc[0].order_match))


class TestLowWeight(unittest.TestCase):
    def test_low_weight_species_are_not_scored_by_default(self):
        def walk(steps):
            return ",".join("{c}:{s}-{e}:{h}:{st}".format(c=x.chrom, s=x.start, e=x.end, h=x.hap, st=x.strand) for x in steps)

        rows = [
            {"decomposition_rank": 1, "objective": 1.0, "species": "ec1", "status": "species", "weight": 12, "ordering_rank": 1,
             "walk": walk(TRUTH)},
            {"decomposition_rank": 1, "objective": 1.0, "species": "ec2", "status": "low_weight", "weight": 1, "ordering_rank": 1,
             "walk": walk(TRUTH[:1])},
        ]
        with temp_dir() as tmp:
            path = os.path.join(tmp, "decompositions.tsv")
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            frame, _ = evaluate(path, [Cycle("truth1", 12.0, [TRUTH])], tol=0)
            frame_all, _ = evaluate(path, [Cycle("truth1", 12.0, [TRUTH])], tol=0, include_low_weight=True)
        self.assertEqual(int(frame.iloc[0].n_species_inferred), 1)
        self.assertEqual(int(frame.iloc[0].n_low_weight), 1)
        self.assertEqual(int(frame_all.iloc[0].n_species_inferred), 2)


class TestSimulatorTruth(unittest.TestCase):
    """demo3_clone1_ecdna: one episome species, 14 segments, 16 steps, phi 16."""

    @classmethod
    def setUpClass(cls):
        cls.cycles = truth_from_simulation(os.path.join(DATA_DIR, "demo3_truth"))

    def test_species(self):
        self.assertEqual(len(self.cycles), 1)
        c = self.cycles[0]
        self.assertEqual(c.name, "ecdna1")
        self.assertEqual(c.weight, 16.0)
        self.assertEqual(len(c.steps), 16)
        self.assertEqual(sum(s.end - s.start + 1 for s in c.steps), 2_343_000)

    def test_haplotypes(self):
        steps = self.cycles[0].steps
        seg4 = [s for s in steps if (s.chrom, s.start) == ("chr4", 6344000)]
        self.assertEqual(sorted(s.hap for s in seg4), ["A", "B"])
        seg5 = [s for s in steps if (s.chrom, s.start) == ("chr4", 6345000)]
        self.assertEqual([s.hap for s in seg5], ["A", "A"])
        seg10 = [s for s in steps if (s.chrom, s.start) == ("chr6", 70090000)]
        self.assertEqual([s.hap for s in seg10], ["B"])
        others = [s for s in steps if s not in seg4 + seg10]
        self.assertTrue(all(s.hap == "A" for s in others))

    def test_truth_matches_itself(self):
        c = self.cycles[0]
        self.assertTrue(order_match(c.steps, rotate(c.steps, 5), tol=0))
        # seg1 / seg2 and seg12 / seg13 are joined by reference joins, so 16 steps form 14 runs
        self.assertEqual(len(runs(c.steps)), 14)


DEMO4_TRUTH = os.path.join(DATA_DIR, "demo4_truth")


def walk_text(steps):
    return ",".join("{c}:{s}-{e}:{h}:{st}".format(c=x.chrom, s=x.start, e=x.end, h=x.hap, st=x.strand) for x in steps)


class TestSimulatorTruthAA(unittest.TestCase):
    """demo4_clone3_ecdna: two species from truth/aa cycles; ecdna1 lost 3 of its 11 segments to a deletion."""

    def test_two_species_with_per_sample_weights(self):
        t1 = {c.name: c for c in truth_from_simulation(DEMO4_TRUTH, sample="tumor1")}
        t2 = {c.name: c for c in truth_from_simulation(DEMO4_TRUTH, sample="sim_tumor2")}  # windows-style name
        self.assertEqual(sorted(t1), ["ecdna1", "ecdna2"])
        self.assertAlmostEqual(t1["ecdna1"].weight, 6.2821, places=3)
        self.assertAlmostEqual(t2["ecdna1"].weight, 7.4376, places=3)
        e1, e2 = t1["ecdna1"], t1["ecdna2"]
        self.assertEqual(len(e1.steps), 8)
        self.assertEqual(sum(s.end - s.start + 1 for s in e1.steps), 1_774_000)
        self.assertEqual({s.hap for s in e1.steps}, {"B"})
        self.assertEqual(len(e2.steps), 4)
        self.assertEqual(sum(s.end - s.start + 1 for s in e2.steps), 1_998_000)
        self.assertEqual({s.hap for s in e2.steps}, {"A"})
        # step coordinates are 0-based inclusive, matching obs/adjacency.tsv extremity names
        first = e2.steps[0]
        self.assertEqual((first.chrom, first.start, first.end, first.strand), ("chr16", 25455000, 25637999, "+"))

    def test_read_aa_cycles(self):
        parsed = read_aa_cycles(os.path.join(DEMO4_TRUTH, "truth", "aa", "ecdna1_cycles.txt"))
        self.assertEqual(len(parsed), 1)
        rows, walk = parsed[0]
        self.assertEqual(len(walk), 8)
        self.assertEqual(len(rows), 8)  # only the segments the cycle uses

    def test_missing_truth_dir_is_a_clear_error(self):
        with temp_dir() as tmp:
            with self.assertRaises(ValueError):
                truth_from_simulation(tmp)


class TestScope(unittest.TestCase):
    """Each demo4 subgraph should be scored only against the species it contains."""

    @classmethod
    def setUpClass(cls):
        from rck.ecdna.subgraph_io import read_subgraph

        cls.truth = {c.name: c for c in truth_from_simulation(DEMO4_TRUTH, sample="tumor1")}
        cls.regions = {a: read_subgraph(os.path.join(DATA_DIR, "demo4", a)).regions() for a in ("amplicon_0", "amplicon_1")}

    def test_fractions(self):
        r0, r1 = self.regions["amplicon_0"], self.regions["amplicon_1"]
        self.assertGreater(fraction_in_regions(self.truth["ecdna1"].steps, r0, pad=10_000), 0.95)
        self.assertLess(fraction_in_regions(self.truth["ecdna2"].steps, r0, pad=10_000), 0.5)
        self.assertGreater(fraction_in_regions(self.truth["ecdna2"].steps, r1, pad=10_000), 0.95)
        self.assertLess(fraction_in_regions(self.truth["ecdna1"].steps, r1, pad=10_000), 0.05)

    def test_out_of_scope_species_are_not_counted_as_missed(self):
        e1 = self.truth["ecdna1"]
        rows = [{"decomposition_rank": 1, "objective": 1.0, "species": "ec1", "status": "species", "weight": 6,
                 "ordering_rank": 1, "walk": walk_text(e1.steps)}]
        with temp_dir() as tmp:
            path = os.path.join(tmp, "decompositions.tsv")
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            frame, summary = evaluate(path, list(self.truth.values()), regions=self.regions["amplicon_0"])
            _, unscoped = evaluate(path, list(self.truth.values()))
        self.assertEqual(summary["n_truth_species"], 1)
        self.assertEqual(summary["n_truth_species_total"], 2)
        self.assertEqual(summary["first_rank_with_all_orders_matched"], 1)
        self.assertEqual(list(frame["truth"]), ["ecdna1"])
        scope = {sp["species"]: sp["in_scope"] for sp in summary["truth_species"]}
        self.assertEqual(scope, {"ecdna1": True, "ecdna2": False})
        # without regions both species are scored, so the missing ecdna2 blocks a full match
        self.assertEqual(unscoped["n_truth_species"], 2)
        self.assertIsNone(unscoped["first_rank_with_all_orders_matched"])

    def test_no_species_in_scope(self):
        rows = [{"decomposition_rank": 1, "objective": 1.0, "species": "ec1", "status": "species", "weight": 6,
                 "ordering_rank": 1, "walk": walk_text(self.truth["ecdna1"].steps)}]
        with temp_dir() as tmp:
            path = os.path.join(tmp, "decompositions.tsv")
            pd.DataFrame(rows).to_csv(path, sep="\t", index=False)
            frame, summary = evaluate(path, list(self.truth.values()), regions=[("chr1", 1, 1000)])
        self.assertEqual(summary["n_truth_species"], 0)
        self.assertIsNone(summary["first_rank_with_all_orders_matched"])
        self.assertIn("note", summary)
        self.assertEqual(len(frame), 0)


if __name__ == "__main__":
    unittest.main()

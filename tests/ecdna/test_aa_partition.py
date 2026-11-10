"""Confusion-matrix metrics of rck.ecdna.aa_eval on small hand-built cases: M1 (species
detection) and M2 (pairwise co-membership of novel adjacencies). Definitions:
docs/pipeline/evaluation.md, section 4.7.

The cases: a merge, a split, a nested pair that shares a novel adjacency (blind and phased), a
background SV left unused (graph universe) and one routed into a walk, the 10 kb against 100 bp
tolerance, a linear walk, undefined rates, pooling, and the AA-file path end to end.

Run from the RCK root: python -m unittest tests.ecdna.test_aa_partition
"""
import math
import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna import aa_eval, aa_format
from rck.ecdna.evaluate import CStep, Cycle

# chr1 segments (0-based inclusive, the walk frame)
S1, S2, S3, S4, S5 = (0, 999), (1000, 1999), (2000, 2999), (3000, 3999), (5000, 5999)


def walk(*segs, hap="A", chrom="chr1", shift=0):
    return [CStep(chrom, s + shift, e + shift, hap, "+") for s, e in segs]


def cyc(name, weight, steps, circular=True):
    return Cycle(name, weight, [steps], circular=circular)


# Outer O = S1 S2 S3 S5: novel adjacencies j1 = S3|S5 and j2 = S5|S1 (wrap-around).
# Inner I = S2 S4: j3 = S2|S4 and j4 = S4|S2 (wrap-around). O and I share reference segments
# (S2) but no novel adjacency, as on the bench. Their Euler union is one walk with j1..j4.
OUTER = walk(S1, S2, S3, S5)
INNER = walk(S2, S4)
MERGED = walk(S1, S2, S4, S2, S3, S5)


def truth_pair():
    return [cyc("outer", 8.0, OUTER), cyc("inner", 4.0, INNER)]


def m2(truth, predicted, tol=10_000, phased=False, extra=()):
    items, ts, ps, _ = aa_eval.adjacency_items(truth, predicted, tol, phased, extra)
    counts = aa_eval.pair_confusion(len(items), ts, ps)
    counts.update(aa_eval.confusion_rates(counts["tp"], counts["fp"], counts["fn"], counts["tn"]))
    counts["n_items"] = len(items)
    return counts


def nan(x):
    return isinstance(x, float) and math.isnan(x)


class TestRates(unittest.TestCase):
    def test_undefined_is_nan_not_one(self):
        r = aa_eval.confusion_rates(0, 0, 0, 0)
        for k in ("precision", "recall", "fpr", "f1"):
            self.assertTrue(nan(r[k]), k)
        r = aa_eval.confusion_rates(3, None, 1)  # best-of: no false-positive count
        self.assertAlmostEqual(r["recall"], 0.75)
        self.assertTrue(nan(r["precision"]) and nan(r["f1"]) and nan(r["fpr"]))

    def test_f1_is_dice(self):
        r = aa_eval.confusion_rates(2, 4, 0, 0)
        self.assertAlmostEqual(r["f1"], 2 * 2 / (2 * 2 + 4 + 0))
        self.assertAlmostEqual(r["f1"], 2 * r["precision"] * r["recall"] / (r["precision"] + r["recall"]))


class TestPairs(unittest.TestCase):
    def test_perfect(self):
        c = m2(truth_pair(), truth_pair())
        self.assertEqual((c["n_items"], c["tp"], c["fp"], c["fn"], c["tn"]), (4, 2, 0, 0, 4))
        self.assertEqual((c["precision"], c["recall"], c["fpr"], c["f1"]), (1.0, 1.0, 0.0, 1.0))

    def test_merge_gives_false_positive_pairs(self):
        # one walk carries both molecules: the four cross pairs j1-j3, j1-j4, j2-j3, j2-j4 are FPs
        c = m2(truth_pair(), [cyc("m", 4.0, MERGED)])
        self.assertEqual((c["n_items"], c["tp"], c["fp"], c["fn"], c["tn"]), (4, 2, 4, 0, 0))
        self.assertAlmostEqual(c["precision"], 1 / 3)
        self.assertEqual(c["recall"], 1.0)
        self.assertEqual(c["fpr"], 1.0)  # every truth-negative pair merged: the merge rate
        self.assertAlmostEqual(c["f1"], 0.5)

    def test_split_gives_false_negative_pairs(self):
        # the truth is one molecule (the figure-8 walk); the answer splits it in two
        c = m2([cyc("t", 4.0, MERGED)], truth_pair())
        self.assertEqual((c["n_items"], c["tp"], c["fp"], c["fn"], c["tn"]), (4, 2, 0, 4, 0))
        self.assertAlmostEqual(c["recall"], 1 / 3)  # 1 - recall: the split rate
        self.assertEqual(c["precision"], 1.0)
        self.assertTrue(nan(c["fpr"]))  # no truth-negative pair at all

    def test_nested_shared_adjacency_is_one_item(self):
        # inner2 = S3 S5 shares j1 = S3|S5 with the outer molecule; its own adjacency is S5|S3
        inner2 = walk(S3, S5)
        truth = [cyc("outer", 8.0, OUTER), cyc("inner", 4.0, inner2)]
        c = m2(truth, truth)
        self.assertEqual((c["n_items"], c["tp"], c["fp"], c["fn"], c["tn"]), (3, 2, 0, 0, 1))
        merged = [cyc("m", 4.0, walk(S1, S2, S3, S5, S3, S5))]  # Euler union of the two
        c = m2(truth, merged)
        self.assertEqual((c["tp"], c["fp"], c["fn"], c["tn"]), (2, 1, 0, 0))
        # phased: the outer molecule on A, the inner on B; the shared adjacency is then two items
        truth_ab = [cyc("outer", 8.0, OUTER), cyc("inner", 4.0, walk(S3, S5, hap="B"))]
        c = m2(truth_ab, truth_ab, phased=True)
        self.assertEqual((c["n_items"], c["tp"], c["fp"], c["fn"], c["tn"]), (4, 2, 0, 0, 4))
        c = m2(truth_ab, truth_ab, phased=False)  # blind: one item again
        self.assertEqual(c["n_items"], 3)

    def test_wrong_haplotype_is_a_false_item_when_phased(self):
        truth = truth_pair()
        answer = [cyc("outer", 8.0, OUTER), cyc("inner", 4.0, walk(S2, S4, hap="B"))]
        blind = m2(truth, answer, phased=False)
        self.assertEqual(blind["f1"], 1.0)
        phased = m2(truth, answer, phased=True)
        self.assertEqual((phased["n_items"], phased["tp"], phased["fp"], phased["fn"]), (6, 1, 1, 1))

    def test_background_sv_left_unused_adds_true_negatives_only(self):
        background = [(("chr9", 100, "+"), ("chr9", 50_000, "-"), ".", ".")]
        native = m2(truth_pair(), truth_pair())
        graph = m2(truth_pair(), truth_pair(), extra=background)
        self.assertEqual(graph["n_items"], native["n_items"] + 1)
        for k in ("tp", "fp", "fn", "precision", "recall", "f1"):
            self.assertEqual(graph[k], native[k], k)
        self.assertEqual(graph["tn"], native["tn"] + native["n_items"])
        # a graph edge that is a truth adjacency (AA frame shifted by 1 bp) is not added twice
        j1 = (("chr1", 2999, "+"), ("chr1", 5000, "-"), ".", ".")
        self.assertEqual(m2(truth_pair(), truth_pair(), extra=[j1])["n_items"], native["n_items"])

    def test_background_sv_routed_into_a_walk_is_a_false_item(self):
        # the answer routes the outer molecule through chr9 (two novel adjacencies no truth has)
        detour = [CStep("chr1", 0, 999, "A", "+"), CStep("chr1", 1000, 1999, "A", "+"), CStep("chr1", 2000, 2999, "A", "+"),
                  CStep("chr9", 100, 999, "A", "+"), CStep("chr1", 5000, 5999, "A", "+")]
        c = m2(truth_pair(), [cyc("o", 8.0, detour), cyc("i", 4.0, INNER)])
        self.assertEqual(c["n_items"], 6)  # j1 missing; S3|chr9 and chr9|S5 false
        # the detour walk keeps j2: its pairs x1-x2, x1-j2 and x2-j2 are FPs; j1-j2 is an FN
        self.assertEqual((c["tp"], c["fp"], c["fn"]), (1, 3, 1))

    def test_tolerance(self):
        # S5 moved 5 kb: both of the outer molecule's adjacencies (S3|S5, S5|S1) are 5 kb off,
        # matched at 10 kb, two false items at 100 bp
        off = walk(S1, S2, S3) + [CStep("chr1", 10_000, 10_999, "A", "+")]
        answer = [cyc("o", 8.0, off), cyc("i", 4.0, INNER)]
        self.assertEqual(m2(truth_pair(), answer, tol=10_000)["f1"], 1.0)
        c = m2(truth_pair(), answer, tol=100)
        self.assertEqual((c["n_items"], c["tp"], c["fp"], c["fn"]), (6, 1, 1, 1))

    def test_junction_distance(self):
        a = (("chr1", 2999, "+"), ("chr1", 5000, "-"), "A", "A")
        near = (("chr1", 3049, "+"), ("chr1", 5000, "-"), "A", "A")
        self.assertEqual(aa_eval.junction_distance(near, a, True), 50)
        self.assertEqual(aa_eval.junction_distance(a, (a[1], a[0], "A", "A"), True), 0)  # either direction
        self.assertIsNone(aa_eval.junction_distance(a, (a[0], a[1], "B", "B"), True))
        self.assertEqual(aa_eval.junction_distance(a, (a[0], a[1], "B", "B"), False), 0)

    def test_one_to_one_matching(self):
        # two answer adjacencies within tolerance of one truth adjacency (50 bp and 1 kb from
        # S3|S5): only the nearer one takes it, the other is a false item. 1.5 kb, because at
        # 10 kb these toy coordinates would also reach the inner molecule's adjacencies
        o_near = walk(S1, S2, S3) + [CStep("chr1", 5050, 5999, "A", "+")]
        o_far = walk(S1, S2, S3) + [CStep("chr1", 6000, 6999, "A", "+")]
        items, ts, ps, n_truth = aa_eval.adjacency_items(truth_pair(), [cyc("a", 1, o_near), cyc("b", 1, o_far)], 1_500, False)
        self.assertEqual(n_truth, 4)
        self.assertIn(0, ps[0])  # j1 (S3|S5) went to the nearer adjacency
        self.assertNotIn(0, ps[1])
        self.assertGreaterEqual(max(ps[1]), n_truth)

    def test_linear_walk_has_no_wrap_around_adjacency(self):
        c = m2(truth_pair(), [cyc("o", 8.0, OUTER, circular=False), cyc("i", 4.0, INNER)])
        self.assertEqual((c["tp"], c["fp"], c["fn"]), (1, 0, 1))  # j2 (wrap-around) is missing


def species_frame(mode, rows):
    return pd.DataFrame([dict(mode=mode, **r) for r in rows])


class TestDetection(unittest.TestCase):
    def frame(self):
        strict = [dict(exact_cycle=1, composite_h0=0.99, exact_cycle_phased=1, composite=0.98),
                  dict(exact_cycle=0, composite_h0=0.70, exact_cycle_phased=0, composite=0.60)]
        best = [dict(exact_cycle=1, composite_h0=0.99, exact_cycle_phased=1, composite=0.98),
                dict(exact_cycle=1, composite_h0=0.97, exact_cycle_phased=1, composite=0.96)]
        return pd.concat([species_frame("strict", strict), species_frame("best_of", best)], ignore_index=True)

    def test_matched_but_failing_partner_is_a_false_positive(self):
        rows = pd.DataFrame(aa_eval.detection_rows(self.frame(), n_truth=2, n_strict=3, phased=True))
        r = rows.set_index(["mode", "haplotype", "criterion"])
        e = r.loc[("strict", "blind", "exact")]
        self.assertEqual((e["tp"], e["fp"], e["fn"]), (1, 2, 1))  # the failing partner and the unmatched one
        self.assertAlmostEqual(e["precision"], 1 / 3)
        self.assertAlmostEqual(e["recall"], 0.5)
        self.assertTrue(pd.isna(e["tn"]) and nan(e["fpr"]))
        self.assertTrue(e["headline"])
        self.assertEqual(e["tol_bp"], aa_eval.GRAPH_TOL)
        c = r.loc[("strict", "blind", "composite_h0>=0.6")]
        self.assertEqual((c["tp"], c["fp"], c["fn"], c["tol_bp"]), (2, 1, 0, aa_eval.SPECIES_TOL))
        self.assertEqual(r.loc[("strict", "phased", "composite>=0.6")]["tp"], 2)
        self.assertEqual(r.loc[("strict", "phased", "composite>=0.7")]["tp"], 1)

    def test_curve_is_monotone_and_best_of_has_recall_only(self):
        rows = pd.DataFrame(aa_eval.detection_rows(self.frame(), n_truth=2, n_strict=3, phased=True))
        curve = rows[(rows["mode"] == "strict") & (rows["haplotype"] == "blind") & (rows["criterion"] != "exact")]
        self.assertEqual(list(curve["criterion"]), ["composite_h0>=%g" % t for t in aa_eval.DETECTION_TAUS])
        self.assertTrue((curve["tp"].diff().dropna() <= 0).all())
        best = rows[(rows["mode"] == "best_of") & (rows["criterion"] == "exact")].iloc[0]
        self.assertEqual(best["tp"], 2)
        self.assertTrue(pd.isna(best["fp"]) and nan(best["precision"]) and nan(best["f1"]))
        self.assertEqual(best["recall"], 1.0)

    def test_unphased_answer_has_blind_rows_only(self):
        rows = pd.DataFrame(aa_eval.detection_rows(self.frame(), n_truth=2, n_strict=2, phased=False))
        self.assertEqual(set(rows["haplotype"]), {"blind"})


class TestPooling(unittest.TestCase):
    def test_micro_and_macro(self):
        base = dict(metric="adjacency_pairs", mode="strict", haplotype="blind", criterion="co_membership", tol_bp=10_000,
                    universe="native", headline=True)
        big = dict(base, tp=90, fp=0, fn=0, tn=10, **aa_eval.confusion_rates(90, 0, 0, 10))
        small = dict(base, tp=1, fp=1, fn=1, tn=0, **aa_eval.confusion_rates(1, 1, 1, 0))
        pooled = aa_eval.pool_partition(pd.DataFrame([big, small])).iloc[0]
        self.assertEqual(pooled["headline_pooling"], "macro")
        self.assertEqual(pooled["n_cases"], 2)
        self.assertAlmostEqual(pooled["micro_f1"], 2 * 91 / (2 * 91 + 1 + 1))
        self.assertAlmostEqual(pooled["macro_f1"], (1.0 + 0.5) / 2)
        self.assertAlmostEqual(pooled["macro_fpr"], (0.0 + 1.0) / 2)
        self.assertEqual(pooled["n_cases_fpr"], 2)

    def test_undefined_rates_are_left_out_of_the_macro_mean(self):
        base = dict(metric="species_detection", mode="best_of", haplotype="blind", criterion="exact", tol_bp=100,
                    universe="species", headline=False, fp=None, tn=None)
        rows = [dict(base, tp=1, fn=1, **aa_eval.confusion_rates(1, None, 1)), dict(base, tp=2, fn=0, **aa_eval.confusion_rates(2, None, 0))]
        pooled = aa_eval.pool_partition(pd.DataFrame(rows)).iloc[0]
        self.assertEqual(pooled["headline_pooling"], "micro")
        self.assertAlmostEqual(pooled["micro_recall"], 0.75)
        self.assertAlmostEqual(pooled["macro_recall"], 0.75)
        self.assertTrue(pooled["fp"] is None and nan(pooled["micro_precision"]))
        self.assertEqual(pooled["n_cases_precision"], 0)


def cycles_text(intervals, segments, cycles):
    lines = ["Interval\t{i}\t{c}\t{s}\t{e}".format(i=i, c=c, s=s, e=e) for i, (c, s, e) in enumerate(intervals, start=1)]
    lines.append("List of cycle segments")
    lines += ["Segment\t{i}\t{c}\t{s}\t{e}".format(i=i, c=c, s=s, e=e) for i, (c, s, e) in enumerate(segments, start=1)]
    lines += ["Cycle={k};".format(k=k) + c for k, c in enumerate(cycles, start=1)]
    return "\n".join(lines) + "\n"


class TestEndToEnd(unittest.TestCase):
    """The AA-file path: a merged answer of the nested pair, phased, with a graph that also holds
    an unused background SV. score() writes aa_partition_metrics.tsv; the existing tables keep
    their columns."""

    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def test_score_writes_partition_metrics(self):
        # AA segments (1-based): 1 = S1, 2 = S2, 3 = S3, 4 = S4, 5 = S5
        segs = [("chr1", s + 1, e + 1) for s, e in (S1, S2, S3, S4, S5)]
        merged = "Copy_count=4.0;Segments=1+,2+,4+,2+,3+,5+;Haplotypes=A,A,A,A,A,A;Status=species;Species=m"
        cycles = os.path.join(self.dir, "a_cycles.txt")
        with open(cycles, "w") as handle:
            handle.write(cycles_text([("chr1", 1, 6000)], segs, [merged]))
        graph = os.path.join(self.dir, "a_graph.txt")
        edges = ["chr1:3000+->chr1:5001-", "chr1:6000+->chr1:1-", "chr1:2000+->chr1:3001-", "chr1:4000+->chr1:1001-",
                 "chr9:101+->chr9:50001-"]
        with open(graph, "w") as handle:
            handle.write("SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped\n")
            for c, s, e in segs:
                handle.write("sequence\t{c}:{s}-\t{c}:{e}+\t4\t4\t{n}\t0\n".format(c=c, s=s, e=e, n=e - s + 1))
            handle.write("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, x, y\n")
            for e in edges:
                handle.write("discordant\t{e}\t4\t0\tNone\tNone\n".format(e=e))
        result = aa_eval.score(cycles, truth_pair(), graph)
        p = result.partition
        self.assertEqual(list(p.columns), aa_eval.PARTITION_COLUMNS)
        det = p[(p["metric"] == "species_detection") & (p["mode"] == "strict") & (p["criterion"] == "exact")].set_index("haplotype")
        for hap in ("blind", "phased"):
            self.assertEqual((det.loc[hap, "tp"], det.loc[hap, "fp"], det.loc[hap, "fn"]), (0, 1, 2))
        pairs = p[(p["metric"] == "adjacency_pairs") & p["headline"].astype(bool)]
        self.assertEqual(set(pairs["haplotype"]), {"blind", "phased"})
        for r in pairs.itertuples():
            self.assertEqual((r.tp, r.fp, r.fn, r.tn), (2, 4, 0, 0))
        g = p[(p["metric"] == "adjacency_pairs") & (p["universe"] == "graph") & (p["tol_bp"] == aa_eval.SPECIES_TOL)].iloc[0]
        self.assertEqual((g["n_items"], g["tp"], g["fp"], g["tn"]), (5, 2, 4, 4))  # the chr9 SV: 4 more TN pairs
        self.assertFalse(((p["universe"] == "graph") & (p["haplotype"] == "phased")).any())
        out = os.path.join(self.dir, "out")
        paths = aa_eval.write(out, result)
        self.assertTrue(os.path.exists(paths["partition"]))
        back = pd.read_csv(paths["partition"], sep="\t")
        self.assertEqual(len(back), len(p))
        self.assertNotIn("tp", pd.read_csv(paths["sample"], sep="\t").columns)  # existing tables unchanged
        self.assertNotIn("tp", pd.read_csv(paths["species"], sep="\t").columns)


if __name__ == "__main__":
    unittest.main()

"""Method-agnostic scorer (rck.ecdna.aa_eval, aa_format, aa_export) on small hand-built cases.

Every metric has a perfect match, a missed edge, a false-positive segment (which fails
error-free), a rotated and reversed cycle for the cyclic LCS and a nested pair for best-of
against strict, as far as the metric applies.

Run from the RCK root: python -m unittest tests.ecdna.test_aa_eval
"""
import math
import os
import tempfile
import unittest

import pandas as pd

from rck.ecdna import aa_eval, aa_export, aa_format
from rck.ecdna import evaluate as ev
from rck.ecdna.evaluate import CStep, Cycle


def graph_text(sequences, discordant, concordant=()):
    lines = ["SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped"]
    for c, s, e, cn in sequences:
        lines.append("sequence\t{c}:{s}-\t{c}:{e}+\t{cn}\t{cn}\t{n}\t0".format(c=c, s=s, e=e, cn=cn, n=e - s + 1))
    lines.append("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, x, y")
    for a, b, cn in discordant:
        lines.append("discordant\t{a}->{b}\t{cn}\t0\tNone\tNone".format(a=a, b=b, cn=cn))
    for a, b, cn in concordant:
        lines.append("concordant\t{a}->{b}\t{cn}\t0\tNone\tNone".format(a=a, b=b, cn=cn))
    return "\n".join(lines) + "\n"


def cycles_text(intervals, segments, cycles):
    lines = ["Interval\t{i}\t{c}\t{s}\t{e}".format(i=i, c=c, s=s, e=e) for i, (c, s, e) in enumerate(intervals, start=1)]
    lines.append("List of cycle segments")
    lines += ["Segment\t{i}\t{c}\t{s}\t{e}".format(i=i, c=c, s=s, e=e) for i, (c, s, e) in enumerate(segments, start=1)]
    lines += ["Cycle={k};".format(k=k) + c for k, c in enumerate(cycles, start=1)]
    return "\n".join(lines) + "\n"


# a three-segment amplicon on chr1 (AA coordinates): 1-1000, 1001-2000, 2001-3000, circle 1+,3-,2+
SEGS = [("chr1", 1, 1000), ("chr1", 1001, 2000), ("chr1", 2001, 3000)]
EDGES = [("chr1:1000+", "chr1:3000+", 10.0), ("chr1:2001-", "chr1:1001-", 10.0), ("chr1:2000+", "chr1:1-", 10.0)]


class Tmp(unittest.TestCase):
    def setUp(self):
        self._tmp = tempfile.TemporaryDirectory()
        self.dir = self._tmp.name

    def tearDown(self):
        self._tmp.cleanup()

    def write(self, name, text):
        path = os.path.join(self.dir, name)
        with open(path, "w") as handle:
            handle.write(text)
        return path


class TestGraphMetrics(Tmp):
    def amplicon(self):
        truth = aa_format.read_graph(self.write("truth_graph.txt", graph_text([s + (12.0,) for s in SEGS], EDGES)))
        segs = pd.DataFrame([{"chrom": c, "start": s, "end": e, "cn": 12.0, "cn_a": 2.0, "cn_b": 10.0} for c, s, e in SEGS])
        return aa_eval.MixtureAmplicon(1, ["ecdna1"], [("chr1", 1, 3000)], truth, segs)

    def score(self, sequences, edges, cn=None, tol=100, concordant=()):
        g = aa_format.read_graph(self.write("inf_graph.txt", graph_text(sequences, edges, concordant)))
        return aa_eval.graph_metrics(self.amplicon(), g, cn, tol)

    def test_perfect_match(self):
        # breakpoints shifted by 1 and 99 bp: within the 100 bp tolerance
        edges = [("chr1:1001+", "chr1:3000+", 9.0), ("chr1:1002-", "chr1:2001-", 9.0), ("chr1:1999+", "chr1:100-", 9.0)]
        seqs = [("chr1", 100, 1000, 12.0), ("chr1", 1001, 2000, 12.0), ("chr1", 2001, 3000, 12.0)]
        r = self.score(seqs, edges)
        self.assertEqual(r["error_free"], 1)
        self.assertEqual(r["error_free_merged"], 1)
        self.assertEqual(r["bp_accuracy"], 1.0)
        self.assertEqual(r["novel_adjacency_fnr"], 0.0)
        self.assertEqual(r["novel_adjacency_fdr"], 0.0)
        self.assertEqual(r["is_violations"], 0)
        self.assertEqual(r["cn_error_mode"], "total")
        self.assertAlmostEqual(r["segment_cn_mae"], 99 * 12 / 3000)  # bases 1-99 have no segment

    def test_one_missed_edge(self):
        r = self.score([s + (12.0,) for s in SEGS], EDGES[:2])
        self.assertEqual(r["error_free"], 0)
        self.assertAlmostEqual(r["bp_accuracy"], 2 / 3)
        self.assertAlmostEqual(r["novel_adjacency_fnr"], 1 / 3)
        self.assertEqual(r["novel_adjacency_fdr"], 0.0)

    def test_tolerance_and_orientation(self):
        off = [("chr1:1101+", "chr1:3000+", 1.0)] + EDGES[1:]  # 101 bp off
        self.assertAlmostEqual(self.score([s + (12.0,) for s in SEGS], off)["bp_accuracy"], 2 / 3)
        self.assertEqual(self.score([s + (12.0,) for s in SEGS], off, tol=101)["bp_accuracy"], 1.0)
        flipped = [("chr1:1000-", "chr1:3000+", 1.0)] + EDGES[1:]
        r = self.score([s + (12.0,) for s in SEGS], flipped)
        self.assertAlmostEqual(r["bp_accuracy"], 2 / 3)
        self.assertAlmostEqual(r["novel_adjacency_fdr"], 1 / 3)

    def test_false_positive_segment_fails_error_free(self):
        # every edge right, one true segment split on a 10 kb-style grid point
        seqs = [("chr1", 1, 1000, 12.0), ("chr1", 1001, 1500, 12.0), ("chr1", 1501, 2000, 12.0), ("chr1", 2001, 3000, 12.0)]
        r = self.score(seqs, EDGES)
        self.assertEqual(r["bp_accuracy"], 1.0)
        self.assertEqual(r["n_false_segments"], 2)
        self.assertEqual(r["error_free"], 0)
        self.assertEqual(r["error_free_merged"], 0)  # no concordant edge at 1500: not a benign split

    def test_benign_copy_number_split_passes_error_free_merged(self):
        # 1001-2000 split at 1500 by a copy-number change (12 -> 14) with a concordant edge only:
        # the strict column counts two false segments, the merged column merges them back
        seqs = [("chr1", 1, 1000, 12.0), ("chr1", 1001, 1500, 12.0), ("chr1", 1501, 2000, 14.0), ("chr1", 2001, 3000, 12.0)]
        conc = [("chr1:1500+", "chr1:1501-", 12.0)]
        r = self.score(seqs, EDGES, concordant=conc)
        self.assertEqual(r["error_free"], 0)
        self.assertEqual(r["n_false_segments"], 2)
        self.assertEqual(r["error_free_merged"], 1)
        self.assertEqual(r["n_false_segments_merged"], 0)
        self.assertEqual(r["n_inferred_segments_merged"], 3)
        # the same split with a discordant edge at the boundary is a breakpoint, not merged
        r = self.score(seqs, EDGES + [("chr1:1500+", "chr1:2500-", 1.0)], concordant=conc)
        self.assertEqual(r["error_free_merged"], 0)
        self.assertEqual(r["n_false_segments_merged"], 2)

    def test_merge_benign_splits(self):
        # a concordant edge next to a discordant one (the truth's 1000/1001 boundary) is kept split;
        # a gap between segments is never merged
        g = aa_format.read_graph(self.write("g.txt", graph_text(
            [("chr1", 1, 1000, 5.0), ("chr1", 1001, 2000, 7.0), ("chr1", 2001, 3000, 9.0), ("chr1", 3101, 4000, 9.0)],
            [("chr1:1000+", "chr1:3000+", 1.0)],
            [("chr1:1000+", "chr1:1001-", 4.0), ("chr1:2001-", "chr1:2000+", 6.0), ("chr1:3000+", "chr1:3101-", 1.0)])))
        merged = aa_eval.merge_benign_splits(g)
        self.assertEqual([(m.start, m.end) for m in merged], [(1, 1000), (1001, 3000), (3101, 4000)])
        self.assertAlmostEqual(merged[1].cn, 8.0)

    def test_flank_segment_is_not_counted(self):
        seqs = [("chr1", 1, 1000, 12.0), ("chr1", 1001, 2000, 12.0), ("chr1", 2001, 3000, 12.0), ("chr1", 2960, 9000, 2.0)]
        self.assertEqual(self.score(seqs, EDGES)["error_free"], 1)  # overlaps the amplicon by only 41 bp

    def test_false_positive_edge(self):
        r = self.score([s + (12.0,) for s in SEGS], EDGES + [("chr1:1500+", "chr1:2500-", 1.0)])
        self.assertEqual(r["bp_accuracy"], 1.0)
        self.assertEqual(r["error_free"], 0)
        self.assertAlmostEqual(r["novel_adjacency_fdr"], 1 / 4)

    def test_infinite_sites(self):
        # chr1:1000+ used by two novel adjacencies (one 50 bp away): one violation
        r = self.score([s + (12.0,) for s in SEGS], EDGES + [("chr1:1050+", "chr1:2500-", 1.0)])
        self.assertEqual(r["is_violations"], 1)
        self.assertEqual(r["is_violations_truth"], 0)
        e = [aa_format.AAEdge("discordant", ("c", 10, "+"), ("c", 500, "-"), 1), aa_format.AAEdge("discordant", ("c", 10, "+"), ("c", 900, "-"), 1),
             aa_format.AAEdge("discordant", ("c", 500, "-"), ("c", 10, "+"), 1)]  # the third repeats the first
        self.assertEqual(aa_eval.is_violations(e, 100), 1)
        self.assertEqual(aa_eval.is_violations(e[:1], 100), 0)

    def test_allele_copy_number_error(self):
        cn = pd.DataFrame([{"chrom": c, "start": s, "end": e, "cn": 12.0, "cn_a": 10.0, "cn_b": 2.0} for c, s, e in SEGS])  # A/B swapped
        r = self.score([s + (12.0,) for s in SEGS], EDGES, cn)
        self.assertEqual(r["cn_error_mode"], "allele")
        self.assertAlmostEqual(r["segment_cn_mae"], 16.0)
        self.assertAlmostEqual(r["segment_cn_mae_total"], 0.0)
        self.assertAlmostEqual(r["segment_cn_rel_err"], 16 / 12)


def steps(*items, hap="A"):
    return [CStep("chr1", s, e, hap, st) for s, e, st in items]


A, B, C, D = (0, 999), (1000, 1999), (2000, 2999), (5000, 5999)


class TestCycleMetrics(unittest.TestCase):
    truth = steps(A + ("+",), C + ("-",), B + ("+",))

    def test_perfect_match(self):
        self.assertEqual(aa_eval.cycle_overlap(self.truth, self.truth), 1.0)
        self.assertEqual(aa_eval.cycle_overlap_coral_code(self.truth, self.truth), 1.0)
        self.assertEqual(aa_eval.cyclic_lcs(self.truth, self.truth)[:2], (1.0, 1.0))
        self.assertEqual(aa_eval.length_error(self.truth, self.truth), 0.0)

    def test_rotated_and_reversed(self):
        # the same circle written from another start, read backwards: B-, C+, A-
        inferred = steps(B + ("-",), C + ("+",), A + ("-",))
        self.assertEqual(aa_eval.cyclic_lcs(self.truth, inferred)[0], 1.0)
        self.assertEqual(aa_eval.cyclic_lcs(self.truth, inferred[1:] + inferred[:1])[0], 1.0)
        # a different order of the same intervals: A+, B+, C- loses one interval in the LCS
        other = steps(A + ("+",), B + ("+",), C + ("-",))
        self.assertAlmostEqual(aa_eval.cyclic_lcs(self.truth, other)[0], 2 / 3)
        self.assertEqual(aa_eval.cycle_overlap(self.truth, other), 1.0)  # same intervals, no order

    def test_missed_interval(self):
        inferred = steps(A + ("+",), C + ("-",))
        self.assertAlmostEqual(aa_eval.cycle_overlap(self.truth, inferred), 2 / 3)
        lcs, lcs_code, raw = aa_eval.cyclic_lcs(self.truth, inferred)
        self.assertAlmostEqual(lcs, 2 / 3)
        self.assertEqual(raw, 2000)
        self.assertAlmostEqual(aa_eval.length_error(self.truth, inferred), math.log2(2 / 3))  # log2(1 + R), R = -1/3

    def test_extra_interval_and_multiplicity(self):
        inferred = steps(A + ("+",), C + ("-",), B + ("+",), D + ("+",))
        self.assertAlmostEqual(aa_eval.cycle_overlap(self.truth, inferred), 3 / 4)
        self.assertAlmostEqual(aa_eval.cyclic_lcs(self.truth, inferred)[0], 1.0)  # D is in one cycle only
        self.assertAlmostEqual(aa_eval.cyclic_lcs(self.truth, inferred)[1], 3 / 4)  # CoRAL's code: over 4 kb
        dup = steps(A + ("+",), A + ("+",), C + ("-",), B + ("+",))
        self.assertEqual(aa_eval.cycle_overlap(self.truth, dup), 1.0)  # no multiplicity
        self.assertAlmostEqual(aa_eval.cycle_overlap_coral_code(self.truth, dup), 3 / 4)
        self.assertAlmostEqual(aa_eval.length_error(self.truth, dup), math.log2(4 / 3))

    def test_other_chromosome_never_overlaps(self):
        other = [CStep("chr2", s.start, s.end, s.hap, s.strand) for s in self.truth]
        self.assertEqual(aa_eval.cycle_overlap(self.truth, other), 0.0)
        self.assertEqual(aa_eval.cyclic_lcs(self.truth, other)[0], 0.0)

    def test_exact_cycle(self):
        t = Cycle("t", 10, [self.truth])
        near = Cycle("c", 10, [[CStep("chr1", s.start + 50, s.end + 50, "A", s.strand) for s in self.truth]])
        self.assertEqual(aa_eval._exact(t, near, 100), 1)
        self.assertEqual(aa_eval._exact(t, near, 10), 0)
        self.assertEqual(aa_eval._exact(t, Cycle("c", 10, [steps(A + ("+",), B + ("+",), C + ("-",))]), 100), 0)


def nested_truth():
    """An outer circle A+ B+ C+ (8 copies) and a nested circle B+ (4 copies) on the same haplotype."""
    return [Cycle("outer", 8.0, [steps(A + ("+",), (1000, 1499, "+"), (1500, 1999, "+"), C + ("+",))]),
            Cycle("nested", 4.0, [steps((1000, 1499, "+"), (1500, 1999, "+"))])]


class TestModes(Tmp):
    """A nested pair: the method's answer merges the two molecules into one walk; its
    alternative has them apart. Strict scores the answer; best-of finds the alternative."""

    def files(self, phased=True):
        segs = [("chr1", 1, 1000), ("chr1", 1001, 1500), ("chr1", 1501, 2000), ("chr1", 2001, 3000)]
        h = lambda n: (";Haplotypes=" + ",".join(["A"] * n)) if phased else ""  # noqa: E731
        merged = ["Copy_count=4.0;Segments=1+,2+,3+,2+,3+,4+" + h(6) + ";Status=species;Species=ec1",
                  "Copy_count=4.0;Segments=1+,2+,3+,4+" + h(4) + ";Status=species;Species=ec2",
                  "Copy_count=0.5;Segments=2+,3+" + h(2) + ";Status=low_weight;Species=ec3"]
        apart = ["Copy_count=8.0;Segments=1+,2+,3+,4+" + h(4) + ";Species=ec1",
                 "Copy_count=4.0;Segments=2+,3+" + h(2) + ";Species=ec2"]
        iv = [("chr1", 1, 3000)]
        return self.write("a_cycles.txt", cycles_text(iv, segs, merged)), self.write("alt_cycles.txt", cycles_text(iv, segs, apart))

    def test_strict_against_best_of(self):
        main, alt = self.files()
        result = aa_eval.score(main, nested_truth(), alt_cycles=[alt])
        strict = result.species[result.species["mode"] == "strict"].set_index("truth")
        best = result.species[result.species["mode"] == "best_of"].set_index("truth")
        self.assertEqual(len(result.strict), 2)  # ec3 is low_weight
        self.assertLess(strict.loc["nested", "composite"], 0.9)
        self.assertGreater(best.loc["nested", "composite"], 0.99)
        self.assertEqual(best.loc["nested", "inferred"], "alt1.ec2")
        self.assertEqual(best.loc["outer", "exact_cycle"], 1)
        self.assertGreaterEqual(best.loc["outer", "composite"], strict.loc["outer", "composite"])
        sample = result.sample.set_index("mode")
        self.assertEqual(sample.loc["strict", "species_count_error"], 0)
        self.assertEqual(sample.loc["strict", "n_false"], 0)
        self.assertEqual(sample.loc["strict", "n_missed"], 0)
        self.assertGreater(sample.loc["best_of", "mean_composite"], sample.loc["strict", "mean_composite"])

    def test_coral_code_columns_are_the_headline(self):
        # CoRAL's code versions come first; the supplementary-methods versions are *_supp. The
        # strict answer's merged walk contains the nested species: the supplement's LCS (over the
        # true length) forgives it, CoRAL's code (over the joint bins) does not
        main, alt = self.files()
        result = aa_eval.score(main, nested_truth(), alt_cycles=[alt])
        cols = list(result.species.columns)
        self.assertLess(cols.index("cycle_overlap_coral_code"), cols.index("cycle_overlap_supp"))
        self.assertLess(cols.index("cyclic_lcs_coral_code"), cols.index("cyclic_lcs_supp"))
        self.assertNotIn("cyclic_lcs", cols)
        self.assertNotIn("cycle_overlap", cols)
        best = result.species[result.species["mode"] == "best_of"].set_index("truth")
        self.assertLess(best.loc["nested", "cyclic_lcs_coral_code"], best.loc["nested", "cyclic_lcs_supp"] + 1e-12)
        sample = result.sample.set_index("mode")
        for c in ("mean_cycle_overlap_coral_code", "mean_cyclic_lcs_coral_code", "mean_cycle_overlap_supp", "mean_cyclic_lcs_supp"):
            self.assertIn(c, sample.columns)
        self.assertNotIn("mean_cyclic_lcs", sample.columns)

    def test_false_and_missed_species(self):
        main, _ = self.files()
        result = aa_eval.score(main, nested_truth()[:1], strict_rule="all")  # three cycles for one true species
        s = result.sample.set_index("mode").loc["strict"]
        self.assertEqual(s["species_count_error"], 2)
        self.assertEqual(s["n_false"], 2)
        far = [Cycle("far", 5.0, [[CStep("chr9", 0, 999, "A", "+")]])]
        s = aa_eval.score(main, nested_truth() + far, all_truth_species=True).sample.set_index("mode")
        self.assertEqual(s.loc["strict", "n_missed"], 1)
        self.assertEqual(s.loc["best_of", "n_missed"], 1)

    def test_unphased_method(self):
        main, alt = self.files(phased=False)
        result = aa_eval.score(main, nested_truth(), alt_cycles=[alt], strict_rule="cn-frac:0.2")
        self.assertFalse(result.info["phased"])
        strict = result.species[result.species["mode"] == "strict"].set_index("truth")
        self.assertTrue(strict["composite"].isna().all())  # phasing accuracy is undefined
        self.assertTrue(strict["phasing_accuracy"].isna().all())
        self.assertGreater(strict.loc["outer", "composite_h0"], 0.5)
        best = result.species[result.species["mode"] == "best_of"].set_index("truth")
        self.assertGreater(best.loc["nested", "composite_h0"], 0.99)

    def test_composite_h0_ignores_haplotypes(self):
        """A phased answer with every haplotype wrong keeps its composite_h0."""
        segs = [("chr1", 1, 1000), ("chr1", 1001, 2000), ("chr1", 2001, 3000)]
        right = self.write("r_cycles.txt", cycles_text([("chr1", 1, 3000)], segs, ["Copy_count=8.0;Segments=1+,2+,3+;Haplotypes=A,A,A"]))
        wrong = self.write("w_cycles.txt", cycles_text([("chr1", 1, 3000)], segs, ["Copy_count=8.0;Segments=1+,2+,3+;Haplotypes=B,B,B"]))
        truth = [Cycle("t", 8.0, [steps(A + ("+",), B + ("+",), C + ("+",))])]
        r = aa_eval.score(right, truth).species.iloc[0]
        w = aa_eval.score(wrong, truth).species.iloc[0]
        self.assertEqual(r["composite"], 1.0)
        self.assertEqual(r["composite_h0"], 1.0)
        self.assertEqual(w["composite_h0"], 1.0)
        self.assertTrue(pd.isna(w["inferred"]))  # the phased matching finds no overlap on haplotype A

    def test_strict_rules(self):
        cycles = [Cycle("a", 10.0, [steps(A + ("+",))]), Cycle("b", 2.0, [steps(B + ("+",))]), Cycle("c", 0.5, [steps(C + ("+",))])]
        records = aa_format.AACycles("x", [], {}, [])
        names = lambda rule: [c.name for c in aa_eval.strict_select(records, cycles, rule)[0]]  # noqa: E731
        self.assertEqual(names("auto"), ["a", "b"])  # no Status: cn-frac:0.1
        self.assertEqual(names("cn-frac:0.3"), ["a"])
        self.assertEqual(names("top-cn:0.8"), ["a"])
        self.assertEqual(names("top-cn:0.9"), ["a", "b"])
        self.assertEqual(names("all"), ["a", "b", "c"])


class TestSampleMetrics(Tmp):
    def test_cn_explained(self):
        g = aa_format.read_graph(self.write("g.txt", graph_text([("chr1", 1, 1000, 10.0), ("chr1", 1001, 2000, 6.0),
                                                                 ("chr1", 2001, 3000, 2.0)], [])))
        cycles = [Cycle("a", 4.0, [steps(A + ("+",), B + ("+",))]), Cycle("b", 2.0, [steps(A + ("+",))]),
                  Cycle("p", 9.0, [steps(B + ("+",))], circular=False)]
        r = aa_eval.cn_explained(g, cycles)  # total 16,000 (the 2-copy segment is below 5)
        self.assertAlmostEqual(r["cn_explained_top1"], 8000 / 16000)
        self.assertAlmostEqual(r["cn_explained_top3"], 10000 / 16000)  # the path is left out

    def test_subwalks(self):
        # a read through A+ -> C- -> B+ (two novel junctions), then one that is not in the cycle
        paths = self.write("paths.tsv", "PID\tPATH\tNSTEP\t#NOV\n"
                           "p1\tchr1:0-@-,chr1:999+@+,chr1:2999+@+,chr1:2000-@-,chr1:1000-@-,chr1:1999+@+\t6\t3\n"
                           "p2\tchr1:0-@-,chr1:999+@+,chr1:1000-@-,chr1:1999+@+,chr1:2000-@-,chr1:2999+@+\t6\t1\n"
                           "p3\tchr1:1000-@-,chr1:1999+@+,chr1:0-@-,chr1:999+@+,chr1:2999+@+,chr1:2000-@-\t6\t1\n")
        reads = aa_eval.read_paths(paths)
        self.assertEqual([p.pid for p in reads], ["p1", "p3"])  # p2 is a reference walk
        cycle = Cycle("c", 1.0, [steps(A + ("+",), C + ("-",), B + ("+",))])
        r = aa_eval.subwalks_satisfied(reads, [cycle], 100)
        self.assertEqual(r["subwalks_satisfied"], 1.0)  # p3 is the same order across the closing join
        reverse = Cycle("r", 1.0, [steps(B + ("-",), C + ("+",), A + ("-",))])
        self.assertEqual(aa_eval.subwalks_satisfied(reads, [reverse], 100)["subwalks_satisfied"], 1.0)
        other = Cycle("o", 1.0, [steps(A + ("+",), B + ("-",), C + ("+",))])
        self.assertEqual(aa_eval.subwalks_satisfied(reads, [other], 100)["subwalks_satisfied"], 0.0)


class TestFormat(Tmp):
    def test_cycles_fields(self):
        segs = [("chr1", 1, 1000), ("chr1", 1001, 2000)]
        path = self.write("c.txt", cycles_text([("chr1", 1, 2000)], segs, [
            "Copy_count=3.5;Segments=1+,2-;Haplotypes=A,B;Alt_segments=2-,1+;Alt_haplotypes=B,A",
            "Copy_count=1.0;Segments=0+,1+,2+,0-;Haplotypes=.,A,A,."]))
        records = aa_format.read_cycles(path)
        self.assertTrue(records.phased)
        cycles = aa_format.to_cycles(records)
        self.assertEqual(cycles[0].steps, [CStep("chr1", 0, 999, "A", "+"), CStep("chr1", 1000, 1999, "B", "-")])
        self.assertEqual(len(cycles[0].orderings), 2)
        self.assertFalse(cycles[1].circular)
        self.assertEqual([s.hap for s in cycles[1].steps], ["A", "A"])
        unphased = aa_format.read_cycles(self.write("u.txt", cycles_text([], segs, ["Copy_count=2;Segments=1+,2+"])))
        self.assertFalse(unphased.phased)
        self.assertEqual(aa_format.to_cycles(unphased)[0].steps[0].hap, aa_format.UNPHASED)

    def test_truth_reader_unchanged(self):
        """evaluate.read_aa_cycles, now on aa_format, still gives the simulator's rows and walk."""
        path = self.write("t.txt", cycles_text([], [("chr17", 3906001, 4925000), ("chr17", 4925001, 5783000)],
                                               ["Copy_count=1.0;Length=1877000;Segments=2+,1-,2+"]))
        rows, walk = ev.read_aa_cycles(path)[0]
        self.assertEqual(rows, [{"#CHR": "chr17", "START": 4925000, "END": 5783000}, {"#CHR": "chr17", "START": 3906000, "END": 4925000}])
        self.assertEqual(walk, [(1, "+"), (2, "-"), (1, "+")])
        bad = self.write("b.txt", cycles_text([], [("chr1", 1, 10)], ["Copy_count=1.0;Segments=0+,1+,0-"]))
        with self.assertRaises(ValueError):
            ev.read_aa_cycles(bad)


class TestExport(Tmp):
    def test_round_trip_reproduces_rck_ecdna_eval(self):
        """The exported cycles score exactly as rck-ecdna-eval scores decompositions.tsv."""
        run = os.path.join(self.dir, "run")
        os.makedirs(os.path.join(run, "aa"))
        walk1 = "chr1:1-1000:A:+,chr1:2001-3000:A:-,chr1:1001-2000:A:+"
        walk1b = "chr1:1-1000:A:+,chr1:1001-2000:A:+,chr1:2001-3000:A:-"
        rows = [dict(decomposition_rank=1, objective=1.0, species="ec1", topology="circular", status="species", weight=7.3,
                     weight_int=7, ordering_rank=1, walk=walk1),
                dict(decomposition_rank=1, objective=1.0, species="ec1", topology="circular", status="species", weight=7.3,
                     weight_int=7, ordering_rank=2, walk=walk1b),
                dict(decomposition_rank=1, objective=1.0, species="ec2", topology="circular", status="low_weight", weight=0.4,
                     weight_int=1, ordering_rank=1, walk="chr1:1001-2000:A:+")]
        pd.DataFrame(rows).to_csv(os.path.join(run, "decompositions.tsv"), sep="\t", index=False)
        with open(os.path.join(run, "summary.json"), "w") as handle:
            handle.write('{"graph": {"name": "amp", "regions": ["chr1:1-3000"]}, "ambiguity": null}')
        paths = aa_export.export(run, os.path.join(self.dir, "out"))
        truth = [Cycle("t", 7.0, [steps((0, 999, "+"), (1000, 1999, "+"), (2000, 2999, "-"))])]
        old = ev.score_decomposition(truth, ev.read_inferred(os.path.join(run, "decompositions.tsv"))[1]["cycles"], 10_000)[0]
        new = aa_eval.score(paths["cycles"], truth).species.iloc[0]
        self.assertTrue(old.order_match)  # only through the second ordering
        for field, column in (("composite", "composite"), ("jaccard_total", "segment_jaccard"), ("jaccard_hap", "haplotype_jaccard"),
                              ("junction_f1", "junction_f1"), ("phasing_accuracy", "phasing_accuracy"),
                              ("order_match", "order_match"), ("cn_rel_err", "cn_rel_err"), ("inferred", "inferred")):
            self.assertEqual(getattr(old, field), new[column], field)
        records = aa_format.read_cycles(paths["cycles"])
        self.assertEqual([r.fields["Status"] for r in records.cycles], ["species", "low_weight"])
        self.assertEqual(records.intervals, [("chr1", 2, 3001)])


if __name__ == "__main__":
    unittest.main()

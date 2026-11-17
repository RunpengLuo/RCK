"""RCK's simulated-input degradations (rck.ecdna.rck_inputs) on a hand-made DIAG folder.

Run from the RCK root: python -m unittest tests.ecdna.test_rck_inputs
"""
import json
import os
import tempfile
import unittest

import numpy as np

from rck.ecdna import rck_inputs as ri


def seg(start, end, a1, b1, a2=None, b2=None, kind="amplicon"):
    cn = {"clone1": {"A": a1, "B": b1}, "clone2": {"A": a1 if a2 is None else a2, "B": b1 if b2 is None else b2}}
    cnb = {c: {"A": {"l": 0, "u": 99}, "B": {"l": 0, "u": 99}} for c in cn}
    return ri.Segment("chr1", start, end, {"cn": cn, "cnb": cnb, "kind": kind})


def folder(tmp, n_adj=4):
    """chr1:1-400000 in 8 segments of 50 kb; n_adj novel adjacencies at segment ends; 1 labelling + 1 molecule group."""
    os.makedirs(tmp, exist_ok=True)
    segs = [seg(50_000 * k + 1, 50_000 * (k + 1), 2 + k, 1) for k in range(8)]
    ri.write_segments(os.path.join(tmp, "rck.scnt.tsv"), segs)
    adjs = [ri.Adjacency("a{k}".format(k=k), "chr1", 50_000 * (k + 1), "+", "chr1", 50_000 * (k + 4) + 1, "-",
                         {"at": "N", "#NOV": "30", "junction_hap": "A", "hap_obs_end1": "A", "hap_obs_end2": "A"})
            for k in range(n_adj)]
    ri.write_adjacencies(os.path.join(tmp, "rck.adj.tsv"), adjs)
    with open(os.path.join(tmp, "rck.adjgroups.tsv"), "w") as h:
        h.write("gid\taids\textra\nm1\ta0,a1\tagt=M;support=1\nl1\ta0,a2\tagt=L;alignment=a0+,a2+\n")
    with open(os.path.join(tmp, "rck.fragments.tsv"), "w") as h:
        h.write("chr\tstart\tend\textra\nchr1\t1\t400000\tphase_block=chr1:0\n")
    return segs, adjs


class Transforms(unittest.TestCase):
    def test_split_keeps_copy_numbers(self):
        segs = ri.split_at([seg(1, 100, 5, 2)], "chr1", 40, "+")
        self.assertEqual([(s.start, s.end) for s in segs], [(1, 40), (41, 100)])
        self.assertEqual([s.extra["cn"]["clone1"]["A"] for s in segs], [5, 5])
        self.assertEqual(len(ri.split_at(segs, "chr1", 41, "-")), 2)  # already a "-" end at 41

    def test_windows_and_average(self):
        segs, frags = ri.window_fragments([seg(1, 30_000, 4, 0), seg(30_001, 70_000, 10, 2)], width=50_000)
        self.assertEqual([(s.start, s.end) for s in segs], [(1, 30_000), (30_001, 50_000), (50_001, 70_000)])
        self.assertEqual(frags, [[0, 1], [2]])
        ri.average(segs, frags)
        self.assertEqual(segs[0].extra["cn"]["clone1"]["A"], 6.4)  # (4 x 30 kb + 10 x 20 kb) / 50 kb
        self.assertEqual(segs[1].extra["cn"]["clone1"]["B"], 0.8)
        self.assertEqual(segs[2].extra["cn"]["clone1"]["A"], 10)

    def test_long_segment_alone_is_not_cut(self):
        """A 1 Mb segment, then two short ones in one window: only that window is cut and averaged."""
        segs, frags = ri.window_fragments([seg(1, 1_000_000, 1, 1), seg(1_000_001, 1_020_000, 9, 1),
                                           seg(1_020_001, 1_060_000, 3, 1)], width=50_000)
        self.assertEqual([(s.start, s.end) for s in segs],
                         [(1, 1_000_000), (1_000_001, 1_020_000), (1_020_001, 1_050_000), (1_050_001, 1_060_000)])
        self.assertEqual(frags, [[0], [1, 2], [3]])

    def test_flip_rate_and_consistency(self):
        segs = [seg(k * 10 + 1, k * 10 + 10, 3, 1, 7, 0) for k in range(2000)]
        frags = ri.segment_fragments(segs)
        flipped = set(ri.chil_flip(segs, frags, np.random.default_rng(1)))
        self.assertAlmostEqual(len(flipped) / 2000, 0.5, delta=0.04)
        for k in range(2000):
            cn = segs[k].extra["cn"]
            want = ({"A": 1, "B": 3}, {"A": 0, "B": 7}) if k in flipped else ({"A": 3, "B": 1}, {"A": 7, "B": 0})
            self.assertEqual((cn["clone1"], cn["clone2"]), want)  # one swap shared by both clones

    def test_perturb_share_and_step(self):
        segs = [seg(k * 10 + 1, k * 10 + 10, 3, 0) for k in range(500)]
        frags = ri.segment_fragments(segs)
        before = [[s.extra["cn"][c][h] for c in ("clone1", "clone2") for h in "AB"] for s in segs]
        n = ri.perturb(segs, frags, np.random.default_rng(2))
        after = [[s.extra["cn"][c][h] for c in ("clone1", "clone2") for h in "AB"] for s in segs]
        diffs = [a - b for x, y in zip(after, before) for a, b in zip(x, y) if a != b]
        self.assertEqual(n, 100)  # 5% of 500 x 2 clones x 2 alleles
        self.assertEqual(len(diffs), 100)
        self.assertTrue(all(abs(d) == 1 for d in diffs))
        self.assertTrue(all(v >= 0 for row in after for v in row))  # 0 only goes up

    def test_shift_and_spurious(self):
        segs = [seg(1, 400_000, 5, 1)]
        adjs = [ri.Adjacency("a{k}".format(k=k), "chr1", 50_000 * (k + 1), "+", "chr1", 210_000 + 40_000 * k, "-", {"#NOV": "30"})
                for k in range(4)]
        orig = [(a.coord1, a.coord2) for a in adjs]
        moved = ri.shift_breakends(adjs, segs, np.random.default_rng(3))
        self.assertEqual(len(moved), 2)
        for a, (c1, c2) in zip(adjs, orig):
            d = (a.coord1 - c1, a.coord2 - c2)
            if a.aid in moved:
                self.assertTrue(all(1 <= abs(x) <= 50 for x in d))
            else:
                self.assertEqual(d, (0, 0))
        new = ri.spurious_adjacencies(adjs, segs, np.random.default_rng(4))
        self.assertEqual(len(new), 1)  # max(1, round(10% of 4))
        self.assertEqual(new[0].extra["spurious"], "1")


class Degrade(unittest.TestCase):
    def test_variants_write_aligned_label_free_inputs(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "src")
            segs, adjs = folder(src, n_adj=4)
            for variant in ri.VARIANTS:
                dst = os.path.join(tmp, variant)
                rec = ri.degrade(src, dst, variant, seed=[7, 1])
                out = ri.read_segments(os.path.join(dst, "rck.scnt.tsv"))
                ends = {(s.start, "-") for s in out} | {(s.end, "+") for s in out}
                for a in ri.read_adjacencies(os.path.join(dst, "rck.adj.tsv")):
                    self.assertIn((a.coord1, a.strand1), ends)
                    self.assertIn((a.coord2, a.strand2), ends)
                    self.assertFalse(set(ri.LABEL_TAGS) & set(a.extra))
                with open(os.path.join(dst, "rck.adjgroups.tsv")) as h:
                    groups = h.read()
                self.assertIn("agt=M", groups)
                self.assertNotIn("agt=L", groups)
                self.assertEqual(rec["labelling_groups_dropped"], 1)
                with open(os.path.join(dst, "transforms.json")) as h:
                    self.assertEqual(json.load(h)["variant"], variant)
                # total copies per clone over the region are kept by averaging and flips (not by ±1)
                mass = lambda ss: sum((s.extra["cn"]["clone1"]["A"] + s.extra["cn"]["clone1"]["B"]) * s.length for s in ss)
                if variant != "chil_err":
                    self.assertAlmostEqual(mass(out), mass(segs), delta=1e-6 * mass(segs))

    def test_chil_flips_undone_give_the_input(self):
        with tempfile.TemporaryDirectory() as tmp:
            src = os.path.join(tmp, "src")
            segs, _ = folder(src)
            dst = os.path.join(tmp, "chil")
            ri.degrade(src, dst, "chil", seed=[9])
            out = ri.read_segments(os.path.join(dst, "rck.scnt.tsv"))
            with open(os.path.join(dst, "rck.fragments.tsv")) as h:
                frags = h.read().splitlines()[1:]
            self.assertEqual(len(frags), len(out))  # chil: one fragment per segment
            for s, line in zip(out, frags):
                flipped = line.endswith("flipped=1")
                orig = next(x for x in segs if x.start <= s.start <= x.end).extra["cn"]
                for c in orig:
                    want = {"A": orig[c]["B"], "B": orig[c]["A"]} if flipped else orig[c]
                    self.assertEqual(s.extra["cn"][c], want)


if __name__ == "__main__":
    unittest.main()

"""--block-flip: per-phase-block allele observations (segment_blocks.tsv) and A/B flips in the flow MILP.

The fixture is a 10-copy ecDNA on haplotype A over two chr1 segments (``circle_with_blocks``).
Its second segment can sit in a phase block whose observed A/B labels are swapped, as on demo7
when subgraph_construct.py cannot orient a block.
"""

import json
import os
import unittest

import pandas as pd

from rck.core.structures import Haplotype, Phasing
from rck.ecdna.phased import seg_key
from rck.ecdna.plot import read_pieces
from rck.ecdna.scale import ScaleParams, build_observations
from rck.ecdna.subgraph_io import read_segment_blocks, read_subgraph
from tests.ecdna.helpers import requires_gurobi, SubgraphBuilder, temp_dir, write_segment_blocks

SWAPPED = 11 / 12.0  # observed BAF of a 10 A + 1 A + 1 B segment whose block has A and B swapped
TRUE = 1 / 12.0


def circle_with_blocks(directory, seg2_baf=SWAPPED, seg2_orientation=0, split=False, list_unoriented=True):
    """chr1 segments s1 (1001-101000), s2 (101001-201000) and s3 (201001-301000).

    10 ecDNA copies on A run through s1 and s2 and close with junction ``circ`` (A); every segment
    also has 1 + 1 background copies. An imputed edge ``imp`` would close s1 alone. Phase blocks:
    s1 in chr1:0 (oriented), s2 in chr1:1 (observed BAF ``seg2_baf``, orientation
    ``seg2_orientation``) or, with ``split``, half in chr1:1 (oriented, correct) and half in chr1:2
    (unoriented, swapped), and s3 (2 copies, balanced) in chr1:3 (unoriented). With
    ``list_unoriented`` False, phase_blocks.tsv lists only the oriented blocks (the format of
    subgraphs built before segment_blocks.tsv). Returns (directory, segment node pairs).
    """
    pooled = (TRUE + SWAPPED) / 2 if split else seg2_baf
    b = SubgraphBuilder()
    segs = b.chain("chr1", [1000, 101000, 201000, 301000], [180.0, 180.0, 30.0], bafs=[TRUE, pooled, 0.5],
                   ref_reads=[12, 72, 12, 12])
    b.variant(segs[1][1], segs[0][0], 60, hap="A", aid="circ")
    b.imputed(segs[0][1], segs[0][0], 50, aid="imp")
    b.write(directory)
    pieces = [{"edge": segs[0], "block": "chr1:0", "start": 1001, "end": 101000, "total": 12.0, "baf": TRUE}]
    orientation = {"chr1:0": 1}
    if split:
        pieces += [{"edge": segs[1], "block": "chr1:1", "start": 101001, "end": 151001, "total": 12.0, "baf": TRUE},
                   {"edge": segs[1], "block": "chr1:2", "start": 151001, "end": 201000, "total": 12.0, "baf": SWAPPED}]
        orientation.update({"chr1:1": 1, "chr1:2": 0})
    else:
        pieces.append({"edge": segs[1], "block": "chr1:1", "start": 101001, "end": 201000, "total": 12.0, "baf": seg2_baf})
        orientation["chr1:1"] = seg2_orientation
    pieces.append({"edge": segs[2], "block": "chr1:3", "start": 201001, "end": 301000, "total": 2.0, "baf": 0.5})
    orientation["chr1:3"] = 0
    if not list_unoriented:
        orientation = {k: v for k, v in orientation.items() if v != 0}
    write_segment_blocks(directory, pieces, orientation)
    return directory, segs


def observations(directory, mode="unoriented", min_snps_allele=5):
    g = read_subgraph(directory)
    blocks = read_segment_blocks(directory, g) if mode != "off" else None
    params = ScaleParams(depth_per_copy=15.0, junction_reads_per_copy=6.0, block_flip=mode, min_snps_allele=min_snps_allele)
    return g, build_observations(g, params, segment_blocks=blocks)


def solve(directory, mode="unoriented", **params):
    from rck.ecdna.ilp_circular import CircularParams, solve_flow

    g, obs = observations(directory, mode)
    return g, obs, solve_flow(g, obs, CircularParams(**params))


def sid_of(g, nodes):
    return g.nodes[nodes[0]].segment


def row(sol, block):
    return [b for b in sol.blocks if b["phase_block"] == block][0]


def status(sol):
    return {b["phase_block"]: (b["flip"], b["status"]) for b in sol.blocks}


class TestPieceInputs(unittest.TestCase):
    """Reading segment_blocks.tsv and building piece observations; no Gurobi."""

    def test_pieces_and_frames(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp)
            g, obs = observations(directory)
        s1, s2, s3 = (sid_of(g, s) for s in segs)
        self.assertEqual(sorted(obs.pieces), sorted([s1, s2, s3]))
        p1, p2, p3 = obs.pieces[s1][0], obs.pieces[s2][0], obs.pieces[s3][0]
        self.assertTrue(p1.fixed)
        self.assertFalse(p2.fixed)
        self.assertAlmostEqual(p2.a, 1.0)  # unoriented: kept in the observed (swapped) frame
        self.assertAlmostEqual(p2.b, 11.0)
        self.assertEqual(p2.alleles(1), (p2.b, p2.a))
        self.assertAlmostEqual(p3.share, 1.0)
        self.assertEqual(obs.block_orientation(), {"chr1:0": 1, "chr1:1": 0, "chr1:3": 0})
        # the segment-level observation pools the pieces after the flips
        self.assertAlmostEqual(obs.segment_alleles(s2)[0], 1.0)
        self.assertAlmostEqual(obs.segment_alleles(s2, {"chr1:1": 1})[0], 11.0)

    def test_orientation_minus_one_is_swapped_and_fixed(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp, seg2_baf=SWAPPED, seg2_orientation=-1)
            g, obs = observations(directory)
        p2 = obs.pieces[sid_of(g, segs[1])][0]
        self.assertTrue(p2.fixed)
        self.assertAlmostEqual(p2.a, 11.0)
        self.assertEqual(p2.count_a, 440)  # 20 SNPs x 24 reads x 11 / 12

    def test_block_missing_from_phase_blocks_is_unoriented(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp, list_unoriented=False)
            g, obs = observations(directory)
        self.assertEqual(obs.block_orientation()["chr1:3"], 0)
        self.assertFalse(obs.pieces[sid_of(g, segs[1])][0].fixed)

    def test_low_snp_piece_is_total_only(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp)
            frame = pd.read_csv(os.path.join(directory, "segment_blocks.tsv"), sep="\t")
            frame.loc[frame["phase_block"] == "chr1:1", "n_snps"] = 3
            frame.to_csv(os.path.join(directory, "segment_blocks.tsv"), sep="\t", index=False)
            g, obs = observations(directory)
            _, obs_strict = observations(directory, min_snps_allele=30)
        p2 = obs.pieces[sid_of(g, segs[1])][0]
        self.assertFalse(p2.has_alleles)
        self.assertAlmostEqual(p2.total, 12.0)
        self.assertTrue(obs.pieces[sid_of(g, segs[0])][0].has_alleles)  # 20 SNPs
        self.assertFalse(any(p.has_alleles for ps in obs_strict.pieces.values() for p in ps))

    def test_missing_file_falls_back_with_warning(self):
        with temp_dir() as tmp:
            directory, _ = circle_with_blocks(tmp)
            os.remove(os.path.join(directory, "segment_blocks.tsv"))
            g = read_subgraph(directory)
            with self.assertLogs("rck.ecdna.subgraph_io", level="WARNING") as captured:
                self.assertIsNone(read_segment_blocks(directory, g))
            self.assertIn("segment_blocks.tsv", captured.output[0])
            _, obs = observations(directory)
        self.assertEqual(obs.pieces, {})

    def test_mapping_through_phantom_contraction(self):
        """A piece of a contracted phantom goes to the segment that absorbed its base; unknown edges are dropped."""
        b = SubgraphBuilder()
        c = "chr2"
        for nid, pos, o in (("chr2:1000-_+", 1000, "+"), ("chr2:1000-_-", 1001, "-"), ("chr2:5000+_+", 5000, "+"),
                            ("chr2:5000+_-", 5001, "-"), ("chr2:5001-_+", 5002, "+"), ("chr2:5001-_-", 5003, "-"),
                            ("chr2:9000+_+", 9000, "+"), ("chr2:9000+_-", 9001, "-"), ("x_+", 20000, "+"), ("x_-", 30000, "-")):
            b.node(nid, c, pos, o)
        b.edge("chr2:1000-_+", "chr2:1000-_-", "extremity", 10)
        b.edge("chr2:5000+_+", "chr2:5000+_-", "extremity", 8)
        b.edge("chr2:5001-_+", "chr2:5001-_-", "extremity", 12)
        b.edge("chr2:9000+_+", "chr2:9000+_-", "extremity", 10)
        b.edge("chr2:1000-_-", "chr2:5000+_+", "reference", 100, n_snps=10, BAF=0.1)
        b.edge("chr2:5000+_-", "chr2:5001-_+", "reference", 100, n_snps=0)  # a 2 bp phantom with a span
        b.edge("chr2:5001-_-", "chr2:9000+_+", "reference", 100, n_snps=10, BAF=0.1)
        b.edge("x_+", "x_-", "reference", 100, n_snps=10, BAF=0.1)
        with temp_dir() as tmp:
            directory = b.write(tmp)
            write_segment_blocks(directory, [
                {"edge": ("chr2:1000-_-", "chr2:5000+_+"), "block": "chr2:0", "start": 1001, "end": 5000, "total": 6.0, "baf": 0.1},
                {"edge": ("chr2:5000+_-", "chr2:5001-_+"), "block": "chr2:0", "start": 5001, "end": 5002, "total": 6.0, "baf": None,
                 "n_snps": 0},
                {"edge": ("chr2:5001-_-", "chr2:9000+_+"), "block": "chr2:0", "start": 5003, "end": 9000, "total": 6.0, "baf": 0.1},
                {"edge": ("x_+", "x_-"), "block": "chr2:0", "start": 20000, "end": 30000, "total": 6.0, "baf": 0.1}])
            g = read_subgraph(directory)
            os.remove(os.path.join(directory, "edges.tsv"))  # only the graph and segment_blocks.tsv are used from here
            edges_kept = [s for s in g.segments.values() if s.chrom == "chr2" and s.start < 10000]
            self.assertEqual(len(edges_kept), 2)  # the phantom was contracted
            right = [s for s in edges_kept if s.tail == "chr2:5001-_-"][0]
            blocks = read_segment_blocks(directory, g)
        by_start = {p.start: p.sid for p in blocks.pieces}
        self.assertEqual(by_start[5001], right.sid)
        self.assertEqual(by_start[5003], right.sid)
        self.assertNotEqual(by_start[1001], right.sid)
        self.assertEqual(blocks.unmapped, 0)


@requires_gurobi
class TestBlockFlipFlow(unittest.TestCase):
    """Small MILPs; within the size-limited Gurobi licence."""

    def test_swapped_unoriented_block_is_flipped(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp)
            g, obs, sol = solve(directory)
        with temp_dir() as tmp:
            reference = solve(circle_with_blocks(tmp, seg2_baf=TRUE)[0])[2]
        s1, s2 = sid_of(g, segs[0]), sid_of(g, segs[1])
        self.assertEqual(status(sol)["chr1:1"], (1, "flipped"))
        self.assertEqual(status(sol)["chr1:0"], (0, "fixed"))
        for sid in (s1, s2):
            self.assertEqual(sol.seg_ec[sid, Haplotype.A], 10, sid)
            self.assertEqual(sol.seg_ec[sid, Haplotype.B], 0, sid)
        self.assertEqual(sol.adj_ec["circ", Phasing.AA], 10)
        self.assertEqual(sol.imputed_used, [])
        self.assertAlmostEqual(sol.objective, reference.objective, places=4)  # only the 1e-6 flip tie-break differs
        self.assertEqual(status(reference)["chr1:1"][0], 0)
        self.assertGreater(row(sol, "chr1:1")["delta_if_reversed"], 1.0)  # reversing the flip costs the swapped misfit

    def test_off_is_unchanged(self):
        """--block-flip off ignores segment_blocks.tsv: the same flow as without the file."""
        with temp_dir() as tmp:
            directory, _ = circle_with_blocks(tmp)
            g, obs, off = solve(directory, mode="off")
            os.remove(os.path.join(directory, "segment_blocks.tsv"))
            _, _, without = solve(directory, mode="off")
            _, _, on = solve(circle_with_blocks(os.path.join(tmp, "on"))[0])
        self.assertEqual(obs.pieces, {})
        self.assertEqual(off.objective, without.objective)
        self.assertEqual(off.seg_ec, without.seg_ec)
        self.assertEqual(off.adj_ec, without.adj_ec)
        self.assertEqual((off.block_flips, off.blocks), ({}, []))
        self.assertGreater(off.objective, on.objective + 1.0)  # the swapped segment costs the off run

    def test_oriented_block_is_never_flipped(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp, seg2_orientation=1)
            g, obs, sol = solve(directory)
            _, _, opt_in = solve(directory, flip_oriented_penalty=0.5)
        self.assertEqual(status(sol)["chr1:1"], (0, "fixed"))
        self.assertEqual(status(sol)["chr1:0"], (0, "fixed"))
        # with --flip-oriented-penalty an oriented block whose votes are wrong may flip
        self.assertEqual(status(opt_in)["chr1:1"], (1, "flipped"))

    def test_block_without_evidence_is_unresolved(self):
        with temp_dir() as tmp:
            directory, _ = circle_with_blocks(tmp)
            _, _, sol = solve(directory)
        flip, state = status(sol)["chr1:3"]
        self.assertEqual((flip, state), (0, "unresolved"))
        self.assertLess(abs(row(sol, "chr1:3")["delta_if_reversed"]), 0.01)
        self.assertEqual(sol.block_flips["chr1:3"], 0)

    def test_segment_over_oriented_and_unoriented_block(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp, split=True)
            g, obs, sol = solve(directory)
        s2 = sid_of(g, segs[1])
        self.assertEqual([p.block for p in obs.pieces[s2]], ["chr1:1", "chr1:2"])
        self.assertEqual(status(sol)["chr1:1"], (0, "fixed"))
        self.assertEqual(status(sol)["chr1:2"], (1, "flipped"))
        self.assertEqual(sol.seg_ec[s2, Haplotype.A], 10)
        self.assertEqual(sol.seg_ec[s2, Haplotype.B], 0)
        self.assertEqual(sol.imputed_used, [])
        a, b = obs.segment_alleles(s2, sol.block_flips)
        self.assertAlmostEqual(a, 11.0)
        self.assertAlmostEqual(b, 1.0)

    def test_refit_uses_flipped_observations(self):
        from rck.ecdna.ilp_circular import CircularParams
        from rck.ecdna.refit import observation_rows

        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(tmp)
            g, obs, sol = solve(directory)
        s2 = sid_of(g, segs[1])
        rows = observation_rows(g, obs, sol, CircularParams())
        a_rows = [r for r in rows if r.keys == (seg_key(s2, Haplotype.A),)]
        self.assertEqual(len(a_rows), 1)
        self.assertAlmostEqual(a_rows[0].target, 11.0 - sol.seg_bg[s2, Haplotype.A])
        # without the flip the target would be the swapped observation
        sol.block_flips = {}
        unflipped = [r for r in observation_rows(g, obs, sol, CircularParams()) if r.keys == (seg_key(s2, Haplotype.A),)]
        self.assertAlmostEqual(unflipped[0].target, 1.0 - sol.seg_bg[s2, Haplotype.A])


@requires_gurobi
class TestBlockFlipRun(unittest.TestCase):
    def run_cli(self, directory, outdir, extra):
        from rck.ecdna_run import build_parser, run

        argv = ["--subgraph-dir", directory, "--outdir", outdir, "--depth-per-copy", "15", "--junction-reads-per-copy", "6",
                "--no-plots", "--num-alternatives", "2"] + extra
        return run(build_parser().parse_args(argv))

    def test_outputs(self):
        with temp_dir() as tmp:
            directory, segs = circle_with_blocks(os.path.join(tmp, "amp"))
            out = os.path.join(tmp, "out")
            result = self.run_cli(directory, out, ["--block-flip", "unoriented"])
            blocks = pd.read_csv(os.path.join(out, "block_flips.tsv"), sep="\t")
            with open(os.path.join(out, "summary.json")) as handle:
                summary = json.load(handle)
            with open(os.path.join(out, "observations.json")) as handle:
                observations_json = json.load(handle)
            segments = pd.read_csv(os.path.join(out, "flow_segments.tsv"), sep="\t")
            pieces = read_pieces(os.path.join(out, "observations.json"), segments)
            off = os.path.join(tmp, "off")
            self.run_cli(directory, off, [])
            off_files = set(os.listdir(off))
            with open(os.path.join(off, "observations.json")) as handle:
                off_json = json.load(handle)
        self.assertEqual(list(blocks.columns), ["phase_block", "orientation", "flip", "status", "n_pieces", "n_snps", "delta_if_reversed"])
        self.assertEqual(dict(zip(blocks["phase_block"], blocks["status"])),
                         {"chr1:0": "fixed", "chr1:1": "flipped", "chr1:3": "unresolved"})
        self.assertEqual(summary["block_flip"]["mode"], "unoriented")
        self.assertEqual((summary["block_flip"]["flipped"], summary["block_flip"]["unresolved"]), (1, 1))
        flipped = {p["block"]: p["flipped"] for p in observations_json["pieces"]}
        self.assertEqual(flipped, {"chr1:0": False, "chr1:1": True, "chr1:3": False})
        s2 = result["graph"].nodes[segs[1][0]].segment
        row = segments[segments["sid"] == s2].iloc[0]
        self.assertAlmostEqual(row["obs_A"], 11.0, places=3)  # oriented observation after the flip
        self.assertEqual(len(segments), 3)  # still one row per segment
        top = result["decompositions"][0]
        self.assertEqual([s.weight for s in top.species], [10.0])
        piece = pieces[pieces["sid"] == s2].iloc[0]  # plot points: the flipped piece observation
        self.assertEqual((piece["start"], piece["end"], round(piece["obs_A"], 6)), (101001, 200999, 11.0))
        self.assertNotIn("block_flips.tsv", off_files)
        self.assertNotIn("pieces", off_json)

    def test_missing_segment_blocks_falls_back(self):
        with temp_dir() as tmp:
            directory, _ = circle_with_blocks(os.path.join(tmp, "amp"))
            os.remove(os.path.join(directory, "segment_blocks.tsv"))
            with self.assertLogs("rck.ecdna.subgraph_io", level="WARNING"):
                result = self.run_cli(directory, os.path.join(tmp, "out"), ["--block-flip", "unoriented"])
        self.assertEqual(result["observations"].pieces, {})
        self.assertEqual(result["flow"].blocks, [])


if __name__ == "__main__":
    unittest.main()

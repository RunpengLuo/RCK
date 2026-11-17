"""rck-ecdna-plot on multi-clone runs (rck.ecdna.plot.multiclone_segments, read_species psi / status).

Run from the RCK root: python -m unittest tests.ecdna.test_plot_multiclone
"""
import json
import os
import tempfile
import unittest

from rck.ecdna import plot


def _segments(path, rows):
    with open(path, "w") as h:
        h.write("chr\tstart\tend\textra\n")
        for chrom, a, b, extra in rows:
            h.write("{c}\t{a}\t{b}\t{e}\n".format(c=chrom, a=a, b=b, e=extra))


class MultiCloneSegments(unittest.TestCase):
    def test_bulk_terms_and_fragment_orientation(self):
        """u = 0.5 / 0.3 (normal 0.2). Fragment f1's input is in the swapped phase frame relative to the inferred
        karyotype, so its observed A / B are drawn swapped; f2's is not."""
        with tempfile.TemporaryDirectory() as tmp:
            run, inp = os.path.join(tmp, "run"), os.path.join(tmp, "inp")
            os.makedirs(run)
            os.makedirs(inp)
            with open(os.path.join(run, "summary.json"), "w") as h:
                json.dump({"rck_input": inp, "clone_proportions": {"clone1": 0.5, "clone2": 0.3}}, h)
            _segments(os.path.join(run, "rck.scnt.tsv"), [
                ("chr1", 1, 100, "cn={'clone1': {'A': 10, 'B': 1}, 'clone2': {'A': 6, 'B': 1}}"),
                ("chr1", 101, 200, "cn={'clone1': {'A': 2, 'B': 4}, 'clone2': {'A': 2, 'B': 4}}")])
            _segments(os.path.join(inp, "rck.scnt.tsv"), [
                ("chr1", 1, 100, "cn={'clone1': {'A': 1, 'B': 10}, 'clone2': {'A': 1, 'B': 6}}"),
                ("chr1", 101, 200, "cn={'clone1': {'A': 2, 'B': 4}, 'clone2': {'A': 2, 'B': 4}}")])
            _segments(os.path.join(inp, "rck.fragments.tsv"), [("chr1", 1, 100, "x=f1"), ("chr1", 101, 200, "x=f2")])
            seg = plot.multiclone_segments(run)
        self.assertEqual(list(seg["start"]), [1, 101])
        self.assertAlmostEqual(seg["ecdna_A"].iloc[0], 0.5 * 10 + 0.3 * 6 + 0.2)   # inferred bulk A
        self.assertAlmostEqual(seg["obs_A"].iloc[0], 0.5 * 10 + 0.3 * 6 + 0.2)     # input swapped back
        self.assertAlmostEqual(seg["obs_B"].iloc[0], 0.5 * 1 + 0.3 * 1 + 0.2)
        self.assertAlmostEqual(seg["obs_A"].iloc[1], 0.5 * 2 + 0.3 * 2 + 0.2)      # f2 kept as given
        self.assertTrue((seg["bg_A"] == 0).all())

    def test_read_species_psi_and_chromosomal(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "decompositions.tsv")
            with open(path, "w") as h:
                h.write("decomposition_rank\tspecies\ttopology\tstatus\tweight\tweight_int\tlength_bp\tordering_rank\twalk\tjoins"
                        "\tpsi_clone1\tpsi_clone2\n")
                h.write("1\tec1\tcircular\tspecies\t27.6\t\t100\t1\tchr1:0-99:A:+\treference:r1:AA\t42\t22\n")
                h.write("1\tpath1\tlinear\tchromosomal\t0.8\t\t100\t1\tchr1:0-99:B:+\tend:x:B\t1\t1\n")
            got = plot.read_species(path)
            self.assertEqual([s.name for s in got], ["ec1"])
            self.assertEqual(got[0].psi, {"clone1": 42.0, "clone2": 22.0})
            self.assertEqual(len(plot.read_species(path, include_chromosomal=True)), 2)


if __name__ == "__main__":
    unittest.main()

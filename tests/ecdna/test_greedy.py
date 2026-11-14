"""Heaviest-first decomposition of per-clone RCK karyotypes (``rck.ecdna.greedy``, multi-clone plan
step 4) on tiny synthetic DIAGs."""

import os
import shutil
import tempfile
import unittest

import pandas as pd

from tests.ecdna.helpers import requires_gurobi
from tests.ecdna.test_multiclone import write_diag

CHROM = "chr1"


def write_two_species(path, copies):
    """Two separate ecDNAs on one chromosome: S1 closed by j1 and S2 closed by j2, with a normal
    segment M between them and 100 kb flanks. ``copies``: clone -> (e1, e2) extra copies on A."""
    os.makedirs(path, exist_ok=True)
    clones = sorted(copies)
    segs = [((1, 100_000), None), ((100_001, 110_000), 0), ((110_001, 210_000), None), ((210_001, 220_000), 1), ((220_001, 320_000), None)]
    ub = 2 * (1 + max(max(v) for v in copies.values())) + 10
    cnb = "{" + ", ".join("'{c}': {{'A': {{'l': 0, 'u': {u}}}, 'B': {{'l': 0, 'u': {u}}}}}".format(c=c, u=ub) for c in clones) + "}"
    with open(os.path.join(path, "rck.scnt.tsv"), "w") as h:
        h.write("chr\tstart\tend\textra\n")
        for (start, end), amp in segs:
            cn = "{" + ", ".join("'{c}': {{'A': {a}, 'B': 1}}".format(c=c, a=1 + (copies[c][amp] if amp is not None else 0)) for c in clones) + "}"
            h.write("{ch}\t{s}\t{e}\tcn={cn};cnb={b}\n".format(ch=CHROM, s=start, e=end, cn=cn, b=cnb))
    with open(os.path.join(path, "rck.adj.tsv"), "w") as h:
        h.write("aid\tchr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\textra\n")
        h.write("j1\t{c}\t110000\t+\t{c}\t100001\t-\tat=N;#NOV=50\n".format(c=CHROM))
        h.write("j2\t{c}\t220000\t+\t{c}\t210001\t-\tat=N;#NOV=50\n".format(c=CHROM))
    with open(os.path.join(path, "rck.telomeres.tsv"), "w") as h:
        h.write("chr\tcoord\tstrand\n{c}\t1\t-\n{c}\t320000\t+\n".format(c=CHROM))
    with open(os.path.join(path, "rck.fragments.tsv"), "w") as h:
        h.write("chr\tstart\tend\textra\n" + "".join("{c}\t{s}\t{e}\tphase_block={i}\n".format(c=CHROM, s=s, e=e, i=i)
                                                     for i, ((s, e), _) in enumerate(segs)))
    with open(os.path.join(path, "clone_proportions.tsv"), "w") as h:
        h.write("clone\tu\n" + "".join("{c}\t{u}\n".format(c=c, u=round(0.8 / len(clones), 4)) for c in clones))
    return path


@requires_gurobi
class GreedyTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="rck_greedy_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def run_greedy(self, path, **greedy):
        from rck.ecdna.gurobi_env import close_model
        from rck.ecdna.greedy import GreedyParams
        from rck.ecdna.multiclone import MultiCloneParams, decompose_karyotypes, load_rck_input, solve

        inp = load_rck_input(path)
        model = solve(inp, MultiCloneParams(threads=1, time_limit=30))
        self.assertAlmostEqual(model.gm.ObjVal, 0.0)
        scnt, acnt = model.get_scnt_from_model(), model.get_acnt_from_model()
        close_model(model.gm)
        params = GreedyParams(threads=1, time_limit=30, **greedy)
        result, decomposition = decompose_karyotypes(inp, scnt, acnt, params)
        return inp, result, decomposition

    def test_single_cycle_two_clones(self):
        inp, result, d = self.run_greedy(write_diag(os.path.join(self.tmp, "d"), {"clone1": 3, "clone2": 8}))
        self.assertEqual(result.stopped, "exhausted")
        self.assertEqual(result.remaining_weight, 0)
        first = result.walks[0]
        self.assertEqual(first.kind, "cycle")
        self.assertEqual(first.psi, {"clone1": 3, "clone2": 8})
        self.assertEqual(sum(v for k, v in first.counts.items() if k[0] == "S"), 2)  # S1 and S2 once each
        paths = result.walks[1:]
        self.assertEqual([w.kind for w in paths], ["path", "path"])  # one chromosome per haplotype, shared by the clones
        for w in paths:
            self.assertEqual(w.psi, {"clone1": 1, "clone2": 1})
            self.assertEqual(sum(v for k, v in w.counts.items() if k[0] == "S"), 4)  # telomere to telomere, once
        self.assertEqual(d.call, "ecDNA")

    def test_two_species_per_clone_psi(self):
        inp, result, d = self.run_greedy(write_two_species(os.path.join(self.tmp, "d"), {"clone1": (3, 8), "clone2": (6, 2)}))
        cycles = [w for w in result.walks if w.kind == "cycle"]
        self.assertEqual(len(cycles), 2)
        self.assertEqual(cycles[0].psi, {"clone1": 8, "clone2": 2})  # 10 copies in all: extracted first
        self.assertEqual(cycles[1].psi, {"clone1": 3, "clone2": 6})
        self.assertEqual(result.remaining_weight, 0)
        self.assertEqual(sorted(w.kind for w in result.walks).count("path"), 2)

    def test_copies_rule_takes_the_cycle_before_a_looping_path(self):
        """With 100 kb flanks and a 20 kb ecDNA of 5 copies, CoRAL's length-weighted rule prefers a
        chromosome path that loops through the amplicon (1 x 260 kb > 5 x 20 kb); AA's copy rule
        takes the cycle with all 5 copies first."""
        diag = write_diag(os.path.join(self.tmp, "d"), {"clone1": 5})
        inp, result, d = self.run_greedy(diag, heaviest="copies")
        self.assertEqual(result.walks[0].kind, "cycle")
        self.assertEqual(result.walks[0].psi, {"clone1": 5})
        shutil.rmtree(diag)
        diag = write_diag(os.path.join(self.tmp, "d"), {"clone1": 5})
        inp, result, d = self.run_greedy(diag, heaviest="length")
        self.assertEqual(result.walks[0].kind, "path")
        looped = max(v for k, v in result.walks[0].counts.items() if k[0] == "S")
        self.assertGreater(looped, 1)
        cycles = [w for w in result.walks if w.kind == "cycle"]
        self.assertLess(sum(w.psi["clone1"] for w in cycles), 5)

    def test_decompositions_table_and_evaluation_reader(self):
        from rck.ecdna.evaluate import read_inferred
        from rck.ecdna.multiclone import write_decompositions

        inp, result, d = self.run_greedy(write_two_species(os.path.join(self.tmp, "d"), {"clone1": (3, 8), "clone2": (6, 2)}))
        path = write_decompositions(self.tmp, inp, result, d)
        frame = pd.read_csv(path, sep="\t")
        top = frame[frame["ordering_rank"] == 1]
        self.assertEqual(sorted(top["species"]), ["ec1", "ec2", "path1", "path2"])
        self.assertEqual(set(top[top["kind"] == "path"]["status"]), {"chromosomal"})
        ec1 = top[top["species"] == "ec1"].iloc[0]
        self.assertEqual((ec1["psi_clone1"], ec1["psi_clone2"]), (8, 2))
        self.assertAlmostEqual(ec1["proportion_clone1"], 8 / 11)
        self.assertAlmostEqual(ec1["proportion_clone2"], 2 / 8)
        self.assertAlmostEqual(ec1["proportion_bulk"], (0.4 * 8 + 0.4 * 2) / (0.4 * 11 + 0.4 * 8))
        self.assertAlmostEqual(ec1["weight"], 0.4 * 8 + 0.4 * 2)  # bulk copies per cell
        inferred = read_inferred(path)
        self.assertEqual(sorted(c.name for c in inferred[1]["cycles"]), ["ec1", "ec2"])  # paths are not species


if __name__ == "__main__":
    unittest.main()

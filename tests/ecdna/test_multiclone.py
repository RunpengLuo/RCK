"""The multi-clone route (``OptModelMultiCloneEcDNA``, ``rck-ecdna --model multiclone``) on tiny
synthetic DIAGs in the formats of amplicon-allele-analysis ``simulations/rck_diag.py``.

The basic DIAG is one chromosome: a 100 kb flank (1/1), two 10 kb amplified segments S1, S2 in
separate phase blocks, a 100 kb flank (1/1), and one novel adjacency S2+ -> S1- that closes an
ecDNA over S1-S2. Each clone carries 1 + e copies of haplotype A on S1 and S2.
"""

import json
import os
import shutil
import tempfile
import unittest

from tests.ecdna.helpers import requires_gurobi

FLANK, AMP = 100_000, 10_000
F1 = (1, FLANK)
S1 = (FLANK + 1, FLANK + AMP)
S2 = (FLANK + AMP + 1, FLANK + 2 * AMP)
F2 = (FLANK + 2 * AMP + 1, 2 * FLANK + 2 * AMP)
CHROM = "chr1"


def _cn(per_clone):
    return "{" + ", ".join("'{c}': {{'A': {a}, 'B': {b}}}".format(c=c, a=a, b=b) for c, (a, b) in sorted(per_clone.items())) + "}"


def _cnb(clones, ub):
    return "{" + ", ".join("'{c}': {{'A': {{'l': 0, 'u': {u}}}, 'B': {{'l': 0, 'u': {u}}}}}".format(c=c, u=ub) for c in sorted(clones)) + "}"


def write_diag(path, extra_copies, swap_s2=False, junction=True, imputed_junction=False, nov_reads=None, proportions=None, ub=None):
    """Write a DIAG folder. ``extra_copies``: clone -> e (ecDNA copies on haplotype A).
    ``swap_s2`` gives S2's input in the other allele frame (as an unoriented phase block would)."""
    os.makedirs(path, exist_ok=True)
    clones = sorted(extra_copies)
    ub = ub or max(10, 2 * (1 + max(extra_copies.values())))
    rows = []
    for (start, end), kind, amp, swap in ((F1, "flank", False, False), (S1, "subgraph", True, False), (S2, "subgraph", True, swap_s2),
                                           (F2, "flank", False, False)):
        per_clone = {}
        for c in clones:
            a, b = (1 + extra_copies[c], 1) if amp else (1, 1)
            per_clone[c] = (b, a) if swap else (a, b)
        rows.append("{ch}\t{s}\t{e}\tcn={cn};cnb={cnb};kind={k}".format(ch=CHROM, s=start, e=end, cn=_cn(per_clone), cnb=_cnb(clones, ub), k=kind))
    with open(os.path.join(path, "rck.scnt.tsv"), "w") as h:
        h.write("chr\tstart\tend\textra\n" + "\n".join(rows) + "\n")
    adj = "j1\t{ch}\t{e}\t+\t{ch}\t{s}\t-\tat=N;#NOV={r};#REF=0,0;junction_hap=A;imputed={i};edge_type={t}\n"
    reads = nov_reads if nov_reads is not None else 100
    header = "aid\tchr1\tcoord1\tstrand1\tchr2\tcoord2\tstrand2\textra\n"
    with open(os.path.join(path, "rck.adj.tsv"), "w") as h:
        h.write(header)
        if junction and not imputed_junction:
            h.write(adj.format(ch=CHROM, e=S2[1], s=S1[0], r=reads, i=0, t="variant"))
    with open(os.path.join(path, "rck.adj.imputed.tsv"), "w") as h:
        h.write(header)
        if imputed_junction:
            h.write(adj.format(ch=CHROM, e=S2[1], s=S1[0], r=reads, i=1, t="imputed_variant").replace("j1", "imp1"))
    with open(os.path.join(path, "rck.telomeres.tsv"), "w") as h:
        h.write("chr\tcoord\tstrand\n{ch}\t{a}\t-\n{ch}\t{b}\t+\n".format(ch=CHROM, a=F1[0], b=F2[1]))
    with open(os.path.join(path, "rck.fragments.tsv"), "w") as h:
        h.write("chr\tstart\tend\textra\n")
        for i, (start, end) in enumerate((F1, S1, S2, F2)):
            h.write("{ch}\t{s}\t{e}\tphase_block={ch}:{i}\n".format(ch=CHROM, s=start, e=end, i=i))
    with open(os.path.join(path, "rck.adjgroups.tsv"), "w") as h:
        h.write("gid\taids\textra\n")
    with open(os.path.join(path, "clone_proportions.tsv"), "w") as h:
        h.write("clone\tu\n" + "".join("{c}\t{u}\n".format(c=c, u=(proportions or {}).get(c, round(0.8 / len(clones), 4))) for c in clones))
    return path


def sid(interval):
    return "{ch}:{s}-{e}".format(ch=CHROM, s=interval[0], e=interval[1])


@requires_gurobi
class MultiCloneTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="rck_multiclone_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def solve(self, extra_copies, params=None, theta=None, **diag):
        from rck.ecdna.multiclone import MultiCloneParams, load_rck_input, solve

        params = params or MultiCloneParams(threads=1, time_limit=30)
        path = write_diag(os.path.join(self.tmp, "rck_input"), extra_copies, **diag)
        inp = load_rck_input(path, read_imputed=params.imputed_penalty is not None)
        model = solve(inp, params, theta)
        self.addCleanup(model.gm.dispose)
        return inp, model

    @staticmethod
    def cn(model, layer, clone, interval):
        from rck.core.structures import Haplotype

        if layer is None:
            scnt = model.get_scnt_from_model()
        else:
            scnt = model.get_layer_scnt_from_model(layer)
        return tuple(scnt[clone].get_cn(sid=sid(interval), haplotype=h) for h in (Haplotype.A, Haplotype.B))

    def assert_layer_invariants(self, inp, model):
        """Totals = chr + ec; the ecDNA layer is balanced at every extremity (a circulation); the
        chromosome layer is balanced except at telomeres."""
        from rck.core.ilp_gurobi_ecdna import CHR, EC

        for clone in inp.clone_ids:
            for s in inp.segments:
                interval = (s.start_coordinate, s.end_coordinate)
                total, chr_, ec = (self.cn(model, layer, clone, interval) for layer in (None, CHR, EC))
                self.assertEqual(total, tuple(c + e for c, e in zip(chr_, ec)))
            for layer, telomeres_free in ((EC, False), (CHR, True)):
                for node in model.iag.nodes(data=False):
                    for h in ("A", "B"):
                        from rck.core.structures import Haplotype

                        value = model.get_lin_expr_for_node_haplotype_layer(clone, node, Haplotype(h), layer).getValue()
                        if telomeres_free and node in model.hapl_telomeres:
                            self.assertGreaterEqual(value, -1e-6)
                        else:
                            self.assertAlmostEqual(value, 0.0, places=6, msg="{l} {c} {v} {h}".format(l=layer, c=clone, v=node, h=h))

    def test_single_clone_single_cycle(self):
        from rck.core.ilp_gurobi_ecdna import CHR

        inp, model = self.solve({"clone1": 5})
        self.assertAlmostEqual(model.gm.ObjVal, 0.0)
        self.assertEqual(sorted(self.cn(model, None, "clone1", S1)), [1, 6])
        self.assertEqual(self.cn(model, None, "clone1", S2), self.cn(model, None, "clone1", S1))
        for flank in (F1, F2):
            self.assertEqual(self.cn(model, CHR, "clone1", flank), (1, 1))
        self.assert_layer_invariants(inp, model)

    def test_two_clones_per_clone_copies(self):
        inp, model = self.solve({"clone1": 3, "clone2": 8})
        self.assertAlmostEqual(model.gm.ObjVal, 0.0)
        c1, c2 = self.cn(model, None, "clone1", S1), self.cn(model, None, "clone2", S1)
        self.assertEqual(sorted(c1), [1, 4])
        self.assertEqual(sorted(c2), [1, 9])
        self.assertEqual(c1.index(4), c2.index(9))  # the fragment flip is shared by the clones
        self.assert_layer_invariants(inp, model)
        acnt = model.get_acnt_from_model()
        junction = [a for a in inp.adjacencies if a.extra.get("aid") == "j1"][0]
        from rck.core.structures import Phasing

        used = {c: {ph: acnt[c].get_cn(aid=junction.stable_id_non_phased, phasing=ph) for ph in (Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB)}
                for c in inp.clone_ids}
        ph = Phasing.AA if c1[0] == 4 else Phasing.BB
        self.assertEqual(used["clone1"][ph], 3)
        self.assertEqual(used["clone2"][ph], 8)  # one labelling shared by the clones
        self.assertEqual(sum(used["clone1"].values()) + sum(used["clone2"].values()), 11)

    def test_flipped_phase_block_is_resolved(self):
        from rck.core.ilp_gurobi import FRAGMENT_ALLELE

        inp, model = self.solve({"clone1": 5}, swap_s2=True)
        self.assertAlmostEqual(model.gm.ObjVal, 0.0)
        a1, a2 = self.cn(model, None, "clone1", S1), self.cn(model, None, "clone1", S2)
        self.assertEqual(a1, a2)  # amplified on one haplotype across both blocks
        self.assertEqual(sorted(a1), [1, 6])
        flips = {fid: round(v.X) for fid, v in model.variables[FRAGMENT_ALLELE].items()}
        b1, b2 = flips[inp.segments_to_fragments[sid(S1)]], flips[inp.segments_to_fragments[sid(S2)]]
        self.assertNotEqual(b1, b2)  # exactly one of the two blocks is flipped against the input frame

    # ---------------------------------------------------------------- risk 1 (plan section 5)
    def _probe(self, model, objective, sense):
        import gurobipy as g

        gm = model.gm
        if not getattr(model, "_probe_fixed", False):
            gm.addLConstr(gm.getObjective(), g.GRB.LESS_EQUAL, gm.ObjVal + 1e-6)
            model._probe_fixed = True
        gm.setObjective(objective, sense)
        gm.optimize()
        self.assertEqual(gm.Status, g.GRB.OPTIMAL)
        return gm.ObjVal

    def test_risk1_linear_paths_pinned_by_flanks(self):
        """Among all optimal solutions, chromosome-layer copies leaving through the telomeres never
        exceed the flanks' copy number: the extra copies cannot run telomere to telomere."""
        import gurobipy as g
        from rck.core.ilp_gurobi_ecdna import CHR
        from rck.core.structures import Haplotype

        inp, model = self.solve({"clone1": 5})
        excess = g.quicksum(model.get_lin_expr_for_node_haplotype_layer("clone1", node, h, CHR)
                            for node in set(model.hapl_telomeres) if node in model.iag.graph for h in (Haplotype.A, Haplotype.B))
        self.assertLessEqual(self._probe(model, excess, g.GRB.MAXIMIZE), 4 + 1e-6)  # 2 telomeres x 2 haplotypes x 1 copy

    @unittest.expectedFailure
    def test_risk1_ecdna_layer_identified(self):
        """KNOWN FAILURE (2026-11-13): the optimum does not determine the layer split. The ecDNA
        copies can sit in the chromosome layer as a closed cycle (a circulation is balanced in
        either layer), so the ecDNA layer ranges from 0 to all e copies at the same objective.
        Expected to pass once a rule separates the layers (docs/multiclone_NEXT_STEPS.md)."""
        import gurobipy as g
        from rck.core.ilp_gurobi_ecdna import EC

        inp, model = self.solve({"clone1": 5})
        ec = model.layer_total(EC, sids=[sid(S1)])
        low = self._probe(model, ec, g.GRB.MINIMIZE)
        high = self._probe(model, ec, g.GRB.MAXIMIZE)
        self.assertAlmostEqual(low, high)

    def test_risk1_ecdna_layer_range(self):
        """The current range, pinned so a change in it is noticed: 0 .. e copies on S1."""
        import gurobipy as g
        from rck.core.ilp_gurobi_ecdna import EC

        inp, model = self.solve({"clone1": 5})
        weight = getattr(inp.segments[1], "length_10")
        ec = model.layer_total(EC, sids=[sid(S1)])
        self.assertAlmostEqual(self._probe(model, ec, g.GRB.MINIMIZE), 0.0)
        self.assertAlmostEqual(self._probe(model, ec, g.GRB.MAXIMIZE), 5 * weight)

    # ------------------------------------------------------------- ecDNA additions (off by default)
    def test_imputed_adjacencies_off_by_default(self):
        from rck.ecdna.multiclone import MultiCloneParams

        inp, model = self.solve({"clone1": 5}, imputed_junction=True)
        self.assertEqual(inp.imputed_aids, [])
        self.assertGreater(model.gm.ObjVal, 0.0)  # no closing junction: the amplification cannot be explained
        inp, model = self.solve({"clone1": 5}, MultiCloneParams(threads=1, time_limit=30, imputed_penalty=2.0), imputed_junction=True)
        self.assertEqual(len(inp.imputed_aids), 1)
        self.assertAlmostEqual(model.gm.ObjVal, 2.0)
        self.assertEqual(sorted(self.cn(model, None, "clone1", S1)), [1, 6])

    def test_novel_penalty(self):
        from rck.ecdna.multiclone import MultiCloneParams

        inp, model = self.solve({"clone1": 5}, MultiCloneParams(threads=1, time_limit=30, novel_penalty=0.5))
        self.assertAlmostEqual(model.gm.ObjVal, 0.5)

    def test_junction_term(self):
        from rck.ecdna.multiclone import MultiCloneParams

        params = MultiCloneParams(threads=1, time_limit=30, junction_weight=1.0)
        # u = 0.8, 5 junction copies -> 4 bulk copies; theta 10 -> 40 reads is an exact fit
        inp, model = self.solve({"clone1": 5}, params, theta=10.0, nov_reads=40, proportions={"clone1": 0.8})
        self.assertAlmostEqual(model.gm.ObjVal, 0.0)
        inp, model = self.solve({"clone1": 5}, params, theta=10.0, nov_reads=60, proportions={"clone1": 0.8})
        self.assertAlmostEqual(model.gm.ObjVal, (6.0 - 4.0) / 6.0 ** 0.5)

    def test_linear_ends(self):
        from rck.ecdna.multiclone import MultiCloneParams

        inp, model = self.solve({"clone1": 5}, junction=False)
        self.assertGreater(model.gm.ObjVal, 1000)  # strictly RCK: no way to explain 5 extra copies without a junction
        params = MultiCloneParams(threads=1, time_limit=30, linear_ends="boundaries", linear_end_penalty=0.25)
        inp, model = self.solve({"clone1": 5}, params, junction=False)
        self.assertAlmostEqual(model.gm.ObjVal, 0.25 * 5 * 2)  # 5 copies, a start and an end each
        self.assertEqual(sorted(self.cn(model, None, "clone1", S1)), [1, 6])


@requires_gurobi
class MultiCloneCliTest(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.mkdtemp(prefix="rck_multiclone_cli_")

    def tearDown(self):
        shutil.rmtree(self.tmp, ignore_errors=True)

    def test_route_writes_native_outputs_with_layers(self):
        import ast

        from rck.ecdna_run import main
        from rck.core.io import read_segments_from_file

        diag = write_diag(os.path.join(self.tmp, "sub", "rck_input"), {"clone1": 2, "clone2": 6})
        out = os.path.join(self.tmp, "out")
        self.assertEqual(main(["--model", "multiclone", "--rck-input", diag, "--outdir", out, "--threads", "1", "--time-limit", "30"]), 0)
        with open(os.path.join(out, "summary.json")) as h:
            summary = json.load(h)
        self.assertEqual(summary["model"], "multiclone")
        self.assertEqual(summary["clones"], ["clone1", "clone2"])
        self.assertAlmostEqual(summary["objective"], 0.0)
        segments = {s.stable_id_non_hap: s for s in read_segments_from_file(os.path.join(out, "rck.scnt.tsv"))}
        s1 = segments[sid(S1)]
        layer = s1.extra["layer"] if isinstance(s1.extra["layer"], dict) else ast.literal_eval(s1.extra["layer"])
        for clone, e in (("clone1", 2), ("clone2", 6)):
            total = s1.extra["cn"][clone]
            self.assertEqual(sorted(layer), ["chr", "ec"])
            self.assertEqual({h: layer["chr"][clone][h] + layer["ec"][clone][h] for h in ("A", "B")},
                             {h.value if hasattr(h, "value") else h: v for h, v in total.items()})
            self.assertEqual(sorted(total.values()), [1, 1 + e])
        with open(os.path.join(out, "rck.acnt.tsv")) as h:
            self.assertIn("layer=", h.read())

    def test_model_flags(self):
        from rck.ecdna_run import build_parser

        args = build_parser().parse_args(["--subgraph-dir", "x", "--outdir", "y"])
        self.assertEqual(args.model, "bulk")
        self.assertIsNone(args.mc_imputed_penalty)
        self.assertEqual((args.mc_novel_penalty, args.mc_junction_weight, args.mc_linear_ends), (0.0, 0.0, "none"))
        from rck.ecdna_run import main

        with self.assertRaises(SystemExit):
            main(["--outdir", "y"])  # bulk needs --subgraph-dir
        with self.assertRaises(SystemExit):
            main(["--model", "multiclone", "--outdir", "y"])  # multiclone needs --rck-input


if __name__ == "__main__":
    unittest.main()

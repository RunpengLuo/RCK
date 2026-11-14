"""Multi-clone RCK with a chromosome and an ecDNA layer: ``OptModelMultiCloneEcDNA``.

A minimal subclass of :class:`rck.core.ilp_gurobi.OptModelMultiClone` (Aganezov & Raphael 2020),
for the multi-clone RCK-ecDNA route (``rck-ecdna --model multiclone``; plan
``docs/multiclone_plan.md`` section 4.1, model description and deviation table in
``docs/ecDNA.md``, section "Multi-clone model").

Every variable, constraint and objective term of ``OptModelMultiClone`` is kept, on the same
names, except the balance constraints. On top of them:

* **Two layers.** Every clone's segment copy number c_{k,s,h} (``SEGMENT_COPY_NUMBER``) and phased
  adjacency copy number N_{k,a,ph} (``PROD``/``PY``) is the sum of a chromosome-layer and an
  ecDNA-layer integer variable, one equality each. RCK's own variables stay the totals, so the
  deltas, the objective, the fragment flips, the bounds, the labelling, the 1 - fp presence rule
  and the adjacency groups read the totals exactly as in RCK.
* **Balance per layer** (``define_constraints_on_nodes`` is replaced). The chromosome layer gets
  RCK's rule (excess allowed at telomere positions, balance elsewhere); the ecDNA layer is
  balanced at every extremity, i.e. a circulation. RCK's balance on the totals follows from the
  two and is not written again.

Optional ecDNA additions, all off by default (``EcDNAOptions``):

* ``imputed_aids`` + ``imputed_penalty``: imputed candidate adjacencies join the graph. They are
  left out of RCK's 1 - fp presence rule and cost ``imputed_penalty`` when used (any labelling,
  any clone).
* ``novel_penalty``: a cost per observed novel adjacency used.
* ``junction_weight`` (with ``junction_copies`` and ``clone_proportions``): a read-support term
  for observed novel adjacencies, ``junction_weight / sqrt(max(v, 1)) * |v - sum_k u_k N_k,a|``
  with v = #NOV / theta in bulk copies, as in the bulk route's ``alpha_novel`` term.
* ``linear_end_nodes`` + ``linear_end_penalty``: extra extremities where the chromosome layer
  may have copy-number excess, at this cost per copy (RCK's telomere positions stay free).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

import gurobipy as g

from rck.core.ilp_gurobi import OptModelMultiClone, SEGMENT_COPY_NUMBER, SEGMENT_LENGTH_ATTRIBUTE, PROD, PY, PP
from rck.core.structures import AdjacencyType, Haplotype, Phasing, SegmentCopyNumberProfile, AdjacencyCopyNumberProfile, CNBoundaries, \
    get_aabb_for_ra, get_abba_for_na_and_position

CHR = "chr"
EC = "ec"
LAYERS = (CHR, EC)
SEGMENT_LAYER = "segment_layer"
ADJACENCY_LAYER = "adjacency_layer"
LINEAR_END = "linear_end"
PHASINGS = (Phasing.AA, Phasing.AB, Phasing.BA, Phasing.BB)
HAPLOTYPES = (Haplotype.A, Haplotype.B)


@dataclass
class EcDNAOptions:
    """The ecDNA additions to ``OptModelMultiClone``; the defaults switch every one off."""
    imputed_aids: Iterable[str] = ()  # stable_id_non_phased of imputed adjacencies in hapl_adjacencies
    imputed_penalty: float = 0.0
    novel_penalty: float = 0.0
    junction_weight: float = 0.0
    junction_copies: Dict[str, float] = field(default_factory=dict)  # aid -> #NOV / theta (bulk copies)
    clone_proportions: Dict[str, float] = field(default_factory=dict)  # clone -> u_k
    linear_end_nodes: Iterable[object] = ()  # IAG nodes (Positions) with penalised chromosome-layer excess
    linear_end_penalty: float = 0.0


class OptModelMultiCloneEcDNA(OptModelMultiClone):
    def __init__(self, *args, ecdna: Optional[EcDNAOptions] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.ecdna = ecdna if ecdna is not None else EcDNAOptions()
        self.imputed_aids = set(self.ecdna.imputed_aids)
        telomeres = set(self.hapl_telomeres)
        self.linear_end_nodes = [v for v in self.ecdna.linear_end_nodes if v not in telomeres]
        self.layer_variables = {
            SEGMENT_LAYER: {layer: {cid: {} for cid in self.clone_ids} for layer in LAYERS},
            ADJACENCY_LAYER: {layer: {cid: {} for cid in self.clone_ids} for layer in LAYERS},
            LINEAR_END: {},
        }
        self.extra_objective_terms = []

    # ------------------------------------------------------------------ variables
    def define_variables(self):
        super().define_variables()
        for clone_id in self.clone_ids:
            for s in self.hapl_segments:
                sid = s.stable_id_non_hap
                ub = max(self.scnb[clone_id].get_cnb(sid=sid, hap=h, boundary_type=CNBoundaries.UPPER) for h in HAPLOTYPES)
                for layer in LAYERS:
                    self.layer_variables[SEGMENT_LAYER][layer][clone_id][sid] = {
                        h: self.gm.addVar(lb=0, ub=ub, vtype=g.GRB.INTEGER, name="x_{{{layer},{cid},{sid},{hap}}}".format(layer=layer, cid=clone_id, sid=sid, hap=h))
                        for h in HAPLOTYPES}
            for adjacency in self.hapl_adjacencies:
                aid = adjacency.stable_id_non_phased
                for layer in LAYERS:
                    self.layer_variables[ADJACENCY_LAYER][layer][clone_id][aid] = {
                        ph: self.gm.addVar(lb=0, vtype=g.GRB.INTEGER, name="y_{{{layer},{cid},{aid},{ph}}}".format(layer=layer, cid=clone_id, aid=aid, ph=ph))
                        for ph in PHASINGS}
            for node in self.linear_end_nodes:
                for h in HAPLOTYPES:
                    self.layer_variables[LINEAR_END][clone_id, node, h] = self.gm.addVar(
                        lb=0, vtype=g.GRB.INTEGER, name="e_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=h))

    # ---------------------------------------------------------------- constraints
    def define_constraints(self):
        super().define_constraints()
        self.define_layer_sum_constraints()

    def define_layer_sum_constraints(self):
        """c_{k,s,h} = x^chr + x^ec and N_{k,a,ph} = y^chr + y^ec: RCK's variables are the totals."""
        seg, adj = self.layer_variables[SEGMENT_LAYER], self.layer_variables[ADJACENCY_LAYER]
        for clone_id in self.clone_ids:
            for s in self.hapl_segments:
                sid = s.stable_id_non_hap
                for h in HAPLOTYPES:
                    self.gm.addLConstr(self.variables[SEGMENT_COPY_NUMBER][clone_id][sid][h], g.GRB.EQUAL,
                                       seg[CHR][clone_id][sid][h] + seg[EC][clone_id][sid][h],
                                       name="layers-seg_{{{cid},{sid},{hap}}}".format(cid=clone_id, sid=sid, hap=h))
            for adjacency in self.hapl_adjacencies:
                aid = adjacency.stable_id_non_phased
                for ph in PHASINGS:
                    self.gm.addLConstr(self.variables[PROD][PY][clone_id][aid][ph], g.GRB.EQUAL,
                                       adj[CHR][clone_id][aid][ph] + adj[EC][clone_id][aid][ph],
                                       name="layers-adj_{{{cid},{aid},{ph}}}".format(cid=clone_id, aid=aid, ph=ph))

    def define_constraints_on_nodes(self):
        """Balance per layer. Chromosome layer: RCK's rule (excess >= 0 at telomere positions, == 0
        elsewhere; with linear ends, excess == e >= 0 at those extremities). ecDNA layer: == 0 at
        every extremity."""
        linear_ends = set(self.linear_end_nodes)
        for clone_id in self.clone_ids:
            for node, data in self.iag.nodes(data=True):
                for haplotype in HAPLOTYPES:
                    chr_exp = self.get_lin_expr_for_node_haplotype_layer(clone_id=clone_id, node=node, haplotype=haplotype, layer=CHR)
                    if node in self.hapl_telomeres:
                        self.gm.addLConstr(chr_exp, g.GRB.GREATER_EQUAL, 0, name="cne-chr_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))
                    elif node in linear_ends:
                        self.gm.addLConstr(chr_exp, g.GRB.EQUAL, self.layer_variables[LINEAR_END][clone_id, node, haplotype],
                                           name="cne-chr-linear_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))
                    else:
                        self.gm.addLConstr(chr_exp, g.GRB.EQUAL, 0, name="cnb-chr_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))
                    ec_exp = self.get_lin_expr_for_node_haplotype_layer(clone_id=clone_id, node=node, haplotype=haplotype, layer=EC)
                    self.gm.addLConstr(ec_exp, g.GRB.EQUAL, 0, name="cnb-ec_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))

    def get_lin_expr_for_node_haplotype_layer(self, clone_id, node, haplotype, layer):
        """``get_lin_expr_for_node_haplotype`` on one layer's variables (same terms, self-loops twice)."""
        seg, adj = self.layer_variables[SEGMENT_LAYER][layer][clone_id], self.layer_variables[ADJACENCY_LAYER][layer][clone_id]
        result = g.LinExpr()
        s_u, s_v, data = self.iag.get_segment_edge(node=node, data=True)
        result.add(seg[data["object"].stable_id_non_hap][haplotype])
        ref_adjs = [d["object"] for _, __, d in self.iag.ref_adjacency_edges(data=True, nbunch=node)]
        assert len(ref_adjs) <= 1
        for adjacency in ref_adjs:
            result.add(adj[adjacency.stable_id_non_phased][get_aabb_for_ra(haplotype=haplotype)], -1)
        for adjacency in [d["object"] for _, __, d in self.iag.nov_adjacency_edges(data=True, nbunch=node)]:
            variables = adj[adjacency.stable_id_non_phased]
            if adjacency.is_self_loop_hapl:
                result.add(variables[get_aabb_for_ra(haplotype=haplotype)], -2)
            else:
                result.add(variables[get_aabb_for_ra(haplotype=haplotype)], -1)
                result.add(variables[get_abba_for_na_and_position(novel_adjacency=adjacency, position=node, haplotype=haplotype)], -1)
        return result

    def define_constraints_nov_adjacency_overall_presence(self):
        """RCK's 1 - fp rule over the observed novel adjacencies only; imputed ones are not SV calls.
        Without observed novel adjacencies the rule is skipped: RCK would write 0 >= 1 - fp, which is
        infeasible (an RCK edge case, e.g. an amplicon whose junctions are all imputed)."""
        observed = [a for a in self.hapl_adjacencies if a.adjacency_type == AdjacencyType.NOVEL and a.stable_id_non_phased not in self.imputed_aids]
        if not observed:
            return
        fp_lin_expr = g.LinExpr()
        for adjacency in observed:
            aid = adjacency.stable_id_non_phased
            fp_lin_expr.add(g.quicksum([self.variables[PROD][PP][aid][ph] for ph in PHASINGS]), mult=(1.0 / len(observed)))
        self.gm.addLConstr(fp_lin_expr, g.GRB.GREATER_EQUAL, g.LinExpr(1 - self.hapl_nov_adjacencies_fp), name="nov-adj-fp")

    # ------------------------------------------------------------------ objective
    def define_objective(self):
        super().define_objective()
        terms = self.ecdna_objective_terms()
        if terms:
            self.gm.update()
            self.gm.setObjective(self.gm.getObjective() + g.quicksum(terms), g.GRB.MINIMIZE)

    def ecdna_objective_terms(self) -> List[object]:
        opt, terms = self.ecdna, []
        pp = self.variables[PROD][PP]
        for adjacency in filter(lambda a: a.adjacency_type == AdjacencyType.NOVEL, self.hapl_adjacencies):
            aid = adjacency.stable_id_non_phased
            used = g.quicksum(pp[aid][ph] for ph in PHASINGS)
            if aid in self.imputed_aids:
                if opt.imputed_penalty:
                    terms.append(opt.imputed_penalty * used)
                continue
            if opt.novel_penalty:
                terms.append(opt.novel_penalty * used)
            value = opt.junction_copies.get(aid)
            if opt.junction_weight and value is not None:
                flow = g.quicksum(opt.clone_proportions.get(cid, 0.0) * self.variables[PROD][PY][cid][aid][ph]
                                  for cid in self.clone_ids for ph in PHASINGS)
                d = self.gm.addVar(lb=0, name="d-junction_{{{aid}}}".format(aid=aid))
                self.gm.addLConstr(d, g.GRB.GREATER_EQUAL, flow - value)
                self.gm.addLConstr(d, g.GRB.GREATER_EQUAL, value - flow)
                terms.append(opt.junction_weight / math.sqrt(max(value, 1.0)) * d)
        if self.layer_variables[LINEAR_END] and opt.linear_end_penalty:
            terms.append(opt.linear_end_penalty * g.quicksum(self.layer_variables[LINEAR_END].values()))
        self.extra_objective_terms = terms
        return terms

    # -------------------------------------------------------------------- results
    def get_layer_scnt_from_model(self, layer) -> Dict[str, SegmentCopyNumberProfile]:
        result = {}
        for clone_id in self.clone_ids:
            scnp = SegmentCopyNumberProfile()
            result[clone_id] = scnp
            for s in self.hapl_segments:
                variables = self.layer_variables[SEGMENT_LAYER][layer][clone_id][s.stable_id_non_hap]
                for h in HAPLOTYPES:
                    scnp.set_cn_record_for_segment(segment=s, cn=int(round(variables[h].X)), haplotype=h)
        return result

    def get_layer_acnt_from_model(self, layer) -> Dict[str, AdjacencyCopyNumberProfile]:
        result = {}
        for clone_id in self.clone_ids:
            acnp = AdjacencyCopyNumberProfile()
            result[clone_id] = acnp
            for adj in self.hapl_adjacencies:
                variables = self.layer_variables[ADJACENCY_LAYER][layer][clone_id][adj.stable_id_non_phased]
                for ph in PHASINGS:
                    acnp.set_cn_record_for_adjacency(adjacency=adj, cn=int(round(variables[ph].X)), phasing=ph)
        return result

    def get_linear_ends_from_model(self) -> Dict[tuple, int]:
        return {key: int(round(v.X)) for key, v in self.layer_variables[LINEAR_END].items() if v.X > 0.5}

    def layer_total(self, layer, sids=None, clone_ids=None):
        """Linear expression: copy number x length (``SEGMENT_LENGTH_ATTRIBUTE``) summed over a layer;
        used to probe how the solution splits copies between the layers."""
        attr = self.extra.get(SEGMENT_LENGTH_ATTRIBUTE, "length_100")
        sids = None if sids is None else set(sids)
        expr = g.LinExpr()
        for clone_id in (clone_ids or self.clone_ids):
            for s in self.hapl_segments:
                sid = s.stable_id_non_hap
                if sids is not None and sid not in sids:
                    continue
                for h in HAPLOTYPES:
                    expr.add(self.layer_variables[SEGMENT_LAYER][layer][clone_id][sid][h], max(getattr(s, attr), 1))
        return expr

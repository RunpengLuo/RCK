"""RCK's multi-clone model with the ecDNA route's optional additions: ``OptModelMultiCloneEcDNA``.

A thin subclass of :class:`rck.core.ilp_gurobi.OptModelMultiClone` (Aganezov & Raphael 2020) for
the multi-clone RCK-ecDNA route (``rck-ecdna --model multiclone``; plan ``docs/multiclone_plan.md``
section 4.1, deviation table in ``docs/ecDNA.md`` section 11.1).

With the default ``EcDNAOptions`` it builds exactly RCK's model (same variables, constraints and
objective), with one exception: RCK's 1 - fp rule is skipped when there are no observed novel
adjacencies, where RCK would write 0 >= 1 - fp, which is infeasible.

The karyotype is not split into a chromosome and an ecDNA layer: a circulation balances in either
layer, so RCK's objective could not tell them apart (risk 1, plan section 5). Which copies form
circles and which linear chromosomes is decided by the decomposition instead (plan 4.3), as in
AmpliconArchitect and CoRAL.

Optional additions, all off by default:

* ``imputed_aids`` + ``imputed_penalty``: imputed candidate adjacencies join the graph. They are
  left out of RCK's 1 - fp rule and cost ``imputed_penalty`` when used (any labelling, any clone).
* ``novel_penalty``: a cost per observed novel adjacency used.
* ``junction_weight`` (with ``junction_copies`` and ``clone_proportions``): a read-support term
  for observed novel adjacencies, ``junction_weight / sqrt(max(v, 1)) * |v - sum_k u_k N_k,a|``
  with v = #NOV / theta in bulk copies, as in the bulk route's ``alpha_novel`` term.
* ``linear_end_nodes`` + ``linear_end_penalty``: extra extremities where copy-number excess is
  allowed, like RCK's telomere positions but at this cost per copy (RCK's telomeres stay free).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import Dict, Iterable, List, Optional

import gurobipy as g

from rck.core.ilp_gurobi import OptModelMultiClone, PROD, PY, PP
from rck.core.structures import AdjacencyType, Haplotype, Phasing

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
    linear_end_nodes: Iterable[object] = ()  # IAG nodes (Positions) with penalised excess
    linear_end_penalty: float = 0.0


class OptModelMultiCloneEcDNA(OptModelMultiClone):
    def __init__(self, *args, ecdna: Optional[EcDNAOptions] = None, **kwargs):
        super().__init__(*args, **kwargs)
        self.ecdna = ecdna if ecdna is not None else EcDNAOptions()
        self.imputed_aids = set(self.ecdna.imputed_aids)
        telomeres = set(self.hapl_telomeres)
        self.linear_end_nodes = [v for v in self.ecdna.linear_end_nodes if v not in telomeres]
        self.linear_end_variables = {}

    def define_variables(self):
        super().define_variables()
        for clone_id in self.clone_ids:
            for node in self.linear_end_nodes:
                for h in HAPLOTYPES:
                    self.linear_end_variables[clone_id, node, h] = self.gm.addVar(
                        lb=0, vtype=g.GRB.INTEGER, name="e_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=h))

    def define_constraints_on_nodes(self):
        """RCK's balance; at a linear-end extremity the excess equals a penalised e >= 0."""
        if not self.linear_end_nodes:
            return super().define_constraints_on_nodes()
        linear_ends = set(self.linear_end_nodes)
        for clone_id in self.clone_ids:
            for node, data in self.iag.nodes(data=True):
                for haplotype in HAPLOTYPES:
                    lin_exp = self.get_lin_expr_for_node_haplotype(clone_id=clone_id, node=node, haplotype=haplotype)
                    if node in self.hapl_telomeres:
                        self.gm.addLConstr(lin_exp, g.GRB.GREATER_EQUAL, 0, name="cne_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))
                    elif node in linear_ends:
                        self.gm.addLConstr(lin_exp, g.GRB.EQUAL, self.linear_end_variables[clone_id, node, haplotype],
                                           name="cne-linear_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))
                    else:
                        self.gm.addLConstr(lin_exp, g.GRB.EQUAL, 0, name="cnb_{{{cid},{v},{hap}}}".format(cid=clone_id, v=node, hap=haplotype))

    def define_constraints_nov_adjacency_overall_presence(self):
        """RCK's 1 - fp rule over the observed novel adjacencies only; imputed ones are not SV calls.
        Skipped without observed novel adjacencies: RCK would write 0 >= 1 - fp, which is
        infeasible (e.g. an amplicon whose junctions are all imputed)."""
        observed = [a for a in self.hapl_adjacencies if a.adjacency_type == AdjacencyType.NOVEL and a.stable_id_non_phased not in self.imputed_aids]
        if not observed:
            return
        if not self.imputed_aids:
            return super().define_constraints_nov_adjacency_overall_presence()
        fp_lin_expr = g.LinExpr()
        for adjacency in observed:
            fp_lin_expr.add(g.quicksum([self.variables[PROD][PP][adjacency.stable_id_non_phased][ph] for ph in PHASINGS]), mult=(1.0 / len(observed)))
        self.gm.addLConstr(fp_lin_expr, g.GRB.GREATER_EQUAL, g.LinExpr(1 - self.hapl_nov_adjacencies_fp), name="nov-adj-fp")

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
        if self.linear_end_variables and opt.linear_end_penalty:
            terms.append(opt.linear_end_penalty * g.quicksum(self.linear_end_variables.values()))
        return terms

    def get_linear_ends_from_model(self) -> Dict[tuple, int]:
        return {key: int(round(v.X)) for key, v in self.linear_end_variables.items() if v.X > 0.5}

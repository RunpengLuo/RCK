"""Translate an AA breakpoint graph back into a subgraph folder that ``rck-ecdna`` reads.

    python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate]
                                  [--background-cn auto|X] [--cn-profile FILE] [--sample NAME]

The inverse of ``python -m rck.ecdna.aa_export subgraph`` (results/graph_route/): it writes
``nodes.tsv``, ``edges.tsv``, ``paths.tsv``, ``scale.json``, ``cn_profile.tsv``,
``summary.txt`` and ``import.json`` (what was mapped, dropped and fitted), which
``rck-ecdna --subgraph-dir OUTDIR --inputs subgraph`` reads. Any AA graph works (AA, CoRAL,
the truth mixture); a graph written by ``aa_export subgraph`` also carries ``#rck_*`` comment
lines, from which the subgraph's node, segment, adjacency and path ids and orders are restored,
so ``read_subgraph`` gives back the original normalised graph (:func:`compare_subgraphs`).

* positions: AA - ``SHIFT`` (the exporter's shift). A breakpoint between two contiguous sequence
  edges is an EID ``chr:pos±`` (``pos`` = left segment end; the sign is the strand of a
  discordant edge there, else ``+``), with nodes ``<EID>_+`` at ``pos`` and ``<EID>_-`` at
  ``pos + 1``, as subgraph_construct.py names them;
* ``sequence`` -> ``reference`` edge (``WEIGHT`` = AverageCoverage, ``CN_TOTAL`` = CN);
  ``concordant`` -> ``extremity`` edge (``WEIGHT`` = reads); ``discordant`` -> ``variant``
  edge (``WEIGHT`` = ``NOV_U`` = reads, no haplotype label, so every phasing is allowed);
  ``source`` edges are dropped (CoRAL adds its own s/t ends) and listed in import.json;
* ``#rck_outside`` lines (joins to a region's outer flank) are always restored; imputed
  candidates follow ``--imputed``: ``comment`` restores the ``#rck_imputed`` lines, ``none``
  leaves them out (what CoRAL and CycleExtractor see), ``regenerate`` proposes RCK-ecDNA's own
  closing candidates on the graph (:func:`regenerate_imputed`);
* ``paths.tsv`` from the ``path_constraint`` lines (and ``#rck_path`` lines);
* ``scale.json`` is fitted from the graph itself (:func:`fit_scale`). ``rck-ecdna`` takes
  segment CN from ``CN_TOTAL`` but junction CN as reads / ``junction_reads_per_copy``, so the
  junction scale matters; a graph without read counts gets ``WEIGHT`` = CN and a scale of 1.
  A graph rewritten by CoRAL (``#CoRAL`` header) is not used for the depth scale: CoRAL's
  writer puts ``lr_nc / size`` in the AverageCoverage column.
"""

from __future__ import annotations

import argparse
import json
import logging
import math
import os
import shutil
import sys
from collections import Counter, OrderedDict
from typing import Dict, List, Optional, Sequence, Tuple

import numpy as np
import pandas as pd

from rck.ecdna import aa_format as aaf
from rck.ecdna.aa_export import SHIFT

log = logging.getLogger(__name__)

IMPORT_FORMAT = "rck-ecdna-aa-import/1"
FLANK_BP = 1_000_000  # rck.ecdna.scale.ScaleParams.flank_bp
EDGE_COLUMNS = ["NODE_1", "NODE_2", "EDGE_TYPE", "WEIGHT", "CHROM_1", "POS_1", "CHROM_2", "POS_2", "AID", "ORIENT_1", "ORIENT_2",
                "n_snps", "NOV_A", "NOV_B", "NOV_U", "seg_hap", "junction_hap", "CN_TOTAL", "CN_A", "CN_B", "CN"]


def _f(text) -> Optional[float]:
    if text is None or text in ("NA", "", "."):
        return None
    return float(text)


def _snap(value: float, ratios: Sequence[float], rel: float = 1e-9) -> float:
    """``value`` rounded to 12 significant digits when every ratio agrees with that to ``rel``:
    reads / (reads / theta) gives theta back only to the last bit, and a scale written as 6.5
    should come back as 6.5."""
    snapped = float("{v:.12g}".format(v=value))
    if all(abs(r - snapped) <= rel * abs(snapped) for r in ratios):
        return snapped
    return value


def weighted_median(values: Sequence[float], weights: Sequence[float]) -> float:
    order = np.argsort(values)
    v, w = np.asarray(values, float)[order], np.asarray(weights, float)[order]
    cum = np.cumsum(w)
    return float(v[np.searchsorted(cum, cum[-1] / 2.0)])


# ------------------------------------------------------------------------------------------ #
# scale
# ------------------------------------------------------------------------------------------ #

def fit_scale(graph: aaf.AAGraph, coral_rewritten: bool, depth_per_copy: Optional[float] = None) -> Dict[str, object]:
    """Depth and junction scales fitted from the graph's own columns.

    ``junction_reads_per_copy``: median reads / CN over discordant edges with reads and CN
    (concordant edges when no discordant edge has reads); no read counts at all -> 1.0, and the
    edges get ``WEIGHT`` = CN (``reads_mode`` ``cn``). ``depth_per_copy``: length-weighted median
    AverageCoverage / CN over sequence edges with both; not fitted from a CoRAL-rewritten graph."""
    out: Dict[str, object] = {}
    edges = [e for e in graph.edges if e.kind in ("concordant", "discordant")]
    have_reads = any((e.reads or 0) > 0 for e in edges)
    if have_reads:
        ratios = [e.reads / e.cn for e in graph.discordant if (e.reads or 0) > 0 and e.cn > 0]
        used = "discordant"
        if not ratios:
            ratios = [e.reads / e.cn for e in graph.concordant if (e.reads or 0) > 0 and e.cn > 0]
            used = "concordant"
        theta = _snap(float(np.median(ratios)), ratios)
        out.update(junction_reads_per_copy=theta, junction_reads_method="median reads / CN over {n} {k} edges".format(n=len(ratios), k=used),
                   reads_mode="reads")
    else:
        out.update(junction_reads_per_copy=1.0, junction_reads_method="fallback: no reads (edge WEIGHT = CN)", reads_mode="cn")
    if depth_per_copy is not None:
        out.update(depth_per_copy=float(depth_per_copy), depth_method="override")
    elif coral_rewritten:
        out.update(depth_per_copy=1.0, depth_method="not fitted: CoRAL-rewritten graph (its AverageCoverage column is lr_nc / size)")
    else:
        pairs = [(s.coverage / s.cn, s.end - s.start + 1) for s in graph.sequences if s.coverage and s.coverage > 0 and s.cn > 0]
        if pairs:
            ratios = [p for p, _ in pairs]
            out.update(depth_per_copy=_snap(weighted_median(ratios, [w for _, w in pairs]), ratios),
                       depth_method="length-weighted median AverageCoverage / CN over {n} sequence edges".format(n=len(pairs)))
        else:
            out.update(depth_per_copy=1.0, depth_method="fallback: no coverage")
    return out


# ------------------------------------------------------------------------------------------ #
# regenerated imputed candidates
# ------------------------------------------------------------------------------------------ #

def amplified_threshold(background: float, profile: Optional[pd.DataFrame], chrom: str, cn_mult: float = 1.5, mad_mult: float = 1.5) -> float:
    """subgraph_construct.py's amplification threshold, ``min(median + MAD_MULT * MAD,
    median * CN_MULT)``, over the chromosome's copy-number profile bins (the window depths it
    is computed on there), else ``background * CN_MULT``."""
    if profile is not None:
        rows = profile[profile["#CHR"].astype(str) == chrom]
        values = rows["CN_TOTAL"].astype(float).values
        values = values[np.isfinite(values)]
        if len(values):
            med = float(np.median(values))
            mad = 1.4826 * float(np.median(np.abs(values - med)))
            return min(med + mad_mult * mad, med * cn_mult)
    return background * cn_mult


def regenerate_imputed(sequences: Sequence[Tuple[str, int, int, float, float, str, str]], background: float,
                       profile: Optional[pd.DataFrame] = None, boundary_min_copies: float = 2.0, boundary_frac: float = 0.25,
                       cn_mult: float = 1.5, mad_mult: float = 1.5) -> List[Dict[str, object]]:
    """RCK-ecDNA's own closing candidates on an AA graph, by subgraph_construct.py's rule
    (``build_full_graph``): every amplified depth segment gets an ``intra`` edge from its end
    (``+``) back to its start (``-``), and consecutive amplified segments an ``inter`` edge.

    subgraph_construct.py finds amplified segments by segmenting window depth; on an AA graph
    the analogue is a maximal run of contiguous sequence edges above the amplification threshold
    (:func:`amplified_threshold`) without a copy-number step inside, a step being rck-ecdna's
    ``--boundary-min-copies`` / ``--boundary-frac`` rule. Its third kind, ``anchor`` (an
    imputed depth breakpoint tied to the nearest observed one), is not made: here every segment
    end is already a breakpoint, so an anchor would be the reference join that read_subgraph
    merges it into. ``inter`` edges between adjacent segments are merged the same way.

    ``sequences``: (chrom, start, end, cn, coverage, tail nid, head nid), sorted."""
    out: List[Dict[str, object]] = []
    groups: List[List[Tuple]] = []
    thresholds: Dict[str, float] = {}
    for s in sequences:
        chrom, start, end, cn = s[0], s[1], s[2], s[3]
        thr = thresholds.setdefault(chrom, amplified_threshold(background, profile, chrom, cn_mult, mad_mult))
        if cn <= thr:
            groups.append([])
            continue
        if groups and groups[-1]:
            last = groups[-1][-1]
            step = abs(cn - last[3]) >= max(boundary_min_copies, boundary_frac * max(cn, last[3]))
            if last[0] == chrom and start == last[2] + 1 and not step:
                groups[-1].append(s)
                continue
        groups.append([s])
    segs = [g for g in groups if g]
    k = 0
    for i, g in enumerate(segs):
        chrom = g[0][0]
        length = sum(x[2] - x[1] + 1 for x in g)
        depth = sum((x[4] if x[4] is not None else x[3]) * (x[2] - x[1] + 1) for x in g) / max(length, 1)
        k += 1
        out.append({"aid": "imputed_adj_intra_{c}_r{k}".format(c=chrom, k=k), "u": g[-1][6], "v": g[0][5], "support": depth, "kind": "intra"})
        if i + 1 < len(segs) and segs[i + 1][0][0] == chrom:
            nxt = segs[i + 1]
            k += 1
            nlen = sum(x[2] - x[1] + 1 for x in nxt)
            ndepth = sum((x[4] if x[4] is not None else x[3]) * (x[2] - x[1] + 1) for x in nxt) / max(nlen, 1)
            out.append({"aid": "imputed_adj_inter_{c}_r{k}".format(c=chrom, k=k), "u": g[-1][6], "v": nxt[0][5],
                        "support": (depth + ndepth) / 2.0, "kind": "inter"})
    return out


# ------------------------------------------------------------------------------------------ #
# import
# ------------------------------------------------------------------------------------------ #

class _Nodes:
    """Node ids of the imported graph, by coordinate (subgraph frame) and id."""

    def __init__(self):
        self.rows: "OrderedDict[str, Dict[str, object]]" = OrderedDict()
        self.at: Dict[Tuple[str, int, str], str] = {}

    def add(self, nid: str, chrom: str, pos: int, side: Optional[str], source: str) -> str:
        if nid not in self.rows:
            self.rows[nid] = {"NODE_ID": nid, "EID": nid.rpartition("_")[0] if nid[-2:] in ("_+", "_-") else nid, "#CHR": chrom,
                              "POS": int(pos), "ORIENT": side if side in ("+", "-") else "None", "SOURCE": source}
            self.at.setdefault((chrom, int(pos), side or "+"), nid)
        return nid


def _eid(chrom: str, q: int, sign: str) -> str:
    return "{c}:{q}{s}".format(c=chrom, q=q, s=sign)


def import_graph(graph_path: str, outdir: str, imputed: str = "comment", background_cn: object = "auto", cn_profile: Optional[str] = None,
                 sample: Optional[str] = None, depth_per_copy: Optional[float] = None, name: Optional[str] = None,
                 boundary_min_copies: float = 2.0, boundary_frac: float = 0.25) -> Dict[str, object]:
    """Write the subgraph folder ``outdir`` for the AA graph ``graph_path`` (module docstring)."""
    if imputed not in ("comment", "none", "regenerate"):
        raise ValueError("imputed must be comment, none or regenerate")
    g = aaf.read_graph(graph_path)
    os.makedirs(outdir, exist_ok=True)
    name = name or os.path.basename(graph_path).replace("_graph.txt", "")
    coral_rewritten = any(c.startswith("#CoRAL") for c in g.comments)
    rck_node, rck_seq, rck_edge = g.rck_lines("rck_node"), g.rck_lines("rck_seq"), g.rck_lines("rck_edge")
    ids = bool(rck_node) and len(rck_seq) == len(g.sequences)
    info: Dict[str, object] = {"format": IMPORT_FORMAT, "created_by": "rck.ecdna.aa_import", "graph": os.path.abspath(graph_path),
                               "name": name, "shift": SHIFT, "imputed": imputed, "rck_ids": ids, "coral_rewritten": coral_rewritten,
                               "warnings": []}
    warn = info["warnings"].append

    def sub(b: aaf.Breakend) -> Tuple[str, int, str]:
        return b[0], b[1] - SHIFT, b[2]

    nodes = _Nodes()
    seqs = []  # (chrom, start, end, cn, coverage, tail, head)
    if ids:
        for f in rck_node:
            nodes.add(f[0], f[1], int(f[2]) - SHIFT, None if f[3] == "." else f[3], "" if f[4] == "." else f[4])
        by_index = {int(f[0]): f for f in rck_seq}
        for i, s in enumerate(g.sequences, start=1):
            _k, _sid, tail, head = by_index[i]
            seqs.append((s.chrom, s.start - SHIFT, s.end - SHIFT, s.cn, s.coverage, tail, head))
        seq_order = sorted(range(len(seqs)), key=lambda i: by_index[i + 1][1])  # the subgraph's sid order
    else:
        sign: Dict[Tuple[str, int], str] = {}
        for e in sorted(g.discordant, key=lambda e: (e.a, e.b)):
            for c, p, st in (sub(e.a), sub(e.b)):
                sign.setdefault((c, p if st == "+" else p - 1), st)
        for s in g.sequences:
            c, a, b = s.chrom, s.start - SHIFT, s.end - SHIFT
            tail = nodes.add(_eid(c, a - 1, sign.get((c, a - 1), "+")) + "_-", c, a, "-", "extremity")
            head = nodes.add(_eid(c, b, sign.get((c, b), "+")) + "_+", c, b, "+", "extremity")
            seqs.append((c, a, b, s.cn, s.coverage, tail, head))
        seq_order = list(range(len(seqs)))

    def node_at(b: aaf.Breakend, create: bool = True) -> Optional[str]:
        c, p, st = sub(b)
        nid = nodes.at.get((c, p, st))
        if nid is None and create:  # a node off the sequence edges: a region flank
            q = p if st == "+" else p - 1
            nid = nodes.add(_eid(c, q, st) + "_" + st, c, p, st, "extremity")
        return nid

    scale = fit_scale(g, coral_rewritten, depth_per_copy)
    cn_reads = scale["reads_mode"] == "cn"
    theta = float(scale["junction_reads_per_copy"])

    def weight(reads, cn):
        if cn_reads:
            return float(cn) * theta
        return 0.0 if reads is None else float(reads)

    rows: List[Tuple[float, Dict[str, object]]] = []  # (order, adjacency row)

    def row(kind, u, v, w, aid=None, cn=None):
        nu, nv = nodes.rows[u], nodes.rows[v]
        r = {"NODE_1": u, "NODE_2": v, "EDGE_TYPE": kind, "WEIGHT": w, "CHROM_1": nu["#CHR"], "POS_1": nu["POS"], "CHROM_2": nv["#CHR"],
             "POS_2": nv["POS"], "AID": aid, "ORIENT_1": nu["ORIENT"], "ORIENT_2": nv["ORIENT"], "CN": cn}
        if kind == "variant":
            r.update(NOV_A=0.0, NOV_B=0.0, NOV_U=w)
        return r

    edge_ids = {(f[0], int(f[1])): f for f in rck_edge}
    counts = Counter()
    for kind in ("concordant", "discordant"):
        for i, e in enumerate((x for x in g.edges if x.kind == kind), start=1):
            meta = edge_ids.get((kind, i)) if ids else None
            if meta is not None:
                aid, u, v, reads, order = meta[2], meta[3], meta[4], _f(meta[5]), float(meta[6])
            else:
                aid, u, v, reads, order = None, node_at(e.a), node_at(e.b), e.reads, len(rows) + 1e6
                if kind == "discordant":
                    aid = "aa_d{k}".format(k=i)
            if kind == "concordant":
                lo, hi = aaf.lower_first(sub(e.a), sub(e.b))
                if not (lo[0] == hi[0] and lo[2] == "+" and hi[2] == "-" and hi[1] == lo[1] + 1):
                    raise ValueError("{p}: concordant edge {a}->{b} does not join neighbouring segments".format(
                        p=graph_path, a=aaf.format_breakend(e.a), b=aaf.format_breakend(e.b)))
                if reads is None and not (u.startswith("imputed_ext") or v.startswith("imputed_ext")) and not cn_reads:
                    warn("unobserved reference join {a} written with WEIGHT 0 (edges.tsv reads it as observed)".format(a=aid))
                rows.append((order, row("extremity", u, v, weight(reads, e.cn), aid, e.cn)))
            else:
                rows.append((order, row("variant", u, v, weight(reads, e.cn), aid, e.cn)))
            counts[kind] += 1
    info["source_edges_dropped"] = [aaf.format_breakend(e.a) + "->" + aaf.format_breakend(e.b) for e in g.edges if e.kind == "source"]

    for f in g.rck_lines("rck_outside"):
        a, b = (aaf.parse_breakend(x) for x in f[0].split("->"))
        cn, reads, kind = _f(f[1]), _f(f[2]), f[3]
        if len(f) >= 8:
            aid, u, v, order = f[4], f[5], f[6], float(f[7])
        else:
            aid, u, v, order = None, node_at(a), node_at(b), len(rows) + 1e6
        etype = "extremity" if kind == "reference" else "variant"
        w = weight(reads, cn) if cn is not None or not cn_reads else 0.0
        rows.append((order, row(etype, u, v, w, aid, cn)))
        counts["outside"] += 1

    profile = pd.read_csv(cn_profile, sep="\t") if cn_profile else None
    if background_cn in (None, "auto"):
        if profile is not None:
            lengths = (profile["END"] - profile["START"]).clip(lower=1).values
            background, bg_method = weighted_median(profile["CN_TOTAL"].astype(float).values, lengths), "auto: length-weighted median of --cn-profile"
        else:
            background, bg_method = 2.0, "auto: 2 (no --cn-profile)"
    else:
        background, bg_method = float(background_cn), "given"

    rck_imputed = g.rck_lines("rck_imputed")
    info["imputed_in_graph"] = len(rck_imputed)
    if imputed == "comment":
        for k, f in enumerate(rck_imputed, start=1):
            a, b = (aaf.parse_breakend(x) for x in f[0].split("->"))
            support = _f(f[1]) if len(f) > 1 else None
            if len(f) >= 6:
                aid, u, v, order = f[2], f[3], f[4], float(f[5])
            else:
                aid, u, v, order = "aa_imputed{k}".format(k=k), node_at(a), node_at(b), len(rows) + 1e6
            rows.append((order, row("imputed_variant", u, v, 0.0 if support is None else support, aid)))
            counts["imputed"] += 1
    elif imputed == "regenerate":
        for k, cand in enumerate(regenerate_imputed(seqs, background, profile, boundary_min_copies, boundary_frac)):
            rows.append((len(rows) + 2e6, row("imputed_variant", cand["u"], cand["v"], cand["support"], cand["aid"])))
            counts["imputed"] += 1
            counts["regenerated_" + cand["kind"]] += 1

    # ---- tables ---- #
    ref_rows = []
    for i in seq_order:
        c, a, b, cn, cov, tail, head = seqs[i]
        w = cov if (cov is not None and not coral_rewritten) else cn * float(scale["depth_per_copy"])
        ref_rows.append({"NODE_1": tail, "NODE_2": head, "EDGE_TYPE": "reference", "WEIGHT": w, "CHROM_1": c, "POS_1": a, "CHROM_2": c,
                         "POS_2": b, "n_snps": 0, "seg_hap": "unknown", "CN_TOTAL": cn})
    adj_rows = [r for _, r in sorted(rows, key=lambda x: x[0])]
    pd.DataFrame(ref_rows + adj_rows, columns=EDGE_COLUMNS).to_csv(os.path.join(outdir, "edges.tsv"), sep="\t", index=False)
    pd.DataFrame(list(nodes.rows.values()), columns=["NODE_ID", "EID", "#CHR", "POS", "ORIENT", "SOURCE"]).to_csv(
        os.path.join(outdir, "nodes.tsv"), sep="\t", index=False)

    # ---- paths ---- #
    edge_nodes: Dict[str, List[Tuple[str, str]]] = {"e": [(s[5], s[6]) for s in seqs]}
    conc_list, disc_list = [], []  # breakpoint edges in file order, as path tokens index them
    for kind, target in (("concordant", conc_list), ("discordant", disc_list)):
        for i, e in enumerate((x for x in g.edges if x.kind == kind), start=1):
            meta = edge_ids.get((kind, i)) if ids else None
            target.append((meta[3], meta[4]) if meta is not None else (node_at(e.a), node_at(e.b)))
    edge_nodes.update(c=conc_list, d=disc_list)
    pc_meta = {int(f[0]): f for f in g.rck_lines("rck_pc")}
    path_rows: List[Tuple[float, Dict[str, object]]] = []
    pstats = Counter()
    for i, (tokens, support) in enumerate(g.path_constraints, start=1):
        steps = _path_steps(aaf.parse_path_tokens(tokens), edge_nodes)
        if steps is None:
            pstats["unreadable"] += 1
            warn("path constraint {i} ({t}) does not follow the graph's edges; skipped".format(i=i, t=tokens))
            continue
        meta = pc_meta.get(i)
        pid, order, nov = (meta[2], float(meta[1]), float(meta[3])) if meta else ("pc{i:05d}".format(i=i), 1e6 + i, float(support))
        if meta and len(meta) > 4 and meta[4] == "1":  # CoRAL's token order reversed the read; restore its direction
            steps = steps[::-1]
        path_rows.append((order, {"PID": pid, "PATH": ",".join(_step(n) for n in steps), "NSTEP": len(steps), "#NOV": nov}))
        pstats["path_constraint"] += 1
    for f in g.rck_lines("rck_path"):
        reason, order, pid, path, nstep, nov = f[0], float(f[1]), f[2], f[3], f[4], float(f[5])
        if reason == "imputed" and imputed != "comment":
            pstats["rck_path_imputed_dropped"] += 1
            continue
        path_rows.append((order, {"PID": pid, "PATH": path, "NSTEP": int(float(nstep)) if nstep else len(path.split(",")), "#NOV": nov}))
        pstats["rck_path_" + reason] += 1
    pd.DataFrame([r for _, r in sorted(path_rows, key=lambda x: x[0])], columns=["PID", "PATH", "NSTEP", "#NOV"]).to_csv(
        os.path.join(outdir, "paths.tsv"), sep="\t", index=False)

    # ---- scale, profile, summary ---- #
    runs = []
    for c, a, b, *_r in seqs:
        if runs and runs[-1][0] == c and a == runs[-1][2] + 1:
            runs[-1] = (c, runs[-1][1], b)
        else:
            runs.append((c, a, b))
    if cn_profile:
        shutil.copyfile(cn_profile, os.path.join(outdir, "cn_profile.tsv"))
        profile_method = "copied from " + os.path.abspath(cn_profile)
    else:
        _write_profile(os.path.join(outdir, "cn_profile.tsv"), seqs, runs, background)
        profile_method = "sequence edges at their CN, {f} bp flanks at the background CN".format(f=FLANK_BP)
    scale_json = {"format": "rck-ecdna-subgraph-scale/1", "created_by": "rck.ecdna.aa_import (fitted from {p})".format(p=os.path.basename(graph_path)),
                  "sample": sample or name, "depth_per_copy": scale["depth_per_copy"], "depth_method": scale["depth_method"],
                  "junction_reads_per_copy": theta, "junction_reads_method": scale["junction_reads_method"],
                  "cell_ploidy": background, "background_cn": background, "background_method": bg_method,
                  "method": "rck.ecdna.aa_import", "files": {"cn_profile": "cn_profile.tsv"}}
    with open(os.path.join(outdir, "scale.json"), "w") as handle:
        json.dump(scale_json, handle, indent=2)
    with open(os.path.join(outdir, "summary.txt"), "w") as handle:
        handle.write("subgraph imported from AA graph {p} (rck.ecdna.aa_import, --imputed {i})\n\nCandidate regions:\n".format(
            p=os.path.abspath(graph_path), i=imputed))
        for c, a, b in runs:
            handle.write("  {c}:{a}-{b}\n".format(c=c, a=a, b=b))
    info.update(counts=dict(counts, sequence=len(seqs)), paths=dict(pstats), scale=scale, background_cn=background,
                background_method=bg_method, cn_profile=profile_method, regions=["{c}:{a}-{b}".format(c=c, a=a, b=b) for c, a, b in runs])
    with open(os.path.join(outdir, "import.json"), "w") as handle:
        json.dump(info, handle, indent=2)
    return info


def _step(nid: str) -> str:
    eid, _, orient = nid.rpartition("_")
    return "{e}@{o}".format(e=eid, o=orient)


def _path_steps(tokens, edge_nodes) -> Optional[List[str]]:
    """Node steps of a path constraint: entry and exit node of every sequence edge, so steps
    (2i, 2i+1) run inside one sequence edge and (2i+1, 2i+2) are the breakpoint edge between."""
    if len(tokens) < 3 or len(tokens) % 2 == 0:
        return None
    steps: List[str] = []
    for j, (kind, idx, strand, _count) in enumerate(tokens):
        if (j % 2 == 0) != (kind == "e"):
            return None
        table = edge_nodes.get(kind)
        if table is None or not 0 <= idx < len(table):
            return None
        if kind == "e":
            tail, head = table[idx]
            entry, out = (tail, head) if strand == "+" else (head, tail)
            if steps and {steps[-1], entry} != set(edge_nodes[tokens[j - 1][0]][tokens[j - 1][1]]):
                return None
            steps += [entry, out]
    return steps


def _write_profile(path: str, seqs, runs, background: float) -> None:
    rows = [{"#CHR": c, "START": a, "END": b + 1, "CN_TOTAL": cn, "N_BINS": max((b - a + 1) / 1000.0, 1e-3)} for c, a, b, cn, *_r in seqs]
    covered = sorted((c, a, b + 1) for c, a, b, *_r in seqs)
    for c, a, b in runs:
        for lo, hi in ((max(a - FLANK_BP, 0), a), (b + 1, b + 1 + FLANK_BP)):
            pieces = [(lo, hi)]
            for cc, s, e in covered:  # keep flanks off other runs' sequence edges
                if cc != c:
                    continue
                pieces = [q for x, y in pieces for q in ((x, min(y, s)), (max(x, e), y)) if q[1] > q[0]]
            rows += [{"#CHR": c, "START": x, "END": y, "CN_TOTAL": background, "N_BINS": (y - x) / 1000.0} for x, y in pieces]
    frame = pd.DataFrame(rows, columns=["#CHR", "START", "END", "CN_TOTAL", "N_BINS"]).drop_duplicates()
    frame.sort_values(["#CHR", "START"]).to_csv(path, sep="\t", index=False)


# ------------------------------------------------------------------------------------------ #
# round-trip comparison
# ------------------------------------------------------------------------------------------ #

def _graph_state(subgraph_dir: str, soft_junction_hap: bool = True, **scale) -> Dict[str, object]:
    from rck.ecdna.aa_export import subgraph_observations
    from rck.ecdna.subgraph_io import read_subgraph
    from rck.ecdna.walk import read_paths

    g = read_subgraph(subgraph_dir, soft_junction_hap=soft_junction_hap)
    obs, _sub = subgraph_observations(subgraph_dir, g, **scale)

    def xy(nid):
        n = g.nodes[nid]
        s = g.segment_at(nid)
        if s is None:
            return (n.chrom, n.pos, n.side)
        return (s.chrom, s.start, "-") if nid == s.tail else (s.chrom, s.end, "+")

    sids = sorted(g.segments)
    p = os.path.join(subgraph_dir, "paths.tsv")
    paths = read_paths(p, g) if os.path.exists(p) and os.path.getsize(p) > 0 else []
    return {
        "nodes": {nid: (xy(nid), n.terminal, n.segment is not None) for nid, n in g.nodes.items()},
        "node_order": list(g.nodes),
        "segments": [(g.segments[k].chrom, g.segments[k].start, g.segments[k].end, g.segments[k].tail, g.segments[k].head) for k in sids],
        "segment_insertion": [(g.segments[k].chrom, g.segments[k].start) for k in g.segments],
        "seg_total": [obs.segments[k].total for k in sids],
        "background": [obs.background[k] for k in sids],
        "adjacencies": {aid: (a.u, a.v, a.kind, a.reads, tuple(sorted(str(x) for x in a.allowed))) for aid, a in g.adjacencies.items()},
        "adj_order": list(g.adjacencies),
        "adj_obs": {aid: obs.adjacencies.get(aid) for aid in g.adjacencies},
        "paths": [([(xy(x), xy(y)) for x, y in r.junctions], r.support) for r in paths],
        "path_pairs": [{frozenset(j) for j in r.junctions} for r in paths],
        "regions": g.regions(),
        "theta": obs.scale.junction_reads_per_copy,
    }


def _close(a, b, rel=1e-12) -> bool:
    if a is None or b is None:
        return a is b
    return abs(a - b) <= rel * max(1.0, abs(a), abs(b))


def compare_subgraphs(original: str, imported: str, rel: float = 1e-12, ignore_imputed: bool = False,
                      original_scale: Optional[Dict[str, float]] = None) -> List[str]:
    """Differences between the normalised graphs ``rck-ecdna`` builds from two subgraph folders
    (total-only, junction labels ignored): nodes (ids, coordinates, terminals and order),
    segments, adjacencies (ids, ends, kind, reads, allowed phasings, order), the observations
    (segment CN, junction CN, flank background), read paths (as junctions in coordinates) and
    regions. ``ignore_imputed``: leave imputed candidates (and the paths through them) out of
    the comparison. An empty list means the round trip is exact."""
    from rck.ecdna.subgraph_io import IMPUTED

    a = _graph_state(original, **(original_scale or {}))
    b = _graph_state(imported)
    out: List[str] = []
    if ignore_imputed:
        drop = {k for k, v in a["adjacencies"].items() if v[2] == IMPUTED}
        for s in (a, b):
            s["adjacencies"] = {k: v for k, v in s["adjacencies"].items() if k not in drop}
            s["adj_order"] = [k for k in s["adj_order"] if k not in drop]
            s["adj_obs"] = {k: v for k, v in s["adj_obs"].items() if k not in drop}
            used = {n for v in s["adjacencies"].values() for n in v[:2]} | {n for seg in s["segments"] for n in seg[3:]}
            s["nodes"] = {k: v for k, v in s["nodes"].items() if k in used}
            s["node_order"] = [k for k in s["node_order"] if k in used]
    for key in ("nodes", "segments", "regions", "node_order", "adj_order", "segment_insertion"):
        if a[key] != b[key]:
            out.append("{k} differ".format(k=key))
    if set(a["adjacencies"]) != set(b["adjacencies"]):
        out.append("adjacency ids differ: only original {x}, only imported {y}".format(
            x=sorted(set(a["adjacencies"]) - set(b["adjacencies"])), y=sorted(set(b["adjacencies"]) - set(a["adjacencies"]))))
    for k in sorted(set(a["adjacencies"]) & set(b["adjacencies"])):
        x, y = a["adjacencies"][k], b["adjacencies"][k]
        if x[:3] != y[:3] or x[4] != y[4] or not _close(x[3], y[3], rel):
            out.append("adjacency {k}: {x} vs {y}".format(k=k, x=x, y=y))
        if not _close(a["adj_obs"][k], b["adj_obs"][k], rel):
            out.append("adjacency {k} observed CN {x} vs {y}".format(k=k, x=a["adj_obs"][k], y=b["adj_obs"][k]))
    for key in ("seg_total", "background"):
        if len(a[key]) != len(b[key]) or not all(_close(x, y, rel) for x, y in zip(a[key], b[key])):
            out.append("{k} differ".format(k=key))
    pa, pb = a["paths"], b["paths"]
    if ignore_imputed:  # paths through an imputed candidate are left out with it
        pairs = {frozenset(v[:2]) for k, v in _graph_state(original, **(original_scale or {}))["adjacencies"].items() if v[2] == IMPUTED}
        pa = [p for p, jp in zip(a["paths"], a["path_pairs"]) if not (jp & pairs)]
    if len(pa) != len(pb) or any(x[0] != y[0] or not _close(x[1], y[1], rel) for x, y in zip(pa, pb)):
        out.append("read paths differ ({n} vs {m})".format(n=len(pa), m=len(pb)))
    if a["theta"] != b["theta"]:
        out.append("junction reads per copy {x!r} vs {y!r}".format(x=a["theta"], y=b["theta"]))
    return out


def main(argv=None) -> int:
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("graph", help="an AA graph (*_graph.txt)")
    parser.add_argument("outdir", help="subgraph folder to write")
    parser.add_argument("--imputed", choices=["comment", "none", "regenerate"], default="comment",
                        help="imputed candidates: restore #rck_imputed lines (default), none (what CoRAL / CycleExtractor see), "
                             "or regenerate RCK-ecDNA's own closing candidates")
    parser.add_argument("--background-cn", default="auto", help="flank / baseline CN (default auto: median of --cn-profile, else 2)")
    parser.add_argument("--cn-profile", default=None, help="copy this cn_profile.tsv (round trip: the original subgraph's)")
    parser.add_argument("--sample", default=None, help="sample name for scale.json (default: the graph's name)")
    parser.add_argument("--depth-per-copy", type=float, default=None, help="override the fitted depth per copy")
    parser.add_argument("--name", default=None)
    parser.add_argument("--boundary-min-copies", type=float, default=2.0, help="copy-number step for --imputed regenerate (as rck-ecdna)")
    parser.add_argument("--boundary-frac", type=float, default=0.25, help="... or this fraction of the larger CN (as rck-ecdna)")
    parser.add_argument("--compare", default=None, metavar="SUBGRAPH_DIR",
                        help="after importing, compare with this subgraph folder (compare_subgraphs) and print the differences")
    args = parser.parse_args(argv)
    info = import_graph(args.graph, args.outdir, imputed=args.imputed, background_cn=args.background_cn, cn_profile=args.cn_profile,
                        sample=args.sample, depth_per_copy=args.depth_per_copy, name=args.name,
                        boundary_min_copies=args.boundary_min_copies, boundary_frac=args.boundary_frac)
    print(json.dumps({k: info[k] for k in ("counts", "paths", "scale", "background_cn", "warnings")}, indent=2))
    if args.compare:
        diffs = compare_subgraphs(args.compare, args.outdir, ignore_imputed=args.imputed != "comment")
        print("\n".join(diffs) if diffs else "round trip exact")
        return 1 if diffs else 0
    return 0


if __name__ == "__main__":
    sys.exit(main())

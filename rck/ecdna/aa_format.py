"""Read AmpliconArchitect (AA) graph and cycles files, from any method.

The AA formats (as ecSimulator's ``ecSim_IO`` writes them, and AA, CoRAL and ``rck-ecdna`` too):

``*_graph.txt``
    ``sequence  chr:start-  chr:end+  PredictedCopyCount  AverageCoverage  Size  NumberReadsMapped``
    and ``concordant`` / ``discordant`` / ``source`` lines
    ``kind  chr:pos±->chr:pos±  PredictedCopyCount  NumberOfReadPairs  ...``.
``*_cycles.txt``
    ``Interval  i  chrom  start  end``, ``Segment  id  chrom  start  end`` and
    ``Cycle=k;Copy_count=w;Segments=1+,2-,...`` lines. A path that is not a cycle starts and ends
    with segment ``0`` (AA's source: ``0+,...,0-``; ecSimulator writes ``0+,...,0+``).

Coordinates are 1-based and inclusive. Cycles become :class:`rck.ecdna.evaluate.Cycle` objects
with 0-based inclusive steps (``start - 1``, ``end - 1``), the convention under which
:func:`rck.ecdna.evaluate.read_aa_cycles` reads the simulator's truth, so truth and inferred
cycles read from AA files are compared in one frame.

Optional keys on a ``Cycle=`` line, which AA and CoRAL ignore:

* ``Haplotypes=A,B,...``: one haplotype per step (``rck-ecdna``; ``.`` for the source). Without
  it, the cycle is unphased and its steps carry haplotype ``.``;
* ``Status=low_weight``: the method itself does not call this cycle a species (``rck-ecdna``);
* ``Alt_segments=...|...`` and ``Alt_haplotypes=...|...``: other orderings of the same species
  (``rck-ecdna`` reports every distinct Eulerian ordering of a species' edges);
* ``Species=name``: a name for the cycle (the truth mixture's cycles files).

:func:`write_graph` writes a graph that CoRAL 3.1.1 (``coral cycle``) and CycleExtractor 1.0.0
both read (results/graph_route/NEXT_STEPS.md lists the rules taken from their parsers):
tab-separated, no blank line, 7-field ``sequence`` lines, 6-field breakpoint lines,
``path_constraint`` lines with 1-based ``e``/``c``/``d`` edge indices in file order, and
``interval`` lines that both tools ignore. :func:`read_graph` reads those back, together with
every ``#`` comment line (``AAGraph.comments``); ``#rck_*`` comment lines carry what the
AA format cannot hold (see :mod:`rck.ecdna.aa_export` and :mod:`rck.ecdna.aa_import`).
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Dict, List, Optional, Sequence, Tuple

from rck.ecdna.evaluate import CStep, Cycle

UNPHASED = "."

Breakend = Tuple[str, int, str]  # (chrom, AA position, strand)


@dataclass
class AAEdge:
    kind: str  # concordant, discordant or source
    a: Breakend
    b: Breakend
    cn: float
    reads: Optional[float] = None


@dataclass
class AASequence:
    chrom: str
    start: int  # 1-based inclusive
    end: int
    cn: float
    coverage: Optional[float] = None  # AverageCoverage (CoRAL's writer puts lr_nc / size here)
    size: Optional[int] = None
    reads: Optional[float] = None  # NumberReadsMapped / NumberOfLongReads


@dataclass
class AAGraph:
    path: str
    sequences: List[AASequence] = field(default_factory=list)
    edges: List[AAEdge] = field(default_factory=list)
    path_constraints: List[Tuple[str, str]] = field(default_factory=list)  # (tokens, support) as written
    intervals: List[Tuple[str, int, int]] = field(default_factory=list)  # ``interval`` lines (ignored by CoRAL 3.1.1)
    comments: List[str] = field(default_factory=list)  # every line starting with "#", without the newline

    @property
    def discordant(self) -> List[AAEdge]:
        return [e for e in self.edges if e.kind == "discordant"]

    @property
    def concordant(self) -> List[AAEdge]:
        return [e for e in self.edges if e.kind == "concordant"]

    def rck_lines(self, tag: str) -> List[List[str]]:
        """Tab fields of the ``#<tag>`` comment lines (e.g. ``rck_imputed``), the tag left out."""
        out = []
        for line in self.comments:
            f = line.split("\t")
            if f[0] == "#" + tag:
                out.append(f[1:])
        return out


@dataclass
class AACycleRecord:
    cycle_id: str
    copy_count: float
    tokens: List[Tuple[str, str]]  # (segment id, strand), the source "0" left out
    linear: bool
    haplotypes: Optional[List[str]]
    fields: Dict[str, str]
    alternatives: List[Tuple[List[Tuple[str, str]], Optional[List[str]]]] = field(default_factory=list)


@dataclass
class AACycles:
    path: str
    intervals: List[Tuple[str, int, int]]  # 1-based inclusive
    segments: Dict[str, Tuple[str, int, int]]
    cycles: List[AACycleRecord]

    @property
    def phased(self) -> bool:
        return bool(self.cycles) and all(c.haplotypes is not None for c in self.cycles)


def parse_breakend(text: str) -> Breakend:
    chrom, rest = text.strip().rsplit(":", 1)
    return chrom, int(rest[:-1]), rest[-1]


def _float(text: str) -> Optional[float]:
    try:
        return float(text)
    except (TypeError, ValueError):
        return None


def read_graph(path: str) -> AAGraph:
    graph = AAGraph(path)
    with open(path) as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith("#"):
                graph.comments.append(line)
                continue
            f = line.split("\t")
            if f[0] == "sequence" and len(f) >= 4:
                a, b = parse_breakend(f[1]), parse_breakend(f[2])
                size = _float(f[5]) if len(f) > 5 else None
                graph.sequences.append(AASequence(a[0], a[1], b[1], float(f[3]), coverage=_float(f[4]) if len(f) > 4 else None,
                                                  size=int(size) if size is not None else None,
                                                  reads=_float(f[6]) if len(f) > 6 else None))
            elif f[0] == "path_constraint" and len(f) >= 3:
                graph.path_constraints.append((f[1], f[2]))
            elif f[0] == "interval" and len(f) >= 4:
                graph.intervals.append((f[1], int(f[2]), int(f[3])))
            elif f[0] in ("concordant", "discordant", "source") and len(f) >= 3:
                left, right = f[1].split("->")
                reads = None
                if len(f) > 3:
                    try:
                        reads = float(f[3])
                    except ValueError:
                        reads = None
                graph.edges.append(AAEdge(f[0], parse_breakend(left), parse_breakend(right), float(f[2]), reads))
    return graph


def chrom_key(name: str) -> Tuple[int, int, str]:
    """CoRAL 3.1.1's ``core_types.chr_sort_key`` (chr1 < chr2 < ... < chrX < chrY < other contigs < chrM).

    ``coral cycle`` sorts sequence and concordant edges with it after parsing
    (``infer_amplicon_intervals`` -> ``sort_edges``) without renumbering the path constraints,
    which index edges in file order; a graph with path constraints must already be in this order."""
    stem = name[3:] if name[:3].lower() == "chr" else name
    if stem.upper() in ("M", "MT"):
        return (2, 0, "")
    if stem.upper() == "X":
        return (1, 0, "")
    if stem.upper() == "Y":
        return (1, 1, "")
    try:
        return (0, int(stem), "")
    except ValueError:
        return (1, 2, stem)


def sequence_key(s: "AASequence"):
    """CoRAL's ``SequenceEdge.__lt__`` order."""
    return chrom_key(s.chrom), s.start, s.end


def edge_key(a: Breakend, b: Breakend):
    """CoRAL's ``BreakpointEdge.__lt__`` order (node1 then node2, natural chromosome order)."""
    return chrom_key(a[0]), a[1], chrom_key(b[0]), b[1]


def lower_first(a: Breakend, b: Breakend) -> Tuple[Breakend, Breakend]:
    """The breakend pair with the lower-coordinate end first (CoRAL's writer convention)."""
    return (a, b) if (chrom_key(a[0]), a[1], a[2]) <= (chrom_key(b[0]), b[1], b[2]) else (b, a)


def format_breakend(b: Breakend) -> str:
    return "{c}:{p}{s}".format(c=b[0], p=b[1], s=b[2])


SEQUENCE_HEADER = "SequenceEdge: StartPosition, EndPosition, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped"
BREAKPOINT_HEADER = ("BreakpointEdge: StartPosition->EndPosition, PredictedCopyCount, NumberOfReadPairs, "
                     "HomologySizeIfAvailable(<0ForInsertions), Homology/InsertionSequence")
PATH_HEADER = "PathConstraint: Path, Support"
INTERVAL_HEADER = "AmpliconIntervals: chr, start, end"


def _num(x: float) -> str:
    return repr(float(x))


def write_graph(path: str, sequences: Sequence[AASequence], edges: Sequence[AAEdge], path_constraints: Sequence[Tuple[str, object]] = (),
                intervals: Sequence[Tuple[str, int, int]] = (), header: Optional[str] = None, comments: Sequence[str] = ()) -> None:
    """Write an AA graph in the order given (callers sort; path-constraint indices refer to it).

    Layout: ``header`` (a ``#`` line), ``SequenceEdge:`` + 7-field ``sequence`` lines,
    ``BreakpointEdge:`` + 6-field ``concordant`` / ``discordant`` lines (``CN reads None None``),
    ``PathConstraint:`` + ``path_constraint`` lines, ``AmpliconIntervals:`` + ``interval`` lines,
    then ``comments``. Floats are written with ``repr``, so they read back exactly; no blank
    line (CoRAL 3.1.1's parser crashes on one)."""
    lines = []
    if header is not None:
        if not header.startswith("#"):
            raise ValueError("the header must be a # comment line")
        lines.append(header)
    lines.append(SEQUENCE_HEADER)
    for s in sequences:
        size = s.size if s.size is not None else s.end - s.start + 1
        lines.append("\t".join(["sequence", "{c}:{p}-".format(c=s.chrom, p=s.start), "{c}:{p}+".format(c=s.chrom, p=s.end), _num(s.cn),
                                _num(s.coverage if s.coverage is not None else s.cn), str(int(size)), str(int(s.reads or 0))]))
    lines.append(BREAKPOINT_HEADER)
    for e in edges:
        lines.append("\t".join([e.kind, format_breakend(e.a) + "->" + format_breakend(e.b), _num(e.cn),
                                _num(e.reads if e.reads is not None else 0.0), "None", "None"]))
    if path_constraints:
        lines.append(PATH_HEADER)
        for tokens, support in path_constraints:
            lines.append("\t".join(["path_constraint", tokens, str(support)]))
    if intervals:
        lines.append(INTERVAL_HEADER)
        for chrom, start, end in intervals:
            lines.append("\t".join(["interval", chrom, str(int(start)), str(int(end))]))
    for c in comments:
        if not c.startswith("#") or "\n" in c:
            raise ValueError("comment lines must start with # and hold no newline: {c!r}".format(c=c))
        lines.append(c)
    if any(not line.strip() for line in lines):
        raise ValueError("refusing to write a blank line (CoRAL 3.1.1 cannot parse it)")
    with open(path, "w") as handle:
        handle.write("\n".join(lines) + "\n")


def coral_path_string(path: Sequence[object]) -> str:
    """CoRAL 3.1.1's ``core_utils.path_to_str`` over a plain path, so the tokens match CoRAL's own
    writer without importing CoRAL (checked against it by results/graph_route/check_coral_format.py).

    ``path`` alternates edges and node strands and starts and ends with a sequence edge:
    ``[(type, idx), strand, (type, idx), strand, ..., (type, idx)]``, type ``e`` / ``c`` / ``d``,
    ``idx`` 0-based in file order, ``strand`` the strand of the node between two edges
    (``+`` = the right end of a sequence edge, ``-`` its left end). Each token's count is the
    number of times its edge occurs in the path (``path_to_edge_count``)."""
    if len(path) < 3 or len(path) % 2 == 0:
        raise ValueError("a path constraint needs at least two sequence edges and one breakpoint edge")
    if path[0][1] > path[-1][1]:  # CoRAL: reverse if the first edge's index is larger than the last's
        path = list(path)[::-1]
    counts: Dict[Tuple[str, int], int] = {}
    for e in path[0::2]:
        counts[e] = counts.get(e, 0) + 1
    out = []
    for i in range(0, len(path), 2):
        kind, idx = path[i]
        if i < len(path) - 2:
            strand = path[i + 1]
        else:
            strand = "+" if path[i - 1] == "-" else "-"
        out.append("{t}{n}{s}:{k}".format(t=kind, n=idx + 1, s=strand, k=counts[(kind, idx)]))
    return ",".join(out)


def parse_path_tokens(tokens: str) -> List[Tuple[str, int, str, int]]:
    """``e2+:1,c2-:1,...`` -> [(type, 0-based index, strand, count)]."""
    out = []
    for t in tokens.split(","):
        t = t.strip()
        if not t:
            continue
        head, _, count = t.partition(":")
        out.append((head[0], int(head[1:-1]) - 1, head[-1], int(count) if count else 1))
    return out


def _tokens(text: str) -> Tuple[List[Tuple[str, str]], bool, List[int]]:
    """(steps, linear, positions of the kept steps among all tokens) of a ``Segments=`` value."""
    raw = [t.strip() for t in text.split(",") if t.strip()]
    steps, kept, linear = [], [], False
    for i, t in enumerate(raw):
        sid, strand = t[:-1], t[-1]
        if strand not in "+-":
            raise ValueError("cycle step {t} has no strand".format(t=t))
        if sid == "0":
            linear = True
            continue
        steps.append((sid, strand))
        kept.append(i)
    return steps, linear, kept


def _haplotypes(text: Optional[str], kept: Sequence[int]) -> Optional[List[str]]:
    if text is None:
        return None
    raw = [h.strip() for h in text.split(",")]
    return [raw[i] if i < len(raw) else UNPHASED for i in kept]


def read_cycles(path: str) -> AACycles:
    intervals, segments, cycles = [], {}, []
    with open(path) as handle:
        for line in handle:
            f = line.rstrip("\n").split("\t")
            if f[0] == "Interval" and len(f) >= 5:
                intervals.append((f[2], int(f[3]), int(f[4])))
            elif f[0] == "Segment" and len(f) >= 5:
                segments[f[1]] = (f[2], int(f[3]), int(f[4]))
            elif line.startswith("Cycle="):
                parts = dict(p.split("=", 1) for p in line.strip().split(";") if "=" in p)
                steps, linear, kept = _tokens(parts.get("Segments", ""))
                haps = _haplotypes(parts.get("Haplotypes"), kept)
                record = AACycleRecord(cycle_id=parts["Cycle"], copy_count=float(parts.get("Copy_count", "nan")), tokens=steps,
                                       linear=linear, haplotypes=haps, fields=parts)
                if parts.get("Alt_segments"):
                    alt_haps = parts.get("Alt_haplotypes", "").split("|") if parts.get("Alt_haplotypes") else []
                    for k, text in enumerate(parts["Alt_segments"].split("|")):
                        alt_steps, _, alt_kept = _tokens(text)
                        record.alternatives.append((alt_steps, _haplotypes(alt_haps[k], alt_kept) if k < len(alt_haps) else None))
                for sid, _ in steps:
                    if sid not in segments:
                        raise ValueError("{p}: cycle {c} uses segment {s}, which is not listed".format(p=path, c=record.cycle_id, s=sid))
                cycles.append(record)
    return AACycles(path, intervals, segments, cycles)


def _steps(cycles: AACycles, tokens, haps) -> List[CStep]:
    out = []
    for i, (sid, strand) in enumerate(tokens):
        chrom, start, end = cycles.segments[sid]
        out.append(CStep(chrom, start - 1, end - 1, haps[i] if haps is not None else UNPHASED, strand))
    return out


def to_cycles(cycles: AACycles, prefix: str = "") -> List[Cycle]:
    """Every cycle (and path) of a cycles file as a scoring :class:`Cycle`, named ``<prefix>``
    + ``Species=`` when given, else ``<prefix>cycle<k>``."""
    out = []
    for rec in cycles.cycles:
        if not rec.tokens:
            continue
        orderings = [_steps(cycles, rec.tokens, rec.haplotypes)]
        orderings += [_steps(cycles, t, h if h is not None else rec.haplotypes) for t, h in rec.alternatives]
        name = prefix + (rec.fields.get("Species") or "cycle{k}".format(k=rec.cycle_id))
        out.append(Cycle(name=name, weight=rec.copy_count, orderings=orderings, circular=not rec.linear))
    return out


def raw_walks(cycles: AACycles) -> List[Tuple[List[dict], List[Tuple[int, str]]]]:
    """[(segment rows, walk)] in the layout of :func:`rck.ecdna.evaluate.read_aa_cycles` (which
    reads the simulator's truth this way): rows ``{#CHR, START = start - 1, END = end}`` in order
    of first use, walk as (1-based row index, strand)."""
    result = []
    for rec in cycles.cycles:
        if rec.linear:
            raise ValueError("{p}: cycle {c} is a path (segment 0); truth species must be cycles".format(p=cycles.path, c=rec.cycle_id))
        rows, walk, index = [], [], {}
        for sid, strand in rec.tokens:
            if sid not in index:
                chrom, start, end = cycles.segments[sid]
                rows.append({"#CHR": chrom, "START": start - 1, "END": end})
                index[sid] = len(rows)
            walk.append((index[sid], strand))
        if walk:
            result.append((rows, walk))
    return result

# Graph route (Task 2): plan and state

Goal: every decomposer gets the same breakpoint graph, so differences come from the decomposition
alone. Two graphs per case: the **truth mixture graph** (oracle) and our **inferred union
subgraph** (realistic). Methods: CoRAL 3.1.1 `cycle`, CycleExtractor (CE), original RCK,
RCK-ecDNA. Cases: demo5–demo8, then the 32 bench cases (demo3: no truth/aa, skipped).

## State (2026-11-03)

- Part A (Task 1 follow-ups) committed: RCK 41686bc, amplicon-allele-analysis 965d9eb (not pushed).
- Step 1 done: environments installed in home (section "Environments" below; lock files in
  `envs/`). Home quota is now tight (288 MB free): see there. Nothing else of Part B is written yet.
- Next session: steps 2 and 3, as planned below. Steps 4–6 are sketched at the end.

## Environments (step 1, installed 2026-11-03)

| env | path | contents | on disk |
|---|---|---|---|
| `coral-3.1.1` | `~/.conda/envs/coral-3.1.1` | conda-forge python 3.11; CoRAL 3.1.1 = commit `7f3bf317` (`pip install --no-deps git+…@7f3bf317`) and the runtime packages CoRAL imports, at its `poetry.lock` versions (gurobipy 12.0.3, pyomo 6.10.0, highspy 1.13.1, pandas 3.0.1, numpy 2.4.3, …). Left out: torch, cnvkit, pomegranate, numba, mypy, ruff, pytest (only `seed` / development use them; CoRAL never imports them) | 1.8 GB |
| `cycle-extractor` | `~/.conda/envs/cycle-extractor` | conda-forge python 3.11, numpy 2.4.3, networkx 3.6.1, gurobipy 12.0.3, highspy 1.13.1 | 0.64 GB |
| CycleExtractor code | `~/tools/CycleExtractor` | git checkout of commit `1b29b33` (`CE.py`, version 1.0.0); run as `$ENV/bin/python ~/tools/CycleExtractor/CE.py` | 2 MB |

Lock files: `envs/<env>.conda-explicit.txt` (`conda list --explicit`) and
`envs/<env>.pip-freeze.txt`. Rebuild: `conda create -n <env> --file envs/<env>.conda-explicit.txt`
then `pip install -r envs/<env>.pip-freeze.txt` (CoRAL's line pins the commit; install it with `--no-deps`).

Checks done: `coral --help` lists `cycle` and `cycle_all`; `coral cycle --help` has `--bp-graph`,
`--output-prefix`, `--solver` (gurobi_direct | scip | bonmin | couenne), `--solver-time-limit`;
`CE.py --version` = 1.0.0; CE with `--solver highs` on its sample graph reproduces its README
cycle (2+,4+,6+, CN 85.19). `coral cycle` itself has not been run (Gurobi: Slurm only).

**Home quota.** On this NFS share small files cost far more than their bytes (CE env: 0.26 GB of
packages, 0.64 GB on disk). After the install home has **288 MB free** of 5 GB; `conda clean -a`
found nothing (every cached package in `~/.conda/pkgs`, 575 MB, is linked by `gurobi-env`).
Pending user decision: drop the `cycle-extractor` env and run CE from `coral-3.1.1` (same pins of
numpy, networkx, gurobipy, highspy; frees ~0.64 GB), or move envs to project space.

**CoRAL's own writer** (`coral_rewrite.py`, run with the coral env's python): parses a graph with
`parse_breakpoint_graph`, infers intervals, writes it back with `output_breakpoint_graph_lr`.
`fixtures/coral311_amplicon1_graph.txt` (CoRAL's test graph) and
`fixtures/coral311_GBM39_amplicon1_graph.txt` (CE's sample, with path constraints) are the 3.1.1
reference layout: `#CoRAL 3.1.1` header, 4-field edge lines, `PathConstraint:` block, no
`AmpliconIntervals` block, tokens like `e2+:1,c2-:1,e3+:1`. The demo7 truth mixture graph parses
as written (2 intervals, 9 concordant, 10 discordant, no constraints).
Quirk: CoRAL's writer puts `lr_nc / size` in AverageCoverage, while its parser reads that column
as `lr_nc`, so coverage shrinks on every CoRAL rewrite. `cycle` does not use it, but the reverse
translator must not fit `depth_per_copy` from a graph CoRAL has rewritten (use reads, or the CN only).

## What the formats require (read from the pinned sources)

**CoRAL 3.1.1** (`coral/breakpoint/parse_graph.py`, `breakpoint_utilities.output_breakpoint_graph_lr`, `cli.py`):

- File name must be `<anything>_amplicon<N>_graph.txt`: `coral cycle` parses N from
  `name.split("_")[-2]`.
- Lines are split on whitespace and dispatched on the first token: `sequence`, `concordant`,
  `discordant`, `path_constraint`; everything else (headers, `source`, `#` comments,
  `interval`) is ignored. **A blank line crashes the parser** (`s[0]` on an empty list).
- `sequence  chr:start-  chr:end+  CN  AverageCoverage  Size  NumberOfLongReads` (7 fields; field 5
  is read as `lr_nc`, field 7 as an int).
- `concordant|discordant  chr:p±->chr:p±  CN  reads [...]`: field 4 must parse as a number (int(float())).
- `path_constraint  <tokens>  <support>` with tokens `e<i>±:<count>`, `c<j>±:<count>`,
  `d<k>±:<count>`: 1-based indices into the file's sequence / concordant / discordant lines, in
  file order (CoRAL's test `test_path_constraint_discordant_index.py`). The path alternates
  sequence and breakpoint edges, starting and ending on a sequence edge.
- **Amplicon intervals are not read from the file in 3.1.1** (the `AmpliconIntervals` block is
  commented out in the writer; old test graphs still carry `interval` lines, ignored).
  `infer_amplicon_intervals` builds them from runs of contiguous sequence edges
  (`next.start == prev.end + 1`), and adds an s/t end node at each run end. So the exported
  sequence edges must be exactly contiguous inside a candidate region and broken between regions.
  Writing an `AmpliconIntervals:` / `interval` block anyway is harmless and documents them.
- `cycle` uses the CN columns; read counts and `lr_nc` only feed `reconstruct`'s CN estimation
  (verify by grepping `cycle_decomposition.py` / `models/` again after the install).
- Solvers for `cycle`: `gurobi_direct` (default), scip, bonmin, couenne; no HiGHS. Defaults:
  `--solver-time-limit 7200`, `--global-time-limit 21600`, `--threads -1`.

**CycleExtractor 1b29b33** (`CE.py`, version 1.0.0):

- Sections by header lines `SequenceEdge:`, `BreakpointEdge:`, `PathConstraint:`; lines are
  `strip()`ped and blank lines skipped.
- `sequence` lines are unpacked into exactly 7 whitespace fields.
- Edge lines are split on **tabs**; with fewer than 4 fields the CN is silently 0.
- `path_constraint` lines must have exactly 3 tab fields; tokens as CoRAL's.
- **Drops discordant self-loops** (`if start != end`): a fold-back whose two breakends are the
  same node (ctrl_dup) is invisible to CE. Report it; do not work around it.
- `--solver highs` needs no licence (used for local smoke tests); default gurobi.
- No time-limit flag.

**Subgraph folder** (amplicon-allele-analysis `subgraph_construct.py`):

- `edges.tsv` types: `reference` (segments; `WEIGHT` = mean depth, `CN_TOTAL`, `CN_A`, `CN_B`),
  `extremity` (joins at a breakpoint between neighbouring segments; `WEIGHT` = spanning reads,
  `CN`), `variant` (junctions; `WEIGHT` = reads = `NOV_A + NOV_B + NOV_U`, `CN`, `junction_hap`),
  `imputed_variant` (`AID` `imputed_adj_{inter,anchor,intra}_*`; no reads).
- `read_subgraph` merges `inter` and `anchor` imputed edges into reference adjacencies
  (`_merge_parallel`); what is left as kind `IMPUTED` is a candidate closing junction.
- Node IDs: EID `chr:pos±` is a breakpoint between `pos` and `pos + 1`; `<EID>_+` is at `pos`
  (right end of the left segment), `<EID>_-` at `pos + 1`. Imputed grid nodes are
  `imputed_ext_<chr>_<k>_±`.
- `paths.tsv`: `PID PATH NSTEP #NOV`, PATH `EID@ORIENT` steps; steps (2i, 2i+1) traverse a run of
  segments, steps (2i+1, 2i+2) are one adjacency (`rck.ecdna.walk.read_paths`).
- Candidate regions: `summary.txt` "Candidate regions" / "Extended along amplified runs", and
  `ECGraph.regions()` (segments merged across ≤ 1 Mb gaps).
- RCK-ecDNA in subgraph mode reads `CN_TOTAL`/`CN_A`/`CN_B`/`CN` directly (they win over
  `scale.json`); `scale.json` gives the scales, `cell_ploidy` (baseline CN) and the sample;
  `cn_profile.tsv` the flank background.

**Coordinate frame.** `rck.ecdna.aa_export` writes every position as subgraph coordinate +
`SHIFT` (1), so that `aa_eval` reads back exactly the steps `rck-ecdna-eval` reads. The graph
exporter uses the same shift, and the reverse translator subtracts it, so a round trip is exact
and aa_eval's graph metrics stay in the frame Task 1 validated (edges within 1 bp of the truth).

## Step 2: AA-format export readable by CoRAL 3.1.1 and CE

### 2a. A shared AA graph writer

Add `write_graph(path, sequences, edges, path_constraints=(), intervals=(), header=None,
comments=())` to `rck/ecdna/aa_format.py`, next to `read_graph`, and a `path_constraint` /
`interval` reader to `read_graph` (new `AAGraph.path_constraints`, `AAGraph.intervals`, and
`AAGraph.comments` for `#rck` lines). Rules, all from the section above:

- tab-separated, no blank line, a `#` header line (`#rck-ecdna aa graph <version> <source>`);
- `SequenceEdge:` header, sequence lines sorted by (chrom, start), 7 fields;
- `BreakpointEdge:` header, then concordant lines sorted, then discordant lines sorted with the
  lower-coordinate end first (CoRAL's convention), each `CN reads None None` (6 fields: AA's
  layout, which CoRAL and CE both accept);
- `PathConstraint: Path, Support` then `path_constraint` lines;
- `AmpliconIntervals: chr, start, end` then `interval` lines (ignored by both tools);
- full-precision floats (`repr`), so the round trip is exact.

`output.write_aa_graph` (RCK-ecDNA's `aa/` answer graph) is refactored to call it, with its
bytes unchanged except the new header line; check with a before/after diff of one demo run's
`aa/*_graph.txt` written from its saved solution, or keep the old header-less layout behind a flag
if bytes must be identical (decide by `test_decompose` / regression expectations).

### 2b. Subgraph → AA (`rck.ecdna.aa_export`, new subcommand)

`python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name NAME] [--amplicon 1]
[--imputed comment|drop] [--min-path-junctions 1]` writes
`OUTDIR/<name>_amplicon<N>_graph.txt` (+ `<name>_amplicon<N>_segment_cn.tsv`,
`<name>_amplicon<N>_export.json`). Existing `aa_export RUN OUT` behaviour stays the default
(subcommand dispatch must keep `python -m rck.ecdna.aa_export RUN OUT` working).

1. `graph = subgraph_io.read_subgraph(dir)`: the normalised graph RCK-ecDNA solves on (phantoms
   contracted, parallel imputed edges merged), so every method sees the same segments.
2. Observations with `scale.build_observations` in subgraph mode (as `ecdna_run` does,
   `--inputs subgraph`), to take the exact CN RCK-ecDNA uses: segment total CN, junction CN,
   reference-adjacency CN. Check the attribute names in `scale.py` (`SegmentObservation`) first;
   fall back to the `edges.tsv` columns when an observation is missing.
3. Sequence edges: one per `graph.segments`, `start + SHIFT`, `end + SHIFT`, CN = total CN,
   AverageCoverage = the reference edge's `WEIGHT` (depth), Size = length,
   NumberReadsMapped = 0 (unused by `cycle`).
4. Concordant edges: every REFERENCE adjacency; CN = its observed CN (extremity `CN`; for
   joins at imputed grid points, where nothing is observed, the min of the two segments' CN,
   recorded in `export.json`); reads = `WEIGHT` or 0.
5. Discordant edges: every NOVEL adjacency; CN = `CN`, reads = `WEIGHT`. Self-loops (fold-backs)
   are written as they are and counted in `export.json` (CE drops them).
6. IMPUTED adjacencies are **not** CoRAL/CE edges (no read evidence; RCK-ecDNA's own
   candidates). `--imputed comment` (default) writes each as
   `#rck_imputed<TAB>chr:p±->chr:p±<TAB>support`, which both parsers skip and the reverse
   translator can restore; `drop` omits them.
7. Path constraints from `paths.tsv`: map each read path's node steps to edges with
   `graph.resolve`; a segment-run traversal (steps 2i → 2i+1) becomes `e<i>±` tokens joined by
   `c<j>±` tokens; each adjacency (2i+1 → 2i+2) a `c` or `d` token. A path that leaves the graph,
   uses an imputed edge, or has fewer than `--min-path-junctions` novel junctions is skipped
   (counted in `export.json`). Support = `#NOV`. Edge counts = occurrences in the path. CoRAL
   keeps only the longest constraints itself, so all are written.
   **Token strand and order**: do not guess. Build the expected tokens with CoRAL's own code
   (2c) and make the unit test match them.
8. Intervals: `graph.regions()` (+ SHIFT) as `interval` lines; assert they equal the runs of
   contiguous sequence edges (what CoRAL infers), else fail loudly.
9. `segment_cn.tsv` (allele CN, `chrom start end cn cn_a cn_b`) for `aa_eval --segment-cn`.

### 2c. Validate against CoRAL 3.1.1 itself

`results/graph_route/check_coral_format.py`, run with the `coral-3.1.1` env's python, no
Gurobi (the parser and writer only):

1. Done 2026-11-03 (`coral_rewrite.py`, `fixtures/`). Write a graph with CoRAL's own writer, as asked: parse CoRAL's `tests/data/amplicon1_graph.txt`
   with `parse_breakpoint_graph`, call `infer_amplicon_intervals`, and write it back with
   `output_breakpoint_graph_lr(g, out, OutputPCOptions.ALL)`. Keep that file as the reference
   layout (`results/graph_route/fixtures/coral311_amplicon1_graph.txt`), and diff its line layout
   with ours.
2. For each exported graph (demo5–demo8 subgraphs, the truth mixtures, then the 32 bench
   subgraphs): `parse_breakpoint_graph` succeeds; counts of sequence / concordant / discordant
   edges and path constraints equal `export.json`'s; `infer_amplicon_intervals` gives our
   `interval` lines; every parsed path constraint, walked through CoRAL's `path` (node list), visits
   the coordinates of the original read path. Then write it back with
   `output_breakpoint_graph_lr` and re-read it with our `aa_format.read_graph`: same edges,
   same CN, same constraints.
3. Golden tokens: for the `amplicon_0` fixture (`tests/ecdna/data`), build a CoRAL
   `BreakpointGraph` from our export, add one read path as a `PathConstraint`, and write it with
   CoRAL's `path_to_str`. Save the resulting `path_constraint` lines as a fixture; the RCK unit
   test (no CoRAL import) compares our writer's tokens with them.
4. CE smoke test, locally with HiGHS (small graphs, no licence):
   `nice python CE.py --graph demo5_amplicon1_graph.txt --output ... --solver highs`
   for demo5 and demo7: it parses (log shows the edge counts and path constraints), and writes
   cycles that `aa_format.read_cycles` reads. Note its `--s-t-strategy` default `all_nodes`.

### 2d. Truth mixture → AA for CoRAL/CE (amplicon-allele-analysis `truth_mixture.py`)

The mixture graph already parses as written (comment header, 6-field edges, filename
`mixture_amplicon<N>_graph.txt`); check it in 2c. Opt-in additions, defaults byte-identical
(test: rewrite demo6/demo7 without the new flags, compare with the committed files):

- `--path-constraints`: the simulator's long-read paths (`obs/path.tsv`) mapped onto the mixture
  graph's edges (breakends within 1 bp; a path is kept only if all its junctions are mixture
  discordant edges). This gives the oracle graph the same read information the inferred one has.
- `--reads-per-copy R` (default off): NumberOfReadPairs = round(CN × R) on discordant and
  concordant edges and AverageCoverage = CN × depth-per-copy, so the reverse translator can fit a
  scale (with 0 reads it falls back, see 3). Use the bench's 6.5 and the run's depth per copy.
- `--outdir` per method folder as today; the graph route writes into
  `/n/fs/scratch/$USER/graph_route/<case>/truth/graph/`, not into `truth/aa/`.

Tests: `simulations/tests/test_truth_mixture.py` gains a path-constraint case (demo7) and the
byte-identical default check.

## Step 3: reverse translator, AA graph → subgraph (`rck/ecdna/aa_import.py`)

`python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed none|comment|regenerate]
[--background-cn auto|X] [--cn-profile FILE] [--sample NAME]`. It writes a subgraph folder that
`rck-ecdna --inputs subgraph` reads: `nodes.tsv`, `edges.tsv`, `paths.tsv`, `scale.json`,
`cn_profile.tsv`, `summary.txt` (candidate regions), and `import.json` (what was mapped,
dropped, fitted).

1. Positions: AA − SHIFT. Breakpoints are the boundaries between contiguous sequence edges plus
   the run ends. EID `chr:pos±` with `pos` = left segment end; the sign is the strand of a
   discordant edge at that boundary when there is one, else `+` (EIDs only need to be unique and
   consistent with `paths.tsv`). Nodes `<EID>_+` at `pos`, `<EID>_-` at `pos + 1`, SOURCE
   `extremity`. Run ends get nodes like the subgraph's terminals (check what `ECGraph.terminals`
   keys on: SOURCE or degree).
2. `reference` edges from sequence edges: `WEIGHT` = AverageCoverage, `CN_TOTAL` = CN, `CN_A` /
   `CN_B` empty, `n_snps` 0, `seg_hap` `unknown`.
3. `extremity` edges from concordant edges: `WEIGHT` = reads, `CN` = CN.
4. `variant` edges from discordant edges: `WEIGHT` = reads, `NOV_U` = reads, `NOV_A` = `NOV_B` =
   0, `junction_hap` empty (every phasing allowed), `CN` = CN, `AID` `aa_d<k>`.
5. `source` edges are dropped and listed in `import.json` (AA's unexplained-flow edges; CoRAL
   adds its own s/t ends). Concordant edges between non-contiguous segments are an error.
6. Imputed candidates: `none` (no candidate edges: exactly what CoRAL/CE see), `comment`
   (restore `#rck_imputed` lines: exact inverse of 2b), `regenerate` (re-impute closing
   candidates at region ends as `subgraph_construct.py` does; needed for Task 3's CoRAL graphs,
   which carry no `#rck` lines). Find the imputation code in amplicon-allele-analysis
   (`subgraph_construct.py`, "imputed extremities and adjacencies") and reuse its rule rather
   than re-implementing it.
7. `paths.tsv` from `path_constraint` lines: tokens → node steps (`EID@ORIENT`), `#NOV` = support,
   `NSTEP`. Inverse of 2b.7.
8. `scale.json` fitted from the graph itself: `depth_per_copy` = length-weighted median of
   AverageCoverage / CN over sequence edges with CN > 0; `junction_reads_per_copy` = median of
   reads / CN over discordant edges with reads > 0; with no reads (truth without
   `--reads-per-copy`) fall back to 1.0 and mark `"junction_reads_method": "fallback: no reads"`.
   `cell_ploidy` / `background_cn` = `--background-cn` (default `auto`: the length-weighted median
   CN of the sequence edges below the amplified level, else 2). `format`
   `rck-ecdna-subgraph-scale/1`, `created_by` aa_import. Since the CN columns win over the scale,
   the fit only matters where a CN is missing; record it anyway.
9. `cn_profile.tsv`: `--cn-profile` copies a given profile (round trip: the original subgraph's);
   otherwise rows for the graph's own sequence edges, and the flanks at the background CN.

### Round-trip test (the acceptance check)

subgraph → AA (2b) → subgraph (3) → RCK-ecDNA total-only must give the same answer as
RCK-ecDNA total-only on the original subgraph.

- **Total-only on both sides.** Original: `--min-snps-allele 100000000` (as `demo7_totalonly`)
  **and** no junction haplotype constraints. Check whether `--soft-junction-hap` with a zero
  penalty is equivalent to ignoring the labels; if not, add an opt-in flag
  `--junction-hap ignore` to `rck-ecdna` (no default change) and use it on the original side
  only. `--block-flip` off on both sides (no phase blocks after import).
- **Offline part (no Gurobi), unit tests** `tests/ecdna/test_aa_import.py`: on the `amplicon_0`
  fixture and a hand-built graph, with `--imputed comment` and `--cn-profile` = the original's:
  the imported `read_subgraph` graph equals the original normalised graph (segments, adjacency
  kinds and endpoints in coordinates, CN observations, read paths as junction lists, regions),
  with node IDs compared through coordinates. Also: `source` edges dropped; non-contiguous
  concordant edge rejected; scale fit and its no-read fallback; `--imputed none` has no IMPUTED
  adjacency; paths survive both directions.
- **Solver part (Gurobi, Slurm, user submits)**: demo5–demo8 and the 32 bench subgraphs, three
  RCK-ecDNA total-only runs each: original; round trip `--imputed comment` (+ original
  cn_profile); round trip `--imputed none` (what CoRAL/CE see). Compare rank-1 decompositions
  (species walks in coordinates, weights) and aa_eval scores. Expected: comment = identical;
  none differs only where an imputed edge was used or a linear end replaces it (known users:
  ctrl_dup_s710, nest_samehap_6to1_s804) and where the flank background changes. Explain every
  other difference before step 5.

## Decisions for the user (before step 5)

1. RCK-ecDNA on the round-tripped inferred graph: `--imputed none` (strictly what CoRAL/CE see;
   recommended as the headline) or `regenerate` (RCK-ecDNA's own candidates, as it would run on
   any AA graph, e.g. Task 3)? Both can be run; the licence time doubles for that arm.
2. Truth graphs: write read counts (`--reads-per-copy 6.5`) and long-read path constraints, so
   the oracle carries the same information as the inferred graph? Recommended yes.

## Steps 4–6 (sketch, for after 2–3)

- **4. Original RCK**: inputs `--scnt` (segments × clone `c1`, allele CN; total-only graphs need
  a haploid split, e.g. `--run-haploid` or A = CN, B = 0; read docs/Segments.md and
  docs/Adjacencies.md) and `--adjacencies` (novel = discordant edges); the simulator already
  writes `truth/rck/rck.acnt.tsv` as a format reference. Score its `rck.acnt.tsv` / `rck.scnt.tsv`
  at graph level (novel adjacency FNR/FDR, allele CN), "no ecDNA species".
- **5. Runs**: one Slurm array `graph_route.sbatch` over cases × {truth, inferred} × methods,
  `%2`, skip-if-done, a `run.json` per run (method, version/commit, solver, time limit, input
  graph, command). Outputs `/n/fs/scratch/$USER/graph_route/<case>/<graph>/<method>/`. STOP
  with the sbatch commands and a licence-time estimate.
- **6. Score/collect**: `aa_eval` strict and best-of, composite-H0 for every comparison with
  CoRAL/CE; graph metrics once per graph; head-to-head table and per-scenario summary in
  `results/graph_route/README.md`.

## Constraints to keep

- Run only the test modules changed or added (name them); no `unittest discover` (Gurobi).
- No default changes in RCK-ecDNA or the simulator; new behaviour opt-in.
- Large outputs on `/n/fs/scratch/$USER/graph_route/` (purged after two weeks).
- Update docs/ecDNA.md, both CHANGELOGs and this file. No commit without showing a diff summary.
- Use `$HOME/.conda/envs/<env>/bin/python` directly; sourcing `~/.bashrc` takes minutes.

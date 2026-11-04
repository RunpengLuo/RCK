# Graph route (Task 2): plan and state

Goal: every decomposer gets the same breakpoint graph, so differences come from the decomposition
alone. Two graphs per case: the **truth mixture graph** (oracle) and our **inferred union
subgraph** (realistic). Methods: CoRAL 3.1.1 `cycle`, CycleExtractor (CE), original RCK,
RCK-ecDNA. Cases: demo5–demo8, then the 32 bench cases (demo3: no truth/aa, skipped).
How to run things and what was validated: [README.md](README.md).

## State (2026-11-03, evening)

- Part A (Task 1 follow-ups) committed: RCK 41686bc, amplicon-allele-analysis 965d9eb (not pushed);
  step 1 and this plan committed as RCK c2ecf50 (not pushed).
- Step 1 finished: the `cycle-extractor` env is removed; CE runs with the coral-3.1.1 python (same
  numpy / networkx / gurobipy / highspy) and reproduces its README cycle with HiGHS. Home: 642 MB free.
- **Steps 2 and 3 built** (sections below): `aa_format.write_graph`, `aa_export subgraph`,
  `aa_import`, truth_mixture.py `--path-constraints` / `--reads-per-copy` / `--depth-per-copy`,
  `check_coral_format.py`, `roundtrip.py`. Offline validation all passes (README "Validation").
- **Round-trip Slurm array submitted 2026-11-03: job 32096382** (`--array=0-143%2`, logs in
  `/n/fs/scratch/$USER/graph_route/roundtrip/logs/`). Next session: check it, run `compare`, explain
  differences (section "Round-trip test"), then step 4.
- **Decisions taken 2026-11-03** (section "Decisions"): head-to-head on `--imputed none`; oracle
  headline on `graph_reads_paths`, `graph_plain` secondary.
- Storage: home and `/n/fs/ugrad/ug28/ym3292` are the same 5 GB share (~600 MB free); large outputs
  are in project space, `/n/fs/yubi-iw/runs/graph_route/` (moved from scratch 2026-11-04; Slurm logs stay on scratch).

## Environments (step 1)

| env | path | contents |
|---|---|---|
| `coral-3.1.1` | `/n/fs/yubi-iw/conda/envs/coral-3.1.1` | conda-forge python 3.11; CoRAL 3.1.1 = commit `7f3bf317` and the runtime packages it imports at its `poetry.lock` versions (gurobipy 12.0.3, pyomo 6.10.0, highspy 1.13.1, pandas 3.0.1, numpy 2.4.3, networkx 3.6.1, …); 1.8 GB |
| CycleExtractor code | `/n/fs/yubi-iw/tools/CycleExtractor` | git checkout of commit `1b29b33` (`CE.py`, version 1.0.0), run as `/n/fs/yubi-iw/conda/envs/coral-3.1.1/bin/python -I /n/fs/yubi-iw/tools/CycleExtractor/CE.py` |

Lock files and rebuild: `envs/README.md`. `coral cycle` itself has not been run (Gurobi: Slurm only).
The host's conda refuses the default Anaconda channels until their terms of service are accepted;
`--override-channels -c conda-forge` works (used to remove the CE env).

**CoRAL's own writer** (`coral_rewrite.py`): `fixtures/coral311_amplicon1_graph.txt` and
`fixtures/coral311_GBM39_amplicon1_graph.txt` are the 3.1.1 reference layout (`#CoRAL 3.1.1`
header, 4-field edge lines, `PathConstraint:` block, no `AmpliconIntervals` block). Its writer puts
`lr_nc / size` in AverageCoverage (the parser reads that column as `lr_nc`), so coverage shrinks on
every rewrite, and it truncates read counts to integers.

## What the formats require (read from the pinned sources)

**CoRAL 3.1.1** (`breakpoint/parse_graph.py`, `breakpoint_utilities.output_breakpoint_graph_lr`, `cli.py`):

- File name `<anything>_amplicon<N>_graph.txt` (`coral cycle` takes N from `name.split("_")[-2]`).
- Lines split on whitespace, dispatched on the first token (`sequence`, `concordant`, `discordant`,
  `path_constraint`); everything else is ignored. **A blank line crashes the parser.**
- `sequence chr:start- chr:end+ CN AverageCoverage Size NumberOfLongReads` (field 5 read as `lr_nc`);
  breakpoint lines `kind chr:p±->chr:p± CN reads [...]` (reads `int(float())`).
- `path_constraint <tokens> <support>`, tokens `e<i>±:<n>` / `c<j>±:<n>` / `d<k>±:<n>`, 1-based
  indices in file order; written by `core_utils.path_to_str` (reverses a path whose first edge index
  exceeds its last; the strand of a token is that of the node after it; counts = occurrences).
- **Not in the original plan, found while building:** `coral cycle` calls `infer_amplicon_intervals`,
  which first calls `sort_edges`: sequence and concordant edges are re-sorted (natural chromosome
  order, `chr_sort_key`) **without renumbering the parsed path constraints**. A graph with
  constraints must already be in that order (the exporter and truth_mixture.py write it; plain string
  order would have corrupted ctrl_separate_s603, chr15 + chr7).
- Amplicon intervals are inferred from runs of contiguous sequence edges (`interval` lines ignored),
  with s/t end nodes at each run end. A self-loop discordant edge is listed twice at its node, so the
  cycle model counts a fold-back traversal twice there.
- Solvers for `cycle`: `gurobi_direct` (default), scip, bonmin, couenne; no HiGHS.

**CycleExtractor 1b29b33**: sections by `SequenceEdge:` / `BreakpointEdge:` / `PathConstraint:`;
`sequence` lines unpacked into exactly 7 fields (NumberReadsMapped must be an integer); edge lines
split on tabs (CN 0 with fewer than 4 fields); `path_constraint` lines need exactly 3 tab fields;
discordant read counts feed its multiplicity bound. **Correction of the earlier note:** CE does not
drop fold-backs; `Create_graph` turns each self-loop into a gadget (`node -> i`, 1 bp sequence edge
`i -> i+1`, `i+1 -> node`), and on ctrl_dup it uses them (README "Fold-back junctions").
`--solver highs` needs no licence; no time-limit flag; `--s-t-strategy` default `all_nodes`.

**RCK-ecDNA** (subgraph mode), what the round trip must preserve: segment CN from `CN_TOTAL`;
**junction and reference-join CN = reads / `junction_reads_per_copy`** (not the `CN` column);
flank background from `cn_profile.tsv` (else `background_cn`); and, since the flow MILP builds its
variables over `sorted(segments)`, `sorted(adjacencies)` and `sorted(nodes)`, the ids and their
order (Gurobi's path among tied optima depends on them). `--soft-junction-hap --soft-junction-penalty 0`
is exactly "junction labels ignored" (`allowed` = every labelling, no cost), so no new rck-ecdna
flag was needed for total-only runs.

## Step 2 as built: AA graphs readable by CoRAL 3.1.1 and CE

- `rck/ecdna/aa_format.py`: `write_graph` (tab-separated, no blank line, 7-field sequence and
  6-field breakpoint lines, `PathConstraint:` and `AmpliconIntervals:` blocks, `#` comments, `repr`
  floats); `read_graph` now also returns path constraints, intervals, comments (`rck_lines(tag)`)
  and sequence coverage / size / reads; `coral_path_string` (CoRAL's `path_to_str`), `chrom_key`
  (CoRAL's `chr_sort_key`), `sequence_key`, `edge_key`, `lower_first`, `parse_path_tokens`.
- `python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name] [--amplicon 1]
  [--imputed comment|drop] [--min-path-junctions 1]` (`export_subgraph`): the normalised graph with
  RCK-ecDNA's observations; table of what goes where in README "What the AA graph holds".
  `python -m rck.ecdna.aa_export RUN OUT` is unchanged.
  Deviations from the plan: junction CN = reads / θ (the plan assumed the `CN` column); flank joins
  go to `#rck_outside` (they end on a node without a segment); ids and orders in `#rck_node` /
  `#rck_seq` / `#rck_edge` / `#rck_pc` lines; a read whose outer step lies outside the trimmed graph
  (327 outer ends over the 418 constraints; `read_paths` keeps such reads) keeps the segment it leaves / enters by;
  `#rck_pc` also records whether CoRAL's token order reversed the read.
- `output.write_aa_graph` was **not** refactored: it would add a header line to every default
  rck-ecdna run's `aa/` file (a default change).
- One behaviour-preserving change in `subgraph_io.build_graph`: an `extremity` row may carry an
  `AID` (aa_import fills it, to restore reference-join ids); subgraph_construct.py never fills it
  (checked on all 181 edges.tsv files under simulations/ and tests/), so existing inputs read as before.
- amplicon-allele-analysis `truth_mixture.py`: `--path-constraints` (obs/path.tsv reads whose
  junctions are mixture discordant edges, steps placed within 1 bp; written in CoRAL's order, plus
  `mixture_amplicon<N>_path_constraints.tsv`), `--reads-per-copy R` (reads = round(CN x R)),
  `--depth-per-copy D` (AverageCoverage = CN x D). Defaults byte-identical (tested on demo6, demo7).
  Both variants for all 36 cases are in `/n/fs/yubi-iw/runs/graph_route` (`$G`): `$G/<case>/truth/graph_plain/` and
  `$G/<case>/truth/graph_reads_paths/` (R = 6.5 and D from the case's scale.json).
- `check_coral_format.py` (coral env): step 2c's checks, plus the independent rebuild of every path
  constraint on CoRAL's adjacencies, the truth sidecar check, `--golden` and `--ce`. Results: README.

## Step 3 as built: AA graph → subgraph (`rck/ecdna/aa_import.py`)

`python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate] [--background-cn auto|X]
[--cn-profile FILE] [--sample NAME] [--depth-per-copy X] [--compare ORIGINAL]` writes `nodes.tsv`,
`edges.tsv`, `paths.tsv`, `scale.json`, `cn_profile.tsv`, `summary.txt`, `import.json`.

- With `#rck` lines: the original ids, orders and flank nodes; without (any AA / CoRAL graph): EIDs
  from coordinates as planned (`chr:pos±`, sign from a discordant edge there, nodes `_+` at pos,
  `_-` at pos + 1). `source` edges dropped and listed; a concordant edge between non-neighbouring
  segments is an error.
- `scale.json` fitted from the graph: θ = median reads / CN over discordant edges (snapped to 12
  significant digits when every ratio agrees, so 6.5 comes back as 6.5); **no reads at all → WEIGHT
  = CN and θ = 1** (the plan's "fallback 1.0" alone would have made every junction observed at CN 0);
  depth = length-weighted median coverage / CN, **not fitted from a `#CoRAL` graph**. `background_cn`
  `auto` = length-weighted median of `--cn-profile`, else 2 (the plan's "median CN below the amplified
  level" has no amplified level to refer to on an amplicon-only graph).
- `cn_profile.tsv`: `--cn-profile` copied, else the sequence edges at their CN and 1 Mb flanks at the
  background CN.
- `--imputed regenerate` (`regenerate_imputed`): subgraph_construct.py's rule on sequence edges:
  threshold `min(median + 1.5 MAD, 1.5 median)` of the profile's CN bins (else 1.5 x background), an
  amplified "depth segment" = a run of contiguous amplified sequence edges without a copy-number step
  (rck-ecdna's `--boundary-min-copies 2` / `--boundary-frac 0.25`), an `intra` edge per segment and
  an `inter` edge between consecutive ones; no `anchor` edges (every segment end is a breakpoint
  already, so an anchor would be the reference join). demo7: 10 candidates after normalisation (the
  subgraph has 19), including the subgraph's cross-gap `inter` edge.
- `compare_subgraphs(original, imported)`: the comparison the round-trip test uses.

## Round-trip test

**Offline part, done** (`python roundtrip.py prepare`, `roundtrip_prepare.tsv`; unit tests
`tests/ecdna/test_aa_import.py` on amplicon_0 + synthetic reads, a hand-built graph, demo5–demo8,
ctrl_dup_s710, nest_samehap_6to1_s804, ctrl_separate_s603, nest_cross_1to1_s505): for all 36 cases,
`rt_comment` is exact, `rt_none_profile` is exact apart from the imputed candidates, `rt_none`
differs from the original only in the flank background.

**Solver part, submitted 2026-11-03 (job 32096382).** 144 runs = 36 cases x 4 arms, one per line of
`roundtrip.list` (columns: case, arm, subgraph folder, output folder, truth dir):

| arm | input | expected vs the original |
|---|---|---|
| `original` | the original subgraph | (reference) |
| `rt_comment` | import `--imputed comment`, original cn_profile and background | identical |
| `rt_none_profile` | import `--imputed none`, original cn_profile and background | different only if the original's flow uses an imputed edge (or a tie among optima is broken otherwise) |
| `rt_none` | import `--imputed none`, graph only (flank background 2) | as above, plus where the flank background changes the answer |

Every run: `rck-ecdna --inputs subgraph --min-snps-allele 100000000 --soft-junction-hap
--soft-junction-penalty 0 --step-penalty 0 --no-plots --threads 4 --time-limit 1800 --truth-dir ...`
(total-only on both sides, labels ignored, block flips off).

```bash
cd ~/RCK/results/graph_route && L=/n/fs/scratch/$USER/graph_route/roundtrip/logs; mkdir -p $L; sbatch --output=$L/%x_%A_%a.out roundtrip.sbatch   # submitted: job 32096382
python roundtrip.py compare                                                  # afterwards, offline
```

Licence estimate: the 40 bench λ = 0 runs took 1394 s in all (max 358 s), so ~4 x 23 min ≈ 1.5 h of
licence time, ~45 min wall at `%2`. `compare` writes `roundtrip_compare.tsv`: per case and arm
`identical` / `same up to haplotype labels` / `differs`, the flow and decomposition objective deltas,
and the imputed edges each flow used. Every `differs` in `rt_comment`, and every `differs` in
`rt_none_profile` without an imputed edge in the original's flow and with a different flow
objective, must be explained before step 5. Known imputed users (allele-specific bench runs):
ctrl_dup_s710 (`imputed_adj_intra_chr2_10`), nest_samehap_6to1_s804 (`imputed_adj_intra_chr3_10`);
the total-only flows may differ.

## Decisions (taken by the user, 2026-11-03)

1. **Candidate edges.** The graph-route head-to-head uses `aa_import --imputed none`: RCK-ecDNA sees
   exactly the graph CoRAL and CycleExtractor see. RCK-ecDNA's **native** run (its own subgraph,
   allele data, its own imputed edges; the bench settings) is reported as a separate `native` arm.
   `--imputed regenerate` is kept for Task 3's realistic part (graphs AA or CoRAL build from reads),
   as its own arm, **never in the head-to-head**.
2. **Truth graphs.** The oracle headline uses `graph_reads_paths` (synthetic reads, R = 6.5, D from
   the case's scale.json, path constraints). `graph_plain` runs as a secondary arm, to measure what
   path constraints (and reads) add for each method.

## Next steps, in order

1. **Check the round-trip array** (job 32096382; commands below), then
   `cd ~/RCK/results/graph_route && python roundtrip.py compare`; read `roundtrip_compare.tsv`.
   Pass: `rt_comment` `identical` for all 36. For every `differs` in `rt_none_profile` / `rt_none`:
   imputed edge in the original flow (`imputed_used`) or flank background (`rt_none` only), or a tie
   (flow objective delta 0). Anything else is a bug to explain before step 5. Record results in
   README.md and commit `roundtrip_compare.tsv`. Failed tasks: read the log, then resubmit the same
   sbatch (finished runs are skipped).
2. **Step 4, original RCK** on the AA graphs (sketch below): an `rck_input` converter (AA graph →
   `--scnt` / `--adjacencies`), a smoke test on demo7 (Slurm, licence).
3. **Step 5, runs** (`graph_route.sbatch`, one array, `%2`, skip-if-done, `run.json` per run). Arms
   per case:

   | graph | methods |
   |---|---|
   | truth `graph_reads_paths` (headline) | CoRAL `cycle`, CE, original RCK, RCK-ecDNA (`aa_import --imputed none`) |
   | truth `graph_plain` (secondary) | the same four |
   | inferred, `aa_export subgraph` (head-to-head) | CoRAL `cycle`, CE, original RCK, RCK-ecDNA (`aa_import --imputed none` + original cn_profile) |
   | `native` | RCK-ecDNA on its own subgraph, allele data, own imputed edges (bench settings, λ = 0 and the pool) |

   RCK-ecDNA on truth graphs: import with `--imputed none`; total-only (no alleles in an AA graph).
   STOP with the sbatch command and a licence-time estimate before submitting.
4. **Step 6, score/collect**: `aa_eval` strict and best-of, composite-H0 for every comparison with
   CoRAL / CE; graph metrics once per graph; head-to-head table and per-scenario summary in README.md.
5. **Task 3 (later)**: graphs AA / CoRAL build from reads; RCK-ecDNA there as its own
   `--imputed regenerate` arm (and `none`), outside the head-to-head; never fit a depth scale from
   a CoRAL-rewritten graph.

### Checking the round-trip array

```bash
squeue -u $USER -n graph-route-rt                 # running / pending tasks (at most 2 at a time)
sacct -j 32096382 --format=JobID%20,State,Elapsed,MaxRSS,ExitCode | grep -v '\.' | sort | uniq -c -f1 | head   # states
L=/n/fs/scratch/$USER/graph_route/roundtrip/logs
tail -n 3 $L/graph-route-rt_32096382_0.out        # one task's log: "== task", gurobi version, "<case> <arm>: N s"
grep -l "attempt .* failed\|Traceback" $L/*.out  # failed tasks
ls /n/fs/yubi-iw/runs/graph_route/roundtrip/*/runs/*/summary.json | wc -l    # finished runs (of 144)
seff 32096382_0                                   # memory / time of a task
```

Poll at most once a minute (the controller rate-limits). Afterwards: `python roundtrip.py compare`
prints and writes `roundtrip_compare.tsv` (here) and `compare.tsv` (`/n/fs/yubi-iw/runs/graph_route/roundtrip/`). A single run's answer
is in `runs/<arm>/decompositions.tsv` and `summary.json` under each case in `/n/fs/yubi-iw/runs/graph_route/roundtrip/`.

## Steps 4–6 (sketch, for after the round trip)

- **4. Original RCK**: inputs `--scnt` (segments x clone `c1`, allele CN; total-only graphs need a
  haploid split, e.g. `--run-haploid` or A = CN, B = 0; read docs/Segments.md and docs/Adjacencies.md)
  and `--adjacencies` (novel = discordant edges); the simulator writes `truth/rck/rck.acnt.tsv` as a
  format reference. Score its `rck.acnt.tsv` / `rck.scnt.tsv` at graph level.
- **5. Runs**: one Slurm array `graph_route.sbatch` over cases x {truth, inferred} x methods, `%2`,
  skip-if-done, a `run.json` per run (method, version/commit, solver, time limit, input graph,
  command). Outputs `/n/fs/yubi-iw/runs/graph_route/<case>/<graph>/<method>/`. CoRAL: `coral cycle
  --bp-graph G --output-prefix ... --solver gurobi_direct`; CE: coral env python, `--solver gurobi`.
  STOP with the sbatch commands and a licence-time estimate.
- **6. Score/collect**: `aa_eval` strict and best-of, composite-H0 for every comparison with CoRAL/CE;
  graph metrics once per graph; head-to-head table and per-scenario summary in README.md.

## Constraints to keep

- Run only the test modules changed or added (name them); no `unittest discover` (Gurobi).
- No default changes in RCK-ecDNA or the simulator; new behaviour opt-in.
- Large outputs on `/n/fs/yubi-iw/runs/graph_route/` (project space); temporary files and Slurm logs on scratch.
- Update docs/ecDNA.md, both CHANGELOGs, README.md and this file. No commit without showing a diff summary.
- Use `/n/fs/yubi-iw/conda/envs/<env>/bin/python` directly; sourcing `~/.bashrc` takes minutes.

# Changelog

All notable changes to this project are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/).
This file covers the `yubi_dev` branch relative to `rck/version-fix`.

## [Unreleased]

### Added

- **`rck-ecdna`** (`python -m rck.ecdna_run`, package `rck.ecdna`): allele-specific
  reconstruction of ecDNA from one amplicon breakpoint subgraph, built on RCK's diploid
  interval adjacency graph and haplotype-aware balance constraints.
  - `subgraph_io`: reads `subgraph_construct.py` folders (`nodes.tsv`, `edges.tsv`,
    `phase_blocks.tsv`) into a normalized graph:
    - contracts phantom 1 bp segments;
    - merges parallel imputed edges;
    - marks terminals;
    - attaches allowed phasings from junction haplotype labels;
    - splits regions at gaps over 1 Mb.
  - `scale`: converts depth and junction reads to copy numbers. Depth per copy comes from an
    override, HATCHet `best.seg.ucn`, the mode of balanced background blocks, or a
    tumor/normal ratio; junction reads per copy from background extremities. It also
    gives allele-specific segment copy numbers and flank background copy numbers.
  - `ilp_circular`: a two-layer flow MILP adapted from `OptModelMultiClone`:
    - a strictly balanced ecDNA layer;
    - a background layer that ends only at terminals;
    - imputed candidate edges with a usage penalty;
    - optional infinite-sites constraint.
  - `decompose`: minimum-species decomposition of the ecDNA flow:
    - integer weights and edge multiplicities, with exact linearization;
    - lazy connectivity cuts;
    - alternatives from the solution pool or no-good cuts;
    - `low_weight` flags;
    - ranking by residual and number of significant species.
  - `walk`: alternating Hierholzer walks, enumeration of distinct orderings up to rotation
    and reversal, and ranking of orderings by long-read support.
  - `output`: `decompositions.tsv`, `flow_segments.tsv`, `flow_adjacencies.tsv`,
    AmpliconArchitect-style graph and cycles files, `observations.json` and `summary.json`.
  - `evaluate` and `rck-ecdna-eval`: scoring against simulator truth:
    - truth read from `truth/aa`, with per-sample copy numbers from `proportions.tsv`;
    - truth species scoped to the subgraph's regions;
    - haplotype Jaccard, junction recall and precision, order match and copy-number error;
    - errors are reported instead of failing the run.
  - `gurobi_env`: thread defaults from Slurm, a quiet `gurobipy` logger, and model disposal.
- **Linear structures in the ecDNA layer**:
  - penalized linear ends at copy-number boundaries and region flanks
    (`--linear-ends`, `--linear-end-penalty`);
  - linear species and a `topology` column;
  - an amplicon call (`ecDNA` / `linear` / `mixed` / `none`) with `circular_fraction`;
  - linear walks written with AA's `0+`/`0-` source;
  - linear walks are never order-matched to a circular truth.
- **Long reads in the decomposition**: multi-junction paths become read constraints, with a
  capped reward (`--path-weight`) or a required explained share (`--path-fraction`, as in
  CoRAL). After solving, each read is checked in order against the orderings and
  `reads_explained` is reported.
- **Fractional copy numbers** (`refit`): after the integer structure search, species weights
  are refitted against the observations with a linear program (scipy, HiGHS). Species below
  `--refit-min-copies` are dropped. Outputs gain `weight_int`, `data_objective` and `pruned`,
  and the evaluation gains `cn_rel_err_int`.
- **Subgraph-only input mode** (`--inputs auto|subgraph|windows`): reads `scale.json`, the
  `CN_*` edge columns, `cn_profile.tsv` and `paths.tsv` written by amplicon-allele-analysis
  `subgraph_construct.py`, so no windows, extremities or HATCHet input is needed.
  `--baseline-cn` defaults to the subgraph's cell ploidy.
- **Recovery matrix** (`recovery_matrix.tsv`, and `recovery_matrix` / `mean_composite` in
  `evaluation.json`): one row per truth species.
  - Scores against the rank-1 decomposition:
    - segment recall, precision and Jaccard;
    - haplotype Jaccard;
    - junction recall, precision and F1, unphased and phased;
    - phasing accuracy, as given and A/B swap-invariant;
    - order and topology match;
    - copy-number error;
    - an equal-weight `composite` in [0, 1].
  - `best_rank`, `best_inferred` and `best_composite`: the best alternative decomposition.
  - Out-of-scope species are listed but not scored. A missed species scores 0.
  - `evaluation.tsv` gains the new per-pair columns.
- **Figures** (`rck.ecdna.plot`, `rck-ecdna-plot`): at the end of a run, `plots/<name>_ecdna.pdf`
  and one PNG per page.
  - One page per species of the top decomposition. On a cnplot `GenomeAxis` (amplicon-allele-analysis
    `third_party/cnplot`) restricted to the species' segments, it shows:
    - the haplotype-coloured topology (CoRAL-style);
    - the matched truth species;
    - the junction arcs;
    - A/B copy number: observed points, the inferred line and this species' share.
  - Beside the tracks, a ring diagram of the walk.
  - A recovery-matrix heatmap page after an evaluation.
  - Options: `--no-plots`, `--plot-ranks`, `--plot-low-weight` and `--cnplot-path` / `$CNPLOT_PATH`.
  - matplotlib and cnplot are optional. Without them, or on a plotting error, the run logs a
    warning and still succeeds.
- `bg_balance_weight`: a small tie-break toward a heterozygous background.
- **Phase-block flips** (`--block-flip unoriented`, default `off`): allele observations per
  phase-block piece, with the A/B labelling of every block the subgraph left unoriented chosen
  by the flow MILP.
  - `subgraph_io`: reads `segment_blocks.tsv` (written by amplicon-allele-analysis
    `subgraph_construct.py`) and `phase_blocks.tsv`, and maps each piece to its segment by the
    parent edge's key, through phantom contraction. Without `segment_blocks.tsv` the run logs a
    warning and keeps the segment observations.
  - `scale`: one observation per piece. Pieces below `--min-snps-allele` SNPs are total-only.
    Pieces in blocks with orientation 1 or −1 are fixed; pieces in blocks with orientation 0 or
    missing from `phase_blocks.tsv` follow the block's flip.
  - `ilp_circular`:
    - one binary flip per unoriented block, shared by its pieces across segments;
    - piece targets `a + f (b − a)` and `b + f (a − b)`, linear in the flip;
    - the segment's data weight split over its pieces by length, so the objective scale is
      unchanged;
    - a 10⁻⁶ cost per flip, so a block without evidence stays as observed;
    - after solving, `delta_if_reversed` per block, and `unresolved` blocks (reversal within
      `--flip-tol`, default 0.01) kept as observed;
    - optional `--flip-oriented-penalty` lets oriented blocks flip at a cost.
  - The flips are fixed after the flow: the refit, `flow_segments.tsv` (segment observations
    pooled over the pieces after the flips) and the CN-A / CN-B figure points use them.
  - New `block_flips.tsv`; `observations.json` gains `pieces`, and `summary.json` gains
    `block_flip` (mode, counts per status, flipped blocks). With `--block-flip off` no output
    changes.
  - demo7, on the union subgraph rebuilt with `segment_blocks.tsv` (run on ionic; outputs in
    amplicon-allele-analysis `simulations/runs/demo7_nestedcycles_ecdna/rck_out_union_blockflip/`):
    - the flow flips chr4:185 and chr4:186 and keeps chr4:171, chr4:172 and chr4:187 as
      oriented; the 12 blocks of the unamplified 172.4–185.5 Mb gap are `unresolved`;
    - the flow objective is 9.13, with no imputed edge (39.26 with `--block-flip off`, which
      reproduces `rck_out_union/`, and 8.8 with total-only observations);
    - both nested species are recovered as separate circles in the rank-1 decomposition: 16.43
      copies against a true 16.90 and 9.08 against 9.10. Composites are 0.994 and 0.9996 (0.43
      and 0 before), with segment, haplotype, junction and phasing scores of 1.0 and matching
      order. The refit drops one 0.39-copy leftover.
  - Tests:
    - `test_block_flip.py`: inputs without Gurobi; small MILPs and runs;
    - `test_demo7.py`: needs a full licence and `RCK_ECDNA_DEMO7_SUBGRAPH`, and is skipped
      otherwise. It passes on ionic;
    - `tests/ecdna/helpers.py` gains `write_segment_blocks` and `have_full_gurobi`.
  - Documented in `docs/ecDNA.md`: the subgraph folder files, inputs, the model, outputs,
    figures and the tests table.
- Tests under `tests/ecdna/`, with fixtures from demo3 and demo4: subgraph reading, scale,
  flow, decomposition, walks, refit, subgraph inputs, evaluation, recovery matrix
  (`test_recovery.py`), figures (`test_plot.py`), Gurobi environment, and end-to-end demo3 and
  demo4.
- Documentation:
  - `docs/ecDNA.md`: install, inputs, model, outputs, figures, evaluation and the recovery
    matrix, ionic, tests and subgraph-only mode;
  - a link from `README.md`.
- Princeton ionic examples in `docs/examples/`:
  - `ecdna_ionic.sbatch` and `ecdna_demo3_tests.sbatch`;
  - `ecdna_bench.sbatch`: a benchmark array, with at most 2 concurrent tasks for the academic
    Gurobi licence;
  - `collect_bench.py`: summarizes the benchmark and tallies log errors.
- `setup.py`: the `rck.ecdna` package and the `rck-ecdna`, `rck-ecdna-eval` and
  `rck-ecdna-plot` console scripts.
- Result and log files:
  - `results/plan2/`: demo3, and demo4 normal, tumor1 and tumor2;
  - `bench_summary.txt`, `test_run.log` and `test_demo3.log`.
- Validation runs of `rck-ecdna` on new simulations. Each was run from the subgraph folder
  alone, with `--truth-dir`, on a Gurobi academic WLS license. The outputs are kept in
  amplicon-allele-analysis under `simulations/runs/<run>/`, and that repository's
  `CHANGELOG.md` gives the full scores. No RCK code changed.
  - demo5, one simple cycle (`demo5_cycle1_ecdna/rck_out/`): one circular species recovered
    exactly (16.82 copies against 16.0), composite 0.99.
  - demo6, two unnested cycles on chr17 and chr15, given as one union subgraph
    (`demo6_cycle2_ecdna/rck_out_union/`): exactly two circular species, each matching its
    truth species one to one, composites 0.984 and 0.994.
  - demo7, a 3 Mb cycle with a second cycle nested inside it on chr4
    (`demo7_nestedcycles_ecdna/rck_out_union/` and `rck_out_union_totalonly/`): not
    recovered. With default settings the composites are 0.43 and 0. With total-only
    observations (`--min-snps-allele 100000000`) they are 0.51 and 0.33. With
    `--block-flip unoriented` (above; `rck_out_union_blockflip/`) both species are recovered,
    with composites 0.994 and 0.9996.
- Known limitations found by demo7:
  - **Nested species merge** (when allele data cannot separate them). When two species share
    reference segments, the minimum-species decomposition prefers one walk that traverses
    both. On demo7 with total-only observations it returned a 7.54 Mb circle at 8 copies that
    runs through the outer cycle twice and the nested cycle once. The flow itself was close to
    the truth (16–17 and 8–10 junction copies against 16.9 and 9.1). Long-read paths only
    break the tie when some read separates the species, and demo7's 3 multi-junction reads
    each lie inside one species.
    - With `--block-flip unoriented` the demo7 species no longer merge. The two species share
      9 reference segments but only one segment haplotype, because `ecdna1` is on B and
      `ecdna2` on A where they overlap. Once the allele observations are correctly phased, a
      single walk cannot reproduce the haplotype-specific flow more cheaply than two circles.
    - The merge is still expected for nested species on the same haplotype; no simulation
      tests that case yet. Possible remedies: penalize walks that repeat a sub-cycle, use
      copy-number differences between species, or require more read support.
  - **Allele observations are trusted as phased** (by default). Segment allele copy numbers
    are taken on the true haplotype scale. When `subgraph_construct.py` cannot orient a phase
    block, here because the only pure junctions there attach to segments shared by both
    species, A and B are swapped against the junction labels. On demo7 this made the flow
    objective 39.3 instead of 8.8 and brought in imputed edges. `--block-flip unoriented`
    (above) lets the flow choose a per-block A/B flip for blocks whose `phase_blocks.tsv`
    orientation is 0 or missing. It needs a subgraph with `segment_blocks.tsv`. On demo7 it
    brings the flow objective to 9.13 with no imputed edge. Blocks the subgraph oriented are
    still trusted, unless `--flip-oriented-penalty` is set. It remains opt-in: `off` is the
    default.

### Changed

- `rck_run.py`: `MIPGapAbs` is set only from the new `--run-g-mip-gap-abs` option. It was
  previously set to the relative `--run-g-mip-gap` value, which forced near-exact
  optimality.

### Fixed

- `rck_run.py`: `--scnb` copy-number boundaries are read from the `--scnb` file. They were
  previously read from the `--scnt` path, so `--scnb` was silently ignored.
- ecDNA depth-scale mode estimate: wrong results when the smoothing kernel was longer than the
  histogram.
- A circle traversed several times is reported as one circle at the summed weight (gcd
  normalization).
- Decompositions are ranked by significant species, so one molecule plus small leftovers
  ranks above the molecule split in two.
- Gurobi models are disposed after solving, and test temp folders work on a full `/tmp` and
  on NFS.
- `evaluation.json` is strict JSON: missing scores in `top_decomposition` are `null`, where they
  were previously `NaN`.
- `--truth-dir` evaluation works for any simulated run and subgraph: several species,
  deleted segments and per-sample copy numbers.

[Unreleased]: https://github.com/RunpengLuo/RCK/compare/rck/version-fix...yubi_dev

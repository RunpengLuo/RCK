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

# bench_IW: status and next steps (branch `yubi_CNP_ionic`)

*Written 2026-11-14, at the end of the session that built bench_IW, while its methods array was
still running. Changes and their reasons are in the `[yubi_CNP_ionic]` sections of both repos'
`CHANGELOG.md` (entries dated 2026-11-14).*

## Resume here (next session)

**State.**
- **Simulation side done and committed** (amplicon-allele-analysis, `yubi_CNP_ionic`): 24 accepted
  runs, union subgraphs (bulk route) and multi-clone DIAGs. All 24 share demo7's clones
  (`configs/bench_IW/shared_clones.tsv`).
- **RCK side committed:** the metrics (`rck/ecdna/bench_metrics.py`), the driver
  `results/bench_IW/bench_iw.py` and `bench_iw.sbatch`.
- **Methods array submitted:** Slurm job **32267060** (96 tasks = 24 cases x 4 methods, `%2`
  for the Gurobi WLS licence, 4 CPUs, 4G, 3 h each).
  - When this file was written (15:32, 2026-11-14), 13 of 96 runs were done: the 4 local pilot
    runs on `simple_S1_s11001` plus array tasks up to 12.
  - Median run time so far is 3 s, the longest 365 s. Expected total: 1–3 h, mostly CoRAL and
    the bulk pool runs on the S=4 cases.
- **Not done yet:** scoring and collecting over all 96 runs, the results write-up, and
  Yubi's decisions on CoRAL / CycleExtractor (below).
- Nothing pushed. Yubi pushes.

**Rules.** Small local runs are fine. New Slurm arrays, sweeps or long limits need Yubi's
approval (the licence allows 2 Gurobi sessions, so never run local Gurobi jobs while the array
holds both). Wait at least 60 s between `squeue` / `sacct` calls. Do not touch `~/RCK` or
`~/amplicon-allele-analysis` (`yubi_dev`).

**Next steps, in order:**
1. **Check the array.**
   ```bash
   cd ~/RCK-cnp/results/bench_IW
   export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
   squeue -u $USER                                   # anything left of 32267060?
   sacct -j 32267060 --format=JobID,State,Elapsed,MaxRSS -X | sort -k2 | uniq -c -f1 | head   # states
   python bench_iw.py status                         # done / missing / failed per (case, method)
   ls /n/fs/yubi-iw/runs/bench_IW/logs/slurm/        # bench-iw_32267060_<task>.out per task
   ```
   - A run is `done` when its outputs exist; `failed` when its `run.json` says so.
   - Look at the task log, and at `<run>/run.json` (command, attempts, elapsed). CoRAL / CE also
     write `<method>.stdout` in the run folder.
   - **If some failed or never ran:** diagnose first (licence error? time limit? a bug?). Only
     then resubmit with `sbatch --array=<idx,...>%2 --output=... bench_iw.sbatch`; finished runs
     are skipped, so the full array can also be resubmitted. Ask Yubi before resubmitting more
     than a handful.
   - `sacct` MaxRSS: if every task is far below 4G, note it for later arrays (`--mem`).
2. **Score and collect (offline, no Gurobi, local).**
   ```bash
   python bench_iw.py score --jobs 8    # export + aa_eval + bench_metrics per finished run
   python bench_iw.py collect           # bench_iw_table.tsv, bench_iw_pooled.tsv, bench_iw_proportions.tsv here
   ```
   - `score` rescores every finished run and overwrites `<run>/score/`. Runs that fail to score
     leave `<run>/score/score_error.txt`, a traceback; `collect` marks them `score_status=error`.
   - The three TSVs are written into `results/bench_IW/`. Commit them with a short README of the
     results (step 3).
3. **Report to Yubi:**
   - the per-(method, class, S) pooled table;
   - failures and rejections;
   - the CoRAL / CycleExtractor behaviour below, across all 24 cases.

   Do not change the metrics or methods before Yubi decides.
4. **Decisions waiting for the full results (Yubi):**
   - CoRAL and CycleExtractor answers include linear walks (`0+ ... 0-`) through non-ecDNA
     segments. They count as ecDNA mass, which lowers copy-number precision and proportion
     scores. Should precision (or everything) drop linear walks, as `score_h2h.normalise_cycles(...,
     drop_paths=True)` can? Recall is unaffected.
   - Whether the bulk route should also be scored with its `lambda_0` answer (not run here; only
     the pool run with lambda 0.2 re-ranking is).

**Pilot (`simple_S1_s11001`, one species, 4.97 copies; run locally before submitting):**

| method | time | cn_recall | cn_precision | proportion_tv (strict / all) | missed / false |
|---|---|---|---|---|---|
| rck_bulk | 1 s | 1.00 | 0.93 | 0 / 0 | 0 / 0 |
| rck_mc | 5 s | 1.00 | 0.96 | 0 / 0 | 0 / 0 |
| coral | 58 s | 1.00 | 0.32 | 0.11 / 0.22 | 0 / 1 |
| ce | 1 s | 1.00 | 0.16 | 0.58 / 0.58 | 0 / 2 |

What the answer files show (no further tests were run):
- **CoRAL** finds the species as a closed cycle (3.96 copies). It then adds three low-copy linear
  walks (0.50, 0.35, 0.27) over segments the species does not use, including the 10 kb flanks
  `graph_route.run_coral` adds.
  - `max_weight` mode keeps decomposing the remaining copy number, which is mostly chromosomal
    background: the exported graph's copy numbers are total copy number.
  - The long-read path constraints are credited to one of the linear walks, not to the cycle. Why
    is not clear from the files.
- **CycleExtractor** returns no closed cycle. It returns three linear walks (2.73, 2.57, 1.21):
  - two carry ecdna1's segments, each repeating part of it;
  - one is a single non-ecDNA segment.
  - Their sum, 6.51, is close to ecDNA plus about 2 copies of background. So it seems to spread
    the ecDNA copies over background-carrying walks (an inference, not checked against its code).

## What bench_IW is

- **Classes.**
  - `simple_S{1,2,4}`: species share no reference base. `dup` is kept, at most 2 traversals of a
    segment.
  - `hard_S{1,2,4}`: S >= 2 nests every species inside the first one's intervals (outer
    0.6–3 Mb). S = 1 is one ctrl_dup-like species with internal repeats or a fold-back; at least
    one segment traversed twice, at most 3.
  - 4 accepted seeds per scenario: 24 runs.
- **Shared ground truth.** Every run has `clone_seed: 7` and demo7's realised clone proportions
  `[0.2, 0.49981378091514816, 0.3001862190848519]`. These reproduce demo7's grid, phasing, clone
  tree, clone CNPs and proportions exactly. Also shared:
  - `ecdna_positive_clones` null (both tumour clones), `ecdna_clone_weights` [1, 0.5];
  - long reads, library 5e6 (about 6.5 junction reads per copy);
  - `hap_labels: inherit`, `p_cross` 0.2.
- **Drawn per seed** (`make_scenarios.iw_draw`, from `default_rng(random_seed)`, written into
  each config and its `bench_iw:` record block):
  - **size per species:** lognormal, median 1 Mb, sigma 1.3, within 0.3–3 Mb; simple sets
    redrawn until they sum to <= 8 Mb;
  - **hard nesting:** `nest_size_fraction` U(0.3, 0.7) per nested species;
  - **copies per species:** lognormal, median 8, sigma 1.2, within 4–60; sets redrawn until
    max/min <= 8 (`abundance_mode: per_sample`, `ecdna_burden` = sum, `ecdna_proportions` =
    shares);
  - **`num_intervals`:** simple U{1,2}, hard U{2,3};
  - **breakpoints per species:** U{2 .. min(8 or 10, ceil(size / 150 kb))}.

  The bounds are AmpliconRepository statistics (2,963 human ecDNA features); the reasons are in
  `make_scenarios.py`.
- **Seed checks** (`bench_nested.check_seed`, intent from `make_scenarios.IW_INTENT`):
  - bulk ratio within 10%;
  - every ecDNA junction has >= 20 reads;
  - traversal cap;
  - footprint <= 8 Mb;
  - every species placed with all its intervals;
  - separate / nested (>= 80% inside) / dup.

  Accepted 24 of 30 candidates. Rejections are in `configs/bench_IW/seeds.tsv`.
- **Subgraphs.** bench_nested's flags plus `--min-component-nodes 2 --min-sv-support 20`. The
  default (>= 5 breakpoint clusters) dropped a one-interval species in `simple_S4_s13004`. Every
  union subgraph has all truth segments and junctions (`subgraphs/bench_IW/subgraph_checks.tsv`).
- **Methods** (`bench_iw.py` docstring has the exact flags):
  - `rck_bulk`: pool run, `--block-flip unoriented`, answer re-ranked at lambda 0.2;
  - `rck_mc`: `--model multiclone` on the DIAG built with `--clone-profiles truth
    --diag-labeling-groups block`;
  - `coral`: CoRAL 3.1.1 on the exported union subgraph, with 10 kb flanks;
  - `ce`: CycleExtractor 1.0.0 on the same graph.
- **Metrics,** per (case, method), in `<run>/score/`:
  - **`aa_eval`:** `aa_species_matrix.tsv` (J_seg, J_adj; blind rows for every method, phased
    rows meaningful for `rck_bulk` only), M1 / M2 (`aa_partition_metrics.tsv`), strict and
    best-of species scores. No truth mixture graph was built, so there are no graph-level
    metrics.
  - **`bench_metrics`** (`rck/ecdna/bench_metrics.py`, approved definitions):
    - `cn_recall` = sum_x min(c, t) / sum_x t, over every base. Ecdna copies only,
      haplotype-blind, bp-weighted, strict cycles. t = sum_s PHI_s x traversals; c = sum of
      cycle bulk weight x traversals. `rck_mc` weights are sum_k u_k psi_k, with u from the DIAG.
    - `cn_precision` = the same over sum_x c.
    - `cn_within_tol`: the share of truth bases with |c - t| <= max(1, 0.1 t).
    - `_all` variants include low_weight cycles.
    - Proportions: truth PHI_s / sum vs inferred w / sum over the strict set and over all
      cycles, greedy blind-Jaccard matching. `proportion_tv` is the total-variation distance
      (unmatched mass counted); `proportion_mae_matched`.
- **Table columns** (`bench_iw_table.tsv`, one row per case x method):
  - case fields: `name`, `scenario`, `class`, `S`, `replicate`, `seed`;
  - run fields: `method`, `run_status`, `elapsed_s`;
  - the bench_metrics fields;
  - `n_inferred_strict`, `n_missed`, `n_false`, `mean_composite_h0`;
  - M1 / M2 headline precision / recall / F1 (`species_detection_*`, `adjacency_pairs_*`);
  - `{blind,phased}_jseg_best_mean`, `_jadj_best_mean`: the mean over truth species of the
    best strict cell; `_exact_cycles`.
- **`bench_iw_pooled.tsv`:** mean and median per (method, class, S), (method, class) and method;
  the `level` column says which.

## Where things are

| what | where |
|---|---|
| configs, seeds.tsv, shared_clones.tsv | amplicon-allele-analysis-cnp `simulations/configs/bench_IW/` (tracked) |
| simulator runs | `/n/fs/yubi-iw/runs/bench_IW/runs/` (linked as `simulations/runs/bench_IW`, ignored) |
| union subgraphs / multi-clone DIAGs | `simulations/subgraphs/bench_IW/`, `subgraphs/bench_IW_cnp/` (tracked) |
| simulation and subgraph logs | `/n/fs/yubi-iw/runs/bench_IW/logs/` |
| method runs and scores | `/n/fs/yubi-iw/runs/bench_IW/methods/<case>/<method>[/amplicon1]/{run.json,score/}` |
| exported graphs (CoRAL, CE input) | `/n/fs/yubi-iw/runs/bench_IW/methods/<case>/graph/` |
| Slurm logs | `/n/fs/yubi-iw/runs/bench_IW/logs/slurm/bench-iw_32267060_<task>.out` |
| driver, sbatch, results TSVs | RCK-cnp `results/bench_IW/` |

Commands on the simulation side (amplicon-allele-analysis-cnp `simulations/`):
- `python bench_IW.py seeds | report | check-shared`;
- `python bench_IW.py [--cnp] rebuild [RUN ...]`;
- `python bench_IW.py scenario <name> --accept 4 --candidates 20` (resumable; existing runs are
  reused).

**Session commits (2026-11-14), all on `yubi_CNP_ionic`, none pushed.**
- **amplicon-allele-analysis:**
  - `fead61a` `clone_seed`, `6283c0e` its test fix;
  - `0e05619` old bench configs and subgraphs removed;
  - `0c3a727` per-species structure settings;
  - `2c8f810` bench_IW code;
  - `21cce8c` bench_IW configs and subgraphs;
  - plus the docs commit that wrote this session's CHANGELOG / README entries.
- **RCK:**
  - `1b5bfd3` bench result data removed (scripts kept);
  - `eef70e6` bench_metrics, driver, sbatch;
  - plus the docs commit that wrote this file.
- **Local only:** in amplicon-allele-analysis-cnp, the `third_party/ecSimulator` (pinned
  `07a0e3c`) and `third_party/cnplot` (pinned `1e8bfef`) submodules were initialised. Without
  ecSimulator the simulator writes no `truth/aa/`, and the seed checker's traversal count
  reads `truth/aa/*_cycles.txt`.

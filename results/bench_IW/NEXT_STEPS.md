# bench_IW: status and next steps (branch `yubi_CNP_ionic`)

*First written 2026-11-14 at the end of the session that built bench_IW; "Resume here" rewritten
the same evening, after the methods array was checked and Yubi chose whole-number psi. Changes
and their reasons are in the `[yubi_CNP_ionic]` sections of both repos' `CHANGELOG.md` (entries
dated 2026-11-14, evening).*

## Resume here (next session)

**State (2026-11-14, evening).**
- **bench_IW's methods array (32267060) is finished and checked.** It is not the benchmark's
  result: bench_IW has a mismatch between the multi-clone input and the scoring truth (below).
  - 92 of 96 tasks completed; no time-outs, no out-of-memory kills (peak 754 MB of 4 GB).
  - Failed:
    - CoRAL on simple_S1_s11004 and hard_S2_s15004 (tasks 14 and 74): the Gurobi licence was
      busy on all 5 attempts;
    - CycleExtractor on simple_S4_s13002 and hard_S2_s15005 (tasks 39 and 79): crashes inside
      CycleExtractor 1.0.0 itself (`NameError: Final_closed_walk`, `CE.py:3757`;
      `NameError: segments`, `CE.py:3954`). A resubmission fails the same way.
  - Stopped at the time limit: `rck_bulk` decomposition on hard_S4_s16003 (gap 85%); CoRAL on
    hard_S4_s16002 and simple_S4_s13002 (gap 9%).
  - `rck_mc`: every solve optimal, longest 3.95 s (hard_S4_s16003). This is the timing Yubi
    asked for before fixing the multi-clone time limit: **60 s stays**.
  - CycleExtractor starts Gurobi with 32 threads on a 4-CPU allocation.
  - Compute per method (summed run.json elapsed, licence waits included): `rck_bulk` 1.09 h,
    `coral` 1.87 h, `ce` 0.97 h, `rck_mc` a few minutes.
  - Not scored or collected.
- **Found: the multi-clone input and the scoring truth disagreed.**
  - psi per clone is real-valued, and the reads and `bench_metrics`' truth use it
    (PHI_s = sum_k u_k psi_k).
  - The clone karyotypes in `truth/rck/`, which are `--clone-profiles truth`, fold in
    round(psi).
  - So `rck_mc` could at best report sum_k u_k round(psi_k), up to 4.7% from PHI_s.
- **Decision (Yubi): whole-number psi by construction (option A).** The rejected alternatives
  and the reasons are in both CHANGELOGs.
- **Done for A (amplicon-allele-analysis, `yubi_CNP_ionic`):**
  - simulator flag `l3_mixture.integer_psi` (commit `2363469`, 7 tests);
  - set `iw_intpsi`;
  - **bench_IW_intpsi**: 24 runs re-simulated from bench_IW's own seeds, both subgraphs built
    (commits `0a164dc` data, `0dc2d69` docs).
  - All 48 subgraphs pass their checks; `check-shared`: 24 runs, 0 differ.
  - `truth/rck/` copy numbers, which are the multi-clone input, are byte-identical to bench_IW's:
    round(psi) did not change. What moved to match them is the reads and the scoring truth.
  - One 1-bp DIAG sliver differs in hard_S2_s15005 (chr5:152,315,000; CHANGELOG). Not changed.
  - The same 24 seeds were accepted. Structures, clones, `segment_max_traversals` and species
    sizes are identical to bench_IW; only psi, PHI and the reads changed.
  - psi is whole (3–85) and PHI_s = sum_k u_k psi_k in every run.
  - Details: that repo's CHANGELOG and `simulations/README.md` 5.5.
- Nothing pushed. Yubi pushes.

**Decisions that apply from here** (Yubi's; they override the earlier text below where they differ):
1. Clone-level ecDNA copies are whole numbers by construction (bench_IW_intpsi). The scoring
   truth stays PHI_s = sum_k u_k psi_k, which is now also what the multi-clone input carries.
2. CoRAL and CycleExtractor run on **RCK bulk's exported graph** (`<case>/rck_bulk/aa/*_graph.txt`),
   not on the union subgraph exported directly as array 32267060 did.
   - Label them *decomposition-only baselines*.
   - Each case's graph-route tasks depend on that case's bulk task.
3. The results table has a **`strict_rule` column** (`Status` or `cn-frac:0.1`).
   `bench_metrics.json` already records it; check that `collect` carries it.
4. **Scores:**
   - Headline: ecDNA-only, haplotype-blind, bp-weighted copy-number recall.
   - For `rck_bulk` and `rck_mc` also **total copy-number recall** (background plus ecDNA).
   - Proportions over strict cycles, with the `low_weight`-inclusive version too.
5. **Per run, record** `segment_max_traversals`, species sizes and the strict rule.
6. **Multi-clone time limit:** 60 s, from the timing above.

**Rules.**
- Small local runs are fine. New Slurm arrays need Yubi's go-ahead.
- The licence allows 2 Gurobi sessions: never run local Gurobi jobs while arrays hold both.
- Wait at least 60 s between `squeue` / `sacct` calls.
- Do not modify bench_IW's configs, seeds, runs or array outputs (kept for the record).
- Do not touch `~/RCK` or `~/amplicon-allele-analysis` (`yubi_dev`).

**Next steps, in order:**
1. **Done (2026-11-14 night, `9620a8a`): driver points at bench_IW_intpsi.** `bench_iw.py --set
   iw_intpsi` switches the case list, simulator runs, subgraphs, method outputs
   (`/n/fs/yubi-iw/runs/bench_IW_intpsi/methods`) and the results TSVs (`bench_iw_intpsi_*.tsv`).
   Default `iw` unchanged.
2. **Done (2026-11-14 night): CoRAL / CE on the bulk route's graph, as decomposition-only
   baselines.** Details and checks: CHANGELOG, "bench_IW_intpsi driver; graph route on the bulk
   route's graph".
   - Graph: `rck.ecdna.aa_export.baseline_graph`, the bulk run's `aa/*_graph.txt` (+1 shift,
     fitted copy numbers, every edge it used) plus the direct subgraph export's long-read path
     constraints, renumbered by coordinate.
   - *Why the path constraints* (Yubi chose to keep them, to stay close to the original tools):
     CoRAL uses them by default as MILP constraints, CE reads them, and the bulk route reads the
     same paths.
   - Written once per case to `<case>/graph_bulk/`. A coral / ce run exits 2 without running
     until the case's `rck_bulk` is done.
   - `run.json`: `baseline: decomposition-only`, `graph_source`, the graph's counts.
   - CE's Gurobi threads are capped by a `gurobi.env` in its run folder; **to confirm in the
     pilot.**
   - Licence errors: up to 15 retries, 2 min apart.
   - `bench_iw_method.sbatch` + `submit_intpsi.sh` (prints the four commands; `--submit` submits).
   - Offline check on bench_IW's 24 bulk runs: all 389 path constraints carried, none dropped;
     3,391 tokens point at the same coordinates.
3. **Total copy-number recall for the RCK routes, and the per-run fields** (decisions 4–5; local,
   about 1 h with tests).
   - `bench_metrics`: total recall = sum_x min(c, t) / sum_x t over the subgraph's segments,
     haplotype-blind, bp-weighted.
     - t = tumour total copy number (background plus ecDNA, `truth/segments.tsv` CN_A + CN_B).
     - `rck_bulk` c: fitted segment copy numbers (`flow_segments.tsv`).
     - `rck_mc` c: sum_k u_k x the clone's inferred copy number (its `rck.scnt.tsv`).
   - `collect` joins `segment_max_traversals` and species sizes (`truth/ecdna.tsv` SIZE) and
     keeps `strict_rule`.
   - Tests for each.
4. **Pilot locally** on simple_S1_s11001 and one hard S=4 case: all four methods plus scoring. It
   uses 1–2 Gurobi sessions for a few minutes; no array may be running.
5. **Methods arrays (needs Yubi's go-ahead):** `./submit_intpsi.sh` prints the commands,
   `./submit_intpsi.sh --submit` submits. Four 24-task arrays, one task per case, 4 CPUs and 2 G
   each (peak 754 MB on bench_IW):

   | array | method | depends on | time limit | wall time (array 32267060 timings) |
   |---|---|---|---|---|
   | A | `rck_bulk` | – | 90 min (max 1908 s, licence waits included) | about 1.1 h at %1 |
   | B | `rck_mc` | – | 15 min (solves <= 4 s) | a few minutes |
   | C | `coral` | `aftercorr:A` (case i after bulk task i succeeds) | 90 min (max 2002 s) | about 1.9 h |
   | D | `ce` | `aftercorr:A` | 60 min (max 374 s; CE's own limit is 7200 s per MILP) | about 1 h |

   - **Gurobi licence (2 sessions):** A and B hold at most one each. C and D join as bulk tasks
     finish, so up to 4 may ask at once; a task waits up to 15 x 2 min.
   - **A failed bulk task** leaves that case's C / D pending (DependencyNeverSatisfied): cancel
     them, diagnose, and resubmit (finished runs are skipped).
   - **Expected total:** about 2–2.5 h wall time, about 4 h compute.
   - CoRAL / CE timings on the bulk route's graph may differ from these.
   - CE crashed on two bench_IW graphs; report crashes as method failures, do not patch CE.
6. **Score and collect (local, no Gurobi):** `bench_iw.py --set iw_intpsi score --jobs 8`, then
   `collect`.
   - Scoring time not yet measured; expect well under 1 h.
   - Commit the TSVs with a short results README.
7. **Report to Yubi:**
   - the pooled table per (method, class, S);
   - failures;
   - CoRAL / CE behaviour.
   - Open decision then: should CoRAL / CE linear walks (`0+ ... 0-`) count in precision and
     proportions (see "Decisions waiting" below)?

## Earlier state (2026-11-14 afternoon, kept for the record)

> Superseded by "Resume here" above. Array 32267060 ran on bench_IW, whose multi-clone input and
> scoring truth disagree; its outputs are kept but are not the benchmark's results.

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
| **bench_IW_intpsi** (the benchmark from 2026-11-14 evening): configs, seeds.tsv | `simulations/configs/bench_IW_intpsi/` (tracked) |
| bench_IW_intpsi simulator runs | `/n/fs/yubi-iw/runs/bench_IW_intpsi/runs/` (linked as `simulations/runs/bench_IW_intpsi`, ignored) |
| bench_IW_intpsi subgraphs | `simulations/subgraphs/bench_IW_intpsi/`, `subgraphs/bench_IW_intpsi_cnp/` (tracked) |
| bench_IW_intpsi simulation and subgraph logs | `/n/fs/yubi-iw/runs/bench_IW_intpsi/logs/` |
| bench_IW_intpsi method runs (planned) | `/n/fs/yubi-iw/runs/bench_IW_intpsi/methods/` |

Commands on the simulation side (amplicon-allele-analysis-cnp `simulations/`; add `--set iw_intpsi`
as the first argument for bench_IW_intpsi):
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

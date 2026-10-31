# bench_nested: next steps (written 2026-10-27)

State at the end of the 2026-10-27 session, and a plan for what comes next; progress of the
2026-10-29 session is marked **[done 2026-10-29]** and in "Progress (2026-10-29)" at the end. The
evidence is in `README.md`: the "bench_nested_v2" and "Follow-up checks" sections. The full
old-vs-new comparison is in `comparison.md` and `comparison.tsv` here.

Scratch purges files older than two weeks. Only run outputs remain there (`bench_nested_v2/runs`,
`rck_runs`, `diag`, `subgraphs_paths`); configs, subgraphs and tables are in the repositories, and a
copy of the v2 configs, subgraphs and `rck_results` is in `~/backups/bench_nested_v2_2026-10-29.tar.gz`.

> **Done 2026-10-31:** array 32025976 finished (all 22 tasks COMPLETED). `diag/rck_noise/SUMMARY.md`
> is complete, except the s403 w0.25 row, which hit Gurobi's `--time-limit` in its pool search and is
> described rather than filled. The 2026-10-31 assessment of the plan is "Assessment (2026-10-31)"
> at the end. Its revised plan replaces "Next steps (revised 2026-10-29; ...)", which is kept for the
> record. Scratch still purges the full run outputs around 2026-11-12. The small outputs are
> archived in `diag/rck_noise/runs/`.

## Where things stand

- **Two benchmarks exist.**
  - The first (committed tables in this folder) used per-step haplotype labels.
  - **bench_nested_v2** (inherited labels, `p_cross` 0.5) lives only on scratch:
    `/n/fs/scratch/$USER/bench_nested_v2/{configs,runs,subgraphs,rck_runs,rck_results,diag}`.
    Its configs, `seeds.tsv`, subgraphs and tables are **not** in either repo.
  - Scratch is not backed up.
- **v2 replaced the first benchmark on 2026-10-29** (step 1, done).
- Recommended `--step-penalty`: **0.2**, not 0.1. On v2, `nest_cross_3sp_s402` switches to the right
  answer only at λ ≥ 0.11, and λ = 0.1 is an exact tie. No control or demo case loses at any λ up
  to 0.5. RCK's default stays 0.
- Remaining v2 failures at λ = 0.2, and what each one is:

| case | failure | cause (from the follow-up checks) |
|---|---|---|
| nest_cross_1to1_s502 | 3 species for 2 (merged) | junction reads per copy: 6.5 is used, the true value is ≈ 5.0. The fitted 5.11 resolves it (0.963) |
| nest_cross_1to1_s505 | 3 species for 2, 0.765 | decomposition: no good solution in the pool; θ 6.5 is right here (fitted 4.93 makes it 0.600) |
| nest_samehap_3to1_s301 | 3 species for 2, 0.693 | decomposition ranking: a 0.973 answer is in the pool but ranked below |
| ctrl_dup_s703 | 2 species for 1, 0.712 | decomposition ranking: a 0.971 answer is in the pool. A −16% depth bias on one segment (simulator `sigma` 0.1, by design) makes the exact fit split the species |
| nest_cross_3sp_s403 | merged, 0.592 | (a) subgraph: `rck_inputs` fits ploidy 6 (scores 0.114 vs 0.122 for ploidy 2), so segment CN is 2× truth; (b) with `--rck-ploidy 2` the scale is right, but RCK still gives 0.545 and the best pool solution is 0.729 |

## Root causes to fix

1. **RCK decomposition: exact cover of a noisy integer flow** (`rck/ecdna/decompose.py`).
   - The flow MILP rounds noisy observations to integers.
   - The decomposition then covers that flow with a residual weight of 1 per copy: `_omega` is 1 on
     segments ≥ `seg_length_cap` (50 kb), and `adj_residual_weight` is 1.0 on junctions. The species
     penalty is also 1.
   - So a 1–3 copy dip from depth bias (±10% per 200 kb, ±2 copies at 20) costs as much as one or
     more extra species. Exact-fit alternatives with leftover or merged species win.
   - This affects s301, s703, s505 and s403.
2. **`rck_inputs.py` purity/ploidy selection** (amplicon-allele-analysis): ploidy classes are
   chosen on near-tied scores. It was wrong for v1 s501, s502 and s406 and for v2 s403; v2 s502
   was right by 0.002.
3. **Junction reads per copy (θ)**: the bench fixes 6.5 for every case. The truth-implied values
   range from about 5.0 to 6.5, and the fitted value is unreliable: right for s502, wrong for s505
   (4.93) and s403 (2.12).
4. **The benchmark's stage rule** (`collect.py`, `pipeline.tsv`) only checks coverage and missing
   junctions. It labels scale errors "RCK". `verify_subgraph.py` reports false `seg_calls_wrong` where
   nested species on different haplotypes overlap.

## Plan, in order

### Step 1. Decide on the benchmark (you; no compute) [done 2026-10-29: v2 adopted]

Read `comparison.md`. If v2 replaces v1:

- copy `/n/fs/scratch/$USER/bench_nested_v2/configs/` into amplicon-allele-analysis
  `simulations/configs/bench_nested/`;
- copy the subgraph folders into `simulations/subgraphs/bench_nested/` (check their size first:
  `du -sh` on that one known folder; home quota is 5 GB);
- copy `rck_results/{bench_sweep,pipeline,rerank_validation}.tsv` over this folder's tables;
- point `runs` at `bench_nested_v2/rck_runs`.

Until then, back up the v2 inputs: `tar czf` the `configs` and `subgraphs` folders to a second
location. The simulations can be regenerated from the configs with `bench_nested.py`.

### Step 2. Test a noise-tolerant decomposition (RCK; small code change, ≈ 1 h of runs) [done 2026-10-29: no setting passes; see diag/rck_noise/SUMMARY.md]

1. Expose two options in `rck/ecdna_run.py` / `DecompositionParams`: `--adj-residual-weight`
   (exists as a param, not on the CLI), and a segment residual scale. The segment scale should be
   either a multiplier on `_omega` or a relative mode, where residual is divided by an expected
   noise of σ × max(flow, 1) with σ ≈ 0.1.
2. Run the four decomposition failures (s301, s703, s505, s403 with `--rck-ploidy 2`) at a few
   settings, for example weights 1, 0.5, 0.25 and relative σ = 0.1. Use `diag/rck_diag.sbatch`
   with a LIST file. Lines are `NAME SUBGRAPH TRUTH OUTDIR [RCK OPTIONS]`, and the array is run
   with `%1`, because the WLS licence allows 2 sessions.
3. Summarize with `diag/summarize.sh <outdirs>` and `diag/pool_summary.py`.
4. Success means s301 and s703 pick their 0.97 answers. s505 and s403 then show whether the
   right answer enters the pool.
5. Rerun the 8 regression and 10 passing bench cases with the best setting. They must not drop
   by more than 0.005.

If the residual weighting is not enough, the larger change is to decompose against the
observations directly: put the species weights into the flow MILP, or refit with a noise model
before choosing among the pool's exact fits.

### Step 3. Long-read paths (cheap check first) [done 2026-10-29: a bug, fixed; paths do not help s301 or s703]

- The v2 `ctrl_dup_s702` and `s703` subgraphs load **0 paths**; the v1 subgraphs loaded 3 and 4.
  Find out why `rck_inputs.py`'s `paths.tsv` came out empty. Compare `union/paths.tsv`, and the
  path filter in `rck_inputs.py`, between `simulations/subgraphs/bench_nested/ctrl_dup_s702` and
  the v2 one.
- Paths are the one data source that can tell exact-fit alternatives apart, as in s703. Then try
  `--path-weight` > 0.5 or `--path-fraction` on s703 and s301.

### Step 4. Ploidy selection in `rck_inputs.py` (amplicon-allele-analysis)

- Log the full score table for every bench case (already in `scale.json` → `alternatives`).
  Look at what separates the wrong picks from the right ones: s403, and v1 s501, s502 and s406.
- Options:
  - prefer ploidy 2 unless a higher class wins by a margin;
  - add a check that segment CN at depth per copy agrees with junction CN at θ;
  - with paired normal depth, use it directly.
- A segment-vs-junction ratio of ≈ 2 (`diag/scale_ratio.py`) flags the error at once.
- Rebuild only the affected subgraphs, then rerun RCK on them.

### Step 5. Junction reads per copy

- Make `rck_inputs.py`'s θ fit robust before using it: it fits from background extremities.
  Check its inputs on s505 and s403, where it is wrong.
- Alternatively, let RCK fit θ jointly in the flow: one free scale factor on junction
  observations.
- Test on s502 (needs ≈ 5.0) and s505 (needs ≈ 6.5).

### Step 6. Benchmark bookkeeping [done 2026-10-29]

- `collect.py`: add a `scale` failure stage. Flag it when `depth_per_copy` deviates more than 30%
  from the median of the cases, or when segment and junction scale disagree.
- `collect.py`: use λ = 0.2 as the reported λ in `pipeline.tsv` (now 0.5 by the plateau rule). Or
  change the rule so that one non-improving case at the grid edge cannot pick the λ.
- `verify_subgraph.py`: score segment calls against the summed haplotype composition of every
  species covering a segment, not the nested species' label alone.
- `simulations/plot_rck_solution.py` still defaults to the removed `runs/demo3_clone1_ecdna`.
  Point it at demo5 or take the run folder as an argument. `subgraph_construct.md` still has a
  "Results on demo3" section.

### Step 7. After the fixes

- Rerun all 22 cases (`sbatch --array=8-21%2 bench_nested.sbatch` with the v2 environment, and the
  regression cases too if RCK changed).
- `python collect.py`, then update `README.md`, both CHANGELOGs, and the λ recommendation.

## Housekeeping

- amplicon-allele-analysis `.git` holds about 380 MB of loose objects. Most of it is the removed
  demo3/demo4 fixtures, which remain in history. `git gc` packs them and frees home quota.
  Removing them from history needs a rewrite and a force-push of `yubi_dev`; decide before doing
  it.
- RCK's untracked `results/bench/` (179 MB), `logs/` and `ecbench_*.out` in the repo root are run
  outputs. Move them to scratch if home quota gets tight.
- Diagnostic outputs of this session: `/n/fs/scratch/$USER/bench_nested_v2/diag/rck/` (s403 with
  ploidy 2; s502, s505, s301 and s703 with fitted θ; Slurm 32005476 and 32005494). Their scripts
  are in `diag/` here.

## Progress (2026-10-29)

- **Step 1.** v2 configs (with `seeds.tsv`) and subgraphs replace the tracked ones in
  amplicon-allele-analysis (v1-only files removed with `git rm`); v2's tables, `comparison.md`
  and `comparison.tsv` are here; `runs` → `bench_nested_v2/rck_runs`, amplicon-allele-analysis
  `simulations/runs/bench_nested` → `bench_nested_v2/runs`. `git gc` in both repositories (no
  history rewrite): amplicon-allele-analysis 380 MB loose → 212 MB packed, RCK 115 MB → 3.5 MB.
- **Step 2, code.** `rck-ecdna --adj-residual-weight`, `--seg-residual-scale`,
  `--residual-mode relative --residual-sigma`, and `--rank-residual sum`, all off by default
  (`docs/ecDNA.md`, "Noise-tolerant residual"). The ranking puts residual first, so lower weights
  alone change only the pool (and `--pool-gap` pruning can then drop the exact fit), not the
  choice among what is kept: the diagnostic runs use `--rank-residual sum`.
  - What the existing pools say: in s301 and s703 **every** pool solution has residual 0 (exact
    fits of the integer flow). s703's 0.971 answer (the walk at 16 + a 2-copy leftover) and its
    rank 1 (two circles) both have 2 significant species and residual 0; they differ only in
    steps (31 against 17) and the tie-break. s703's flow has the B insertion and both its
    junctions at 16 against 18–19 elsewhere, so a single walk leaves about 9–11 weighted copies of
    residual: it needs weight ≲ 0.1 to beat one species. Expect 0.5, 0.25 and relative σ 0.1 not
    to fix s703; s301 may differ. (Confirmed by the runs, below.)
  - **Runs** (Slurm 32025976, 22 tasks; `diag/rck_noise/`: `SUMMARY.md`, `noise_summary.tsv`, the
    LIST files and each run's small outputs). **No setting passes.**
    - s301 (0.693) and s703 (0.712) keep their rank 1 at every residual setting and every path
      setting.
    - s505 stays at 0.765 at every setting (w1, w0.5, w0.25, rel0.1), the best in its pool: no
      good answer enters.
    - s403 (ploidy 2) stays at 0.545 at every setting. The pool's best is 0.729 (4 significant
      species), never picked. The w0.25 run hit Gurobi's 1800 s limit in its pool search, so its
      row is not filled; its archived outputs say the same.
    - Tasks 15 and 18–21 were recorded on 2026-10-31.
  - **What was learnt:**
    - The 0.97 answers lose on **leftovers counted as significant species**, not on residual.
      s703's walk plus a 2-copy leftover counts 2 species, like the two circles; s301's answer
      counts 4 against 3.
    - **The observations cannot separate them.** Refitted to the observations, every s703 pool
      solution has the same misfit (`data_objective` 5.097), and s301's 0.97 answer is better by
      0.002 out of 5.1. A refit before choosing, or species weights in the flow MILP, therefore
      cannot fix s703, so the larger change NEXT_STEPS proposed is not the next step for these
      cases.
    - Raising `--min-species-weight-frac` to 0.15, checked offline on the saved pools, fixes s703
      (0.971) but not s301, and makes s505 worse (0.765 → 0.53): a global threshold is not it
      either.
    - **Completed record (2026-10-31).** No residual weight changes rank 1 or adds a better pool
      answer in s403 or s505. The reading above holds for s301 and s703. **s505 and s403 fail
      before the ranking:** their pools hold no good answer (s505 best 0.765 = rank 1; s403 best
      0.729).
      - s505's flow is wrong where its observations are: depth on about 300 kb of haplotype B
        reads 12–27% under the truth (the simulator's correlated depth bias).
    - **The leftover prior proposed as step 2b does not work** (offline recount of all 22 saved
      pools; "Assessment (2026-10-31)"). Counting leftovers against the heaviest *overlapping*
      species gives the same answers as the global threshold, because every leftover in these
      pools overlaps the heaviest species. s505 and s502 break because the default threshold is
      max(1, 0.1 × heaviest) = 1 copy whenever the heaviest species has ≤ 10 copies: there their
      structurally needed 1-copy species only count because of that floor.
    - **RCK already finds the two right answers that lose:** s703's 0.971 answer is the first
      same-fit alternative of its rank 1, and s301's 0.973 answer is the second. Its ambiguity
      report records both, but it reports only counts, not the alternatives' walks.
- **Step 3.** The empty `paths.tsv` was a bug in amplicon-allele-analysis `rck_inputs.py`: a path
  was dropped unless every simulator segment boundary on it, read ends included, was a subgraph
  node. Under inherited labels some former haplotype-switch junctions are plain reference
  continuations without a node, which emptied s702 and s703; most paths were dropped in every
  case. Fixed (`clip_path`). Rebuilt s702, s703 and s301 into
  `/n/fs/scratch/$USER/bench_nested_v2/subgraphs_paths/` (graphs byte-identical; paths with ≥ 2
  junctions: s702 0, s703 1 (3 reads), s301 3 (5 reads)), so paths can only weakly help.
  **The committed subgraphs predate the fix**, and paths are on by default (`--path-weight 0.5`):
  rebuilding them changes RCK's input for every case. Do it before step 7's rerun.
- **Step 6.** `collect.py`: `scale` failure stage (v2 s403 is now `scale`; s104 passes but is
  flagged, junction scale 0.675) and the 2-case λ rule (v2: 0.2; the first benchmark would give
  0.05). `verify_subgraph.py` scores segment calls against summed haplotype copies (v2: 34 wrong
  → 0; the committed `verify_report.txt` and `subgraph_checks.tsv` still hold the old counts,
  refreshed by the next subgraph rebuild). `plot_rck_solution.py` takes the run folder.
- `--rck-ploidy` is a `subgraph_construct.py` / `rck_inputs.py` option, not an `rck-ecdna` one:
  s403 at ploidy 2 uses the rebuilt subgraph `bench_nested_v2/diag/s403_ploidy2/union`.

## Next steps (revised 2026-10-29; superseded by "Assessment (2026-10-31)", kept for the record)

In order. Steps 4, 5 and 7 are as planned above; step 2 continues as 2b.

1. **Finish the step 2 record** (array 32025976 tasks 15 and 18–21, see "Check first"). Fill the
   pending s403 and s505 rows of `diag/rck_noise/SUMMARY.md`
   (`python diag/noise_summary.py diag/rck_noise/ordered.list --archive diag/rck_noise/runs`).
   If s403 or s505 at some setting gets a better answer into the pool, note it there.
2. **Step 2b: a noise-aware leftover prior** (RCK, opt-in). Count a species as significant only if
   its weight exceeds the copy-number noise of the species it overlaps (σ × the heaviest overlapping
   species' weight), not 10% of the heaviest species overall.
   - Test it offline first: recount significant species in the saved pools and re-rank them with
     `rank_pool`, as in `SUMMARY.md`.
   - First find out why s505 loses at a higher global threshold.
   - Then run it on the four cases, and on the 8 regression and 10 passing bench cases (no drop
     above 0.005).
3. **Rebuild the 14 committed subgraphs with the `rck_inputs.py` path fix** (amplicon-allele-analysis
   `bench_nested.sbatch`). Paths are on by default, so RCK's input changes for every case, and the
   rebuild refreshes `verify_report.txt` and `subgraph_checks.tsv` to the new segment-call scoring.
   Do this before step 7.
4. **Step 4** (ploidy selection in `rck_inputs.py`) for s403, which `pipeline.tsv` now labels
   `scale`, and **step 5** (θ) for s502 and s505.
5. Optional: raise the long-read count in the s703 and s301 configs, to see whether paths can
   separate their answers when there are more of them.
6. **Step 7**: rerun all 22 cases (and the regression cases if RCK changed), `python collect.py`,
   then update `README.md`, both CHANGELOGs and the λ recommendation.


## Assessment (2026-10-31)

An evaluation of "Next steps (revised 2026-10-29)" against the completed step 2 diagnostic and the
saved pools. It uses no new runs: only offline recounts of the 22 bench pools on scratch
(`bench_nested_v2/rck_runs/*/pool`) and the four `diag/rck_noise/*_w1` pools. The bench mean
composite is 0.928 at λ = 0.2 (22 cases), and a case's change moves the mean by Δ/22. Each fixed
failure is worth +0.012 to +0.017. All four together are worth at most about +0.05.

New offline scripts (in `diag/`, opt-in; no RCK or amplicon-allele-analysis default changes):
- `leftover_recount.py` re-flags every pool solution's species under a leftover rule and
  re-ranks with RCK's own `rank_pool`. With the default rule (`global 0.1`) it reproduces every
  case's recorded rank 1 (mean 0.9283) and SUMMARY.md's global-threshold table. Output:
  `diag/rck_noise/leftover_recount.tsv`.
- `obs_vs_truth.py` compares segment depth and junction reads with the truth copy numbers.
  Output: `diag/rck_noise/obs_vs_truth.tsv`.

### 1. Step 2b, the overlap-relative leftover prior: tested offline, drop it

**Rule tested.** A species is low weight (not counted in the ranking) when its integer weight is
below max(1, σ × w), where w is the heaviest *other* species it overlaps. A species that overlaps
nothing stays significant. **Overlap**, as used for the headline numbers: the two species share at
least one *phased* segment, i.e. the same chromosome, interval and haplotype label in their walks
(`hap / 0` in the TSV). Three other definitions were also run:
- haplotype-blind (`total`: interval only);
- shared segments covering more than 50% of the lighter species' length (`/ 0.5`);
- the rule on top of the default 10% (`overlap+global`).

Species that the weight refit dropped have no saved walk. They are assumed to overlap every
species (`--pruned all`). Rerunning with `--pruned none` changes 3 cells (noted below), and no
conclusion.

Rank-1 mean composite at λ = 0.2:

| case | default (global 0.1) | overlap σ 0.1 | σ 0.15 | σ 0.2 | best in pool |
|---|---|---|---|---|---|
| **s301** | 0.693 | 0.693 | 0.693 | 0.703 | 0.973 |
| **s703** | 0.712 | 0.712 | **0.971** | **0.971** | 0.971 |
| **s505** | 0.765 | 0.765 | 0.528 | 0.532 | 0.765 |
| **s403** (bench, ploidy 6) | 0.592 | 0.592 | 0.592 | 0.592 | 0.592 |
| **s403** (ploidy 2, `rck_noise` w1) | 0.545 | 0.545 | 0.506 | 0.507 | 0.729 |
| s502 | 0.955 | 0.955 | 0.513 | 0.513 | 0.960 |
| s402 | 0.978 | 0.978 | 0.978 | 0.801 | 0.978 |
| s104 | 0.962 | 0.962 | 0.962 | 0.949 | 0.962 |
| demo8 | 0.987 | 0.987 | 0.987 | 0.836 | 0.992 |
| demo7_totalonly | 0.976 | 0.978 | 0.969 | 0.969 | 0.988 |
| the other 13 passing cases (s702, s601, s603, s103, s201, s202, s304, demo3, demo4 ×2, demo5, demo6, demo7_blockflip) | unchanged | unchanged | unchanged | unchanged | |
| **mean, 22 cases** | **0.928** | 0.928 | 0.909 | 0.894 | |
| cases > 0.005 below default | | 0 | 3 | 6 | |

- The haplotype-blind rule gives 0.910 at σ 0.15 and 0.895 at σ 0.2, case for case the same as
  the global threshold at 0.15 and 0.2 in the same recount.
- The 50% overlap variant is no better (0.908 at σ 0.15).
- `--pruned none` changes demo7_totalonly at σ 0.1 (0.864), and s301 (0.665) and s104 (0.962)
  at σ 0.2.

**Why the overlap rule is no better than the global one.** In every pool tested, each leftover
overlaps the heaviest species: nested species share the outer species' segments, and leftovers
are pieces of the heavier walks. So "σ × heaviest overlapping" = "σ × heaviest". The haplotype-aware
rule differs only where a species sits wholly on the other haplotype (demo7_totalonly).

**Why s505 (and s502) get worse at 0.15.** It is not overlap. It is the 1-copy floor:
- The default threshold is max(1, 0.1 × heaviest). For a heaviest species of ≤ 10 copies it is
  1, so then *every* species counts as significant.
- s505's pool #56 (10 / 4 / 1 / 1) therefore has 4 significant species at 0.1 and loses. At
  0.15 (threshold 1.5) its two 1-copy species drop out and it wins, with 2 against rank 1's
  3 (12 / 7 / 3 / [1]).
- Those 1-copy species are not noise: they are 7.0 and 7.5 Mb walks with 7 and 9 junctions,
  and the 10-copy species merges both true species (0.528).
- s502 is the same with pool #7 (10 / 1 / 1; 0.955 → 0.513).

Any rule that makes a 1-copy species light enough to drop does this, overlap-relative or not.

**A structural variant was also tested** (`novel σ`): a species counts as a leftover only if it
adds no junction that heavier species lack, and weighs < σ × the heaviest.
- It separates s301's and s703's leftovers (both are subsets of the heavy walk's junctions) from
  real nested species. It fixes s301 (0.973) and s703 (0.971) at σ ≥ 0.3, and moves s403 (ploidy
  2) to 0.595.
- But it breaks s502 (0.514), s505 (0.532) and s402 (0.824). At σ 0.5 it also breaks s201 and
  demo4_amp0, and at σ 1 also s103, s202 and both demo7 cases. Mean 0.916 at σ 0.3.
- It loses for the same reason as the overlap rule. In s502 and s505, solutions with one heavy
  walk that carries every junction, plus 1-copy pieces, become cheap.

| step 2b | |
|---|---|
| fixes | s703 only (σ ≥ 0.15); s301 not below σ 0.3 and then only with the structural variant |
| bench mean | 0.928 → 0.909 (σ 0.15) / 0.894 (σ 0.2); best variant 0.916 |
| risk to passing | s502 drops at every σ that helps; demo8, s402 and s104 at σ 0.2 |
| cost | small code; offline test done; a run would cost ≈ 1 h of licence time for nothing |
| **verdict** | **drop.** A prior on weights alone cannot tell s703's 2-copy leftover from s502's structurally needed 1-copy species |

### 2. The shared root cause the plan misses

The four failures are two different problems, and neither is a ranking-weight problem.

**(a) s301 and s703: the data cannot tell the answers apart. Report them as ambiguous; don't
add a prior.**
- Both right answers fit the integer flow exactly. After the refit they fit the observations
  as well as rank 1, to 0.002 or better (s703: 5.097 for every pool solution).
- The long reads cannot help either.
  - Simulated reads average 15 kb (N50 20 kb, `l5_reads` in the configs).
  - s703's two answers (circles at 18 and 16, or one walk at 16 plus a 2-copy copy of the
    A part) differ only in how entries and exits pair across the shared 2 kb + 61 kb stretch
    (chr1:9838001–9900999). Only a read longer than ≈ 63 kb would separate them.
  - So step 5's "more long reads" cannot help as stated: more reads of the same length add
    nothing. Only an ultra-long-read experiment (reads > 65 kb) would show whether long reads
    can resolve it.
- **The prior.** The only prior that would be defensible here is "a light species that adds no
  junction is coverage unevenness of the heavier walk". It is the structural variant above, and
  it costs s502, s505 and s402. With 2 seeds per scenario there is no evidence that any prior
  is right more often than it is wrong.
- **What RCK already has.** The ambiguity report finds both right answers:
  - s703's 0.971 answer (16 / 2) is the **first** same-fit alternative of rank 1 (4 in all);
  - s301's 0.973 answer (11 / 4 / 3) is the **second** (19 in all).
  - Both decisions are `decided_by = step_penalty`.
- **What is missing is reporting.**
  - `ambiguous` is true for 17 of 22 cases, so it says little.
  - `decided_by = step_penalty` marks 6 of the 22 bench cases: the decomposition failures s301,
    s505 and s703, and three passing cases (s603, demo8 and s202). It also marks s403 at ploidy 2
    (`rck_noise` w1).
  - The alternatives' walks are not written out.
  - A "top-k of the ranked list" would not help: at λ = 0.2, `pool_keeps` prunes s301's answer
    from the ranked list, and s703's is ranked 4th.

**(b) s505 and s403: the pool has no good answer, because of the flow and the subgraph.** Ranking
cannot fix these.
- **s505.** On chr3:121.68–121.98 Mb (haplotype B), depth reads 22–27 copies against a truth of
  30.5, and the flow follows it (`flow_vs_truth.py`: 14 of 21 segments off by > 1.5 copies, the
  most of any case).
- **Would junction reads fix it?** In s703 the B insertion's junctions read 18.6 and 18.2 copies
  (truth 19.0) while the flow sets 16 from depth. But over the bench, junction reads are not
  more accurate than depth (`obs_vs_truth.tsv`). Mean relative error, junctions against
  segment depth:
  - s703 3% against 12%; s502 8% against 12%;
  - s505 7.9% against 7.2%; s402 10.5% against 7.0%.
  - So raising `--alpha-novel` (already on the CLI) would trade cases, like the leftover prior.
- **s403 at ploidy 2.** The best pool answer is 0.729, and its flow error (9% of segment copies,
  like s301's) is the bench's largest.

What would show whether (a) and (b) come from the simulated depth bias (`sigma` 0.1, 200 kb
correlation) or from a model error: rerun the failing seeds with the bias off (run A below). If
they pass at σ = 0:
- they are identifiability limits of this data;
- the right response is reporting (step A1 below) and data requirements, not more priors.

### 3. Order of steps 3–5

| step | fixes | Δ bench mean | risk to the 18 passing | cost | verdict |
|---|---|---|---|---|---|
| 3. rebuild the 14 subgraphs with the path fix | none expected: paths did not move s301 or s703, and the reads are too short to separate their answers | ≈ 0, unknown sign | changes every case's RCK input (paths are on by default) | CPU only (`bench_nested.sbatch`, 7 tasks ≤ 45 min); then the 22-case rerun, ≈ 4350 s of RCK wall time (≈ 1.2 licence-hours) | **proceed, merged with 4 into one re-baseline** |
| 4. ploidy selection | s403's scale (`pipeline.tsv` stage `scale`) | **−0.002 alone** (0.592 → 0.545); +0.006 if s403 then reached its pool best | low: see the offline check below | small code in `rck_inputs.py`; rebuild with 3 | **proceed, with the paired-normal rule** (superseded: an explicit `--rck-ploidy 2`; see the correction at the end) |
| 5. θ fit | s502 (0.955 → 0.963) | +0.0004 | s505: the fitted 4.93 gives 0.600 (−0.0075); s403: fitted 2.12 | code + rebuild + reruns | **defer**: the bench fixes the library size, so θ = 6.5 is right by construction; it matters for real data only |
| 5 (optional). more long reads for s703 and s301 | none (reads too short; section 2a) | 0 | none | simulation + ≈ 2 min of RCK | **drop as stated**. Replace with an ultra-long-read experiment only if the data-requirement question matters |

**Offline check for step 4 (`scale.json` of all 14 tracked subgraphs).** *(Superseded: this rule is a hidden diploid prior; see "Agreed order … and a correction" at the end.)* The simulation has a
paired normal sample (`all_samples.sim_normal`).
- At the right class, the tumour's depth per copy matches the normal's within 5% in every case
  (0.98–1.05).
- s403's chosen class 6 is at 0.48× the normal. Its class 2 would be at 1.02×.
- "Choose the ploidy class whose depth per copy is closest to the paired normal's" picks the
  right class in 14 of 14 cases.
- A margin rule ("prefer 2 unless beaten by m") is fragile: s502's right pick wins by 0.002 and
  s403's wrong one by 0.008.
- For real data the normal's depth per copy must be scaled by the tumour/normal library-size
  ratio; the simulation uses equal library sizes (`normal_library_size: null`). Without a
  normal, keep `scale_ratio.py`'s segment/junction check as the flag.

**Does any of them change step 2b's inputs?** Yes:
- step 3 changes every case's input (paths);
- step 4 changes s403's segment copy numbers by 2×;
- step 5 would change every junction observation.

So no ranking prior should be tuned before steps 3 and 4 are rebuilt. With 2b dropped, the point
is moot for 2b, but it holds for any later ranking change (A2 below).

**A metric problem shown by step 4.** The composite ranks s403's wrong-scale answer (segment CN
2× truth, 0.592) above the right-scale one (0.545). The CN term is too weak to see a 2× scale
error. Report it separately (`pipeline.tsv` already has the `scale` stage), and do not read
step 4's composite change as a loss.

### 4. Other things that look more valuable

- **Ambiguity reporting (RCK output only; rank 1 unchanged).**
  - Write the walks and weights of the first few same-fit alternatives when `decided_by` is a
    prior, so the user gets "11 / 3 + leftover, *or* 11 / 4 / 3".
  - In `collect.py`, add a set metric next to the rank-1 one: the best composite over rank 1
    and its first 2 same-fit alternatives.
  - On the saved pools this gives s301 0.973 and s703 0.971, a set-metric mean of at least
    0.952 against 0.928. It cannot lower any case, since rank 1 is unchanged.
  - This is the honest answer for cases the data cannot decide, and it can be tested offline from
    `pool.json` without Gurobi.
- **Too few seeds.**
  - With 2 seeds per scenario, each failing case is one seed. Every candidate prior here fixes 1–2
    cases and breaks 2–3, so the bench cannot tell a good prior from a bad one.
  - λ = 0.2 itself rests on one case (s402).
  - The per-case no-drop limit of 0.005 is fine as a guard. The problem is that the gains
    being judged come from single cases.
  - Before tuning any prior, add seeds to the scenarios with failures: nest_samehap_3to1,
    nest_cross_1to1 and ctrl_dup, to about 6 each.
- **A missing scenario.** No scenario has a light true nested species (6:1 or more on the same
  haplotype). Every leftover-type prior would kill such a species, and the bench cannot see that
  risk. Add `nest_samehap_6to1`.
- **A missing control.** No case has the depth bias off, so noise and model errors cannot be told
  apart (run A).

### Revised plan (2026-10-31), in order

1. **[done] Step 2 record** (`diag/rck_noise/SUMMARY.md`).
2. **A1. Ambiguity reporting** (RCK, opt-in output plus a `collect.py` column; no Gurobi to test).
   - Write the top same-fit alternatives' walks when the choice is `decided_by` a prior.
   - Score the "rank 1 + first 2 alternatives" set in `collect.py` from the saved pools.
   - Expected: rank-1 mean unchanged (0.928), set metric ≥ 0.952, no risk to passing cases.
3. **Run A: σ = 0 control** (simulation and subgraphs on CPU, then about 45 min of licence time).
   It decides whether s301, s703 and s505 are noise limits or model errors, and so whether any
   further RCK model work (A2) is worth it. Commands below.
4. **Re-baseline: steps 3 + 4 together.**
   - Paired-normal ploidy rule in `rck_inputs.py` (opt-in flag first; it becomes the bench's flag).
     *Superseded: an explicit `--rck-ploidy 2` for the bench; see the correction at the end.*
   - Rebuild all 14 subgraphs with the path fix and the rule.
   - Then step 7: rerun all 22 cases and `collect.py`.
   - About 1.2 licence-hours.
5. **Bench growth before any tuning.** About 4 more seeds for nest_samehap_3to1, nest_cross_1to1
   and ctrl_dup, and a `nest_samehap_6to1` scenario. Simulation on CPU, then RCK, about 2–3
   licence-hours (nest_cross_1to1 pools take 5–25 min each).
6. **A2 (only if run A shows model errors).** Depth-bias-aware flow weighting: segments within one
   correlation length count as one observation, or a probe of `--alpha-novel` (run B). Judge it
   on the grown bench.
7. **Deferred or dropped.** Step 2b (dropped; offline test above). Step 5 θ fit (deferred to
   real-data work). More long reads at 15 kb (dropped). Species weights in the flow MILP or a
   noise-model refit (dropped for s301/s703 by the identical refit misfit; revisit for s505 only
   after run A).

### Runs to answer the open questions (not submitted)

**Run A, σ = 0 control.** The same seeds, with the simulator's correlated depth bias off.
- With `sigma: 0` the bias multiplier is exp(0) = 1, and the random draws are still made
  (`src/emit_counts.py::_coverage_bias`), so walks and reads are unchanged.
- Check that `truth/` matches the v2 run before reading the results.
- Step 1 runs on CPU only (no Gurobi), from amplicon-allele-analysis `simulations/`:

```bash
export BENCH_NESTED_ROOT=/n/fs/scratch/$USER/bench_nested_sigma0
mkdir -p $BENCH_NESTED_ROOT/configs logs
for c in nest_samehap_3to1_s301 ctrl_dup_s703 nest_cross_1to1_s505 nest_cross_1to1_s502; do
  sed 's/^  sigma: 0.1$/  sigma: 0.0/' configs/bench_nested/$c.yaml > $BENCH_NESTED_ROOT/configs/$c.yaml
done
grep -c '^  sigma: 0.0$' $BENCH_NESTED_ROOT/configs/*.yaml          # 1 each
sbatch --job-name=sigma0-sim --cpus-per-task=1 --mem=4G --time=03:00:00 --output=logs/%x_%j.out \
  --export=ALL,BENCH_NESTED_ROOT --wrap 'export PATH=$HOME/.conda/envs/gurobi-env/bin:$PATH; python -c "
import os, bench_nested as b
for sc, seed in [(\"nest_samehap_3to1\", 301), (\"ctrl_dup\", 703), (\"nest_cross_1to1\", 505), (\"nest_cross_1to1\", 502)]:
    b.build_subgraph(b.simulate(sc, seed), os.path.join(b.SUBGRAPHS, b.run_name(sc, seed)))"'
```

Step 2 needs Gurobi: 4 runs, s505 up to 25 min. s502 is the passing control from the same
scenario. LIST file `/n/fs/scratch/$USER/bench_nested_v2/diag/rck_sigma0/sigma0.list`:

```
nest_samehap_3to1_s301_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/nest_samehap_3to1_s301/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/nest_samehap_3to1_s301 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_sigma0/nest_samehap_3to1_s301_sigma0
ctrl_dup_s703_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/ctrl_dup_s703/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/ctrl_dup_s703 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_sigma0/ctrl_dup_s703_sigma0
nest_cross_1to1_s505_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/nest_cross_1to1_s505/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/nest_cross_1to1_s505 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_sigma0/nest_cross_1to1_s505_sigma0
nest_cross_1to1_s502_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/nest_cross_1to1_s502/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/nest_cross_1to1_s502 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_sigma0/nest_cross_1to1_s502_sigma0
```

```bash
cd /n/fs/scratch/$USER/bench_nested_v2/diag/rck_sigma0 && mkdir -p logs
LIST=$PWD/sigma0.list DIAG_POOL_LAMBDA=0.2 DIAG_LAMBDA0=0 sbatch --array=0-3%2 --time=02:00:00 ~/RCK/results/bench_nested/diag/rck_diag.sbatch
python ~/RCK/results/bench_nested/diag/noise_summary.py sigma0.list   # when done
```

Read it as follows:
- s301 or s703 picking its 0.97 answer at σ = 0 means the failure is the depth bias (a noise
  limit).
- Still failing means a prior or model problem in RCK.
- s505's pool gaining a good answer means its failure is the biased flow.

**Run B (optional, only after run A shows model errors): flow junction weight.** The existing v2
subgraphs at `--alpha-novel 1.0` (default 0.5). s402, s201 and s502 are guards: their junctions are
noisier than, or about as noisy as, their depth. LIST file
`/n/fs/scratch/$USER/bench_nested_v2/diag/rck_flow/alpha.list`, submitted with:

```bash
cd /n/fs/scratch/$USER/bench_nested_v2/diag/rck_flow && mkdir -p logs
LIST=$PWD/alpha.list DIAG_POOL_LAMBDA=0.2 DIAG_LAMBDA0=0 sbatch --array=0-5%2 --time=02:00:00 ~/RCK/results/bench_nested/diag/rck_diag.sbatch
```

The lines are:

```
ctrl_dup_s703_an1 /n/fs/scratch/ym3292/bench_nested_v2/subgraphs/ctrl_dup_s703/union /n/fs/scratch/ym3292/bench_nested_v2/runs/ctrl_dup_s703 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_flow/ctrl_dup_s703_an1 --alpha-novel 1.0
nest_cross_1to1_s505_an1 /n/fs/scratch/ym3292/bench_nested_v2/subgraphs/nest_cross_1to1_s505/union /n/fs/scratch/ym3292/bench_nested_v2/runs/nest_cross_1to1_s505 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_flow/nest_cross_1to1_s505_an1 --alpha-novel 1.0
nest_samehap_3to1_s301_an1 /n/fs/scratch/ym3292/bench_nested_v2/subgraphs/nest_samehap_3to1_s301/union /n/fs/scratch/ym3292/bench_nested_v2/runs/nest_samehap_3to1_s301 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_flow/nest_samehap_3to1_s301_an1 --alpha-novel 1.0
nest_cross_1to1_s502_an1 /n/fs/scratch/ym3292/bench_nested_v2/subgraphs/nest_cross_1to1_s502/union /n/fs/scratch/ym3292/bench_nested_v2/runs/nest_cross_1to1_s502 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_flow/nest_cross_1to1_s502_an1 --alpha-novel 1.0
nest_cross_3sp_s402_an1 /n/fs/scratch/ym3292/bench_nested_v2/subgraphs/nest_cross_3sp_s402/union /n/fs/scratch/ym3292/bench_nested_v2/runs/nest_cross_3sp_s402 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_flow/nest_cross_3sp_s402_an1 --alpha-novel 1.0
nest_samehap_2to1_s201_an1 /n/fs/scratch/ym3292/bench_nested_v2/subgraphs/nest_samehap_2to1_s201/union /n/fs/scratch/ym3292/bench_nested_v2/runs/nest_samehap_2to1_s201 /n/fs/scratch/ym3292/bench_nested_v2/diag/rck_flow/nest_samehap_2to1_s201_an1 --alpha-novel 1.0
```

### Agreed order (2026-10-31, later) and a correction to the ploidy recommendation

You chose to skip ahead to more seeds, with the re-baseline (old plan steps 3 + 4) done alongside
it rather than skipped. The order is now:

1. **c. Re-baseline the existing 22 cases (run now).** Rebuild the 14 bench subgraphs with the
   current `rck_inputs.py`, which has the `clip_path` fix, and the ploidy handling below. Build
   into a scratch root, so the tracked subgraphs and the v2 runs stay as they are until reviewed:
   - root `/n/fs/scratch/$USER/bench_nested_v2/rebuild/`, with `configs/` (a copy of the tracked
     ones), `runs` (→ the v2 simulations, unchanged) and `subgraphs/`;
   - RCK outputs in `/n/fs/scratch/$USER/bench_nested_v2/rebuild/rck_runs/`.
   Then run RCK on all 22 cases. The 8 regression cases rerun too, so that the new baseline is one
   self-consistent set of runs.
2. **d. More seeds.**
   - Simulate about 4 more accepted seeds each for nest_samehap_3to1, nest_cross_1to1 and
     ctrl_dup, and add a `nest_samehap_6to1` scenario. Simulation needs no Gurobi.
   - Build their subgraphs with the same settings as in c, and run RCK on them.
3. **b. σ = 0 control.** Optional; run it alongside d if licence time allows (commands above).
4. **a. Ambiguity reporting.** Deferred until after d.
5. Tune nothing until c and d are both in.

**Correction: the paired normal cannot choose the ploidy class.** The "Assessment" section above
recommends choosing the class whose depth per copy is closest to the paired normal's (14 of 14
right). That rule is a hidden diploid prior:
- In `rck_inputs.py::_fit`, a class's depth per copy is `med / cn1` and its cell ploidy is
  `Σ rdr × cn1 × w`. Their product is the tumour's mean depth, the same for every class.
- So with matched libraries, depth per copy only says how many genome copies the reads are
  spread over. It cannot tell a pure diploid tumour from an impure tetraploid one: that is
  exactly what the lattice and BAF scores are for.
- "Closest to the normal" therefore means "closest to cell ploidy 2". It is right 14 of 14 times
  only because every simulated tumour is near-diploid (tumour ploidy 1.85–2.0). On a
  whole-genome-doubled tumour it would be wrong.
- The bench has no such tumour, so it cannot test any ploidy rule.

**Ploidy handling for c and d: an explicit `--rck-ploidy 2`.**
- All bench subgraphs are built with `--rck-ploidy 2`, set through the new opt-in
  `BENCH_NESTED_RCK_PLOIDY=2` of amplicon-allele-analysis `bench_nested.py`. Without it,
  `bench_nested.py` is unchanged and still lets `rck_inputs.py` fit the class.
- `scale.json` still records every class's fitted score (`alternatives`), so how often the auto
  fit would have been wrong stays measurable. In v2 it was 1 case of 14 (s403).
- This makes the bench measure the subgraph-to-RCK decomposition, with the copy-number scale
  given, which is what steps 2b and later tune.
- **Ploidy inference becomes a separate item.** It needs its own scenarios: tetraploid and
  impure tetraploid tumours, and low purity. Without them no rule can be judged, so
  `rck_inputs.py`'s default auto fit stays as it is.
- `collect.py`'s `scale` stage still flags any depth-per-copy or segment/junction-scale outlier.

**Step c, status (2026-10-31, 14:5x):**
- **Subgraphs rebuilt.** All 14 bench subgraphs are rebuilt in
  `/n/fs/scratch/$USER/bench_nested_v2/rebuild/subgraphs/`, from
  amplicon-allele-analysis `simulations/` with:
  ```bash
  BENCH_NESTED_ROOT=/n/fs/scratch/$USER/bench_nested_v2/rebuild BENCH_NESTED_RCK_PLOIDY=2 python bench_nested.py rebuild
  ```
  - Subgraph-stage failures: none.
  - Against the tracked v2 subgraphs, `nodes`, `edges`, `segment_blocks`, `phase_blocks` and
    `cn_profile` are byte-identical in 13 cases.
  - `paths.tsv` grows in all 14: from 0–10 to 6–23 paths.
  - s403 changes copy numbers only (auto fit class 6 → set 2). It is identical to the 2026-10-27
    `diag/s403_ploidy2` subgraph except for paths.
  - The auto fit would have chosen class 2 in the other 13 cases.
- **RCK submitted.** Slurm array **32041913** (`--array=0-21%2 --time=02:00:00`), submitted from
  here with:
  - `BENCH_NESTED_ROOT` as above, and
    `BENCH_NESTED_OUT=/n/fs/scratch/$USER/bench_nested_v2/rebuild/rck_runs`;
  - logs in `/n/fs/scratch/$USER/bench_nested_v2/rebuild/logs/`.
- **When it finishes:**
  - check it with `sacct -j 32041913 -X -o JobID%14,State,Elapsed`;
  - then run `collect.py` with the same two variables (writing to a scratch folder, not over this
    folder's v2 tables) and compare the result with v2;
  - only then copy the rebuilt subgraphs into amplicon-allele-analysis
    `simulations/subgraphs/bench_nested/`, and their RCK tables here.

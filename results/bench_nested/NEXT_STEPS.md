# bench_nested: next steps (written 2026-10-27)

> **Resume here (2026-11-12):** the current state and the ordered steps to continue are in
> "Resume here (2026-11-12): state and what to do next" at the end of this file. Read it first.

State at the end of the 2026-10-27 session, and a plan for what comes next; progress of the
2026-10-29 session is marked **[done 2026-10-29]** and in "Progress (2026-10-29)" at the end. The
evidence is in `README.md`: the "bench_nested_v2" and "Follow-up checks" sections. The full
old-vs-new comparison is in `comparison.md` and `comparison.tsv` here.

Scratch purges files older than two weeks. Only run outputs remain there (`bench_nested_v2/runs`,
`rck_runs`, `diag`, `subgraphs_paths`); configs, subgraphs and tables are in the repositories, and a
copy of the v2 configs, subgraphs and `rck_results` is in `~/backups/bench_nested_v2_2026-10-29.tar.gz`.

> **2026-11-01:** steps c and d are done and the grown bench (40 cases) is adopted. The current
> plan is "Grown bench: failure analysis and revised plan (2026-11-01)" at the end; it replaces the
> "Agreed order" for everything after d.

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
    `/n/fs/yubi-iw/runs/bench_nested_v2/{configs,runs,subgraphs,rck_runs,rck_results,diag}`.
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

- copy `/n/fs/yubi-iw/runs/bench_nested_v2/configs/` into amplicon-allele-analysis
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
- Diagnostic outputs of this session: `/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck/` (s403 with
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
  `/n/fs/yubi-iw/runs/bench_nested_v2/subgraphs_paths/` (graphs byte-identical; paths with ≥ 2
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
  --export=ALL,BENCH_NESTED_ROOT --wrap 'export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH; python -c "
import os, bench_nested as b
for sc, seed in [(\"nest_samehap_3to1\", 301), (\"ctrl_dup\", 703), (\"nest_cross_1to1\", 505), (\"nest_cross_1to1\", 502)]:
    b.build_subgraph(b.simulate(sc, seed), os.path.join(b.SUBGRAPHS, b.run_name(sc, seed)))"'
```

Step 2 needs Gurobi: 4 runs, s505 up to 25 min. s502 is the passing control from the same
scenario. LIST file `/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0/sigma0.list`:

```
nest_samehap_3to1_s301_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/nest_samehap_3to1_s301/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/nest_samehap_3to1_s301 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0/nest_samehap_3to1_s301_sigma0
ctrl_dup_s703_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/ctrl_dup_s703/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/ctrl_dup_s703 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0/ctrl_dup_s703_sigma0
nest_cross_1to1_s505_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/nest_cross_1to1_s505/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/nest_cross_1to1_s505 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0/nest_cross_1to1_s505_sigma0
nest_cross_1to1_s502_sigma0 /n/fs/scratch/ym3292/bench_nested_sigma0/subgraphs/nest_cross_1to1_s502/union /n/fs/scratch/ym3292/bench_nested_sigma0/runs/nest_cross_1to1_s502 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0/nest_cross_1to1_s502_sigma0
```

```bash
cd /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0 && mkdir -p logs
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
`/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/alpha.list`, submitted with:

```bash
cd /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow && mkdir -p logs
LIST=$PWD/alpha.list DIAG_POOL_LAMBDA=0.2 DIAG_LAMBDA0=0 sbatch --array=0-5%2 --time=02:00:00 ~/RCK/results/bench_nested/diag/rck_diag.sbatch
```

The lines are:

```
ctrl_dup_s703_an1 /n/fs/yubi-iw/runs/bench_nested_v2/subgraphs/ctrl_dup_s703/union /n/fs/yubi-iw/runs/bench_nested_v2/runs/ctrl_dup_s703 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/ctrl_dup_s703_an1 --alpha-novel 1.0
nest_cross_1to1_s505_an1 /n/fs/yubi-iw/runs/bench_nested_v2/subgraphs/nest_cross_1to1_s505/union /n/fs/yubi-iw/runs/bench_nested_v2/runs/nest_cross_1to1_s505 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/nest_cross_1to1_s505_an1 --alpha-novel 1.0
nest_samehap_3to1_s301_an1 /n/fs/yubi-iw/runs/bench_nested_v2/subgraphs/nest_samehap_3to1_s301/union /n/fs/yubi-iw/runs/bench_nested_v2/runs/nest_samehap_3to1_s301 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/nest_samehap_3to1_s301_an1 --alpha-novel 1.0
nest_cross_1to1_s502_an1 /n/fs/yubi-iw/runs/bench_nested_v2/subgraphs/nest_cross_1to1_s502/union /n/fs/yubi-iw/runs/bench_nested_v2/runs/nest_cross_1to1_s502 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/nest_cross_1to1_s502_an1 --alpha-novel 1.0
nest_cross_3sp_s402_an1 /n/fs/yubi-iw/runs/bench_nested_v2/subgraphs/nest_cross_3sp_s402/union /n/fs/yubi-iw/runs/bench_nested_v2/runs/nest_cross_3sp_s402 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/nest_cross_3sp_s402_an1 --alpha-novel 1.0
nest_samehap_2to1_s201_an1 /n/fs/yubi-iw/runs/bench_nested_v2/subgraphs/nest_samehap_2to1_s201/union /n/fs/yubi-iw/runs/bench_nested_v2/runs/nest_samehap_2to1_s201 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_flow/nest_samehap_2to1_s201_an1 --alpha-novel 1.0
```

### Agreed order (2026-10-31, later) and a correction to the ploidy recommendation

You chose to skip ahead to more seeds, with the re-baseline (old plan steps 3 + 4) done alongside
it rather than skipped. The order is now:

1. **c. Re-baseline the existing 22 cases. [done 2026-11-01: adopted; see "Step c, result" at the end]** Rebuild the 14 bench subgraphs with the
   current `rck_inputs.py`, which has the `clip_path` fix, and the ploidy handling below. Build
   into a scratch root, so the tracked subgraphs and the v2 runs stay as they are until reviewed:
   - root `/n/fs/yubi-iw/runs/bench_nested_v2/rebuild/`, with `configs/` (a copy of the tracked
     ones), `runs` (→ the v2 simulations, unchanged) and `subgraphs/`;
   - RCK outputs in `/n/fs/yubi-iw/runs/bench_nested_v2/rebuild/rck_runs/`.
   Then run RCK on all 22 cases. The 8 regression cases rerun too, so that the new baseline is one
   self-consistent set of runs.
2. **d. More seeds. [done 2026-11-01: d1 Slurm 32049944, d2 Slurm 32050408; adopted; see "Step d2, result" and "Grown bench: failure analysis and revised plan (2026-11-01)" at the end]**
   - Simulate about 4 more accepted seeds each for nest_samehap_3to1, nest_cross_1to1 and
     ctrl_dup, and add a `nest_samehap_6to1` scenario. Simulation needs no Gurobi.
   - Build their subgraphs with the same settings as in c, and run RCK on them.
3. **b. σ = 0 control.** Optional; run it alongside d if licence time allows (commands updated for the re-baseline in "Step b, commands" at the end).
4. **a. Ambiguity reporting.** Deferred until after d.
5. Tune nothing until c and d are both in. [both in since 2026-11-01; the plan after d is "Grown bench: failure analysis and revised plan (2026-11-01)" at the end]

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
  `/n/fs/yubi-iw/runs/bench_nested_v2/rebuild/subgraphs/`, from
  amplicon-allele-analysis `simulations/` with:
  ```bash
  BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2/rebuild BENCH_NESTED_RCK_PLOIDY=2 python bench_nested.py rebuild
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
    `BENCH_NESTED_OUT=/n/fs/yubi-iw/runs/bench_nested_v2/rebuild/rck_runs`;
  - logs in `/n/fs/scratch/$USER/bench_nested_v2/rebuild/logs/`.
- **When it finishes:**
  - check it with `sacct -j 32041913 -X -o JobID%14,State,Elapsed`;
  - then run `collect.py` with the same two variables (writing to a scratch folder, not over this
    folder's v2 tables) and compare the result with v2;
  - only then copy the rebuilt subgraphs into amplicon-allele-analysis
    `simulations/subgraphs/bench_nested/`, and their RCK tables here.

**Step c, result (2026-11-01): done, and adopted as the baseline.**
- Slurm 32041913: 22 of 22 tasks COMPLETED; the longest is s403 (25:41). 2762 s of RCK wall time.
- Against v2 (`rebuild_2026-10-31/comparison_vs_v2.txt`, which now holds the real row-by-row
  comparison; the first commit's file was a copy of `collect.log`):
  - rank 1 at λ = 0.2 is unchanged in 21 of 22 cases. The path fix moves no composite.
  - `nest_cross_3sp_s403`: 0.592 → 0.545. Mean at λ = 0.2: 0.9283 → 0.9262. λ = 0.2 stays the
    recommendation.
  - Small side effects: s402's λ = 0 run (0.5221 → 0.5206), s703's second weight (16.22 → 15.53,
    same walks and composite), and pool sizes ±1–2 in four cases.
- **s403 against the 2026-10-27 ploidy-2 diagnostic** (`/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck/nest_cross_3sp_s403_ploidy2`):
  - `observations.json`, `flow_segments.tsv` and `flow_adjacencies.tsv` are byte-identical in
    both the λ = 0 and the pool run.
  - The 32 rank-1 rows of `decompositions.tsv` are identical except the path columns
    (`reads_explained`, `path_support`, `paths_matched`).
  - Rank 1 is 5.75 / 4.06 / 3.33 (0.519 / 0.593 / 0.523), merged (ec1 carries ecdna1 + ecdna2,
    ec2 carries ecdna1 + ecdna3).
  - The pool's best is 0.729 (#15 of 83), as in the step 2 w1 pool (#16 of 85), and it is not
    picked.
  - The scale is right: depth per copy 13.92 (+1% from the median); segment and junction CN
    1.03× and 1.04× truth (ratio 0.99). The flow objective is 11.2 (41.2 at class 6).
    `pipeline.tsv` stage: `RCK`.
  - **So s403's remaining failure is in the decomposition, not the scale.** It is the type (b)
    failure of "Assessment" section 2: the pool has no good answer (best 0.729).
- **Adopted:**
  - the 14 rebuilt subgraphs were copied into amplicon-allele-analysis
    `simulations/subgraphs/bench_nested/` (now identical to `rebuild/subgraphs/`);
  - `rebuild_2026-10-31/{bench_sweep,pipeline,rerank_validation}.tsv` were copied here
    (`rerank_validation.tsv` is unchanged);
  - `python collect.py` with `BENCH_NESTED_OUT=/n/fs/yubi-iw/runs/bench_nested_v2/rebuild/rck_runs`
    and the tracked subgraphs reproduces all three tables byte for byte;
  - README.md has the re-baseline's results, with v2's kept under "v2 results before the
    re-baseline".
- The local `runs` link still points at v2's `rck_runs`. Pass `BENCH_NESTED_OUT` as above, or
  repoint it.

**Step d, prepared (2026-11-01; d1 ran the same day as Slurm 32049944, d2 submitted as Slurm 32050408; see "Job list" below).**

Code (amplicon-allele-analysis `simulations/`; no default changed):
- `make_scenarios.py`: new scenario `nest_samehap_6to1`, last in `NESTED_SCENARIOS`, so the others
  keep their seed blocks; its block is 801–810.
  - It is nest_samehap_3to1 with `ecdna_proportions` [6, 1] and `num_ecdna_units` 60 (demo7: 30).
  - At 30 units the light species would get 2.5–4.5 copies (bulk totals of the samehap seeds are
    17.7–30.9 at 30 units). At about 6 reads per copy, its junctions would mostly fail the
    ≥ 20-read seed check. At 60 units it gets about 5–9 copies and the heavy one 30–53.
  - Intent: `nested`, `samehap`.
  - The 24 tracked configs regenerate byte-identically.
- `bench_nested.py scenario --candidates N`: the first N seeds of the scenario's block of 100
  (default 10 = `NESTED_SEEDS`, as before).
  - It is needed because the 10-seed blocks are too small for 6 accepted seeds. nest_cross_1to1
    accepted 2 of 5 and has 5 candidates left; nest_samehap_3to1 accepted 2 of 4 and has 6 left.
- `bench_nested.sbatch`:
  - `BENCH_NESTED_ARGS` is passed to `scenario` (unset: nothing);
  - `--array` is 0-7 (the 8th scenario).
- Checked locally:
  - re-checking the 12 existing seeds of the three scenarios gives their recorded `seeds.tsv`
    fields exactly. So `--accept 6` re-derives the accepted seeds in slots 1–2 without
    re-simulating, and adds slots 3–6.

Root `/n/fs/yubi-iw/runs/bench_nested_v2/grow/` (set up; nothing run in it):
- `configs/`: copies of the tracked configs and `seeds.tsv`;
- `runs/`: one link per existing v2 simulation, plus copies of the per-scenario `*.seeds.tsv`.
  New seeds are simulated into this folder, not into v2's;
- `subgraphs/`: copies of the 14 adopted subgraphs (so they are not rebuilt);
- `rck_runs/`: one link per case to `rebuild/rck_runs/<case>`. `run_case.py` sees all 22 as
  done, so only new seeds run;
- `logs/`.

Step d1, simulation and subgraphs (CPU only, no Gurobi). Array tasks 2, 4, 6 and 7 are
nest_samehap_3to1, nest_cross_1to1, ctrl_dup and nest_samehap_6to1.
- Each task wants 6 accepted seeds out of up to 30 candidates.
- About 2 min per simulation and 1–3 min per subgraph. Peak 2.4 GB in the v1 array 31981888.
- It also re-checks the existing seeds, which takes seconds.

```bash
cd ~/amplicon-allele-analysis/simulations
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2/grow BENCH_NESTED_RCK_PLOIDY=2
BENCH_NESTED_ARGS="--accept 6 --candidates 30" sbatch --array=2,4,6,7 --mem=3G --time=01:30:00 \
  --output=$BENCH_NESTED_ROOT/logs/%x_%A_%a.out bench_nested.sbatch
# when done (one sacct, not a polling loop):
sacct -j <JOBID> -X -o JobID%14,State,Elapsed,MaxRSS
python bench_nested.py seeds && python bench_nested.py report     # re-merge after concurrent tasks
awk -F'\t' 'NR>1 && $6=="True" {print $1}' $BENCH_NESTED_ROOT/configs/seeds.tsv | sort | uniq -c
```

Before step d2, check:
- 6 accepted per grown scenario;
- `subgraph_checks.tsv` has no `subgraph_fail`;
- every new `scale.json` says `ploidy class 2 (set with --rck-ploidy)`.

Step d2, RCK (Gurobi), on the new cases only:
- 16–22 cases: 4 each for nest_samehap_3to1, nest_cross_1to1 and ctrl_dup, and up to 6 for 6to1.
- Estimated cost: nest_cross_1to1 5–25 min each, the others 1–3 min. About 1–2 licence-hours,
  about 1 h of wall time at 2 concurrent tasks.

```bash
cd ~/RCK/results/bench_nested
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2/grow BENCH_NESTED_OUT=/n/fs/yubi-iw/runs/bench_nested_v2/grow/rck_runs
IDX=$(python run_case.py list | while IFS=$'\t' read i name kind sub; do
        [ -e "$BENCH_NESTED_OUT/$name/pool/pool/pool.json" ] || echo $i; done | paste -sd,)
echo "$IDX"                                   # only the new cases
sbatch --array="$IDX"%2 --output=$BENCH_NESTED_ROOT/logs/%x_%A_%a.out bench_nested.sbatch
# when done:
mkdir -p $BENCH_NESTED_ROOT/rck_results
BENCH_NESTED_RESULTS=$BENCH_NESTED_ROOT/rck_results python collect.py
```

- Compare `rck_results/` with this folder's tables: the 22 existing rows must be identical.
- Then decide whether to adopt, as for c: configs, `seeds.tsv` and the new subgraphs into
  amplicon-allele-analysis; the tables here.
- `collect.py`'s per-scenario table and λ rule then see 6 seeds per grown scenario. Its
  "7 scenarios × 2 seeds" wording in README.md will need updating.

**Step b, commands (updated 2026-11-01 for the re-baseline; optional, not submitted).**

"Run A" above, with two changes:
- the subgraphs are built with `BENCH_NESTED_RCK_PLOIDY=2`, as in c and d;
- the comparison baseline is `rebuild/rck_runs/<case>`, not v2's `rck_runs`.

The `clip_path` fix is already in `rck_inputs.py`.

```bash
cd ~/amplicon-allele-analysis/simulations
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export BENCH_NESTED_ROOT=/n/fs/scratch/$USER/bench_nested_sigma0 BENCH_NESTED_RCK_PLOIDY=2
mkdir -p $BENCH_NESTED_ROOT/configs $BENCH_NESTED_ROOT/logs
for c in nest_samehap_3to1_s301 ctrl_dup_s703 nest_cross_1to1_s505 nest_cross_1to1_s502; do
  sed 's/^  sigma: 0.1$/  sigma: 0.0/' configs/bench_nested/$c.yaml > $BENCH_NESTED_ROOT/configs/$c.yaml
done
grep -c '^  sigma: 0.0$' $BENCH_NESTED_ROOT/configs/*.yaml          # 1 each
sbatch --job-name=sigma0-sim --cpus-per-task=1 --mem=3G --time=01:00:00 \
  --output=$BENCH_NESTED_ROOT/logs/%x_%j.out --export=ALL --wrap 'python -c "
import os, bench_nested as b
for sc, seed in [(\"nest_samehap_3to1\", 301), (\"ctrl_dup\", 703), (\"nest_cross_1to1\", 505), (\"nest_cross_1to1\", 502)]:
    b.build_subgraph(b.simulate(sc, seed), os.path.join(b.SUBGRAPHS, b.run_name(sc, seed)))"'
# truth must match v2 (sigma only changes the counts):
for c in nest_samehap_3to1_s301 ctrl_dup_s703 nest_cross_1to1_s505 nest_cross_1to1_s502; do
  diff -rq $BENCH_NESTED_ROOT/runs/$c/truth /n/fs/yubi-iw/runs/bench_nested_v2/runs/$c/truth && echo "$c truth ok"; done
```

The RCK step is the LIST and `rck_diag.sbatch` command of "Run A" unchanged: 4 runs, s505 up to
25 min, about 45 min of licence time. **The licence allows 2 concurrent sessions**, so do not run it
next to the step d2 array. Chain it with `--dependency=afterany:<d2 JOBID>`, or submit it before d2.

**Job list (2026-11-01): what is submitted and what is left.** All commands are in "Step d,
prepared" and "Step b, commands" above.

| step | job | Gurobi | status | needs |
|---|---|---|---|---|
| d1: simulate and build subgraphs, 6 accepted seeds each for nest_samehap_3to1, nest_cross_1to1, ctrl_dup and nest_samehap_6to1 (array tasks 2, 4, 6, 7) | Slurm **32049944** | no | done 2026-11-01: 4 of 4 COMPLETED (5:53–9:16); all checks pass ("Checking d1, result" below) | — |
| d2: RCK on the 18 new cases only (`--array=14,15,16,17,22,23,24,25,30,31,32,33,34,35,36,37,38,39%2`) | Slurm **32050408** | yes, about 1–2 licence-hours | done 2026-11-01: 18 of 18 COMPLETED ("Step d2, result" below) | d1 |
| d2 collect: `collect.py` into `grow/rck_results/`, then compare the 22 existing rows with this folder's tables | — (local, CPU, about 1–2 min) | no | done 2026-11-01; adopted ("Step d2, result") | d2 |
| b: σ = 0 simulation and subgraphs (`sigma0-sim`) | not submitted | no | optional; can run any time | — |
| b: σ = 0 RCK (`rck_diag.sbatch`, 4 runs) | not submitted | yes, about 45 min | optional; d2 is done, so both licence sessions are free | b sim |
| a: ambiguity reporting; any tuning | — | — | d is in: see "Grown bench: failure analysis and revised plan (2026-11-01)" | — |

**Checking d1 (Slurm 32049944).** Estimated 30–45 min from the start; worst case about 75 min,
within the 1:30 limit. Logs: `/n/fs/scratch/$USER/bench_nested_v2/grow/logs/bench-nested-sim_32049944_<task>.out`.

```bash
sacct -j 32049944 -X -o JobID%14,JobName%18,State,Elapsed,MaxRSS     # COMPLETED x4? (do not loop it)
tail -n 3 /n/fs/scratch/$USER/bench_nested_v2/grow/logs/bench-nested-sim_32049944_*.out
# a task exits 1 when a scenario got fewer than 6 accepted seeds; its log says which slot failed
cd ~/amplicon-allele-analysis/simulations && export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2/grow BENCH_NESTED_RCK_PLOIDY=2
python bench_nested.py seeds && python bench_nested.py report
awk -F'\t' 'NR>1 && $6=="True" {print $1}' $BENCH_NESTED_ROOT/configs/seeds.tsv | sort | uniq -c   # 6 for the 4 scenarios, 2 for the others
grep -c "set with --rck-ploidy" $BENCH_NESTED_ROOT/subgraphs/*/union/scale.json | grep -v ':1$'  # empty
```

Also check that `subgraph_checks.tsv` has no `subgraph_fail` True. If a task timed out or failed
for a transient reason, resubmit the same command with `--array=<that task>`: finished steps are
skipped.

**Checking d1, result (2026-11-01): all pass.**
- `sacct -j 32049944`: tasks 2, 4, 6 and 7 COMPLETED, in 6:17, 5:53, 9:16 and 9:01.
- `bench_nested.py seeds` and `report` re-merge (exit 0); 32 runs, subgraph-stage failures: none.
- Accepted seeds in `grow/configs/seeds.tsv`: 6 each for nest_samehap_3to1, nest_cross_1to1,
  ctrl_dup and nest_samehap_6to1; 2 for the other four scenarios.
- `grow/subgraphs/subgraph_checks.tsv`: `subgraph_fail` is False in all 32 rows; the 18 new runs
  have `verify_verdict` PASS and depth per copy 12.7–15.1.
- All 32 `scale.json` say `set with --rck-ploidy`.
- New accepted seeds: nest_samehap_3to1 s305–s308, nest_cross_1to1 s506–s509, ctrl_dup s706,
  s708, s709, s710, nest_samehap_6to1 s801–s806 (18 runs).

**Step d2, submitted (2026-11-01): Slurm 32050408.**
- Submitted with the "Step d2" block unchanged (`--array="$IDX"%2`, the sbatch's defaults of 4
  CPUs, 2 GB and 1:30), from `~/RCK/results/bench_nested`, with `BENCH_NESTED_ROOT` and
  `BENCH_NESTED_OUT` = `grow/` and `grow/rck_runs`. Logs:
  `/n/fs/scratch/$USER/bench_nested_v2/grow/logs/bench-nested_32050408_<index>.out`.
- `IDX` = 14,15,16,17,22,23,24,25,30,31,32,33,34,35,36,37,38,39. **The indices of the existing
  cases moved** (`run_case.py list` now has 40 cases), so the index numbers here are not those of
  the c array 32041913. Case behind each index:

  | index | case | index | case |
  |---|---|---|---|
  | 14 | nest_samehap_3to1_s305 | 30 | ctrl_dup_s706 |
  | 15 | nest_samehap_3to1_s306 | 31 | ctrl_dup_s708 |
  | 16 | nest_samehap_3to1_s307 | 32 | ctrl_dup_s709 |
  | 17 | nest_samehap_3to1_s308 | 33 | ctrl_dup_s710 |
  | 22 | nest_cross_1to1_s506 | 34 | nest_samehap_6to1_s801 |
  | 23 | nest_cross_1to1_s507 | 35 | nest_samehap_6to1_s802 |
  | 24 | nest_cross_1to1_s508 | 36 | nest_samehap_6to1_s803 |
  | 25 | nest_cross_1to1_s509 | 37 | nest_samehap_6to1_s804 |
  |  |  | 38 | nest_samehap_6to1_s805 |
  |  |  | 39 | nest_samehap_6to1_s806 |

- Checked before submitting: the 22 indices left out (0–13, 18–21, 26–29) are exactly the 8
  regression cases and the 14 existing bench cases. Each is a link to `rebuild/rck_runs/<case>`,
  and `run_case.done()` is True for both its `lambda_0` and `pool` runs. So none of them reruns,
  even if submitted by mistake. All 18 new cases have `edges.tsv`. No other job of ours was
  queued, so the array has both licence sessions.
- The task-to-case mapping is recomputed at run time from `grow/configs/seeds.tsv`. Do not
  re-run `bench_nested.py scenario` in `grow/` while 32050408 is queued.

When it finishes:

```bash
sacct -j 32050408 -X -o JobID%14,JobName%14,State,Elapsed,MaxRSS     # 18 x COMPLETED? (once, no loop)
grep -L "pool in" /n/fs/scratch/$USER/bench_nested_v2/grow/logs/bench-nested_32050408_*.out  # any task without a finished pool run
cd ~/RCK/results/bench_nested && export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export BENCH_NESTED_ROOT=/n/fs/yubi-iw/runs/bench_nested_v2/grow BENCH_NESTED_OUT=/n/fs/yubi-iw/runs/bench_nested_v2/grow/rck_runs
mkdir -p $BENCH_NESTED_ROOT/rck_results
BENCH_NESTED_RESULTS=$BENCH_NESTED_ROOT/rck_results python collect.py 2>&1 | tee $BENCH_NESTED_ROOT/rck_results/collect.log
# the 22 existing rows must be identical to this folder's tables:
for t in bench_sweep pipeline rerank_validation; do
  python - "$t" <<'PY'
import sys, pandas as pd
t = sys.argv[1]; key = "run" if t == "pipeline" else None
old = pd.read_csv(f"{t}.tsv", sep="\t", dtype=str)
new = pd.read_csv(f"/n/fs/yubi-iw/runs/bench_nested_v2/grow/rck_results/{t}.tsv", sep="\t", dtype=str)
merged = old.merge(new, how="left", indicator=True)
print(t, "old rows", len(old), "new rows", len(new), "old rows missing or changed in new:", (merged["_merge"] != "both").sum())
PY
done
```

**Step d2, result (2026-11-01): done, and adopted.**
- `sacct -j 32050408`: 18 of 18 tasks COMPLETED (elapsed 6:56–27:07, including licence waits).
  RCK wall time over all 40 cases is 7030 s (2762 s for the 22 existing ones).
- `collect.py` into `grow/rck_results/` (copied to `grow_2026-11-01/` here): mean composite at
  λ = 0.2 is 0.8531 over 40 cases (0.9262 over 22). λ = 0.2 stays the recommendation.
- **The 22 existing cases are unchanged:**
  - `bench_sweep.tsv` (132 rows) and `rerank_validation.tsv` (48 rows) are identical row for row.
  - `compare_existing.txt` checked `bench_sweep.tsv` only. In `pipeline.tsv` all 14 existing rows
    differ, in `depth_vs_median` alone (by 0.000–0.005). That column is relative to the bench
    cases' median depth per copy, which moves from 13.762 to 13.7065 with the 18 new cases.
    `depth_per_copy` and `scale_fail` are unchanged in every row, and recomputing the column
    from the new median reproduces it exactly.
- **Adopted** as for c:
  - into amplicon-allele-analysis `simulations/`: the 21 new configs (including the rejected
    s704, s705 and s707, as the existing scenarios keep theirs), `seeds.tsv` (rows added, none
    changed), the 18 new subgraph folders and `subgraph_checks.tsv`, now identical to
    `grow/configs` and `grow/subgraphs` (23 MB in all);
  - the three tables here, from `grow_2026-11-01/`;
  - the untracked links: `runs` → `grow/rck_runs`, and amplicon-allele-analysis
    `simulations/runs/bench_nested` → `grow/runs`;
  - `python collect.py` with no `BENCH_NESTED_*` variables then reproduces the three tables and
    the README results block byte for byte.
- **Two traps for the next full rebuild or rerun** (README "How to run" now says so; no file
  changed):
  - RCK `bench_nested.sbatch` still has `#SBATCH --array=0-21%2`, so a plain `sbatch` misses
    the 18 new cases: pass `--array=0-39%2`.
  - amplicon-allele-analysis `bench_nested.sbatch` defaults to `--accept 2`. Rerun on the grown
    scenarios, it rewrites their seed records to 2 slots and re-merges `seeds.tsv` without
    s305–s308, s506–s509, s706–s710 and s801–s806. Pass `BENCH_NESTED_ARGS="--accept 6 --candidates 30"`
    for tasks 2, 4, 6 and 7.

## Grown bench: failure analysis and revised plan (2026-11-01)

Steps c and d are in, so the "tune nothing until c and d are in" hold is lifted. This section reads
the 17 failures of the 40-case bench from the saved pools alone: no new runs and no Gurobi.

**Scripts** (offline, in `diag/`):
- `grow_failures.py` gives, per case at λ = 0.2:
  - rank 1, with its integer weights and low_weight flags and its refit drops;
  - its composites with and without low_weight species, and which rank-1 species overlap
    which truth species;
  - the pool's best solution: its structure cost against rank 1's, whether the λ pruning keeps
    it, and whether it is a recorded same-fit alternative.
- `grow_rules.py` re-ranks all 40 pools with RCK's own `rank_pool` under candidate rules. Its
  `default` column reproduces `bench_sweep.tsv` at λ = 0.2 (0.8531) and the 17 failures.

**Outputs** are in `grow_2026-11-01/`: `grow_failures.tsv`, `grow_failures_species.tsv` and
`grow_rules_lambda{0.2,0}.{tsv,txt}`. λ = 0 there is emulated from the pool (0.7714 against 0.7713
run).

A correction to the counts in the brief: nest_cross_1to1 fails in **4** of 6 seeds (s502, s505,
s507, s508), not 3. nest_samehap_6to1 fails in all 6: in 5 the nested species scores 0, and s804
scores 0.841. ctrl_dup fails in 3 of 6 (s703, s708, s710), and nest_samehap_3to1 in 3 of 6 (s301,
s307, s308).

### 1. nest_samehap_6to1: the nested species is split into low_weight pieces

It is not absent from the flow, not dropped by the refit, and not absorbed into the outer species.

| seed | truth | rank 1 (integer; L = low_weight) | threshold max(1, 0.1 × heaviest) | ecdna2 with low_weight included | pool solutions with ecdna2 ≥ 0.9 | pool best (rank at λ 0.2) | refit misfit, best − rank 1 |
|---|---|---|---|---|---|---|---|
| s801 | 36.4 / 6.1 | 33, 3L, 2L, 2L | 3.3 | 0.675 (`ec2`, 3 copies, Jaccard 0.97 to ecdna2) | 4 of 87 | 0.972: 32 / 5 (#32) | +0.007 |
| s802 | 36.9 / 6.1 | 36, 2L, 1L | 3.6 | 0.520 | 10 of 61 | 0.971: 35 / 4 (pruned at λ 0.2) | +0.005 |
| s803 | 37.6 / 6.3 | 35, 3L, 2L, 2L | 3.5 | 0.591 | **0** of 70 | 0.755: 33 / 4 (#62) | +0.007 |
| s804 | 40.8 / 6.8 | 38, 6, 2L, 1L | 3.8 | 0.695 (nested found, wrong walk) | **0** of 19 | 0.842 = rank 1 | +0.015 |
| s805 | 41.9 / 7.0 | 40, 3L, 2L, 2L | 4.0 | 0.584 (`ec2` Jaccard 0.78) | 4 of 55 | 0.955: 40 / 4 (#35) | +0.013 |
| s806 | 43.0 / 7.2 | 42, 3L, 2L, 1L | 4.2 | 0.578 | 16 of 72 | 0.978: 42 / 6 (#56) | −0.166 |

- **Flagged low_weight.** In the 5 seeds that score 0, rank 1 has one significant species (33–42
  copies, composite 0.98–0.996). The nested species' copies are spread over 2–3 walks of 1–3
  copies, each below max(1, 0.1 × heaviest) = 3.3–4.2. In s801 one of them *is* the nested walk
  (3 copies, Jaccard 0.97); the others are chimeras of outer and nested segments. The evaluation
  leaves low_weight species out, so ecdna2 gets 0. Scored with them, it gets 0.52–0.675.
- **Not dropped by the refit.** No 6to1 rank 1 has a refit drop; the only drops are 1-copy species
  at ranks 3–4 of s803.
- **Not absorbed.** The outer species is at or below its true weight (33–42 against 36–43).
- **The data cannot decide.** The 2-species answers fit the observations as well as rank 1:
  refit misfit within 0.015 in 5 seeds, and s806's is better by 0.17 of 12.
- **Why rank 1 wins.**
  - The ranking puts residual first (all 0), then species_penalty × significant species
    + λ × steps.
  - Splitting the nested species below the threshold makes it free: 1 significant species
    against 2.
  - At λ = 0.2 the right answer also pays for its 7–16 extra significant steps (1.4–3.2), more
    than a species.
  - At λ = 0 the leftover rule alone already decides: the 6to1 mean at λ = 0 is the same 0.553.
- **Two seeds also fail before the ranking.** s803 and s804 have no ≥ 0.9 answer for the nested
  species in their pools. s804's flow has chr3:41.12–41.52 Mb at 62 copies against a truth of
  47.65 (depth reads 64), and it uses an imputed edge (`imputed_adj_intra_chr3_10`). This is
  the s505 type (b) failure.
- **Same mechanism in nest_samehap_3to1_s307 (new).** Rank 1 is 15, 1L, 1L, with the nested species
  (5.3) gone. Its 0.977 answer (16 / 4 / 1) is rank 1's first same-fit alternative.

So 6to1 is s703/s301's identifiability problem, made systematic. It is the risk "Assessment"
section 4 predicted for a leftover-type prior, and the *default* rule (10% of the heaviest
species) already has it.

### 2. ctrl_dup and nest_cross_1to1 (and the new 3to1 failures): old causes or new?

| case | new? | composite | cause | same as s703 / s505 / s502? |
|---|---|---|---|---|
| ctrl_dup_s703 | old | 0.712 | indistinguishable exact fits; 2-copy leftover counted (unchanged) | — |
| **ctrl_dup_s708** | new | 0.732 | **multiplicity cap**: the truth traverses a segment 4 times, and `--max-multiplicity` is 3. No pool solution has 1 species (pool: 15 with 2, 54 with 3) | **no: new cause** |
| **ctrl_dup_s710** | new | 0.723 | **multiplicity cap**, as s708 (pool: none below 3 species). Also an imputed edge in the flow and depth +40% on chr2:103.76–104.03 Mb (163.8 against 117.0) | **no: new cause** |
| nest_cross_1to1_s502 | old | 0.955 | θ (junction scale 0.77); fails only on the label (see the refit re-flag below) | — |
| nest_cross_1to1_s505 | old | 0.765 | biased flow; no good answer in the pool (unchanged) | — |
| **nest_cross_1to1_s507** | new | **0.961** | composites 0.969 / 0.952, above the 0.95 line. Fails only on the label: two 1-copy leftovers, low_weight in the ranking, refit to 1.23 and 1.09, which is ≥ max(1, 0.1 × 10.77), so the output counts them, "4 species for 2", and `ec3` carries both junction sets | **leftover**, but in the report, not the ranking |
| **nest_cross_1to1_s508** | new | 0.672 | indistinguishable exact fits: the 0.988 answer has the same significant weights (13 / 12) and steps (41) as rank 1. It is a same-fit alternative (within the first 5), but outside `pool_gap` at λ 0.2; the refit misfit differs by 0.105 of 5.5. θ is right (junction scale 0.98) | **s703 type** |
| nest_samehap_3to1_s301 | old | 0.693 | indistinguishable; 0.973 answer is the 2nd same-fit alternative (unchanged) | — |
| **nest_samehap_3to1_s307** | new | 0.493 | the 6to1 mechanism (section 1) | s301 / 6to1 type |
| **nest_samehap_3to1_s308** | new | 0.949 | nested weight 3.6 against 5.8; junction scale 0.75, like s502's 0.77, so θ may contribute; the 0.980 answer is ranked 2nd | marginal; possibly θ |

- **θ** explains none of the new failures outright. The junction scale is 0.98–1.06 in s507, s508
  and every 6to1 seed. s708 (0.77) and s710 (0.85) are decided by the multiplicity cap.
- **Multiplicity.** s708 and s710 are the only accepted cases with `segment_max_traversals` 4
  (`seeds.tsv`); every other case has ≤ 3. The cap is per edge and per species (`decompose.py`,
  `mult <= max_multiplicity`), so the true single species cannot be expressed.
  - The bench's seed check does not compare traversals with RCK's cap.
  - The 22-case bench never hit it: s703 has 3.
- **The refit re-flag.** The ranking flags low_weight on integer weights. `refit_decompositions`
  re-flags on refitted weights, and the output, `collect.py` and the failure rule read the
  latter. Rank 1's significant-species count differs between the two in 8 of 40 cases, both ways:
  - the output counts *more* in s502, s507, s508 and s301: a 1-copy leftover refits to 1.2–1.5;
  - it counts *fewer* in s603, s506 and s509: 1-copy species refit to < 1 or are dropped. Those
    three pass because of this.
  - Aligning the report with the ranking would swap 2 failures for 3.
  - Ranking on the refitted flags (`grow_rules.py` `refit`) gives 0.828, with 2 cases up and 5
    down (s502 0.955 → 0.572).
  - So neither direction is a fix. Read s502 and s507 as composite passes with a reported 1-copy
    leftover.

### 3. Candidate rules, offline (40 pools; λ = 0.2, λ = 0 in brackets)

| rule | mean composite | up / down (> 0.005) | breaks |
|---|---|---|---|
| default | **0.8531** (0.7714) | — | — |
| leftover mass (low_weight copies that add up to the threshold count as species) | 0.8442 (0.8081) | 0 / 2 (4 / 1) | demo3, demo7_totalonly. At λ 0 it fixes s801, s806 and s307, and λ 0.2's step cost undoes it |
| rank on refitted flags | 0.8281 (0.7125) | 2 / 5 (0 / 7) | s502, s403, s301, s305, s804 |
| global 0.05 | 0.8546 | 5 / 4 | s706, demo3, s507, s305 |
| global 0.075 | 0.8308 | 1 / 3 | s706, demo3, s507 |
| global 0.15 | 0.8914 | 7 / 5 | s502, s505, s403, s308, s804 |
| **set metric, rank 1 + first k same-fit alternatives** (not a ranking) | k = 1: 0.9113; **k = 2: 0.9229**; k = 5: 0.9373 | 6 / 0; 8 / 0; 10 / 0 | none (rank 1 is unchanged) |
| pool best (ceiling) | 0.9470 | | |

Every ranking rule trades cases, as step 2b did. With 6 seeds per scenario the trade is now seen
in more than one seed: each rule that fixes 6to1 seeds (global 0.15: s802, s805, s806) breaks
s502, s505 or s403. The set metric is the only change with no loss.

### Revised plan (2026-11-01), in order

Δ is the change in the 40-case mean composite at λ = 0.2. One case moving by d moves it by d/40.
Licence time is Gurobi time at the WLS limit of 2 sessions.

| # | step | cases it targets | expected Δ | risk to passing cases | licence time |
|---|---|---|---|---|---|
| 1 | **a. Ambiguity reporting** (RCK output, opt-in; `collect.py` set metric) | s703, s301, s307, s802, s805, s806 (≥ 0.955 at k = 2); s801 (0.813) and s403 (0.606) partly; k = 5 adds s508 (0.987) and s803 (0.755) | rank-1 mean **0** (unchanged); set metric **+0.070** (k = 2, 0.923) or **+0.084** (k = 5, 0.937) | **none**: rank 1 unchanged, 0 cases down offline | **0**: offline from `pool.json`; code + tests only |
| 2 | **Multiplicity check** for ctrl_dup s708 and s710: `rck_diag.sbatch` at `--max-multiplicity 4` (opt-in option, no default change) | s708, s710 | if both reach about 0.97: **+0.012** rank-1 and set metric | none as a 2-case diagnostic. Adopting it bench-wide changes every pool (3 bits per edge and species instead of 2) and needs the 40-case rerun (about 2 licence-hours). A bench seed rule "traversals ≤ `--max-multiplicity`" instead needs no licence, but changes the case set | **< 10 min** (their runs take 2 s each now; allow 10× for the larger MILP) |
| 3 | **b. σ = 0 control**, extended | as planned: s301, s703, s505, s502; added: **s801** (does the 6to1 split survive noise-free data? If so the prior alone causes it), **s803** and **s804** (no good answer in the pool: noise or model?), **s508** | **0** (diagnostic). It decides whether A2 (depth-bias-aware flow) or a nested-species prior is worth building | none | about **45–50 min** (s505 is 13 min of the 40-case total; the 4 added cases about 5 min) |
| 4 | **6to1 decision (yours; no compute).** After 1–3, decide whether RCK should prefer "outer + a light nested species" over "outer + leftovers" when the data tie. Lead to test then: the leftover-mass rule at λ 0 (fixes s801, s806, s307; costs demo3), which means revisiting λ = 0.2. λ = 0.2 was chosen on a bench without a light nested species, and its step cost is what blocks the fix | s801, s802, s805, s806, s307 | up to +0.06 if all 5 were fixed, but no rule tested gets there without losses (section 3) | high: every rule tested breaks 2–5 cases | offline first (0); a rerun only after a rule passes offline |
| 5 | Bookkeeping (no licence; independent) | — | 0 | none | 0 |

Step 5's items:
- extend `compare_existing.txt`'s check to all three tables;
- decide whether RCK `bench_nested.sbatch`'s `--array` becomes 0-39 (your call; README says
  pass it);
- note in `collect.py`'s failure rule that the species count includes refitted 1-copy leftovers
  (s502 and s507 fail only on that). Do not change the rule to hide them: the same re-flag
  makes s506, s509 and s603 pass.

**Order rationale.**
- Step 1 is the largest gain, with no risk and no licence time, and the evidence for it got
  stronger: 6 of the new failures are same-fit ties the data cannot decide.
- Step 2 is the one rank-1 fix with a clear mechanism and almost no cost.
- Step 3 rises in value because it can now separate prior from noise in 6to1, too.
- Step 4 waits for 1–3, because every ranking rule tested trades cases.

**Deferred or dropped.**
- Dropped:
  - global threshold changes (all trade cases);
  - ranking on refitted flags;
  - leftover mass at λ = 0.2;
  - step 2b (as before).
- Deferred:
  - θ fit, for s502 and possibly s308;
  - A2, until step 3 reports.

### How steps 1–3 would be carried out (2026-11-01; nothing submitted)

**Step 1, ambiguity reporting.** Code and tests only; no Gurobi needed to judge it.
1. RCK `rck/ecdna/decompose.py`:
   - `same_fit_structures` keeps only each alternative's weights. Make it keep the walks too.
   - Add an opt-in `rck-ecdna --report-alternatives K` (default 0, so outputs are unchanged).
     When rank 1's choice is `decided_by` a prior (`species_penalty` or `step_penalty`), it writes
     the first K same-fit alternatives to `alternatives.tsv`, in the `decompositions.tsv` layout:
     walks, integer and refitted weights, status.
   - Document it in `docs/ecDNA.md` and add a test in `tests/ecdna/` (a small graph with two
     same-fit decompositions).
2. `collect.py`: add `set_composite` (best of rank 1 and its first k = 2 alternatives, as
   `diag/grow_rules.py` `set 2` computes it from `pool.json`) to `bench_sweep.tsv`, `pipeline.tsv`
   and the README scenario table, next to the rank-1 composite. The rank-1 columns, the λ rule and
   the failure rule stay as they are.
3. Check: `python collect.py` must leave every existing column byte-identical, and give a
   set-metric mean of 0.9229 at λ = 0.2 (0.8944 at λ = 0, emulated).
   - `--report-alternatives` itself is first exercised by the step 2 runs below (add
     `--report-alternatives 2` to their lines once it exists).
   - Its walks must match the pool entries `grow_rules.py` matched.

**Step 2, multiplicity check** (Gurobi, < 10 min). ctrl_dup s708 and s710, plus two guards:
- s703, which has 3 traversals and fails;
- s709, which has 3 traversals and passes.

The option changes only the decomposition MILP's per-edge cap. LIST file
`/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_mult/mult4.list`:

```
ctrl_dup_s708_m4 /n/fs/yubi-iw/runs/bench_nested_v2/grow/subgraphs/ctrl_dup_s708/union /n/fs/yubi-iw/runs/bench_nested_v2/grow/runs/ctrl_dup_s708 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_mult/ctrl_dup_s708_m4 --max-multiplicity 4
ctrl_dup_s710_m4 /n/fs/yubi-iw/runs/bench_nested_v2/grow/subgraphs/ctrl_dup_s710/union /n/fs/yubi-iw/runs/bench_nested_v2/grow/runs/ctrl_dup_s710 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_mult/ctrl_dup_s710_m4 --max-multiplicity 4
ctrl_dup_s703_m4 /n/fs/yubi-iw/runs/bench_nested_v2/grow/subgraphs/ctrl_dup_s703/union /n/fs/yubi-iw/runs/bench_nested_v2/grow/runs/ctrl_dup_s703 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_mult/ctrl_dup_s703_m4 --max-multiplicity 4
ctrl_dup_s709_m4 /n/fs/yubi-iw/runs/bench_nested_v2/grow/subgraphs/ctrl_dup_s709/union /n/fs/yubi-iw/runs/bench_nested_v2/grow/runs/ctrl_dup_s709 /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_mult/ctrl_dup_s709_m4 --max-multiplicity 4
```

```bash
cd /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_mult && mkdir -p logs   # write mult4.list here first
LIST=$PWD/mult4.list DIAG_POOL_LAMBDA=0.2 sbatch --array=0-3%2 --time=00:30:00 ~/RCK/results/bench_nested/diag/rck_diag.sbatch
# when done (one sacct, no loop):
python ~/RCK/results/bench_nested/diag/noise_summary.py mult4.list
python ~/RCK/results/bench_nested/diag/grow_failures.py . ctrl_dup_s708_m4 ctrl_dup_s710_m4 ctrl_dup_s703_m4 ctrl_dup_s709_m4
```

Reading the result:
- **s708 and s710 reach about 0.97 with 1 species.** The cap is the whole cause. Then choose
  between two fixes:
  - (i) a seed rule in amplicon-allele-analysis `check_seed`: reject `segment_max_traversals`
    above 3. No licence time, but it changes the case set: s708 and s710 are replaced by new
    seeds, and the mean is not comparable;
  - (ii) `--max-multiplicity 4` for every bench case. Every pool changes (3 bits per edge and
    species instead of 2), so it needs the 40-case rerun, about 2 licence-hours.
  - RCK's default stays 3 either way.
- **They gain a 1-species answer that still loses.** It is then a ranking question, like s703.
- **s703 or s709 change.** The larger MILP alone moves answers, and (ii) needs the full rerun
  before it can be judged.

**Step 3, σ = 0 control, extended** (simulation on CPU, then about 45–50 min of Gurobi). It is
"Step b, commands" with 4 more cases:
- nest_samehap_6to1 s801 (the prior on noise-free data);
- s803 and s804 (no good answer in the pool);
- nest_cross_1to1 s508 (a same-fit tie).

Their configs are now tracked, so the `sed` loop takes all 8:

```bash
cd ~/amplicon-allele-analysis/simulations
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
export BENCH_NESTED_ROOT=/n/fs/scratch/$USER/bench_nested_sigma0 BENCH_NESTED_RCK_PLOIDY=2
mkdir -p $BENCH_NESTED_ROOT/configs $BENCH_NESTED_ROOT/logs
CASES="nest_samehap_3to1_s301 ctrl_dup_s703 nest_cross_1to1_s505 nest_cross_1to1_s502 nest_samehap_6to1_s801 nest_samehap_6to1_s803 nest_samehap_6to1_s804 nest_cross_1to1_s508"
for c in $CASES; do sed 's/^  sigma: 0.1$/  sigma: 0.0/' configs/bench_nested/$c.yaml > $BENCH_NESTED_ROOT/configs/$c.yaml; done
grep -c '^  sigma: 0.0$' $BENCH_NESTED_ROOT/configs/*.yaml          # 1 each, 8 files
sbatch --job-name=sigma0-sim --cpus-per-task=1 --mem=3G --time=01:30:00 \
  --output=$BENCH_NESTED_ROOT/logs/%x_%j.out --export=ALL --wrap 'python -c "
import os, bench_nested as b
for sc, seed in [(\"nest_samehap_3to1\", 301), (\"ctrl_dup\", 703), (\"nest_cross_1to1\", 505), (\"nest_cross_1to1\", 502),
                 (\"nest_samehap_6to1\", 801), (\"nest_samehap_6to1\", 803), (\"nest_samehap_6to1\", 804), (\"nest_cross_1to1\", 508)]:
    b.build_subgraph(b.simulate(sc, seed), os.path.join(b.SUBGRAPHS, b.run_name(sc, seed)))"'
# truth must match the bench's (sigma only changes the counts):
for c in $CASES; do diff -rq $BENCH_NESTED_ROOT/runs/$c/truth /n/fs/yubi-iw/runs/bench_nested_v2/grow/runs/$c/truth && echo "$c truth ok"; done
```

The RCK step uses the "Run A" LIST, extended to 8 lines of the same form
(`<case>_sigma0 $ROOT/subgraphs/<case>/union $ROOT/runs/<case> .../diag/rck_sigma0/<case>_sigma0`):

```bash
cd /n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_sigma0 && mkdir -p logs
LIST=$PWD/sigma0.list DIAG_POOL_LAMBDA=0.2 DIAG_LAMBDA0=0 sbatch --array=0-7%2 --time=02:00:00 ~/RCK/results/bench_nested/diag/rck_diag.sbatch
```

- Do not run it next to the step 2 array (licence: 2 sessions). Chain it with
  `--dependency=afterany:<step 2 JOBID>`.
- Reading:
  - s801 still splitting its nested species at σ = 0 means the prior alone causes it, which
    feeds step 4;
  - s803 and s804 gaining a good pool answer means their failure is the simulated depth bias;
  - s301, s703, s505 and s502 are read as in "Run A".

## Offline test: junction-incidence and long-read-path priors (2026-11-12): negative

Two more ranking priors were tested offline on the 40 saved pools with
`diag/junction_rules.py` (no Gurobi; gurobipy blocked at import). Outputs are in
`diag/junction_rules_2026-11-12/`. The `default` rule reproduces `bench_sweep.tsv` (0.8531 at
λ = 0.2) and the 17 failures. Failures are scored as in `collect.py`: composite < 0.95, merged, or a
wrong reported species count.

- **Junction incidence** (`jinc μ`): μ × Σ over novel junctions of (species using it − 1), all
  species counted, low_weight included, added to the structure cost. At λ = 0.2 the best is
  16 failures (μ ≥ 0.5: s307 fixed). It costs demo3 (a control, down > 0.005) and s308. μ = 0.1 is
  harmless but fixes nothing. Without the λ pruning (`np`) the results are the same.
  - **Why it fails.** The good answers also carry 1–2 copy leftovers, and those leftovers repeat
    junctions too. For example, s801's 0.972 answer `[32, 5, 2, 1]` scores 13 against rank 1's 12.
    Every pool also has bad answers with fewer repeats (minimum 2–8).
- **Long-read paths** (`paths first`, `paths sum ν`): no change in any case. `reads_explained`
  is 1.0 for almost every pool solution, and 15 of the 40 pools have no multi-junction reads at all
  (among them 6to1 s802 and s803). This bench has no path signal to rank on, so the test says
  nothing about paths on data that has them.

Reading: this is the same conclusion as section 3 above. The right answers need noise leftovers
too, so no prior on leftovers separates them. Next are the two Gurobi diagnostics (plan steps 2
and 3: `--max-multiplicity 4`, and σ = 0), which address the cause rather than the ranking.

## Plan steps 2 and 3 run: multiplicity 4 and σ = 0 (2026-11-12)

Both ran with `diag/rck_diag.sbatch` (pool run at λ = 0.2, never more than 2 licence sessions).
The LIST files and summaries (`diag/noise_summary.py`) are in `diag/rck_mult/` and `diag/rck_sigma0/`.
The run outputs are in `/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_{mult,sigma0}/`.

### Step 2, `--max-multiplicity 4` (Slurm 32226591, 4 tasks, 4–17 s each)

| case | bench (cap 3) | cap 4, λ = 0.2 | best in the cap-4 pool |
|---|---|---|---|
| ctrl_dup_s708 (4 traversals) | 0.732, 2 species | 0.732, 2 species | 0.994 (`21, 4, 1`), rank 4 at λ 0.2, a same-fit alternative |
| ctrl_dup_s710 (4 traversals) | 0.723, 3 species | 0.723, 2 species | 0.976 (`49, 26, 5`), rank 7 |
| ctrl_dup_s703 (3 traversals, fails) | 0.712, 2 species | **0.971, 1 species** | the same |
| ctrl_dup_s709 (3 traversals, passes) | 0.991 | 0.991 | the same |

λ = 0 gives the same rank 1 in all four.

- **The cap was real but not sufficient.** With cap 4 the true walk enters s708's and s710's pools
  (with cap 3, none had it), but it loses the ranking.
- **Why it loses (s708).** The good answer also carries a 4-copy noise leftover, which is
  significant (cutoff 2.1). So it counts 2 species, like the split rank 1 (`25, 21`). Its 4-traversal
  walk has 57 steps against 41, so λ prefers the split.
- **s703.** It has only 3 traversals, so the cap cannot explain its change. Read it as the larger
  MILP returning a different same-fit answer (tie-breaking), not as a fix.
- **Do not adopt cap 4 bench-wide** (the 40-case rerun) on this evidence.

### Step 3, σ = 0 control (simulation Slurm 32226593, RCK Slurm 32226645)

- Simulations and subgraphs are in `/n/fs/yubi-iw/runs/bench_nested_sigma0/`, in project space,
  not scratch. They are built at `BENCH_NESTED_RCK_PLOIDY=2`, and each config differs from the
  tracked one only in `sigma: 0.0`.
- **Truth check.** The species truth is identical to the bench (`ecdna.tsv`, `events.tsv`, the AA
  truth cycles up to the timestamp line).
  - The long-read observations differ: path supports change, a few low-support paths appear or
    disappear (`walks.walk`), and so do `rck/rck.adjgroups.tsv` groups.
  - So "Run A"'s statement that reads are unchanged at σ = 0 holds for the walks, not for the
    sampled reads.

Scored as `collect.py` scores a bench case (rank 1 at λ = 0.2):

| case | bench (σ = 0.1) | σ = 0 | fails at σ = 0 |
|---|---|---|---|
| ctrl_dup_s703 | 0.712, 2 species | 0.988, 1 species | no |
| nest_cross_1to1_s502 | 0.955, 3 species | 0.972, 2 species | no |
| nest_cross_1to1_s505 | 0.765, 3 species | 0.984, 2 species | no |
| nest_cross_1to1_s508 | 0.672, 3 species | 0.991, 2 species | no |
| nest_samehap_6to1_s801 | 0.492, 1 species | 0.995, 2 species (after the refit; integer rank 1 is `34, 3L, 3L, 1L`) | no |
| nest_samehap_3to1_s301 | 0.693 | 0.492: `12, 1L, 1L`, the 0.992 answer `12, 5, 1` ranks below | **yes** |
| nest_samehap_6to1_s803 | 0.493 | 0.496: `36, 2L, 1L`; best 0.959 is pruned at λ 0.2 | **yes** |
| nest_samehap_6to1_s804 | 0.841 | 0.493: `38, 3L, 2L`; best 0.989 (`40, 6, 1`) ranks below | **yes** |

None of these rank-1 answers is merged.

Reading:

- **Noise: 5 of 8** (s703, s502, s505, s508, s801). Without the correlated depth bias, the current
  ranking picks the right answer. The step-2 noise-tolerant residual on the decomposition did not
  help these, so the fix belongs in the flow stage (A2, depth-bias-aware flow). It is no longer
  deferred.
  - ctrl_dup s708 and s710 probably belong here too: their remaining obstacle at cap 4 is a noise
    leftover.
- **Prior (identifiability): 3 of 8** (s301, s803, s804). Even on noise-free depth, "outer species
  + 1–3 copy leftovers" ties the true nested answer on the fit and wins on species count.
  - Total copy number cannot decide, and these are same-haplotype cases, so alleles cannot either.
    Only long reads spanning the nested species' own junctions could.
  - This feeds step 4 (the 6to1 decision). Expect most of nest_samehap_6to1 to be in this class.
- **Caveat.** These are single runs. s703's change at cap 4 shows that answers move with model size.

## Flow-stage noise, step 0 (2026-11-12, offline): option 1 chosen

`diag/flow_noise.py` (no Gurobi) reads the 32 bench runs and the 8 σ = 0 runs. Its outputs are in
`diag/flow_noise_2026-11-12/`.

- **The noise is in segment depth and is shared by neighbours.**
  - Depth error is 8.7% at σ = 0.1 against 1.9% at σ = 0.
  - The lag-1 correlation between neighbouring segments is 0.60 against 0.26.
  - Junction copy-number error is about 9% at both σ, so it is not the depth wave. It is most
    likely θ = 6.5 fixed (root cause 3) and read sampling.
- **The flow follows the depth.** Where it misses the truth, it is on the observed side in 85% of
  the length. In s703 it lowers the two crossing junctions from 18 to 16 (reads 18.6 and 18.2) to
  match a segment 17% low.
- **The noise also adds subgraph breakpoints** (s703: 8 ecDNA segments against 5). Balance keeps
  copy numbers equal across them.
- **Pricing.** The σ = 0 flow, priced on the σ = 0.1 observations against RCK's σ = 0.1 flow
  (`flow_noise_sigma0flow.tsv`):
  - option 2 (relative-error weights) cuts the noisy flow's advantage by about 10–20% and flips
    no case;
  - option 1 (one bias per 200 kb window at cost κ|β|) cuts it to 0.07–2.3 at κ = 1 (s703 −0.15,
    s803 0.07, s804 0.58, s505 2.25 from 6.93). At κ ≥ 5 it does almost nothing;
  - weighting junction reads ×2 or ×4 makes it worse in 7 of 8.
- **The remaining gap is in the junction term**, so per-case θ is option 1's companion.
- **Caveats.** Totals only; reference-adjacency and background terms left out; windows independent
  (no smoothness). So the gains are an upper bound, and passing cases are not tested.

**Decision:** implement option 1 behind opt-in flags (`--depth-bias-window`, `--depth-bias-cost`),
default off. Then validate on the 8 diagnostic cases, s708 and s710 at cap 4, and passing guards
(Gurobi, about 1 licence-hour; ask first).

## Depth-bias validation (2026-11-12, Slurm 32227536)

`rck-ecdna --depth-bias-window 200000` (commit 13a308a) was run on 16 cases at bias cost and
smoothness 1 (`b1`) and 2 (`b2`): 32 `--keep-pool` runs at λ = 0.2.
- Each case uses its exact bench configuration (`diag/rck_valid.sbatch`, which takes the case from
  `run_case.py`, so demo3 and demo8 keep their options).
- s708 and s710 also use `--max-multiplicity 4`.
- The list and per-run scores are in `diag/rck_depthbias/` (`valid.list`, `valid_scores.tsv`).
- Runs are in `/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_depthbias/`.
- Scoring is `collect.py`'s rule (composite, merged, reported species count), against the bench
  rank 1 at λ = 0.2.

| case | bench | b1 | b2 |
|---|---|---|---|
| ctrl_dup_s703 | 0.712 ✗ | **0.990 ✓** | **0.979 ✓** |
| ctrl_dup_s708 (cap 4) | 0.732 ✗ | **0.999 ✓** | 0.732 ✗ |
| ctrl_dup_s710 (cap 4) | 0.723 ✗ | **0.992 ✓** | 0.731 ✗ |
| nest_cross_1to1_s502 | 0.955 ✗ | **0.974 ✓** | 0.319 ✗ |
| nest_cross_1to1_s508 | 0.672 ✗ | **0.995 ✓** | **0.995 ✓** |
| nest_cross_1to1_s505 | 0.765 ✗ | 0.962 ✗ (3 species for 2: `15, 10, 2`) | **0.968 ✓** |
| nest_samehap_3to1_s301 | 0.693 ✗ | 0.483 ✗ | 0.741 ✗ |
| nest_samehap_6to1_s801 | 0.492 ✗ | 0.493 ✗ | 0.491 ✗ |
| nest_samehap_6to1_s803 | 0.493 ✗ | 0.491 ✗ | 0.494 ✗ |
| nest_samehap_6to1_s804 | 0.841 ✗ | 0.493 ✗ | 0.495 ✗ |
| guards: s709, s506, s201, s306, demo3, demo8 | all pass | all pass (demo3 0.959 → 0.975) | all pass |

Reading:

- **b1 fixes 5 of the 10 targeted failures and breaks no guard.** The fitted biases are about
  −0.23 to +0.28 (s710 reaches the 0.3 bound in one window).
- **s708 and s710 need both changes.** They pass only with the bias and cap 4 together.
- **b2 is unstable.** It fixes s505, but loses s708, s710 and s502.
- **The prior class does not move** (s301, s801, s803), as the σ = 0 control predicted. s801 passed
  at σ = 0 but not with the bias.
- **s804 got worse** (0.841 → about 0.49, nested species lost) in both settings. This is **not** the
  bias absorbing a real copy-number change (checked offline, below). It moved s804 into the prior
  class.
- **Projection.** At b1 with cap 4, the bench would go from 17 to 12 failures, if the 22 untested
  cases hold. That needs the full 40-case rerun.
- **Licence.** Every run succeeded, but Gurobi reported "5 active sessions for a baseline of 2" and
  "over the baseline for 30 minutes". Short tasks (2–20 s) leave sessions counted after they exit.
  Run arrays of short tasks at `%1`, or pause between tasks.

### s804 under the depth bias (offline, 2026-11-12)

Truth: an outer species at 40.8 copies and a nested one at 6.8, so chr3:40.99–42.21 Mb carries 47.7
ecDNA copies plus 3 background.

| chr3 | truth | depth | β (b1) | corrected | bench flow | b1 flow | σ = 0 flow |
|---|---|---|---|---|---|---|---|
| 40.99–41.12 Mb | 47.7 | 56.1 | +0.11 | 49.9 | 48 | 47 | 47 |
| 41.12–41.52 Mb | 47.7 | 65.7 / 66.6 | +0.24 | 50.0 / 50.8 | **62** | **47** | **47** |
| 41.52–42.21 Mb | 47.7 | 45–60 | −0.11 to +0.16 | ≈ 50 | 46–48 | 47 | 47 |

- **The +32% bump is the simulator's depth bias.** The same region reads 49.8 at σ = 0.
- **Without the bias, the bench flow followed it.** It put 62 copies there and used
  `imputed_adj_intra_chr3_10` to balance them. That distorted flow happened to give a 6-copy walk
  over the nested region, so s804 scored 0.841 with a wrong walk.
- **With the bias,** the flow equals the σ = 0 flow (flat 47, no imputed edge), and the nested
  junctions carry 7–8 copies (truth 6.8).
- **Then the prior loses the nested species,** as at σ = 0 (0.493). Rank 1 is `38, 3L, 1L`. The
  pool has 11 answers with the nested species ≥ 0.9. The best is 0.969 (`39, 5, 2, 1`), a same-fit
  alternative, but the λ = 0.2 pruning drops it (11 more steps).
- **So s804 now belongs to the prior class,** with s801, s803 and s301. That class is step 4.

## Full 40-case rerun with the depth bias (Slurm 32228969, submitted 2026-11-12)

The rerun is the bench protocol (`bench_nested.sbatch`: a λ = 0 run and a `--keep-pool` run at
λ = 0.1 per case), with every run given
`--depth-bias-window 200000 --depth-bias-cost 1 --depth-bias-smooth 1 --max-multiplicity 4`.
`BENCH_NESTED_EXTRA` (run_case.py, commit b0302ae) passes these, and each `case.json` records them.

- Array `0-39%1`: one task at a time, because of the licence-session overage seen in the
  validation. Time limit 2 h per task.
- Runs: `/n/fs/yubi-iw/runs/bench_nested_v2/depthbias/rck_runs/<case>/`.
- Logs: `/n/fs/yubi-iw/runs/bench_nested_v2/depthbias/logs/bench-nested_32228969_<task>.out`.
- The bench's own runs and tables are not touched.
- Expected wall time: the bench took 7030 s of RCK time, so about 2–3 h at `%1`.

### How to check it when it has finished

All of this is offline (no Gurobi). Do not loop on `sacct`.

```bash
R=/n/fs/yubi-iw/runs/bench_nested_v2/depthbias
sacct -j 32228969 -X -o JobID%14,State,Elapsed,ExitCode | awk 'NR>2{print $2}' | sort | uniq -c   # 40 COMPLETED?
grep -l 'attempt .* failed' $R/logs/*.out | wc -l               # licence retries (fine if the run then finished)
grep -h 'GurobiError' $R/logs/*.out | sed 's/.*GurobiError/GurobiError/' | sort | uniq -c
# every case has both runs, and the options reached them:
for d in $R/rck_runs/*/; do [ -f $d/lambda_0/summary.json ] && [ -f $d/pool/pool/pool.json ] || echo "MISSING $d"; done
grep -L 'depth-bias-window' $R/rck_runs/*/case.json              # should print nothing
# tables into their own folder (does not touch the repo's tables or README):
cd ~/RCK/results/bench_nested
BENCH_NESTED_OUT=$R/rck_runs BENCH_NESTED_RESULTS=$R/rck_results \
  /n/fs/yubi-iw/conda/envs/gurobi-env/bin/python collect.py | tee $R/collect.log
```

- `collect.py` first prints the re-ranking check against `results/step_sweep/`. The demo cases ran
  with different options, so mismatches there are expected and are not an error.
- To compare per case with the bench at λ = 0.2:

```bash
/n/fs/yubi-iw/conda/envs/gurobi-env/bin/python - <<'PY'
import pandas as pd
R = "/n/fs/yubi-iw/runs/bench_nested_v2/depthbias/rck_results"
old = pd.read_csv("pipeline.tsv", sep="\t").set_index("case")
new = pd.read_csv(R + "/pipeline.tsv", sep="\t").set_index("case")
cols = ["mean_composite", "n_species", "n_truth", "merged", "failure_stage"]
t = old[cols].join(new[cols], lsuffix="_bench", rsuffix="_bias")
t["fixed"] = t.failure_stage_bench.notna() & t.failure_stage_bias.isna()
t["broken"] = t.failure_stage_bench.isna() & t.failure_stage_bias.notna()
print(t.to_string()); print("fixed:", list(t.index[t.fixed])); print("broken:", list(t.index[t.broken]))
PY
```

(The `pipeline.tsv` in the first line is the repo's bench table, read from `results/bench_nested/`.)

What to look for, in this order:

1. **Broken cases (the main risk).** Look for any bench case that passes in the bench and fails
   here, or any control or demo more than 0.005 below its λ = 0 score (the README block's "control
   / demo cases below λ = 0" row).
   - The validation covered only s709, s506, s201, s306, demo3 and demo8 as guards. The other 22
     cases are untested.
   - For each broken case, read `flow_segments.tsv` (`depth_bias`, `obs_corrected`): is a window at
     or near ±0.3, and does the flow now miss the truth where it did not before?
     `diag/flow_noise.py`'s per-segment comparison does this.
   - Note that `--max-multiplicity 4` also changes every pool, so a change can come from the cap,
     not the bias. If it is unclear, rerun that case at cap 3 with `diag/rck_valid.sbatch`.
2. **The failure count at λ = 0.2.** Expected 17 → 12: s703, s708, s710, s502 and s508 fixed.
   s505 still fails on the species count at cost 1. s804 falls to about 0.49 (explained above).
   The prior class (6to1, s301, s307) is unchanged.
3. **The λ rule.** Is 0.2 still the recommended `--step-penalty`, and does λ = 0's mean go up? The
   bias removes noise leftovers, so the step penalty may matter less.
4. **The bias itself.** In each `summary.json` (`depth_bias`: windows, non-zero, min, max), count
   the windows at the ±0.3 bound. A cluster at the bound in one case means the bias is doing more
   than noise correction there.
5. **The scale stage.** `collect.py`'s `scale` check uses the uncorrected `obs_total`. Read
   s104's flag as before.
6. **Runtime.** `collect.py`'s "Compute" line, against 7030 s.

Decision after the check:
- **If no passing case breaks and the failures fall as expected,** propose adopting the options for
  the bench. That means rerunning `collect.py` into the repo tables and README, which needs your
  decision; RCK's defaults stay off.
- **If cases break,** diagnose them as in point 1 before any default or bench change.

## Step 4, offline: an "unexplained junction" prior for the nested species (2026-11-12)

`diag/nested_prior.py` (no Gurobi) re-ranks saved pools. Outputs are in `diag/nested_prior_2026-11-12/`.

**Signature of the prior class.** The wrong rank 1 ("outer + 1–3 copy leftovers") leaves the nested
species' own novel junctions, which carry 6–8 copies of flow, to low_weight walks. The ranking counts
significant species only, so it does not see this.

**Rule `unexpl μ`.** Count the *unexplained junctions* of a decomposition: novel junctions whose
ecDNA flow is at least its low_weight threshold max(1, 0.1 × heaviest), but which no significant
species uses. Add μ per unexplained junction to the structure cost. `np` = without the λ pool pruning.

| pools, λ | rule | mean | failures | fixed | broken / down > 0.005 |
|---|---|---|---|---|---|
| bench (40), 0.2 | default | 0.8531 | 17 | – | – |
| | **unexpl 1** | **0.9113** | **12** | s307, s801, s802, s805, s806 | **none** |
| | unexpl 2 | 0.9179 | 12 | the same | none |
| | unexpl 0.25, 0.5 | 0.8531 | 17 | none | none |
| | mass | 0.8442 | 17 | none | demo3, demo7_totalonly |
| bench (40), 0 | default | 0.7714 | 21 | – | – |
| | unexpl 1 | 0.8325 | 17 | s307, s802, s805, s806 | none |
| | unexpl 1 np | 0.8687 | 16 | + s305 | none |
| depth bias b1 (16), 0.2 | default | 0.8621 | 5 | – | – |
| | **unexpl 1** | **0.9517** | **2** | s301, s801, s803 | **none** (guards s709, s506, s201, s306, demo3, demo8 unchanged) |
| | unexpl 1 np | 0.9815 | 1 | + s804 | none |

Reading:

- **This is the first ranking rule here that fixes cases without breaking any**, on both the noisy
  bench pools and the bias pools, at λ = 0.2 and λ = 0. Every earlier rule traded cases (§3 of the
  failure analysis, the junction-incidence test).
- **μ = 1 is the natural setting.** It makes an unexplained junction cost one species. μ ≤ 0.5 is too
  weak: the right answer pays a species plus about 10 steps at λ 0.2.
- **Fixed cases** are those whose unexplained count at rank 1 is 2–6 (6to1, s307, s301 on the bias
  pools). Cases with 0 (the noise class, s403, s308, s505) are untouched, as intended.
- **s804 needs `np`:** its right answer is pruned at λ = 0.2. An implementation in RCK should prune
  with the new cost (`pool_keeps` reads the structure cost), not after it.
- **Caveats.**
  - These are pools generated under the old ranking, so a solve can only do as well as what the pool
    holds.
  - The bias pools are 16 cases, not 40. s802, s805, s806 and s307 are fixed on the bench pools but
    not yet tested on bias pools.
  - The threshold reuses `--min-species-weight-frac` (0.1).

**Next (offline):** when the full rerun (Slurm 32228969) is checked, run the same script on its 40
pools. That is the real test of depth bias + cap 4 + this rule together:

```bash
cd ~/RCK/results/bench_nested
for L in 0.2 0; do /n/fs/yubi-iw/conda/envs/gurobi-env/bin/python diag/nested_prior.py \
  /n/fs/yubi-iw/runs/bench_nested_v2/depthbias/rck_runs --lambda $L -o /n/fs/yubi-iw/runs/bench_nested_v2/depthbias/nested_prior_lambda$L.tsv; done
```

What to look for: no case broken or down > 0.005 (controls and demos especially), and the count of
remaining failures. Projected: s505, s403, s308, s507 and perhaps s804. If it holds, implement it in
`rck/ecdna/decompose.py` (`ranking_key` / `pool_keeps`) as an opt-in option, e.g.
`--unexplained-junction-penalty 1`, and validate it with a solve.

## Resume here (2026-11-12): state and what to do next

Written at the end of the 2026-11-12 session, for a fresh session. Everything above is the record;
this section is the plan.

### State when the session ended

- **RCK `yubi_dev`.** 9 local commits not pushed: 7c7e4d7 … 17a01d3, plus this note.
  - `rck-ecdna --depth-bias-window` is implemented (13a308a; docs `docs/ecDNA.md` §4.8, tests
    `tests/ecdna/test_depth_bias.py`). It is opt-in, and the defaults are unchanged.
  - Untracked, and left alone on purpose: old `ecbench_*.out` files at the repo root, `logs/`, and
    `results/bench/`.
- **amplicon-allele-analysis.** Unchanged (5447a99). `git fetch` over HTTPS needs credentials.
- **Running.** Slurm **32228969**, the full 40-case rerun with the depth bias
  (`--depth-bias-window 200000 --depth-bias-cost 1 --depth-bias-smooth 1 --max-multiplicity 4`).
  - Array `0-39%1`: the bench protocol, a λ = 0 run and a pool run per case.
  - Started about 14:50; expected 2–3 h.
  - Runs: `/n/fs/yubi-iw/runs/bench_nested_v2/depthbias/rck_runs/`.
  - Logs: `/n/fs/yubi-iw/runs/bench_nested_v2/depthbias/logs/`.
- **Results so far.**
  1. Depth bias at cost 1 fixes s703, s708 and s710 (with cap 4), s502 and s508 in the 16-case
     validation, and breaks no guard ("Depth-bias validation").
  2. The "unexplained junction" prior at μ = 1 fixes s307, s801, s802, s805 and s806 on the bench
     pools, and s301, s801 and s803 on the bias pools, offline, breaking nothing ("Step 4,
     offline").
  - Together they project roughly 17 → about 5 failures. This is not yet tested together.
- **Open item outside this plan.** 380 graph-route Slurm logs exist only on scratch
  (`/n/fs/scratch/$USER/graph_route/h2h/logs`, 13 MB). Scratch was being purged on 2026-11-12.
  - The proposed copy (not run; it needs your approval):
    `tar -czf /n/fs/yubi-iw/runs/graph_route/h2h_logs_2026-11-04.tar.gz -C /n/fs/scratch/$USER/graph_route h2h`.
  - They may already be gone. Check with `ls /n/fs/scratch/$USER/graph_route/h2h/logs | wc -l`.

### Rules that apply to every step

- No Gurobi without an explicit yes from the user (`rck-ecdna`, `sbatch` of any RCK script, the
  Gurobi unit tests). The WLS licence allows 2 sessions.
  - Short tasks leave sessions counted after they exit. The 16-case validation hit "5 active
    sessions" and "over the baseline for 30 minutes".
  - So run arrays at `%1` and check that nothing else is running (`squeue -u $USER`).
- Do not loop `sacct` or `squeue`. Check once, or wait at least 60 s between checks.
- Offline scripts block gurobipy themselves (`sys.modules["gurobipy"] = None`), so they are safe.
  Run them with `/n/fs/yubi-iw/conda/envs/gurobi-env/bin/python`.
- Unit tests without Gurobi:

  ```bash
  cd ~/RCK && python -c "import sys; sys.modules['gurobipy']=None; import unittest; \
  r=unittest.TextTestRunner().run(unittest.defaultTestLoader.discover('tests/ecdna', top_level_dir='.')); sys.exit(not r.wasSuccessful())"
  ```

  At the end of this session: 235 tests, OK, 59 skipped.
- Do not run `collect.py` into the repo folder, or anything else that rewrites the repo's tables or
  README, without the user's decision. Use `BENCH_NESTED_RESULTS=<elsewhere>`.

### Steps, in order

**1. Check the rerun (offline).**
- Is it finished? `sacct -j 32228969 -X -o JobID%14,State,Elapsed | awk 'NR>2{print $2}' | sort | uniq -c`
- Then follow "Full 40-case rerun with the depth bias" → "How to check it when it has finished"
  above. It covers states, licence errors, missing outputs, `collect.py` into
  `/n/fs/yubi-iw/runs/bench_nested_v2/depthbias/rck_results`, the per-case comparison with the
  bench, and what to look for.
- **The main question is whether any passing case breaks.** 22 of the 40 cases were never tested
  with the bias. Expected: 17 → 12 failures at λ = 0.2.
- If a task failed: read its log. Licence errors are retried 5 times inside the task. Resubmit only
  the missing tasks (`--array=<ids>%1`, same environment variables as in the submission above; runs
  that finished are skipped), and only after the user agrees.

**2. Score the rerun's pools with the nested-species prior (offline).**

```bash
cd ~/RCK/results/bench_nested
for L in 0.2 0; do /n/fs/yubi-iw/conda/envs/gurobi-env/bin/python diag/nested_prior.py \
  /n/fs/yubi-iw/runs/bench_nested_v2/depthbias/rck_runs --lambda $L \
  -o /n/fs/yubi-iw/runs/bench_nested_v2/depthbias/nested_prior_lambda$L.tsv \
  | tee /n/fs/yubi-iw/runs/bench_nested_v2/depthbias/nested_prior_lambda$L.txt; done
```

- First check that the `default` row reproduces the rerun's own numbers (its failures and mean at
  λ = 0.2 should match `depthbias/rck_results/pipeline.tsv` and `bench_sweep.tsv`).
- Then read `unexpl 1` and `unexpl 1 np`:
  - no case "broken" and no case "down > 0.005", above all the controls and demos;
  - how many failures remain (projected: s505, s403, s308, s507, perhaps s804);
  - whether `np` (no λ pruning) adds fixes. That says the RCK implementation must apply the rule in
    the pruning too.
- Record the result in this file, under a new heading, and commit.
- **If anything breaks,** look at why (rank-1 weights, unexplained counts in the `-o` TSV, and
  `diag/grow_failures.py <runs> <case>`) before going on.

**3. Implement the prior in RCK as an opt-in option (code and tests; no Gurobi).** Only if step 2
breaks nothing.
- **Parameter.** `DecompositionParams.unexplained_junction_penalty` (default 0 = off), and an
  `rck-ecdna --unexplained-junction-penalty` flag.
- **Definition (as in `diag/nested_prior.py::unexplained`).**
  - An unexplained junction is a NOVEL adjacency whose ecDNA flow (the sum of the flow over its
    `ADJ` edge keys) is at least the decomposition's low_weight threshold (`_flag_low_weight`'s:
    `min_species_weight` if set, else max(1, `min_species_weight_frac` × heaviest weight), on
    integer weights).
  - It is unexplained when no significant (not `low_weight`) species uses it.
  - Species usage: `Species.counts` keys `(ADJ, aid, phasing)`.
- **Ranking.** In `decompose.py`, add penalty × unexplained to the structure cost:
  - in `ranking_key` (a new argument), in `decompose_with_pool`'s candidate loop, and in
    `structure_cost`;
  - **in the pruning too:** `pool_keeps` / `prune_pool` must compare the full structure cost,
    including the new term, not only steps. s804's right answer is otherwise pruned at λ = 0.2.
- **Pools.** Store `unexplained` in `_pool_record`'s entries, so `rank_pool` and
  `rck.ecdna.rerank` can re-rank offline. Old pools without the field count as 0.
- **Ambiguity report.** Check that `same_fit_structures` / `ambiguity_report` (`decided_by`) still
  read sensibly with the new term.
- **Tests without Gurobi** (`tests/ecdna/`, new file):
  - the unexplained count on a hand-made flow and species;
  - `ranking_key` and `rank_pool` order with the penalty on and off;
  - the defaults are unchanged, and the CLI flag parses.
  - Add one `requires_gurobi` solve test, and run it only with permission.
- **Docs.** `docs/ecDNA.md` §5 (ranking), the CHANGELOG, and this file.
- **Check offline.** With the implemented `rank_pool`, re-rank the rerun's pools. The result must
  equal step 2's `unexpl 1` (pruned) or `unexpl 1 np`, depending on the pruning choice. This
  needs no Gurobi, because the pools hold all the entries.

**4. Validation solve (Gurobi; ask first).**
- Use `diag/rck_valid.sbatch`, which runs each case in its exact bench configuration plus extra
  options, as a pool run at λ = 0.2.
- Options on every line:
  `--depth-bias-window 200000 --depth-bias-cost 1 --depth-bias-smooth 1 --max-multiplicity 4 --unexplained-junction-penalty 1`.
- Cases: the prior class (s301, s307, s801–s806), the noise class fixed by the bias (s703, s708,
  s710, s502, s508), and guards (s709, s506, s201, s306, s402, s104, demo3, demo7_totalonly, demo8).
- Run at `%1`. The format is in `diag/rck_depthbias/valid.list`; write the new list to
  `/n/fs/yubi-iw/runs/bench_nested_v2/diag/rck_unexpl/`.
- Score as in "Depth-bias validation". The scoring snippet used `diag/junction_rules.py::Case`; it
  is also how `diag/rck_depthbias/valid_scores.tsv` was made.

**5. Whole-bench run with all options, and the adoption decision (Gurobi; the user decides).**
- Repeat step 1 with `BENCH_NESTED_EXTRA` also holding `--unexplained-junction-penalty 1`, into a
  new `BENCH_NESTED_OUT` (for example `/n/fs/yubi-iw/runs/bench_nested_v2/depthbias_unexpl/rck_runs`).
- Check it the same way.
- **The user decides** whether the bench adopts these options. If so, rerun `collect.py` into the
  repo folder, which rewrites `bench_sweep.tsv`, `pipeline.tsv` and the README results block.
- RCK's defaults stay off unless the user decides otherwise.

**Not in this plan, still open** (see the sections above):
- per-case θ (junction reads per copy), the companion of the depth bias for the junction-term gap
  ("Flow-stage noise, step 0");
- s505 at cost 1 (one 2-copy leftover);
- s403 and s308;
- seed repeats for tie-breaking noise (graph-route step e).

# bench_nested: next steps (written 2026-10-27)

State at the end of the 2026-10-27 session, and a plan for what comes next; progress of the
2026-10-29 session is marked **[done 2026-10-29]** and in "Progress (2026-10-29)" at the end. The
evidence is in `README.md`: the "bench_nested_v2" and "Follow-up checks" sections. The full
old-vs-new comparison is in `comparison.md` and `comparison.tsv` here.

Scratch purges files older than two weeks. Only run outputs remain there (`bench_nested_v2/runs`,
`rck_runs`, `diag`, `subgraphs_paths`); configs, subgraphs and tables are in the repositories, and a
copy of the v2 configs, subgraphs and `rck_results` is in `~/backups/bench_nested_v2_2026-10-29.tar.gz`.

> **Check first (left running on 2026-10-29 at 15:55):** Slurm array **32025976** on ionic, the
> step 2 diagnostic. Tasks 0–14, 16 and 17 are done and recorded. Still running or pending:
> - 15 `nest_cross_3sp_s403_w0.25` and 18 `nest_cross_1to1_s505_w0.5`, running;
> - 19 `nest_cross_3sp_s403_w0.5`, 20 `nest_cross_1to1_s505_w1` and 21 `nest_cross_3sp_s403_w1`,
>   pending.
>
> Expected to finish at about 17:00 on 2026-10-29. Check them with
> `sacct -j 32025976 -X -o JobID%14,State,Elapsed`, then finish `diag/rck_noise/SUMMARY.md` as its
> status section says (`diag/noise_summary.py ... --archive`) and commit. Outputs:
> `/n/fs/scratch/$USER/bench_nested_v2/diag/rck_noise/<name>/`, which scratch deletes after two
> weeks (about 2026-11-12).

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
    - s505 at w0.25 and rel0.1 stays at 0.765, the best in its pool: no good answer enters.
    - s403 at rel0.1 stays at 0.545 (pool best 0.729, not picked).
    - s403 at w0.25, w0.5 and w1, and s505 at w0.5 and w1, were still running or pending when the
      session ended (array 32025976, tasks 15 and 18–21; see "Check first" at the top).
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

## Next steps (revised 2026-10-29)

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


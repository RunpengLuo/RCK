# Step 2 diagnostic: noise-tolerant decomposition, and long-read paths

NEXT_STEPS.md step 2.3 (and the step 3 path runs). Slurm array 32025976 on ionic, submitted
2026-10-29 15:12, 22 tasks at 2 concurrent tasks (WLS licence).

**Status: complete, 22 of 22 runs (recorded 2026-10-31).** `sacct -j 32025976` shows every task
COMPLETED, none near the 2 h Slurm limit (longest 31 min, task 15), and every run has `summary.json`,
`evaluation.tsv` and the kept pool. One run hit a time limit inside RCK: **task 15,
`nest_cross_3sp_s403_w0.25`**. Its first decomposition solve (the pool search) stopped at Gurobi's
`--time-limit 1800` (`Time limit reached`, best bound = best objective, gap 0%; the run took 1866 s),
so the optimum is proven but the pool may be incomplete. Its row below is not filled; its outputs are
archived and described under it. The other 21 runs solved every model to optimality.

`noise_summary.tsv` and `runs/` were refreshed with
`python noise_summary.py rck_noise/ordered.list --archive rck_noise/runs > rck_noise/noise_summary.tsv`
(from `diag/`). The 17 runs archived before are unchanged. Scratch deletes the full outputs (the pool's
`decompositions.tsv` and the Gurobi logs) after two weeks, about 2026-11-12.

## Setup

- **Cases**: the four decomposition failures of bench_nested_v2.
  - `nest_samehap_3to1_s301`, `ctrl_dup_s703` and `nest_cross_1to1_s505` use their committed
    subgraph.
  - `nest_cross_3sp_s403` uses the ploidy-2 subgraph from the 2026-10-27 diagnostics
    (`subgraph_construct.py --rck-ploidy 2`, in `/n/fs/scratch/$USER/bench_nested_v2/diag/s403_ploidy2/`).
  - The path runs use s703 and s301 rebuilt with the `rck_inputs.py` path fix
    (`/n/fs/scratch/$USER/bench_nested_v2/subgraphs_paths/`; same graph, only `paths.tsv` differs).
- **Settings** (`rck-ecdna` options added on 2026-10-29, all off by default):

  | setting | options |
  |---|---|
  | w1 (baseline) | defaults |
  | w0.5 | `--adj-residual-weight 0.5 --seg-residual-scale 0.5 --rank-residual sum` |
  | w0.25 | `--adj-residual-weight 0.25 --seg-residual-scale 0.25 --rank-residual sum` |
  | rel0.1 | `--residual-mode relative --residual-sigma 0.1 --rank-residual sum` |
  | paths | defaults on the path-fixed subgraph (`--path-weight 0.5`) |
  | pw2 | `--path-weight 2` on the path-fixed subgraph |
  | pf0.9 | `--path-fraction 0.9` on the path-fixed subgraph |

- **Every run**:
  - has `--block-flip unoriented`, plus `run_case.py`'s `--no-plots --threads 4 --time-limit 1800`
    and `--truth-dir`;
  - is one `--keep-pool` run at `--step-penalty 0.2`
    (`DIAG_POOL_LAMBDA=0.2 DIAG_LAMBDA0=0 sbatch --array=0-21%2 rck_diag.sbatch`), so rank 1 is the
    λ = 0.2 answer.
- **Files**:
  - `ordered.list` is the submitted list; `noise.list` is the same lines in case order.
  - `noise_summary.tsv` is the full table, from `python ../noise_summary.py ordered.list --archive runs`.
  - `runs/<name>/` holds each run's small outputs: summary, evaluation, flow, the ranked
    decompositions (gzipped), and the pool table and evaluation.
  - `leftover_recount.tsv` (from `../leftover_recount.py`) and `obs_vs_truth.tsv` (from
    `../obs_vs_truth.py`) are the offline checks of 2026-10-31 on the saved pools and flows. No
    new runs; see NEXT_STEPS.md "Assessment (2026-10-31)".
  - The pool's full `decompositions.tsv` (up to 9 MB a run) and the Gurobi logs stay in
    `/n/fs/scratch/$USER/bench_nested_v2/diag/rck_noise/`, which scratch purges after two weeks.

## Results

How to read the columns:
- **species**: significant species at rank 1 / true species.
- **best in pool**: the best mean composite among all the pool's solutions (pool rank and
  significant species).
- **picked**: rank 1 is within 0.005 of the pool's best.
- **inexact**: pool solutions that do not fit the integer flow exactly (residual > 0).

| case | setting | species | composite | best in pool | picked | pool / inexact | runtime |
|---|---|---|---|---|---|---|---|
| s301 | w1 | 3 / 2 | 0.693 | 0.973 (#3, 3 sig.) | no | 60 / 0 | 15 s |
| s301 | w0.5 | 3 / 2 | 0.693 | 0.973 (#3) | no | 62 / 3 | 20 s |
| s301 | w0.25 | 3 / 2 | 0.693 | 0.973 (#6) | no | 64 / 10 | 19 s |
| s301 | rel0.1 | 3 / 2 | 0.693 | 0.973 (#3) | no | 63 / 0 | 14 s |
| s301 | paths | 3 / 2 | 0.693 | 0.973 (#3) | no | 62 / 0 | 13 s |
| s301 | pw2 | 3 / 2 | 0.693 | 0.973 (#3) | no | 60 / 0 | 16 s |
| s301 | pf0.9 | 3 / 2 | 0.693 | 0.973 (#3) | no | 61 / 0 | 13 s |
| s703 | w1 | 2 / 1 | 0.712 | 0.971 (#4, 2 sig.) | no | 37 / 0 | 2 s |
| s703 | w0.5 | 2 / 1 | 0.712 | 0.971 (#4) | no | 37 / 0 | 12 s |
| s703 | w0.25 | 2 / 1 | 0.712 | 0.971 (#4) | no | 37 / 0 | 2 s |
| s703 | rel0.1 | 2 / 1 | 0.712 | 0.971 (#4) | no | 37 / 0 | 2 s |
| s703 | paths | 2 / 1 | 0.712 | 0.971 (#3) | no | 48 / 0 | 3 s |
| s703 | pw2 | 2 / 1 | 0.712 | 0.971 (#3) | no | 47 / 0 | 2 s |
| s703 | pf0.9 | 2 / 1 | 0.712 | 0.971 (#3) | no | 45 / 0 | 2 s |
| s505 | w0.25 | 3 / 2 | 0.765 | 0.765 (#1) | yes | 88 / 11 | 1244 s |
| s505 | rel0.1 | 3 / 2 | 0.765 | 0.765 (#1) | yes | 88 / 11 | 313 s |
| s505 | w0.5 | 3 / 2 | 0.765 | 0.765 (#1) | yes | 88 / 11 | 791 s |
| s505 | w1 | 3 / 2 | 0.765 | 0.765 (#1) | yes | 77 / 0 | 551 s |
| s403 | w0.25 | not filled: Gurobi time limit in the pool search (see below the table) | | | | | 1865 s |
| s403 | rel0.1 | 3 / 3 | 0.545 | 0.729 (#15, 4 sig.) | no | 81 / 0 | 432 s |
| s403 | w0.5 | 3 / 3 | 0.545 | 0.729 (#15, 4 sig.) | no | 82 / 0 | 1376 s |
| s403 | w1 | 3 / 3 | 0.545 | 0.729 (#16, 4 sig.) | no | 85 / 0 | 1104 s |

The s403 w0.25 row is not filled because of the time limit. For the record, its archived outputs
give the same rank 1 as the other three s403 settings (3 / 3, 5.75 / 4.06 / 3.33, 0.545), pool best
0.729 (#14, 4 sig., not picked), 96 / 8 pool solutions / inexact. A rerun without the limit could
only add pool solutions, and none of the other settings' pools holds a better answer.

Before this diagnostic, the composites at λ = 0.2 were:
- s301 0.693, s703 0.712 and s505 0.765 (v2 bench runs);
- s403 0.545 on the ploidy-2 subgraph (2026-10-27 diagnostics).

## Reading

1. **The success test fails: no setting makes s301 or s703 pick its 0.97 answer.** Rank 1 is
   the same at every setting, down to the weights.
2. **The residual weight is not what those two answers lose on.** Every pool solution that wins,
   and both 0.97 answers, fit the integer flow exactly (residual 0).
   - **s703.** The 0.971 answer is the full walk at 16 plus a 2-copy leftover. Rank 1 is two
     circles at 18 and 16. Both count 2 significant species, because the leftover is above 10% of
     the heaviest species. The step penalty then prefers the two circles (17 steps against 31).
     Even at weight 0.25 no inexact fit enters the pool: every pool objective is 3.005–4.008, so
     a single walk's weighted residual is above the window of 2.0.
   - **s301.** The 0.973 answer has 4 significant species (two leftovers of 2.1 and 1.2 copies),
     against 3 at rank 1. (Correction, 2026-10-31: these are the counts after the weight refit.
     The ranking counts significant species from the integer weights: rank 1, pool #0, is
     11 / 3 / [1] with 2 significant species and 30 steps, and the 0.973 answer, pool #2, is
     11 / 4 / 3 / [1] with 3 significant species and 30 steps. It loses by exactly one species.) Lower weights bring 3–10 inexact fits into the pool, but none ranks
     above the exact 3-species fit.
3. **s403 (ploidy-2 subgraph).** Rank 1 is the same at w1, w0.5, rel0.1, and at the time-limited
   w0.25: 0.545, three species at 5.75 / 4.06 / 3.33, as in the 2026-10-27 ploidy-2 run. The
   pool's best is 0.729 at every setting (pool #14–16, 4 significant species), and rank 1 never
   picks it. Only the time-limited w0.25 run lets inexact fits in (8), and none of them is better.
4. **s505.** At every setting (w1, w0.5, w0.25, rel0.1) the answer is the baseline's 0.765
   (12.06 / 7.02 / 2.84), and it is the best in the pool. w0.5 and w0.25 let 11 inexact fits in
   (rel0.1 too); none is better. No good answer is in any s505 pool.
5. **Long-read paths (step 3) change nothing.** With the path fix, s703 loads 10 paths (one with
   ≥ 2 junctions, 3 reads) and s301 15 (three, 5 reads). Rank 1 is unchanged at
   `--path-weight` 0.5 and 2 and at `--path-fraction 0.9`. The wrong rank-1 answers already
   explain 100% of that read support, so these reads cannot tell the answers apart.
6. **Leftover threshold, an offline check (no new runs).** I re-ranked three saved pools at
   λ = 0.2 with `rank_pool`, recounting significant species from the integer weights at higher
   `--min-species-weight-frac`:

   | pool | best in pool | 0.10 (default) | 0.15 | 0.20 | 0.25 |
   |---|---|---|---|---|---|
   | s703 w1 | 0.971 | 0.712 | **0.971** | **0.971** | **0.971** |
   | s301 w1 | 0.973 | 0.693 | 0.693 | 0.703 | 0.703 |
   | s505 w0.25 | 0.765 | 0.765 | 0.528 | 0.532 | 0.532 |

   - A threshold of 0.15 fixes s703, does not fix s301, and breaks s505.
   - Why s505 gets worse (checked 2026-10-31): the default threshold is max(1, 0.1 × heaviest),
     which is 1 copy whenever the heaviest species has at most 10 copies, so then every 1-copy
     species counts as significant. At 0.15, pool #56 (10 / 4 / 1 / 1, four significant species at
     0.10) drops its two 1-copy species and wins with 2 significant species against rank 1's 3
     (12 / 7 / 3 / [1]). Its 1-copy species are not noise: they are 7.0 and 7.5 Mb walks with 7 and
     9 junctions, and its 10-copy species merges both true species (0.528). s502 breaks the same
     way (pool #7, 10 / 1 / 1, 0.955 → 0.513). See NEXT_STEPS.md "Assessment (2026-10-31)".
   - So one global threshold is not the fix either. The pools were collected at 0.10, so a run
     at another threshold could also differ in what reaches the pool.

## Open questions, answered (2026-10-31)

- **Does any weight setting put a better answer into the pool for s403 or s505?** No.
  - s403 (ploidy 2): the best pool composite is 0.729 at every setting (w1 #16, w0.5 #15, rel0.1
    #15, and #14 in the time-limited w0.25 run). Rank 1 never picks it: it stays 0.545.
  - s505: the best pool composite is 0.765 at every setting, and it is rank 1. Nothing better
    enters, not even among the 11 inexact fits at w0.5, w0.25 and rel0.1.
- **Does any setting change rank 1 for either case?** No. The rank-1 weights and composites are
  identical at every setting: s403 5.75 / 4.06 / 3.33 (0.545), s505 12.06 / 7.02 / 2.84 (0.765).
- **Does any result change the reading that the right answers lose on how leftover species are
  counted, not on residual?** It holds for s301 and s703, and the residual part holds for all four.
  It needs one qualification for s505 and s403: they do not lose in the ranking at all.
  - s505's pool has no good answer. Its flow is wrong where the observations are: over about
    300 kb of haplotype B (chr3:121.68–121.98 Mb), depth reads 22–27 copies against a truth of
    30.5 (12–27% low, the simulator's correlated depth bias), and the flow follows it
    (`obs_vs_truth.tsv`; `flow_vs_truth.py`: 14 of 21 segments off by more than 1.5 copies, the
    most of any bench case).
  - In s403 (ploidy 2) the best pool answer is itself only 0.729.
  - And the offline recount (NEXT_STEPS.md "Assessment (2026-10-31)", `leftover_recount.tsv`)
    shows that no recount of leftovers picks the right s301 and s703 answers without breaking
    s502 and s505. For s301 and s703 the choice is between exact fits that the data cannot tell
    apart, so it is a prior-only choice.

## Conclusion

No setting meets the success test, and none changes rank 1 in any of the four cases. Per
NEXT_STEPS step 2.4 there is no regression LIST. The noise-tolerant options stay in RCK, off by
default. In s505 they let inexact fits into the pool, but no better answer. The proposed next
change below was tested offline on 2026-10-31 and is not recommended. See NEXT_STEPS.md
"Assessment (2026-10-31)".

## Proposed next change (not implemented; option 1 tested offline 2026-10-31 and not recommended)

NEXT_STEPS names two larger changes: refit with a noise model before choosing among exact fits,
or put the species weights into the flow MILP. **An offline check of the saved pools argues
against both for these cases.** `data_objective`, the misfit to the observations after each pool
solution's weights are refitted as fractions (`rck/ecdna/refit.py`, already in the pool's
`decompositions.tsv`), gives:

| pool | rank 1 | 0.97 answer | spread over the pool | refit-ranked choice (data + structure cost) |
|---|---|---|---|---|
| s703 w1 | 5.097 | 5.097 | 5.097–5.097 | rank 1 (0.712) |
| s301 w1 | 5.133 | 5.131 | 5.01–5.18 | pool #4 (0.592) |
| s505 w0.25 | 8.024 | (rank 1 is best) | 7.98–8.41 | rank 1 (0.765) |

- **s703.** Every pool solution fits the observations identically. The two circles and the walk
  plus a leftover span the same copy-number space, so the depth and junction data cannot separate
  them, and neither a refit nor decomposing against the observations can.
- **s301.** The 0.973 answer fits better by 0.002 out of 5.1, which is noise, and the refit ranking
  picks a worse answer.
- So the choice in s301 and s703 is made by the priors, and the reads are too few to help (point
  5 above).

Options, cheapest first:

1. **A noise-aware leftover prior.** Count a species as significant only if its weight exceeds
   the copy-number noise of the species it overlaps (for example σ × the heaviest overlapping
   species' weight, σ ≈ 0.15), instead of 10% of the heaviest species overall. s703's
   2-copy leftover inside a 16-copy walk would drop out. What s505 loses at a higher global
   threshold has to be understood first. The offline
   threshold check above shows the risk, so test it offline on the saved pools first, then on the
   18 regression and passing cases.
2. **More long-read evidence.** The bench simulates few long reads (s703: one path with ≥ 2
   junctions). Raising the long-read count in the s703 and s301 configs would show whether paths
   can separate the answers, as CoRAL's do.
3. **Species weights in the flow MILP or a noise-model refit**, as in NEXT_STEPS. They cannot fix
   s703 (identical fits). They may still matter for s505 and s403, where the flow itself may be
   wrong; the pending s403 and s505 runs bear on that.

Both should be judged on the four cases and then on the 8 regression and 10 passing bench cases
(no drop above 0.005), as step 2.5 plans.

# Step 2 diagnostic: noise-tolerant decomposition, and long-read paths

NEXT_STEPS.md step 2.3 (and the step 3 path runs). Slurm array 32025976 on ionic, submitted
2026-10-29 15:12, 22 tasks at 2 concurrent tasks (WLS licence).

**Status: partial, 17 of 22 runs (2026-10-29 15:55); the array was left running when the session ended.**

| task | run | state at 15:55 |
|---|---|---|
| 15 | `nest_cross_3sp_s403_w0.25` | running since 15:34 |
| 18 | `nest_cross_1to1_s505_w0.5` | running since 15:48 |
| 19 | `nest_cross_3sp_s403_w0.5` | pending |
| 20 | `nest_cross_1to1_s505_w1` | pending |
| 21 | `nest_cross_3sp_s403_w1` | pending |

The array should finish at about 17:00 (s403 takes 7–20 min a run, s505 5–25 min). **To finish
this file:**

```bash
sacct -j 32025976 -X -o JobID%14,State,Elapsed           # all 22 COMPLETED?
cd ~/RCK/results/bench_nested/diag
python noise_summary.py rck_noise/ordered.list --archive rck_noise/runs > rck_noise/noise_summary.tsv
```

Then fill the pending rows below from `noise_summary.tsv`, and commit `rck_noise/`. Scratch
deletes the full outputs after two weeks (about 2026-11-12). A task that failed or timed out
(2 h limit) can be rerun alone with `LIST=$PWD/ordered.list DIAG_POOL_LAMBDA=0.2 DIAG_LAMBDA0=0
sbatch --array=<task> --time=02:00:00 ~/RCK/results/bench_nested/diag/rck_diag.sbatch`, submitted
from `/n/fs/scratch/$USER/bench_nested_v2/diag/rck_noise/`.

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
| s505 | w0.5 | pending | | | | | |
| s505 | w1 | pending | | | | | |
| s403 | w0.25 | pending | | | | | |
| s403 | rel0.1 | 3 / 3 | 0.545 | 0.729 (#15, 4 sig.) | no | 81 / 0 | 432 s |
| s403 | w0.5 | pending | | | | | |
| s403 | w1 | pending | | | | | |

Before this diagnostic, the composites at λ = 0.2 were:
- s301 0.693, s703 0.712 and s505 0.765 (v2 bench runs);
- s403 0.545 on the ploidy-2 subgraph (2026-10-27 diagnostics).

## Reading (partial)

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
     against 3 at rank 1. Lower weights bring 3–10 inexact fits into the pool, but none ranks
     above the exact 3-species fit.
3. **s403 (rel0.1).** Rank 1 is unchanged: 0.545, three species at 5.75 / 4.06 / 3.33, the same
   as the 2026-10-27 ploidy-2 run. The pool's best is still 0.729 (4 significant species), and no
   inexact fit enters. The three weight settings are pending.
4. **s505.** At w0.25 and rel0.1 the answer is still the baseline's 0.765, and it is the best in
   the pool. The 11 inexact fits that now enter are no better: no good answer enters the pool.
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
   - s505's rank 1 changes to a worse pool solution at 0.15; why has not been checked yet. Its
     2.84-copy species is the one that carries both true species' junctions (2026-10-27 reading).
   - So one global threshold is not the fix either. The pools were collected at 0.10, so a run
     at another threshold could also differ in what reaches the pool.

## Conclusion so far

No setting meets the success test. Per NEXT_STEPS step 2.4 there is no regression LIST; the
next change is proposed below, not implemented. The pending s403 and s505 runs can only change the reading
for those cases. The noise-tolerant options stay in RCK, off by default. For s505
they let inexact fits into the pool, but no better answer so far.

## Proposed next change (not implemented)

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

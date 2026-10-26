# bench_nested: `--step-penalty` choice and nested-ecDNA performance

A minimum benchmark of amplicon-allele-analysis `simulations/subgraph_construct.py` → `rck-ecdna`
with two uses: choosing a default for `--step-penalty` (λ), and measuring end to end how nested
ecDNA species and haplotype-crossing junctions are recovered. It sits alongside the older
`results/bench/` and `results/step_sweep/`, which are unchanged.

## Cases

Bench cases: 7 scenarios × 2 accepted seeds = 14 simulations, from amplicon-allele-analysis
`simulations/make_scenarios.py --set nested` (demo7's config with one or two keys changed: 3
clones, 1 tumor sample, long-read counts, no reads; library size, read model and bin size as
demo7, so junction reads per copy stay 6.5).

| scenario | species | nesting | haplotypes | bulk ratio | purpose |
|---|---|---|---|---|---|
| `nest_cross_2to1` | 2 | nested | `p_cross` 0.2 | 2:1 | demo7-like, new seeds |
| `nest_samehap_2to1` | 2 | nested | `nest_same_haplotype`, `p_cross` 0 | 2:1 | merging saves a species (untested before) |
| `nest_samehap_3to1` | 2 | nested | `nest_same_haplotype`, `p_cross` 0 | 3:1 | demo8-like, a tie decided by λ |
| `nest_cross_3sp` | 3 | 2 inside 1 | `p_cross` 0.2 | 0.5 / 0.3 / 0.2 | several nested species with crossings |
| `nest_cross_1to1` | 2 | nested | `p_cross` 0.2 | 0.55 / 0.45 (1.22:1) | close, non-integer copy numbers |
| `ctrl_separate` | 2 | none (no shared base) | `p_cross` 0.2 | 1:1 | demo6-like control: λ must not hurt |
| `ctrl_dup` | 1 | – | `p_cross` 0.2 | – | true duplications and fold-backs (`sv_probs` dup 0.8, fback 0.3, del 0): λ must not split it |

Each scenario draws its seeds from its own block (101–110, 201–210, …, 701–710), because the
walks depend on the seed alone: with shared seeds, every nested scenario would have the same
geometry. A seed is rejected, recorded, and the next one tried (at most 5 tries per accepted
seed) when it misses the scenario's intent: nested species ≥ 80% of their bases inside the
outer species; ≥ 1 crossing ecDNA junction (`cross`); no crossing and every species on one
haplotype (`samehap`); no shared base (`ctrl_separate`); a segment traversed twice
(`ctrl_dup`); and for every scenario the bulk copies within 10% of the target ratio and every
ecDNA junction with ≥ 20 reads. Records: amplicon-allele-analysis
`simulations/configs/bench_nested/seeds.tsv`.

Regression cases: the 8 cases of `results/step_sweep/` (demo3, demo4 amplicon_0 / amplicon_1,
demo5, demo6 union, demo7 with `--block-flip`, demo7 total-only, demo8), with the same
subgraphs and options.

## Runs

Every case gets two `rck-ecdna` runs (`run_case.py`), each from the subgraph folder alone with
`--truth-dir`, `--no-plots`, `--threads 4` and `--time-limit 1800`; bench cases add
`--block-flip unoriented`:

- `runs/<case>/lambda_0/`: λ = 0, the current default;
- `runs/<case>/pool/`: λ = 0.1 with `--keep-pool`. λ only enters the ranking, and every λ > 0
  widens the decomposition pool to the same `--ambiguity-window`, so this one pool holds every
  answer λ > 0 can give. `python -m rck.ecdna.rerank runs/<case>/pool` re-ranks it for any λ
  with the run's own pruning and ranking functions, and reports the composites (from
  `pool/evaluation.tsv`) and the ambiguity fields (from the same-fit alternatives recorded for
  each rank-1 solution of the λ grid).

λ = 0 is run rather than re-ranked: a λ = 0 run solves a narrower pool (10% relative gap, 40
solutions), where ties can be broken differently. demo3 shows it in `results/step_sweep/`:
rank 1 at 12.77 copies at λ = 0, and an equivalent solution at 12.70 for every λ > 0. The λ = 0
runs cost about a quarter of the pool runs (demo8: 8 s against 31 s).

λ grid: 0, 0.02, 0.05, 0.1, 0.2, 0.5.

## How to run (Princeton ionic)

```bash
# 1. simulations, seed checks, subgraphs and their checks (no Gurobi; 7 tasks, one per scenario)
cd ~/amplicon-allele-analysis/simulations && mkdir -p logs && sbatch bench_nested.sbatch

# 2. first check the re-ranking on the step sweep's hardest cases (demo3, demo7_totalonly, demo8)
cd ~/RCK/results/bench_nested && mkdir -p logs && sbatch --array=0,6,7%2 bench_nested.sbatch
python collect.py --validate          # rerank_validation.tsv: every "pool" row should match

# 3. every case (finished runs are skipped)
sbatch bench_nested.sbatch            # 22 tasks, 2 at a time (academic WLS licence)
python collect.py                     # bench_sweep.tsv, pipeline.tsv, results below
```

Both arrays skip every step whose outputs exist, so a failed or interrupted array can be
resubmitted as is. Estimated cost: simulation about 1.5 min and subgraph about 1 min per seed
(≈ 14–25 simulations with rejections, ≈ 0.5–1 CPU-hour in all); RCK about 10–60 s per λ = 0
run and 30 s–3 min per pool run with 4 threads (22 cases, ≈ 1–3 CPU-hours, about 30 min of wall
time at 2 concurrent tasks).

## Outputs

- `rerank_validation.tsv`: rank-1 composites and ambiguity fields, re-ranked (and λ = 0 run),
  against `results/step_sweep/sweep.tsv`, per case and λ.
- `bench_sweep.tsv`: one row per case and λ: rank-1 significant species against the truth
  count, weights, steps, residual, composite per truth species and their mean, `merged`,
  `ambiguous`, `decided_by`, `gap`, `n_alternatives`, runtime (`source` is `run` for λ = 0 and
  `pool` for re-ranked λ, whose runtime is the pool run's). `merged`: an inferred species carries
  two or more truth species, or a truth species is missed while another inferred species carries
  it; "carries" means junction recall ≥ 0.5, because nested species share their segments by
  construction but no junction.
- `pipeline.tsv`: one row per bench case at the chosen λ: subgraph scores (segments and
  junctions recovered, ecDNA phase blocks oriented right, wrong or not at all, the junction reads
  per copy `rck_inputs.py` would have fitted against the 6.5 used), the flow objective, imputed
  edges used, block flips against the truth orientation, the RCK scores, and for a failing case
  the stage that explains it (`subgraph` when the subgraph misses a true junction or covers under
  80% of a true segment, else `RCK`).
- `runs/`: every run's full output (not tracked).

The chosen λ is the smallest on the plateau (within 0.001) of best mean composite over all
cases, with no control (`ctrl_*`) or demo case more than 0.005 below its λ = 0 score. It is a
recommendation; RCK's default `--step-penalty` stays 0.

## Status and reading (2026-10-25)

Run on ionic (Slurm 31981878/31981889 validation, 31981902/31981984 RCK). 20 of 22 cases are
complete. `nest_cross_3sp_s406` and `nest_cross_1to1_s502` were cancelled at the session's time
limit. The first was still at an 81% MIP gap after 800 s in its λ = 0 decomposition; the second
finished λ = 0 (488 s) but not its pool run. Both need a longer Slurm limit (`--time=05:00:00`)
and are excluded from the choice below until then. The RCK runs live on scratch
(`runs` → `/n/fs/scratch/$USER/bench_nested/rck_runs`; a 5 GB home quota cannot hold them).

The rule below picks 0.5 because of a single case, `nest_samehap_3to1_s301`, which stays merged
and wrong at every λ (0.693 → 0.740). It is also the edge of the grid. Every other case gives
identical answers from λ = 0.05 to 0.5. **λ = 0.05 is the smallest value that gives the large
gain** (mean composite 0.800 → 0.896, merges fixed on `nest_cross_2to1`, `nest_samehap_2to1`
and demo7 total-only), with no control or demo case worse than at λ = 0. The earlier
recommendation of 0.1 sits on the same plateau.

## Results

<!-- results:start -->
**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.5**, the smallest λ on the plateau of best mean composite over all 20 cases (0.5), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.5. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.8002 | 0.8026 | 0.8958 | 0.8958 | 0.8958 | 0.8981 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.5 | merged at 0 | merged at 0.5 |
|---|---|---|---|---|---|
| ctrl_dup | 2 | 0.870 | 0.870 | 0 | 0 |
| ctrl_separate | 2 | 0.670 | 0.670 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 1 | 0.716 | 0.729 | 0 | 0 |
| nest_cross_2to1 | 2 | 0.564 | 0.972 | 2 | 0 |
| nest_cross_3sp | 1 | 0.596 | 0.632 | 1 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 2 | 0.836 | 0.859 | 1 | 1 |

Remaining failures at λ = 0.5 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 6 of 12 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 4 species for 2 true; composites ecdna1=0.962,ecdna2=0.517 |
| nest_cross_3sp_s402 | RCK | merged: ec3 carries ecdna1+ecdna3; composites ecdna1=0.628,ecdna2=0.691,ecdna3=0.577 |
| nest_cross_1to1_s501 | RCK | composites ecdna1=0.893,ecdna2=0.564 |
| ctrl_separate_s601 | subgraph | missing junctions none; min segment coverage 0.0% |
| ctrl_separate_s602 | subgraph | missing junctions none; min segment coverage 0.0% |
| ctrl_dup_s702 | RCK | 1 ecDNA phase block(s) wrong after the flow; imputed edges imputed_adj_intra_chr16_6; composites ecdna1=0.744 |

Subgraph stage: 2 of 12 subgraphs miss a true junction or cover under 80% of a true segment (ctrl_separate_s601, ctrl_separate_s602).

Compute: 754 s of RCK wall time over 20 cases (one λ = 0 run and one pool run each, 4 threads).
<!-- results:end -->

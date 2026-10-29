# bench_nested: `--step-penalty` choice and nested-ecDNA performance

A minimum benchmark of amplicon-allele-analysis `simulations/subgraph_construct.py` → `rck-ecdna`
with two uses: choosing a default for `--step-penalty` (λ), and measuring end to end how nested
ecDNA species and haplotype-crossing junctions are recovered. It sits alongside the older
`results/bench/` and `results/step_sweep/`, which are unchanged.

**Benchmark in use: bench_nested_v2** (inherited haplotype labels; see "Inherited haplotype labels"
and "bench_nested_v2" below). It replaced the first benchmark on 2026-10-29. The tracked configs
(with `seeds.tsv`) and subgraphs in amplicon-allele-analysis, and this folder's `bench_sweep.tsv`,
`pipeline.tsv` and `rerank_validation.tsv`, are v2's. `comparison.md` and `comparison.tsv` compare
the two benchmarks. The first benchmark's results are kept below ("First benchmark"), and its
tables in git history (9b74b87).

Run outputs stay on scratch, which **purges files older than two weeks**: `runs` →
`/n/fs/scratch/$USER/bench_nested_v2/rck_runs` (the 8 regression cases there link to
`/n/fs/scratch/$USER/bench_nested/rck_runs`), amplicon-allele-analysis `simulations/runs/bench_nested`
→ `/n/fs/scratch/$USER/bench_nested_v2/runs`, and the diagnostics in
`/n/fs/scratch/$USER/bench_nested_v2/diag/`. Everything this README or `NEXT_STEPS.md` needs to
keep is in the repositories; the simulations can be regenerated from the configs with
`bench_nested.py`, and the RCK runs with `bench_nested.sbatch`.

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

### Inherited haplotype labels (bench_nested_v2)

The cases above were simulated with per-step haplotype labels, which give some junctions
several phasings and some fold-backs a haplotype crossing. Neither arises from one physical
join, and they caused part of the `ctrl_dup_s702` and `nest_cross_3sp_s402` failures. See
amplicon-allele-analysis `CHANGELOG.md`, `hap_labels`. The benchmark was rebuilt with
inherited labels (these are now the tracked configs and subgraphs, so steps 1–3 above run it; the
commands below are how it was first built and run, outside the tracked folders):
`/n/fs/scratch/$USER/bench_nested_v2/{configs,runs,subgraphs}`, 14 accepted seeds. Its RCK runs
go to their own folders. The regression cases are unchanged, so their existing runs are linked
in rather than rerun:

```bash
cd ~/RCK/results/bench_nested && mkdir -p logs
export BENCH_NESTED_ROOT=/n/fs/scratch/$USER/bench_nested_v2      # cases: its seeds.tsv and subgraphs
export BENCH_NESTED_OUT=$BENCH_NESTED_ROOT/rck_runs                # runs
export BENCH_NESTED_RESULTS=$BENCH_NESTED_ROOT/rck_results         # collect.py tables and README
mkdir -p $BENCH_NESTED_OUT
for c in demo3 demo4_amp0 demo4_amp1 demo5 demo6_union demo7_blockflip demo7_totalonly demo8; do
  ln -sfn /n/fs/scratch/$USER/bench_nested/rck_runs/$c $BENCH_NESTED_OUT/$c
done
sbatch --array=8-21%2 bench_nested.sbatch                          # the 14 bench cases (sbatch passes the environment)
python collect.py                                                  # afterwards, with the same three variables set
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
  the stage that explains it, the first that applies: `subgraph` when the subgraph misses a true
  junction or covers under 80% of a true segment; `scale` when its depth per copy is more than 30%
  from the cases' median, or the flow's segment and junction copy numbers against the truth differ
  by more than 1.5× (`diag/scale_ratio.py`); else `RCK`. The scale columns are in the table.
- `runs/`: every run's full output (not tracked).

The chosen λ is the smallest on the plateau of best mean composite over all cases, with no
control (`ctrl_*`) or demo case more than 0.005 below its λ = 0 score. A λ is on the plateau when
its mean is within 0.001 of the best λ's, or when fewer than 2 cases gain more than 0.005 from it
to the best λ, so that one case at the edge of the grid cannot pick λ (2026-10-29; before, the
0.001 test alone picked 0.5 on both benchmarks because of `nest_samehap_3to1_s301`). It is a
recommendation; RCK's default `--step-penalty` stays 0. `python collect.py --lambda X` reports
another λ.

## First benchmark (2026-10-25; replaced by bench_nested_v2)

### Status and reading

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

### Results of the first benchmark

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

## Results (bench_nested_v2)

Written by `python collect.py` (the λ rule and failure stages of its docstring).

<!-- results:start -->
**Re-ranking check** (`rerank_validation.tsv`, against `results/step_sweep/sweep.tsv`): re-ranked λ > 0 composites match at 32/32 (case, λ), ambiguity fields at 32/32; the λ = 0 runs match at 8/8; λ = 0 emulated from the pool matches at 7/8.

**Recommended `--step-penalty`: 0.2**, the smallest λ on the plateau of best mean composite over all 22 cases (0.2, 0.5; a λ whose shortfall to the best comes from fewer than 2 cases is on it), with no control or demo case more than 0.005 below its λ = 0 score. λ giving the same rank-1 answers on every case: 0.2. The default stays 0.

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all cases | 0.8269 | 0.8600 | 0.9076 | 0.9076 | 0.9283 | 0.9304 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |
| cases gaining from λ to the best λ | demo7_totalonly, nest_cross_2to1_s104, nest_cross_3sp_s402, nest_cross_3sp_s403, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | demo7_totalonly, nest_cross_3sp_s402, nest_samehap_2to1_s202, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_cross_3sp_s402, nest_samehap_3to1_s301 | nest_samehap_3to1_s301 | none |

Mean composite per scenario (mean over seeds; merged = rank-1 decompositions flagged merged):

| scenario | cases | λ = 0 | λ = 0.2 | merged at 0 | merged at 0.2 |
|---|---|---|---|---|---|
| ctrl_dup | 2 | 0.855 | 0.855 | 0 | 0 |
| ctrl_separate | 2 | 0.979 | 0.979 | 0 | 0 |
| demo3 | 1 | 0.960 | 0.959 | 0 | 0 |
| demo4_amp0 | 1 | 0.983 | 0.983 | 0 | 0 |
| demo4_amp1 | 1 | 0.997 | 0.997 | 0 | 0 |
| demo5 | 1 | 0.990 | 0.990 | 0 | 0 |
| demo6_union | 1 | 0.990 | 0.990 | 0 | 0 |
| demo7_blockflip | 1 | 0.997 | 0.997 | 0 | 0 |
| demo7_totalonly | 1 | 0.422 | 0.976 | 1 | 0 |
| demo8 | 1 | 0.987 | 0.987 | 0 | 0 |
| nest_cross_1to1 | 2 | 0.860 | 0.860 | 2 | 2 |
| nest_cross_2to1 | 2 | 0.609 | 0.967 | 1 | 0 |
| nest_cross_3sp | 2 | 0.550 | 0.785 | 2 | 1 |
| nest_samehap_2to1 | 2 | 0.745 | 0.991 | 1 | 0 |
| nest_samehap_3to1 | 2 | 0.836 | 0.836 | 1 | 1 |

Remaining failures at λ = 0.2 (mean composite < 0.95, merged, or a wrong species count; `pipeline.tsv`): 5 of 14 bench cases.

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | merged: ec2 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.964,ecdna2=0.423 |
| nest_cross_3sp_s403 | scale | depth per copy 6.51 is -53% from the median 13.74; segment / junction scale 2.12 (segment CN 2.20x, junction CN 1.04x truth); composites ecdna1=0.635,ecdna2=0.483,ecdna3=0.657 |
| nest_cross_1to1_s502 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.967,ecdna2=0.943 |
| nest_cross_1to1_s505 | RCK | merged: ec3 carries ecdna1+ecdna2; 3 species for 2 true; composites ecdna1=0.944,ecdna2=0.586 |
| ctrl_dup_s703 | RCK | 2 species for 1 true; composites ecdna1=0.712 |

Scale: 2 of 14 cases have depth per copy more than 30% from the median or segment / junction scales apart by more than 1.5x (nest_cross_2to1_s104, nest_cross_3sp_s403).
Subgraph stage: 0 of 14 subgraphs miss a true junction or cover under 80% of a true segment.

Compute: 4349 s of RCK wall time over 22 cases (one λ = 0 run and one pool run each, 4 threads).
<!-- results:end -->

## bench_nested_v2 (2026-10-27; the benchmark since 2026-10-29)

These are the results of the rebuilt benchmark (inherited haplotype labels, see above), run on ionic
(Slurm 32000320, the 14 bench cases; the 8 regression cases are links to the runs above). They to `/n/fs/scratch/$USER/bench_nested_v2/rck_results/` and are now this folder's tables,
with `comparison.md` and `comparison.tsv` against the first benchmark's. The rest of this section
is the 2026-10-27 reading, under the old λ rule and failure stages; "Results (bench_nested_v2)"
above is the current one.

The two benchmarks share scenarios and seed blocks. 9 bench cases keep their name: the 4 samehap
cases are identical, and the 5 cross and control cases have the same walks and copy numbers but
new haplotype labels. 5 are replacement seeds, because s101, s501, s602 and s701 have no crossing
under inherited labels and s403 now passes its read check: s104, s505, s603, s703 and s403
(in place of s406).

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, all 22 cases | 0.8269 | 0.8600 | 0.9076 | 0.9076 | 0.9283 | 0.9304 |
| control / demo cases below λ = 0 | none | none | none | none | none | none |

- `ctrl_dup_s702` is fixed: 0.744 → 0.997 at every λ. It has no wrong junction label, no imputed
  edge and no ecDNA block wrong after the flow.
- `ctrl_separate_s601` is fixed by `--max-extend`: 0.759 → 0.982. Its replacement pair, s603, scores
  0.976.
- `nest_cross_3sp_s402` has clean labels, flips and imputed edges. It is merged at λ ≤ 0.1
  (ecdna2 + ecdna3, 2 species for 3, 0.523), and recovered at λ ≥ 0.2 (3 species, 0.978).
- **The λ = 0.05–0.5 plateau no longer holds**, and 0.1 is no longer on it. λ = 0.2 is the smallest
  grid value with every real gain. The rule picks 0.5 again, for `nest_samehap_3to1_s301` alone
  (0.693 → 0.740, still merged). The choice of 0.2 rests on one case. Re-ranking the s402 pool on a
  0.01 grid puts its switch between 0.10 and 0.11: at 0.1, one extra species (penalty 1.0) exactly
  ties 10 fewer steps, so 0.1 sits on the tie.
- `nest_cross_1to1_s502` finishes (13 s + 32 s). `nest_cross_3sp_s406` is replaced by s403.

Remaining failures at λ = 0.1, 6 of 14 bench cases (5 at λ = 0.2):

| case | stage | why |
|---|---|---|
| nest_samehap_3to1_s301 | RCK | unchanged: ec2 carries ecdna1+ecdna2; 3 species for 2; 0.693 |
| nest_cross_3sp_s402 | RCK (λ) | ec2 carries ecdna2+ecdna3; 2 species for 3; 0.523; fixed at λ ≥ 0.2 |
| nest_cross_3sp_s403 | subgraph (scale), then RCK | `rck_inputs` ploidy fit 5.81 (depth per copy 6.51, about 13.8 elsewhere): segment CN 2× truth, flow objective 41, 181 steps; 0.592; `pipeline.tsv` says RCK. Rebuilt with `--rck-ploidy 2`, the scale is right (segment and junction CN 1.03 / 1.04× truth), but it still fails (0.545); see Follow-up checks |
| nest_cross_1to1_s502 | RCK | 3 species for 2: a 1.48-copy species carries both species' junctions; 0.955 |
| nest_cross_1to1_s505 | RCK | the same, with 2.84 copies; ecdna2 0.586 |
| ctrl_dup_s703 | RCK | 2 species for 1: the B insertion's observed CN (16.1) is below the truth (19.2) |

The old `nest_cross_1to1_s501`, `s502` and `nest_cross_3sp_s406` subgraphs had the same scale error
(ploidy ≈ 4, depth per copy 7.4–7.9). s501's old failure therefore comes from the subgraph stage
(its old segment CN was 1.81× the truth), not from RCK. `pipeline.tsv`'s stage rule checks coverage and
missing junctions only.

RCK wall time: 4349 s over 22 cases, 2860 s of it in `nest_cross_3sp_s403`.

### Follow-up checks (bench_nested_v2, 2026-10-27)

These are diagnostic runs on the existing outputs and subgraphs, with no new simulation, in
`/n/fs/scratch/$USER/bench_nested_v2/diag/` (Slurm 32005476 and 32005494; scripts `rck_diag.sbatch`,
`pool_summary.py`, `flow_vs_truth.py`, `scale_ratio.py`, `summarize.sh`, kept in `diag/` here). The plan
that follows from them is in `NEXT_STEPS.md`. "Fitted θ" is the
junction-reads-per-copy value `rck_inputs` fits; the bench fixes it at 6.5.

| case | check | result |
|---|---|---|
| nest_cross_3sp_s402 | re-rank the pool, λ = 0.10–0.20 in steps of 0.01 | 3 species, 0.978, from λ = 0.11 on |
| nest_cross_3sp_s403 | subgraph with `--rck-ploidy 2`, then RCK | scale fixed, still 0.545 at every λ; decomposition optimal; best pool solution 0.729 |
| nest_cross_1to1_s502 | RCK with fitted θ 5.11 | **resolved**: 2 species, 0.963 at λ ≥ 0.02 (junctions observed at 0.77× truth with θ 6.5) |
| nest_cross_1to1_s505 | RCK with fitted θ 4.93 | worse, 0.600 (θ 6.5 was right: junctions at 1.00× truth); no good pool solution with θ 6.5 (best = rank 1, 0.765) |
| nest_samehap_3to1_s301 | pool, and RCK with fitted θ 5.09 | the pool holds a 0.973 answer (4 species, 2 of them leftovers ≤ 2.1 copies), ranked below the 3-species rank 1; fitted θ: 0.691 |
| ctrl_dup_s703 | window depth; pool; RCK with fitted θ 4.53 | the B insertion's depth per copy is 12.0 against 12.9–16.4 elsewhere (the simulator's `sigma` 0.1 correlated depth bias); the pool holds a 0.971 answer (the full walk at 16.2 plus a 1.7 leftover), ranked below rank 1; fitted θ: 0.718 |

Reading:

- s402 is a λ choice: any λ ≥ 0.11 resolves it.
- s502 is a junction-scale input error. The fitted θ fixes it here, but it is wrong in s505 (4.93
  against a true ≈ 6.5) and s403 (2.12), so it cannot replace 6.5 as is.
- s301, s703, s505 and s403 (with the right scale) are decomposition failures. The ecDNA flow is
  integer and fitted to noisy observations, and the decomposition must cover it exactly: residual
  weight 1 per copy on a segment ≥ 50 kb (`adj_residual_weight` 1.0 on junctions), against
  `species_penalty` 1. A 1–3 copy dip from the simulated depth bias (±10% per 200 kb, so ±2 copies
  at 20) therefore costs as much as an extra species. Exact-fit alternatives then win: extra
  leftover species (s703, s301) or wrong merges (s505, s403). In s301 and s703 the near-correct answer
  is in the pool but loses on the ranking. In s505 and s403 it is not in the pool at all.

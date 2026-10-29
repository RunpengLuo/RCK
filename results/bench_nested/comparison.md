# bench_nested: old (committed) vs new (bench_nested_v2) benchmark

- **Old**: RCK `results/bench_nested/{bench_sweep,pipeline}.tsv` (2026-10-25). Its simulations used per-step
  haplotype labels and `p_cross` 0.2. Its subgraphs are in amplicon-allele-analysis
  `simulations/subgraphs/bench_nested/`. Its RCK runs are in `/n/fs/scratch/ym3292/bench_nested/rck_runs`.
- **New**: `/n/fs/scratch/ym3292/bench_nested_v2` (2026-10-27). It uses `hap_labels: inherit`, `p_cross` 0.5 (0
  for samehap), `--max-extend` 5 Mb (the default) and the new `evaluate.py`, which scores a junction against
  every truth phasing. Its 14 bench cases were run with `sbatch --array=8-21%2` (job 32000320). The 8
  regression cases (demo3–demo8) are links to the old runs.
- Per-case table: `comparison.tsv` (one row per matched pair and λ ∈ {0, 0.1}). It includes the
  simulation change, species counts, composites, the merged flag, pass/fail, the pipeline stage and reason,
  junction-label errors, imputed edges, phase blocks wrong after the flow, the fitted depth per copy, and
  runtimes.

A case **fails** when its mean composite is < 0.95, it is merged, or its species count is wrong (the
`pipeline.tsv` rule).

## Summary

- Fixing the haplotype labels resolves `ctrl_dup_s702` (0.744 → 0.997). It also removes the junction-label
  error in `nest_cross_3sp_s402`.
- `nest_cross_3sp_s402` now succeeds at λ ≥ 0.2 (0.978). At λ = 0.05–0.1 it is merged (2 species for 3;
  0.523). **The λ = 0.05–0.5 plateau no longer holds, and 0.1 is no longer on it.**
- **New finding, not in the old analysis:** the depth-per-copy (purity/ploidy) fit in `rck_inputs` sometimes
  picks ploidy ≈ 4 or 6 instead of 2. That doubles every segment copy number against the junction reads.
  This happened in old `nest_cross_1to1_s501`, `s502` and `nest_cross_3sp_s406`, and it happens in new
  `nest_cross_3sp_s403`. It is a subgraph-stage failure, but `pipeline.tsv` labels it "RCK" because its
  stage rule only checks coverage and missing junctions.
- `ctrl_separate_s601` is fixed by `--max-extend` (0.759 → 0.982).
- At λ = 0.1 the new benchmark has 6 failures out of 14 bench cases. At λ = 0.2 and 0.5 it has 5.

## 1. Case matching

The new benchmark has the same 7 scenarios and seed blocks (101–110 … 701–710), 2 accepted seeds per
scenario, and the same acceptance rules. With `inherit` the walks draw the same random numbers, so **a
given seed has the same walks, segments, total copy numbers and clone proportions as before**. This was
checked on the `ecdna` records of `walks.walk`, on `segments.tsv` totals, on `proportions.tsv` and on
`ecdna.tsv`. What changes for crossing scenarios is the haplotype labels: allele copy numbers, junction
phasings and therefore reads. That changes which seeds pass the "≥ 1 crossing junction" and "≥ 20
junction reads" checks.

| old case | new case | match | simulation |
|---|---|---|---|
| demo3 … demo8 (8) | same | name | unchanged; old runs linked in, results identical |
| nest_samehap_2to1_s201, s202 | same | name | identical (truth, observations, subgraph) |
| nest_samehap_3to1_s301, s304 | same | name | identical |
| nest_cross_2to1_s103 | same | name | relabelled (9 of 87 segments change allele CN); crossings 6 → 4 |
| nest_cross_3sp_s402 | same | name | relabelled (8/82); crossings 5 → 2; ecdna1=A, ecdna2=mixed, ecdna3=B |
| nest_cross_1to1_s502 | same | name | relabelled (13/73); crossings 4 → 6 |
| ctrl_separate_s601 | same | name | relabelled (7/84); crossings 4 → 6 |
| ctrl_dup_s702 | same | name | relabelled (5/70); the AA+BB junction (chr16 58.78 → 60.12 Mb) now has one phasing |
| nest_cross_2to1_s101 | nest_cross_2to1_s104 | scenario | s101 rejected in v2: no crossing junction |
| nest_cross_3sp_s406 | nest_cross_3sp_s403 | scenario | s403 now passes (old: a junction with 17 reads); s406 not reached |
| nest_cross_1to1_s501 | nest_cross_1to1_s505 | scenario | s501 rejected in v2: no crossing (s503, s504 also rejected) |
| ctrl_separate_s602 | ctrl_separate_s603 | scenario | s602 rejected in v2: no crossing |
| ctrl_dup_s701 | ctrl_dup_s703 | scenario | s701 rejected in v2: no crossing |

In the old truth, two junctions had several phasings: `ctrl_dup_s702` (AA+BB) and `nest_cross_3sp_s402`
(AA+BA). The new truth has none (`truth/rck/rck.acnt.tsv`, all cases).

Of the 14 bench cases, 9 keep their names (4 identical, 5 relabelled) and 5 are new seeds. Bench-level
means are therefore not directly comparable. Over the 8 same-name bench cases that both benchmarks ran at
λ = 0.1 (s502 excluded), the mean composite is old 0.846 → new 0.891. At λ = 0 it is 0.732 → 0.829.

## 2. Per case, old vs new

Mean composite (rank 1), species found / truth, merged flag (M). `—`: not run.

| old → new | λ = 0 old | λ = 0 new | λ = 0.1 old | λ = 0.1 new |
|---|---|---|---|---|
| demo3–demo8 (8 cases) | — | identical | — | identical (demo7_totalonly 0.422 M at λ = 0, 0.976 at 0.1) |
| nest_samehap_2to1_s201 | 0.994 2/2 | 0.994 2/2 | 0.994 2/2 | 0.994 2/2 |
| nest_samehap_2to1_s202 | 0.495 2/2 M | 0.495 2/2 M | 0.989 2/2 | 0.989 2/2 |
| nest_samehap_3to1_s301 | 0.693 3/2 M | 0.693 3/2 M | 0.693 3/2 M | 0.693 3/2 M |
| nest_samehap_3to1_s304 | 0.978 2/2 | 0.978 2/2 | 0.978 2/2 | 0.978 2/2 |
| nest_cross_2to1_s103 | 0.601 3/2 M | **0.972 2/2** | 0.979 2/2 | 0.972 2/2 |
| nest_cross_3sp_s402 | 0.596 3/3 M | 0.522 2/3 M | 0.632 3/3 M | **0.523 2/3 M** |
| nest_cross_1to1_s502 | 0.626 3/2 | 0.955 3/2 M | — | 0.955 3/2 M |
| ctrl_separate_s601 | 0.759 2/2 | **0.982 2/2** | 0.759 2/2 | **0.982 2/2** |
| ctrl_dup_s702 | 0.744 1/1 | **0.997 1/1** | 0.744 1/1 | **0.997 1/1** |
| nest_cross_2to1_s101 → s104 | 0.526 2/2 M | 0.246 1/2 M | 0.964 2/2 | 0.962 2/2 |
| nest_cross_3sp_s406 → s403 | — | 0.578 3/3 M | — | 0.592 3/3 M |
| nest_cross_1to1_s501 → s505 | 0.716 2/2 | 0.765 3/2 M | 0.729 2/2 | 0.765 3/2 M |
| ctrl_separate_s602 → s603 | 0.581 3/1¹ | 0.976 2/2 | 0.581 3/1¹ | 0.976 2/2 |
| ctrl_dup_s701 → s703 | 0.997 1/1 | 0.712 2/1 | 0.997 1/1 | 0.712 2/1 |

¹ The old s602 subgraph (from before `--max-extend`) missed one species entirely, so only 1 truth species was scored.

Per-species composites at λ = 0.1 (new, failures only; all values are in `comparison.tsv`):

- s402: ecdna1 0.973, ecdna2 0.000, ecdna3 0.596;
- s403: 0.635 / 0.483 / 0.657;
- s502: 0.967 / 0.943;
- s505: 0.944 / 0.586;
- s301: 0.964 / 0.423;
- s703: 0.712.

**Failure stage and reason** (`pipeline.tsv`, which is reported at λ = 0.5 in both benchmarks):

| case | old | new |
|---|---|---|
| nest_samehap_3to1_s301 | RCK: ec2 carries ecdna1+ecdna2; 4 species for 2 | same |
| nest_cross_3sp_s402 | RCK: ec3 carries ecdna1+ecdna3 | ok (0.978, 3/3); at λ = 0.1 merged, 2/3 |
| nest_cross_1to1_s501 → s505 | RCK: composites 0.893 / 0.564 | RCK: ec3 carries ecdna1+ecdna2; 3 for 2 |
| nest_cross_1to1_s502 | not finished | RCK: ec3 carries ecdna1+ecdna2; 3 for 2 |
| nest_cross_3sp_s406 → s403 | not finished | RCK: ec2 carries ecdna2+ecdna3 (but see §4: scale) |
| ctrl_separate_s601 | subgraph: min segment coverage 0% | ok |
| ctrl_separate_s602 → s603 | subgraph: min segment coverage 0% | ok |
| ctrl_dup_s702 | RCK: 1 block wrong after flow; imputed edge `imputed_adj_intra_chr16_6` | ok |
| ctrl_dup_s701 → s703 | ok | RCK: 2 species for 1 |

Subgraph stage, new: 0 of 14 subgraphs miss a junction or cover < 80% of a segment (old: 2 of 12).
RCK wall time: 4349 s over 22 cases, against 754 s over 20 before. 65% of it is s403 (965 s at λ = 0 plus
1895 s for the pool) and 23% is s505 (447 + 565 s).

## 3. `ctrl_dup_s702` and `nest_cross_3sp_s402`

**ctrl_dup_s702: now succeeds.** Truth: 1 species, φ 27.57, 13 steps (fold-back).

| | old | new |
|---|---|---|
| junction labels (verify_subgraph) | 5 ok, 1 wrong | 4 ok, 0 wrong |
| imputed edges used | `imputed_adj_intra_chr16_6` | none |
| ecDNA blocks after flow (ok / wrong) | 4 / 1 (2 flipped) | 5 / 0 (3 flipped) |
| flow objective | 15.98 | 3.26 |
| rank 1, every λ | 1 species, 24.91, 26 steps, 0.744 | 1 species, 27.21, 12 steps, **0.997** |

Nothing is left for this case. The re-scored old run (new `evaluate.py`, old labels) was 0.781. The rest of
the gain comes from the simulation fix.

**nest_cross_3sp_s402: fixed upstream; fails at λ ≤ 0.1, succeeds at λ ≥ 0.2.**

| | old | new |
|---|---|---|
| junction labels | 16 ok, 1 wrong (verdict CHECK) | 16 ok, 0 wrong |
| imputed edges | none | none |
| ecDNA blocks after flow (ok / wrong) | 5 / 0 | 5 / 0 (1 flipped) |
| flow objective | 16.49 | 11.81 |
| λ = 0 | 3 sp, 171 steps, M (ec2: e1+e2; ec3: e1+e3), 0.596 | 2 sp (11.97, 4.29), 82 steps, M (ec2: e2+e3), 0.522 |
| λ = 0.02–0.1 | 3 sp, 160 steps, M (ec3: e1+e3), 0.632 | 2 sp (11.75, 4.44), 82 steps, M, 0.523 |
| λ = 0.2–0.5 | same as 0.1 | **3 sp (12.12, 7.46, 4.67), 72 steps, 0.978** |

Junction labels, block flips and imputed edges are all clean now. What remains is in the decomposition
ranking. The two inner species, ecdna2 (mixed, 8.16 bulk copies) and ecdna3 (B, 5.44), are merged into
one walk. That saves one species penalty (1.0) at the cost of 10 more steps. So it is preferred until
λ × 10 steps exceeds the residual difference. In the λ = 0.1 log, the closest alternative (3 species,
86 steps) sits +0.398 behind. (Follow-up: re-ranking on a 0.01 grid puts the switch between 0.10 and 0.11;
at λ = 0.1 the two answers tie on structure cost, 2 + 8.2 = 3 + 7.2.)

Verifier note: in the new `verify_report.txt` for s402, 12 of 33 segment allele calls are "wrong". These
are ecdna3 (B) segments that lie inside ecdna1 (A) and are called AB or A. The verifier's expected
per-species label ignores overlap between nested species on different haplotypes. This is a reporting
artefact, not a pipeline error. The same likely explains the other `seg_calls_wrong` counts: verdict
CHECK on 6 of 14 new subgraphs (s103, s104, s402, s403, s502, s505), all from segment calls.

## 4. Attribution of the 6 old failures, and the 2 unfinished cases

| old failure | attribution | evidence | status now |
|---|---|---|---|
| ctrl_dup_s702 | **simulator artifact** | AA+BB junction, 1 wrong label, imputed edge, 1 block wrong | fixed (0.997) |
| nest_cross_3sp_s402 | **simulator artifact + RCK** | 1 wrong label, AA+BA junction | labels fixed; RCK merge remains at λ ≤ 0.1, fixed at λ ≥ 0.2 |
| ctrl_separate_s601 | **subgraph** (coverage, fixed by `--max-extend`) | min coverage 0% → 100% | fixed (0.982) |
| ctrl_separate_s602 | **subgraph** (coverage, fixed by `--max-extend`) | rebuilt committed subgraph covers 100% | seed rejected in v2; its replacement s603 passes (0.976). The rebuilt s602 subgraph was never run through RCK |
| nest_cross_1to1_s501 | **subgraph (scale)**, newly attributed; labelled RCK before | depth per copy 7.74 (ploidy fit 3.88; others 13.4–14.3); segment CN 1.81× truth in the old flow | seed rejected in v2; replacement s505 has a correct scale and fails in RCK |
| nest_samehap_3to1_s301 | **RCK** | identical simulation, subgraph and result | unchanged (0.693, merged) |

Totals: 2 involve simulator artifacts (1 fully, 1 partly), 3 are subgraph stage (2 coverage, 1 scale),
and 2 remain RCK failures (s301, and s402 at λ ≤ 0.1). s402 is counted in both its rows.

`--max-extend`: the new subgraphs were built with the default `--max-extend` 5 Mb. `ctrl_separate_s601`'s
new subgraph has the same nodes as the committed rebuilt one (70 nodes, same coordinates, 39.31 Mb), and
covers every segment at 100%.

Cases that never finished before:

- **nest_cross_1to1_s502** now finishes: 13 s at λ = 0, 32 s for the pool, against 488 s at λ = 0 and an
  unfinished pool before. Its old subgraph had the same scale error (depth per copy 7.45, segment CN
  1.83× truth). The relabelled simulation fits ploidy 2 (14.13). The result is 0.955 but merged: a third
  species of weight 1.48 carries junctions of both.
- **nest_cross_3sp_s406** is not in the new benchmark. Its old subgraph also had the scale error (depth per
  copy 7.92, ploidy 3.97). Its replacement **s403 has the same problem**: ploidy fit 5.81, depth per copy
  6.51 against about 13.8, segment CN 1.99× truth, flow objective 41.2 (others 3–12), and 181 steps. It
  costs 2860 s in total. Its failure is subgraph stage (scale), not RCK. The `subgraph_construct.log`
  ploidy scores were 0.122 (ploidy 2, purity 1.00) against 0.114 (ploidy 6, purity 0.57). Truth purity
  is 0.80.

## 5. λ

| λ | 0 | 0.02 | 0.05 | 0.1 | 0.2 | 0.5 |
|---|---|---|---|---|---|---|
| mean composite, old (20 cases) | 0.8002 | 0.8026 | 0.8958 | 0.8958 | 0.8958 | 0.8981 |
| mean composite, new (22 cases) | 0.8269 | 0.8600 | 0.9076 | 0.9076 | **0.9283** | 0.9304 |

- **Plateau: no longer holds.** Old: every case except `nest_samehap_3to1_s301` had the same rank-1 answer
  from 0.05 to 0.5. New: `nest_cross_3sp_s402` changes between 0.1 and 0.2 (0.523 → 0.978), and s301
  changes at 0.5 as before (0.693 → 0.740, still merged, 4 species: not a real fix). Every other new case
  is constant from 0.05 to 0.5. The 0.02 → 0.05 step fixes `nest_samehap_2to1_s202` and
  `demo7_totalonly`. `nest_cross_2to1_s104` is already fixed at 0.02.
- **Controls and demos below λ = 0:** none, at any λ, in either benchmark. The largest drop is demo3,
  0.9596 → 0.9587 (−0.0009). The control results do not depend on λ at all (s702, s703, s601 and s603
  are constant).
- **Recommendation:** 0.1 is no longer on the plateau. It misses the s402 fix, and at 0.1 s402 scores below its old result
  (0.523 against 0.632). **λ = 0.2 is the smallest grid value that captures every real gain**,
  with no control or demo cost. collect.py's rule picks 0.5 again, but only because of s301's
  0.693 → 0.740 at the edge of the grid. The case for 0.2 over 0.1 rests on one case (s402), with an
  switch point between 0.10 and 0.11 (follow-up re-rank). That is thin, so confirm it with the re-rank below before
  changing the recommendation.

## 6. Remaining failures (new benchmark, λ = 0.1) and cheapest next diagnostic

None of these needs a new simulation or RCK run. Items marked "re-rank" use the existing pools.

1. **nest_cross_3sp_s403: subgraph scale (ploidy mis-fit).** Segment CN 1.99× truth against
   junction-read CN. Next: `subgraph_construct.py … --rck-ploidy 2` for this one case (about 1 min, no
   Gurobi), then compare `rck.scnt.tsv` totals with truth. Longer term, `pipeline.tsv`'s stage rule should
   flag `depth_per_copy` far from the other cases, or segment CN / truth ≠ 1. Re-check the
   `rck_inputs` ploidy scoring too: it chose ploidy 6 by 0.008.
2. **nest_cross_3sp_s402: λ-dependent merge of ecdna2 + ecdna3.** Next: re-rank
   `python -m rck.ecdna.rerank $BENCH_NESTED_OUT/nest_cross_3sp_s402/pool` at λ = 0.12, 0.15 and 0.18 to
   locate the switch. It runs in seconds and needs no Gurobi solve.
3. **nest_cross_1to1_s502 / s505: a low-weight third species (1.48 / 2.84 copies) carries junctions of
   both species**, so the case is merged with 3 species for 2. The s502 composites are high
   (0.967 / 0.943); s505's ecdna2 is 0.586. The 1.22:1 copy ratio is non-integer, so integer flow
   leaves a remainder for a third walk to absorb. Next: the λ = 0 runs report a near-tie decided by
   species_penalty (s502 gap 0.004, s505 0.002). Read the 2-species alternatives in
   `pool/decompositions.tsv` and their scores in `pool/evaluation.tsv`.
4. **nest_samehap_3to1_s301: unchanged RCK merge** (ec2 carries ecdna1+ecdna2; 3 for 2). This is the
   old open case. Next: compare the truth 2-species decomposition's objective with rank 1 (the truth's
   integer CNs against the flow). Is it in the pool, or does it fit the flow worse?
5. **ctrl_dup_s703: 2 species for 1** (ec1 17.93 with both fold-backs; ec2 16.22, the 3-step piece
   through the B insertion chr1:182.61–183.29 Mb). The flow puts 16 copies on the crossing junctions
   sim00002 and sim00003, against 18–19 on the fold-backs. The B segment's observed CN is 16.1 against
   a truth of 19.2, while the shared segments match the truth. So one weight cannot fit all of them.
   Next: look at that segment's window depth in `obs/windows` against the truth (a simulator depth
   artifact, or a real noise outlier?). No RCK change is suggested yet.

## Caveats

- The old `pipeline.tsv` stage rule cannot see scale errors. s501 and s406 were labelled RCK or
  unfinished, but their subgraph inputs carried twice the copy numbers.
- The new README in `rck_results/` says "re-ranking check … 32/32". `rerank_validation.tsv` is
  byte-identical to the committed one (it comes from the linked demo runs). It was not re-validated.
- Gurobi WLS refused sessions ("5 active sessions for a baseline of 2") on several tasks of job 32000320.
  Every task succeeded on a retry, and every case has results.

## Follow-up checks (2026-10-27)

See RCK `results/bench_nested/README.md`, "Follow-up checks". Diagnostic runs and scripts are in
`/n/fs/scratch/ym3292/bench_nested_v2/diag/`. Corrections to the text above: the s402 switch is at λ 0.10 → 0.11, not
0.14–0.2. s403 with the correct scale (`--rck-ploidy 2`) still fails (0.545), so its failure is
subgraph scale and then RCK.

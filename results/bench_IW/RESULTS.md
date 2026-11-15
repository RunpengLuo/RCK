# RCK-ecDNA vs CoRAL and CycleExtractor: manuscript figures and Table 1

*Written 2026-11-15. Benchmark: bench_IW_intpsi (arrays 32280369–32280372, all 96 runs finished).
Run-level details, failures and the open items are in `RESULTS_intpsi.md` and `NEXT_STEPS.md`.*

## How to regenerate

```bash
cd ~/RCK-cnp/results/bench_IW
export PATH=/n/fs/yubi-iw/conda/envs/gurobi-env/bin:$PATH
python make_figures.py                                    # CoRAL / CE scored with linear walks kept
python make_figures.py --baseline-variant walks_dropped --out figures_intpsi_walks_dropped
```

- `make_figures.py` reads only `bench_iw_intpsi_table.tsv` and `bench_iw_intpsi_proportions.tsv`
  (the output of `bench_iw.py --set iw_intpsi collect`). It takes a few seconds and uses no Gurobi.
- Output goes to `figures_intpsi/`: each figure as PDF (vector, fonts embedded) and PNG (300 dpi),
  `table1.md` / `table1.tsv`, `table1_paired_tests.tsv`, and `figure_data/*.tsv` with the exact
  numbers each figure plots (useful for source-data files).
- `figures_intpsi_walks_dropped/` is the same set with CoRAL / CE scored with linear walks
  dropped. Use it as the supplementary robustness check; Table 1 lists both variants either way.

## The benchmark in one paragraph (for the Methods)

- **Design:** 24 simulated tumours: 2 classes (`simple` and `nested`, where one ecDNA species
  contains another's segments) x S in {1, 2, 4} ecDNA species x 4 replicates.
- **Clones:** all cases share one set of clones. Each clone carries a whole number of copies (psi)
  of each species, so the scoring truth, PHI_s = sum_k u_k psi_k, is exactly what the multi-clone
  input encodes.
- **Methods:**
  - **RCK-ecDNA bulk** (`rck_bulk`): builds and fits its own allele-specific graph, then
    decomposes it.
  - **RCK-ecDNA multi-clone** (`rck_mc`): given the true clone profiles.
  - **CoRAL and CycleExtractor (CE):** run as *decomposition-only* baselines on RCK bulk's
    exported graph (fitted copy numbers) plus the long-read path constraints.
- **Scoring** is by `rck.ecdna.aa_eval` / `bench_metrics` against the simulator's truth
  (`docs/pipeline/evaluation.md`).

## Conventions used in every figure

- Colour = method: RCK bulk blue, RCK multi-clone orange, CoRAL green, CE yellow.
- Marker = class: circle simple, triangle nested.
- One marker = one case (24 per method), unless stated.
- \* on the multi-clone route: given the true clone profiles. Present it as an upper bound for that
  route, not as a like-for-like competitor.
- † on CoRAL / CE: decomposition only, on RCK bulk's graph.
- Main figures show CoRAL / CE with **linear walks kept**, which is the tools' unmodified output.
- The four colours are categorical slots 1–4 of the reference palette. Adjacent pairs pass the
  colour-vision-deficiency check, but four colours cannot all be told apart in one scatter. So no
  panel relies on colour alone: methods are also separated by panel or x position and named in
  titles or tick labels.

---

## Figure 2: ecDNA copy-number recall vs precision

`fig2_cn_recall_precision.pdf`

**What it shows.** One panel per method; each point is a case.
- *x:* the share of the true ecDNA copy-number mass (bp x copies, haplotype-blind) that the
  method's strict ecDNA cycles recover.
- *y:* the share of the method's ecDNA mass that is true.
- *+:* the mean.

**Reading.**
- Both RCK routes sit in the top-right corner: bulk mean 0.96 / 0.98, multi-clone 0.99 / 1.00.
- CoRAL (0.85 / 0.48) and CE (0.96 / 0.40) recover most of the true ecDNA but report a lot that
  is not ecDNA. Their precision spreads from about 0 to 1 across cases.
- This is the headline figure: the same graph, decomposed by RCK, gives ecDNA copy numbers that are
  both complete and clean.

**Limitations.**
- Copy-number recall and precision are haplotype-blind and say nothing about cycle structure;
  Figure 3 covers that.
- The one RCK multi-clone outlier (recall 0.71, hard_S4_s16003) is unexplained (`NEXT_STEPS.md`,
  "After step 7").
- CE's recall equals RCK bulk's. Only precision separates them (paired tests: recall n.s.,
  precision p < 0.001).

## Figure 3: structure accuracy as the number of species grows

`fig3_structure_vs_S.pdf`

**What it shows.**
- *Columns:* simple and nested cases. *Rows:* best-match segment Jaccard (J_seg) and adjacency
  Jaccard (J_adj).
- Both Jaccards are haplotype-blind, per true species, averaged within a case.
- *Small points:* cases. *Large points and lines:* means per S.

**Reading.**
- On simple cases RCK stays near 1 at every S; CE declines with S (0.97 to 0.74) and CoRAL is
  lower from S = 1 (0.73).
- On nested cases the methods separate. At S = 4, RCK multi-clone keeps J_seg 0.94 and J_adj
  0.81; RCK bulk drops to 0.61 / 0.56, CE to 0.45 / 0.29 and CoRAL to 0.31 / 0.28.
- This is the "harder instances" figure: nested ecDNA with several species is where
  decomposition is ambiguous, and where the allele-specific, clone-aware model pays off.

**Limitations.**
- Only 4 cases per point. The means are noisy (see the spread of the small points), and S takes
  three values only.
- J_adj for RCK multi-clone on simple S = 2 (0.76) is below the other methods. That comes from
  the extra cycles it sometimes adds (Figure 6; `RESULTS_intpsi.md`).
- 3 CoRAL runs with linear walks dropped return no cycle at all; they count as 0 here. The pooled
  TSV from `collect` skips them, which overstates CoRAL (0.75 instead of 0.65).

## Figure 4: how much non-ecDNA each method reports as ecDNA

`fig4_mass_ratio.pdf`

**What it shows.** Inferred ecDNA copy-number mass divided by true ecDNA mass, per case, on a log
scale. The bar is the median, also given under each method; the dashed line marks 1 (exact).

**Reading.**
- Both RCK routes sit at about 1 in every case: median 0.99 and 1.00, none above 1.5.
- CoRAL (median 3.30) and CE (3.25) report more than 1.5x the true ecDNA mass in 15 and 19 of 24
  cases, and over 10x in 4 cases each.
- In those cases they return the chromosomal background as circular ecDNA. This explains the low
  precision in Figure 2, and it is a failure a reader can picture.

**Limitations.**
- With linear walks dropped, CoRAL's median falls to 0.77 (2 cases above 1.5x) but its recall
  falls too: Table 1 and `figures_intpsi_walks_dropped/`. Show both, or say which you chose and
  why.
- CE is unchanged either way, because it returns the background as closed cycles.
- The baselines get RCK bulk's *total* copy-number graph, with no ecDNA vs background separation.
  Part of what the figure shows is that RCK's model separates the two and theirs does not.
  Reviewers may argue that this is not CoRAL's intended input; the end-to-end comparison in
  `NEXT_STEPS.md` addresses that.

## Figure 5: species proportions

`fig5_proportions.pdf`

**What it shows.**
- *Top:* each true species' share of the case's ecDNA copies against the share the method
  assigns its matched cycle, pooled over all cases (strict cycles).
  - × on the y axis: a false species (truth 0).
  - × on the x axis: a missed species (inferred 0).
  - Titles give the totals.
- *Bottom:* total-variation distance between the true and inferred proportion vectors, per case
  (0 = exact, lower is better); the bar is the mean.

**Reading.**
- RCK points hug the diagonal: mean TV 0.07 (bulk) and 0.06 (multi-clone).
- CoRAL (0.23) and CE (0.33) scatter widely, and CE adds 48 false species over 24 cases.
- RCK multi-clone misses no species.

**Limitations.**
- Proportions depend on the strict-cycle rule, which differs: RCK uses its own species /
  `low_weight` labels, CoRAL / CE a 10% copy-number fraction cut-off (`cn-frac:0.1`). The
  `low_weight`-inclusive version is in the TSVs (`proportion_tv_all`) and tells the same story.
- Species are matched one-to-one by best blind Jaccard. A method that splits one species into
  two is charged a false species even when the copy numbers are right. In the RCK multi-clone
  case checked by hand (simple_S1_s11004) its false species are of this kind.

## Figure 6: allele-specific results (RCK only)

`fig6_allele_specific.pdf`

CoRAL and CE report no alleles, so this figure has no baseline. It shows what RCK-ecDNA adds.

**a: phased vs blind structure.**
- J_seg with alleles required to match against J_seg ignoring alleles, per case.
- 37 of 48 runs sit at (1, 1). Points on the diagonal mean the cycles are not only right but on
  the right haplotype.
- RCK bulk loses phase only in one simple case (blind 1.00, phased 0.75).

**b: phase-block identifiability (multi-clone).**
- Missed + false species per case, scored with absolute A/B labels (open circles) and within
  phase blocks (filled; one A/B swap per block, shared by all species).
- Most of the absolute errors vanish within blocks: 30 vs 8 in all.
- This is how the paper should justify within-block scoring: A and B are identified only within
  a phase block, as in the original RCK's per-fragment allele flip.

**c: exact cycles.**
- Exact cycles per true species, blind (bars) and phased (black ticks), by class.
- On nested cases, multi-clone recovers 0.68 against 0.46 for bulk.
- On simple cases bulk is higher (0.96 vs 0.75), because of multi-clone's extra split cycles.

**Limitations.**
- Within-block scoring chooses the swap that agrees best with the truth, which is favourable by
  construction. Always show the absolute numbers next to it, as panel b does.
- Panel a overplots: most points are at (1, 1); the count is printed in the panel.

## Table 1: per method variant, all cases and per class

`table1.md` (paste-ready), `table1.tsv` (full precision)

**What it shows.**
- Mean ± SD over cases, for all 24 cases and for the 12 simple and 12 nested cases.
- Columns: ecDNA copy-number recall, precision, proportion TV, species-detection F1 (blind,
  strict cycles), J_seg, J_adj, and the exact-cycle rate (exact cycles / true species, summed).
- All six method variants are listed: CoRAL / CE with linear walks kept and dropped. The best
  mean per column is in bold.

**Reading.**
- RCK is best in every column. The gap is largest on nested cases, where multi-clone leads every
  column.
- Neither CoRAL nor CE reaches RCK's precision and recall together.

**`table1_paired_tests.tsv`.** Two-sided Wilcoxon signed-rank tests over the 24 paired cases, for
RCK bulk and multi-clone against each baseline variant, with a Holm correction within each pair of
methods. Columns:
- the median difference;
- the number of cases where RCK is better, and ties;
- the raw and Holm-corrected p-values.

Main results:
- Copy-number precision: p < 0.002 against every walks-kept baseline and against CE either way.
  Against CoRAL with walks dropped it is not significant after correction (Holm p 0.38 for bulk,
  0.05 for multi-clone), because that variant has high precision and low recall.
- Proportion TV: p < 0.03 against every baseline variant.
- CE's copy-number recall is not significantly different from RCK's (CE is higher in 20 of 24
  cases for bulk, by a median of 0.004).

**Limitations.**
- n = 24 (12 per class) is small for per-class tests; report per-class differences descriptively.
- Species F1 is lower for multi-clone (0.70) than for bulk (0.76) because of its extra split
  cycles. Don't hide it.
- The strict rule differs between RCK and the baselines (see Figure 5). The rule-free columns
  (copy-number recall, precision, Jaccards) are the fairest comparison.
- Copy-number precision is undefined for the 3 CoRAL walks-dropped runs with no cycle. Those
  pairs are left out of that test (n = 21).

---

## Limitations that apply to the whole comparison

Say these in the paper; reviewers will raise them otherwise. How to address each is in
`NEXT_STEPS.md`, "Additions reviewers will likely ask for".

1. **The baselines only decompose.** CoRAL and CE run on RCK bulk's graph, so this benchmark
   compares decomposition given the same graph. It does not compare full pipelines from reads.
   Frame the claim as "given the same graph, RCK's decomposition is more accurate", or add the
   end-to-end comparison.
2. **The multi-clone route is given the true clone profiles.** Its numbers are an upper bound for
   that route. RCK bulk is the like-for-like comparison with CoRAL / CE.
3. **Simulated data only**, from one simulator (ecSimulator-based) with one read model. There is
   no sweep over depth, purity, read length or breakpoint noise.
4. **Small benchmark:** 24 cases, 4 per (class, S), S <= 4, one shared set of clones.
5. **Different strict rules** for RCK and the baselines; the rule-free metrics are the fair ones.
6. **Scoring choices that favour RCK multi-clone,** both disclosed and shown both ways:
   - within-phase-block counts use the best-agreeing swap;
   - linear walks are scored both kept and dropped for the baselines.
7. **Timing is not comparable as recorded:** `elapsed_s` includes Gurobi licence waits. Use the
   solver times in each run's logs / `summary.json` if runtime is reported.
8. **Two solves stopped at their time limit** and are scored on the incumbent: RCK bulk on
   hard_S4_s16003 (gap 79%) and CoRAL on simple_S4_s13002 (gap 9%). There were no other failures
   in this benchmark. CE crashed on 2 of 24 graphs on the earlier bench_IW and is otherwise not
   tested here.

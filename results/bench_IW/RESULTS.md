# RCK MC vs CoRAL and CE: manuscript figures and Table 1

*Written 2026-11-15; rewritten 2026-11-21 for the three-method comparison (Yubi: compare only RCK MC, CoRAL
and CE, labelled so, with no footnote marks on the figures). Benchmark: bench_IW_intpsi (arrays
32280369–32280372, all 96 runs finished). Run-level details, failures and the open items are in
`RESULTS_intpsi.md` and `NEXT_STEPS.md`.*

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
- **Methods compared:**
  - **RCK MC** (`rck_mc`, RCK-ecDNA's multi-clone route): decomposes the allele-specific multi-clone
    graph, **given the true clone copy-number profiles and proportions** (as RCK's own simulations
    used the true copy numbers).
  - **CoRAL and CE (CycleExtractor):** run as *decomposition-only* baselines. Their input is a breakpoint
    graph with fitted total copy numbers, built by RCK-ecDNA's bulk route (`rck_bulk`) and exported, plus
    the long-read path constraints. `rck_bulk` is used only to build this graph; its own answers are not
    shown or scored here.
- **Scoring** is by `rck.ecdna.aa_eval` / `bench_metrics` against the simulator's truth
  (`docs/pipeline/evaluation.md`).
- **Say in the text that the inputs differ** (Limitations, item 2): RCK MC is given the true clone
  profiles, and CoRAL / CE get a total copy-number graph. The figures carry no marks for this.

## Conventions used in every figure

- Labels: RCK MC, CoRAL, CE. No footnote marks.
- Colour = method: RCK MC orange, CoRAL green, CE yellow (categorical slots 2–4 of the reference
  palette, in their validated order; the same colours as before RCK bulk was dropped).
- Marker = class: circle simple, triangle nested.
- One marker = one case (24 per method), unless stated.
- Main figures show CoRAL / CE with **linear walks kept**, which is the tools' unmodified output.
- No panel relies on colour alone: methods are also separated by panel or x position and named in
  titles, legends or tick labels.

---

## Figure 2: ecDNA copy-number recall vs precision

`fig2_cn_recall_precision.pdf`

**What it shows.** One panel per method; each point is a case.
- *x:* the share of the true ecDNA copy-number mass (bp x copies, haplotype-blind) that the
  method's strict ecDNA cycles recover.
- *y:* the share of the method's ecDNA mass that is true.
- *+:* the mean.

**Reading.**
- RCK MC sits in the top-right corner: mean recall 0.99, precision 1.00.
- CoRAL (0.85 / 0.48) and CE (0.96 / 0.40) recover most of the true ecDNA but report a lot that
  is not ecDNA. Their precision spreads from about 0 to 1 across cases.
- The headline figure: RCK MC's ecDNA copy numbers are both complete and clean.

**Limitations.**
- Copy-number recall and precision are haplotype-blind and say nothing about cycle structure;
  Figure 3 covers that.
- RCK MC has one outlier (recall 0.71, hard_S4_s16003), unexplained (`NEXT_STEPS.md`, "After step 7").
- CE's recall is not different from RCK MC's (paired test p 0.53; CE is higher in 15 of 24 cases, by tiny
  margins). Only precision separates them (Holm p < 0.001).

## Figure 3: structure accuracy as the number of species grows

`fig3_structure_vs_S.pdf`

**What it shows.**
- *Columns:* simple and nested cases. *Rows:* best-match segment Jaccard (J_seg) and adjacency
  Jaccard (J_adj).
- Both Jaccards are haplotype-blind, per true species, averaged within a case.
- *Small points:* cases. *Large points and lines:* means per S.

**Reading.**
- On simple cases RCK MC's J_seg stays at 0.96–1.00 for every S; CE declines with S (0.97 to 0.74) and
  CoRAL is lower from S = 1 (0.73).
- On nested cases the methods separate. At S = 4, RCK MC keeps J_seg 0.94 and J_adj 0.81; CE drops to
  0.45 / 0.29 and CoRAL to 0.31 / 0.28.
- This is the "harder instances" figure: nested ecDNA with several species is where decomposition is
  ambiguous, and where the allele-specific, clone-aware model pays off.

**Limitations.**
- Only 4 cases per point. The means are noisy (see the spread of the small points), and S takes
  three values only.
- RCK MC's J_adj is below the baselines' on simple S = 2 (0.76 vs 1.00) and S = 4 (0.85 vs 0.88–0.94).
  That comes from the extra cycles it sometimes adds (Figure 6; `RESULTS_intpsi.md`). Over all cases
  RCK MC's J_adj is not significantly different from either baseline with walks kept (Table 1 tests).
- 3 CoRAL runs with linear walks dropped return no cycle at all; they count as 0 here. The pooled
  TSV from `collect` skips them, which overstates CoRAL (0.75 instead of 0.65).

## Figure 4: how much non-ecDNA each method reports as ecDNA

`fig4_mass_ratio.pdf`

**What it shows.** Inferred ecDNA copy-number mass divided by true ecDNA mass, per case, on a log
scale. The bar is the median, also given under each method; the dashed line marks 1 (exact).

**Reading.**
- RCK MC sits at about 1 in every case: median 1.00, none above 1.5.
- CoRAL (median 3.30) and CE (3.25) report more than 1.5x the true ecDNA mass in 15 and 19 of 24
  cases, and over 10x in 4 cases each.
- In those cases they return the chromosomal background as circular ecDNA. This explains their low
  precision in Figure 2, and it is a failure a reader can picture.

**Limitations.**
- With linear walks dropped, CoRAL's median falls to 0.77 (2 cases above 1.5x), but its recall
  falls too (0.69): Table 1 and `figures_intpsi_walks_dropped/`. Show both, or say which you chose and
  why.
- CE is nearly unchanged either way (median 3.22), because it returns the background as closed cycles.
- CoRAL / CE get a *total* copy-number graph, with no ecDNA vs background separation; RCK MC's input
  separates the two. Part of what the figure shows is this difference in input and model. Reviewers may
  argue that this is not CoRAL's intended input; the end-to-end comparison in `NEXT_STEPS.md` addresses
  that.

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
- RCK MC's points hug the diagonal: mean TV 0.06, 0 missed species, 8 false over 24 cases.
- CoRAL (TV 0.23; 11 missed, 7 false) and CE (0.33; 6 missed, 48 false) scatter widely.

**Limitations.**
- Proportions depend on the strict-cycle rule, which differs: RCK MC uses its own species /
  `low_weight` labels, CoRAL / CE a 10% copy-number fraction cut-off (`cn-frac:0.1`). The
  `low_weight`-inclusive version is in the TSVs (`proportion_tv_all`) and tells the same story.
- Species are matched one-to-one by best blind Jaccard. A method that splits one species into
  two is charged a false species even when the copy numbers are right. In the RCK MC case checked by
  hand (simple_S1_s11004) its false species are of this kind.
- With linear walks dropped, CoRAL misses 21 species (1 false) and CE 7 (21 false).

## Figure 6: allele-specific results (RCK MC only)

`fig6_allele_specific.pdf`

CoRAL and CE report no alleles, so this figure has no baseline (the caption says so). It shows what the
allele-specific model adds.

**a: phased vs blind structure.**
- J_seg with alleles required to match against J_seg ignoring alleles, per case.
- 20 of 24 runs sit at (1, 1). Points on the diagonal mean the cycles are not only right but on the right
  haplotype.
- Phasing costs almost nothing: only 2 cases lose any J_seg to phase (hard_S2_s15003 0.82 → 0.80,
  hard_S4_s16003 0.76 → 0.74).

**b: phase-block identifiability.**
- Missed + false species per case, scored with absolute A/B labels (open circles) and within
  phase blocks (filled; one A/B swap per block, shared by all species).
- Most of the absolute errors vanish within blocks: 30 vs 8 in all.
- This is how the paper should justify within-block scoring: A and B are identified only within
  a phase block, as in the original RCK's per-fragment allele flip.

**c: exact cycles.**
- Exact cycles per true species, blind (bars) and phased (black ticks), by class: 0.75 simple, 0.68
  nested, the same blind and phased.
- For comparison (Table 1, blind): CoRAL 0.68 / 0.18 and CE 0.64 / 0.25 (simple / nested).

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
- Five method variants: RCK MC, and CoRAL / CE with linear walks kept and dropped. The best mean per
  column is in bold.

**Reading.**
- Over all cases RCK MC has the best mean in every column (recall 0.99, precision 1.00, TV 0.06, F1 0.70,
  J_seg 0.97, J_adj 0.87, exact-cycle rate 0.71 against at most 0.45).
- The gap is largest on nested cases: exact-cycle rate 0.68 against 0.18 (CoRAL) and 0.25 (CE); J_adj
  0.87 against 0.52–0.65.
- On simple cases the baselines' J_adj (0.96, walks kept) is higher than RCK MC's (0.87); every other
  simple-case column favours RCK MC.
- Neither CoRAL nor CE reaches RCK MC's precision and recall together.

**`table1_paired_tests.tsv`.** Two-sided Wilcoxon signed-rank tests over the 24 paired cases, RCK MC
against each baseline variant, with a Holm correction over the 7 metrics within each pair. Columns:
- the median difference;
- the number of cases where RCK MC is better, and ties;
- the raw and Holm-corrected p-values.

Main results (Holm-corrected):
- Copy-number precision: p < 0.001 against CoRAL with walks kept and against CE either way. Against
  CoRAL with walks dropped p = 0.05 (n = 21; 14 ties), because that variant has high precision and low
  recall.
- Copy-number recall: p < 0.004 against CoRAL either way; not significant against CE (walks kept p 0.53,
  dropped p 0.054).
- Proportion TV: p < 0.025 against every baseline variant.
- J_seg: p 0.03–0.054 against every variant. J_adj: not significant against any walks-kept variant
  (p 0.18–0.53); p 0.04 against CoRAL with walks dropped.
- Exact cycles and species F1: p 0.008–0.07, around the 0.05 line; report as a trend.

**Limitations.**
- n = 24 (12 per class) is small for per-class tests; report per-class differences descriptively.
- RCK MC's species F1 (0.70) is limited by its extra split cycles. Don't hide it.
- The strict rule differs between RCK MC and the baselines (see Figure 5). The rule-free columns
  (copy-number recall, precision, Jaccards) are the fairest comparison.
- Copy-number precision is undefined for the 3 CoRAL walks-dropped runs with no cycle. Those
  pairs are left out of that test (n = 21).

---

## Limitations that apply to the whole comparison

Say these in the paper; reviewers will raise them otherwise. How to address each is in
`NEXT_STEPS.md`, "Additions reviewers will likely ask for".

1. **The baselines only decompose.** CoRAL and CE run on a graph built by RCK-ecDNA's bulk route,
   so this benchmark compares decomposition. It does not compare full pipelines from reads.
2. **The inputs are not the same, and no figure shows a like-for-like comparison.** RCK MC is given the
   true clone copy-number profiles and proportions; CoRAL / CE get a fitted total copy-number graph. RCK
   MC's numbers are therefore an upper bound for that route. Earlier versions of these figures included
   RCK-ecDNA's bulk route as the like-for-like comparison and marked the difference with * / †; both were
   removed on 2026-11-21 (Yubi), so the text must state it. Frame the claim accordingly, or add the
   end-to-end comparison.
3. **Simulated data only**, from one simulator (ecSimulator-based) with one read model. There is
   no sweep over depth, purity, read length or breakpoint noise.
4. **Small benchmark:** 24 cases, 4 per (class, S), S <= 4, one shared set of clones.
5. **Different strict rules** for RCK MC and the baselines; the rule-free metrics are the fair ones.
6. **Scoring choices that favour RCK MC,** both disclosed and shown both ways:
   - within-phase-block counts use the best-agreeing swap;
   - linear walks are scored both kept and dropped for the baselines.
7. **Timing is not comparable as recorded:** `elapsed_s` includes Gurobi licence waits. Use the
   solver times in each run's logs / `summary.json` if runtime is reported.
8. **One shown solve stopped at its time limit** and is scored on the incumbent: CoRAL on
   simple_S4_s13002 (gap 9%). (The bulk-route run that built the baselines' graph for hard_S4_s16003 also
   stopped at its limit, gap 79%, so CoRAL / CE decompose that incumbent's graph there.) There were no
   other failures in this benchmark. CE crashed on 2 of 24 graphs on the earlier bench_IW and is otherwise
   not tested here.

# Evaluation: scorers, metrics and AA-format translators

This document defines every number the pipeline uses to judge a reconstruction against the
simulator's truth, and the translators that put every method's output into one format. It
covers stage 4 of the pipeline (see [the overview](README.md)):

1. `rck-ecdna-eval` (`rck/ecdna/evaluate.py`): `rck-ecdna`'s own scorer. It produces
   `evaluation.tsv`, `evaluation.json` and the recovery matrix with its composite score.
2. The method-agnostic scorer `python -m rck.ecdna.aa_eval` (`rck/ecdna/aa_eval.py`,
   `rck/ecdna/aa_format.py`). It reads AmpliconArchitect-format (AA-format) files only, so
   AmpliconArchitect (AA) [Deshpande et al. 2019](references.md#deshpande2019), CoRAL [Zhu et al. 2024](references.md#zhu2024), CycleExtractor
   [Faizrahnemoon et al. 2026](references.md#cycleextractor), the original RCK [Aganezov & Raphael 2020](references.md#aganezov2020) and `rck-ecdna` all go through
   the same code.
3. The translators between our subgraphs and the AA graph format: `aa_export` (answer → AA,
   subgraph → AA), `aa_import` (AA → subgraph), and `coral_flank` (the CoRAL flank fix). They are
   the basis of the graph route (`results/graph_route/`).

Every step below carries a provenance line with its category: **(a)** taken from an existing
method, **(b)** a patch that changes no result, or **(c)** new or modified for this pipeline. It
also names the code. The terms follow RCK [Aganezov & Raphael 2020](references.md#aganezov2020):

- *novel adjacency*: the code's "junction", AA's "discordant edge";
- *reference adjacency*: the code's "reference join", AA's "concordant edge";
- *extremity*: a segment end, AA's "breakend" `chr:pos±`;
- *haplotype labelling* (phasing) of an adjacency: AA, AB, BA or BB.

A *species* is one distinct ecDNA molecule type; its *walk* (AA: *cycle*, or *path* when it is
linear) is the ordered list of oriented, haplotype-labelled segments it traverses.

---

## Contents

1. [Common objects and notation](#1-common-objects-and-notation)
2. [Reading the truth](#2-reading-the-truth)
3. [`rck-ecdna-eval`: species metrics, matching, recovery matrix and composite](#3-rck-ecdna-eval-species-metrics-matching-recovery-matrix-and-composite)
4. [The method-agnostic scorer `aa_eval`](#4-the-method-agnostic-scorer-aa_eval)
5. [Translators: `aa_export`, `aa_import`, `coral_flank`](#5-translators-aa_export-aa_import-coral_flank)
6. [Benchmark-level summaries built on these metrics](#6-benchmark-level-summaries-built-on-these-metrics)
7. [Known limitations](#7-known-limitations)
- [Appendix A. How to run](#appendix-a-how-to-run)
- [Appendix B. Output files and columns](#appendix-b-output-files-and-columns)
- [Appendix C. Tests](#appendix-c-tests)

---

## 1. Common objects and notation

> **Provenance:** (c) new, RCK yubi_dev@c144584. **Code:** `rck/ecdna/evaluate.py::CStep`,
> `::Cycle`, `::runs`, `::junctions`.

A step is $s = (c, a, b, h, o)$:

- $c$ is the chromosome;
- $[a, b]$ is a 0-based inclusive interval;
- $h \in \lbrace \mathrm{A}, \mathrm{B}, . \rbrace$ is the haplotype, where `.` means unphased;
- $o \in \lbrace +, - \rbrace$ is the orientation.

Its length is $\ell(s) = b - a + 1$. A walk $C = (s_1, \dots, s_n)$ is circular or linear, has a
weight $w_C$ (copies per cell), and may carry several *orderings*. RCK-ecDNA reports every
distinct Euler ordering of a species' edges, and the first one is $C$ itself. The length of a
walk counts multiplicity: $\lvert C \rvert = \sum_i \ell(s_i)$.

**Multiplicity.** For a base $x$ on chromosome $c$ and haplotype $h$:

```math
m_C(c, x, h) = \#\lbrace i : c_i = c,\ h_i = h,\ a_i \le x \le b_i \rbrace,
\qquad
m_C(c, x) = \sum_{h} m_C(c, x, h).
```

**Reference contiguity.** Consecutive steps $s, t$ are contiguous when they agree in $c$, $h$
and $o$ and the second continues the first within 2 bp:

- $\lvert a_t - (b_s + 1) \rvert \le 2$ for $o = +$;
- $\lvert a_s - (b_t + 1) \rvert \le 2$ for $o = -$.

`runs(C)` merges maximal contiguous stretches into one step, cyclically for a circle. A different
segmentation of the same molecule therefore gives the same runs. Because contiguity requires the
same haplotype, a haplotype switch between reference neighbours is a join, not a continuation.

**Joins (novel adjacencies of a walk).** `junctions(C)` lists every non-contiguous consecutive
pair as $(e^{\mathrm{exit}}, e^{\mathrm{entry}}, h_s, h_t)$. Here $e = (c, \mathrm{pos}, \pm)$ is
an extremity: the exit end is $(c, b_s, +)$ when $o_s = +$ and $(c, a_s, -)$ otherwise. For a
circle the last step joins the first; a linear walk has no such wrap-around join.

**Extremity match.** Two extremities match within a tolerance $\tau$ when they have the same
chromosome, the same strand and $\lvert \mathrm{pos} - \mathrm{pos}' \rvert \le \tau$. Two joins
match when both ends match, in either direction. A *phased* match also requires the haplotype
pair to agree, swapped when the direction is.

Tolerances used below:

| Symbol | Option | Default | Used for |
|---|---|---|---|
| $\tau_s$ | `rck-ecdna --eval-tol-bp`, `rck-ecdna-eval --tol-bp`, `aa_eval --species-tol` | 10 000 bp | species metrics (the subgraph's imputed breakpoints sit on a 10 kb grid) |
| $\tau_g$ | `aa_eval --graph-tol` | 100 bp | graph metrics, `exact_cycle`, subwalks (the AA and CoRAL convention) |

> **Inclusive tolerance.** All matches here use $\le \tau$. CoRAL 3.1.1 uses a strict $< \tau$
> (`coral/scoring/scoring_utils.py::is_oriented_breakpoint_match`) [CoRAL 3.1.1 code](references.md#coral_code).

---

## 2. Reading the truth

### 2.1 Species walks and coordinates

> **Provenance:** (c) new, RCK yubi_dev@c144584, generalised in yubi_dev@3ba7627. Parsing has
> been shared with `aa_format` since yubi_dev@e4c2409. The AA cycles format itself is (a)
> AmpliconArchitect's [Deshpande et al. 2019](references.md#deshpande2019), as ecSimulator writes it [ecSimulator](references.md#ecsimulator).
> **Code:** `rck/ecdna/evaluate.py::truth_from_simulation`, `::read_aa_cycles`,
> `rck/ecdna/aa_format.py::raw_walks`.

For each species listed in `truth/ecdna.tsv` (or else in `truth/aa/`):

1. **Walk and coordinates.** They come from `truth/aa/<species>_cycles.txt`, which the
   simulator writes with ecSimulator's own AA writer. This works whatever the species went
   through: deletions, shared chromosomes, several species. Older runs without `truth/aa/` fall
   back to rebuilding the walk from `truth/walks.walk` and `truth/segments.tsv`
   (`_segments_from_profiles`). There, the species' segments are the rows whose copy number
   exceeds the clone-mixture background, and the interval order is the one under which every
   join of the walk is a reference join or a junction of that species. AA coordinates are
   1-based inclusive; a segment row is read as `START = start − 1`, `END = end`, and a step as
   $[\mathrm{start} - 1, \mathrm{end} - 1]$.
2. **Copy number.** $w_T$ is the species' copies per cell in the evaluated sample: the
   `ecdna_burden` rows of `truth/proportions.tsv`, column `PHI` (`species_weights`). The sample
   is `--truth-sample`, else `--sample` (`sim_tumor1` maps to `tumor1`, `resolve_sample`), else
   the first tumor. Without proportions it falls back to `ecdna.tsv` `PHI`.

The rebuilt truth is written to `truth_cycles.tsv`. For real data, write that table yourself and
pass it as `--truth-cycles`: columns `species order chrom start end strand hap weight`, 0-based
inclusive coordinates (`read_truth_table`).

### 2.2 Truth haplotypes, including junctions with several phasings

> **Provenance:** (c) new, RCK yubi_dev@c144584. The several-phasings rule is (c),
> yubi_dev@184ec95. **Code:** `rck/ecdna/evaluate.py::_label_haplotypes`, `::truth_phasings`,
> `::species_junctions`, `::_majority_hap`.

Haplotypes come from the species' junction phasings (`truth/adjacency_clone.tsv` plus
`obs/adjacency.tsv`) and are carried across reference joins with a union-find over the walk's
steps. A step that no junction labels takes the species' majority haplotype.

`truth/adjacency_clone.tsv` keeps one `PHASING` per junction, its most common. Runs simulated with
per-step haplotype labels (amplicon-allele-analysis `hap_labels: per_step`, every run before that
key) can traverse one junction with several phasings, for example AA on one copy and BB on its
duplicate. `truth/rck/rck.acnt.tsv` lists them all, and `truth_phasings` reads every phasing with
copies in some clone. A junction between reference neighbours may also be traversed as a plain
reference join on one haplotype, so AA and BB are added to its options.

Let $Q$ be the traversals that have more than one option, and $\Phi_q$ the options of traversal
$q$. The chosen assignment is the lexicographic minimum:

```math
\pi^\star = \mathop{\mathrm{arg\,lexmin}}\limits_{\pi \in \prod_{q \in Q} \Phi_q}
\Big( \mathrm{conflicts}(\pi),\ -\sum_{j} \big\lvert \lbrace \pi_q : q \text{ traverses } j,\ \pi_q \in \Phi^{\mathrm{listed}}_j \rbrace \big\rvert \Big)
```

where

- $\mathrm{conflicts}(\pi)$ is the number of steps whose union-find class receives two different
  haplotypes from the fixed and assigned phasings;
- the second term is the number of distinct listed phasings of each junction $j$ that the
  assignment uses, summed over junctions;
- the search stops at the first assignment with no conflict that uses every listed phasing, and
  in any case after 50 000 assignments (`MAX_PHASING_ASSIGNMENTS`).

**Why.** Applying one `PHASING` to every traversal mislabelled 0.995 Mb of the 5.37 Mb
bench_nested `ctrl_dup_s702` truth walk. Re-scored, its committed λ = 0 result went from phasing
accuracy 0.815 to 1.0, haplotype Jaccard from 0.687 to 1.0, and composite from 0.744 to 0.781
(RCK `CHANGELOG.md`, Fixed). Of the 22 bench_nested and step-sweep truth walks of that time, only
`ctrl_dup_s702` changed. The simulator's own fix, inherited haplotype labels (amplicon-allele-analysis
`hap_labels: inherit`), removes such junctions from new runs. This rule keeps old runs scorable.

### 2.3 Which truth species are scored (scope)

> **Provenance:** (c) new, RCK yubi_dev@3ba7627. **Code:** `rck/ecdna/evaluate.py::fraction_in_regions`,
> `::evaluate`; `rck/ecdna/aa_eval.py::regions_of`, `::score`.

A run covers one amplicon subgraph, but a simulation can hold several species. A truth species
$T$ is in scope when

```math
f(T) = \frac{\sum_{s \in T} \big\lvert [a_s, b_s] \cap \bigcup_{(c_s, r_0, r_1) \in R} [r_0 - \tau_s,\ r_1 + \tau_s] \big\rvert}{\lvert T \rvert} \ \ge\ \phi
```

where

- $R$ is the run's regions:
  - for `rck-ecdna`, `summary.json` `graph.regions` (`rck/ecdna/subgraph_io.py::ECGraph.regions`): the subgraph's
    segments split wherever segments on one chromosome are more than 1 Mb apart;
  - for `aa_eval`, the answer's `Interval` lines, or else the graph's segments merged across
    gaps up to 1 Mb;
- each region is padded by $\tau_s$;
- $\phi$ is `--min-truth-fraction` (default 0.5);
- $f(T)$ is capped at 1.

Out-of-scope species are listed in `evaluation.json` `truth_species` with `in_scope: false` and
`fraction_in_subgraph`. They are not scored and do not count as missed. `--all-truth-species`
scores every species.

---

## 3. `rck-ecdna-eval`: species metrics, matching, recovery matrix and composite

### 3.1 Base-pair overlap: segment recall, precision and Jaccard; haplotype Jaccard

> **Provenance:** (c) new, RCK yubi_dev@c144584. A base-pair Jaccard index between reconstructed
> and true intervals is a standard set-similarity score [Jaccard 1912](references.md#jaccard1912); the
> multiplicity-weighted, haplotype-resolved form is ours. **Code:**
> `rck/ecdna/evaluate.py::_overlap_bp`, `::jaccard`, `::score_pair`.

For a truth species $T$ and an inferred species $C$:

```math
O(T, C) = \sum_{c} \sum_{x} \min\big(m_T(c, x),\ m_C(c, x)\big),
\qquad
O_H(T, C) = \sum_{c} \sum_{h} \sum_{x} \min\big(m_T(c, x, h),\ m_C(c, x, h)\big)
```

```math
\mathrm{segment\_recall} = \frac{O}{\lvert T \rvert},\qquad
\mathrm{segment\_precision} = \frac{O}{\lvert C \rvert},\qquad
\mathrm{segment\_jaccard} = \frac{O}{\lvert T \rvert + \lvert C \rvert - O},\qquad
\mathrm{haplotype\_jaccard} = \frac{O_H}{\lvert T \rvert + \lvert C \rvert - O_H}
```

where

- $O$ counts bases both walks cover, with multiplicity and ignoring haplotype;
- $O_H$ counts them per haplotype;
- $\lvert T \rvert$ and $\lvert C \rvert$ are the walk lengths with multiplicity (§1);
- a Jaccard index with an empty union is 0.

In `evaluate.py` the segment Jaccard is the field `jaccard_total` and the haplotype Jaccard is
`jaccard_hap`; the recovery matrix renames them (`_RENAMED`). Coverage does not depend on
segmentation, so these scores are the same on steps and on runs.

### 3.2 Phasing accuracy

> **Provenance:** (c) new, RCK yubi_dev@4e53747. **Code:** `rck/ecdna/evaluate.py::score_pair`,
> `::_swap_hap`.

```math
\mathrm{phasing\_accuracy} = \frac{O_H(T, C)}{O(T, C)},
\qquad
\mathrm{phasing\_accuracy\_swap\_invariant} = \frac{\max\big(O_H(T, C),\ O_H(T, \bar C)\big)}{O(T, C)}
```

where $\bar C$ is $C$ with A and B exchanged on every step. Both scores are 0 when $O = 0$. The
swap-invariant form is for real data, whose phase labels need not match the truth's.

### 3.3 Junction (novel adjacency) recall, precision and F1, unphased and phased

> **Provenance:** (c) new, RCK yubi_dev@c144584. Matching adjacencies by both extremities within a
> tolerance follows CoRAL's breakpoint matching [CoRAL 3.1.1 code](references.md#coral_code)
> (`coral/scoring/scoring_utils.py::is_breakpoint_match`). Here the tolerance is inclusive and the
> phased variant is added. **Code:** `rck/ecdna/evaluate.py::junction_scores`, `::_junction_match`,
> `::_unique`, `::_f1`.

Let $\mathcal J_T$ and $\mathcal J_C$ be the *distinct* joins of the two walks (§1;
duplicates are removed by exact phased equality, `_unique`). For a match relation $\sim$ (unphased
or phased, tolerance $\tau_s$):

```math
\mathrm{recall} = \frac{\lvert \lbrace j \in \mathcal J_T : \exists\, k \in \mathcal J_C,\ j \sim k \rbrace \rvert}{\lvert \mathcal J_T \rvert},
\quad
\mathrm{precision} = \frac{\lvert \lbrace k \in \mathcal J_C : \exists\, j \in \mathcal J_T,\ k \sim j \rbrace \rvert}{\lvert \mathcal J_C \rvert},
\quad
F_1 = \frac{2\,\mathrm{recall}\cdot\mathrm{precision}}{\mathrm{recall} + \mathrm{precision}}
```

Conventions:

- recall is 1 when $\mathcal J_T = \emptyset$, and precision is 1 when $\mathcal J_C = \emptyset$;
- $F_1 = 0$ when recall + precision = 0;
- a linear inferred walk has no wrap-around join.

Columns: `junction_recall`, `junction_precision`, `junction_f1` (unphased) and
`phased_junction_*` (phased).

### 3.4 Order match and topology match

> **Provenance:** (c) new, RCK yubi_dev@c144584. The linear case is (c), yubi_dev@10416aa.
> **Code:** `rck/ecdna/evaluate.py::order_match`, `::_runs_equal`, `::_flip`.

```math
\mathrm{order\_match}(T, C) = \mathbb 1\Big[\ \mathrm{circ}(T) = \mathrm{circ}(C)\ \wedge\ \exists\, k,\ \exists\, \rho \in \mathcal R :\ \mathrm{runs}(T) \approx_{\tau} \rho\big(\mathrm{runs}(C^{(k)})\big)\ \Big]
```

where

- $C^{(k)}$ ranges over $C$'s orderings; `order_match_ordering` records the first $k$ that
  matches;
- $\mathcal R$ is every rotation of the run sequence and of its reverse. Reversal reverses the
  order and flips every orientation. A linear walk is not rotated;
- $\approx_\tau$ means the same number of runs and, position by position, the same chromosome,
  haplotype and orientation, with start and end each within $\tau$;
- $\mathrm{circ}(\cdot)$ is the topology.

`topology_match` is $\mathbb 1[\mathrm{circ}(T) = \mathrm{circ}(C)]$, and `inferred_topology`
names the inferred one. Truth species are circular, so a linear inferred walk never counts as an
order match.

### 3.5 Copy-number error

> **Provenance:** (c) new, RCK yubi_dev@c144584. The integer variant is (c), yubi_dev@c5e0b9c
> (with the refit). **Code:** `rck/ecdna/evaluate.py::score_pair`.

```math
\mathrm{cn\_rel\_err} = \frac{\lvert w_C - w_T \rvert}{w_T},
\qquad
\mathrm{cn\_rel\_err\_int} = \frac{\lvert w^{\mathrm{int}}_C - w_T \rvert}{w_T}
```

where

- $w_C$ is the inferred species' weight after the copy-number refit (`weight` in
  `decompositions.tsv`);
- $w^{\mathrm{int}}_C$ is its integer weight before the refit (`weight_int`);
- $w_T$ is the truth copies per cell (§2.1).

Both are undefined when $w_T = 0$. Reporting both makes the effect of the refit visible per
species.

### 3.6 One-to-one matching of truth to inferred species

> **Provenance:** (c) new, RCK yubi_dev@c144584. **Code:** `rck/ecdna/evaluate.py::score_decomposition`.

All pairs $(T, C)$ are sorted by haplotype Jaccard, highest first. The pairs are then accepted
greedily: a pair is skipped if either side is already used or if its haplotype Jaccard is 0. A
truth species without a partner is *missed*; its row has no inferred species and composite 0. An
inferred species without a partner is *false*. Matching on haplotype Jaccard means an unphased
answer (every haplotype `.`) matches nothing. That is why the method-agnostic scorer also matches
haplotype-stripped walks (§4.4).

### 3.7 Composite score

> **Provenance:** (c) new, RCK yubi_dev@4e53747. **Code:** `rck/ecdna/evaluate.py::composite_score`.

```math
\mathrm{composite}(T, C) = \frac{1}{5}\Big( J_{\mathrm{seg}} + F_1^{\mathrm{junc}} + \mathrm{PA} + \mathbb 1[\mathrm{order}] + \max(0,\ 1 - \varepsilon) \Big),
\qquad \mathrm{composite}(T, \varnothing) = 0
```

where

- $J_{\mathrm{seg}}$ is `segment_jaccard` (§3.1);
- $F_1^{\mathrm{junc}}$ is the *unphased* `junction_f1` (§3.3);
- $\mathrm{PA}$ is `phasing_accuracy` (§3.2, as given, not swap-invariant);
- $\mathbb 1[\mathrm{order}]$ is `order_match` (§3.4);
- $\varepsilon$ is `cn_rel_err` (§3.5); a missing $\varepsilon$ gives 0 for that term.

The weights are equal and the score lies in [0, 1].

**Why.** The species metrics measure five different things: which bases, which novel
adjacencies, which haplotype, in what order, and how many copies. Choosing a step penalty, or
comparing two decompositions, needs a single number per truth species, so that a better answer
on one axis and a worse one on another can be ranked. No published ecDNA benchmark gives one:
AA [Deshpande et al. 2019](references.md#deshpande2019) and CoRAL [Zhu et al. 2024](references.md#zhu2024) report each metric separately. Equal weights
are a choice, not a fitted value. `composite_h0` (§4.4) is the same score without phasing, for
methods that do not phase.

**Effect.** The composite is the score behind every `--step-penalty` choice and every pass/fail
in `results/bench_nested/`: its λ rule and failure threshold are defined on it (that README). It
also separated the demo7 diagnosis into its parts: with `--block-flip unoriented` both nested
species reach 0.994 and 0.9996, against 0.43 and 0 before (RCK `CHANGELOG.md`).

**A known blind spot.** On `nest_cross_3sp_s403` the composite scored the answer with a 2× scale
error (0.592) above the answer at the right scale (0.545) (RCK `CHANGELOG.md`, "Assessment
(2026-10-31)"). The copy-number term is capped at 0, so a large copy-number error weighs no more
than a moderate one.

### 3.8 The recovery matrix, `best_*` and `mean_composite`

> **Provenance:** (c) new, RCK yubi_dev@4e53747. **Code:** `rck/ecdna/evaluate.py::recovery_matrix`,
> `::evaluate`.

`decompositions.tsv` holds ranked decompositions $D_1, D_2, \dots$, with rank 1 being RCK's
answer. Unless `--include-low-weight` is given, species flagged `low_weight` are left out
(`read_inferred`). For each in-scope truth species $T$:

```math
\mathrm{composite}(T) = \mathrm{composite}\big(T,\ \mu_{1}(T)\big),
\qquad
\mathrm{best\_composite}(T) = \max_{r} \mathrm{composite}\big(T,\ \mu_{r}(T)\big)
```

```math
\mathrm{mean\_composite} = \frac{1}{\lvert \mathcal T_{\mathrm{in}} \rvert} \sum_{T \in \mathcal T_{\mathrm{in}}} \mathrm{composite}(T)
```

where

- $\mu_r(T)$ is $T$'s partner in decomposition $r$ under the matching of §3.6, or $\varnothing$;
- the maximum runs over every reported decomposition. `best_rank` is the arg max, with ties going
  to the lower rank, and `best_inferred` is its matched species;
- $\mathcal T_{\mathrm{in}}$ is the in-scope truth species.

`best_composite` is an upper bound on what the ranking could have picked from the reported
alternatives. `evaluation.json` also reports:

- `first_rank_with_all_orders_matched`: the first rank in which every in-scope truth species
  matches in order;
- `best_mean_jaccard_hap` and its rank.

**Failures.** If the truth files cannot be read, `rck-ecdna` still finishes: the reconstruction
outputs are written, and `evaluation.json` records the error and traceback.

---

## 4. The method-agnostic scorer `aa_eval`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409; `error_free_merged` and the CoRAL-code headline
> are from yubi_dev@41686bc. **Code:** `rck/ecdna/aa_eval.py::score`, `::main`.

The comparison with AA, CoRAL and the original RCK needs one scorer that reads every method's
output the same way. `aa_eval` reads AA-format cycles and graph files only, and none of
`rck-ecdna`'s own files. Species matching and the species scores are §3's functions, called
unchanged (`evaluate.score_decomposition` and `evaluate.score_pair`). Truth species are read and
scoped as in §2. `aa_eval` adds the AA-file reading, the graph- and sample-level metrics, CoRAL's
cycle metrics, and two scoring modes.

### 4.1 Reading AA files

> **Provenance:** (a) the AA graph and cycles formats of AmpliconArchitect [Deshpande et al. 2019](references.md#deshpande2019),
> as written by ecSimulator's `ecSim_IO.py` [ecSimulator](references.md#ecsimulator) and by CoRAL 3.1.1
> [CoRAL 3.1.1 code](references.md#coral_code). The optional `Cycle=` keys are (c) new, RCK yubi_dev@e4c2409.
> **Code:** `rck/ecdna/aa_format.py::read_cycles`, `::to_cycles`, `::read_graph`.

- `Interval` lines give the region the method looked at; `Segment` lines its segments. A
  `Cycle=` line is one walk. A walk whose `Segments=` contains segment `0` (AA's source, `0+ … 0-`)
  is a linear path.
- Optional keys, which AA and CoRAL ignore:
  - `Haplotypes=`: one haplotype per step. It makes the answer phased; without it every
    haplotype is `.` and the scores are haplotype-blind. A file is phased only if every cycle,
    and every `--alt-cycles` file, carries it;
  - `Status=low_weight`: the method's own non-species;
  - `Alt_segments=` / `Alt_haplotypes=`: a species' other orderings, `|`-separated;
  - `Species=`: a name.
- Coordinates: AA's are 1-based inclusive, and segments are read as $[\mathrm{start} - 1,\ \mathrm{end} - 1]$,
  the frame in which the simulator's AA truth is read (§2.1).
- `read_cycles` reads `Cycle=` lines only. CoRAL 3.1.1's and CycleExtractor's linear walks are
  written as `Path=` lines and are invisible to it unless rewritten as `Cycle=`, which
  `results/graph_route/score_h2h.py` does (§7).

### 4.2 Truth mixture graph

> **Provenance:** (c) new, amplicon-allele-analysis yubi_dev@40c6c73. Grouping into amplicons, one
> graph and cycles file per amplicon, is (a) ecSimulator's output layout [ecSimulator](references.md#ecsimulator).
> **Code:**
> [`simulations/truth_mixture.py::build`](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/truth_mixture.py),
> `::group_amplicons`, `::mixture_segments`, `::merge_silent`, `::background_junctions`.

The species truth (§2) is one walk per species. The graph truth is a different object: the
amplicon as a method sees it, with every species at its copies per cell plus the background
genome. The simulator writes only one AA graph per species, so `truth_mixture.py` builds the
mixture from the truth tables:

- `segments.tsv`, `clones.tsv` and `proportions.tsv`;
- the species' `truth/aa` files;
- the background junctions in `adjacency_clone.tsv` and `rck.acnt.tsv`.

It writes `truth/aa/mixture_amplicon<N>_graph.txt`, `_cycles.txt` (one cycle per species,
`Copy_count` = PHI, `Species=` naming it) and `_segments.tsv` (allele copy numbers, which an AA
graph cannot hold), plus an index `mixture_amplicons.tsv`.

- **Amplicons follow ecSimulator.** ecSimulator simulates each amplicon on its own and writes one
  `<prefix>_amplicon<N>_graph.txt` / `_cycles.txt` per amplicon; it has no joint multi-amplicon
  graph. Here an amplicon is a group of species whose intervals overlap or abut (`group_amplicons`),
  and each gets its own pair of files. demo6 has two amplicons; demo7's nested pair has one.
- **Segments.** The merged intervals are cut at every species segment boundary, every
  `segments.tsv` boundary and every background breakend inside them. A boundary with no junction
  and no copy-number change is dropped (`merge_silent`), since no method can see it. ecSimulator
  keeps such a boundary for an SV that a later event removed.
- **Copy numbers.** For a mixture segment $g$:

  ```math
  \mathrm{cn}(g) = \sum_{k} U_k\,(a_{k,g} + b_{k,g}) \;+\; \sum_{s} \phi_s\, \kappa_s(g),
  \qquad
  \mathrm{cn}(d) = \sum_s \phi_s\, \nu_s(d) \;+\; \sum_k U_k\, \beta_k(d)
  ```

  where

  - $U_k$ is clone $k$'s share of cells in the sample;
  - $a_{k,g}, b_{k,g}$ are its background allele copy numbers over $g$;
  - $\phi_s$ is species $s$'s copies per cell (PHI);
  - $\kappa_s(g)$ is how many times $s$'s walk covers $g$;
  - for a discordant edge $d$, $\nu_s(d)$ is how many times $s$ traverses $d$ and $\beta_k(d)$ is
    clone $k$'s copies of a background junction with an end in the amplicon.

  Concordant edges carry what node balance leaves.
- **Alleles** are known for the first tumor sample only, which `segments.tsv` describes.
- **Where the files are.** The tracked demo5–demo8 runs carry theirs in amplicon-allele-analysis
  (`simulations/runs/demo*/truth/aa/mixture_*`). The bench simulations (`/n/fs/yubi-iw/runs/`)
  carry theirs in place. `results/bench_nested/diag/aa_regression.py` uses
  `truth/aa/mixture_amplicons.tsv` wherever it exists and otherwise builds one (demo4: under
  `/n/fs/yubi-iw/runs/bench_nested_aa/mixtures/`).

Options for the graph route, opt-in, with byte-identical defaults: `--reads-per-copy R` (read
counts `round(CN × R)` on breakpoint edges), `--depth-per-copy D` (AverageCoverage = CN × D) and
`--path-constraints`. The last writes the run's long-read paths whose junctions are mixture
discordant edges, as CoRAL 3.1.1 `path_constraint` lines in CoRAL's sort order.

### 4.3 Graph-level metrics (`aa_graph_metrics.tsv`)

> **Code:** `rck/ecdna/aa_eval.py::graph_metrics`, `::edges_match`, `::scoped_edges`,
> `::scoped_sequences`, `::_segment_match`.

The metrics are computed once per truth amplicon $A$ that contains an in-scope species. The
following symbols are used throughout this section:

| Symbol | Meaning |
|---|---|
| $I_A$ | the amplicon's intervals (1-based, from the mixture cycles file) |
| $E_T$ | the mixture graph's discordant edges (novel adjacencies) |
| $E_C$ | the method's discordant edges with at least one breakend within $I_A$ padded by $\tau_g$ |
| $S_T$ | the mixture graph's sequence edges (segments) |
| $S_C$ | the method's sequence edges overlapping $I_A$ by more than $\tau_g$ bp (a flank that only touches the amplicon is not part of it) |
| $e \sim f$ | same orientations, both breakpoints within $\tau_g$, in either direction (`edges_match`; CoRAL's `is_breakpoint_match` with $\le$ instead of $<$) |
| $g \simeq g'$ | same chromosome, starts within $\tau_g$ and ends within $\tau_g$ |

#### `bp_accuracy`: breakpoint-graph accuracy

> **Provenance:** (a) CoRAL's breakpoint-graph accuracy [Zhu et al. 2024](references.md#zhu2024). In CoRAL 3.1.1 it is the
> count `n_breakpoint_edges_covered` of `coral/scoring/score_simulation.py::score_reconstruction`
> over `n_breakpoint_edges`, counted by `coral/scoring/scoring_utils.py::find_overlapping_bp_edges`
> [CoRAL 3.1.1 code](references.md#coral_code).

```math
\mathrm{bp\_accuracy} = \frac{\lvert \lbrace e \in E_T : \exists f \in E_C,\ e \sim f \rbrace \rvert}{\lvert E_T \rvert}
```

It is undefined when $E_T = \emptyset$. Differences from CoRAL's code:

1. CoRAL counts against each reconstructed amplicon graph and keeps the amplicon with the most
   matches. Here the method's whole graph is scoped to $I_A$.
2. CoRAL uses a strict $<$.

#### `novel_adjacency_fnr`, `novel_adjacency_fdr`

> **Provenance:** (a) RCK's novel-adjacency recovery [Aganezov & Raphael 2020](references.md#aganezov2020). The rates are written
> here as FNR = 1 − recall and FDR. The FDR is (c), chosen instead of an FPR, which would need a
> negative set.

```math
\mathrm{FNR} = 1 - \mathrm{bp\_accuracy},
\qquad
\mathrm{FDR} = \frac{\lvert \lbrace f \in E_C : \nexists\, e \in E_T,\ f \sim e \rbrace \rvert}{\lvert E_C \rvert}
```

#### `error_free`: AA's error-free reconstruction

> **Provenance:** (a) AmpliconArchitect's simulation benchmark [Deshpande et al. 2019](references.md#deshpande2019).

```math
\mathrm{error\_free} = \mathbb 1\Big[ \forall g \in S_T\, \exists g' \in S_C : g \simeq g' \ \wedge\ \forall g' \in S_C\, \exists g \in S_T : g' \simeq g \ \wedge\ \forall e \in E_T\, \exists f \in E_C : e \sim f \ \wedge\ \forall f \in E_C\, \exists e \in E_T : f \sim e \Big]
```

In words: every true segment is found with both ends within tolerance, there is no false
segment, and every true discordant edge is found with no false one. `rck-ecdna`'s imputed
breakpoints (the subgraph's 10 kb grid) count as false segments.

#### `error_free_merged`: the same rule after merging benign splits

> **Provenance:** (c) modified from AA's rule, RCK yubi_dev@41686bc. **Code:**
> `rck/ecdna/aa_eval.py::merge_benign_splits`.

```math
\mathrm{error\_free\_merged} = \mathrm{error\_free}\ \text{with}\ S_T, S_C\ \text{replaced by}\ M(S_T), M(S_C)
```

where $M(S)$ is the graph's sequence edges with every benign split removed. Two segments $g, g'$
on one chromosome with $\mathrm{start}(g') = \mathrm{end}(g) + 1$ become one segment when:

- a concordant edge joins $\mathrm{end}(g)+$ to $\mathrm{start}(g')-$; and
- no discordant or source edge has a breakend within 1 bp of either side of that boundary.

Merged segments are then scoped to $I_A$ (overlap > $\tau_g$). The discordant-edge part of the
rule is unchanged. Columns `n_{true,inferred,false}_segments_merged` and
`n_segments_found_merged` count the merged segments. `error_free` stays the headline.

**Why.** AA's rule asks for every segment boundary. AA and CoRAL graphs, like `rck-ecdna`'s, also
cut segments where the copy number changes without a breakpoint, and `rck-ecdna`'s subgraphs add
a 10 kb imputed grid. The strict rule therefore scores every method near 0 for boundaries that no
breakpoint explains. On the grown bench at λ = 0.2, 1 of 42 truth amplicons is `error_free`
(demo6's chr17 amplicon), and all 532 unmatched inferred segment boundaries are imputed
extremities of the subgraph. `error_free_merged` merges those boundaries first, in both graphs,
and so tests the breakpoints and the segments they delimit. A boundary that carries a discordant
or source edge is never merged, so a missed or false breakpoint still fails it.

**Effect.** 38 of 42 truth amplicons pass `error_free_merged` (RCK `results/bench_nested/README.md`,
"Method-agnostic metrics (2026-11-03)"). Two of the four failures use an imputed edge
(`ctrl_dup_s710`, `nest_samehap_6to1_s804`). The other two miss a copy-number-2 flank segment
inside the truth interval that the subgraph shortens (`demo4_amp0`, `nest_samehap_6to1_s803`).

#### `segment_cn_mae`, `segment_cn_rel_err`, `segment_cn_mae_total`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. A length-weighted absolute copy-number error in
> the spirit of RCK's segment copy-number distance [Aganezov & Raphael 2020](references.md#aganezov2020). **Code:**
> `rck/ecdna/aa_eval.py::segment_cn_error`, `::_piecewise`.

Over the truth amplicon's mixture segments $k$, with span $[s_k, e_k]$, length $\ell_k$, total
copy number $n_k$ and alleles $(\alpha_k, \beta_k)$:

```math
\mathrm{err}_{\mathrm{tot}} = \sum_k \sum_{x = s_k}^{e_k} \big\lvert \hat n(x) - n_k \big\rvert,
\qquad
\mathrm{err}_{\mathrm{al}} = \sum_k \sum_{x = s_k}^{e_k} \Big( \big\lvert \hat\alpha(x) - \alpha_k \big\rvert + \big\lvert \hat\beta(x) - \beta_k \big\rvert \Big)
```

```math
\mathrm{segment\_cn\_mae} = \frac{\mathrm{err}}{\sum_k \ell_k},
\qquad
\mathrm{segment\_cn\_rel\_err} = \frac{\mathrm{err}}{\sum_k n_k\, \ell_k},
\qquad
\mathrm{segment\_cn\_mae\_total} = \frac{\mathrm{err}_{\mathrm{tot}}}{\sum_k \ell_k}
```

where

- $\hat n(x)$, $\hat\alpha(x)$ and $\hat\beta(x)$ are the method's total and allele copy numbers
  at base $x$. Totals come from `--segment-cn` (`cn`, or `cn_a + cn_b`), else from the graph's
  sequence edges; alleles come from `--segment-cn`;
- bases the method gives no segment count as copy number 0;
- $\mathrm{err} = \mathrm{err}_{\mathrm{al}}$ (`cn_error_mode` = `allele`) when both sides have
  allele copy numbers, else $\mathrm{err}_{\mathrm{tot}}$ (`total`).

#### `is_violations`, `is_violations_truth`

> **Provenance:** (a) the infinite-sites assumption of RCK [Aganezov & Raphael 2020](references.md#aganezov2020), as enforced by
> RCK's constraint that an extremity takes part in at most one novel adjacency. This count of
> violations is (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_eval.py::is_violations`.

Distinct discordant edges are taken first. Their breakends are grouped per (chromosome, strand)
by single linkage at gap $\le \tau_g$, and each group is one extremity. The count is the number
of extremities used by more than one distinct edge. A fold-back (both ends on one node) counts as
one edge. `is_violations_truth` is the same count on $E_T$. Simulated duplications reuse
extremities, so the truth itself violates the assumption: 25 violations on the bench, the same 25
in RCK's graphs (RCK `results/bench_nested/README.md`).

### 4.4 Species-level metrics (`aa_species_metrics.tsv`)

> **Code:** `rck/ecdna/aa_eval.py::strict_rows`, `::best_of_rows`, `::phased_part`, `::blind_part`.

Each in-scope truth species gets one row per mode (§4.6). There are two families of columns:

- **The `rck-ecdna-eval` columns** of §3, at $\tau_s$, under the recovery-matrix names, plus
  `exact_cycle_phased`. They are computed for phased answers only. For an unphased answer the
  haplotype columns and `composite` are empty.
- **Haplotype-blind columns**, from a *separate* one-to-one matching (§3.6) of
  haplotype-stripped walks (`strip_haplotypes` sets every haplotype to `.`). The matched species
  is `inferred_h0`.

The cycle metrics are computed on haplotype-stripped walks: the CoRAL metrics on the walks
merged into reference runs (§1), the other columns on the walks' steps.

#### `composite_h0`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_eval.py::composite_h0`.

```math
\mathrm{composite\_h0}(T, C) = \frac{1}{4}\Big( J_{\mathrm{seg}}^{0} + F_1^{0} + \mathbb 1[\mathrm{order}^{0}] + \max(0,\ 1 - \varepsilon^{0}) \Big),
\qquad \mathrm{composite\_h0}(T, \varnothing) = 0
```

The superscript 0 marks scores computed on the haplotype-stripped pair matched by the
haplotype-blind matching: `segment_jaccard_h0`, `junction_f1_h0`, `order_match_h0` and
`cn_rel_err_h0`.

**Why.** AA, CoRAL, CycleExtractor and the original RCK on an AA graph do not phase, and the
phased matching of §3.6 matches nothing to an unphased answer. `composite_h0` is the composite
(§3.7) without phasing accuracy, and with its order match haplotype-blind. It makes `rck-ecdna`
comparable with them on equal terms. **Effect.** On the grown bench at λ = 0.2, strict
`composite_h0` is 0.835 against a composite of 0.853. Nothing in the composite rests on phasing
alone (RCK `results/bench_nested/README.md`).

#### CoRAL's cycle-interval overlap: `cycle_overlap_coral_code` (headline) and `cycle_overlap_supp` (secondary)

> **Provenance:** (a) CoRAL [Zhu et al. 2024](references.md#zhu2024). The code version is CoRAL 3.1.1
> `coral/scoring/scoring_utils.py::get_fragments_similarity_unweighted` [CoRAL 3.1.1 code](references.md#coral_code). The
> supplementary-methods version follows the supplement's description, as recorded in the
> previous version of this documentation [Zhu et al. 2024](references.md#zhu2024). **Code:**
> `rck/ecdna/aa_eval.py::cycle_overlap_coral_code`, `::cycle_overlap`, `::_bins`, `::_counts`.

Joint bins are taken per chromosome: the consecutive points of the set
$\lbrace a_s, b_s + 1 \rbrace$ over every step of both walks. For a bin $\beta$ with length
$\lvert\beta\rvert$, $e_T(\beta)$ and $e_C(\beta)$ are the numbers of steps of $T$ and $C$ that
cover it (multiplicity).

```math
\mathrm{cycle\_overlap\_coral\_code}(T, C) = 1 - \frac{\sum_\beta \lvert e_T(\beta) - e_C(\beta) \rvert\, \lvert\beta\rvert}{\sum_\beta \max\big(e_T(\beta),\ e_C(\beta)\big)\, \lvert\beta\rvert}
```

```math
\mathrm{cycle\_overlap\_supp}(T, C) = \frac{\lvert U_T \cap U_C \rvert}{\lvert U_T \cup U_C \rvert}
```

where $U_X$ is the union of $X$'s step intervals, without multiplicity. Intervals on different
chromosomes never overlap. The two are equal when no interval is used twice.

The CoRAL code differs from the supplement because it counts multiplicity. The code's version is
the headline because it produced CoRAL's published numbers. In CoRAL's code the true side is
every true cycle of the amplicon, and the minimum distance over the reconstructed cycles is kept.
Here the score is per truth species against its matched cycle (strict mode), or the best over
candidates (best-of).

#### CoRAL's cyclic longest common subsequence: `cyclic_lcs_coral_code` (headline) and `cyclic_lcs_supp` (secondary)

> **Provenance:** (a) CoRAL [Zhu et al. 2024](references.md#zhu2024). The code version is CoRAL 3.1.1
> `coral/scoring/scoring_utils.py::compute_normalized_longest_cycle_subsequence` with
> `::find_longest_subsequence` [CoRAL 3.1.1 code](references.md#coral_code). The LCS itself is the classical dynamic
> programme [Wagner & Fischer 1974](references.md#wagnerfischer1974). **Code:** `rck/ecdna/aa_eval.py::cyclic_lcs`,
> `::_weighted_lcs`, `::_labels`.

The bins are those defined for the overlap. Let $K = \lbrace \beta : e_T(\beta) > 0,\ e_C(\beta) > 0 \rbrace$
be the bins both walks cover. Write each walk as the sequence $\sigma_X$ of $(\beta, o)$ over the
bins of $K$ in traversal order; a step on `−` lists its bins in reverse. Then:

```math
L^\star = \max_{\rho \in \mathcal R(\sigma_C)} \mathrm{WLCS}\big(\sigma_T,\ \rho\big),
\qquad
\mathrm{cyclic\_lcs\_supp} = \frac{L^\star}{\lvert T \rvert},
\qquad
\mathrm{cyclic\_lcs\_coral\_code} = \frac{L^\star}{\sum_\beta \max\big(e_T(\beta), e_C(\beta)\big)\, \lvert\beta\rvert}
```

where

- $\mathrm{WLCS}$ is the longest common subsequence in which a match of $(\beta, o)$ scores
  $\lvert\beta\rvert$;
- $\mathcal R(\sigma_C)$ is every rotation of $\sigma_C$ and of its reverse, with orientations
  flipped (only $\sigma_C$ itself for a linear walk);
- the code's denominator sums over *all* joint bins, including those only one walk covers.

The two normalisations differ: the code normalises by the bins' total with multiplicity, the
supplement by the true cycle's length. The code's version is the headline. It produced CoRAL's
published numbers, and it charges a merged walk for its extra bases. The supplement's version
gives 1 to a merged walk that contains the true cycle, which flatters exactly the nested-species
failure the bench tests. **Effect** at λ = 0.2 on the grown bench: code version 0.871 strict and
0.958 best-of; supplement version 0.901 strict and 0.986 best-of (RCK `CHANGELOG.md`, Changed;
`results/bench_nested/README.md`).

#### `length_error`

> **Provenance:** (a) CoRAL [Zhu et al. 2024](references.md#zhu2024). The ratio $R$ is in CoRAL 3.1.1's code: the second
> return value of `get_fragments_similarity_unweighted`, stored as `reconstruction_length_ratio`.
> The $\log_2(1 + R)$ form follows the supplementary methods, as recorded in the previous version
> of this documentation. **Code:** `rck/ecdna/aa_eval.py::length_error`.

```math
R = \frac{\lvert C \rvert - \lvert T \rvert}{\lvert T \rvert},
\qquad
\mathrm{length\_error} = \log_2(1 + R) = \log_2 \frac{\lvert C \rvert}{\lvert T \rvert}
```

The error is signed and is 0 at the right length. Best-of takes the value of smallest absolute
size.

#### `exact_cycle`, `exact_cycle_phased`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409, after AA's and CoRAL's exact-reconstruction
> criterion [Deshpande et al. 2019](references.md#deshpande2019) [Zhu et al. 2024](references.md#zhu2024). **Code:** `rck/ecdna/aa_eval.py::_exact`.

```math
\mathrm{exact\_cycle}(T, C) = \mathrm{order\_match}_{\tau_g}(T^{0}, C^{0}),
\qquad
\mathrm{exact\_cycle\_phased}(T, C) = \mathrm{order\_match}_{\tau_g}(T, C)
```

This is §3.4 at the graph tolerance $\tau_g$, on haplotype-stripped walks or with haplotypes,
over every ordering of $C$. Every run end lies within 100 bp, so every junction does too.

### 4.5 Sample-level metrics (`aa_sample_metrics.tsv`)

> **Code:** `rck/ecdna/aa_eval.py::score` (sample rows), `::cn_explained`, `::subwalks_satisfied`.

There is one row per mode. $\mathcal S$ is the strict answer, $\mathcal T_{\mathrm{in}}$ the
in-scope truth species, and $\mathcal C$ the candidates (§4.6).

| Column | Definition | Provenance |
|---|---|---|
| `species_count_error` | $\lvert \mathcal S \rvert - \lvert \mathcal T_{\mathrm{in}} \rvert$ (strict only) | (c) new |
| `n_missed` | strict: truth species without a partner in the matching (phased, or haplotype-blind for an unphased answer); best-of: truth species with `composite_h0` $\le 0$ | (c) new |
| `n_false` | strict species without a partner (strict only) | (c) new |
| `mean_composite`, `mean_composite_h0`, `mean_cycle_overlap_coral_code`, `mean_cyclic_lcs_coral_code`, `mean_cycle_overlap_supp`, `mean_cyclic_lcs_supp` | means over $\mathcal T_{\mathrm{in}}$ of the species columns (`mean_composite` for phased answers only) | (c) new |
| `n_exact_cycles` | $\sum_T$ `exact_cycle` | (c) new |
| `cn_explained_top1`, `cn_explained_top3` | below | (a) CoRAL, modified |
| `subwalks_satisfied`, `subwalks_satisfied_support`, `n_subwalks` | below | (a) CoRAL, modified |
| `error_free`, `error_free_merged` | 1 only if every scored amplicon passes | (c) aggregation of §4.3 |
| `bp_accuracy`, `novel_adjacency_fnr`, `novel_adjacency_fdr`, `segment_cn_mae` | means over the scored amplicons | (c) aggregation of §4.3 |
| `is_violations` | sum over the scored amplicons | (c) aggregation of §4.3 |

#### `cn_explained_topk`: copy number explained by the heaviest cycles

> **Provenance:** (a) CoRAL's copy-number ratios `best_copy_number_ratio` and
> `top_three_copy_number_ratio` [Zhu et al. 2024](references.md#zhu2024), CoRAL 3.1.1
> `coral/scoring/scoring_utils.py::get_cycle_copy_number_ratio` [CoRAL 3.1.1 code](references.md#coral_code). (c) modified
> as described below.

```math
\mathrm{cn\_explained\_top}k = \frac{\sum_{C \in \mathrm{top}_k} w_C\, \lvert C \rvert}{\sum_{g \in S_C,\ \mathrm{cn}(g) \ge \gamma} \mathrm{cn}(g)\, \ell(g)}
```

where

- $\mathrm{top}_k$ is the $k$ circular walks of the method's cycles file with the largest
  $w_C \lvert C \rvert$ (copy count × length with multiplicity; paths are left out, as CoRAL
  scores cyclic walks only);
- the denominator sums over the method's own graph's sequence edges $g$ with copy number at
  least $\gamma$ = `--cn-explained-min-cn` (5, CoRAL's constant).

No truth is needed. Difference from CoRAL's code: CoRAL returns 1 when the number of cycles is at
most $k$ (`if n >= len(cycle_to_ratio): return 1.0`), whereas this computes the share whatever the
number of cycles. CoRAL's lengths are half-open BED lengths ($\mathrm{end} - \mathrm{start}$),
and these are inclusive.

#### `subwalks_satisfied`: long-read subwalk constraints

> **Provenance:** (a) CoRAL's long-read subwalk (path) constraints [Zhu et al. 2024](references.md#zhu2024). The
> satisfaction share is (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_eval.py::read_paths`,
> `::path_in_cycle`, `::subwalks_satisfied`.

The paths come from `--paths`, a subgraph `paths.tsv`. Only reads with at least two *novel*
joins are kept; a reference join is the same chromosome, opposite strands, and the `−` position
equals the `+` position + 1 within 2 bp. This matches the paths `rck-ecdna` uses. A read $p$, with
novel joins $(j_1, \dots, j_m)$, is satisfied by a walk whose novel joins in traversal order
(haplotype-blind) contain $(j_1, \dots, j_m)$ consecutively, in either direction, cyclically for a
circle, every extremity within $\tau_g$.

```math
\mathrm{subwalks\_satisfied} = \frac{\lvert \lbrace p : \exists C \in \mathcal X,\ p \preceq C \rbrace \rvert}{\lvert P \rvert},
\qquad
\mathrm{subwalks\_satisfied\_support} = \frac{\sum_{p \text{ satisfied}} \mathrm{NOV}_p}{\sum_p \mathrm{NOV}_p}
```

where

- $\mathcal X$ is the mode's walks: the strict answer, or every candidate in best-of;
- $P$ is the kept reads;
- $\mathrm{NOV}_p$ is the read's support (`#NOV`).

### 4.6 Modes: `strict` and `best_of`, and `--strict-rule`

> **Provenance:** `strict` is (c) new, RCK yubi_dev@e4c2409. `best_of` is (a) CoRAL's scoring,
> which keeps each metric's best value over the reconstructed cycles independently
> (`coral/scoring/score_simulation.py::score_reconstruction`: fragment overlap and LCS are
> maximised separately) [CoRAL 3.1.1 code](references.md#coral_code), applied per truth species here. **Code:**
> `rck/ecdna/aa_eval.py::strict_select`, `::strict_rows`, `::best_of_rows`.

**Strict.** This is the method's own answer $\mathcal S \subseteq$ its cycles, chosen by
`--strict-rule`:

| Rule | $\mathcal S$ |
|---|---|
| `status` | cycles without `Status=low_weight`; for `rck-ecdna`, rank 1 without `low_weight` |
| `cn-frac:F` | $\lbrace C : w_C \ge F \cdot \max_{C'} w_{C'} \rbrace$ |
| `top-cn:X` | the fewest heaviest cycles by $w_C \lvert C \rvert$ whose cumulative mass reaches $X \sum_C w_C \lvert C \rvert$ |
| `all` | every cycle |
| `auto` (default) | `status` when the file carries `Status=`, else `cn-frac:0.1` |

AA's and CoRAL's cycles files list many low-copy cycles that the methods themselves do not report
as species. `cn-frac:0.1` approximates `rck-ecdna`'s own significance rule
(`--min-species-weight-frac` 0.1), so the three methods face a comparable threshold. It is not
identical: see §7. Strict mode matches one to one as in §3.6 and counts missed and false species.

**Best-of.** For each truth species $T$, the candidates $\mathcal C(T)$ are every cycle of
`--cycles` and of `--alt-cycles`, `low_weight` included, that overlaps $T$: haplotype Jaccard
$> 0$ when phased, haplotype-blind Jaccard $> 0$ for the blind columns. Each metric $m$ takes its
best value on its own:

```math
m^{\mathrm{best}}(T) = \max_{C \in \mathcal C(T)} m(T, C)
\quad\text{(min for cn\_rel\_err, cn\_rel\_err\_h0; smallest}\ \lvert\cdot\rvert\ \text{for length\_error)}
```

`inferred` (`inferred_h0`) names the candidate with the best `composite` (`composite_h0`). For
`rck-ecdna` the alternatives are rank 1's first $k$ same-fit alternatives (`aa_export`, $k = 2$
by default), the candidates of the set metric (§6.2). Best-of picks per truth species, not per
decomposition, so it is never below that set metric.

**Why two modes.** Strict mode is what a user of the method gets. Best-of is what CoRAL's own
benchmark reports, so the numbers are comparable with its paper, and for `rck-ecdna` it measures
whether the ambiguity report contains the right answer. **Effect** at λ = 0.2 on the grown bench:
best-of mean composite 0.931 against 0.853 strict; 55 of 72 truth species are exact cycles in
strict mode and 61 in best-of (RCK `results/bench_nested/README.md`).

### 4.7 Regression against `rck-ecdna-eval`

> **Code:** `results/bench_nested/diag/aa_regression.py`; `results/bench_nested/collect.py --aa-metrics`.

The AA-format path must reproduce §3 exactly for `rck-ecdna`'s answers. On the 40 bench_nested
cases at λ = 0 and 0.2 (80 pairs), the following are identical at full precision:

- every part of the composite;
- every composite;
- the matched species, `n_species`, `n_truth` and `merged`.

With `--aa-metrics` all 240 (case, λ) match, and the existing columns of `bench_sweep.tsv`,
`pipeline.tsv` and `rerank_validation.tsv` are unchanged. Two things were needed:

- a species' other orderings: `nest_samehap_3to1_s304` matches its order only through ordering
  2, which RCK's own `aa/` files leave out;
- reading the old tables with `float_precision="round_trip"`.

---

## 5. Translators: `aa_export`, `aa_import`, `coral_flank`

### 5.1 `rck-ecdna` answers in AA format (`aa_export RUN OUT`)

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_export.py::export`,
> `::write_cycles`, `::run_alternatives`, `::pool_alternatives`, `::shift_graph`, `::write_segment_cn`.

The files `rck-ecdna` writes itself (`aa/<name>_decomp<k>_cycles.txt`) lack the `low_weight`
status, the other orderings, and the solutions of a kept pool, which the λ > 0 results come from.
`python -m rck.ecdna.aa_export RUN OUT [--lambda L | --pool-index I]` writes these from
`decompositions.tsv` (or `pool/decompositions.tsv`):

- `<name>_cycles.txt`: the answer, one cycle per species, with `Copy_count` = the refitted weight,
  `Haplotypes=`, `Status=`, `Species=` and the other orderings. Its `Interval` lines are the run's
  regions;
- `<name>_alt<k>_cycles.txt`: its first `--alternatives` (2) same-fit alternatives. They come from
  `summary.json` `ambiguity`, or from `pool.json` `same_fit` matched to pool entries by objective
  (within $10^{-3}$), steps and sorted significant weights, as `grow_rules.py` matches them;
- `<name>_graph.txt`: the run's AA graph;
- `<name>_segment_cn.tsv`: allele copy numbers from `flow_segments.tsv`.

Positions are the walk's coordinates + 1 (`SHIFT`). The scorer subtracts 1 when it reads AA
files (§4.1), so it reads back exactly the steps `rck-ecdna-eval` reads. Graph edges match the
truth within 1 bp either way, because a subgraph segment's start and end are not in one
convention.

### 5.2 Subgraph → AA graph (`aa_export subgraph`)

> **Provenance:** (c) new, RCK yubi_dev@2f921a5. The graph layout and sort order are (a) CoRAL 3.1.1's
> and CycleExtractor 1.0.0's parsers [CoRAL 3.1.1 code](references.md#coral_code) [Faizrahnemoon et al. 2026](references.md#cycleextractor); the path-constraint
> tokens are CoRAL's `coral/core_utils.py::path_to_str`. **Code:**
> `rck/ecdna/aa_export.py::export_subgraph`, `::subgraph_observations`, `::_path_edges`;
> `rck/ecdna/aa_format.py::write_graph`, `::coral_path_string`, `::chrom_key`.

The graph route gives every decomposer the same breakpoint graph, so that decompositions are
compared on equal input (`results/graph_route/`).
`python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name N] [--amplicon 1] [--imputed comment|drop]`
writes `<name>_amplicon<N>_graph.txt`, plus `_segment_cn.tsv` and `_export.json`. The graph is the
normalised graph `rck-ecdna` solves on (`subgraph_io.read_subgraph`: phantom segments contracted,
parallel imputed edges merged), with the copy numbers it uses (`scale.build_observations`):

```math
\mathrm{CN}(\text{segment}) = \mathrm{CN\_TOTAL},
\qquad
\mathrm{CN}(\text{novel adjacency}) = \frac{\text{reads}}{\theta}
```

where $\theta$ is the junction reads per copy in `scale.json`. This is the *observed* adjacency
copy number, not the rounded `CN` column of `edges.tsv`.

| Line | Content | Seen by CoRAL / CycleExtractor |
|---|---|---|
| `sequence` | one per segment, CoRAL's sort order; AverageCoverage = depth | yes |
| `concordant` | reference adjacencies between segments; an unobserved join (imputed grid point) gets the smaller segment CN | yes |
| `discordant` | observed novel adjacencies; fold-backs (both ends on one node) as self-loops, listed in `_export.json` | yes |
| `path_constraint` | read paths of `paths.tsv` over segments, reference joins and observed junctions, as CoRAL 3.1.1 tokens; support = `#NOV` rounded | yes |
| `interval` | runs of contiguous segments (what CoRAL's `infer_amplicon_intervals` infers) | ignored |
| `#rck_imputed` | RCK-ecDNA's imputed candidates (no read evidence); `--imputed drop` leaves them out | no |
| `#rck_outside` | joins to a region's outer flank (a node without a segment); CoRAL and CycleExtractor add their own s/t ends there | no |
| `#rck_node`, `#rck_seq`, `#rck_edge`, `#rck_pc`, `#rck_path` | the subgraph's node, segment, adjacency and path ids and orders, so the import is exact | no |

Path constraints are written in CoRAL's sort order because `coral cycle` re-sorts sequence and
concordant edges (`infer_amplicon_intervals`) without renumbering constraints that were already
parsed by file index.

### 5.3 AA graph → subgraph (`aa_import`)

> **Provenance:** (c) new, RCK yubi_dev@2f921a5. The `regenerate` mode reuses (c)
> `subgraph_construct.py`'s imputation rule and amplification threshold (amplicon-allele-analysis).
> **Code:** `rck/ecdna/aa_import.py::import_graph`, `::fit_scale`, `::regenerate_imputed`,
> `::amplified_threshold`, `::weighted_median`.

`python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate] [--cn-profile F] [--background-cn X]`
turns any AA graph into a subgraph folder for `rck-ecdna --inputs subgraph`. The mapping:

- `sequence` becomes a `reference` edge (`WEIGHT` = AverageCoverage, `CN_TOTAL` = CN);
- `concordant` becomes an `extremity` edge;
- `discordant` becomes a `variant` edge with no haplotype label, so every phasing is allowed;
- `source` edges are dropped and listed in `import.json`.

The scale is fitted from the graph itself:

```math
\hat\theta = \operatorname{median}_{d \in D,\ r_d > 0,\ \mathrm{cn}_d > 0} \frac{r_d}{\mathrm{cn}_d},
\qquad
\widehat{\mathrm{dpc}} = \operatorname{wmedian}_{g:\ \mathrm{cov}_g > 0,\ \mathrm{cn}_g > 0}\Big( \frac{\mathrm{cov}_g}{\mathrm{cn}_g};\ \ell_g \Big)
```

where

- $D$ is the discordant edges, or the concordant edges if no discordant edge has reads;
- $r_d$ and $\mathrm{cn}_d$ are an edge's read count and copy number;
- $\mathrm{cov}_g$ is a sequence edge's AverageCoverage, and wmedian is the length-weighted median;
- a graph without reads gets $\hat\theta = 1$ and `WEIGHT` = CN;
- a CoRAL-rewritten graph (`#CoRAL` header) gets no depth fit, because CoRAL's writer puts
  `lr_nc / size` in AverageCoverage;
- a ratio that agrees to $10^{-9}$ with its 12-significant-digit rounding is snapped to it, so
  6.5 comes back as 6.5 (`_snap`).

Imputed candidates, `--imputed`:

- `comment` restores the `#rck_imputed` lines;
- `none` leaves them out, which is what CoRAL and CycleExtractor see;
- `regenerate` proposes closing candidates by `subgraph_construct.py`'s rule:
  - an `intra` edge from each amplified run's end back to its start;
  - an `inter` edge between consecutive amplified runs;
  - "amplified" means above $\min(\mathrm{med} + 1.5\,\mathrm{MAD},\ 1.5\,\mathrm{med})$ of the
    chromosome's `cn_profile.tsv`, or else 1.5 × background, where MAD is the scaled median
    absolute deviation ($1.4826 \times$);
  - a run is split at a copy-number step of at least $\max(2,\ 0.25 \times \text{larger side})$
    (`--boundary-min-copies`, `--boundary-frac`).

`#rck_outside` joins are always restored.

**Round trip.** `aa_import.compare_subgraphs` compares the normalised graphs `rck-ecdna` builds
from two folders: nodes, segments, adjacencies (ids, ends, kind, reads, allowed phasings,
order), observations, read paths and regions. It compares them total-only with junction labels
ignored, to a relative $10^{-12}$. Export then import with `--imputed comment` and the original
`cn_profile.tsv` is exact on 36 of 36 cases (demo5–demo8 and the 32 bench cases). `--imputed
none` is exact apart from the imputed candidates. A graph-only `none` differs only in the flank
background (RCK `results/graph_route/README.md`, "Validation done offline"). Total-only runs on
both sides use `--min-snps-allele 100000000 --soft-junction-hap --soft-junction-penalty 0`, which
ignores junction labels exactly.

The one reader change: an `extremity` row of `edges.tsv` may carry an `AID` (the import fills it
to keep reference-join ids); `subgraph_construct.py` leaves it empty, so its subgraphs read as
before (`rck/ecdna/subgraph_io.py::build_graph`).

### 5.4 The CoRAL flank fix (`coral_flank`)

> **Provenance:** (c) new, RCK yubi_dev@bb15c21. **Code:** `rck/ecdna/coral_flank.py::add_flanks`,
> `::flank_graph`, `::_runs`, `::_renumber`.

**The failure.** CoRAL 3.1.1 infers amplicon intervals from runs of contiguous sequence edges
and makes the two outer nodes of each run *end nodes*
(`coral/breakpoint/breakpoint_graph.py::BreakpointGraph.infer_amplicon_intervals`). At an end node
$v$, its Eulerian constraint (`coral/models/constraints.py::get_eulerian_node_constraints`)
[CoRAL 3.1.1 code](references.md#coral_code) is

```math
x_{(s, v)} + x_{(v, t)} = x_{\mathrm{seq}(v)}
```

where $x$ is the cycle's edge multiplicity, $s, t$ are CoRAL's source and sink, and
$\mathrm{seq}(v)$ is the sequence edge at $v$. Discordant and concordant edges at an end node are
left out of the balance, so a cycle through a breakpoint at an interval end is infeasible.
CoRAL's own graphs extend each amplicon beyond its breakpoints. Our graphs, the truth mixtures
and `aa_export subgraph`, end at breakpoints in 69 of 75 cases. On them CoRAL returned
zero-weight paths and explained 0 copies (RCK `results/graph_route/NEXT_STEPS.md`; RCK
`CHANGELOG.md`).

**The fix.** At both ends of every run, a flank sequence edge of length $F$ is added (`--flank`,
default 10 000 bp), joined to the run by a concordant edge:

```math
\mathrm{cn}(\text{flank at } v) = \max\Big(0,\ \mathrm{cn}(\mathrm{seq}(v)) - \sum_{d \ni v} \mathrm{cn}(d)\Big)
```

where the sum runs over the discordant edges with an end at $v$, a self-loop counted twice (as
CoRAL lists it twice). Coverage and read counts of the new edges are 0. A flank takes at most
$\lfloor(\text{gap} - 2)/2\rfloor$ bp of the gap to the next run on the chromosome, so two flanks
never touch, which would make CoRAL merge them. Sequence and concordant edges are written in
CoRAL's sort order, and path-constraint `e`/`c` indices are renumbered. Discordant edges and their
order are unchanged. `#` and `interval` lines are not copied.

**Effect.** With CoRAL's own parser and model offline (balance, copy number and multiplicity, not
the spanning-tree constraints), truth cycles feasible in CoRAL's greedy model go from 28 to 66 of
67. The remaining one (`nest_samehap_6to1_s804`, CN 40.8) exceeds CoRAL's own read-clustered
discordant-multiplicity bound. Breakpoints on end nodes go from 157 to 0. Slurm test 32104711
recovered the demo5 truth cycle at CN 16 (RCK `results/graph_route/NEXT_STEPS.md`). The graph route
gives the flanked graph to CoRAL only: CycleExtractor joins s/t at every node, and RCK-ecDNA has
its own region-flank joins and flank background, so the fix gives CoRAL the same picture, not an
advantage.

---

## 6. Benchmark-level summaries built on these metrics

The λ choice rule and the failure-stage rule of `results/bench_nested/` are defined in that
README. Two summaries defined on the composite are used across documents.

### 6.1 Mean composite of a benchmark

> **Code:** `results/bench_nested/collect.py`.

This is the mean over cases of each case's `mean_composite` (§3.8) at the reported λ.

### 6.2 Set metric (rank 1 plus its first $k$ same-fit alternatives)

> **Provenance:** (c) new, RCK yubi_dev@de5a805. **Code:** `results/bench_nested/diag/grow_rules.py::main`
> (`set k`).

```math
\mathrm{set}_k = \frac{1}{\lvert \mathcal K \rvert} \sum_{\kappa \in \mathcal K} \max\Big( \bar c_\kappa(D_1),\ \max_{j \le k} \bar c_\kappa(A_{\kappa, j}) \Big)
```

where

- $\mathcal K$ is the cases;
- $\bar c_\kappa(D)$ is case $\kappa$'s mean composite for decomposition $D$;
- $A_{\kappa, j}$ is rank 1's $j$-th same-fit alternative (`pool.json` `same_fit`), matched to its
  pool entry by objective, steps and significant weights; unmatched alternatives are skipped.

**Why.** The grown bench's failure analysis found the right answer among rank 1's same-fit
alternatives in several failing cases (s703, s301, and others), and every re-ranking rule tried
broke as many cases as it fixed. **Effect** at λ = 0.2 (`results/bench_nested/NEXT_STEPS.md`,
"Grown bench: failure analysis"): rank-1 mean 0.8531; $\mathrm{set}_1 = 0.9113$,
$\mathrm{set}_2 = 0.9229$, $\mathrm{set}_5 = 0.9373$; pool best 0.9470. The set metric changes no
case downward. It is not a ranking: it scores what an ambiguity report would let a user see.

---

## 7. Known limitations

- **`aa_format.read_cycles` reads `Cycle=` lines only.** CoRAL 3.1.1 and CycleExtractor write
  linear walks as `Path=` lines, which the scorer skips unless they are rewritten as `Cycle=`
  (`results/graph_route/score_h2h.py` does this; its `_cycles_only` tables are the sensitivity
  check).
- **`aa_eval.read_paths` cannot parse the subgraph's `imputed_ext_*` extremity ids**, so
  `subwalks_satisfied` was not computed for the graph-route arms (RCK `results/graph_route/NEXT_STEPS.md`).
- **`cn-frac:0.1` is not exactly `rck-ecdna`'s significance rule.** `rck-ecdna` flags a species
  `low_weight` below $\max(1,\ 0.1 \times \text{heaviest})$ on *integer* weights
  (`rck/ecdna/decompose.py`, threshold line). `cn-frac:0.1` has no 1-copy floor and reads
  `Copy_count`, which for `rck-ecdna` is the refitted weight. `rck-ecdna` answers use `status`,
  so this affects only AA and CoRAL answers.
- **The composite is equal-weight and copy-number-capped**: see §3.7 (`nest_cross_3sp_s403`).
- **The phased matching needs haplotype overlap.** An unphased answer matches nothing in the
  phased columns, which is why `composite_h0` and the blind matching exist.
- **CoRAL metrics are per truth species, not per amplicon.** CoRAL's code scores all true cycles
  of an amplicon against each reconstructed cycle. Here each truth species is scored against its
  matched cycle (strict) or the best candidate (best-of).
- **The log2 length error and the supplementary-methods overlap and LCS** follow the CoRAL
  supplement as recorded in the previous version of this document. The supplement itself is not
  in either repository (amplicon-allele-analysis `data/COLO320/coral_paper_supplement/` holds its
  README only).
- **Truth alleles exist for the first tumor sample only** (§4.2).
- **Unobserved concordant joins** (at imputed grid points, 10 in demo7) get the smaller
  neighbouring segment CN in the exported AA graph; RCK-ecDNA has no observation there.
- **Truth and inferred read paths differ.** The truth graph gets every read of `obs/path.tsv`
  whose junctions are mixture edges (demo7: 12). The subgraph's `paths.tsv` was filtered by
  `rck_inputs.clip_path` (demo7: 3). Synthetic truth reads are `round(CN × R)`, so the fitted
  junction scale is close to, not exactly, R (demo7: 6.497 for 6.5).

---

## Appendix A. How to run

`rck-ecdna-eval` on a finished run:

```bash
rck-ecdna-eval --decompositions rck_ecdna_out/demo3/amplicon_0/decompositions.tsv \
               --truth-dir simulations/runs/demo3_clone1_ecdna
```

Or pass `--truth-dir` (and optionally `--truth-sample`) to `rck-ecdna` itself, which evaluates
at the end of the run. Options:

| Option | Default | Meaning |
|---|---|---|
| `--decompositions` | required | `decompositions.tsv` written by `rck-ecdna` |
| `--truth-dir` | | simulator run directory (with `truth/` and `obs/`) |
| `--truth-cycles` | | explicit truth table (`species order chrom start end strand hap weight`, 0-based inclusive) |
| `--truth-sample` | `--sample`, else first tumor | sample in `truth/proportions.tsv` for the truth copy numbers |
| `--tol-bp` (`rck-ecdna --eval-tol-bp`) | 10000 | species tolerance $\tau_s$ |
| `--include-low-weight` | off | also score species flagged `low_weight` |
| `--min-truth-fraction` | 0.5 | scope threshold $\phi$ |
| `--all-truth-species` | off | score every truth species |
| `--outdir` | the directory of `--decompositions` | |

The method-agnostic scorer:

```bash
python -m rck.ecdna.aa_eval --cycles METHOD_cycles.txt --graph METHOD_graph.txt --truth-dir SIM_RUN \
    [--alt-cycles ALT_cycles.txt ...] [--segment-cn cn.tsv] [--paths SUBGRAPH/paths.tsv] \
    [--strict-rule auto] [--graph-tol 100] [--species-tol 10000] --outdir OUT
```

| Option | Default | Meaning |
|---|---|---|
| `--cycles` | required | the method's AA cycles file (its answer) |
| `--alt-cycles` | none | other answers of the method, for best-of only |
| `--graph` | | the method's AA graph (graph-level metrics, `cn_explained`) |
| `--segment-cn` | | allele copy numbers per segment: `chrom start end [cn] cn_a cn_b`, AA coordinates |
| `--paths` | | long-read paths (subgraph `paths.tsv`), for `subwalks_satisfied` |
| `--truth-dir`, `--truth-cycles`, `--truth-sample` | | as `rck-ecdna-eval` |
| `--mixture-dir` | `TRUTH_DIR/truth/aa` when it holds `mixture_amplicons.tsv` | truth mixture files; without them graph metrics are skipped with a warning |
| `--strict-rule` | `auto` | `auto`, `status`, `all`, `cn-frac:F`, `top-cn:X` |
| `--species-tol` | 10000 | $\tau_s$ |
| `--graph-tol` | 100 | $\tau_g$ |
| `--min-truth-fraction`, `--all-truth-species` | 0.5, off | scope |
| `--cn-explained-min-cn` | 5 | $\gamma$ |
| `--outdir` | required | |

Building the truth mixture (amplicon-allele-analysis):

```bash
python simulations/truth_mixture.py RUN_DIR [--sample tumor1] [--outdir DIR] \
    [--reads-per-copy R] [--depth-per-copy D] [--path-constraints]
```

Exporting an `rck-ecdna` answer for the scorer:

```bash
python -m rck.ecdna.aa_export RUN_DIR OUT [--pool-index I | --lambda L] [--alternatives 2]
```

The graph route translators:

```bash
# inferred subgraph -> AA graph (+ _segment_cn.tsv, _export.json)
python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name N] [--amplicon 1] \
    [--imputed comment|drop] [--min-path-junctions 1]
# AA graph -> subgraph folder for rck-ecdna --inputs subgraph
python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate] \
    [--background-cn auto|X] [--cn-profile FILE] [--sample NAME]
# flank sequence edges so CoRAL 3.1.1 can use breakpoints at interval ends
python -m rck.ecdna.coral_flank GRAPH OUT_GRAPH [--flank 10000]
```

The bench re-scoring through the AA path: `python results/bench_nested/collect.py --aa-metrics`
(and `results/bench_nested/diag/aa_regression.py`). Head-to-head scoring of the graph route:
`results/graph_route/score_h2h.py`, with commands and caveats in `results/graph_route/README.md`.

## Appendix B. Output files and columns

**`evaluation.tsv`**: one row per (decomposition rank, truth species). It has every field of
`PairScore` (`jaccard_hap`, `jaccard_total`, `segment_*`, `junction_*`, `phased_junction_*`,
`phasing_accuracy*`, `order_match`, `order_match_ordering`, `truth_weight`, `inferred_weight`,
`cn_rel_err`, `inferred_weight_int`, `cn_rel_err_int`, `inferred_topology`, `topology_match`,
`composite`), plus `decomposition_rank`, `objective`, `n_species_inferred`, `n_low_weight` and
`n_species_truth`.

**`evaluation.json`**: `n_decompositions`, `n_truth_species`, `n_truth_species_total`,
`truth_species` (per species: `weight`, `length_bp`, `n_steps`, `n_junctions`,
`fraction_in_subgraph`, `in_scope`), `first_rank_with_all_orders_matched`,
`best_mean_jaccard_hap`, `best_mean_jaccard_rank`, `top_decomposition`, `tol_bp`,
`min_truth_fraction`, `mean_composite`, `recovery_matrix`, and `note` when no species is in
scope. It is strict JSON: missing scores are `null`.

**`recovery_matrix.tsv`**: one row per truth species.

| Column | Meaning |
|---|---|
| `truth`, `in_scope`, `fraction_in_subgraph` | truth species and whether it is scored here (out-of-scope rows have no scores) |
| `truth_length_bp`, `truth_n_junctions`, `truth_weight` | its size, its distinct junctions and its copies per cell |
| `inferred`, `inferred_topology`, `inferred_weight` | the matched rank-1 species (empty when missed) |
| `segment_recall`, `segment_precision`, `segment_jaccard` | §3.1 |
| `haplotype_jaccard` | §3.1 |
| `junction_recall`, `junction_precision`, `junction_f1` | §3.3, unphased |
| `phased_junction_recall`, `phased_junction_precision`, `phased_junction_f1` | §3.3, phased |
| `phasing_accuracy`, `phasing_accuracy_swap_invariant` | §3.2 |
| `order_match`, `topology_match` | §3.4 |
| `cn_rel_err` | §3.5 |
| `composite` | §3.7; 0 for a missed species |
| `best_rank`, `best_inferred`, `best_composite` | §3.8 |

**`truth_cycles.tsv`**: the rebuilt truth (`species order chrom start end strand hap weight`).

**`aa_graph_metrics.tsv`**: one row per scored truth amplicon:

- `amplicon`, `species`, `intervals`, `graph_tol_bp`;
- `n_true_segments`, `n_segments_found`, `n_inferred_segments`, `n_false_segments`;
- `n_true_edges`, `n_edges_found`, `n_inferred_edges`, `n_false_edges`;
- `error_free`, the `*_merged` counts, `error_free_merged`;
- `bp_accuracy`, `novel_adjacency_fnr`, `novel_adjacency_fdr`;
- `is_violations`, `is_violations_truth`;
- `cn_error_mode`, `segment_cn_mae`, `segment_cn_rel_err`, `segment_cn_mae_total`.

**`aa_species_metrics.tsv`**: one row per (mode, truth species):

- `mode`, `truth`, `truth_weight`, `truth_length_bp`, `inferred`, `inferred_weight`,
  `inferred_topology`;
- the §3 columns and `exact_cycle_phased`;
- `inferred_h0`, `segment_jaccard_h0`, `junction_f1_h0`, `order_match_h0`, `cn_rel_err_h0`,
  `composite_h0`;
- `cycle_overlap_coral_code`, `cyclic_lcs_coral_code`, `cycle_overlap_supp`, `cyclic_lcs_supp`,
  `length_error`, `exact_cycle`.

The supplementary-methods versions were named `cycle_overlap` and `cyclic_lcs` before 2026-11-03.

**`aa_sample_metrics.tsv`**: one row per mode: `mode`, `strict_rule`, `phased`, `n_truth`,
`n_inferred`, and the columns of §4.5.

**`aa_eval.json`**: `info` (inputs, rule, tolerances, `truth_in_scope`, `regions`,
`mixture_graphs`), `sample`, `graph`.

## Appendix C. Tests

| Test (RCK `tests/ecdna/`) | Covers | Gurobi |
|---|---|---|
| `test_evaluate.py` | truth reading, several-phasing assignment, species metrics | no |
| `test_recovery.py` | recovery matrix, composite, `best_*` | no |
| `test_aa_eval.py` | AA parsing, graph metrics incl. `error_free_merged` and `merge_benign_splits`, CoRAL metrics (code and supplement), modes, sample metrics | no |
| `test_aa_import.py` | subgraph export, import, scale fit, round trip; reference tokens written by CoRAL in `tests/ecdna/data/graph_route/` | no |
| `test_coral_flank.py` | flank lengths, residual copy numbers, renumbered constraints | no |

amplicon-allele-analysis `simulations/tests/test_truth_mixture.py` covers the truth mixture
(demo6 two amplicons; demo7 nested amplicon, node balance, silent boundary merged). Run named
modules only (`python -m unittest tests.ecdna.test_aa_eval ...`). `unittest discover` also runs
the Gurobi tests.

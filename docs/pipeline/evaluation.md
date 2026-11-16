# Evaluation: scorers, metrics and AA-format translators

This document defines every number the pipeline uses to judge a reconstruction against the
simulator's truth, and the translators that put every method's output into one format. It
covers stage 4 of the pipeline (see [the overview](README.md)):

1. `rck-ecdna-eval` (`rck/ecdna/evaluate.py`): `rck-ecdna`'s own scorer. It produces
   `evaluation.tsv`, `evaluation.json` and the recovery matrix with its composite score. The same
   code runs inside `rck-ecdna` when it is given `--truth-dir` or `--truth-cycles`.
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
also names the code. The terms follow RCK [Aganezov & Raphael 2020](references.md#aganezov2020)
and `docs/ecDNA.md` (the [glossary](glossary.md) maps the code's names):

- *segment*: AA's and CoRAL's "sequence edge";
- *extremity*: a segment end, AA's "breakend" `chr:pos±`;
- *novel adjacency*: the code's "junction", AA's and CoRAL's "discordant edge";
- *reference adjacency*: the code's "reference join", AA's and CoRAL's "concordant edge";
- *haplotype labelling* (phasing) of an adjacency: AA, AB, BA or BB.

A *species* is one distinct ecDNA molecule type. Its *walk* (AA and CoRAL: a *cycle*, or a *path*
when it is linear) is the ordered list of oriented, haplotype-labelled segments it traverses.

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
> `::Cycle`, `::runs`, `::_contiguous`, `::junctions`, `::_unique`, `::_be_match`,
> `::_junction_match`. Where each object comes from:
>
> - **Step $s$ and weight $w_C$**: (a) one segment of a cycle in AmpliconArchitect's cycles
>   file (chromosome, interval, orientation) and the cycle's `Copy_count`
>   [Deshpande et al. 2019](references.md#deshpande2019). CoRAL writes the same format
>   [Zhu et al. 2024](references.md#zhu2024).
> - **Haplotype $H$**: (a) RCK's haplotype labels A and B of a segment, whose adjacencies are
>   labelled AA, AB, BA or BB [Aganezov & Raphael 2020](references.md#aganezov2020). AA and CoRAL
>   have no haplotype field. `.` (unphased) is (c).
> - **Extremity and join**: (a) RCK's segment extremity and adjacency. A join that is not
>   reference-contiguous is RCK's *novel adjacency*, the *discordant edge* of AA and CoRAL
>   [Aganezov & Raphael 2020](references.md#aganezov2020) [Deshpande et al. 2019](references.md#deshpande2019).
> - **Per-base multiplicity $m_q$, reference runs and the 2 bp contiguity slack**: (c). No tool
>   defines them. The code records no reason for 2 bp. The slack is wide enough to absorb the 1 bp
>   coordinate-convention differences of §5.1.
> - **Join match within $\tau$**: (a) CoRAL 3.1.1's `is_breakpoint_match`, which tries both
>   directions of `is_oriented_breakpoint_match` (same chromosome and strand at both ends,
>   positions within the tolerance), with $\le$ instead of $<$ (box below). $\tau_{\mathrm{gr}}$ =
>   100 bp is (a), the AA and CoRAL breakpoint tolerance (CoRAL's default `tolerance`). $\tau_{\mathrm{sp}}$ = 10 kb is (c), set by
>   the subgraph's 10 kb imputation grid.

**Steps.** A step (code `CStep`) is
$s = (\mathrm{chr}_s,\ \mathrm{start}_s,\ \mathrm{end}_s,\ H_s,\ \mathrm{str}_s)$:

- $\mathrm{chr}_s$ is the chromosome;
- $[\mathrm{start}_s, \mathrm{end}_s]$ is a 0-based inclusive interval of integers;
- $H_s \in \lbrace \mathrm{A}, \mathrm{B}, \cdot \rbrace$ is the haplotype, where $\cdot$ (code `.`)
  means unphased;
- $\mathrm{str}_s \in \lbrace +, - \rbrace$ is the orientation (strand) of the traversal.

Its length is $L(s) = \mathrm{end}_s - \mathrm{start}_s + 1$, RCK's segment length $L(j)$.

**Walks.** A walk (code `Cycle`, whatever its topology) $q = (s_1, \dots, s_n)$ has:

- a topology $\mathrm{circ}(q) \in \lbrace 0, 1 \rbrace$ (code `circular`): 1 for a circle, whose
  last step joins its first; 0 for a linear walk, which has no such wrap-around join. Truth walks
  are always circular (§2.1);
- a weight (code `weight`), copies per cell: $w_T$ for a truth walk $T$ (§2.1), $w_C$ for an
  inferred walk $C$ (§3.5); an inferred walk may also have an integer weight $w^{\mathrm{int}}_C$;
- one or more *orderings* $q^{(1)}, \dots, q^{(K_q)}$ (code `orderings`). A truth walk has one.
  For an inferred species, `rck-ecdna` reports up to `--max-orderings` (default 20) distinct
  Euler orderings of the species' edges (`docs/ecDNA.md` §7); $q = q^{(1)}$ is the top one, and
  every score except order match (§3.4) uses only it.

The length of a walk counts multiplicity: $\lvert q \rvert = \sum_{i=1}^{n} L(s_i)$.

Inferred walks are read from `decompositions.tsv` (`read_inferred`): one row per (decomposition
rank, species, ordering); `walk` lists the steps as `chr:start-end:hap:strand` (0-based inclusive),
rows are ordered by `ordering_rank`, `topology` = `linear` gives $\mathrm{circ} = 0$, and `weight`,
`weight_int` give $w_C$, $w^{\mathrm{int}}_C$.

**Multiplicity.** For a base at position $\mathrm{pos}$ on chromosome $\mathrm{chr}$ and haplotype
$H$, with $q = (s_1, \dots, s_n)$:

```math
m_q(\mathrm{chr}, \mathrm{pos}, H) = \#\lbrace i : \mathrm{chr}_{s_i} = \mathrm{chr},\ H_{s_i} = H,\ \mathrm{start}_{s_i} \le \mathrm{pos} \le \mathrm{end}_{s_i} \rbrace,
\qquad
m_q(\mathrm{chr}, \mathrm{pos}) = \sum_{H} m_q(\mathrm{chr}, \mathrm{pos}, H).
```

The sum runs over $H \in \lbrace \mathrm{A}, \mathrm{B}, \cdot \rbrace$, so $\cdot$ is a third
label, not a wildcard.

**Reference contiguity.** Consecutive steps $s, s'$ are contiguous (`_contiguous`) when they agree
in $\mathrm{chr}$, $H$ and $\mathrm{str}$ and the second continues the first within 2 bp:

- $\lvert \mathrm{start}_{s'} - (\mathrm{end}_s + 1) \rvert \le 2$ for $\mathrm{str} = +$;
- $\lvert \mathrm{start}_s - (\mathrm{end}_{s'} + 1) \rvert \le 2$ for $\mathrm{str} = -$.

`runs(q)` merges maximal contiguous stretches into one step with the union span
$[\min \mathrm{start}, \max \mathrm{end}]$. For a circle the scan starts after the first
non-contiguous pair, so a stretch that wraps around is one run; a circle whose every pair is
contiguous becomes one run. A different segmentation of the same molecule therefore gives the same
runs. Because contiguity requires the same haplotype, a haplotype switch between reference
neighbours is a join, not a continuation.

**Extremities.** An extremity is $v = (\mathrm{chr}, \mathrm{pos}, \pm)$ (code `Breakend`):
$+$ is the right end of a segment (RCK's head $j_h$, at $\mathrm{end}$), $-$ its left end (RCK's
tail $j_t$, at $\mathrm{start}$), as in AA's `chr:pos±`. A step leaves by its exit extremity and
arrives by its entry extremity:

```math
v^{\mathrm{out}}_s = \begin{cases} (\mathrm{chr}_s, \mathrm{end}_s, +) & \mathrm{str}_s = + \\ (\mathrm{chr}_s, \mathrm{start}_s, -) & \mathrm{str}_s = - \end{cases}
\qquad
v^{\mathrm{in}}_s = \begin{cases} (\mathrm{chr}_s, \mathrm{start}_s, -) & \mathrm{str}_s = + \\ (\mathrm{chr}_s, \mathrm{end}_s, +) & \mathrm{str}_s = - \end{cases}
```

**Joins (novel adjacencies of a walk).** `junctions(q)` lists, in walk order, every
non-contiguous consecutive pair $(s, s')$ as a join

```math
a = \big(v^{\mathrm{out}}_s,\ v^{\mathrm{in}}_{s'},\ H_s,\ H_{s'}\big).
```

For a circle the last step joins the first; a linear walk has no such wrap-around join. Joins are
taken on the steps, with the 2 bp slack; on runs they would be the same.
$\mathcal A_N(q)$ is the list of *distinct* joins (`_unique`): a join is dropped when an earlier one
is identical, in the same direction, or reversed with its haplotype pair swapped. Distinctness is
phased: one novel adjacency traversed once as AA and once as BB gives two elements of
$\mathcal A_N(q)$, for the unphased scores too.

**Matching.** Two extremities match within $\tau$, $v \simeq_\tau v'$, when they have the same
chromosome, the same sign and $\lvert \mathrm{pos} - \mathrm{pos}' \rvert \le \tau$ (`_be_match`).
For joins $a = (v_1, v_2, H_1, H_2)$ and $a' = (v'_1, v'_2, H'_1, H'_2)$ (`_junction_match`):

- $a \sim_\tau a'$ (unphased) when $v_1 \simeq_\tau v'_1 \wedge v_2 \simeq_\tau v'_2$ (same
  direction) or $v_1 \simeq_\tau v'_2 \wedge v_2 \simeq_\tau v'_1$ (reversed);
- $a \sim^{\mathrm{ph}}_\tau a'$ (phased) also requires $(H_1, H_2) = (H'_1, H'_2)$ in the same
  direction, or $(H_1, H_2) = (H'_2, H'_1)$ reversed. The same direction is tried first: when its
  extremities match, its haplotype test decides, and the reversed direction is not tried.

Tolerances used below:

| Symbol | Option | Default | Used for |
|---|---|---|---|
| $\tau_{\mathrm{sp}}$ | `rck-ecdna --eval-tol-bp`, `rck-ecdna-eval --tol-bp`, `aa_eval --species-tol` | 10 000 bp | species metrics and the scope padding (§2.3) (the subgraph's imputed breakpoints sit on a 10 kb grid) |
| $\tau_{\mathrm{gr}}$ | `aa_eval --graph-tol` | 100 bp | graph metrics, `exact_cycle`, subwalks (the AA and CoRAL convention) |

> **Inclusive tolerance.** All matches here use $\le \tau$. CoRAL 3.1.1 uses a strict $< \tau$
> (`coral/scoring/scoring_utils.py::is_oriented_breakpoint_match`) [CoRAL 3.1.1 code](references.md#coral_code).

**Relation to the DIAG of `docs/ecDNA.md`.** An inferred walk $C$ is one walk $q$ of the
decomposition (`docs/ecDNA.md` §5 and §7; $q \in \mathcal W$ in §11.3), written out by coordinates.
A step $s$ is one traversal of a segment edge $j^H$, $H = H_s$: from $j^H_t$ to $j^H_h$ when
$\mathrm{str}_s = +$, and from $j^H_h$ to $j^H_t$ when $\mathrm{str}_s = -$; so
$v^{\mathrm{out}}_s$ is $j_h$ or $j_t$ without its haplotype, and the haplotype travels in the
join. A contiguous pair that `runs` merges is a reference adjacency of $E_R$ labelled AA or BB. A
join $a$ is a labelled adjacency $a^\varphi$ with $\varphi = H_s H_{s'}$, read in walk order
(RCK orders $\varphi$ by the adjacency's first end, so a join read backwards has $\varphi$
reversed). It is usually a novel adjacency of $E_N$; it can also be a reference adjacency
traversed with a haplotype switch (AB or BA), which `docs/ecDNA.md` §2.5 does not allow, or a
reference adjacency across a gap of more than 2 bp. Truth walks $T$ are not walks of the run's
DIAG: their steps are the simulator's intervals, which need not coincide with the DIAG's
segments. Every comparison is therefore by coordinates within a tolerance, never by edge
identity. An unphased step ($H_s = \cdot$, AA and CoRAL answers) has no DIAG counterpart.

**Symbols of §1–§3.**

| Symbol | Meaning | Type / domain | Counterpart (RCK / `docs/ecDNA.md` / AA / CoRAL) |
|---|---|---|---|
| $s = (\mathrm{chr}_s, \mathrm{start}_s, \mathrm{end}_s, H_s, \mathrm{str}_s)$ | step (`CStep`) | chromosome name, integers (0-based inclusive), $\lbrace \mathrm{A}, \mathrm{B}, \cdot \rbrace$, $\lbrace +, - \rbrace$ | traversal of segment edge $j^H$; AA/CoRAL cycle token `k±` of a sequence edge |
| $L(s)$ | step length, $\mathrm{end}_s - \mathrm{start}_s + 1$ | bp, $\ge 1$ | RCK $L(j)$ |
| $H$ | haplotype | $\lbrace \mathrm{A}, \mathrm{B}, \cdot \rbrace$ | RCK $H \in \lbrace A, B \rbrace$; none in AA/CoRAL |
| $q$; $T$; $C$ | a walk; a truth walk; an inferred walk (`Cycle`) | sequence of steps | ecDNA.md walk $q \in \mathcal W$; AA/CoRAL cycle (circular) or path (linear) |
| $\mathrm{circ}(q)$ | topology | $\lbrace 0, 1 \rbrace$ | AA/CoRAL cycle vs path (source edges) |
| $q^{(\ell)}$, $K_q$ | $\ell$-th ordering of $q$, number of orderings | $1 \le \ell \le K_q \le$ `--max-orderings` | ecDNA.md §7 orderings |
| $\lvert q \rvert$ | walk length with multiplicity | bp | AA/CoRAL cycle length |
| $w_T$, $w_C$, $w^{\mathrm{int}}_C$ | truth weight, inferred weight after refit, integer weight before refit | copies per cell, $\ge 0$ | AA/CoRAL `Copy_count`; ecDNA.md §6 refit, §5 slot weight $w_k$; multi-clone $\sum_i u_i \psi_{q,i}$ |
| $m_q(\mathrm{chr}, \mathrm{pos}, H)$, $m_q(\mathrm{chr}, \mathrm{pos})$ | per-base multiplicity, with and without haplotype | $\mathbb N$ | none (new) |
| $\mathrm{runs}(q)$ | $q$ with contiguous steps merged | sequence of steps | none (new) |
| $v = (\mathrm{chr}, \mathrm{pos}, \pm)$ | extremity (`Breakend`) | $+$ right end, $-$ left end | RCK $j_\sigma$, $\sigma \in \lbrace t, h \rbrace$ ($-$ = $t$, $+$ = $h$); AA breakend `chr:pos±` |
| $v^{\mathrm{out}}_s$, $v^{\mathrm{in}}_s$ | exit and entry extremity of a step | extremity | ecDNA.md $j_h$ / $j_t$ by strand |
| $a = (v_1, v_2, H_1, H_2)$ | join between consecutive non-contiguous steps | extremity pair and haplotype pair | ecDNA.md $a^\varphi$, $\varphi = H_1 H_2$; AA/CoRAL discordant edge |
| $\mathcal A_N(q)$ | distinct joins of $q$ | list of joins | RCK $A_N$ restricted to $q$ |
| $\simeq_\tau$, $\sim_\tau$, $\sim^{\mathrm{ph}}_\tau$ | extremity match, join match, phased join match | relations | CoRAL `is_oriented_breakpoint_match`, `is_breakpoint_match` (strict $<$) |
| $\tau_{\mathrm{sp}}$, $\tau_{\mathrm{gr}}$ | species and graph tolerance | bp; 10 000, 100 | CoRAL tolerance 100 |
| $\mathcal W^\ast$, $\mathcal W^\ast_{\mathrm{in}}$ | truth species, in-scope truth species | sets of walks | none |
| $\mathcal Q$, $t$, $\Phi_t$, $\Phi^{\mathrm{listed}}_a$, $\boldsymbol\varphi$ | truth traversals with several options, one traversal, its options, a junction's listed phasings, an assignment | §2.2 | ecDNA.md $\Phi(a)$, $\varphi$ |
| $H_0$ | default haplotype of a truth species | $\lbrace \mathrm{A}, \mathrm{B} \rbrace$ | none |
| $\mathrm{Reg}$, $(\mathrm{chr}, r_0, r_1)$ | regions of a run, one region | 0-based inclusive spans | ecDNA.md §2.6 regions |
| $f(T)$, $f_{\min}$ | fraction of $T$ inside the padded regions, its threshold | $[0, 1]$ | none |
| $O(T, C)$, $O_H(T, C)$ | base-pair overlap without and with haplotype | bp | none |
| $\bar C$ | $C$ with A and B exchanged | walk | none |
| $\mathrm{Rot}$, $\mathrm{rot}$ | rotations and reversals of a run sequence, one of them | set of maps | CoRAL cyclic LCS search space |
| $\approx_\tau$ | run sequences equal within $\tau$ | relation | none |
| $r$, $\mathcal W_r$ | decomposition rank, its inferred walks | $r \ge 1$; set of walks | ecDNA.md §5.5, §5.7 ranked decompositions |
| $\mathrm{match}_r(T)$ | $T$'s partner in $\mathcal W_r$ (§3.6) | walk or $\varnothing$ | none |
| $J_{\mathrm{seg}}$, $F_1^{\mathrm{junc}}$, $\mathrm{PA}$, $\mathrm{err}_{\mathrm{cn}}$ | composite terms (§3.7) | $[0,1]$, $[0,1]$, $[0,1]$, $\ge 0$ | none |
| $i$, $\mathrm{Cl}$, $u_i$ | clone, clones scored, clone proportion | clone name; $u_i \in [0, 1]$ | RCK clone $i$, $u_i$; ecDNA.md §11 |
| $\psi_{q,i}$, $\psi^\ast_{T,i}$ | inferred and truth copies of a walk in clone $i$ | $\mathbb N$ | ecDNA.md §11.3 $\psi_{q,i}$ |
| $\mathcal B$, $b$, $\mathrm{sw}_b$ | phase blocks, one block, its A/B swap | spans; $\lbrace 0, 1 \rbrace$ | ecDNA.md §3.6 block $b$, §11.4 fragment $f$ |
| $q\vert_b$ | $q$'s steps clipped to block $b$ | sequence of steps | none |

**Symbols of §4–§6.** Defined where first used; listed here so that none is read with its
`docs/ecDNA.md` meaning.

| Symbol | Meaning | Section | Counterpart / clash avoided |
|---|---|---|---|
| $\mathrm{Amp}$, $I_{\mathrm{Amp}}$ | truth amplicon (connected group of species), its intervals | §4.2 | AA amplicon; not haplotype $A$ |
| $\mu^\ast(j)$, $\mu^\ast(j^H)$, $\mu^\ast(a)$ | truth mixture copy number of segment $j$ (total, haplotype $H$) and of adjacency $a$ | §4.2 | RCK $\mu_i$ summed over clones $i$ with weights $u_i$ |
| $\hat\mu$, $\hat\mu_A$, $\hat\mu_B$ | a method's segment copy number (total, allele) | §4.3 | AA `PredictedCopyCount` |
| $\mathrm{cov}_T(j)$, $\mathrm{trav}_T(a)$ | times truth walk $T$ covers segment $j$ / traverses adjacency $a$ | §4.2 | |
| $\mathrm{disc}(v)$ | summed copy number of discordant edges at extremity $v$ | §4.2 | not the DIAG $D$ |
| $\mathcal A^\ast_N$, $\hat{\mathcal A}_N$ | truth and answer novel adjacencies (discordant edges) of a graph | §4.3 | RCK $A_N$ |
| $\mathrm{Seg}^\ast$, $\widehat{\mathrm{Seg}}$, $\mathrm{merge}(G)$ | truth and answer segments; graph $G$ with benign splits merged | §4.3 | AA sequence edges |
| $T^0$, $C^0$ | haplotype-stripped walks (`inferred_h0` is the partner of $T^0$) | §4.4 | |
| $\mathrm{Bins}(T,C)$, $\Delta$, $\mathrm{Bins}_\cap$ | joint bins of two walks, one bin, bins both cover | §4.4 | CoRAL `repartition_intervals` |
| $\mathrm{Cov}(X)$ | bases covered by walk $X$, without multiplicity | §4.4 | |
| $\mathrm{seq}_X$, $\mathrm{WLCS}$, $\mathrm{LCS}^\ast$ | bin-label sequence, weighted LCS, its best value over $\mathrm{Rot}$ | §4.4 | CoRAL `find_longest_subsequence` |
| $R_{\mathrm{len}}$ | length ratio | §4.4 | CoRAL length ratio; not reference $R$ |
| $\mathcal S$, $\mathcal C(T)$, $\mathcal X$ | strict answer, best-of candidates of $T$, the mode's walks | §4.5–§4.6 | |
| $\mathcal E_{\mathrm{seq}}$, $\mathrm{CN}(j)$, $\mathrm{cn}_{\min}$ | sequence edges of a graph, their copy number, CoRAL's minimum (5) | §4.5 | CoRAL `get_cycle_copy_number_ratio` |
| $\upsilon$, $\Upsilon$, $N_\upsilon$ | long-read path, set of paths, support | §4.5 | ecDNA.md §5.4 |
| $f_{\mathrm{cn}}$, $f_{\mathrm{mass}}$ | `--strict-rule` thresholds | §4.6 | not $F_e$ |
| $\mathrm{pass}(T,C)$, $\mathrm{thr}$ | M1 detection criterion and its threshold | §4.7 | not $\kappa_\ast$, $\theta$ |
| $\mathcal I$, $\mathcal I(X)$, $\mathrm{dist}(a,a')$ | M2 item universe, items of walk $X$, adjacency distance | §4.7 | not $\iota$, $U_j$ |
| $\mathrm{tog}^\ast$, $\widehat{\mathrm{tog}}$ | pair together in some truth / answer species | §4.7 | not $y_e$ |
| $J_{\mathrm{seg}}$, $J_{\mathrm{adj}}$, $n_\cap$, $n_T$, $n_C$ | species-matrix Jaccards and adjacency counts | §4.7 | |
| $\theta$, $\gamma$, $\hat\theta$, $\hat\gamma$, $\mathrm{rd}_e$, $\mathrm{cn}_e$ | junction reads and depth per copy (given / estimated), an edge's reads and copy number | §5.2–§5.3 | same as ecDNA.md §3.1–§3.2 |
| $L_{\mathrm{fl}}$ | CoRAL flank length | §5.4 | not $F_e$ |
| $\mathcal K$, $\mathfrak c$ | cases of a benchmark, one case | §6 | not slot $k$ |
| $\overline{\mathrm{comp}}_{\mathfrak c}$, $\mathcal W^{\mathrm{alt}}_{\mathfrak c,l}$, $n_{\mathrm{alt}}$ | case mean composite, $l$-th same-fit alternative, number of alternatives | §6.1–§6.2 | |
| $\mathcal G$, $\mathfrak a$, $\mathrm{rate}$ | pooled group of answers, one answer, a confusion rate | §6.3 | |

---

## 2. Reading the truth

### 2.1 Species walks and coordinates

> **Provenance:** (c) new, RCK yubi_dev@c144584, generalised in yubi_dev@3ba7627. Parsing has
> been shared with `aa_format` since yubi_dev@e4c2409. The AA cycles format itself is (a)
> AmpliconArchitect's [Deshpande et al. 2019](references.md#deshpande2019), as ecSimulator writes it [ecSimulator](references.md#ecsimulator).
> **Code:** `rck/ecdna/evaluate.py::truth_from_simulation`, `::read_aa_cycles`,
> `::_segments_from_profiles`, `::species_weights`, `::resolve_sample`, `::read_truth_table`;
> `rck/ecdna/aa_format.py::raw_walks`.

The truth species $\mathcal W^\ast$ are the species listed in `truth/ecdna.tsv` (column
`SPECIES`), or, without that file, one per `truth/aa/<species>_cycles.txt`. For each:

1. **Walk and coordinates.** They come from `truth/aa/<species>_cycles.txt`, which the
   simulator writes with ecSimulator's own AA writer. This works whatever the species went
   through: deletions, shared chromosomes, several species.
   - AA coordinates are 1-based inclusive. A segment row is read as `START = start − 1`,
     `END = end` (0-based half-open, as `truth/segments.tsv`), and a step as
     $[\mathrm{start} - 1, \mathrm{end} - 1]$.
   - A file with several cycles gives several truth walks, named `<species>.<k>`, each with the
     species' weight. A path (segment `0`) in the file is an error: truth species must be
     circular.
   - The file's haplotype field, if any, is not read. Haplotypes come from §2.2.
   - Older runs without `truth/aa/` fall back to rebuilding the walk from `truth/walks.walk`
     (circular walks only; a species without one is skipped with a warning) and
     `truth/segments.tsv` (`_segments_from_profiles`). There, the species' segments are the rows
     whose `CN_A + CN_B` exceeds the clone-mixture background
     $\sum_i u_i (\mathrm{A}_i + \mathrm{B}_i)$ of `truth/clones.tsv` by more than 0.5. Rows that
     touch form intervals. Among the combinations of intervals with `N_SEGMENTS` rows in total
     (else the walk's largest segment index), the first order, tried over every permutation, under
     which every join of the walk is a reference join or a junction of that species is kept. The
     search gives up after 50 000 orders; the run then fails and asks for `--truth-cycles`.
2. **Copy number.** $w_T$ is the species' copies per cell in the evaluated sample: the
   `ecdna_burden` rows of `truth/proportions.tsv`, column `PHI` (`species_weights`). A species
   without such a row keeps `ecdna.tsv` `PHI`; a species in neither has no weight (NaN).
   The sample (`resolve_sample`) is `--truth-sample` (in `rck-ecdna`, else `--sample`). A name that is not a sample of `proportions.tsv` is mapped to the longest sample
   name that it ends with or contains (`sim_tumor1` → `tumor1`). Without a name, or with no hit
   (a warning), it is the first sample of `proportions.tsv` other than `normal`.

The truth is written to `truth_cycles.tsv` (`write_truth_table`). For real data, write that table
yourself and pass it as `--truth-cycles`: columns `species order chrom start end strand hap
weight`, 0-based inclusive coordinates, one row per step, sorted by `order` within a species
(`read_truth_table`). Every walk of the table is circular.

### 2.2 Truth haplotypes, including junctions with several phasings

> **Provenance:** (c) new, RCK yubi_dev@c144584. The several-phasings rule is (c),
> yubi_dev@184ec95. **Code:** `rck/ecdna/evaluate.py::_label_haplotypes`, `::truth_phasings`,
> `::species_junctions`, `::_majority_hap`, `::_concordant`.

Haplotypes come from the species' junction phasings (`truth/adjacency_clone.tsv` plus
`obs/adjacency.tsv`) and are carried across reference joins with a union-find over the walk's
steps.

`truth/adjacency_clone.tsv` keeps one `PHASING` per junction, its most common. Runs simulated with
per-step haplotype labels (amplicon-allele-analysis `hap_labels: per_step`, every run before that
key) can traverse one junction with several phasings, for example AA on one copy and BB on its
duplicate. `truth/rck/rck.acnt.tsv` lists them all, and `truth_phasings` reads every phasing with
copies in some clone. For a junction $a$ of the species, $\Phi^{\mathrm{listed}}_a$ is that set,
in the order AA, AB, BA, BB, or the one `PHASING` when `rck.acnt.tsv` does not list it. A
phasing $\varphi = H H'$ gives $H$ to the step at the junction's `EID-1` end and $H'$ to the
other.

The labelling, for a walk $(s_1, \dots, s_n)$ with joins $t = 1, \dots, n$ (join $t$ goes from
$s_t$ to $s_{t+1}$, cyclically):

1. If join $t$ is a junction $a$ of the species, its options are
   $\Phi_t = \Phi^{\mathrm{listed}}_a$. If its two extremities are also reference neighbours
   (`_concordant`: same chromosome, opposite signs, the $-$ position is the $+$ position plus 1),
   AA and BB are appended when missing: such a crossing may also be traversed as a plain
   reference join on one haplotype. A join with one option is *fixed*: it labels its two steps.
   A join with several is a traversal $t \in \mathcal Q$.
2. Else, if its extremities are reference neighbours, the two steps are merged in the union-find.
3. Else the join is unexplained: it is counted in a warning and labels nothing.
4. An assignment $\boldsymbol\varphi = (\varphi_t)_{t \in \mathcal Q}$ picks one option per
   traversal. The fixed labels, then the assigned ones, each in walk order, are applied to the
   steps' union-find classes. A class keeps the first label it receives.
   $\mathrm{conflicts}(\boldsymbol\varphi)$ is the number of applied labels that differ from the
   label their class already has.
5. The chosen assignment is the first lexicographic minimum, over the product of the $\Phi_t$
   enumerated in walk order and option order:

```math
\boldsymbol\varphi^\star = \mathop{\mathrm{arg\,lexmin}}\limits_{\boldsymbol\varphi \in \prod_{t \in \mathcal Q} \Phi_t}
\Big( \mathrm{conflicts}(\boldsymbol\varphi),\ -\sum_{a} \big\lvert \lbrace \varphi_t : t \in \mathcal Q \text{ traverses } a,\ \varphi_t \in \Phi^{\mathrm{listed}}_a \rbrace \big\rvert \Big)
```

   The second term counts, summed over junctions $a$, the distinct listed phasings of $a$ that the
   assignment uses. The search stops at the first assignment that reaches
   $(0, -\sum_a \min(n_a, \lvert \Phi^{\mathrm{listed}}_a \rvert))$, where $n_a$ is the number of
   traversals of $a$ in $\mathcal Q$, and in any case after 50 000 assignments
   (`MAX_PHASING_ASSIGNMENTS`, with a warning), keeping the best so far.
6. A step whose class has no label takes the species' default haplotype $H_0$
   (`_majority_hap`): the more frequent letter among the listed phasings of all the species'
   junctions (each junction once, every listed phasing counted); A on a tie or without junctions.

**Why.** Applying one `PHASING` to every traversal mislabelled 0.995 Mb of the 5.37 Mb
bench_nested `ctrl_dup_s702` truth walk. Re-scored, its committed λ = 0 result went from phasing
accuracy 0.815 to 1.0, haplotype Jaccard from 0.687 to 1.0, and composite from 0.744 to 0.781
(RCK `CHANGELOG.md`, Fixed). Of the 22 bench_nested and step-sweep truth walks of that time, only
`ctrl_dup_s702` changed. The simulator's own fix, inherited haplotype labels (amplicon-allele-analysis
`hap_labels: inherit`), removes such junctions from new runs. This rule keeps old runs scorable.

### 2.3 Which truth species are scored (scope)

> **Provenance:** (c) new, RCK yubi_dev@3ba7627. **Code:** `rck/ecdna/evaluate.py::fraction_in_regions`,
> `::evaluate`, `::regions_from_summary`; `rck/ecdna/subgraph_io.py::ECGraph.regions`;
> `rck/ecdna/aa_eval.py::regions_of`, `::score`.

A run covers one amplicon subgraph, but a simulation can hold several species. A truth species
$T$ is in scope, $T \in \mathcal W^\ast_{\mathrm{in}}$, when $f(T) \ge f_{\min}$, with

```math
f(T) = \min\left(1,\ \frac{1}{\lvert T \rvert} \sum_{s \in T}\ \sum_{(\mathrm{chr}_s, r_0, r_1) \in \mathrm{Reg}} \big\lvert [\mathrm{start}_s, \mathrm{end}_s] \cap [r_0 - \tau_{\mathrm{sp}},\ r_1 + \tau_{\mathrm{sp}}] \big\rvert \right)
```

and $f(T) = 0$ when $\lvert T \rvert = 0$. Here

- $\mathrm{Reg}$ is the run's regions, each $(\mathrm{chr}, r_0, r_1)$, 0-based inclusive:
  - in `rck-ecdna` (both routes), the run graph's `ECGraph.regions`: its segments sorted by
    position, with a new region wherever a segment starts more than 1 Mb after the current
    region's end (`docs/ecDNA.md` §2.6). The multi-clone route uses its DIAG graph;
  - in `rck-ecdna-eval`, `summary.json` `graph.regions` in the folder of `--decompositions`
    (the same regions, written by `rck-ecdna`). Without that file, every species is in scope;
  - in `aa_eval`, the answer's `Interval` lines, or else the graph's segments merged across
    gaps up to 1 Mb;
- each region is padded by $\tau_{\mathrm{sp}}$ (the scorer's `--tol-bp` / `--eval-tol-bp`);
- the sum runs over regions, not their union: overlapping padded regions count a base twice,
  and the cap at 1 bounds the result;
- $f_{\min}$ is `--min-truth-fraction` (default 0.5).

Out-of-scope species are listed in `evaluation.json` `truth_species` with `in_scope: false` and
`fraction_in_subgraph` ($f(T)$ rounded to 4 decimals). They are not scored and do not count as
missed. `n_truth_species` counts $\mathcal W^\ast_{\mathrm{in}}$ and `n_truth_species_total`
counts $\mathcal W^\ast$. `--all-truth-species` scores every species; `min_truth_fraction` is then
`null`. When no species is in scope, `evaluation.json` has a `note`.

---

## 3. `rck-ecdna-eval`: species metrics, matching, recovery matrix and composite

### 3.1 Base-pair overlap: segment recall, precision and Jaccard; haplotype Jaccard

> **Provenance:** (c) new, RCK yubi_dev@c144584. A base-pair Jaccard index between reconstructed
> and true intervals is a standard set-similarity score [Jaccard 1912](references.md#jaccard1912); the
> multiplicity-weighted, haplotype-resolved form is ours. **Code:**
> `rck/ecdna/evaluate.py::_overlap_bp`, `::jaccard`, `::score_pair`.

For a truth walk $T$ and an inferred walk $C$ (first orderings):

```math
O(T, C) = \sum_{\mathrm{chr}} \sum_{\mathrm{pos}} \min\big(m_T(\mathrm{chr}, \mathrm{pos}),\ m_C(\mathrm{chr}, \mathrm{pos})\big),
\qquad
O_H(T, C) = \sum_{\mathrm{chr}} \sum_{H} \sum_{\mathrm{pos}} \min\big(m_T(\mathrm{chr}, \mathrm{pos}, H),\ m_C(\mathrm{chr}, \mathrm{pos}, H)\big)
```

```math
\mathrm{segment\_recall} = \frac{O}{\lvert T \rvert},\qquad
\mathrm{segment\_precision} = \frac{O}{\lvert C \rvert},\qquad
\mathrm{segment\_jaccard} = \frac{O}{\lvert T \rvert + \lvert C \rvert - O},\qquad
\mathrm{haplotype\_jaccard} = \frac{O_H}{\lvert T \rvert + \lvert C \rvert - O_H}
```

where

- $O$ counts bases both walks cover, with multiplicity and ignoring haplotype;
- $O_H$ counts them per haplotype $H \in \lbrace \mathrm{A}, \mathrm{B}, \cdot \rbrace$, so an
  unphased step overlaps only unphased steps;
- $\lvert T \rvert$ and $\lvert C \rvert$ are the walk lengths with multiplicity (§1);
- recall (precision) is 0 when $\lvert T \rvert = 0$ ($\lvert C \rvert = 0$), and a Jaccard
  index with an empty union is 0. All four lie in [0, 1].

In `evaluate.py` the segment Jaccard is the field `jaccard_total` and the haplotype Jaccard is
`jaccard_hap`; the recovery matrix renames them (`_RENAMED`). Coverage does not depend on
segmentation, so these scores are the same on steps and on runs.

### 3.2 Phasing accuracy

> **Provenance:** (c) new, RCK yubi_dev@4e53747. **Code:** `rck/ecdna/evaluate.py::score_pair`,
> `::_swap_hap`. $O$ and $O_H$ are §3.1's base-pair overlaps without and with haplotype
> agreement. The A/B labels are RCK's [Aganezov & Raphael 2020](references.md#aganezov2020). Neither
> AA nor CoRAL scores phasing. The swap-invariant form follows the convention of haplotype-assembly
> benchmarks, where the two haplotypes are an unordered pair. Their switch error
> [Edge et al. 2017](references.md#edge2017) counts discordant relative phase between adjacent
> variants, and so is unchanged when A and B are exchanged. Here the exchange is global, one swap
> per species, not one per phase block. That makes it coarser than a switch error.

```math
\mathrm{phasing\_accuracy} = \frac{O_H(T, C)}{O(T, C)},
\qquad
\mathrm{phasing\_accuracy\_swap\_invariant} = \frac{\max\big(O_H(T, C),\ O_H(T, \bar C)\big)}{O(T, C)}
```

where $\bar C$ is $C$ with A and B exchanged on every step ($\cdot$ stays $\cdot$). Both scores are
0 when $O = 0$, and lie in [0, 1]. The swap-invariant form is for real data, whose phase labels
need not match the truth's.

### 3.3 Junction (novel adjacency) recall, precision and F1, unphased and phased

> **Provenance:** (c) new, RCK yubi_dev@c144584. Matching adjacencies by both extremities within a
> tolerance follows CoRAL's breakpoint matching [CoRAL 3.1.1 code](references.md#coral_code)
> (`coral/scoring/scoring_utils.py::is_breakpoint_match`). Here the tolerance is inclusive and the
> phased variant is added. **Code:** `rck/ecdna/evaluate.py::junction_scores`, `::_junction_match`,
> `::_unique`, `::_f1`.

Let $\mathcal A_N(T)$ and $\mathcal A_N(C)$ be the distinct joins of the two walks (§1). For a
match relation $\sim$, either $\sim_{\tau_{\mathrm{sp}}}$ (unphased) or
$\sim^{\mathrm{ph}}_{\tau_{\mathrm{sp}}}$ (phased):

```math
\mathrm{recall} = \frac{\lvert \lbrace a \in \mathcal A_N(T) : \exists\, a' \in \mathcal A_N(C),\ a \sim a' \rbrace \rvert}{\lvert \mathcal A_N(T) \rvert},
\quad
\mathrm{precision} = \frac{\lvert \lbrace a' \in \mathcal A_N(C) : \exists\, a \in \mathcal A_N(T),\ a' \sim a \rbrace \rvert}{\lvert \mathcal A_N(C) \rvert},
\quad
F_1 = \frac{2\,\mathrm{recall}\cdot\mathrm{precision}}{\mathrm{recall} + \mathrm{precision}}
```

Conventions:

- recall is 1 when $\mathcal A_N(T) = \emptyset$, and precision is 1 when $\mathcal A_N(C) = \emptyset$;
- $F_1 = 0$ when recall + precision = 0;
- matching is not one-to-one: several joins of one walk may match the same join of the other;
- a linear inferred walk has no wrap-around join.

All lie in [0, 1]. Columns: `junction_recall`, `junction_precision`, `junction_f1` (unphased) and
`phased_junction_recall`, `phased_junction_precision`, `phased_junction_f1` (phased). The recovery
matrix's `truth_n_junctions` is $\lvert \mathcal A_N(T) \rvert$.

### 3.4 Order match and topology match

> **Provenance:** (c) new, RCK yubi_dev@c144584. The linear case is (c), yubi_dev@10416aa.
> **Code:** `rck/ecdna/evaluate.py::order_match`, `::_runs_equal`, `::_flip`, `::score_pair`. The
> criterion is (a) the exact-reconstruction criterion of AA's and CoRAL's simulation benchmarks
> [Deshpande et al. 2019](references.md#deshpande2019) [Zhu et al. 2024](references.md#zhu2024):
> the same segments in the same cyclic order and orientation. Treating rotations, and the reverse
> with every orientation flipped, as the same molecule is the usual equivalence for a circular
> sequence that can be read on either strand. It is the same set $\mathrm{Rot}$ that §4.4's CoRAL
> cyclic LCS searches. The following are (c): the comparison on reference runs rather than
> segments, the haplotype requirement, the search over every reported ordering, and the
> $\tau_{\mathrm{sp}}$ tolerance. §4.4's `exact_cycle` is the same function at $\tau_{\mathrm{gr}}$.

```math
\mathrm{order\_match}(T, C) = \mathbb 1\Big[\ \mathrm{circ}(T) = \mathrm{circ}(C)\ \wedge\ \exists\, \ell \le K_C,\ \exists\, \mathrm{rot} \in \mathrm{Rot} :\ \mathrm{runs}(T) \approx_{\tau_{\mathrm{sp}}} \mathrm{rot}\big(\mathrm{runs}(C^{(\ell)})\big)\ \Big]
```

where

- $C^{(\ell)}$ ranges over $C$'s orderings (§1); `order_match_ordering` records the first $\ell$
  that matches (empty when none does). $T$ uses its one ordering;
- $\mathrm{Rot}$ is every rotation of the run sequence and of its reverse. The reverse reverses the
  order and flips every $\mathrm{str}$. A linear walk is only taken as given or reversed, not
  rotated;
- $\approx_\tau$ means the same number of runs and, position by position, the same
  $\mathrm{chr}$, $H$ and $\mathrm{str}$, with $\mathrm{start}$ and $\mathrm{end}$ each within
  $\tau$.

`topology_match` is $\mathbb 1[\mathrm{circ}(T) = \mathrm{circ}(C)]$, and `inferred_topology`
names $C$'s (`circular` or `linear`). Truth walks are circular, so a linear inferred walk never
counts as an order match.

### 3.5 Copy-number error

> **Provenance:** (c) new, RCK yubi_dev@c144584. The integer variant is (c), yubi_dev@c5e0b9c
> (with the refit). **Code:** `rck/ecdna/evaluate.py::score_pair`. The quantity compared is (a) a
> cycle's copy count, AA's and CoRAL's `Copy_count`
> [Deshpande et al. 2019](references.md#deshpande2019) [Zhu et al. 2024](references.md#zhu2024).
> The per-species relative error against a truth copy number is (c). CoRAL's copy-number scores
> (§4.5 `cn_explained_topk`) compare cycles with the graph, not with the truth. $w_C$ comes from
> the refit LP of `docs/ecDNA.md` §6, and $w^{\mathrm{int}}_C$ from the decomposition MILP
> (`docs/ecDNA.md` §5, after CoRAL's integer cycle weights).

```math
\mathrm{cn\_rel\_err} = \frac{\lvert w_C - w_T \rvert}{w_T},
\qquad
\mathrm{cn\_rel\_err\_int} = \frac{\lvert w^{\mathrm{int}}_C - w_T \rvert}{w_T}
```

where

- $w_C$ is the inferred species' weight after the copy-number refit (`weight` in
  `decompositions.tsv`, reported as `inferred_weight`). In the multi-clone route it is
  $\sum_i u_i \psi_{C,i}$ (§3.9);
- $w^{\mathrm{int}}_C$ is its integer weight before the refit (`weight_int`, reported as
  `inferred_weight_int`); the multi-clone route leaves it empty, and so `cn_rel_err_int`;
- $w_T$ is the truth copies per cell (§2.1, `truth_weight`).

Both lie in $[0, \infty)$. Both are empty when $w_T = 0$; when $w_T$ is NaN (§2.1) they are NaN.
Reporting both makes the effect of the refit visible per species.

### 3.6 One-to-one matching of truth to inferred species

> **Provenance:** (c) new, RCK yubi_dev@c144584. **Code:** `rck/ecdna/evaluate.py::score_decomposition`,
> `::_match`. The procedure is the greedy heuristic for maximum-weight bipartite matching
> [Avis 1983](references.md#avis1983), here on haplotype Jaccard (§3.1). It is not the optimal
> assignment (the Hungarian method [Kuhn 1955](references.md#kuhn1955)). One-to-one matching is the
> choice of the GA4GH benchmarking framework, which pairs each call with at most one truth record
> [Krusche et al. 2019](references.md#krusche2019). It differs from CoRAL 3.1.1, which takes each
> metric's best value over all reconstructed cycles independently, so one cycle can serve several
> truth cycles (§4.6 `best_of`) [CoRAL 3.1.1 code](references.md#coral_code).

For the in-scope truth walks $\mathcal W^\ast_{\mathrm{in}}$ and the inferred walks $\mathcal W_r$
of one decomposition:

1. Score every pair $(T, C)$ by `haplotype_jaccard`.
2. Sort the triples (score, truth name, inferred name) in decreasing order (a Python tuple sort):
   highest score first; ties go to the larger truth name, then the larger inferred name, in
   string order.
3. Go through the list. Accept a pair unless its truth or inferred walk is already used, or its
   score is 0.

$\mathrm{match}_r(T)$ is the inferred walk accepted for $T$, or $\varnothing$. A truth species with
$\mathrm{match}_r(T) = \varnothing$ is *missed*: its row has no inferred species, every score 0 and
composite 0. An inferred species without a partner is *false*. Matching on haplotype Jaccard means
an unphased answer (every haplotype $\cdot$) matches nothing. That is why the method-agnostic scorer
also matches haplotype-stripped walks (§4.4).

### 3.7 Composite score

> **Provenance:** (c) new, RCK yubi_dev@4e53747. **Code:** `rck/ecdna/evaluate.py::composite_score`.

```math
\mathrm{composite}(T, C) = \frac{1}{5}\Big( J_{\mathrm{seg}} + F_1^{\mathrm{junc}} + \mathrm{PA} + \mathbb 1[\mathrm{order}] + \max(0,\ 1 - \mathrm{err}_{\mathrm{cn}}) \Big),
\qquad \mathrm{composite}(T, \varnothing) = 0
```

where

- $J_{\mathrm{seg}}$ is `segment_jaccard` (§3.1);
- $F_1^{\mathrm{junc}}$ is the *unphased* `junction_f1` (§3.3);
- $\mathrm{PA}$ is `phasing_accuracy` (§3.2, as given, not swap-invariant);
- $\mathbb 1[\mathrm{order}]$ is `order_match` (§3.4);
- $\mathrm{err}_{\mathrm{cn}}$ is `cn_rel_err` (§3.5); an empty or NaN $\mathrm{err}_{\mathrm{cn}}$ gives 0
  for that term.

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
(2026-10-31)"). The copy-number term is floored at 0, so a large copy-number error weighs no more
than a moderate one.

### 3.8 The recovery matrix, `best_*` and `mean_composite`

> **Provenance:** (c) new, RCK yubi_dev@4e53747. **Code:** `rck/ecdna/evaluate.py::recovery_matrix`,
> `::evaluate`, `::read_inferred`. The ranked alternatives $\mathcal W_r$ come from the
> decomposition pool and ambiguity report of `docs/ecDNA.md` §5 (c). Scoring every alternative and
> keeping the best is an oracle upper bound. It is the same idea as CoRAL's best-of scoring
> (§4.6), applied over ranked decompositions instead of over the cycles of one answer.
> `mean_composite` is a plain mean over truth species (macro averaging,
> [Sokolova & Lapalme 2009](references.md#sokolova2009)). The `low_weight` filter is
> `rck-ecdna`'s significance rule (`docs/ecDNA.md` §5.6; §7 below).

`decompositions.tsv` holds ranked decompositions with walks $\mathcal W_1, \mathcal W_2, \dots$
(`decomposition_rank`); the smallest rank, 1, is RCK's answer. Unless `--include-low-weight` is
given, species flagged `low_weight` are left out (`read_inferred`; `n_low_weight` counts them). A
rank whose species are all left out stays, with no walks. Every rank is scored (§3.1–§3.7, matching
§3.6) and gives one row per in-scope truth species in `evaluation.tsv`. For each
$T \in \mathcal W^\ast_{\mathrm{in}}$:

```math
\mathrm{composite}(T) = \mathrm{composite}\big(T,\ \mathrm{match}_{1}(T)\big),
\qquad
\mathrm{best\_composite}(T) = \max_{r} \mathrm{composite}\big(T,\ \mathrm{match}_{r}(T)\big)
```

```math
\mathrm{mean\_composite} = \frac{1}{\lvert \mathcal W^\ast_{\mathrm{in}} \rvert} \sum_{T \in \mathcal W^\ast_{\mathrm{in}}} \mathrm{composite}(T)
```

where

- $\mathrm{match}_r(T)$ is $T$'s partner in decomposition $r$ under the matching of §3.6, or $\varnothing$;
- the maximum runs over every reported decomposition. `best_rank` is the arg max, with ties going
  to the lower rank, and `best_inferred` is its matched species (empty when it is $\varnothing$);
- `mean_composite` is empty (`null`) when $\mathcal W^\ast_{\mathrm{in}} = \emptyset$.

`recovery_matrix.tsv` has one row per truth species, in or out of scope: the rank-1 scores (§3.1–§3.7)
and `best_rank`, `best_inferred`, `best_composite`; out-of-scope rows have only the truth columns.
`best_composite` is an upper bound on what the ranking could have picked from the reported
alternatives. `evaluation.json` also reports:

- `first_rank_with_all_orders_matched`: the smallest rank in which every in-scope truth species
  has `order_match`, or `null`;
- `best_mean_jaccard_hap`: the largest, over ranks, mean `haplotype_jaccard` over the in-scope
  truth species (a missed species counts 0), and its rank `best_mean_jaccard_rank` (ties to the
  lower rank).

**Failures.** If the truth files cannot be read or scoring fails, `rck-ecdna` still finishes and
writes the reconstruction outputs. The bulk route writes `evaluation.json` with the `error` and
its `traceback`. The multi-clone route writes no `evaluation.json`; `summary.json` `evaluation`
records the `error`.

### 3.9 Per-clone copies ψ (multi-clone route)

> **Provenance:** (c) new, RCK `yubi_CNP_ionic` fc22a34 (2026-11-13). **Code:**
> `rck/ecdna/evaluate.py::psi_scores`, called by `rck/ecdna/multiclone.py::evaluate_multiclone`
> when `rck-ecdna --model multiclone` has `--truth-dir` (`docs/ecDNA.md` §11.3).

The multi-clone route reports each cycle $q$'s integer copies per clone $i$, $\psi_{q,i}$
(`decompositions.tsv` column `psi_<clone>`).

**Truth.** The `ecdna_psi` rows of `truth/proportions.tsv` give `ID` (species), `CLONE` and
`PSI`: $\psi^\ast_{T,i}$, 0 when there is no row.

**Structure scoring is unchanged.** §3.1–§3.8 run on `decompositions.tsv` as for the bulk route,
with two differences:

- the paths (status `chromosomal`, the telomere-to-telomere chromosome walks) are skipped by
  `read_inferred`;
- a cycle's `weight` is its bulk copies per cell, $w_C = \sum_i u_i \psi_{C,i}$, comparable to the
  truth's bulk weight $w_T$.

**Pairing.** Each $T \in \mathcal W^\ast_{\mathrm{in}}$ is paired with $C = \mathrm{match}_1(T)$
(§3.6, on haplotype Jaccard, not the within-block pairing of §3.10). When
$\mathrm{match}_1(T) = \varnothing$, $\psi_{C,i} = 0$. The inferred $\psi$ are read from the rank-1
decomposition, ordering 1.

**Notation.** $\mathrm{Cl}$ is the union of the inferred clones (`psi_<clone>` columns) and the truth's
`CLONE` values. $\mathcal W^{\mathrm{cyc}}_1$ is every rank-1 row of `kind` `cycle`, *including*
`low_weight` cycles. $u_i$ is clone $i$'s proportion from the DIAG's `clone_proportions.tsv` (0 for
a clone it does not list).

**`evaluation_psi.tsv`** has one row per (truth species $T$, clone $i \in \mathrm{Cl}$), with `truth`,
`clone`, `inferred` and:

```math
\mathrm{truth\_psi} = \psi^\ast_{T,i},\quad
\mathrm{inferred\_psi} = \psi_{C,i},\quad
\mathrm{abs\_err} = \lvert \psi_{C,i} - \psi^\ast_{T,i} \rvert,\quad
\mathrm{truth\_proportion} = \frac{\psi^\ast_{T,i}}{\sum_{T' \in \mathcal W^\ast_{\mathrm{in}}} \psi^\ast_{T',i}},\quad
\mathrm{inferred\_proportion} = \frac{\psi_{C,i}}{\sum_{q \in \mathcal W^{\mathrm{cyc}}_1} \psi_{q,i}}
```

A proportion is empty when its denominator is 0.

**`evaluation.json`** `psi` holds:

- `n_pairs`, the number of rows, $\lvert \mathcal W^\ast_{\mathrm{in}} \rvert \cdot \lvert I \rvert$;
- `psi_exact`, the rows with `abs_err` = 0;
- `psi_mean_abs_err`, the mean `abs_err` (`null` without rows);
- `proportion_mean_abs_err`, the mean of
  $\lvert \mathrm{inferred\_proportion} - \mathrm{truth\_proportion} \rvert$ over the rows where
  both are defined (`null` when either column is empty throughout);
- `bulk`, per truth species: `truth_bulk` $= \sum_{i \in I} u_i \psi^\ast_{T,i}$, `inferred_bulk`
  $= \sum_{i \in I} u_i \psi_{C,i}$, `truth_proportion_bulk` (`truth_bulk` over its sum over
  $\mathcal W^\ast_{\mathrm{in}}$) and `inferred_proportion_bulk` (`inferred_bulk` over
  $\sum_{q \in \mathcal W^{\mathrm{cyc}}_1} \sum_i u_i \psi_{q,i}$);
- `rows`, the table above.

`summary.json` `evaluation` copies `psi_exact`, `psi_pairs` (`n_pairs`) and `psi_mean_abs_err`.

### 3.10 Phasing within phase blocks (`*_blockwise`)

> **Provenance:** (c) new, RCK `yubi_CNP_ionic` 140e938 (2026-11-13). **Code:**
> `rck/ecdna/evaluate.py::read_phase_blocks`, `::block_swaps`, `::relabel`, `::score_blockwise`,
> `::score_decomposition` (`blocks=`). Used by `rck-ecdna --model multiclone` with `--truth-dir` or
> `--truth-cycles` (blocks from the DIAG's `rck.fragments.tsv`, whenever that file exists) and by
> `rck-ecdna-eval --phase-blocks`. The bulk route of `rck-ecdna` does not use it.

**Why.** The haplotype-sensitive scores (phasing accuracy §3.2, phased junctions §3.3, order match
§3.4) compare A/B with the truth's across the whole chromosome.

- A phase block's observed A/B is arbitrary, and without phase information between blocks no
  method can tie neighbouring blocks together.
- The multi-clone route (`docs/ecDNA.md` §11.4) therefore reconstructs haplotypes correctly only
  *within* each phase block, and the absolute scores also count its arbitrary per-block A/B
  choice.

**Phase blocks.** $\mathcal B$ is read from an RCK fragments file (`read_phase_blocks`; columns
`chr`, `start`, `end`, `extra`): one block $b$ per row, the span $[\mathrm{start}, \mathrm{end}]$
taken as inclusive, named by `phase_block=` in `extra`, else `chr:start`. The DIAG writes one
fragment per phase block. $q\vert_b$ is $q$'s steps clipped to $b$'s span.

**Definition.** For one decomposition $\mathcal W_r$:

1. Pair truth and inferred walks one-to-one as in §3.6, but on `segment_jaccard` (haplotype-free).
   Let $\mathcal P$ be these pairs.
2. For every block $b$, swap A/B when the swapped haplotypes overlap the truth's more, in bases,
   summed over the pairs (ties: no swap):

```math
\mathrm{sw}_b = \mathbb 1\Big[ \sum_{(T, C) \in \mathcal P} O_H\big(T\vert_b,\ \bar C\vert_b\big) > \sum_{(T, C) \in \mathcal P} O_H\big(T\vert_b,\ C\vert_b\big) \Big]
```

   The swap is **one per block, shared by every species**, so a phase error *between* two species
   in one block still counts.
3. Relabel every inferred walk, paired or not, in every ordering: a step whose midpoint
   $\lfloor (\mathrm{start}_s + \mathrm{end}_s)/2 \rfloor$ lies in a block $b$ (the first in file
   order) with $\mathrm{sw}_b = 1$ gets A and B exchanged. A step in no block keeps its label.
   DIAG segments never cross a block end, so this agrees with the clipping of step 2.
4. Pair the relabelled walks with the truth again by haplotype Jaccard (§3.6).
5. On each new pair, score `phasing_accuracy_blockwise` (§3.2's $\mathrm{PA}$),
   `phased_junction_f1_blockwise` (§3.3, phased) and `order_match_blockwise` (§3.4, any ordering).
6. `composite_blockwise` is §3.7's composite with these in place of $\mathrm{PA}$ and
   $\mathbb 1[\mathrm{order}]$. When step 4 pairs $T$ with a different walk than §3.6 (or §3.6
   left $T$ missed), `inferred_blockwise` names it, and $J_{\mathrm{seg}}$,
   $F_1^{\mathrm{junc}}$ and $\mathrm{err}_{\mathrm{cn}}$ come from that pair.
7. A truth species left unpaired in step 4 scores 0 (`order_match_blockwise` false).

`mean_composite_blockwise` is the mean of the rank-1 `composite_blockwise` over
$\mathcal W^\ast_{\mathrm{in}}$.

**Outputs.**

- Without phase blocks, nothing changes: no new columns, and the bulk route's evaluation is
  byte-identical (checked on demo5, RCK `CHANGELOG.md`).
- With them, `evaluation.tsv` (every rank) and `recovery_matrix.tsv` (rank 1) gain the
  `*_blockwise` columns, and `evaluation.json` gains `mean_composite_blockwise` and
  `phase_blocks` (the number of blocks read).

**Example.** demo7 with the DIAG labelling groups of one junction label per junction has absolute
composites 0.716 and 0.670 (phasing accuracy 0.58 and 0.35), and 1.0 within blocks for both
species. Without the groups, both scores are 0.892 and 0.601: the inner species' error there is
structural, so the block swaps do not remove it. With per-end junction labels (amplicon-allele-analysis
`54226ed`), demo7 scores 1.0 both absolute and within blocks (RCK `CHANGELOG.md`).

---
## 4. The method-agnostic scorer `aa_eval`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409; `error_free_merged` and the CoRAL-code headline
> are from yubi_dev@41686bc. **Code:** `rck/ecdna/aa_eval.py::score`, `::main`.

The comparison with AA, CoRAL and the original RCK needs one scorer that reads every method's
output the same way. `aa_eval` reads the method's answer from AA-format cycles and graph files
only, and none of `rck-ecdna`'s own files. The optional inputs `--segment-cn` (§4.3),
`--paths` (§4.5) and `--phase-blocks` (§3.10) are plain tables. Species matching and the species
scores are §3's functions, called unchanged (`evaluate.score_decomposition` and
`evaluate.score_pair`). Truth species are read and scoped as in §2. `aa_eval` adds the AA-file
reading, the graph- and sample-level metrics, CoRAL's cycle metrics, two scoring modes (§4.6),
and the confusion-matrix metrics and species matrix of §4.7.

### 4.1 Reading AA files

> **Provenance:** (a) the AA graph and cycles formats of AmpliconArchitect [Deshpande et al. 2019](references.md#deshpande2019),
> as written by ecSimulator's `ecSim_IO.py` [ecSimulator](references.md#ecsimulator) and by CoRAL 3.1.1
> [CoRAL 3.1.1 code](references.md#coral_code). The optional `Cycle=` keys are (c) new, RCK yubi_dev@e4c2409;
> `Status=chromosomal` is (c), RCK eb5c85a. **Code:** `rck/ecdna/aa_format.py::read_cycles`, `::to_cycles`,
> `::read_graph`, `::AAGraph`, `::AAEdge`, `::AASequence`; `rck/ecdna/aa_eval.py::species_cycles`.

Both files are tab-separated. Lines starting with `#` are comments. Coordinates are AA's:
1-based, inclusive.

**Breakend (extremity).** $v = (\mathrm{chr}, \mathrm{pos}, \pm)$, written `chr:pos±`. A sequence
edge $j$ has two: its left end $j_t = (\mathrm{chr}_j, \mathrm{start}_j, -)$ (`chr:start-`) and its
right end $j_h = (\mathrm{chr}_j, \mathrm{end}_j, +)$ (`chr:end+`). These are RCK's extremities
$j_\sigma$, $\sigma \in \lbrace t, h \rbrace$.

**Graph file** (`*_graph.txt`, `read_graph` → `AAGraph`). One edge per line; the first field is
its kind:

| Line | Fields read | Object |
|---|---|---|
| `sequence` | `chr:start-`, `chr:end+`, PredictedCopyCount, AverageCoverage, Size, NumberReadsMapped | a *sequence edge* (segment) $j = (\mathrm{chr}_j, \mathrm{start}_j, \mathrm{end}_j)$ with copy count $\mu(j)$; length $L(j) = \mathrm{end}_j - \mathrm{start}_j + 1$ |
| `concordant` | `v->v'`, PredictedCopyCount, NumberOfReadPairs | a *reference adjacency* $a = \lbrace v, v' \rbrace$ with copy count $\mu(a)$; it joins $j_h$ to $j'_t$ of the next segment |
| `discordant` | as `concordant` | a *novel adjacency* $a = \lbrace v, v' \rbrace$, any two breakends (both may be equal: a fold-back) |
| `source` | as `concordant` | an end of a path: copy number entering or leaving the graph at $v$ (AA's source node) |
| `path_constraint` | tokens, support | kept as text (§5.2) |
| `interval` | chr, start, end | `AAGraph.intervals` (CoRAL 3.1.1 ignores them) |

Only the fields listed are used: the PredictedCopyCount of every edge, and the breakends. The
remaining fields (homology) are ignored. Edges are unordered: $\lbrace v, v' \rbrace$ and
$\lbrace v', v \rbrace$ are the same edge for every match below.

**Cycles file** (`*_cycles.txt`, `read_cycles` → `AACycles`):

- `Interval  i  chr  start  end`: the region the method looked at;
- `Segment  id  chr  start  end`: a segment the walks use; a `Cycle=` step whose id is not
  listed is an error;
- `Cycle=k;Copy_count=w;Segments=id±,id±,...`: one walk, with weight $w_C$ = `Copy_count`
  (`nan` when absent). A walk whose `Segments=` contains segment `0` (AA's source, `0+ … 0-`;
  ecSimulator writes `0+ … 0+`) is a linear path, and every `0` token is dropped. A walk with no
  step left is skipped. The walk is named `Species=` when given, else `cycle<k>`.

Optional keys on a `Cycle=` line, which AA and CoRAL ignore:

- `Haplotypes=`: one haplotype per `Segments=` token, the `0` tokens included and then dropped;
  missing trailing entries are `.`. It makes the answer phased; without it every step has
  $H_s = \cdot$ and the scores are haplotype-blind. A file is phased only if it has at least one
  `Cycle=` line and every one, and every `--alt-cycles` file, carries `Haplotypes=`;
- `Status=low_weight`: the method's own non-species (§4.6);
- `Status=chromosomal`: a telomere-to-telomere path of the multi-clone route; it is dropped
  before scoring (`species_cycles`);
- `Alt_segments=` / `Alt_haplotypes=`: a species' other orderings, `|`-separated. An ordering
  without its own `Alt_haplotypes=` entry takes the walk's `Haplotypes=`;
- `Species=`: a name.

**Steps.** Each step `id±` of segment `id` = $(\mathrm{chr}, \mathrm{start}, \mathrm{end})$
becomes the step $s = (\mathrm{chr},\ \mathrm{start} - 1,\ \mathrm{end} - 1,\ H_s,\ \pm)$ of §1, a
0-based inclusive interval. This is the frame in which the simulator's AA truth is read (§2.1).

**Linear walks of other tools.** `read_cycles` reads `Cycle=` lines only. CoRAL 3.1.1 and
CycleExtractor write linear walks as `Path=` lines, which it skips unless they are rewritten as
`Cycle=`. `results/graph_route/score_h2h.py` does that (§7).

### 4.2 Truth mixture graph

> **Provenance:** (c) new, amplicon-allele-analysis yubi_dev@40c6c73; the graph-route options are
> from yubi_dev@5e2da37. Grouping into amplicons, one graph and cycles file per amplicon, is (a)
> ecSimulator's output layout [ecSimulator](references.md#ecsimulator).
> **Code:**
> [`simulations/truth_mixture.py::build`](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/truth_mixture.py),
> `::group_amplicons`, `::mixture_segments`, `::merge_silent`, `::species_coverage`,
> `::background_junctions`, `::write_amplicon`; RCK `rck/ecdna/aa_eval.py::load_mixture`.

The species truth (§2) is one walk per species. The graph truth is a different object: the
amplicon as a method sees it, with every species at its copies per cell plus the background
genome. The simulator writes only one AA graph per species, each with the copy counts of one
molecule, so `truth_mixture.py` builds the mixture for one tumor sample (`--sample`, default
the first) from:

- `truth/proportions.tsv`: clone shares $u_i$ (rows `KIND = clone`, column `U`) and species
  copies per cell $w_T$ (rows `KIND = ecdna_burden`, column `PHI`; 0 when absent);
- `truth/clones.tsv`: each clone's background allele copy numbers (columns `A_<clone>`,
  `B_<clone>`);
- `truth/segments.tsv`: the sample's allele copy numbers (`CN_A`, `CN_B`);
- each species' `truth/aa/<species>_cycles.txt` and `_graph.txt`;
- the background junctions: `truth/adjacency_clone.tsv` (rows with `IS_ECDNA` not true),
  their coordinates in `obs/adjacency.tsv`, and their per-clone copies in `truth/rck/rck.acnt.tsv`.

The truth tables are 0-based; a simulator position $p$ becomes AA position $p + 1$.

**Amplicons follow ecSimulator.** ecSimulator simulates each amplicon on its own and writes one
`<prefix>_amplicon<N>_graph.txt` / `_cycles.txt` per amplicon; it has no joint multi-amplicon
graph. Here a truth amplicon $\mathrm{Amp}$ is a connected group of species whose intervals
overlap or abut (`group_amplicons`, union–find). Amplicons are numbered $N = 1, 2, \dots$ by
their smallest interval (chromosome name as text, then start). The amplicon's intervals
$I_{\mathrm{Amp}}$, each $(\mathrm{chr}, r_0, r_1)$, 1-based inclusive, are its species' intervals
merged where they overlap or abut. demo6 has two amplicons; demo7's nested pair has one.

**Construction** (`build`), for each amplicon:

1. **Cuts.** Cut each interval of $I_{\mathrm{Amp}}$ at every species segment boundary, every
   `segments.tsv` boundary, and every breakend of a background junction kept in step 4 (a `-`
   breakend at $p$ starts a segment at $p$, a `+` breakend at $p$ ends one at $p$). The pieces are
   the mixture segments $j$ (`mixture_segments`).
2. **Segment copy numbers.** For each mixture segment $j$:

   ```math
   \mu^\ast(j) = \sum_i u_i\,(a_{i,j} + b_{i,j}) \;+\; \sum_{T \in \mathcal W^\ast_{\mathrm{Amp}}} w_T\, \mathrm{cov}_T(j)
   ```

   where $i$ runs over the clones of the sample; $a_{i,j}, b_{i,j}$ are clone $i$'s background
   copy numbers of haplotypes A and B over $j$ (`clones.tsv`, length-weighted mean over $j$);
   $\mathcal W^\ast_{\mathrm{Amp}}$ is the amplicon's truth species; and $\mathrm{cov}_T(j)$ is
   the number of steps of $T$ whose segment contains $j$.
3. **Allele copy numbers.** $\mu^\ast(j^A)$, $\mu^\ast(j^B)$ are the length-weighted means of
   `segments.tsv` `CN_A`, `CN_B` over $j$, for the first tumor sample only (the one
   `segments.tsv` describes); empty for another sample. They are not forced to sum to
   $\mu^\ast(j)$; a mismatch above $10^{-6}$ relative is printed as a warning.
4. **Novel adjacencies.** The discordant edges are every species junction and every background
   junction with at least one breakend inside $I_{\mathrm{Amp}}$:

   ```math
   \mu^\ast(a) = \sum_{T \in \mathcal W^\ast_{\mathrm{Amp}}} w_T\, \mathrm{trav}_T(a) \;+\; \sum_i u_i\, \mu_i(a)
   ```

   where $\mathrm{trav}_T(a)$ is how many times $T$ traverses $a$ (the PredictedCopyCount of the
   species' own `_graph.txt`, or, without it, counted from the walk), and $\mu_i(a)$ is clone
   $i$'s copies of background junction $a$, summed over its labellings (`rck.acnt.tsv`; without
   that file, 1 for every clone in `SHARED_CLONES`). A background junction with
   $\sum_i u_i \mu_i(a) \le 0$ is dropped.
5. **Silent boundaries** (`merge_silent`). Two neighbouring segments $j, j'$
   ($\mathrm{start}_{j'} = \mathrm{end}_j + 1$) become one when no discordant edge has a breakend
   exactly at $j_h$ or $j'_t$ and $\mu^\ast$, $\mu^\ast(\cdot^A)$ and $\mu^\ast(\cdot^B)$ agree
   (relative $10^{-9}$). No method can see such a boundary. ecSimulator keeps one for an SV that a
   later event removed.
6. **Reference adjacencies.** Each pair of neighbouring segments gets a concordant edge
   $\lbrace j_h, j'_t \rbrace$ with what node balance leaves:

   ```math
   \mu^\ast(\lbrace j_h, j'_t \rbrace) = \max\Big(0,\ \min\big(\mu^\ast(j) - \mathrm{disc}(j_h),\ \mu^\ast(j') - \mathrm{disc}(j'_t)\big)\Big)
   ```

   where $\mathrm{disc}(v)$ is the summed $\mu^\ast(a)$ of the discordant edges with a breakend at $v$ (a
   fold-back counted once). There are no source edges, and no edge across a gap between
   intervals.

**Output files** (in `truth/aa/`):

| File | Content |
|---|---|
| `mixture_amplicon<N>_graph.txt` | a `#` header naming the sample; `sequence` lines ($\mu^\ast(j)$; AverageCoverage = $\mu^\ast(j)$; NumberReadsMapped 0); `concordant`, then `discordant` lines ($\mu^\ast(a)$; NumberOfReadPairs 0) |
| `mixture_amplicon<N>_cycles.txt` | `Interval` lines ($I_{\mathrm{Amp}}$); `Segment` lines numbered 1, 2, … in segment order; one `Cycle=k;Copy_count=`$w_T$`;Length=`$\lvert T\rvert$`;Segments=…;Species=<name>` per species, its walk cut into the mixture segments (linear: `0+ … 0-`) |
| `mixture_amplicon<N>_segments.tsv` | `chrom start end cn cn_a cn_b`: $j$, $\mu^\ast(j)$, $\mu^\ast(j^A)$, $\mu^\ast(j^B)$, AA coordinates (alleles, which an AA graph cannot hold) |
| `mixture_amplicons.tsv` | one row per amplicon: `amplicon` ($N$), `species` (comma-separated), `intervals`, `sample`, `n_segments`, `n_discordant`, and the three file names `graph`, `cycles`, `segments` |

`aa_eval` reads them with `load_mixture`: $I_{\mathrm{Amp}}$ from the cycles file's `Interval`
lines, $G^\ast$ from the graph file, the allele table from `_segments.tsv`.

**Where the files are.** The tracked demo5–demo8 runs carry theirs in amplicon-allele-analysis
(`simulations/runs/demo{5,6,7,8}_*/truth/aa/mixture_*`). The bench simulations
(`/n/fs/yubi-iw/runs/`) carry theirs in place. `results/bench_nested/diag/aa_regression.py`
uses `truth/aa/mixture_amplicons.tsv` wherever it exists and otherwise builds one (demo4: under
`/n/fs/yubi-iw/runs/bench_nested_aa/mixtures/`).

**Options for the graph route.** Opt-in; without them every file is byte-identical:

- `--reads-per-copy R`: NumberOfReadPairs = $\mathrm{round}(\mu^\ast(a) \times R)$ on concordant
  and discordant edges;
- `--depth-per-copy D`: AverageCoverage = $\mu^\ast(j) \times D$;
- `--path-constraints`: the run's long-read paths (`obs/path.tsv`) as CoRAL 3.1.1
  `path_constraint` lines (support = `#NOV`). A read is kept only if its breakends lie on the
  amplicon's breakends (within 1 bp) and every junction it crosses is a discordant edge of the
  mixture. The graph's lines are then written in CoRAL's sort order, which the constraints index,
  and `mixture_amplicon<N>_path_constraints.tsv` lists each constraint's read and junctions.

### 4.3 Graph-level metrics (`aa_graph_metrics.tsv`)

> **Code:** `rck/ecdna/aa_eval.py::graph_metrics`, `::edges_match`, `::scoped_edges`,
> `::scoped_sequences`, `::_segment_match`. **Provenance of the shared objects:** sequence and
> discordant edges, and amplicon intervals, are (a) AmpliconArchitect's breakpoint graph and
> `interval` lines [Deshpande et al. 2019](references.md#deshpande2019). $a \sim a'$ is (a) CoRAL
> 3.1.1's `is_breakpoint_match` with $\le$ [CoRAL 3.1.1 code](references.md#coral_code). The
> following are (c): scoping to one truth amplicon ($\hat{\mathcal A}_N$ and
> $\widehat{\mathrm{Seg}}$, `bp_accuracy` difference 1) and the segment match $j \simeq j'$. Each
> metric below has its own provenance line.

The metrics need `--graph` and a truth mixture (`--mixture-dir`, default
`<truth-dir>/truth/aa` when `mixture_amplicons.tsv` is there); without either the table is empty.
They are computed once per truth amplicon $\mathrm{Amp}$ whose `species` include an in-scope
truth species (§2.3; the species name up to its first `.`), one row each. All matches use
$\tau_{\mathrm{gr}}$ (`--graph-tol`, default 100 bp; column `graph_tol_bp`). Symbols:

| Symbol | Meaning |
|---|---|
| $G^\ast$, $\hat G$ | the amplicon's truth mixture graph (§4.2) and the method's graph (`--graph`) |
| $I_{\mathrm{Amp}}$ | the amplicon's intervals $(\mathrm{chr}, r_0, r_1)$, 1-based inclusive (§4.2) |
| $\mathcal A_N^\ast$ | all discordant edges of $G^\ast$ |
| $\hat{\mathcal A}_N$ | the discordant edges of $\hat G$ with at least one breakend $(\mathrm{chr}, \mathrm{pos}, \pm)$ such that $r_0 - \tau_{\mathrm{gr}} \le \mathrm{pos} \le r_1 + \tau_{\mathrm{gr}}$ for some $(\mathrm{chr}, r_0, r_1) \in I_{\mathrm{Amp}}$ |
| $\mathrm{Seg}^\ast$ | all sequence edges of $G^\ast$ |
| $\widehat{\mathrm{Seg}}$ | the sequence edges $j$ of $\hat G$ with $\sum_{(\mathrm{chr}_j, r_0, r_1) \in I_{\mathrm{Amp}}} \lvert [\mathrm{start}_j, \mathrm{end}_j] \cap [r_0, r_1] \rvert > \tau_{\mathrm{gr}}$ (a flank that only touches the amplicon is not part of it) |
| $v \approx v'$ | same chromosome, same strand, $\lvert \mathrm{pos} - \mathrm{pos}' \rvert \le \tau_{\mathrm{gr}}$ (§1's extremity match) |
| $a \sim a'$ | $a = \lbrace v_1, v_2 \rbrace$, $a' = \lbrace v'_1, v'_2 \rbrace$ with $v_1 \approx v'_1 \wedge v_2 \approx v'_2$, or $v_1 \approx v'_2 \wedge v_2 \approx v'_1$ (`edges_match`; CoRAL's `is_breakpoint_match` with $\le$ instead of $<$) |
| $j \simeq j'$ | $\mathrm{chr}_j = \mathrm{chr}_{j'}$, $\lvert \mathrm{start}_j - \mathrm{start}_{j'} \rvert \le \tau_{\mathrm{gr}}$ and $\lvert \mathrm{end}_j - \mathrm{end}_{j'} \rvert \le \tau_{\mathrm{gr}}$ |

Copy counts, read counts and edge kinds other than `sequence` and `discordant` play no part in
these matches. Count columns: `n_true_edges` $= \lvert \mathcal A_N^\ast \rvert$, `n_inferred_edges`
$= \lvert \hat{\mathcal A}_N \rvert$, `n_edges_found` (truth edges matched), `n_false_edges`
(method edges matching none); `n_true_segments`, `n_inferred_segments`, `n_segments_found`,
`n_false_segments` likewise for $\mathrm{Seg}^\ast$ and $\widehat{\mathrm{Seg}}$. A row also
carries `amplicon`, `species` and `intervals`.

#### `bp_accuracy`: breakpoint-graph accuracy

> **Provenance:** (a) CoRAL's breakpoint-graph accuracy [Zhu et al. 2024](references.md#zhu2024). In CoRAL 3.1.1 it is the
> count `n_breakpoint_edges_covered` of `coral/scoring/score_simulation.py::score_reconstruction`
> over `n_breakpoint_edges`, counted by `coral/scoring/scoring_utils.py::find_overlapping_bp_edges`
> [CoRAL 3.1.1 code](references.md#coral_code).

```math
\mathrm{bp\_accuracy} = \frac{\lvert \lbrace a \in \mathcal A_N^\ast : \exists a' \in \hat{\mathcal A}_N,\ a \sim a' \rbrace \rvert}{\lvert \mathcal A_N^\ast \rvert} \in [0, 1]
```

It is empty (undefined) when $\mathcal A_N^\ast = \emptyset$. Differences from CoRAL's code:

1. CoRAL counts the true edges against each reconstructed amplicon graph (`*_graph.txt`) and
   keeps the largest count. Here the method's whole graph is scoped to $I_{\mathrm{Amp}}$.
2. CoRAL uses a strict $<$.

#### `novel_adjacency_fnr`, `novel_adjacency_fdr`

> **Provenance:** (a) RCK's novel-adjacency recovery [Aganezov & Raphael 2020](references.md#aganezov2020). The rates are written
> here as FNR = 1 − recall and FDR. The FDR is (c), chosen instead of an FPR, which would need a
> negative set.

```math
\mathrm{novel\_adjacency\_fnr} = 1 - \mathrm{bp\_accuracy},
\qquad
\mathrm{novel\_adjacency\_fdr} = \frac{\lvert \lbrace a' \in \hat{\mathcal A}_N : \nexists\, a \in \mathcal A_N^\ast,\ a' \sim a \rbrace \rvert}{\lvert \hat{\mathcal A}_N \rvert}
```

Both are in $[0, 1]$. The FNR is empty when $\mathcal A_N^\ast = \emptyset$, the FDR when
$\hat{\mathcal A}_N = \emptyset$. Matching is many-to-many: one method edge can match several
true edges.

#### `error_free`: AA's error-free reconstruction

> **Provenance:** (a) AmpliconArchitect's simulation benchmark [Deshpande et al. 2019](references.md#deshpande2019).

```math
\mathrm{error\_free} = \mathbb 1\Big[ \forall j \in \mathrm{Seg}^\ast\, \exists j' \in \widehat{\mathrm{Seg}} : j \simeq j' \ \wedge\ \forall j' \in \widehat{\mathrm{Seg}}\, \exists j \in \mathrm{Seg}^\ast : j' \simeq j \ \wedge\ \forall a \in \mathcal A_N^\ast\, \exists a' \in \hat{\mathcal A}_N : a \sim a' \ \wedge\ \forall a' \in \hat{\mathcal A}_N\, \exists a \in \mathcal A_N^\ast : a' \sim a \Big]
```

In words: every true segment is found with both ends within tolerance, there is no false
segment, every true discordant edge is found, and there is no false one. Copy numbers are not
compared. Empty sets satisfy their clauses. `rck-ecdna`'s imputed breakpoints (the subgraph's
10 kb grid) count as false segments.

#### `error_free_merged`: the same rule after merging benign splits

> **Provenance:** (c) modified from AA's rule, RCK yubi_dev@41686bc. **Code:**
> `rck/ecdna/aa_eval.py::merge_benign_splits`.

```math
\mathrm{error\_free\_merged} = \mathrm{error\_free}\ \text{with}\ \mathrm{Seg}^\ast, \widehat{\mathrm{Seg}}\ \text{replaced by}\ \mathrm{Seg}^\ast_{\mathrm{m}}, \widehat{\mathrm{Seg}}_{\mathrm{m}}
```

$\mathrm{merge}(G)$ (`merge_benign_splits`) is $G$'s sequence edges with every benign split
removed:

1. Sort the sequence edges by (chromosome, start, end).
2. Walk them in order, keeping a current segment $j$. The next segment $j'$ is merged into $j$
   (new end $\mathrm{end}_{j'}$, copy count the length-weighted mean, unused here) when all hold:
   - $\mathrm{chr}_{j'} = \mathrm{chr}_j$ and $\mathrm{start}_{j'} = \mathrm{end}_j + 1$;
   - a concordant edge matches $\lbrace j_h, j'_t \rbrace$ within 1 bp ($a \sim a'$ at
     tolerance 1);
   - no breakend of a discordant or source edge matches $j_h$ or $j'_t$ within 1 bp (same
     chromosome and strand).

   Otherwise $j'$ becomes the current segment. Merges chain, so a run of benign splits becomes
   one segment.

$\mathrm{Seg}^\ast_{\mathrm{m}}$ and $\widehat{\mathrm{Seg}}_{\mathrm{m}}$ are
$\mathrm{merge}(G^\ast)$ and $\mathrm{merge}(\hat G)$, each kept when it overlaps $I_{\mathrm{Amp}}$
by more than $\tau_{\mathrm{gr}}$ bp (as $\widehat{\mathrm{Seg}}$; here the truth side is scoped
too). The discordant-edge part of the rule is unchanged. Columns `n_true_segments_merged`,
`n_inferred_segments_merged`, `n_segments_found_merged` and `n_false_segments_merged` count the
merged segments. `error_free` stays the headline.

**Why.** AA's rule asks for every segment boundary. Graphs also cut segments where the copy
number changes without a breakpoint, and `rck-ecdna`'s subgraphs add a 10 kb imputed grid. The
strict rule therefore scores a method near 0 for boundaries that no breakpoint explains. On the
grown bench at λ = 0.2 (40 cases, 42 truth amplicons), 1 truth amplicon is `error_free` (demo6's chr17
amplicon), and all 532 inferred segment boundaries that match no true boundary are imputed
extremities of the subgraph. `error_free_merged` merges such boundaries first, in both graphs,
and so tests the breakpoints and the segments they delimit. A boundary that carries a discordant
or source edge is never merged, so a missed or false breakpoint still fails it.

**Effect.** 38 of 42 truth amplicons pass `error_free_merged` (RCK `results/bench_nested/README.md`,
"Method-agnostic metrics (2026-11-03)"). Two of the four failures use an
imputed discordant edge (`ctrl_dup_s710`, `nest_samehap_6to1_s804`). The other two miss a
copy-number-2 flank segment inside the truth interval that the subgraph shortens (`demo4_amp0`,
`nest_samehap_6to1_s803`).

#### `segment_cn_mae`, `segment_cn_rel_err`, `segment_cn_mae_total`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. A length-weighted absolute copy-number error in
> the spirit of RCK's allele-specific segment copy-number distance [Aganezov & Raphael 2020](references.md#aganezov2020)
> (ecDNA.md §4.1). **Code:** `rck/ecdna/aa_eval.py::segment_cn_error`, `::_piecewise`.

Over the rows $j$ of the amplicon's `_segments.tsv` (all mixture segments of $G^\ast$, §4.2),
with span $[\mathrm{start}_j, \mathrm{end}_j]$ (AA coordinates), length $L(j)$, total
$\mu^\ast(j)$ and alleles $\mu^\ast(j^A)$, $\mu^\ast(j^B)$:

```math
\mathrm{err}_{\mathrm{tot}} = \sum_j \sum_{\mathrm{pos} = \mathrm{start}_j}^{\mathrm{end}_j} \big\lvert \hat\mu(\mathrm{pos}) - \mu^\ast(j) \big\rvert,
\qquad
\mathrm{err}_{\mathrm{al}} = \sum_j \sum_{\mathrm{pos} = \mathrm{start}_j}^{\mathrm{end}_j} \Big( \big\lvert \hat\mu_A(\mathrm{pos}) - \mu^\ast(j^A) \big\rvert + \big\lvert \hat\mu_B(\mathrm{pos}) - \mu^\ast(j^B) \big\rvert \Big)
```

```math
\mathrm{segment\_cn\_mae} = \frac{\mathrm{err}}{\sum_j L(j)},
\qquad
\mathrm{segment\_cn\_rel\_err} = \frac{\mathrm{err}}{\sum_j \mu^\ast(j)\, L(j)},
\qquad
\mathrm{segment\_cn\_mae\_total} = \frac{\mathrm{err}_{\mathrm{tot}}}{\sum_j L(j)}
```

where

- $\hat\mu(\mathrm{pos})$, $\hat\mu_A(\mathrm{pos})$ and $\hat\mu_B(\mathrm{pos})$ are the method's total and haplotype A and B
  copy numbers at base $\mathrm{pos}$. With `--segment-cn` (a TSV `chrom start end [cn] cn_a cn_b`, AA
  coordinates), totals are `cn`, or `cn_a + cn_b` where `cn` is missing or empty, and the graph is
  not used. Without it, totals are $\mu(j')$ of $\hat G$'s sequence edges $j'$ (all of them, not
  only $\widehat{\mathrm{Seg}}$). Alleles come from `--segment-cn` only;
- a base no method row covers has copy number 0; where method rows overlap, the one with the
  smaller start wins;
- $\mathrm{err} = \mathrm{err}_{\mathrm{al}}$ (column `cn_error_mode` = `allele`) when
  `--segment-cn` has columns `cn_a` and `cn_b` and every truth `cn_a` is set (first tumor sample,
  §4.2); else $\mathrm{err} = \mathrm{err}_{\mathrm{tot}}$ (`total`).

Units: copies per base (`mae`), a fraction of the truth's copy number (`rel_err`; in `allele` mode
it can exceed 1). Empty when the denominator is 0. Unlike RCK's distance, there is no minimum
over the A/B swap: the method's A and B are compared with the truth's A and B as labelled.

#### `is_violations`, `is_violations_truth`

> **Provenance:** (a) the infinite-sites assumption of RCK [Aganezov & Raphael 2020](references.md#aganezov2020): extremity
> exclusivity, under which an extremity takes part in at most one novel adjacency (a constraint of
> RCK's paper; RCK's code only checks it after solving, ecDNA.md §4.1). This count of violations is
> (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_eval.py::is_violations`.

`is_violations` counts on $\hat{\mathcal A}_N$, `is_violations_truth` on $\mathcal A_N^\ast$:

1. Keep distinct edges: an edge equal to an earlier one, breakends exactly equal in either
   direction, is dropped. Edges that differ by a few bp stay distinct.
2. Give each edge its breakends, grouped by (chromosome, strand).
3. In each group, sort the positions and cut wherever two consecutive positions are more than
   $\tau_{\mathrm{gr}}$ apart (single linkage). Each cluster is one extremity.
4. The count is the number of extremities used by more than one distinct edge. An edge with both
   breakends in one extremity (a fold-back) counts once there.

The count is haplotype-blind: AA graphs carry no haplotype, so two novel adjacencies at one
extremity on different haplotypes, which RCK's per-haplotype exclusivity allows, count as a
violation. Simulated duplications reuse extremities, so the truth itself violates the
assumption: the truth has 25 violations over the 40 cases, and RCK's graphs the same 25 (19 on
the 32 bench cases; RCK `results/bench_nested/README.md`).
### 4.4 Species-level metrics (`aa_species_metrics.tsv`)

> **Code:** `rck/ecdna/aa_eval.py::strict_rows`, `::best_of_rows`, `::phased_part`, `::blind_part`.
> **Provenance:** the column families are §3's functions (c), called unchanged. The haplotype-blind
> matching (`strip_haplotypes`, `inferred_h0`) is (c), RCK yubi_dev@e4c2409. It exists so that
> AA, CoRAL and original RCK, which write unphased cycles, can be matched at all. The CoRAL
> columns are (a), and each has its own provenance line below.

Each in-scope truth species $T \in \mathcal W^\ast_{\mathrm{in}}$ (§2.3) gets one row per mode
(§4.6). An answer is *phased* when every `Cycle=` line of `--cycles` and of every `--alt-cycles`
file carries `Haplotypes=` (`aa_format.py::AACycles.phased`). There are two families of columns:

- **The `rck-ecdna-eval` columns** of §3, at $\tau_{\mathrm{sp}}$, under the recovery-matrix names
  (`segment_recall`, `segment_precision`, `segment_jaccard`, `haplotype_jaccard`,
  `junction_*`, `phased_junction_*`, `phasing_accuracy`, `phasing_accuracy_swap_invariant`,
  `order_match`, `topology_match`, `cn_rel_err`, `composite`), plus `inferred`,
  `inferred_weight`, `inferred_topology` and `exact_cycle_phased`. `cn_rel_err_int`,
  `order_match_ordering` and the `*_blockwise` scores are not written. For a phased answer they
  score the partner of the phased matching (§3.6). For an unphased answer they score the
  partner of the haplotype-blind matching, and the haplotype columns (`haplotype_jaccard`,
  `phased_junction_*`, `phasing_accuracy*`), `composite` and `exact_cycle_phased` are empty.
- **Haplotype-blind columns**, from a *separate* one-to-one matching (§3.6) of
  haplotype-stripped walks $T^0$, $C^0$ (`strip_haplotypes` sets every $H_s$ to `.`). The matched
  species is `inferred_h0`. Every column below except `exact_cycle_phased` scores this partner,
  for phased answers too.

The CoRAL columns use each walk's first ordering. The overlap and `length_error` are computed
on the walks' steps, the cyclic LCS on their reference runs (`evaluate.py::runs`, §1; the truth
always cyclically). The overlap does not depend on this choice up to the runs' 2 bp slack.

#### `composite_h0`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_eval.py::composite_h0`.

```math
\mathrm{composite\_h0}(T, C) = \frac{1}{4}\Big( J_{\mathrm{seg}}^{0} + F_1^{0} + \mathbb 1[\mathrm{order}^{0}] + \max(0,\ 1 - \mathrm{err}_{\mathrm{cn}}^{0}) \Big),
\qquad \mathrm{composite\_h0}(T, \varnothing) = 0
```

The superscript 0 marks the §3 scores of the haplotype-stripped pair $(T^0, C^0)$ matched by the
haplotype-blind matching, at $\tau_{\mathrm{sp}}$:

- $J_{\mathrm{seg}}^{0}$ is `segment_jaccard_h0` (§3.1);
- $F_1^{0}$ is `junction_f1_h0` (§3.3, unphased). It can differ from `junction_f1` of the same
  pair: once haplotypes are stripped, a haplotype switch between reference neighbours is a
  continuation, not a join (§1);
- $\mathbb 1[\mathrm{order}^{0}]$ is `order_match_h0` (§3.4, over every ordering of $C$);
- $\mathrm{err}_{\mathrm{cn}}^{0}$ is `cn_rel_err_h0` (§3.5); a missing value ($w_T = 0$) gives 0 for that term.

The range is $[0, 1]$.

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
> `rck/ecdna/aa_eval.py::cycle_overlap_coral_code`, `::cycle_overlap`, `::_merged`, `::_bins`, `::_counts`.

**Joint bins.** For each chromosome, sort the distinct points of
$\lbrace \mathrm{start}_s,\ \mathrm{end}_s + 1 \rbrace$ over every step $s$ of both walks into
$p_1 < \dots < p_n$. The bins are the intervals $\Delta = [p_i,\ p_{i+1} - 1]$, with length
$\lvert \Delta \rvert = p_{i+1} - p_i$. $\mathrm{Bins}(T, C)$ is the set of bins on every chromosome.
Every step covers whole bins. For $X \in \lbrace T, C \rbrace$, $m_X(\Delta)$ is the number of steps
of $X$ that cover $\Delta$, that is the multiplicity $m_X(\mathrm{chr}, \mathrm{pos})$ of §1 at any base
$\mathrm{pos} \in \Delta$.

```math
\mathrm{cycle\_overlap\_coral\_code}(T, C) = 1 - \frac{\sum_{\Delta \in \mathrm{Bins}} \lvert m_T(\Delta) - m_C(\Delta) \rvert\, \lvert \Delta \rvert}{\sum_{\Delta \in \mathrm{Bins}} \max\big(m_T(\Delta),\ m_C(\Delta)\big)\, \lvert \Delta \rvert}
= \frac{\sum_{\Delta} \min\big(m_T(\Delta),\ m_C(\Delta)\big)\, \lvert \Delta \rvert}{\sum_{\Delta} \max\big(m_T(\Delta),\ m_C(\Delta)\big)\, \lvert \Delta \rvert}
```

```math
\mathrm{cycle\_overlap\_supp}(T, C) = \frac{\lvert \mathrm{Cov}(T) \cap \mathrm{Cov}(C) \rvert}{\lvert \mathrm{Cov}(T) \cup \mathrm{Cov}(C) \rvert}
```

where $\mathrm{Cov}(X)$ is the set of (chromosome, base) pairs covered by at least one step of
$X$, without multiplicity (`_merged`). Intervals on different chromosomes never overlap. Both
lie in $[0, 1]$. Both are 0 for a missed species, and 0 when the denominator is 0.

The second form of the code version is $\max - \lvert\cdot\rvert = \min$ binwise. It is the
multiplicity-weighted Jaccard of §3.1, so `cycle_overlap_coral_code` equals
`segment_jaccard_h0` of the same pair (checked numerically to $10^{-16}$). The two versions are
equal when no base is covered twice by either walk.

The CoRAL code differs from the supplement because it counts multiplicity. The code's version is
the headline because it produced CoRAL's published numbers. In CoRAL's code the true side is
every true cycle of the amplicon (`true_cycles_bed`, all `cycle_id`s), and the minimum distance
over the reconstructed cycles is kept. Here the score is per truth species against its matched
cycle (strict mode), or the best over candidates (best-of). CoRAL reads a cycles file as
half-open BED intervals $[\mathrm{start}, \mathrm{end})$ of AA's 1-based inclusive coordinates,
one base short per run (`coral/scoring/io_utils.py::read_cycles_file_to_bed`); here steps are
inclusive.

#### CoRAL's cyclic longest common subsequence: `cyclic_lcs_coral_code` (headline) and `cyclic_lcs_supp` (secondary)

> **Provenance:** (a) CoRAL [Zhu et al. 2024](references.md#zhu2024). The code version is CoRAL 3.1.1
> `coral/scoring/scoring_utils.py::compute_normalized_longest_cycle_subsequence` with
> `::find_longest_subsequence` and `::repartition_intervals` [CoRAL 3.1.1 code](references.md#coral_code). The LCS itself is the classical dynamic
> programme [Wagner & Fischer 1974](references.md#wagnerfischer1974). **Code:** `rck/ecdna/aa_eval.py::cyclic_lcs`,
> `::_weighted_lcs`, `::_labels`.

The bins $\mathrm{Bins}$ and multiplicities $m_X$ are those of the overlap, here built from the
runs of $T^0$ and $C^0$. Let
$\mathrm{Bins}_\cap = \lbrace \Delta \in \mathrm{Bins} : m_T(\Delta) > 0,\ m_C(\Delta) > 0 \rbrace$ be the bins both walks
cover. Write each walk $X$ as the sequence $\mathrm{seq}_X$ of labels $(\Delta, \mathrm{str}_s)$: its
steps in traversal order, each listing its bins of $\mathrm{Bins}_\cap$ in increasing position, or in
reverse when $\mathrm{str}_s = -$. Then:

```math
\mathrm{LCS}^\ast = \max_{\mathrm{rot} \in \mathrm{Rot}(\mathrm{seq}_C)} \mathrm{WLCS}\big(\mathrm{seq}_T,\ \mathrm{rot}\big),
\qquad
\mathrm{cyclic\_lcs\_supp} = \frac{\mathrm{LCS}^\ast}{\lvert T \rvert},
\qquad
\mathrm{cyclic\_lcs\_coral\_code} = \frac{\mathrm{LCS}^\ast}{\sum_{\Delta \in \mathrm{Bins}} \max\big(m_T(\Delta),\ m_C(\Delta)\big)\, \lvert \Delta \rvert}
```

where

- $\mathrm{WLCS}(x, y)$ is the longest common subsequence of label sequences $x$ and $y$, in which
  a match of two equal labels $(\Delta, \mathrm{str})$ scores $\lvert \Delta \rvert$ (dynamic programme
  $\mathrm{DP}_{i,j} = \mathrm{DP}_{i-1,j-1} + \lvert \Delta_i \rvert$ if $x_i = y_j$, else $\max(\mathrm{DP}_{i-1,j}, \mathrm{DP}_{i,j-1})$);
- $\mathrm{Rot}(\mathrm{seq}_C)$ is every rotation of $\mathrm{seq}_C$ and of its reverse, the
  reverse with every orientation flipped. For a linear $C$ it is $\mathrm{seq}_C$ and its
  reverse, unrotated. $\mathrm{seq}_T$ is not rotated, as in CoRAL: rotating one side is enough
  for a cyclic common subsequence;
- $\lvert T \rvert$ is the truth's length with multiplicity (§1);
- the code's denominator sums over *all* joint bins, including those only one walk covers.

Both lie in $[0, 1]$. Both are 0 when $\mathrm{Bins}_\cap$ is empty and for a missed species.

The two normalisations differ: the code normalises by the bins' total with multiplicity, the
supplement by the true cycle's length. The code's version is the headline. It produced CoRAL's
published numbers, and it charges a merged walk for its extra bases. The supplement's version
gives 1 to a merged walk that contains the true cycle, which flatters exactly the nested-species
failure the bench tests. **Effect** at λ = 0.2 on the grown bench: code version 0.871 strict and
0.958 best-of; supplement version 0.901 strict and 0.986 best-of (RCK `CHANGELOG.md`, Changed;
`results/bench_nested/README.md`).

Differences from CoRAL's code, besides the per-species true side (as for the overlap; CoRAL's
true string is every true cycle's rows concatenated in file order):

- CoRAL rotates every reconstructed walk, linear ones too; this does not rotate a linear $C$;
- CoRAL builds runs by merging exactly adjacent segments (`io_utils.py::read_cycles_file_to_bed`);
  this uses §1's runs, with the 2 bp slack and the haplotype requirement (moot on stripped walks).

#### `length_error`

> **Provenance:** (a) CoRAL [Zhu et al. 2024](references.md#zhu2024). The ratio $R_{\mathrm{len}}$ is in CoRAL 3.1.1's code: the second
> return value of `get_fragments_similarity_unweighted`, stored as `reconstruction_length_ratio`
> (`score_simulation.py::score_reconstruction`).
> The $\log_2(1 + R_{\mathrm{len}})$ form follows the supplementary methods, as recorded in the previous version
> of this documentation. **Code:** `rck/ecdna/aa_eval.py::length_error`.

```math
R_{\mathrm{len}} = \frac{\lvert C \rvert - \lvert T \rvert}{\lvert T \rvert},
\qquad
\mathrm{length\_error} = \log_2(1 + R_{\mathrm{len}}) = \log_2 \frac{\lvert C \rvert}{\lvert T \rvert}
```

$\lvert \cdot \rvert$ is the walk length with multiplicity (§1), on the steps of the haplotype-blind
partner. The error is signed, in $(-\infty, \infty)$, and is 0 at the right length. It is empty
for a missed species. Best-of takes the value of smallest absolute size. CoRAL's code instead
reports the $R_{\mathrm{len}}$ of the reconstructed cycle with the best overlap, against the
length of every true cycle of the amplicon.

#### `exact_cycle`, `exact_cycle_phased`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409, after AA's and CoRAL's exact-reconstruction
> criterion [Deshpande et al. 2019](references.md#deshpande2019) [Zhu et al. 2024](references.md#zhu2024). **Code:** `rck/ecdna/aa_eval.py::_exact`.

```math
\mathrm{exact\_cycle}(T, C) = \mathrm{order\_match}_{\tau_{\mathrm{gr}}}(T^{0}, C^{0}),
\qquad
\mathrm{exact\_cycle\_phased}(T, C) = \mathrm{order\_match}_{\tau_{\mathrm{gr}}}(T, C)
```

This is §3.4 at the graph tolerance $\tau_{\mathrm{gr}}$ (`--graph-tol`, 100 bp), with `0`/`1`
values: the truth's first ordering against every ordering of $C$, topologies equal, on
haplotype-stripped walks or with haplotypes. Every run start and end lies within 100 bp, so every
novel adjacency does too. `exact_cycle` scores the haplotype-blind partner (`inferred_h0`) and
`exact_cycle_phased` the phased one (`inferred`). Both are 0 for a missed species;
`exact_cycle_phased` is empty for an unphased answer.

### 4.5 Sample-level metrics (`aa_sample_metrics.tsv`)

> **Code:** `rck/ecdna/aa_eval.py::score` (sample rows), `::cn_explained`, `::subwalks_satisfied`.

There is one row per mode. $\mathcal S$ is the strict answer, $\mathcal W^\ast_{\mathrm{in}}$ the
in-scope truth species, and $\mathcal C(T)$ the candidates of $T$ (§4.6). Means are over the
mode's rows, one per $T \in \mathcal W^\ast_{\mathrm{in}}$, skipping empty cells; with no
in-scope truth they are empty.

| Column | Definition | Provenance |
|---|---|---|
| `mode`, `strict_rule`, `phased` | the mode, the rule applied (§4.6, `auto` resolved) and whether the answer is phased (§4.4) | (c) new |
| `n_truth`, `n_inferred` | $\lvert \mathcal W^\ast_{\mathrm{in}} \rvert$; $\lvert \mathcal S \rvert$ (strict) or the number of candidate walks of `--cycles` and `--alt-cycles` (best-of) | (c) new |
| `species_count_error` | $\lvert \mathcal S \rvert - \lvert \mathcal W^\ast_{\mathrm{in}} \rvert$ (strict only) | (c) new |
| `n_missed` | strict: truth species without a partner in the matching (phased, or haplotype-blind for an unphased answer); best-of: truth species with `composite_h0` $\le 0$, that is with no candidate of haplotype-blind Jaccard $> 0$ | (c) new |
| `n_false` | strict species without a partner in the same matching (strict only) | (c) new |
| `mean_composite`, `mean_composite_h0`, `mean_cycle_overlap_coral_code`, `mean_cyclic_lcs_coral_code`, `mean_cycle_overlap_supp`, `mean_cyclic_lcs_supp` | means over $\mathcal W^\ast_{\mathrm{in}}$ of the species columns (`mean_composite` for phased answers only) | (c) new |
| `n_exact_cycles` | $\sum_T$ `exact_cycle` (haplotype-blind) | (c) new |
| `cn_explained_top1`, `cn_explained_top3` | below; the same value in both modes | (a) CoRAL, modified |
| `subwalks_satisfied`, `subwalks_satisfied_support`, `n_subwalks` | below | (a) CoRAL, modified |
| `phase_frame`, `n_missed_absolute`, `n_false_absolute`, `n_blocks_swapped` | with `--phase-blocks` only (§3.10) | (c) new |
| `error_free`, `error_free_merged` | 1 only if every scored amplicon passes | (c) aggregation of §4.3 |
| `bp_accuracy`, `novel_adjacency_fnr`, `novel_adjacency_fdr`, `segment_cn_mae` | means over the scored amplicons | (c) aggregation of §4.3 |
| `is_violations` | sum over the scored amplicons | (c) aggregation of §4.3 |

The §4.3 columns are present only with `--graph` and a truth mixture graph.

#### `cn_explained_topk`: copy number explained by the heaviest cycles

> **Provenance:** (a) CoRAL's copy-number ratios `best_copy_number_ratio` and
> `top_three_copy_number_ratio` [Zhu et al. 2024](references.md#zhu2024), CoRAL 3.1.1
> `coral/scoring/scoring_utils.py::get_cycle_copy_number_ratio`, called with `n` = 1 and 3 from
> `coral/scoring/score_simulation.py::score_reconstruction` [CoRAL 3.1.1 code](references.md#coral_code). (c) modified
> as described below.

```math
\mathrm{cn\_explained\_top}k = \frac{\sum_{C \in \mathrm{top}_k} w_C\, \lvert C \rvert}{\sum_{j \in \mathcal E_{\mathrm{seq}},\ \mathrm{CN}(j) \ge \mathrm{cn}_{\min}} \mathrm{CN}(j)\, L(j)},
\qquad k \in \lbrace 1, 3 \rbrace
```

where

- $\mathrm{top}_k$ is the $k$ circular walks with the largest $w_C \lvert C \rvert$ (copy count ×
  length with multiplicity), among every walk of `--cycles` except `Status=chromosomal`,
  `low_weight` included, `--alt-cycles` not used. Paths are left out, as CoRAL scores cyclic
  walks only. With fewer than $k$ circular walks, all of them;
- $\mathcal E_{\mathrm{seq}}$ is the `sequence` lines of the method's own `--graph`, $j$ one of
  them (AA/CoRAL *sequence edge*), $\mathrm{CN}(j)$ its `PredictedCopyCount`, and
  $L(j) = \mathrm{end}_j - \mathrm{start}_j + 1$ its length in AA's inclusive coordinates;
- $\mathrm{cn}_{\min}$ = `--cn-explained-min-cn` (5, CoRAL's constant).

No truth is needed. It is empty without `--graph` or when the denominator is 0. It is not capped:
a walk through segments below $\mathrm{cn}_{\min}$ can push it above 1.

Differences from CoRAL's code:

- CoRAL returns 1 when the number of cyclic walks is at most $k$
  (`if n >= len(cycle_to_ratio): return 1.0`); this computes the share whatever their number;
- CoRAL's lengths are $\mathrm{end} - \mathrm{start}$ of AA's inclusive coordinates, one base
  short per run and per sequence edge; these are inclusive;
- CoRAL computes it per reconstructed amplicon and reports the one of the amplicon whose cycles
  have the best overlap; this uses the whole cycles file and graph.

#### `subwalks_satisfied`: long-read subwalk constraints

> **Provenance:** (a) CoRAL's long-read subwalk (path) constraints [Zhu et al. 2024](references.md#zhu2024). The
> satisfaction share is (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_eval.py::read_paths`,
> `::_reference_join`, `::_cycle_junctions`, `::path_in_cycle`, `::subwalks_satisfied`.

The paths come from `--paths`, a subgraph `paths.tsv` (`PID`, `PATH`, `NSTEP`, `#NOV`). `PATH` is
a comma-separated list of extremity tokens `chr:pos±@orient`; `pos` is in the walks' 0-based
frame, `+` the right end and `−` the left end of a segment, and `@orient` is ignored. With tokens
$v_0, \dots, v_{2n-1}$ (0-based), each pair $(v_{2i-1}, v_{2i})$, $i = 1, \dots, n-1$, is one join.

1. A join $(v, v')$ is a *reference adjacency* when $v$ and $v'$ are on the same chromosome, on
   opposite strands, and the `−` position equals the `+` position + 1 within 2 bp. Every other
   join is a novel adjacency.
2. A read path $\upsilon$ (ecDNA.md's notation) is kept when it has at least two novel adjacencies,
   counted with repeats. It becomes the sequence $\mathcal A_N(\upsilon) = (a_1, \dots, a_m)$ of its
   novel adjacencies, each $a = (v, v')$ ordered as read. This matches `rck-ecdna`'s
   `--min-path-junctions` (2); `rck-ecdna` identifies reference adjacencies on its own graph
   instead (`walk.py::read_paths`) and also requires support $\ge 1$.
3. A walk $C$ gives the sequence $\mathcal A_N(C) = (a'_1, \dots, a'_n)$ of its novel
   adjacencies (exit, entry) in traversal order (`evaluate.py::junctions` on $C^0$, so
   haplotype-blind; the closing join only for a circle).
4. $\upsilon \preceq C$ ($C$ satisfies $\upsilon$) when, for some start $i$, every
   $a_k$ matches $a'_{i+k}$ end by end in the same order (exit with exit, entry with entry), both
   extremities within $\tau_{\mathrm{gr}}$ (§1). Indices are taken mod $n$ for a circle, so the
   read may wrap; for a linear walk $i + m \le n$. The same test is run on the reversed read,
   $\big((v'_m, v_m), \dots, (v'_1, v_1)\big)$. A walk without novel adjacencies satisfies nothing.

```math
\mathrm{subwalks\_satisfied} = \frac{\lvert \lbrace \upsilon \in \Upsilon : \exists C \in \mathcal X,\ \upsilon \preceq C \rbrace \rvert}{\lvert \Upsilon \rvert},
\qquad
\mathrm{subwalks\_satisfied\_support} = \frac{\sum_{\upsilon\ \text{satisfied}} N_\upsilon}{\sum_{\upsilon \in \Upsilon} N_\upsilon}
```

where

- $\mathcal X$ is the mode's walks: $\mathcal S$ in strict mode, every candidate walk in best-of;
- $\Upsilon$ is the kept read paths, and `n_subwalks` $= \lvert \Upsilon \rvert$;
- $N_\upsilon$ is the read path's support, its `#NOV` (a 0 or empty value is read as 1).

Both lie in $[0, 1]$. Both are empty when $\Upsilon$ is empty, including without `--paths`.

### 4.6 Modes: `strict` and `best_of`, and `--strict-rule`

> **Provenance:** `strict` is (c) new, RCK yubi_dev@e4c2409. `best_of` is (a) CoRAL's scoring,
> which keeps each metric's best value over the reconstructed cycles independently
> (`coral/scoring/score_simulation.py::score_reconstruction`: fragment overlap and LCS are
> maximised separately) [CoRAL 3.1.1 code](references.md#coral_code), applied per truth species here. **Code:**
> `rck/ecdna/aa_eval.py::species_cycles`, `::strict_select`, `::strict_rows`, `::best_of_rows`.

Both modes start from the method's walks: every `Cycle=` line of `--cycles` except those flagged
`Status=chromosomal` (`species_cycles`), each with its copy count $w_C$ (`Copy_count`).

**Strict.** This is the method's own answer $\mathcal S \subseteq$ its walks, chosen by
`--strict-rule`:

| Rule | $\mathcal S$ |
|---|---|
| `status` | the walks without `Status=low_weight`; for `rck-ecdna`, whose `--cycles` is rank 1, rank 1 without `low_weight` |
| `cn-frac:`$f_{\mathrm{cn}}$ | $\lbrace C : w_C \ge f_{\mathrm{cn}} \cdot \max_{C'} w_{C'} \rbrace$ |
| `top-cn:`$f_{\mathrm{mass}}$ | the walks in decreasing order of $w_C \lvert C \rvert$, taken until their cumulative mass reaches $f_{\mathrm{mass}} \sum_C w_C \lvert C \rvert$ (every walk if the total is 0) |
| `all` | every walk |
| `auto` (default) | `status` when any `Cycle=` line carries `Status=`, else `cn-frac:0.1` |

Linear walks count in every rule. AA's and CoRAL's cycles files list many low-copy cycles that the
methods themselves do not report as species. `cn-frac:0.1` approximates `rck-ecdna`'s own
significance rule (`--min-species-weight-frac` 0.1: with more than one species, a species is
`low_weight` when $w < \max(1,\ 0.1 \max w)$, `decompose.py`), so the three methods face a
comparable threshold. It is not identical: see §7. Strict mode matches one to one as in §3.6, phased
and haplotype-blind (§4.4), and counts missed and false species.

**Best-of.** For each truth species $T$, the candidates $\mathcal C(T)$ are every walk of
`--cycles` and of `--alt-cycles` (`Status=chromosomal` left out, `low_weight` included) that
overlaps $T$: haplotype Jaccard $> 0$ for the phased columns (phased answers), haplotype-blind
Jaccard $> 0$ for the haplotype-blind columns (and, for an unphased answer, for the other §3
columns). Each metric takes its best value on its own:

```math
\mathrm{score}^{\mathrm{best}}(T) = \max_{C \in \mathcal C(T)} \mathrm{score}(T, C)
\quad\text{(min for cn\_rel\_err, cn\_rel\_err\_h0; smallest}\ \lvert\cdot\rvert\ \text{for length\_error)}
```

Booleans count as 0/1, and empty values are skipped. With no candidate, `composite_h0`,
`exact_cycle` and the four CoRAL columns are 0, `composite` and `exact_cycle_phased` are 0 for a
phased answer, and every other column is empty.

`inferred`, `inferred_weight` and `inferred_topology` name the candidate with the best
`composite` (for an unphased answer, the best `composite_h0`); `inferred_h0` names the candidate
with the best `composite_h0`. Ties go to the first candidate: `--cycles` in file order, then each
`--alt-cycles` file in order (names prefixed `alt1.`, `alt2.`, ...). For `rck-ecdna` the
alternatives are rank 1's first $k$ same-fit alternatives (`aa_export --alternatives`, $k = 2$ by
default), the candidates of the set metric (§6.2). Best-of picks per truth species, not per
decomposition, so it is never below that set metric.

CoRAL's own best-of differs: it does not optimise the length ratio or the copy-number ratios on
their own but reports those of the best-overlap cycle, and its LCS keeps the cycle with the best
normalised LCS within an amplicon and the best raw LCS across amplicons.

**Why two modes.** Strict mode is what a user of the method gets. Best-of is what CoRAL's own
benchmark reports, so the numbers are comparable with its paper, and for `rck-ecdna` it measures
whether the ambiguity report contains the right answer. **Effect** at λ = 0.2 on the grown bench:
best-of mean composite 0.931 against 0.853 strict; 55 of 72 truth species are exact cycles in
strict mode and 61 in best-of (RCK `results/bench_nested/README.md`).
### 4.7 Confusion-matrix metrics: species detection (M1) and adjacency pairs (M2)

> **Provenance:** (c) new, RCK yubi_dev, 2026-11-07. The building blocks are (a): the
> confusion-matrix rates and the F-measure, precision-recall reporting without true negatives
> [Krusche et al. 2019](references.md#krusche2019), pair counting [Rand 1971](references.md#rand1971) generalised to overlapping
> clusters [Collins & Dent 1988](references.md#collins1988). Their combination for ecDNA species is (c). **Code:**
> `rck/ecdna/aa_eval.py::partition_metrics`, `::detection_rows`, `::pair_rows`, `::walk_adjacencies`,
> `::adjacency_items`, `::junction_distance`, `::pair_confusion`, `::confusion_rates`, `::pool_partition`.
> Output: `aa_partition_metrics.tsv` (Appendix B).

**Why these metrics.** The goal is to reconstruct each ecDNA molecule on its own. Two errors defeat
that goal while leaving most other scores high:

- a *merge*: one inferred walk carries two molecules;
- a *split*: one molecule is reported as several walks.

The composite (§3.7) charges them only in part. A merged walk is matched to one of the two truth
species and loses some segment Jaccard, junction precision and its order match. The other truth
species is missed, but the merged walk is never counted as a false species. The bench's own
`merged` flag (`results/bench_nested/collect.py::merged_flag`) detects merges but gives no rate.
Two metrics are defined here:

- **M1** asks, for each molecule: *was it reconstructed?*
- **M2** asks, for each pair of novel adjacencies: *did the answer put them in the same molecule
  exactly when the truth does?*

Both report the counts of true positives (TP), false positives (FP), false negatives (FN) and,
where defined, true negatives (TN), with these rates (`confusion_rates`):

```math
\mathrm{precision} = \frac{TP}{TP + FP},
\qquad
\mathrm{recall} = \frac{TP}{TP + FN},
\qquad
\mathrm{FPR} = \frac{FP}{FP + TN},
\qquad
F_1 = \frac{2\,TP}{2\,TP + FP + FN}
```

All four lie in $[0, 1]$. Conventions:

- $F_1$ is written in the form above, which equals the harmonic mean of precision and recall and is
  Dice's coefficient [Dice 1945](references.md#dice1945). It is defined whenever anything is positive, even
  when precision or recall is not.
- A rate whose denominator is 0, or that needs a count that is undefined (FP or TN written empty), is
  *undefined*: written empty (NaN), never 1. `junction_scores` (§3.3) returns 1 for an empty set;
  that convention would score an empty answer as perfect, so it is not used here.

Each row of `aa_partition_metrics.tsv` has the fields `metric`, `mode`, `haplotype`, `criterion`,
`tol_bp`, `universe`, `headline`, `n_truth`, `n_predicted`, `n_items`, `tp`, `fp`, `fn`, `tn`,
`precision`, `recall`, `fpr`, `f1`. An undefined count or field is empty.

Notation used in both metrics:

| Symbol | Meaning |
|---|---|
| $\mathcal W^\ast_{\mathrm{in}}$ | the in-scope truth species (§2.3); $T \in \mathcal W^\ast_{\mathrm{in}}$ |
| $\mathcal S$ | the method's strict answer (§4.6: `status`, else `cn-frac:0.1`); $C \in \mathcal S$. `Status=chromosomal` paths are never species (`species_cycles`) |
| $X^0$ | walk $X$ with every haplotype set to `.` (`strip_haplotypes`, §4.4) |
| $\mathrm{match}(T)$, $\mathrm{match}^0(T)$ | $T$'s partner in the phased (haplotype Jaccard), resp. haplotype-blind (Jaccard of $T^0$, $C^0$), one-to-one matching (§3.6, §4.4), or $\varnothing$ |
| $\tau_{\mathrm{sp}}$, $\tau_{\mathrm{gr}}$ | the species and graph tolerances, `--species-tol` 10 000 bp and `--graph-tol` 100 bp (§1) |
| $v = (\mathrm{chr}, \mathrm{pos}, \pm)$ | an extremity (§1) |
| $a = (v_1, v_2, H_1, H_2)$ | a novel adjacency of a walk, as `evaluate.junctions` lists it (§1): exit extremity $v_1$ of step $s$, entry extremity $v_2$ of the next step $s'$, and $(H_1, H_2) = (H_s, H_{s'})$ |
| $\mathcal A_N(X)$ | the distinct novel adjacencies of walk $X$ (`walk_adjacencies`): those of $X^0$ for `blind`, of $X$ for `phased`; duplicates at distance 0 (below) removed; no wrap-around adjacency for a linear walk |

A *novel adjacency* is the code's "junction" (§1): RCK's novel adjacency, AA's and CoRAL's
discordant edge [Aganezov & Raphael 2020](references.md#aganezov2020) [Deshpande et al. 2019](references.md#deshpande2019) [Zhu et al. 2024](references.md#zhu2024).

#### M1. Species detection (`metric` = `species_detection`)

**Unit and classes.** The unit is a species. A truth species is *positive in the truth*; a strict
species is *positive in the prediction*. For a criterion with indicator
$\mathrm{pass}(T, C) \in \lbrace 0, 1 \rbrace$:

```math
TP = \big\lvert \lbrace T \in \mathcal W^\ast_{\mathrm{in}} : \mathrm{match}(T) \ne \varnothing,\ \mathrm{pass}\big(T, \mathrm{match}(T)\big) = 1 \rbrace \big\rvert,
\qquad
FN = \lvert \mathcal W^\ast_{\mathrm{in}} \rvert - TP,
\qquad
FP = \lvert \mathcal S \rvert - TP
```

```math
\mathrm{precision} = \frac{TP}{\lvert \mathcal S \rvert},
\qquad
\mathrm{recall} = \frac{TP}{\lvert \mathcal W^\ast_{\mathrm{in}} \rvert},
\qquad
F_1 = \frac{2\,TP}{2\,TP + FP + FN} = \frac{2\,TP}{\lvert \mathcal W^\ast_{\mathrm{in}} \rvert + \lvert \mathcal S \rvert}
```

Precision is undefined for $\mathcal S = \varnothing$, recall for $\mathcal W^\ast_{\mathrm{in}} = \varnothing$,
$F_1$ when both are empty. TN and FPR are **undefined** (`tn` empty). Fields: `n_truth` =
$\lvert \mathcal W^\ast_{\mathrm{in}} \rvert$, `n_predicted` = $\lvert \mathcal S \rvert$,
`universe` = `species`, `n_items` empty.

$TP$ is read from the strict rows of `aa_species_metrics.tsv` (§4.4): one row per truth species,
with its matched partner's scores, and 0 for a missed species. The criteria (`match` is the blind
matching $\mathrm{match}^0$ for the blind rows):

| `haplotype` | `criterion` | $\mathrm{pass}(T, C) = 1$ when | tolerance (`tol_bp`) |
|---|---|---|---|
| `blind` | `exact` (**headline**) | `exact_cycle` $= 1$ (§4.4): same run sequence, haplotypes ignored | $\tau_{\mathrm{gr}}$ |
| `phased` | `exact` | `exact_cycle_phased` $= 1$: the same with haplotypes | $\tau_{\mathrm{gr}}$ |
| `blind` | `composite_h0>=`$\mathrm{thr}$ | `composite_h0` $\ge \mathrm{thr}$ (§4.4) | $\tau_{\mathrm{sp}}$ |
| `phased` | `composite>=`$\mathrm{thr}$ | `composite` $\ge \mathrm{thr}$ (§3.7) | $\tau_{\mathrm{sp}}$ |

for $\mathrm{thr} \in \mathrm{Thr} = \lbrace 0.5, 0.6, 0.7, 0.8, 0.9, 0.95 \rbrace$ (`DETECTION_TAUS`), written
with `%g` (`composite_h0>=0.95`). An empty score counts as 0. Phased rows are written only for a
phased answer, as for `composite`. The `headline` field is true for both strict `exact` rows,
`blind` and `phased`; methods are compared on `blind`.

**Best-of** (`mode` = `best_of`). $\mathcal C(T)$ is $T$'s candidates (§4.6): the cycles of
`--cycles` and `--alt-cycles`, `low_weight` included, that overlap $T$ (haplotype Jaccard $> 0$
for the phased criteria, blind Jaccard $> 0$ for the blind ones). The TP count is the number of
truth species with *some* candidate that passes:

```math
TP^{\mathrm{best}} = \big\lvert \lbrace T \in \mathcal W^\ast_{\mathrm{in}} : \exists\, C \in \mathcal C(T),\ \mathrm{pass}(T, C) = 1 \rbrace \big\rvert,
\qquad
\mathrm{recall}^{\mathrm{best}} = \frac{TP^{\mathrm{best}}}{\lvert \mathcal W^\ast_{\mathrm{in}} \rvert}
```

The code reads it from the best-of rows of §4.6, whose `exact_cycle` and composites are already the
maxima over $\mathcal C(T)$. FP, precision and $F_1$ are undefined in best-of (`fp` and
`n_predicted` empty), and so are TN and FPR.

**Worked example** (truth: two molecules O and I):

| Answer | TP | FP | FN | precision | recall | $F_1$ |
|---|---|---|---|---|---|---|
| O and I, both exact | 2 | 0 | 0 | 1 | 1 | 1 |
| one walk carrying O and I (merge) | 0 | 1 | 2 | 0 | 0 | 0 |
| O exact; I reported as two walks (split) | 1 | 2 | 1 | 1/3 | 1/2 | 2/5 |

**Reasoning, decision by decision.**

- **The unit is the species.** The goal is each molecule, so the molecule is the unit. Gene-structure
  benchmarks score at several levels (nucleotide, exon, protein), the exon level counting the exons
  identified exactly [Burset & Guigó 1996](references.md#burset1996). Transcript assemblers are scored on whole transcripts
  whose intron chain matches exactly [Pertea & Pertea 2020](references.md#pertea2020). Both are the
  analogue of a species: a whole walk through a splice graph or a breakpoint graph.
- **FP $= \lvert\mathcal S\rvert - TP$.** Precision's denominator is every prediction, as in variant
  benchmarking, where every query call that matches no truth call is a false positive
  [Krusche et al. 2019](references.md#krusche2019). So a matched partner that fails the criterion is a false
  positive as well as a false negative. Counting only *unmatched* inferred species as FP would score
  a merge as FP 0 and precision 1: the merged walk is matched, so nothing would be false.
- **The matching is the existing one (§3.6).** Then $TP$ at $\mathrm{thr}$ is the number of strict
  species rows (for `rck-ecdna`, the recovery-matrix rows, §4.8) with composite $\ge \mathrm{thr}$. It
  cannot increase with $\mathrm{thr}$, so the precision-recall points form a curve. An exact pair has
  Jaccard close to 1, so the greedy matching takes it first.
- **Exact is the headline.** The order match at $\tau_{\mathrm{gr}}$ is AA's and CoRAL's criterion for a correct
  cycle [Deshpande et al. 2019](references.md#deshpande2019) [Zhu et al. 2024](references.md#zhu2024). It ignores the copy number,
  which the composite does not. With the composite, a walk of the right structure but 30% copy-number
  error falls below 0.95.
- **The composite thresholds are a curve, not one number.** No published threshold exists for a
  partly correct ecDNA reconstruction [citation needed: a published threshold for a partly correct
  ecDNA reconstruction]. The full precision-recall curve is reported instead, which compares methods
  without choosing a threshold [Davis & Goadrich 2006](references.md#davis2006). The highest point, 0.95, is
  the bench's own failure threshold (`results/bench_nested/collect.py::FAIL_COMPOSITE`); it is a
  project choice, not a published value.
- **No TN, no FPR.** A true negative would be a walk that is not a truth species and is not reported.
  The walks of a breakpoint graph cannot be listed, and with multiplicity there are infinitely many,
  so TN cannot be counted. Variant-calling benchmarks report precision and recall instead of
  specificity for the same reason [Krusche et al. 2019](references.md#krusche2019); precision-recall is also the
  more informative summary when negatives far outnumber positives [Saito & Rehmsmeier 2015](references.md#saito2015).
- **Best-of has recall only.** Best-of picks, per truth species, the best of every cycle the method
  returned (CoRAL's convention, §4.6); those candidates are not a set the method claims, so there is
  nothing to compute a precision on.
- **Pooling over cases: micro (headline).** Counts are summed over cases, then the rates are
  computed, so every molecule weighs the same [Sokolova & Lapalme 2009](references.md#sokolova2009). Cases have
  1 to 3 in-scope species, so a per-case $F_1$ takes few values. Micro pooling also gives more weight
  to the multi-species cases, where merges and splits happen. The macro mean (the mean of the per-case
  rates) is reported beside it (§6.3).

#### M2. Pairwise co-membership of novel adjacencies (`metric` = `adjacency_pairs`)

**Distance between novel adjacencies** (`junction_distance`). For $a = (v_1, v_2, H_1, H_2)$ and
$a' = (v'_1, v'_2, H'_1, H'_2)$, a *direction* pairs $v_1$ with $v'_1$ and $v_2$ with $v'_2$
(forward), or $v_1$ with $v'_2$ and $v_2$ with $v'_1$ (reverse). A direction is admissible when the
paired extremities agree in chromosome and strand and, for `phased`, when the labelling agrees in
that direction: $(H'_1, H'_2) = (H_1, H_2)$ forward, $(H'_2, H'_1) = (H_1, H_2)$ reverse. Then

```math
\mathrm{dist}(a, a') = \min\Big( \underbrace{\max\big(\lvert \mathrm{pos}_1 - \mathrm{pos}'_1 \rvert,\ \lvert \mathrm{pos}_2 - \mathrm{pos}'_2 \rvert\big)}_{\text{forward}},\ \underbrace{\max\big(\lvert \mathrm{pos}_1 - \mathrm{pos}'_2 \rvert,\ \lvert \mathrm{pos}_2 - \mathrm{pos}'_1 \rvert\big)}_{\text{reverse}} \Big)
```

over the admissible directions, where $\mathrm{pos}_k$ is the position of $v_k$. With no admissible
direction, $\mathrm{dist} = \infty$ (code `None`). Two adjacencies are *equal* when
$\mathrm{dist} = 0$, and *match within $\tau$* when $\mathrm{dist} \le \tau$.

**Items.** The unit is an unordered pair of distinct *items*. An item is a novel adjacency
(`haplotype` = `blind`, on haplotype-stripped walks), or a novel adjacency with its haplotype
labelling $(H_1, H_2)$ (`phased`): AA, AB, BA or BB, RCK's labelling $\varphi$
[Aganezov & Raphael 2020](references.md#aganezov2020), or a pair with `.` where a step is unphased. The universe
$\mathcal I$ of items is built per tolerance $\tau \in \lbrace \tau_{\mathrm{sp}}, \tau_{\mathrm{gr}} \rbrace$
(`adjacency_items`):

1. **Truth items:** $\bigcup_{T \in \mathcal W^\ast_{\mathrm{in}}} \mathcal A_N(T)$, with equal adjacencies
   merged. A novel adjacency that two truth species share exactly is one item, so two distinct truth
   adjacencies are never one item, at any tolerance.
2. **Predicted adjacencies** $\bigcup_{C \in \mathcal S} \mathcal A_N(C)$, equal ones merged over the whole
   answer, are matched **one to one** to truth items, nearest first, within $\tau$:
   1. list every (predicted adjacency $a$, truth item $a'$) with $\mathrm{dist}(a, a') \le \tau$;
   2. sort the list by $\mathrm{dist}$, then by the index of $a$, then by the index of $a'$. Indices
      follow first appearance: species in input order (cycles-file order for $\mathcal S$), each
      walk's adjacencies in walk order;
   3. accept each entry whose $a$ and $a'$ are both still free.

   A predicted adjacency left unmatched becomes a new, *false* item.
3. **Universe** (`universe`):
   - `native`: $\mathcal I$ = items 1 and 2;
   - `graph` (blind only, when `--graph` is given): $\mathcal I$ also holds every discordant edge of the
     graph (`graph_adjacencies`, AA positions minus 1, labels `.`) that matches no item already in
     $\mathcal I$ within $\tau$ (not one to one; edges added earlier count). In the graph route, every
     method decomposes the same graph, so they all have the same set of negative pairs.

Each species becomes the set of its items: $\mathcal I(T) \subseteq \mathcal I$ and
$\mathcal I(C) \subseteq \mathcal I$.

**Classes.** For a pair $\lbrace a, a' \rbrace$ of distinct items $a \ne a' \in \mathcal I$, *together*
in the truth and in the answer are:

```math
\mathrm{tog}^\ast(a, a') = \mathbb 1\big[\exists\, T \in \mathcal W^\ast_{\mathrm{in}} : a, a' \in \mathcal I(T)\big],
\qquad
\widehat{\mathrm{tog}}(a, a') = \mathbb 1\big[\exists\, C \in \mathcal S : a, a' \in \mathcal I(C)\big]
```

```math
TP = \sum_{\lbrace a, a' \rbrace} \mathrm{tog}^\ast\, \widehat{\mathrm{tog}},
\quad
FP = \sum_{\lbrace a, a' \rbrace} (1 - \mathrm{tog}^\ast)\, \widehat{\mathrm{tog}},
\quad
FN = \sum_{\lbrace a, a' \rbrace} \mathrm{tog}^\ast\, (1 - \widehat{\mathrm{tog}}),
\quad
TN = \binom{\lvert \mathcal I \rvert}{2} - TP - FP - FN
```

(`pair_confusion`). Fields: `criterion` = `co_membership`, `mode` = `strict`, `n_truth` =
$\lvert \mathcal W^\ast_{\mathrm{in}} \rvert$, `n_predicted` = $\lvert \mathcal S \rvert$, `n_items` =
$\lvert \mathcal I \rvert$. Rows: `blind`, and `phased` for a phased answer; at `tol_bp` =
$\tau_{\mathrm{sp}}$ and $\tau_{\mathrm{gr}}$ (one set of rows if they are equal); `native`, and `graph`
for `blind`. `headline` is true at $\tau_{\mathrm{sp}}$ in the `native` universe, for `blind` and `phased`.

Precision, recall, FPR and $F_1$ are as above. Their meaning here:

- **FPR is the merge rate:** the share of the pairs that the truth keeps apart which the answer puts
  together.
- **$1 - \mathrm{recall}$ is the split rate:** the share of the pairs that the truth keeps together
  which the answer separates.
- A false item contributes false-positive pairs with the items of its walk, and a missed truth item
  contributes false-negative pairs.

**When the rates are undefined.**

| Rate | Undefined when | Example |
|---|---|---|
| recall | no truth pair, $TP + FN = 0$ | every truth species has at most one novel adjacency |
| precision | no predicted pair, $TP + FP = 0$ | every inferred species has at most one novel adjacency |
| FPR | no negative pair, $FP + TN = 0$ | one truth species and no false items (the bench's `ctrl_dup` cases) |
| $F_1$ | nothing positive, $TP + FP + FN = 0$ | |

**Worked example** (`tests/ecdna/test_aa_partition.py`). O has the novel adjacencies $a_1, a_2$, I
has $a_3, a_4$, and they share reference segments only. $\mathcal I = \lbrace a_1, a_2, a_3, a_4 \rbrace$
gives 6 pairs: 2 positive in the truth ($a_1 a_2$, $a_3 a_4$) and 4 negative.

| Answer | TP | FP | FN | TN | precision | recall | FPR | $F_1$ |
|---|---|---|---|---|---|---|---|---|
| O and I | 2 | 0 | 0 | 4 | 1 | 1 | 0 | 1 |
| one walk, the Euler union of O and I (merge) | 2 | 4 | 0 | 0 | 1/3 | 1 | 1 | 1/2 |
| truth is the single figure-8 walk; the answer splits it into O and I (split) | 2 | 0 | 4 | 0 | 1 | 1/3 | undefined | 1/2 |

**Reasoning, decision by decision.**

- **Pairs of items, not species.** Pair counting compares two groupings of the same items through the
  pairs each puts together. It yields all four counts, so a TN and an FPR exist
  [Rand 1971](references.md#rand1971) [Hubert & Arabie 1985](references.md#hubert1985). Pair precision and recall are the
  two factors of the Fowlkes-Mallows index [Fowlkes & Mallows 1983](references.md#fowlkes1983). A merge
  creates pairs that the truth keeps apart (FP); a split loses pairs that the truth keeps together
  (FN). M1 tells whether a molecule is right; M2 tells how far its novel adjacencies are grouped
  wrongly.
- **Novel adjacencies, not bases.** A nested molecule lies inside its outer one in base space, so
  every pair of its bases is also a pair of the outer molecule's bases, and positive in the truth. A
  merged walk therefore adds no false pair of bases, and a base-level co-membership cannot see the
  nested merge, which is the failure the bench tests. The novel adjacencies of the two molecules are
  distinct: in all 40 bench cases, truth species share no novel adjacency, at 100 bp or at 10 kb,
  phased or blind (checked 2026-11-07 on the cases' `truth_cycles.tsv`).
- **"Together in at least one species".** A novel adjacency may lie on two species (nested species
  can share one on real data), so the truth's groups can overlap. The relation "some species holds
  both" is defined for overlapping groups, as the Omega index's generalisation of the Rand index to
  non-disjoint clusterings is [Collins & Dent 1988](references.md#collins1988). The Omega index also counts *how
  many* groups hold both. That count differs from the 0/1 relation only when species share an
  adjacency, which no bench case has, so it is not computed.
- **Phased items.** A labelled item separates two traversals of one adjacency on different
  haplotypes (the bench's `ctrl_dup` duplications), and charges a walk on the wrong haplotype. A
  haplotype switch between reference neighbours is a novel adjacency of the phased walk (§1, §3.3),
  as it is in `junction_scores`. It is absent from the blind walk.
- **One-to-one matching, nearest first.** Without it, two predicted adjacencies near one truth
  adjacency would both be credited to it. One-to-one matching of calls to truth is also the
  variant-benchmarking rule [Krusche et al. 2019](references.md#krusche2019).
- **Tolerance: $\tau_{\mathrm{sp}}$ = 10 kb is the headline (`headline` = true), $\tau_{\mathrm{gr}}$ = 100 bp the sensitivity
  rows.** M2 measures grouping. Breakpoint accuracy is measured by `exact_cycle` and `bp_accuracy`.
  At 100 bp, an RCK-ecDNA imputed breakpoint (on the subgraph's 10 kb grid) would be charged twice:
  once as a wrong adjacency, and again in every pair it takes part in. The composite's junction $F_1$
  (§3.3) uses $\tau_{\mathrm{sp}}$ for the same reason. On the bench, only one pair of truth adjacencies matches at
  10 kb but not at 100 bp, and both lie in one species (demo3). Truth items are merged only when
  equal, so the wider tolerance cannot merge them.
- **The `graph` universe changes TN and FPR only.** The extra items lie in no species, so TP, FP, FN,
  precision, recall and $F_1$ are the same in both universes. The FPR is comparable between methods
  only in the graph route, where every method has the same candidate edges; in the native universe
  each method's false items add their own negative pairs. A background SV the method leaves unused
  adds true-negative pairs only. One that it routes into a walk becomes a false item.
- **Strict mode only.** M2 is a property of a whole answer. A best-of answer assembled from different
  decompositions is not a partition of anything. Scoring rank 1's alternative decompositions as whole
  answers is left for later.
- **Pooling over cases: macro (headline).** Pair counts grow with the square of a case's size: demo3
  has 14 phased items (91 pairs; 11 blind items, 55 pairs), a three-adjacency case has 3 pairs.
  Summed counts (micro) would let demo3 outweigh such a case about 18 to 1 (blind) or 30 to 1
  (phased), so the headline is the mean of the per-case rates over the cases where each is defined
  (macro) [Sokolova & Lapalme 2009](references.md#sokolova2009). The micro rates are reported beside it (§6.3).
- **A known weakness.** Amigó et al. compare clustering metrics against formal constraints and argue
  that pair-counting metrics fail one that B-cubed satisfies [Amigó et al. 2009](references.md#amigo2009) (this
  claim was not checked against the paper's text; see [references](references.md#not-verified)). In pair
  counting, an error in a large cluster costs more pairs than the same error in a small one. B-cubed
  gives averaged precision and recall per item, but no TP, FP and TN counts, so it cannot give the FPR
  asked for here; it is not computed.

#### Per-cell species matrix: (J_seg, J_adj) (`aa_species_matrix.tsv`)

> **Provenance:** (c) new, RCK yubi_dev, 2026-11-13. Its two elements reuse (c) pieces defined
> above: J_seg is §3.1's base-pair Jaccard [Jaccard 1912](references.md#jaccard1912), and J_adj is a
> Jaccard index over M2's items. The tuple, the column layout and the merge and split read-out are
> (c). **Code:** `rck/ecdna/species_matrix.py::species_matrix`, `::column_items`, `::item_jaccard`,
> `::pair_jaccard`, `::weighted_jaccard`, `::plot_species_matrix`; written by `aa_eval.score` /
> `::write` for every answer (heatmaps with `aa_eval --matrix-plot`). Bench run:
> `results/bench_nested/diag/species_matrix_bench.py`.

**Why.** M1 and M2 score a whole answer. They say *that* an answer merges or splits, not *which*
truth species went into *which* inferred walk. The matrix shows that per (truth, inferred) pair.

**Layout.**
- **Rows:** the in-scope truth species $\mathcal W^\ast_{\mathrm{in}}$ (§2.3), by truth copy number $w_T$, high
  to low (`truth_rank` 1, 2, ...).
- **Columns** (`column_rank` 1, 2, ...):
  1. the strict species $\mathcal S$ (§4.6), by copy number $w_C$ (`Copy_count`, the refitted weight
     for `rck-ecdna`), high to low;
  2. then the answer's other cycles, by copy number, flagged `low_weight`: every cycle of `--cycles`
     not in $\mathcal S$ (not `--alt-cycles`, not `Status=chromosomal`). For `rck-ecdna` these are
     its `low_weight` species; for AA and CoRAL under `auto`, the cycles under `cn-frac:0.1`.
- Ties are broken by name.
- Rows and columns are both ordered by copy number, so a correct answer is the diagonal.
- One block of rows per `haplotype`: `blind` always, `phased` for a phased answer.

**Cell.** For truth species $T$ and column $C$, $\mathrm{cell}(T, C) = (J_\mathrm{seg}, J_\mathrm{adj})$:

```math
J_\mathrm{seg}(T, C) = \frac{O(T, C)}{\lvert T \rvert + \lvert C \rvert - O(T, C)},
\qquad
J_\mathrm{adj}(T, C) = \frac{\lvert \mathcal I(T) \cap \mathcal I(C) \rvert}{\lvert \mathcal I(T) \cup \mathcal I(C) \rvert}
```

with $O$ the multiplicity-weighted overlap of §3.1, and $\mathcal I(\cdot)$ the item sets below. Both
lie in $[0, 1]$.

- **$J_\mathrm{seg}$** (`j_seg`) is §3.1's Jaccard with multiplicity (`evaluate.jaccard`):
  - `segment_jaccard` on haplotype-stripped walks in the `blind` rows;
  - `haplotype_jaccard` ($O_H$ in place of $O$) in the `phased` rows;
  - 0 for an empty union.
- **$J_\mathrm{adj}$** (`j_adj`) is a Jaccard index over M2's items, at $\tau_{\mathrm{sp}}$ (`tol_bp`):
  - items are novel adjacencies (`blind`), or novel adjacencies with their labelling (`phased`);
  - a linear walk has no wrap-around adjacency;
  - it is NaN for an empty union, never 1, as in M2.
- `exact_cycle` is §4.4 per cell: `exact_cycle` in the `blind` rows and `exact_cycle_phased` in the
  `phased` rows (the column is named `exact_cycle` in both).
- **Optional columns:**
  - `j_adj_pairs`: M2's pair form on one cell,
    $\binom{n_\cap}{2} / \big(\binom{n_T}{2} + \binom{n_C}{2} - \binom{n_\cap}{2}\big)$, with
    $n_\cap = \lvert\mathcal I(T)\cap\mathcal I(C)\rvert$, $n_T = \lvert\mathcal I(T)\rvert$,
    $n_C = \lvert\mathcal I(C)\rvert$; NaN for a zero denominator;
  - `j_seg_cn`: copy-weighted,
    $\sum_x \min(w_T m_T, w_C m_C) / \sum_x \max(w_T m_T, w_C m_C)$, over every base $x$ per
    chromosome (`blind`) or per (chromosome, haplotype) (`phased`), with $m$ the per-base multiplicity
    (§1); 0 for an empty denominator;
  - `n_items_truth` $= n_T$, `n_items_inferred` $= n_C$ and `n_items_shared` $= n_\cap$, so that either
    form can be recomputed.
- Other fields: `case` (`aa_eval --case`, default the cycles file name), `truth_weight` $= w_T$,
  `inferred`, `copy_number` $= w_C$, `low_weight`, and the flags `merge_column`, `split_row` below.

**Items** (`column_items`), per haplotype block at $\tau_{\mathrm{sp}}$, native universe:
1. **Strict columns** use M2's universe and item sets (`aa_eval.adjacency_items`): truth items, then
   the strict adjacencies matched one to one, nearest first, within $\tau_{\mathrm{sp}}$. The species are
   passed in matrix order (by copy number), not in input order, so a tie in distance can be broken
   differently from the `aa_partition_metrics.tsv` rows.
2. **The adjacencies of `low_weight` columns are added afterwards**, equal ones merged. Each takes, in
   order:
   1. an existing item it equals (any item: truth, or a strict false item). This does not mark a
      truth item as taken;
   2. otherwise a truth item that no strict adjacency and no earlier adjacency of this step took, one
      to one and nearest first within $\tau_{\mathrm{sp}}$ (sorted as in M2);
   3. otherwise a new false item.

   A truth item can therefore hold one `low_weight` adjacency from step 1 and another from step 2.

**Worked example** (`tests/ecdna/test_species_matrix.py`, blind). Two truth species:
- the outer O: 2 Mb, novel adjacencies $a_1, a_2$;
- the nested I: 0.5 Mb inside O in base space, novel adjacencies $a_3, a_4$.

The merged walk M is their Euler union: 2.5 Mb, I's bases covered twice, $a_1$ to $a_4$.

| answer | cells on row O | cells on row I |
|---|---|---|
| correct, columns O, I | (O, O) = (1, 1); (O, I) = (0.25, 0) | (I, O) = (0.25, 0); (I, I) = (1, 1) |
| merged, one column M | (O, M) = (0.8, 0.5) | (I, M) = (0.2, 0.5) |
| truth is M; the answer reports O and I (split) | (M, O) = (0.8, 0.5); (M, I) = (0.2, 0.5) | |

**Reading 0 and 1.**

| cell | meaning |
|---|---|
| (1, 1) | the same molecule, up to walk order (order: `exact_cycle`) |
| (high, low) | C lies on T's DNA but joins it differently: the nested-merge signature |
| (> 0, 0) | a neighbour: shared reference DNA, no shared novel adjacency (a nested species against its outer one) |
| (0, 0) | unrelated |

Two flags, computed per haplotype block on the strict columns only:
- a strict column with $\ge 2$ rows of $J_\mathrm{adj} > 0$ is a **merge** (`merge_column`, on each of
  its cells; always false on `low_weight` columns);
- a row with $\ge 2$ strict columns of $J_\mathrm{adj} > 0$ is a **split** (`split_row`, on each of its
  cells).

NaN counts as 0 in both.

**Bench, λ = 0.2** (`diag/species_matrix_bench.py`, 40 cases, blind):
- The 55 cells with `exact_cycle` = 1 are all $(1, 1)$, up to $J_\mathrm{seg} \ge 0.99996$ (next
  paragraph).
- `nest_samehap_3to1_s301`:
  - the outer species is exact;
  - the two light strict walks score (0.75, 0.38) and (0.32, 0.25) on the outer row, and
    (0.07, 0.50) and (0.09, 0.60) on the nested row;
  - so both are merge columns: walks re-traversing the outer that carry the nested species'
    adjacencies.
- `nest_samehap_6to1_s801`:
  - the outer species is exact;
  - the nested species is spread over the three `low_weight` columns: (0.97, 0.80), (0.55, 0.67)
    and (0.51, 0.44).
- `ctrl_dup_s708`: one truth species, two strict columns, (0.20, 0.40) and (0.80, 0.80): a split.
- The passing `nest_cross` cases (`nest_cross_2to1_s103`, `nest_cross_1to1_s506`): on the strict
  columns, a clean diagonal of (1, 1), and (0.47–0.50, 0) off it.

**Decisions and reasoning.**

- **Items, not pairs, in J_adj.**
  - The pair form is $0/0$ for every cell where both species have at most one novel adjacency: a
    simple circle with only its closing adjacency, or s306's nested species.
  - On one cell, the pair form is a function of the three item counts only, so it adds no
    information.
  - It is still written (`j_adj_pairs`) for comparison with M2.
- **M2's universe, unchanged for strict columns.** Every strict cell then uses the same item sets
  as M2 (up to distance ties, **Items** 1). Adding `low_weight` adjacencies into the nearest-first matching
  could hand a truth item to a leftover and change a strict cell, so they come second and take only
  what the strict answer left.
- **`low_weight` columns shown, but kept out of the flags.**
  - They are not part of the answer (M1 and M2 leave them out).
  - But the 6:1 failure *is* the nested species living in them, so they are shown.
  - Leftovers often reuse the outer species' adjacencies (s801's leftovers score 0.11 and 0.18 on
    the outer row), which would flag a false split.
- **No copy number in the tuple.** The tuple asks *which molecule*; copy number is §3.5's and the
  composite's job. `j_seg_cn` is the copy-weighted variant, for when it matters.
- **Rows and columns by copy number**, so that the diagonal is the right answer, and the
  highest-copy, most confident species comes first.

**Caveats and undefined cases.**
- **$J_\mathrm{adj}$ is undefined (NaN)** when neither species has a novel adjacency. That means two
  linear walks without one; a circle always has its closing adjacency.
- **`exact_cycle` = 1 does not give exactly (1, 1).**
  - `exact_cycle` allows $\tau_{\mathrm{gr}}$ = 100 bp at each run end, so $J_\mathrm{seg}$ can fall just short
    of 1. On the bench every exact cell has $0.99996 \le J_\mathrm{seg} < 1$, the AA export's 1 bp
    coordinate convention (§5.1), the same as `segment_jaccard_h0` in `aa_species_metrics.tsv`.
  - $J_\mathrm{adj}$ can fall below 1 if another strict species has a different adjacency nearer to
    the same truth item, because the one-to-one matching then makes the exact walk's adjacency a
    false item. This never happens on the bench.
- **Multiplicity is invisible to $J_\mathrm{adj}$.** Items are distinct, so a novel adjacency traversed
  twice counts once. Duplications (bench `ctrl_dup`) are separated by $J_\mathrm{seg}$, which counts
  multiplicity, and by `exact_cycle`.
- **The flags are diagnostics, not classes.** A merged walk that re-traverses the outer species'
  adjacencies (s301) makes its column a merge *and* the outer row a split. A truth adjacency shared
  by two species (none on the bench, §4.7 M2) would flag a correct answer.
- **Fold-backs and tolerance.** Items use M2's extremity matching (same chromosome and strand,
  either direction), so a fold-back is an item like any other. Two distinct strict adjacencies
  within $\tau_{\mathrm{sp}}$ of one truth item still match one to one.
- **Summing cells does not reproduce M2.**
  - Inferred species can share items (a merged walk and its outer species), and M2 counts a pair
    once if *some* species holds it.
  - The matrix also holds `low_weight` columns that M2 leaves out.
  - Per cell, $J_\mathrm{adj}$ is a ratio of item counts, not of pairs.

#### How M1 and M2 handle methods, modes, linear walks, `low_weight` and tolerance

| | M1 species detection | M2 adjacency pairs |
|---|---|---|
| **Haplotype-blind methods (AA, CoRAL, CycleExtractor, original RCK)** | `blind` rows only: blind matching $\mathrm{match}^0$, `exact_cycle`, `composite_h0` | `blind` rows only: items are novel adjacencies of the haplotype-stripped walks |
| **Phased answers (`rck-ecdna`)** | `blind` and `phased` rows; compare methods on `blind` | the same |
| **`--phase-blocks`** | the strict and candidate walks are relabelled per phase block first (§3.10), so the `phased` rows are within-block | the same; `blind` rows are unchanged |
| **Strict mode** | TP, FP, FN, precision, recall, $F_1$ | TP, FP, FN, TN, precision, recall, FPR, $F_1$ |
| **Best-of mode** | TP, FN, recall only | not computed |
| **Linear walks** | count in $\lvert\mathcal S\rvert$; never `exact` (truth species are circular, and `order_match` needs the same topology); may pass a composite threshold | no wrap-around adjacency, so a linear walk has one fewer item. CoRAL's and CycleExtractor's `Path=` lines are read only once rewritten (§4.1, §7) |
| **`low_weight` species** | not in $\mathcal S$ (`status` rule; `cn-frac:0.1` for AA and CoRAL), so neither a prediction nor an FP. In strict mode a molecule split into sub-threshold walks (bench `nest_samehap_6to1`) is an FN; in best-of they are candidates | not in $\mathcal S$: their adjacencies are not predicted items. `--strict-rule all` scores every walk |
| **10 kb against 100 bp** | `exact` at $\tau_{\mathrm{gr}}$; the composite thresholds at $\tau_{\mathrm{sp}}$, as the composite itself | both, in rows with `tol_bp` 10000 (headline) and 100 |

### 4.8 Regression against `rck-ecdna-eval`

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. It is a software regression test, not a metric: it
> checks that `aa_eval` reproduces §3 and introduces no new definitions. **Code:**
> `results/bench_nested/diag/aa_regression.py`; `results/bench_nested/collect.py --aa-metrics`.

The AA-format path must reproduce §3 exactly for `rck-ecdna`'s answers. `aa_regression.py` exports
each answer with `aa_export` (rank 1 and its first 2 same-fit alternatives), scores it with
`aa_eval`, and compares the strict species rows field by field with `rck-ecdna-eval`'s
(`lambda_0/recovery_matrix.tsv` at λ = 0, `pool/pool/evaluation.tsv` at the re-ranked pool entry for
λ > 0). On the 40 bench_nested cases at λ = 0 and 0.2 (`--lambdas` default; 80 pairs), the following
are identical at full precision:

- every part of the composite (`inferred`, `segment_jaccard`, `haplotype_jaccard`, `junction_f1`,
  `phased_junction_f1`, `phasing_accuracy`, `order_match`, `cn_rel_err`);
- every composite;
- the matched species, `n_species`, `n_truth` and `merged`.

With `--aa-metrics` all 240 (case, λ) match, and the existing columns of `bench_sweep.tsv`,
`pipeline.tsv` and `rerank_validation.tsv` are unchanged. Two things were needed:

- a species' other orderings: `nest_samehap_3to1_s304` matches its order only through ordering
  2, which RCK's own `aa/` files leave out;
- reading the old tables with `float_precision="round_trip"`.

---

## 5. Translators: `aa_export`, `aa_import`, `coral_flank`

All AA files below use AA's 1-based inclusive coordinates. The walks and subgraphs of `rck-ecdna`
use 0-based inclusive coordinates (§1). The translators convert with `SHIFT` = 1
(`rck/ecdna/aa_export.py`): AA position = walk position + 1, and `aa_eval` reads an AA position
back as position − 1 (§4.1).

### 5.1 `rck-ecdna` answers in AA format (`aa_export RUN OUT`)

> **Provenance:** (c) new, RCK yubi_dev@e4c2409. **Code:** `rck/ecdna/aa_export.py::export`,
> `::write_cycles`, `::run_alternatives`, `::pool_alternatives`, `::shift_graph`, `::write_segment_cn`.

The files `rck-ecdna` writes itself (`aa/<name>_decomp<k>_cycles.txt`) lack the `low_weight`
status, the other orderings, and the solutions of a kept pool (`--keep-pool`), which the λ > 0
results come from.
`python -m rck.ecdna.aa_export RUN OUT [--pool-index I | --lambda L] [--alternatives 2]` writes
these from `RUN/decompositions.tsv`, or from `RUN/pool/decompositions.tsv` with a pool option:

- no option: the answer is the lowest `decomposition_rank` of `decompositions.tsv`;
- `--pool-index I`: pool entry $I$ (0-based; its `decomposition_rank` is $I+1$);
- `--lambda L`: the pool entry `rck.ecdna.rerank` ranks first at `--step-penalty` L.

`--pool-index` and `--lambda` exclude each other. It writes, with `<name>` = `summary.json`
`graph.name`:

- `<name>_cycles.txt`: the answer. `Interval` lines are the run's regions (`summary.json`
  `graph.regions`). `Segment` lines number the distinct segments of all its walks, sorted by
  (chromosome name, start, end). One `Cycle=` line per species, in file order:
  `Copy_count` = the refitted `weight`, `Segments=`, `Haplotypes=` (one $H_s$ per step),
  `Status=` (`species` or `low_weight`), `Species=`, and the other orderings as
  `Alt_segments=` / `Alt_haplotypes=` (`|`-separated). A linear walk is written AA's way,
  `0+,…,0-`, with haplotype `.` on the two source tokens;
- `<name>_alt<k>_cycles.txt`, $k = 1, 2, \dots$: the same for rank 1's first `--alternatives`
  (default 2) same-fit alternatives:
  - without a pool option: the `rank` of the first `--alternatives` entries of `summary.json`
    `ambiguity.alternatives`; an entry without a rank is skipped;
  - with one: the first `--alternatives` entries of `pool.json` `same_fit[I]`, each matched to the
    first pool entry with $\lvert\Delta\,\mathrm{objective}\rvert < 10^{-3}$, equal `steps` and
    equal sorted significant (not `low_weight`) integer weights, as
    `results/bench_nested/diag/grow_rules.py` matches them (§6.2); an unmatched one is skipped;
- `<name>_graph.txt`: the first `RUN/aa/*_graph.txt`, every position + 1 (only if one exists);
- `<name>_segment_cn.tsv`: columns `chrom start end cn cn_a cn_b` from `flow_segments.tsv`, with
  `cn_a` = `ecdna_A` + `bg_A`, `cn_b` = `ecdna_B` + `bg_B`, `cn` = `cn_a` + `cn_b`, positions + 1
  (only if `flow_segments.tsv` exists).

Every position is the walk's coordinate + 1, so the scorer reads back exactly the steps
`rck-ecdna-eval` reads from `decompositions.tsv`, and scores them identically
(`tests/ecdna/test_aa_eval.py::TestExport`). Graph edges match the truth within 1 bp either way,
because a subgraph segment's start and end are not in one convention.

### 5.2 Subgraph → AA graph (`aa_export subgraph`)

> **Provenance:** (c) new, RCK yubi_dev@2f921a5. The graph layout and sort order are (a) CoRAL 3.1.1's
> and CycleExtractor 1.0.0's parsers [CoRAL 3.1.1 code](references.md#coral_code) [Faizrahnemoon et al. 2026](references.md#cycleextractor); the path-constraint
> tokens are CoRAL's `coral/core_utils.py::path_to_str`. **Code:**
> `rck/ecdna/aa_export.py::export_subgraph`, `::subgraph_observations`, `::_path_edges`, `::_runs`;
> `rck/ecdna/aa_format.py::write_graph`, `::coral_path_string`, `::chrom_key`, `::edge_key`.

The graph route gives every decomposer the same breakpoint graph, so that decompositions are
compared on equal input (`results/graph_route/`).
`python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name N] [--amplicon 1] [--imputed comment|drop] [--min-path-junctions 1]`
writes `<name>_amplicon<N>_graph.txt`, `<name>_amplicon<N>_segment_cn.tsv` and
`<name>_amplicon<N>_export.json` (`<name>` defaults to the folder name; CoRAL reads $N$ from the
file name). The graph is the normalised graph `rck-ecdna` solves on (`subgraph_io.read_subgraph`:
phantom segments of at most `--phantom-max-bp` (2) bp contracted, parallel imputed edges merged),
with the copy numbers it uses (`scale.build_observations`):

```math
\mathrm{CN}(j) = \mathrm{CN\_TOTAL}(j),
\qquad
\mathrm{CN}(a) = \frac{\mathrm{reads}(a)}{\theta}
```

where

- $j$ is a segment, and `CN_TOTAL` is its `edges.tsv` column. It is used when the folder has
  `scale.json` and no `--depth-per-copy` is given; otherwise $\mathrm{CN}(j)$ = depth / $\gamma$;
- $a$ is a novel or reference adjacency with an observation, and $\mathrm{reads}(a)$ its read
  count (`WEIGHT`);
- $\theta$ is the junction reads per copy (`docs/ecDNA.md`) of `scale.json`, or
  `--junction-reads-per-copy`; $\gamma$ is the depth per copy, likewise.

This is the *observed* adjacency copy number, not the rounded `CN` column of `edges.tsv`.

| Line | Content | Seen by CoRAL / CycleExtractor |
|---|---|---|
| `sequence` | one per segment, sorted by (`chrom_key`, start) = CoRAL's order; CN = $\mathrm{CN}(j)$, AverageCoverage = the segment's depth (`WEIGHT`), NumberReadsMapped = 0 | yes |
| `concordant` | reference adjacencies between two neighbouring segment ends; CN = $\mathrm{CN}(a)$, or for an unobserved join (an imputed grid point) the smaller of the two segment CNs (listed in `_export.json`) | yes |
| `discordant` | novel adjacencies between two segments; CN = $\mathrm{CN}(a)$, 0 when unobserved; fold-backs (both ends on one node) as self-loops, listed in `_export.json` | yes |
| `path_constraint` | read paths of `paths.tsv` (below), as CoRAL 3.1.1 tokens; support = `#NOV` rounded | yes |
| `interval` | the runs of contiguous segments (next start = previous end + 1), what CoRAL's `infer_amplicon_intervals` infers | ignored |
| `#rck_imputed` | RCK-ecDNA's imputed candidates (no read evidence), with their support; `--imputed drop` leaves them out | no |
| `#rck_outside` | adjacencies with an end on a node without a segment (a region's outer flank); CoRAL and CycleExtractor add their own s/t ends there | no |
| `#rck_node`, `#rck_seq`, `#rck_edge`, `#rck_pc` | the subgraph's node, segment, adjacency and path ids and orders, so the import is exact | no |
| `#rck_path` | read paths `rck-ecdna` reads that cannot be written as a constraint, with the reason | no |

Breakpoint edges are written lower breakend first and sorted by CoRAL's `BreakpointEdge` order
(`edge_key`: chromosome, position of the first end, then of the second); concordant before
discordant. The export stops with an error when a region is not a run of contiguous segments
(CoRAL would infer other intervals), or when a reference adjacency does not join neighbouring
segment ends.

A read path of `paths.tsv` becomes a `path_constraint` by these steps:

1. Resolve its node steps. A read may start or end outside the graph; its outer node is then the
   other end of the segment it enters or leaves by. A read with another unresolved node is dropped,
   as `rck-ecdna` drops it.
2. Walk it in CoRAL's layout: sequence edge, node strand, breakpoint edge, node strand, …,
   sequence edge (`_path_edges`). Steps $(2i, 2i+1)$ run along segments and reference joins;
   steps $(2i+1, 2i+2)$ are one adjacency.
3. Write it when it uses only segments, reference joins and observed novel adjacencies and has at
   least $\max(1,$ `--min-path-junctions`$)$ novel adjacencies. CoRAL reverses a path whose first
   edge index is larger than its last; the reversal is recorded in `#rck_pc`.
4. Otherwise: a path with no novel adjacency is dropped (`rck-ecdna` ignores it too); one with too
   few, through an imputed candidate, or through a join the graph lacks becomes a `#rck_path`
   line (through an imputed candidate with `--imputed drop`: dropped).

Path constraints index sequence and concordant edges in CoRAL's sort order because `coral cycle`
re-sorts them (`infer_amplicon_intervals` → `sort_edges`) without renumbering constraints that were
already parsed by file index.

### 5.3 AA graph → subgraph (`aa_import`)

> **Provenance:** (c) new, RCK yubi_dev@2f921a5. The `regenerate` mode reuses (c)
> `subgraph_construct.py`'s imputation rule and amplification threshold (amplicon-allele-analysis).
> **Code:** `rck/ecdna/aa_import.py::import_graph`, `::fit_scale`, `::regenerate_imputed`,
> `::amplified_threshold`, `::weighted_median`, `::compare_subgraphs`.

`python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate] [--cn-profile F] [--background-cn auto|X]`
turns any AA graph (AA, CoRAL, the truth mixture) into a subgraph folder for
`rck-ecdna --inputs subgraph`: `nodes.tsv`, `edges.tsv`, `paths.tsv`, `scale.json`,
`cn_profile.tsv`, `summary.txt` and `import.json` (what was mapped, dropped and fitted).
Positions are AA − 1. The mapping:

- `sequence` becomes a `reference` edge (`WEIGHT` = AverageCoverage, `CN_TOTAL` = CN);
- `concordant` becomes an `extremity` edge (`WEIGHT` = reads, `CN` = CN);
- `discordant` becomes a `variant` edge (`WEIGHT` = `NOV_U` = reads, `NOV_A` = `NOV_B` = 0) with no
  haplotype label, so every phasing is allowed;
- `source` edges are dropped and listed in `import.json`;
- `path_constraint` lines become `paths.tsv` rows; `#rck_path` lines are copied back.

A graph written by `aa_export subgraph` also restores the subgraph's ids and orders from its
`#rck_*` lines. Any other graph gets node ids `chr:pos±_±` as `subgraph_construct.py` names them.

The scale is fitted from the graph itself (`fit_scale`):

```math
\hat\theta = \operatorname{median}_{e \in \mathcal E_\theta} \frac{\mathrm{rd}_e}{\mathrm{cn}_e},
\qquad
\hat\gamma = \operatorname{wmedian}_{j:\ \mathrm{cov}_j > 0,\ \mathrm{CN}(j) > 0}\Big( \frac{\mathrm{cov}_j}{\mathrm{CN}(j)};\ L(j) \Big)
```

where

- $\mathrm{rd}_e$ and $\mathrm{cn}_e$ are a breakpoint edge's read count and copy number;
- $\mathcal E_\theta$ is the discordant edges with $\mathrm{rd}_e > 0$ and $\mathrm{cn}_e > 0$, or the
  concordant edges with both if no discordant edge qualifies;
- $\mathrm{cov}_j$ is a sequence edge's AverageCoverage, $\mathrm{CN}(j)$ its CN and
  $L(j) = \mathrm{end} - \mathrm{start} + 1$ its length;
- $\operatorname{wmedian}(v; w)$ is the lower weighted median: the smallest $v$ whose cumulative
  weight, in increasing $v$, reaches half the total;
- a graph without any read count gets $\hat\theta = 1$, and every breakpoint edge gets `WEIGHT` = CN;
- a CoRAL-rewritten graph (a `#CoRAL` header) gets $\hat\gamma = 1$ and reference `WEIGHT` = CN,
  because CoRAL's writer puts `lr_nc / size` in AverageCoverage; `--depth-per-copy` overrides
  $\hat\gamma$; a graph without coverage gets $\hat\gamma = 1$;
- a ratio that agrees to $10^{-9}$ (relative) with its 12-significant-digit rounding is snapped to
  it, so 6.5 comes back as 6.5 (`_snap`).

The background copy number $c_{\mathrm{bg}}$ (`scale.json` `background_cn`) is `--background-cn`;
with `auto` (default) it is the length-weighted median of `--cn-profile`'s `CN_TOTAL`, else 2.
`cn_profile.tsv` is a copy of `--cn-profile`, or else the sequence edges at their CN plus 1 Mb flanks
on both sides of every run at $c_{\mathrm{bg}}$, kept off other runs.

Imputed candidates, `--imputed`:

- `comment` (default) restores the `#rck_imputed` lines and the `#rck_path` lines through them;
- `none` leaves both out, which is what CoRAL and CycleExtractor see;
- `regenerate` proposes closing candidates by `subgraph_construct.py`'s rule:
  1. per chromosome, the threshold is $\min(\tilde c + 1.5\,\mathrm{MAD},\ 1.5\,\tilde c)$ over the
     chromosome's `CN_TOTAL` bins of `--cn-profile`, with $\tilde c$ their median and MAD
     $= 1.4826 \times$ the median absolute deviation; without a profile it is $1.5\,c_{\mathrm{bg}}$;
  2. a sequence edge is amplified when its CN is strictly above the threshold;
  3. an amplified run is a maximal run of contiguous amplified sequence edges without a step:
     a step is a CN change of at least $\max(2,\ 0.25 \times$ the larger CN$)$
     (`--boundary-min-copies`, `--boundary-frac`);
  4. each run gets an `intra` edge from its last head back to its first tail, with support = its
     length-weighted mean coverage (CN where there is none);
  5. each pair of consecutive runs on one chromosome gets an `inter` edge from the first's head
     to the second's tail, with support = the mean of the two supports.

  `subgraph_construct.py`'s third kind, `anchor`, is not made: every segment end here is already a
  breakpoint.

`#rck_outside` joins are always restored.

**Round trip.** `aa_import.compare_subgraphs` (`--compare SUBGRAPH_DIR`) compares the normalised
graphs `rck-ecdna` builds from two folders: nodes (ids, coordinates, terminals, order), segments,
adjacencies (ids, ends, kind, reads, allowed phasings, order), observations (segment CN, junction
CN, flank background), read paths and regions. Junction labels are ignored
(`soft_junction_hap`). Numbers agree when $\lvert x - y \rvert \le 10^{-12} \max(1, \lvert x \rvert, \lvert y \rvert)$; $\hat\theta$ must be equal. With `--imputed none` or `regenerate`, the
original's imputed candidates and the paths through them are left out of the comparison. Export
then import with `--imputed comment` and the original `cn_profile.tsv` is exact on 36 of 36 cases
(demo5–demo8 and the 32 bench cases). `--imputed none` is exact apart from the imputed candidates.
A graph-only `none` differs only in the flank background (RCK `results/graph_route/README.md`,
"Validation done offline"). Total-only runs on both sides use
`--min-snps-allele 100000000 --soft-junction-hap --soft-junction-penalty 0`, which ignores junction
labels exactly.

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
x_{(\mathrm s, v)} + x_{(v, \mathrm t)} = x_{\mathrm{seq}(v)}
```

where $x$ is a cycle's edge multiplicity (`docs/ecDNA.md`'s $x_e$), $\mathrm s$ and
$\mathrm t$ are CoRAL's source and sink, and $\mathrm{seq}(v)$ is the sequence edge at $v$.
Discordant and concordant edges at an end node are left out of the balance, so a cycle through a
breakpoint at an interval end is infeasible. CoRAL's own graphs extend each amplicon beyond its
breakpoints. Our graphs, the truth mixtures and `aa_export subgraph`, end at breakpoints in 69 of
75 cases. On them CoRAL returned zero-weight paths and explained 0 copies (RCK
`results/graph_route/README.md`, `NEXT_STEPS.md`; RCK `CHANGELOG.md`).

**The fix.** `python -m rck.ecdna.coral_flank GRAPH OUT_GRAPH [--flank 10000]`. For every run
$[\mathrm{start}_r, \mathrm{end}_r]$ of contiguous sequence edges (sorted in CoRAL's order), a
flank sequence edge is added at each end, joined to the run by a concordant edge. Its length is

```math
L_{\mathrm{left}} = \min\Big(L_{\mathrm{fl}},\ \Big\lfloor \frac{\mathrm{start}_r - \mathrm{end}_{r^-} - 2}{2} \Big\rfloor\Big),
\qquad
L_{\mathrm{right}} = \min\Big(L_{\mathrm{fl}},\ \Big\lfloor \frac{\mathrm{start}_{r^+} - \mathrm{end}_r - 2}{2} \Big\rfloor\Big)
```

where

- $L_{\mathrm{fl}}$ = `--flank` (default 10 000 bp, at least 1);
- $\mathrm{end}_{r^-}$ is the largest end of an earlier run on the chromosome; without one,
  $L_{\mathrm{left}} = \min(L_{\mathrm{fl}},\ \mathrm{start}_r - 1)$;
- $\mathrm{start}_{r^+}$ is the smallest start of a later run on the chromosome; without one,
  $L_{\mathrm{right}} = L_{\mathrm{fl}}$ (no chromosome-end check).

So two flanks never touch, which would make CoRAL merge them. A length below 1 is an error. The
flank's copy number, also that of its concordant edge, is

```math
\mathrm{cn}(\text{flank at } v) = \max\Big(0,\ \mathrm{cn}(\mathrm{seq}(v)) - \sum_{e \in E_{\mathrm d}(v)} \mathrm{cn}(e)\Big)
```

where $v$ is the run's outer node and $E_{\mathrm d}(v)$ is the discordant edges with an end at
$v$, a self-loop counted twice (as CoRAL lists it twice). Negative residuals are listed as
`negative_residuals` in the printed info. Coverage and read counts of the new edges are 0.
Sequence and concordant edges are written in CoRAL's sort order, and path-constraint `e`/`c`
indices are renumbered. Discordant edges and their order are unchanged; `source` edges are kept.
`#` and `interval` lines are not copied.

**Effect.** With CoRAL's own parser and model offline (balance, copy number and multiplicity, not
the spanning-tree constraints), truth cycles feasible in CoRAL's greedy model go from 28 to 66 of
67. The remaining one (`nest_samehap_6to1_s804`, CN 40.8) exceeds CoRAL's own read-clustered
discordant-multiplicity bound. Breakpoints on end nodes go from 157 to 0. Slurm test 32104711
recovered the demo5 truth cycle at CN 16 (RCK `results/graph_route/NEXT_STEPS.md`). The graph route
gives the flanked graph to CoRAL only: CycleExtractor joins s/t at every node, and RCK-ecDNA has
its own region-flank joins and flank background, so the fix gives CoRAL the same picture, not an
advantage.

### 5.5 Bulk graph with the subgraph's read paths (`aa_export baseline`)

> **Provenance:** (c) new, RCK `yubi_CNP_ionic` (2026-11-14). **Code:**
> `rck/ecdna/aa_export.py::baseline_graph`, `::main_baseline`. Used by `results/bench_IW/bench_iw.py`
> (`--set iw_intpsi`, §6.4).

`python -m rck.ecdna.aa_export baseline BULK_RUN_DIR SUBGRAPH_GRAPH DST` writes the graph that
CoRAL and CycleExtractor decompose as *decomposition-only baselines*: `rck-ecdna`'s bulk graph with
the long-read path constraints of the direct export (§5.2). The baselines then decompose the graph
the bulk route fitted, so the comparison isolates the decomposition. Steps:

1. Read `BULK_RUN_DIR/aa/*_graph.txt` (exactly one) with every position + 1. Its copy numbers are
   the bulk route's fitted ones (background plus ecDNA), and every edge it used is present,
   imputed ones included.
2. Sort the sequence edges in CoRAL's order. Keep a concordant or discordant edge when both ends
   are sequence-edge ends; write it lower end first, sorted by `edge_key`. An edge to any other
   node is left out and recorded (`#rck_outside_bulk` lines, `outside_edges`): it continues the
   reference past an interval end, and `coral_flank` adds the same join. Any other edge kind
   (`source`) is an error.
3. For each `path_constraint` of `SUBGRAPH_GRAPH`, renumber every `e`/`c`/`d` index to the edge of
   step 2 with the same coordinates. Token strands and counts, and the support, are kept. A
   constraint with an edge not found is dropped, counted by reason (`no_sequence_edge`,
   `no_concordant_edge`, `no_discordant_edge`). Constraints through an edge with CN ≤ 0 are kept
   and counted.
4. Write `interval` lines for the runs of contiguous sequence edges, and `<DST stem>.json` with the
   counts (`sequence`, `concordant`, `discordant`, `outside_dropped`, `path_constraints_in`,
   `_written`, `_dropped`, `_through_zero_cn`, `interval`).

---

## 6. Benchmark-level summaries built on these metrics

The λ choice rule and the failure-stage rule of `results/bench_nested/` are defined in its
`collect.py` docstring and README. Two summaries defined on the composite are used across
documents (§6.1, §6.2), §6.3 pools the confusion-matrix metrics of §4.7, and §6.4 lists what the
bench_IW benchmark scores.

Notation: $\mathcal K$ is the set of cases of a benchmark, and $\mathfrak c \in \mathcal K$ one
case.

### 6.1 Mean composite of a benchmark

> **Provenance:** (c) new, RCK `results/bench_nested/`. A macro average over cases
> [Sokolova & Lapalme 2009](references.md#sokolova2009). **Code:** `results/bench_nested/collect.py::run_row`,
> `::pool_rows`, `::choose`.

```math
\overline{\mathrm{comp}}(\lambda) = \frac{1}{\lvert \mathcal K \rvert} \sum_{\mathfrak c \in \mathcal K} \overline{\mathrm{comp}}_{\mathfrak c}(\lambda)
```

where

- $\overline{\mathrm{comp}}_{\mathfrak c}(\lambda)$ is case $\mathfrak c$'s `mean_composite` (§3.8,
  the mean `composite` over its in-scope truth species) for the rank-1 decomposition at step
  penalty λ: from the λ = 0 run's `recovery_matrix.tsv`, or from the kept pool re-ranked at
  λ > 0 (`rck.ecdna.rerank`);
- $\mathcal K$ is the cases that have a result at every λ of the grid.

The reported λ is the chosen λ (`collect.py` docstring), or `--lambda`.

### 6.2 Set metric (rank 1 plus its first same-fit alternatives)

> **Provenance:** (c) new, RCK yubi_dev@de5a805. **Code:** `results/bench_nested/diag/grow_rules.py::main`
> (`set k`). It is an oracle over rank 1 and its first $n_{\mathrm{alt}}$ same-fit alternatives from
> the ambiguity report (`docs/ecDNA.md` §5), macro-averaged over cases as in §6.1.

```math
\mathrm{set}_{n_{\mathrm{alt}}} = \frac{1}{\lvert \mathcal K \rvert} \sum_{\mathfrak c \in \mathcal K} \max\Big( \overline{\mathrm{comp}}_{\mathfrak c}(\mathcal W_{\mathfrak c,1}),\ \max_{l \in A_{\mathfrak c}(n_{\mathrm{alt}})} \overline{\mathrm{comp}}_{\mathfrak c}(\mathcal W^{\mathrm{alt}}_{\mathfrak c, l}) \Big)
```

where

- $\mathcal K$ is the cases with a kept pool (`<case>/pool/pool/pool.json`);
- $\mathcal W_{\mathfrak c,1}$ is the walks of case $\mathfrak c$'s rank-1 decomposition at λ
  (`--lambda`, default 0.2), ranked by `decompose.rank_pool`;
- $\overline{\mathrm{comp}}_{\mathfrak c}(\cdot)$ is the mean `composite` of that decomposition's
  rows of `pool/pool/evaluation.tsv` (in-scope truth species; `low_weight` species not scored);
- $\mathcal W^{\mathrm{alt}}_{\mathfrak c, l}$ is rank 1's $l$-th same-fit alternative
  (`pool.json` `same_fit`), matched to the first pool entry with
  $\lvert\Delta\,\mathrm{objective}\rvert < 10^{-3}$, equal `steps` and equal sorted significant
  integer weights;
- $A_{\mathfrak c}(n_{\mathrm{alt}})$ is the matched ones among the first $n_{\mathrm{alt}}$
  alternatives (unmatched ones are skipped; the inner max is omitted when none matches).

The script reports $n_{\mathrm{alt}} \in \lbrace 1, 2, 5 \rbrace$ (columns `set 1`, `set 2`,
`set 5`) and `pool_best`, the best decomposition of the pool. Unlike §3.8's `best_composite`, which
takes each truth species' best over all ranks, the set metric takes the best *decomposition*.

**Why.** The grown bench's failure analysis found the right answer among rank 1's same-fit
alternatives in several failing cases (s703, s301, and others), and every re-ranking rule tried
broke as many cases as it fixed. **Effect** at λ = 0.2 over 40 cases
(`results/bench_nested/NEXT_STEPS.md`, "Grown bench: failure analysis"): rank-1 mean 0.8531;
$\mathrm{set}_1 = 0.9113$, $\mathrm{set}_2 = 0.9229$, $\mathrm{set}_5 = 0.9373$; pool best 0.9470.
The set metric changes no case downward. It is not a ranking: it scores what an ambiguity report
would let a user see.

### 6.3 Pooled confusion-matrix metrics (M1, M2)

> **Provenance:** (c) new, RCK yubi_dev, 2026-11-07. Micro and macro averaging as defined by
> [Sokolova & Lapalme 2009](references.md#sokolova2009). **Code:** `rck/ecdna/aa_eval.py::pool_partition`,
> `::confusion_rates`; `results/bench_nested/collect.py::aa_partition`;
> `results/graph_route/score_h2h.py::partition_rows`.

For a group $\mathcal G$ of scored answers $\mathfrak a$ (cases, or (case, amplicon) runs in the
graph route) that share `metric`, `mode`, `haplotype`, `criterion`, `tol_bp` and `universe` (and
the caller's keys, below), with per-answer counts $TP_{\mathfrak a}, FP_{\mathfrak a}, FN_{\mathfrak a}, TN_{\mathfrak a}$ (§4.7):

```math
\mathrm{micro\_rate} = \mathrm{rate}\Big(\sum_{\mathfrak a \in \mathcal G} TP_{\mathfrak a},\ \sum_{\mathfrak a} FP_{\mathfrak a},\ \sum_{\mathfrak a} FN_{\mathfrak a},\ \sum_{\mathfrak a} TN_{\mathfrak a}\Big),
\qquad
\mathrm{macro\_rate} = \frac{1}{\lvert \mathcal G_{\mathrm{rate}} \rvert} \sum_{\mathfrak a \in \mathcal G_{\mathrm{rate}}} \mathrm{rate}\big(TP_{\mathfrak a}, FP_{\mathfrak a}, FN_{\mathfrak a}, TN_{\mathfrak a}\big)
```

where

- $\mathrm{rate}$ is precision, recall, FPR or $F_1$ (§4.7, `confusion_rates`);
- $\mathcal G_{\mathrm{rate}}$ is the answers whose $\mathrm{rate}$ is defined, and `n_cases_<r>` (`<r>` the rate's name) = $\lvert \mathcal G_{\mathrm{rate}} \rvert$;
  `macro_<r>` is undefined (NaN) when $\mathcal G_{\mathrm{rate}}$ is empty;
- $\sum FP$ is undefined when any answer's FP is (best-of rows), and then micro precision, $F_1$
  and FPR are undefined; $\sum TN$ likewise (M1 has no TN), and then micro FPR is undefined.
  Micro recall needs only TP and FN;
- `headline` of a pooled row is true when every pooled row is a headline row.

`headline_pooling` names the one the docs report: `micro` for M1, so that every molecule weighs the
same, and `macro` for M2, so that one large case does not dominate (pair counts grow
quadratically with a case's size; §4.7).

Files, in addition to the existing tables, which are unchanged:

- `python results/bench_nested/collect.py --aa-metrics` (opt-in): `aa_partition_metrics.tsv` (every
  (case, λ)'s rows with `case`, `kind`, `step_penalty` in front) and `aa_partition_pooled.tsv`
  (pooled per `step_penalty`, over all cases and over the bench cases, `cases` = `all` or
  `bench`).
- `python results/graph_route/score_h2h.py` (always): `h2h_partition.tsv` (every run's rows with
  `case`, `scenario`, `graph`, `method`, `amplicon` in front) and `h2h_partition_pooled.tsv` (per
  graph and method with `scenario` = `all`, and per scenario, graph and method). The
  `--drop-paths` sensitivity run writes the same files with suffix `_cycles_only`. The unit of
  macro pooling is the (case, amplicon) run.

### 6.4 bench_IW: copy-number recall, ecDNA proportions and the per-run table

> **Provenance:** (c) new, RCK `yubi_CNP_ionic` (2026-11-14). **Code:** `rck/ecdna/bench_metrics.py::score`,
> `::cn_scores`, `::proportion_scores`, `::multiclone_bulk_weights`, `::total_cn_all`, `::total_cn_scores`;
> `results/bench_IW/bench_iw.py::score_one`, `::collect`, `::matrix_summary`. Results:
> `results/bench_IW/RESULTS.md`, `RESULTS_intpsi.md`.

bench_IW (`bench_iw.py --set iw`) and bench_IW_intpsi (`--set iw_intpsi`, whole-number ψ per
clone) run four methods per case: `rck_bulk` (bulk route, answer re-ranked at λ = 0.2), `rck_mc`
(multi-clone route), `coral` (CoRAL 3.1.1, flanked graph, §5.4) and `ce` (CycleExtractor 1.0.0).
CoRAL and CE decompose the direct subgraph export (§5.2, `iw`) or `rck_bulk`'s graph with the
subgraph's read paths (§5.5, `iw_intpsi`).

**Scoring one run** (`bench_iw.py score`, into `<run>/score/`):

1. The answer as an AA cycles file: `rck_bulk` via §5.1 at the λ = 0.2 pool entry with 2
   alternatives; `rck_mc` via §5.1, rank 1, no alternatives; `coral` and `ce` with `Path=` read
   as `Cycle=` (`score_h2h.normalise_cycles`). CoRAL and CE are scored twice, linear walks kept
   (`method_variant` `<m>:walks_kept`, `score/`) and dropped (`<m>:walks_dropped`,
   `score_nolinear/`).
2. `aa_eval` (§4) with its defaults, `--case <case>:<method_variant>`, the decomposed graph as
   `--graph` (not for `rck_mc`), and for `rck_mc` `--phase-blocks <DIAG>/rck.fragments.tsv`
   (§3.10). Graph-level metrics need a truth mixture graph (§4.2), which bench_IW does not build.
3. `bench_metrics` (below), with `rck_mc`'s weights converted to bulk copies.

**ecDNA copy-number recall** (`cn_*`; ecDNA only, haplotype-blind, base-pair weighted). For every
reference base $(\mathrm{chr}, x)$:

```math
t(\mathrm{chr},x) = \sum_{T \in \mathcal W^\ast} w_T\, m_T(\mathrm{chr},x),
\qquad
c(\mathrm{chr},x) = \sum_{C \in \mathcal S} w_C\, m_C(\mathrm{chr},x)
```

```math
\mathrm{cn\_recall} = \frac{\sum \min(c, t)}{\sum t},
\qquad
\mathrm{cn\_precision} = \frac{\sum \min(c, t)}{\sum c},
\qquad
\mathrm{cn\_within\_tol} = \frac{\big\lvert \lbrace x : t > 0,\ \lvert c - t \rvert \le \max(1,\ 0.1\,t) \rbrace \big\rvert}{\big\lvert \lbrace x : t > 0 \rbrace \big\rvert}
```

where

- the sums run over all bases;
- $\mathcal W^\ast$ is *every* truth species of the tumour sample (not only the in-scope ones of
  §2.3), with $w_T$ its bulk copies per cell (`PHI`, §2.1);
- $\mathcal S$ is the strict answer (§4.6, `strict_select` with rule `auto`; cycles with
  `Status=chromosomal` are never species);
- $m_X(\mathrm{chr},x) = \sum_H m_X(\mathrm{chr},x,H)$ is the number of times walk $X$ traverses
  the base (§1), on its first ordering;
- $w_C$ is the cycle's `Copy_count`; for `rck_mc` it is $\sum_i u_i\, \psi_{C,i}$, with $u_i$ from
  the DIAG's `clone_proportions.tsv` and $\psi_{C,i}$ the `psi_<i>` columns of rank 1 in
  `decompositions.tsv` (`multiclone_bulk_weights`);
- a ratio with a zero denominator is `null`.

The same over every species cycle (with `low_weight` ones) is `cn_recall_all`, `cn_precision_all`,
`cn_within_tol_all`. `truth_mass`, `inferred_mass`, `overlap_mass` and `truth_bp` are the sums.

**ecDNA proportions.** For a set $X \in \lbrace \mathcal S, \text{all species cycles} \rbrace$:
$p_T = w_T / \sum_{T'} w_{T'}$ over $\mathcal W^\ast$, $p_C = w_C / \sum_{C' \in X} w_{C'}$. Truth
and inferred are matched one to one, greedily by haplotype-blind Jaccard (`evaluate._match`, §3.6),
giving a partner $\mathrm{match}(T)$ or none.

```math
\mathrm{proportion\_tv} = \frac{1}{2}\Big( \sum_{T\ \text{matched}} \lvert p_{\mathrm{match}(T)} - p_T \rvert + \sum_{T\ \text{unmatched}} p_T + \sum_{C\ \text{unmatched}} p_C \Big),
\qquad
\mathrm{proportion\_mae\_matched} = \operatorname{mean}_{T\ \text{matched}} \lvert p_{\mathrm{match}(T)} - p_T \rvert
```

`proportion_tv` ∈ [0, 1], lower is better; `proportion_mae_matched` is `null` with no match. The
`all` set's values carry the suffix `_all`, with `n_matched` / `n_matched_all`.

**Total copy number** (`total_cn_*`; `rck_bulk` and `rck_mc` only; background plus ecDNA,
haplotype-blind). Over $\mathcal D$ = the bases of the method's segments that `truth/segments.tsv`
covers:

- $t(x)$ = `CN_A` + `CN_B` of `truth/segments.tsv` (0-based half-open, read as `END` − 1), the
  tumour sample's cell-averaged copy number;
- $c(x)$ = `inferred_total` of `flow_segments.tsv` (`rck_bulk`), or
  $\sum_i u_i (A_i + B_i) + 2\,(1 - \sum_i u_i)$ over the clones of `rck.scnt.tsv`, its
  `kind=flank` segments left out (`rck_mc`);
- `total_cn_recall`, `total_cn_precision` and `total_cn_within_tol` are as above over $\mathcal D$;
  `total_cn_mae` = $\operatorname{mean}_{\mathcal D} \lvert c - t \rvert$; `total_cn_bases` =
  $\lvert\mathcal D\rvert$, `total_cn_bases_no_truth` = method bases the truth does not cover;
- `total_cn_ecdna_*`: the same over the bases of $\mathcal D$ that a truth species traverses.

$\mathcal D$ is mostly background, so `total_cn_*` is dominated by background and purity;
`total_cn_ecdna_*` is the part where the ecDNA copies matter (`bench_metrics.py` docstring).

**Per-run table** (`bench_iw.py collect`; `bench_iw[_intpsi]_table.tsv`), one row per (case,
method variant): the case fields (`name`, `scenario`, `class`, `S`, `replicate`, `seed`,
`segment_max_traversals`, `species_size_*`, `footprint_bp`), `method`, `method_variant`,
`linear_walks`, `baseline` (`decomposition-only` for CoRAL / CE in `iw_intpsi`), `run_status`,
`elapsed_s` (wall time including licence waits), every `bench_metrics.json` field, and from
`aa_eval`:

- strict sample columns (§4.5): `n_inferred_strict`, `n_missed`, `n_false`, `mean_composite_h0`,
  `mean_composite`, `phase_frame`, `n_missed_absolute`, `n_false_absolute`;
- headline M1 and M2 rows (§4.7): `<metric>_<mode>_<haplotype>_<precision|recall|f1>`;
- from the species matrix (§4.7), per `haplotype` ∈ {`blind`, `phased`}, over its strict columns:
  `<haplotype>_jseg_best_mean` and `<haplotype>_jadj_best_mean`, the mean over truth species of
  the row's largest $J_{\mathrm{seg}}$ and largest $J_{\mathrm{adj}}$ (taken separately, possibly
  from different columns), and `<haplotype>_exact_cycles`, the truth species with an exact strict
  column.

`bench_iw[_intpsi]_pooled.tsv` holds the mean and median (`<column>_mean`, `<column>_median`) of
every numeric column per (`method_variant`, `class`, `S`), (`method_variant`, `class`) and
`method_variant` (`level`), with `n` rows and `n_scored` rows that have `cn_recall`.
`bench_iw[_intpsi]_proportions.tsv` concatenates every run's `bench_proportions.tsv`.

---

## 7. Known limitations

- **`aa_format.read_cycles` reads `Cycle=` lines only.** CoRAL 3.1.1 and CycleExtractor write
  linear walks as `Path=` lines, which the scorer skips unless they are rewritten as `Cycle=`
  (`results/graph_route/score_h2h.py::normalise_cycles` does this; its `--drop-paths` /
  `_cycles_only` tables are the sensitivity check, and bench_IW's `walks_dropped` variant).
- **`aa_eval.read_paths` cannot parse the subgraph's `imputed_ext_*` extremity ids**, so
  `subwalks_satisfied` was not computed for the graph-route arms (RCK `results/graph_route/NEXT_STEPS.md`).
- **`cn-frac:0.1` is not exactly `rck-ecdna`'s significance rule.** `rck-ecdna` flags a species
  `low_weight` below $\max(1,\ 0.1 \times \text{heaviest})$ on *integer* weights, and only when the
  decomposition has more than one species (`rck/ecdna/decompose.py::_flag_low_weight`;
  `--min-species-weight` replaces the threshold). `cn-frac:0.1` keeps cycles with
  `Copy_count` ≥ 0.1 × the heaviest, has no 1-copy floor, and reads `Copy_count`, which for
  `rck-ecdna` is the refitted weight. `rck-ecdna` answers use `status`, so this affects only AA,
  CoRAL and CycleExtractor answers.
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
- **M1 has no true negatives**, so it has no FPR (§4.7). Its composite thresholds include the
  copy-number error, and the 0.95 operating point is the bench's failure threshold, not a published
  value.
- **M2 sees only novel adjacencies.** A molecule with at most one novel adjacency forms no pair, so
  its own grouping is not scored (bench `nest_samehap_3to1_s306`'s nested species has one). Two
  molecules with the same novel adjacencies but different multiplicities (for example a segment
  traversed twice) are the same set of items: M1's `exact` separates them, M2 does not.
- **M2 is strict-only** and is not computed over rank 1's alternative decompositions.
- **The species matrix's J_adj does not see multiplicity**, and its merge and split flags can both
  fire on one merged walk, or on a correct answer whose truth species share a novel adjacency
  (§4.7, "Per-cell species matrix").
- **The species matrix's `phased` mode scores absolute haplotypes** unless `aa_eval` is given the
  phase blocks. The multi-clone route's phase is identified only within phase blocks (§3.10), so a
  correct answer with a flipped block loses J_seg and J_adj in absolute `phased` mode. Since
  2026-11-14, `aa_eval --phase-blocks <DIAG>/rck.fragments.tsv` applies §3.10's per-block swaps
  before every phased score: the strict pairing, `n_missed` / `n_false`, M1, the composites and the
  phased matrix. It keeps the absolute counts as `n_missed_absolute` / `n_false_absolute`, and
  `phase_frame` says which was used. bench_IW_intpsi scores `rck_mc` this way.
- **Multi-clone runs:** `aa_eval` leaves out walks with `Status=chromosomal` (the route's
  telomere-to-telomere paths), so they are neither species nor matrix columns. The route writes no
  `aa/*_graph.txt`, so the graph-level metrics are skipped; M1, M2 and the matrix do not need it.
- **M2's `graph` universe** adds every discordant edge of `--graph`, not only those within the truth
  amplicon. It is meant for the graph route, where that graph is the amplicon's shared input.
- **Unobserved concordant joins** (at imputed grid points, 10 in demo7) get the smaller
  neighbouring segment CN in the exported AA graph; RCK-ecDNA has no observation there.
- **Truth and inferred read paths differ.** The truth graph gets every read of `obs/path.tsv`
  whose junctions are mixture edges (demo7: 12). The subgraph's `paths.tsv` was filtered by
  `rck_inputs.clip_path` (demo7: 3). Synthetic truth reads are `round(CN × R)`, so the fitted
  junction scale is close to, not exactly, R (demo7: 6.497 for 6.5).
- **`bench_metrics` scores every truth species**, not only those in scope (§2.3) as `aa_eval` does,
  so a truth species outside the answer's regions lowers `cn_recall` and `proportion_tv`.

---

## Appendix A. How to run

`rck-ecdna-eval` on a finished run:

```bash
rck-ecdna-eval --decompositions rck_ecdna_out/demo3/amplicon_0/decompositions.tsv \
               --truth-dir simulations/runs/demo3_clone1_ecdna
```

Or pass `--truth-dir` (and optionally `--truth-sample`) to `rck-ecdna` itself, which evaluates
at the end of the run with `--eval-tol-bp`, `--min-truth-fraction` and `--all-truth-species`;
the multi-clone route also scores within phase blocks when the DIAG has `rck.fragments.tsv`, and
writes `evaluation_psi.tsv` (§3.9). Options of `rck-ecdna-eval` (`rck/ecdna/evaluate.py::main`):

| Option | Default | Meaning |
|---|---|---|
| `--decompositions` | required | `decompositions.tsv` written by `rck-ecdna` |
| `--truth-dir` | | simulator run directory (with `truth/` and `obs/`) |
| `--truth-cycles` | | explicit truth table (`species order chrom start end strand hap weight`, 0-based inclusive); used instead of `--truth-dir` |
| `--truth-sample` | first tumor (`rck-ecdna`: `--sample`, else first tumor) | sample in `truth/proportions.tsv` for the truth copy numbers |
| `--tol-bp` (`rck-ecdna --eval-tol-bp`) | 10000 | species tolerance $\tau_{\mathrm{sp}}$ |
| `--include-low-weight` | off | also score species flagged `low_weight` (not in `rck-ecdna`) |
| `--min-truth-fraction` | 0.5 | scope threshold $f_{\min}$ (§2.3) |
| `--all-truth-species` | off | score every truth species |
| `--outdir` | the directory of `--decompositions` | |
| `--phase-blocks` | off | RCK fragments file (`rck_input/rck.fragments.tsv`): also score phasing within phase blocks (§3.10) |

`rck-ecdna-eval` reads the regions for the scope from `summary.json` next to `--decompositions`;
`rck-ecdna` uses its graph's regions.

The method-agnostic scorer:

```bash
python -m rck.ecdna.aa_eval --cycles METHOD_cycles.txt --graph METHOD_graph.txt --truth-dir SIM_RUN \
    [--alt-cycles ALT_cycles.txt ...] [--segment-cn cn.tsv] [--paths SUBGRAPH/paths.tsv] \
    [--strict-rule auto] [--graph-tol 100] [--species-tol 10000] [--case NAME] \
    [--phase-blocks rck.fragments.tsv] [--matrix-plot] --outdir OUT
```

| Option | Default | Meaning |
|---|---|---|
| `--cycles` | required | the method's AA cycles file (its answer) |
| `--alt-cycles` | none | other answers of the method, for best-of only |
| `--graph` | | the method's AA graph (graph-level metrics, `cn_explained`, M2's `graph` universe) |
| `--segment-cn` | | allele copy numbers per segment: `chrom start end [cn] cn_a cn_b`, AA coordinates |
| `--paths` | | long-read paths (subgraph `paths.tsv`), for `subwalks_satisfied` |
| `--truth-dir`, `--truth-cycles`, `--truth-sample` | | as `rck-ecdna-eval` |
| `--mixture-dir` | `TRUTH_DIR/truth/aa` when it holds `mixture_amplicons.tsv` | truth mixture files; without them graph metrics are skipped with a warning |
| `--strict-rule` | `auto` | `auto`, `status`, `all`, `cn-frac:F`, `top-cn:X` (§4.6) |
| `--species-tol` | 10000 | $\tau_{\mathrm{sp}}$ |
| `--graph-tol` | 100 | $\tau_{\mathrm{gr}}$ |
| `--min-truth-fraction`, `--all-truth-species` | 0.5, off | scope (§2.3) |
| `--cn-explained-min-cn` | 5 | $\mathrm{cn}_{\min}$ of `cn_explained_topk` (§4.5) |
| `--case` | the cycles file name | the answer's name in `aa_species_matrix.tsv` |
| `--phase-blocks` | off | RCK fragments file: phased scores within phase blocks (§3.10, §7) |
| `--matrix-plot` | off | also draw `aa_species_matrix_<blind|phased>.png` |
| `--outdir` | required | |

Building the truth mixture (amplicon-allele-analysis):

```bash
python simulations/truth_mixture.py RUN_DIR [RUN_DIR ...] [--sample tumor1] [--outdir DIR] \
    [--reads-per-copy R] [--depth-per-copy D] [--path-constraints]
```

(default `--outdir` `RUN_DIR/truth/aa`; with several runs, `OUTDIR/<run name>`).

Exporting an `rck-ecdna` answer for the scorer (§5.1):

```bash
python -m rck.ecdna.aa_export RUN_DIR OUT [--pool-index I | --lambda L] [--alternatives 2]
```

The graph route translators (§5.2–§5.5):

```bash
# inferred subgraph -> AA graph (+ _segment_cn.tsv, _export.json)
python -m rck.ecdna.aa_export subgraph SUBGRAPH_DIR OUTDIR [--name N] [--amplicon 1] \
    [--imputed comment|drop] [--min-path-junctions 1] \
    [--depth-per-copy D] [--junction-reads-per-copy R] [--baseline-cn C] [--phantom-max-bp 2]
# AA graph -> subgraph folder for rck-ecdna --inputs subgraph
python -m rck.ecdna.aa_import GRAPH OUTDIR [--imputed comment|none|regenerate] \
    [--background-cn auto|X] [--cn-profile FILE] [--sample NAME] [--name NAME] \
    [--depth-per-copy D] [--boundary-min-copies 2] [--boundary-frac 0.25] [--compare SUBGRAPH_DIR]
# flank sequence edges so CoRAL 3.1.1 can use breakpoints at interval ends
python -m rck.ecdna.coral_flank GRAPH OUT_GRAPH [--flank 10000]
# bulk route's graph + the direct export's path constraints (decomposition-only baselines)
python -m rck.ecdna.aa_export baseline BULK_RUN_DIR SUBGRAPH_GRAPH DST
```

`aa_export subgraph`'s scale options override `scale.json` (and are needed without it);
`--baseline-cn` and `--phantom-max-bp` are as in `rck-ecdna`. `aa_import --compare` prints the
differences of §5.3's round trip and exits 1 if there are any.

The bench re-scoring through the AA path: `python results/bench_nested/collect.py --aa-metrics`
(also `--lambda L` for the reported λ, `--validate` for the re-ranking check only; the AA exports go
to `BENCH_NESTED_AA_OUT`), and `results/bench_nested/diag/aa_regression.py`. Head-to-head scoring of
the graph route: `python results/graph_route/score_h2h.py [--jobs 16] [--drop-paths]`, with commands
and caveats in `results/graph_route/README.md`. The set metric (§6.2):
`python results/bench_nested/diag/grow_rules.py RCK_RUNS [--lambda 0.2] [-o rules.tsv]`.

bench_IW (§6.4):

```bash
python -m rck.ecdna.bench_metrics --cycles ANSWER_cycles.txt --truth-dir SIM_RUN --outdir OUT \
    [--truth-sample S] [--strict-rule auto] [--case NAME] \
    [--bulk-run RCK_BULK_RUN | --multiclone-run RUN_DIR --clone-proportions DIAG/clone_proportions.tsv]
python results/bench_IW/bench_iw.py [--set iw|iw_intpsi] score [--jobs 8]
python results/bench_IW/bench_iw.py [--set iw|iw_intpsi] collect
```

`--multiclone-run` needs `--clone-proportions`; `--bulk-run` or `--multiclone-run` adds the
`total_cn_*` fields.

## Appendix B. Output files and columns

**`evaluation.tsv`**: one row per (decomposition rank, truth species in scope). It has every field of
`PairScore` (`truth`, `inferred`, `jaccard_hap`, `jaccard_total`, `segment_recall`,
`segment_precision`, `junction_*`, `phased_junction_*`, `phasing_accuracy`,
`phasing_accuracy_swap_invariant`, `order_match`, `order_match_ordering`, `truth_weight`,
`inferred_weight`, `cn_rel_err`, `inferred_weight_int`, `cn_rel_err_int`, `inferred_topology`,
`topology_match`, `composite`), plus `decomposition_rank`, `objective`, `n_species_inferred`,
`n_low_weight` and `n_species_truth`. With phase blocks it also has the `*_blockwise` columns of
§3.10 (`phasing_accuracy_blockwise`, `phased_junction_f1_blockwise`, `order_match_blockwise`,
`composite_blockwise`, `inferred_blockwise`).

**`evaluation.json`**: `n_decompositions`, `n_truth_species`, `n_truth_species_total`,
`truth_species` (per species: `species`, `weight`, `length_bp`, `n_steps`, `n_junctions`,
`fraction_in_subgraph`, `in_scope`), `first_rank_with_all_orders_matched`,
`best_mean_jaccard_hap`, `best_mean_jaccard_rank`, `top_decomposition`, `tol_bp`,
`min_truth_fraction`, `mean_composite`, `recovery_matrix`, and `note` when no species is in
scope. With phase blocks: `mean_composite_blockwise` and `phase_blocks`; multi-clone route: `psi`
(§3.9). When `rck-ecdna`'s own evaluation fails, it holds only `error` and `traceback`. It is
strict JSON: missing scores are `null`.

**`evaluation_psi.tsv`** (multi-clone route, §3.9): one row per (truth species, clone): `truth`,
`clone`, `inferred`, `truth_psi`, `inferred_psi`, `abs_err`, `truth_proportion`,
`inferred_proportion`.

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
| `*_blockwise` (5 columns) | §3.10, only with phase blocks |

**`truth_cycles.tsv`**: the rebuilt truth (`species order chrom start end strand hap weight`).

**`aa_graph_metrics.tsv`**: one row per scored truth amplicon:

- `amplicon`, `species`, `intervals`, `graph_tol_bp`;
- `n_true_segments`, `n_segments_found`, `n_inferred_segments`, `n_false_segments`;
- `n_true_edges`, `n_edges_found`, `n_inferred_edges`, `n_false_edges`;
- `error_free`, the four `*_merged` segment counts, `error_free_merged`;
- `bp_accuracy`, `novel_adjacency_fnr`, `novel_adjacency_fdr`;
- `is_violations`, `is_violations_truth`;
- `cn_error_mode`, `segment_cn_mae`, `segment_cn_rel_err`, `segment_cn_mae_total`.

**`aa_species_metrics.tsv`**: one row per (mode, truth species):

- `mode`, `truth`, `truth_weight`, `truth_length_bp`, `inferred`, `inferred_weight`,
  `inferred_topology`;
- the §3 columns of `recovery_matrix.tsv` (`segment_recall` … `composite`) and
  `exact_cycle_phased`;
- `inferred_h0`, `segment_jaccard_h0`, `junction_f1_h0`, `order_match_h0`, `cn_rel_err_h0`,
  `composite_h0`;
- `cycle_overlap_coral_code`, `cyclic_lcs_coral_code`, `cycle_overlap_supp`, `cyclic_lcs_supp`,
  `length_error`, `exact_cycle`.

The supplementary-methods versions were named `cycle_overlap` and `cyclic_lcs` before 2026-11-03.

**`aa_sample_metrics.tsv`**: one row per mode (`strict`, `best_of`): `mode`, `strict_rule`,
`phased`, `n_truth`, `n_inferred`, `species_count_error`, `n_missed`, `n_false`,
`mean_composite`, `mean_composite_h0`, `mean_cycle_overlap_coral_code`,
`mean_cyclic_lcs_coral_code`, `mean_cycle_overlap_supp`, `mean_cyclic_lcs_supp`,
`n_exact_cycles`, `phase_frame`, `n_missed_absolute`, `n_false_absolute`, `n_blocks_swapped`,
`cn_explained_top1`, `cn_explained_top3`, `n_subwalks`, `subwalks_satisfied`,
`subwalks_satisfied_support` (§4.5, §4.6), and, with graph metrics, `error_free`,
`error_free_merged`, `bp_accuracy`, `novel_adjacency_fnr`, `novel_adjacency_fdr`,
`segment_cn_mae`, `is_violations` (over the amplicons: all, all, mean, mean, mean, mean, sum).

**`aa_partition_metrics.tsv`** (§4.7): one row per (metric, mode, haplotype, criterion,
tolerance, universe):

| Column | Meaning |
|---|---|
| `metric` | `species_detection` (M1) or `adjacency_pairs` (M2) |
| `mode` | `strict` or `best_of` (M1 only) |
| `haplotype` | `blind`, or `phased` (phased answers only) |
| `criterion` | M1: `exact`, `composite_h0>=x` (blind) or `composite>=x` (phased), x ∈ {0.5, 0.6, 0.7, 0.8, 0.9, 0.95}; M2: `co_membership` |
| `tol_bp` | tolerance: $\tau_{\mathrm{gr}}$ for `exact`, $\tau_{\mathrm{sp}}$ for the composite thresholds; M2 both |
| `universe` | M1: `species`; M2: `native`, or `graph` (blind, with `--graph`) |
| `headline` | the headline rows: M1 strict `exact` (blind and phased); M2 strict at $\tau_{\mathrm{sp}}$ on the `native` universe (blind for comparing methods, phased for `rck-ecdna` against the truth) |
| `n_truth`, `n_predicted` | $\lvert\mathcal W^\ast_{\mathrm{in}}\rvert$ and $\lvert\mathcal S\rvert$ (`n_predicted` empty in best-of) |
| `n_items` | M2: $\lvert U \rvert$ |
| `tp`, `fp`, `fn`, `tn` | the counts; empty where undefined |
| `precision`, `recall`, `fpr`, `f1` | the rates; empty where undefined |

The pooled files of §6.3 have the group keys, `headline`, `headline_pooling`, `n_cases`, the summed
`tp`, `fp`, `fn`, `tn`, and `micro_<rate>`, `macro_<rate>` and `n_cases_<rate>` for each rate.

**`aa_species_matrix.tsv`** (§4.7, "Per-cell species matrix"): one row per (haplotype mode, truth
species, column), long format:

| Column | Meaning |
|---|---|
| `case` | the answer's name (`aa_eval --case`, default the cycles file name) |
| `haplotype` | `blind`, and `phased` for a phased answer |
| `tol_bp` | $\tau_{\mathrm{sp}}$, the tolerance of the items |
| `truth`, `truth_rank`, `truth_weight` | the row: truth species, its rank by copy number, its copy number |
| `column_rank`, `inferred`, `copy_number`, `low_weight` | the column: rank (strict first, then `low_weight`, each by copy number), name, copy number, whether it is outside the strict answer |
| `j_seg`, `j_adj` | the cell |
| `exact_cycle` | §4.4 for this cell (`exact_cycle_phased` in the `phased` rows) |
| `j_adj_pairs`, `j_seg_cn` | optional: M2's pair form on the cell; copy-weighted J_seg |
| `n_items_truth`, `n_items_inferred`, `n_items_shared` | the item counts behind `j_adj` |
| `merge_column`, `split_row` | the flags (strict columns only) |

`aa_eval --matrix-plot` also draws `aa_species_matrix_<blind|phased>.png`: two panels (J_seg,
J_adj), with the `low_weight` columns after a dashed line and NaN in grey
(`python -m rck.ecdna.species_matrix OUTDIR [--haplotype blind|phased] [-o PNG]` draws it from the
table).

**`aa_eval.json`**: `info` (`cycles`, `graph`, `alt_cycles`, `strict_rule`, `phased`,
`species_tol_bp`, `graph_tol_bp`, `n_truth_total`, `truth_in_scope`, `regions`,
`mixture_graphs`), `sample` and `graph` (the rows of the two tables above).

**`bench_metrics.json`** (§6.4): `strict_rule`, `n_truth`, `n_cycles`, `n_strict`; `cn_recall`,
`cn_precision`, `cn_within_tol`, `truth_mass`, `inferred_mass`, `overlap_mass`, `truth_bp`;
`cn_recall_all`, `cn_precision_all`, `cn_within_tol_all`; `proportion_tv`,
`proportion_mae_matched`, `n_matched` and their `_all` versions; with `--bulk-run` or
`--multiclone-run`, `total_cn_recall`, `total_cn_precision`, `total_cn_within_tol`, `total_cn_mae`,
`total_cn_bases`, `total_cn_bases_no_truth` and the same with prefix `total_cn_ecdna_`; `case`
(`bench_iw.py` adds `method`, `variant`).

**`bench_proportions.tsv`** (§6.4): one row per truth species or unmatched inferred cycle and set:
`set` (`strict`, `all`), `truth`, `inferred`, `truth_proportion`, `inferred_proportion`,
`abs_err`, `jaccard_blind` (haplotype-blind Jaccard of the pair), `case` (`bench_iw.py`: also
`method`, `method_variant`).

**Translator side files** (§5): `<name>_amplicon<N>_export.json` (`aa_export subgraph`: counts,
intervals, scale, path statistics, fold-backs, concordant CNs taken from segments, outside and
imputed ids, normalisation log); `import.json` (`aa_import`: counts, path statistics, scale,
background, dropped `source` edges, warnings); `<DST stem>.json` (`aa_export baseline`, §5.5).

## Appendix C. Tests

| Test (RCK `tests/ecdna/`) | Covers | Gurobi |
|---|---|---|
| `test_evaluate.py` | truth reading, several-phasing assignment, species metrics, scope, within-block phasing (§3.10) | no |
| `test_recovery.py` | recovery matrix, composite, `best_*` | no |
| `test_aa_eval.py` | AA parsing, graph metrics incl. `error_free_merged` and `merge_benign_splits`, CoRAL metrics (code and supplement), modes, strict rules, sample metrics, the §5.1 export scoring as `rck-ecdna-eval` | no |
| `test_aa_eval_blocks.py` | `aa_eval --phase-blocks`: absolute by default, within-block pairing and phases, absolute counts kept (§3.10, §7) | no |
| `test_aa_partition.py` | M1 and M2 (§4.7): a merge, a split, a nested pair sharing a novel adjacency (blind and phased), a background SV left unused (graph universe) and one routed into a walk, 10 kb against 100 bp, one-to-one matching, a linear walk, undefined rates, micro and macro pooling, the AA-file path end to end | no |
| `test_species_matrix.py` | the per-cell species matrix (§4.7): the worked example (correct and merged), a split, a single-junction species (pair form NaN), an empty union (NaN), the copy-weighted variant, `low_weight` columns (strict item sets unchanged from M2's; a leftover cannot take a strict item), J_seg equal to the recovery matrix's Jaccard for matched pairs, the AA-file path end to end | no |
| `test_aa_import.py` | AA graph writer, subgraph export (§5.2), import modes incl. `regenerate`, scale fit and fallbacks, fold-backs, round trip; reference tokens written by CoRAL in `tests/ecdna/data/graph_route/` | no |
| `test_coral_flank.py` | flank lengths, short flanks between close runs, residual copy numbers (self-loop twice, negative floored), renumbered constraints, read-back (§5.4) | no |
| `test_aa_export_baseline.py` | `aa_export baseline` (§5.5): shift and sort, outside edges, constraint renumbering and drops, zero-CN counts, intervals and info file | no |
| `test_bench_metrics.py` | `bench_metrics` (§6.4): copy-number recall/precision, traversals, empty answer, proportions (matched, unmatched, strict and all), multi-clone weights, total copy number (half-open truth, bulk, multi-clone with flanks and normal cells, ecDNA restriction) | no |
| `test_demo3.py`, `test_demo4.py` | `rck-ecdna --truth-dir` end to end on shipped subgraphs (demo4: per-subgraph scope, requested sample) | yes |

amplicon-allele-analysis `simulations/tests/test_truth_mixture.py` covers the truth mixture
(demo6 two amplicons; demo7 nested amplicon, node balance, silent boundary merged; path
constraints and CoRAL order). No unit test calls `evaluate.psi_scores` (§3.9). Run named modules
only (`python -m unittest tests.ecdna.test_aa_eval ...`). `unittest discover` also runs the Gurobi
tests.


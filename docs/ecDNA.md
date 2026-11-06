# Allele-specific ecDNA reconstruction (`rck-ecdna`)

`rck-ecdna` is stage 3 of the pipeline described in [docs/pipeline/README.md](pipeline/README.md).
It takes the breakpoint subgraph of one amplicon, built from a sample's observations by
amplicon-allele-analysis
[`simulations/subgraph_construct.py`](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/subgraph_construct.md)
(stage 2), and reports the ecDNA species the data support. Each species is an ordered walk over
haplotype-labelled segments, with a copy number per cell (cell-averaged over the bulk sample) and
alternatives. The model extends RCK, the clone- and haplotype-specific karyotype reconstruction of
Aganezov and Raphael [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020), from linear chromosomes to circular, haplotype-specific
molecules. Its decomposition step is modelled on CoRAL's cycle extraction [Zhu et al. 2024](pipeline/references.md#zhu2024).

```
subgraph folder (nodes.tsv, edges.tsv, scale.json, ...)
   │  §2  normalised diploid interval adjacency graph (DIAG)
   │  §3  copy-number observations (segment, allele-specific, adjacency, flank background)
   ▼
§4  two-layer flow MILP: inferred ecDNA + background copy numbers per segment haplotype
    and per adjacency labelling
   ▼
§5  decomposition MILP: fewest species reproducing the ecDNA flow; pool of alternatives;
    ranking (step penalty), ambiguity report
   ▼
§6  refit LP: fractional species copy numbers against the observations
   ▼
§7  alternating walks; alternative orderings ranked by long-read paths
   ▼
outputs (decompositions.tsv, AA-format cycles, ...) → evaluation (docs/pipeline/evaluation.md)
```

Code lives in `rck/ecdna/`. The entry points are `rck/ecdna_run.py` (`rck-ecdna`),
`rck/ecdna/evaluate.py` (`rck-ecdna-eval`), `rck/ecdna/plot.py` (`rck-ecdna-plot`) and
`rck/ecdna/rerank.py` (`rck-ecdna-rerank`). RCK's own karyotype model (`rck`, `rck/core/*`) is
unchanged apart from the patches and two fixes in §8.

Terms follow RCK [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020): *segment*, *extremity*, *reference adjacency*, *novel
adjacency*, *haplotype labelling* (AA, AB, BA, BB), *segment* and *adjacency copy number*. The code's
names (junction, variant edge, node, ...) are kept where files and columns use them, and are mapped
to these terms in the [glossary](pipeline/glossary.md). Commands, file formats and options are in
the appendices.

**Contents**

1. [Provenance at a glance](#1-provenance-at-a-glance)
2. [Input graph](#2-input-graph-the-normalised-diag)
3. [Copy-number observations](#3-copy-number-observations)
4. [Flow MILP](#4-flow-milp-inferred-copy-numbers)
5. [Decomposition into species](#5-decomposition-into-species)
6. [Refit of species copy numbers](#6-refit-of-species-copy-numbers)
7. [Walks and orderings](#7-walks-and-orderings)
8. [Changes to original RCK](#8-changes-to-original-rck-patches-and-fixes)
9. [Evaluation](#9-evaluation)
10. [Known limitations](#10-known-limitations)
- Appendices: [A install](#appendix-a-install), [B run](#appendix-b-how-to-run),
  [C inputs and options](#appendix-c-inputs-and-options), [D outputs](#appendix-d-outputs),
  [E figures](#appendix-e-figures), [F ionic](#appendix-f-running-on-princeton-ionic),
  [G tests](#appendix-g-tests)

---

## 1. Provenance at a glance

Every step is in exactly one category. **(a)** taken from an existing method; **(b)** a patch
needed to run, which changes no result; **(c)** new or modified for this pipeline. Commits are on
RCK `yubi_dev` unless marked otherwise. The [pipeline overview](pipeline/README.md) has the
same table for all stages.

| Step | Section | Category | Source |
|---|---|---|---|
| DIAG, haplotype labellings, labelling rules, haplotype-aware balance | §2.5, §4.3 | (a) | RCK [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020), `rck/core/ilp_gurobi.py::OptModelMultiClone` @ `4c6c288` |
| Edge-type mapping, phantom contraction, parallel merges, terminals, regions | §2.1–§2.4, §2.6 | (c) new | `c257e54` |
| Junction labellings from read-based `junction_hap`, soft mode | §2.5 | (c) new | `c257e54`, `7d8dff9` |
| Depth per copy γ, junction reads per copy θ | §3.1–§3.2 | (c) new (HATCHet route reads HATCHet's output [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020)) | `c257e54`; mode fix `41f5804` |
| Allele split $c_A = c(1-\mathrm{BAF})$ | §3.3 | (a) | the BAF–allele copy-number relation used by ASCAT [Van Loo et al. 2010](pipeline/references.md#vanloo2010) and HATCHet [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020) |
| Subgraph-only inputs (`scale.json`, `CN_*`, `cn_profile.tsv`) | §3 | (c) new | `a987ffa` |
| Two-layer flow MILP | §4 | (c) modified from RCK `OptModelMultiClone` | `7d8dff9` |
| Background slack at terminals | §4.3 | (a) RCK's telomere relaxation, moved to the background layer | `7d8dff9` |
| Imputed candidate adjacencies with a usage penalty | §4.5 | (c) new | `7d8dff9` |
| Length-capped, $1/\sqrt{c}$-normalised data terms; junction read term | §4.5 | (c) new | `7d8dff9` |
| Linear ends | §4.4 | (c) new, motivated by AA and CoRAL's linear amplicons [Deshpande et al. 2019](pipeline/references.md#deshpande2019) [Zhu et al. 2024](pipeline/references.md#zhu2024) | `10416aa` |
| `bg_balance_weight` | §4.5 | (c) new | `10416aa` |
| Phase-block flips | §4.6 | (c) modified from RCK's fragment allele flip [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020) | `a422c39` |
| Infinite sites (`--enforce-is`) | §4.3 | (a) RCK's extremity exclusivity [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020) | `7d8dff9` |
| Decomposition MILP (fewest species, residual) | §5.1 | (c) modified from CoRAL's `minimize_cycles` [Zhu et al. 2024](pipeline/references.md#zhu2024) [CoRAL 3.1.1 code](pipeline/references.md#coral_code) | `3c82929` |
| Exact linearisation of weight × multiplicity | §5.1 | (a) standard binary-expansion linearisation [Glover 1975](pipeline/references.md#glover1975) | `3c82929` |
| Lazy connectivity cuts | §5.2 | (a) subtour-elimination cuts [Dantzig et al. 1954](pipeline/references.md#dantzig1954), as lazy constraints [Gurobi](pipeline/references.md#gurobi) | `3c82929` |
| Noise-tolerant residual | §5.3 | (c) new | `2fdf73d` |
| Long-read constraints | §5.4 | (c) modified from CoRAL's subpath constraints [Zhu et al. 2024](pipeline/references.md#zhu2024) [CoRAL 3.1.1 code](pipeline/references.md#coral_code) | `c5e0b9c` |
| Alternatives: solution pool, no-good cuts | §5.5 | (a) Gurobi solution pool [Gurobi](pipeline/references.md#gurobi); standard no-good cut | `3c82929` |
| gcd normalisation, `low_weight` flag, ranking by significant species | §5.6–§5.7 | (c) new | `41f5804`, `f84faff` |
| Step penalty, ranking-only | §5.7 | (c) new | `90ac98a` |
| Ambiguity report | §5.8 | (c) new | `90ac98a` |
| Pool re-ranking (`--keep-pool`, `rerank`) | §5.9 | (c) new | `949f796` |
| Amplicon call (`ecDNA` / `linear` / `mixed` / `none`) | §5.10 | (c) new | `10416aa` |
| Refit LP (fractional weights) | §6 | (c) new; fractional weights as in CoRAL [CoRAL 3.1.1 code](pipeline/references.md#coral_code); LP solved with HiGHS [Huangfu & Hall 2018](pipeline/references.md#huangfu2018) | `c5e0b9c` |
| Alternating Euler walk | §7 | (a) Hierholzer's algorithm [Hierholzer & Wiener 1873](pipeline/references.md#hierholzer1873), adapted to alternation `[citation needed: Euler walks without forbidden transitions; Kotzig 1968 proposed, not verified]` | `3c82929` |
| Ordering enumeration and long-read ranking | §7 | (c) new | `3c82929` |
| AA-format graph and cycles output | App. D | (a) AmpliconArchitect's format [Deshpande et al. 2019](pipeline/references.md#deshpande2019); `Haplotypes=` is (c) | `3c82929` |
| `gurobi_env` (threads, logger, disposal), NFS-safe test folders | App. F, G | (b) | `41f5804`, `8f0f810` |
| `rck/version-fix`: current gurobipy, networkx, pyvcf3 | §8 | (b) | `0b41976` (Runpeng Luo) |
| `rck_run.py` `--scnb` path and `MIPGapAbs` | §8 | (c) fix to original RCK; changes results when the option is used | `b26bb96` |

---

## 2. Input graph: the normalised DIAG

RCK represents a genome as a *diploid interval adjacency graph* (DIAG): every segment has two
extremities (its tail and head), each present once per haplotype, and segments are joined by
reference and novel adjacencies, each with a haplotype labelling [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020).
`rck-ecdna` builds one DIAG per amplicon subgraph and normalises it so that every non-terminal
extremity has exactly one segment.

### 2.1 From subgraph edges to DIAG elements

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::build_graph`.

| `edges.tsv` `EDGE_TYPE` | DIAG element | Observation carried |
|---|---|---|
| `reference` | segment (code: segment edge, reference edge) | `WEIGHT` = mean depth; `BAF`, `n_snps`; `CN_TOTAL`, `CN_A`, `CN_B` in subgraph mode |
| `extremity` | reference adjacency between the two extremities at a breakpoint (code: reference join) | `WEIGHT` = `#REF`, reads spanning the reference join |
| `variant` | observed novel adjacency (code: junction, variant edge) | `WEIGHT` = `#NOV`; `NOV_A`, `NOV_B`; `junction_hap` |
| `imputed_variant` | candidate novel adjacency proposed by stage 2's imputation (code: imputed edge) | none (its `WEIGHT` is kept as `imputed_support`, not as evidence) |

A reference adjacency touching an `imputed_ext*` extremity carries no read count.

### 2.2 Phantom segments

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::_contract_phantoms`.

The simulator names both sides of a breakpoint as separate extremities (`x+` and `x+1-`). Between
them sits a segment of 1 bp (`seg_hap = no_span`) whose two ends each have exactly one reference
adjacency, so it only carries reference flow. Every segment of at most `--phantom-max-bp` (default
2) bp with that shape is contracted. Its two reference adjacencies become one, whose read count is
their mean; its base joins the right-hand segment; and any other edge on its tail (head) moves to
the true tail (head) on the same side. `--phantom-max-bp 0` turns this off.

*Why:* without it every breakpoint adds a segment, two extremities and a reference adjacency that
the MILPs must balance, and novel adjacencies attach to the phantom instead of to the segment
they join.

### 2.3 Parallel adjacencies

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::_merge_parallel`.

Edges between the same pair of extremities are merged:

- several reference adjacencies become one, with the mean read count;
- several observed novel adjacencies become one, with the summed read count;
- an imputed candidate parallel to a reference or observed novel adjacency is merged into it
  (it is that adjacency). Its weight is kept only as `imputed_support`. Stage 2 produces such
  edges as `inter` edges at imputed extremities and as anchors that loop back onto their own
  extremity;
- several parallel imputed candidates become one.

### 2.4 Terminals

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::build_graph`.

An extremity is *terminal* when it has no segment (the outer flanks of each region), is a
chromosome-boundary node, has no side (`_+` / `_-`), or has no reference adjacency. Only the
background layer of the flow MILP may start or end at a terminal (§4.3).

### 2.5 Allowed haplotype labellings

> **Provenance:** labelling rules (a) RCK [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020); labels from read evidence (c) new,
> `c257e54`; soft mode (c) new, `7d8dff9`. **Code:** `rck/ecdna/subgraph_io.py::_junction_phasings`.

Each adjacency $j=\{u,v\}$ gets a set $\Phi(j)$ of allowed labellings. Labelling $\varphi=(h_1,h_2)$
joins the copy of $u$ on haplotype $h_1$ to the copy of $v$ on haplotype $h_2$, as in RCK
(`rck/core/structures.py::Phasing`).

- Reference adjacencies: $\Phi=\{AA, BB\}$ (RCK: a reference adjacency never joins homologs).
- Observed novel adjacencies follow stage 2's read-based label `junction_hap`: `A` → $\{AA\}$,
  `B` → $\{BB\}$, `mixed` → $\{AB, BA\}$ (a haplotype crossing), missing → all four.
- Imputed candidates: $\{AA, BB\}$.
- A novel self-loop ($u=v$, a fold-back) can never be AB or BA, as in RCK
  (`OptModelMultiClone.define_constraints_for_novel_labeling`); a `mixed` self-loop falls back to
  $\{AA,BB\}$.

With `--soft-junction-hap`, $\Phi(j)$ is all four labellings for observed novel adjacencies, and the
read-based ones become *preferred*: any other labelling used costs `--soft-junction-penalty`
(§4.5).

### 2.6 Regions

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::ECGraph.regions`.

The *regions* of a graph are the spans its segments cover, split wherever two segments on one
chromosome are more than 1 Mb apart. They define the flank background (§3.5) and which truth
species a run is scored against (evaluation).

---

## 3. Copy-number observations

Two scales turn read counts into copy numbers (`rck/ecdna/scale.py`):

- $\gamma$, **depth per copy**: tumour read depth contributed by one copy of a locus, averaged over
  the bulk sample. Observed segment copy number $\hat c_s = \bar d_s/\gamma$.
- $\theta$, **junction reads per copy**: reads spanning one copy of an adjacency. Observed
  adjacency copy number $\hat c_j = r_j/\theta$, where $r_j$ is `#NOV` (novel) or `#REF`
  (reference).

All copy numbers here are cell-averaged over the sample, so they need not be integers.

### 3.1 Depth per copy γ

> **Provenance:** (c) new, `c257e54` (`_log_mode` fixed in `41f5804`); subgraph route `a987ffa`.
> The HATCHet route reads HATCHet's `best.seg.ucn` [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020).
> **Code:** `rck/ecdna/scale.py::estimate_depth_scale`, `::_gamma_from_hatchet`, `::_log_mode`,
> `::build_observations`.

γ comes from the first available source:

0. **Subgraph folder** (`--inputs subgraph`, or `auto` when `scale.json` exists): γ =
   `scale.json` `depth_per_copy`, fitted by stage 2's `rck_inputs.py` over purity and ploidy from
   genome-wide depth and SNPs (documented in
   [subgraph_construct.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/subgraph_construct.md)).
1. **Override** `--depth-per-copy`.
2. **HATCHet** (`--hatchet-seg best.seg.ucn`): over background segments $g$ at least 1 Mb long and
   more than 1 Mb from every region, the length-weighted median of

   ```math
   \frac{\bar d_g}{\sum_{i} u_i \left(a_{g,i}+b_{g,i}\right)}
   ```

   where $\bar d_g$ is mean tumour depth over $g$, $u_i$ is clone $i$'s proportion and
   $a_{g,i}$, $b_{g,i}$ its allele-specific copy numbers in `best.seg.ucn`. The denominator is the
   cell-averaged total copy number. This is the route without a matched normal (cell lines).
3. **Windows mode** (default without the above): 1 kb windows are aggregated into blocks of
   `--block-bp` (1 Mb). A block is *background* when it has at least half its windows and lies
   more than 1 Mb from every region. It is *balanced* when it also has at least 20 SNPs and mirrored
   BAF ≥ 0.40. Then

   ```math
   \gamma = \frac{\operatorname{mode}_{\log}\{\bar d_B : B \text{ balanced}\}}{c_0}
   ```

   where $c_0$ is `--baseline-cn` (default 2, or `scale.json`'s `background_cn` in subgraph
   mode) and $\operatorname{mode}_{\log}$ is the mode of a histogram of $\log \bar d_B$ (bin width
   0.02) smoothed by a Gaussian kernel of σ = 2 bins. With fewer than 10 balanced blocks every
   background block is used. With fewer than 3 values, or values within one bin width, the
   median is used.

`--depth-mode ratio` instead divides tumour by matched-normal depth window by window. Segment copy
numbers are then $c_0\, \rho_s/\rho_{\mathrm{mode}}$, with $\rho$ the tumour/normal ratio. This
removes coverage bias that both libraries share (GC, mappability). The amplicon-allele-analysis
simulator draws the bias independently per sample, so on simulated data the ratio only adds noise,
and `tumor` is the default.

*Why the mode:* an amplicon sample is mostly unamplified background, and its commonest copy number
on allelically balanced blocks is the cell ploidy. The 2026-10-19 fix replaced `np.convolve(...,
mode="same")`, which returned the wrong length whenever the kernel was longer than the histogram.

### 3.2 Junction reads per copy θ

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/scale.py::estimate_junction_scale`.

1. `scale.json` `junction_reads_per_copy` in subgraph mode, or the override
   `--junction-reads-per-copy`.
2. **Background extremities** (`--extremities` and windows): over extremities $v$ more than 1 Mb
   from every region whose local copy number $\hat c_v = \bar d(v\pm 25\,\mathrm{kb})/\gamma$ lies in
   $[0.5, 8]$,

   ```math
   \theta = \frac{\sum_v \#\mathrm{REF}_v}{\sum_v \hat c_v}
   ```

   It needs at least 5 such extremities.
3. **Fallback** (logged as a warning): the median `#REF` of the subgraph's own reference
   adjacencies divided by $c_0$.

### 3.3 Segment observations and the allele split

> **Provenance:** allele split (a), the BAF relation used by ASCAT [Van Loo et al. 2010](pipeline/references.md#vanloo2010) and HATCHet
> [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020); subgraph-mode columns (c) new, `a987ffa`.
> **Code:** `rck/ecdna/scale.py::build_observations`.

For segment $s$:

```math
\hat c_s = \frac{\bar d_s}{\gamma}, \qquad \hat c_{s,A} = \hat c_s\,(1-\mathrm{BAF}_s), \qquad \hat c_{s,B} = \hat c_s\,\mathrm{BAF}_s
```

where $\bar d_s$ is the mean tumour depth over $s$ (from the windows, or the subgraph's `WEIGHT`
when the windows have no bins there) and $\mathrm{BAF}_s$ is the phase-oriented B-allele fraction
from stage 2. In subgraph mode $\hat c_s$, $\hat c_{s,A}$, $\hat c_{s,B}$ are the `CN_TOTAL`,
`CN_A`, `CN_B` columns. A segment with fewer than `--min-snps-allele` (5) SNPs gets a total-only
observation. `--min-snps-allele 100000000` makes every observation total-only, a diagnostic used
below.

### 3.4 Adjacency observations

$\hat c_j = r_j/\theta$ for every reference or novel adjacency with a read count $r_j$; imputed
candidates have none. This is the *observed adjacency copy number*. The `CN` column of
`edges.tsv` holds the same quantity. The flow's integer value (§4) is the *inferred* adjacency
copy number in RCK's sense.

### 3.5 Flank background

For each region $[a,b]$ on chromosome $c$, the flank copy number is the bin-weighted mean of the
copy number over $[a-1\,\mathrm{Mb},\, a-10\,\mathrm{kb}]$ and $[b+10\,\mathrm{kb},\, b+1\,\mathrm{Mb}]$
(windows, or `cn_profile.tsv` in subgraph mode), else $c_0$. Every segment in the region gets it as
$\hat c^{\mathrm{bg}}_s$, the background prior of §4.5.

### 3.6 Per-phase-block observations (`--block-flip unoriented`)

> **Provenance:** (c) new, `a422c39`. **Code:** `rck/ecdna/scale.py::piece_observations`,
> `rck/ecdna/subgraph_io.py::read_segment_blocks`.

The segment BAF is on the true haplotype scale only where stage 2 could orient every phase block the
segment touches. With `--block-flip unoriented`, each segment with rows in `segment_blocks.tsv` gets
one observation per *piece* $\pi$, its part inside one phase block $b(\pi)$:

- pieces are matched to segments by their parent edge key (`NODE_1 NODE_2 CHROM_1 POS_1 POS_2`),
  through phantom contraction (a contracted phantom's pieces go to the segment that absorbed its
  base);
- the piece total $\hat c_\pi$ is its `CN_TOTAL`, and $(a_\pi,b_\pi)$ its `CN_A`/`CN_B` (or, in
  windows mode or with `--depth-per-copy`, the total split by the piece BAF). A piece with fewer
  than `--min-snps-allele` SNPs is total-only;
- its share $\sigma_\pi = \ell_\pi/\sum_{\pi'\in s}\ell_{\pi'}$ is its fraction of the segment's
  length;
- a piece in a block with orientation 1 or −1 in `phase_blocks.tsv` is *fixed* and taken in the
  oriented frame. A piece in a block with orientation 0, or one missing from `phase_blocks.tsv`, is
  taken in the observed frame, and its A/B may be swapped by the block's flip (§4.6).

Without `segment_blocks.tsv` the run logs a warning and uses the segment observations, as with
`--block-flip off` (the default).

---

## 4. Flow MILP: inferred copy numbers

> **Provenance:** (c) modified from RCK's `OptModelMultiClone` [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020), `7d8dff9`; linear
> ends and `bg_balance_weight` `10416aa`; block flips `a422c39`.
> **Code:** `rck/ecdna/ilp_circular.py::OptModelCircular` (`build`, `_balance`, `_objective`,
> `_infinite_sites`, `_piece_terms`, `block_report`, `linear_end_nodes`). Solved with Gurobi
> [Gurobi](pipeline/references.md#gurobi).

The flow MILP infers integer copy numbers per segment haplotype and per adjacency labelling, split
into two layers: an **ecDNA layer**, strictly balanced so that it decomposes into closed walks, and
a **background layer** for the chromosomal copies under the amplicon.

### 4.1 What RCK solves

For reference, RCK's model (`rck/core/ilp_gurobi.py::OptModelMultiClone` at upstream `4c6c288`,
unchanged by `rck/version-fix`) infers, for each clone $i$, integer copy numbers $c_{i,s,h}$ of
every segment haplotype and $y_{i,j,\varphi}$ of every adjacency labelling. Given input
allele-specific segment copy numbers $\hat c_{i,s,A},\hat c_{i,s,B}$ (from HATCHet, TitanCNA,
Battenberg or ReMixT), it minimises

```math
\min \sum_{i}\sum_{s} L_s\left(\delta_{i,s,A}+\delta_{i,s,B}\right),\quad
\delta_{i,s,A} \ge \left| c_{i,s,A} - \left(b_{f(s)}\,\hat c_{i,s,A} + (1-b_{f(s)})\,\hat c_{i,s,B}\right)\right|
```

with $\delta_{i,s,B}$ defined symmetrically, subject to:

- **balance** at every extremity $v$ and haplotype $h$: $c_{i,s(v),h} = \sum$ (adjacency copy
  numbers occupying $(v,h)$, a self-loop counted twice), relaxed to $\ge$ at allowed telomere
  extremities (`define_constraints_on_nodes`, `get_lin_expr_for_node_haplotype`);
- reference adjacencies only AA or BB; a novel adjacency has at most one labelling over all clones;
  a novel self-loop never AB or BA (`define_constraints_for_*_labeling`);
- at least a fraction $1-\mathrm{fp}$ of the novel adjacencies present (`--run-nas-fp`, default
  0.1; `define_constraints_nov_adjacency_overall_presence`);
- copy-number bounds per segment haplotype (`--scnb`), and optional adjacency groups and
  infinite-sites (extremity-exclusivity) constraints.

where $L_s=\lceil \ell_s/10\rceil$ by default (`--run-segment-length-attr length_10`) and
$b_f\in\{0,1\}$ is the allele flip of the *fragment* $f(s)$ containing $s$. A fragment is a span on
which the input copy-number tool inferred allele-specific copy numbers, so its A/B labelling is
arbitrary and shared by its segments (`docs/Segments.md`, "Segments vs fragments").

### 4.2 Variables

| Symbol | Domain | Meaning | Code |
|---|---|---|---|
| $x^{e}_{s,h}$ | $\mathbb Z\cap[0,U_s]$ | ecDNA-layer copy number of segment $s$ on haplotype $h\in\{A,B\}$ | `xe` |
| $x^{b}_{s,h}$ | $\mathbb Z\cap[0,\beta_{\max}]$ | background copy number | `xb` |
| $y^{e}_{j,\varphi}$ | $\mathbb Z\cap[0,M_j]$ | ecDNA-layer copy number of adjacency $j$ with labelling $\varphi\in\Phi(j)$ | `ye` |
| $y^{b}_{j,\varphi}$ | $\mathbb Z\cap[0,\beta_{\max}]$ | background copy number; only for reference adjacencies and (unless `--no-background-novel`) observed novel ones, never imputed | `yb` |
| $p_{j,\varphi}$ | $\{0,1\}$ | labelling $\varphi$ of novel or imputed $j$ is used | `pv` |
| $t_{v,h}$ | $\mathbb Z\cap[0,4\beta_{\max}]$ | background copies ending at terminal $v$ | `t` |
| $e_{v,h}$ | $\mathbb Z\cap[0,U_{s(v)}]$ | ecDNA-layer copies ending at linear-end extremity $v\in\mathcal L$ (§4.4) | `te` |
| $f_b$ | $\{0,1\}$ | phase block $b$'s A and B are swapped (§4.6) | `flip` |

where $U_s=\lceil 1.5\max(\hat c_s,1)\rceil+3$ (`ub_factor`, `ub_add`); $\beta_{\max}$ is `--bg-max`
(4); $M_j=\max$ of $U_s$ over the segments at $j$'s ends (or $\beta_{\max}$ if none), halved for a
self-loop (`_big_m`). $x_{s,h}=x^e_{s,h}+x^b_{s,h}$ and
$y_j=\sum_{\varphi}(y^e_{j,\varphi}+y^b_{j,\varphi})$ are the inferred segment and adjacency copy
numbers.

### 4.3 Constraints

Let $k_{j,\varphi}(v,h)\in\{0,1,2\}$ be the number of times labelling $\varphi$ of $j$ occupies
extremity $v$ on haplotype $h$ (2 for a self-loop AA or BB at $v$; `Adjacency.hap_at`).

**Balance** (`_balance`). For every extremity $v$ and $h\in\{A,B\}$, if $v$ has a segment $s(v)$:

```math
x^{e}_{s(v),h} = \sum_{j \ni v}\sum_{\varphi\in\Phi(j)} k_{j,\varphi}(v,h)\, y^{e}_{j,\varphi} + e_{v,h},
\qquad
x^{b}_{s(v),h} = \sum_{j \ni v}\sum_{\varphi\in\Phi(j)} k_{j,\varphi}(v,h)\, y^{b}_{j,\varphi} + t_{v,h}
```

and if $v$ has no segment (a flank): $\sum k\,y^e = 0$ and $\sum k\,y^b = t_{v,h}$. Here
$e_{v,h}=0$ unless $v\in\mathcal L$, and $t_{v,h}=0$ unless $v$ is terminal. This is RCK's
haplotype-aware balance, applied per layer. Because the ecDNA layer has no telomeres (apart from
linear ends), it decomposes into closed alternating walks (§7). The terminal slack $t$ plays the
role of RCK's telomere relaxation, for the background layer only.

**Labelling and presence** (`build`). For every novel or imputed $j$ and $\varphi\in\Phi(j)$:

```math
y^{e}_{j,\varphi}+y^{b}_{j,\varphi} \le (M_j+\beta_{\max})\,p_{j,\varphi},\qquad
y^{e}_{j,\varphi}+y^{b}_{j,\varphi} \ge p_{j,\varphi},\qquad
\sum_{\varphi\in\Phi(j)} p_{j,\varphi}\le 1
```

At most one labelling per novel adjacency, as in RCK. RCK's minimum share of novel adjacencies
used is dropped: chromosomal SVs that touch the amplicon may stay unused.

**Infinite sites** (`--enforce-is`, off by default; `_infinite_sites`). For every extremity $v$ and
haplotype $h$: $\sum_{j\ni v,\ j\ \text{novel or imputed}}\ \sum_{\varphi:\,h\in\varphi(v)} p_{j,\varphi}\le 1$.
This is RCK's extremity exclusivity [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020). It is off because chromothripsis-derived
ecDNA violates it.

### 4.4 Linear ends

> **Provenance:** (c) new, `10416aa`. **Code:** `rck/ecdna/ilp_circular.py::OptModelCircular.linear_end_nodes`.

Not every focal amplification is circular. AmpliconArchitect and CoRAL also report linear
amplicons (breakage–fusion–bridge, heavily rearranged or tandem-duplicated), and reach a segment's
end through their source edges [Deshpande et al. 2019](pipeline/references.md#deshpande2019) [Zhu et al. 2024](pipeline/references.md#zhu2024). The ecDNA layer may therefore
start or stop at a set $\mathcal L$ of extremities chosen by `--linear-ends`:

- `boundaries` (default): extremity $v$ of segment $s$ is in $\mathcal L$ when it has no reference
  adjacency, or its reference neighbour has no segment (a region flank), or, for a neighbouring
  segment $s'$, $\lvert \hat c_s-\hat c_{s'}\rvert \ge \max(m, \phi\max(\hat c_s,\hat c_{s'}))$ with
  $m$ = `--boundary-min-copies` (2) and $\phi$ = `--boundary-frac` (0.25). That is where the end of a
  linear amplicon shows in the depth;
- `all`: every extremity of a segment;
- `none`: the strictly circular model.

Each copy-end costs $\kappa_{\mathrm{end}}$ = `--linear-end-penalty` (0.25), so a linear walk costs
0.5 per copy. An imputed closing edge costs `--imputed-penalty` (2.0) once, so with 4 or more copies
an affordable imputed closure wins, while a long gap with no other support is explained as linear.

*Why:* without ends, an amplicon with no closing junction can only be explained by forcing a
circle, an imputed join, or deviations. The first run with linear ends is demo4 in
amplicon-allele-analysis `results/rck_ecdna_plan1/` ("Plan 1"). Neither the commit nor the CHANGELOG
records a before/after score for it.

### 4.5 Objective

> **Provenance:** (c) new, `7d8dff9`; `bg_balance_weight` `10416aa`. **Code:**
> `rck/ecdna/ilp_circular.py::OptModelCircular._objective`, `::_seg_weight`.

```math
\min\ \sum_{s} D_s
+ \beta \sum_s w_s \left|\textstyle\sum_h x^{b}_{s,h} - \hat c^{\mathrm{bg}}_s\right|
+ \varepsilon_{\mathrm{bal}} \sum_{s\in\mathcal S_{AB}} \left|x^{b}_{s,A}-x^{b}_{s,B}\right|
+ \sum_{j:\,\hat c_j\ \text{defined}} \frac{\alpha_j}{\sqrt{\max(\hat c_j,1)}}\left|y_j-\hat c_j\right|
+ \sum_{j\ \text{novel or imputed}}\ \sum_{\varphi\in\Phi(j)} \left(\pi_j+\pi_{\mathrm{soft}}\,[\varphi\notin\Phi^{\ast}(j)]\right) p_{j,\varphi}
+ \kappa_{\mathrm{end}} \sum_{v,h} e_{v,h}
+ \sum_b \kappa_b f_b
```

with the segment data term

```math
D_s =
\begin{cases}
w_s\, n_s \sum_{h} \left|x_{s,h}-\hat c_{s,h}\right| & \text{allele-specific observation}\\[2pt]
w_s\, n_s \left|x_{s,A}+x_{s,B}-\hat c_s\right| & \text{total-only observation}\\[2pt]
\sum_{\pi\in s} w_s\,\sigma_\pi\, n_s\, D_\pi & \text{with pieces (§4.6)}
\end{cases}
\qquad
w_s=\max\!\left(\min\!\left(\tfrac{\ell_s}{\ell_{\mathrm{cap}}},1\right), w_{\min}\right),\quad
n_s=\frac{1}{\sqrt{\max(\hat c_s,1)}}
```

where:

| Symbol | Meaning | Default (option) |
|---|---|---|
| $\ell_s$ | segment length (bp) | |
| $\ell_{\mathrm{cap}}$, $w_{\min}$ | length cap and floor of the segment weight | 50,000; 0.02 (`seg_length_cap`, `min_seg_weight`) |
| $\hat c_s$, $\hat c_{s,h}$, $\hat c_j$, $\hat c^{\mathrm{bg}}_s$ | observations of §3 | |
| $\beta$ | pull of the background towards the flank copy number | 0.1 (`--bg-prior-weight`) |
| $\varepsilon_{\mathrm{bal}}$, $\mathcal S_{AB}$ | tie-break towards a heterozygous background, over segments with allele data | 0.01 (`bg_balance_weight`) |
| $\alpha_j$ | weight of the read term: novel / reference adjacency | 0.5 / 0.25 (`--alpha-novel`, `--alpha-ref`) |
| $\pi_j$ | cost per labelling used: imputed / observed novel adjacency | 2.0 / 0.05 (`--imputed-penalty`, `--novel-penalty`) |
| $\pi_{\mathrm{soft}}$, $\Phi^{\ast}(j)$ | cost of a labelling outside the read-based ones (`--soft-junction-hap` only); the preferred labellings | 1.0 (`--soft-junction-penalty`) |
| $\kappa_{\mathrm{end}}$ | cost per linear copy-end | 0.25 |
| $\kappa_b$ | flip cost (§4.6) | $10^{-6}$ |

Absolute values are linearised with one auxiliary variable each (`_abs_dev`). Gurobi's relative MIP
gap is `--mip-gap` ($10^{-4}$).

*Why these terms.* The background layer and its flank prior let the chromosomal copies under an
amplicon (for example 1A + 1B) stay out of the cycles. Without them those copies would be forced
into the ecDNA layer or charged to the deviation term. Capping the length weight at 50 kb stops a
long background segment from dominating short amplified ones. The $1/\sqrt{\max(\hat c,1)}$
normaliser makes a deviation count relative to the size of the observation, as a Poisson standard
deviation would. The read term lets junction reads, not only depth, decide adjacency copy numbers.
The novel-adjacency cost is small parsimony. The imputed cost makes a candidate edge win only when
it explains several copies. `bg_balance_weight` stops linear copies from swapping haplotype at no
cost (it was added with the linear ends). These weights are heuristic: no published model is the
source, and they were set on demo3/demo4 and kept fixed since.

### 4.6 Phase-block flips (`--block-flip unoriented`)

> **Provenance:** (c) modified from RCK's fragment allele flip $b_f$ (§4.1) [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020),
> `a422c39`. **Code:** `rck/ecdna/ilp_circular.py::OptModelCircular.flippable_blocks`,
> `::_piece_terms`, `::_flip_cost`, `::block_report`.

Stage 2 orients a phase block only from votes of pure-haplotype novel adjacencies at segments that
read as one haplotype. A block without votes keeps an arbitrary A/B labelling, so its allele
observations can contradict the adjacency labellings. RCK has the same problem for its fragments and
gives each fragment a binary flip. Here the flip is per **phase block**, not per copy-number
fragment, and only for blocks the subgraph left unoriented:

- one binary $f_b$ per phase block with orientation 0 (or missing) that has an allele-specific piece,
  shared by its pieces across segments; with `--flip-oriented-penalty P`, also per oriented block;
- a piece with observed $(a_\pi,b_\pi)$ fits the targets

  ```math
  \tau_{\pi,A}=a_\pi+f_{b(\pi)}\,(b_\pi-a_\pi),\qquad \tau_{\pi,B}=b_\pi+f_{b(\pi)}\,(a_\pi-b_\pi),\qquad
  D_\pi=\sum_h \left|x_{s,h}-\tau_{\pi,h}\right|
  ```

  which are linear in $f_b$, so the absolute-deviation constraints stay linear. A total-only piece
  has $D_\pi=\lvert x_{s,A}+x_{s,B}-\hat c_\pi\rvert$. Every piece of a segment shares the segment's
  inferred $x_{s,h}$;
- the segment's weight is split over its pieces by length ($\sigma_\pi$), so each segment's total
  weight and the objective scale are unchanged. $n_s$ uses the segment total, so it does not depend
  on $f_b$;
- $\kappa_b=10^{-6}$ for an unoriented block (`flip_penalty`), `--flip-oriented-penalty` for an
  oriented one. With no evidence either way a block stays as observed, and the adjacency labellings
  decide its orientation through the balance constraints;
- after solving, with every copy number fixed, $\Delta_b$ (`delta_if_reversed`) is the change in the
  objective from reversing $f_b$. An unoriented block with $\Delta_b\le$ `--flip-tol` (0.01) is
  `unresolved` and kept as observed. The statuses are `fixed`, `flipped`, `kept` and `unresolved`.

The flips are fixed after the flow. The refit (§6), `flow_segments.tsv` and the figure points use
the observations with the flips applied.

*Why (demo7).* In demo7 (two nested species on chr4, the inner one on haplotype A and the outer on
B where they overlap), blocks chr4:185 and chr4:186 got no votes, so about 1.6 Mb of `ecdna1` had A
and B swapped. With default settings the flow objective was 39.3, with two imputed edges, and the
rank-1 composites were 0.43 and 0. Total-only observations gave 8.8.

*Effect (demo7 union, `rck_out_union_blockflip/` in amplicon-allele-analysis):* the flow flips
chr4:185 and chr4:186 and keeps chr4:171, chr4:172 and chr4:187 as oriented. The 12 blocks of the
unamplified 172.4–185.5 Mb gap are `unresolved`. The objective is 9.13, with no imputed edge.
Both nested species are recovered as separate circles: 16.43 copies against a true 16.90 and 9.08
against 9.10, composites 0.994 and 0.9996. Correctly phased, the two species share 9 reference
segments but only one segment haplotype, so a single merged walk can no longer reproduce the
haplotype-specific flow. `--block-flip off` reproduces the earlier output exactly.

### 4.7 Differences from RCK, in one table

| | RCK `OptModelMultiClone` | `rck-ecdna` `OptModelCircular` |
|---|---|---|
| copy numbers | per clone; integer | one bulk sample; integer, two layers (ecDNA + background) |
| telomeres | allowed at listed telomere extremities ($\ge$ balance) | ecDNA layer: only at linear ends, each copy-end costs `--linear-end-penalty`; background: at terminals |
| novel adjacency use | at least $1-\mathrm{fp}$ of them used | no requirement; each used costs `--novel-penalty` |
| missing adjacencies | n/a | imputed candidates, each costing `--imputed-penalty`; background may not use them |
| allele flips | per input fragment, free | per unoriented phase block, $10^{-6}$; oriented blocks fixed unless `--flip-oriented-penalty` |
| objective | length-weighted L1 to input allele copy numbers | length-capped, $1/\sqrt{c}$-normalised L1 to allele-specific (or total) segment observations, plus adjacency read term, background prior and penalties |
| infinite sites | on | optional (`--enforce-is`) |

---

## 5. Decomposition into species

> **Provenance:** (c) modified from CoRAL's cycle decomposition, `minimize_cycles` [Zhu et al. 2024](pipeline/references.md#zhu2024)
> [CoRAL 3.1.1 code](pipeline/references.md#coral_code), `3c82929`. **Code:** `rck/ecdna/decompose.py::CycleDecomposer.build`,
> `::_read_terms`, `::_callback`, `::solve`.

The ecDNA layer of the flow is split into the fewest species that reproduce it. Splitting a flow
into the fewest weighted paths or cycles is NP-hard in general [Vatinlen et al. 2008](pipeline/references.md#vatinlen2008). Like CoRAL, the
decomposition is solved as a MILP with a bounded number of slots.

### 5.1 The decomposition MILP

The input is the positive ecDNA-layer flow $F_e$ on **phased edges** $e\in E$
(`FlowSolution.ecdna_edges`, `rck/ecdna/phased.py`): segment haplotypes $(s,h)$ with
$F=x^e_{s,h}$, adjacency labellings $(j,\varphi)$ with $F=y^e_{j,\varphi}$, and linear ends
$(v,h)$ with $F=e_{v,h}$. Their endpoints are *phased extremities* $(v,h)$. Linear ends all lead to
one virtual end node. There are $K$ = `--max-species` (4) slots; $W=\max_e F_e$; $\bar m$ =
`--max-multiplicity` (3); $n_{\mathrm{bit}}=\lceil\log_2(\bar m+1)\rceil$.

| Symbol | Domain | Meaning |
|---|---|---|
| $u_k$ | $\{0,1\}$ | slot $k$ is used |
| $w_k$ | $\mathbb Z\cap[0,W]$ | integer copy weight of slot $k$ |
| $\beta_{keb}$ | $\{0,1\}$ | bit $b$ of the multiplicity of edge $e$ in slot $k$ |
| $m_{ke}=\sum_b 2^b\beta_{keb}$ | $\{0,\dots,\bar m\}$ | multiplicity: times species $k$ traverses $e$ |
| $\zeta_{keb}$ | $[0,W]$ | $=w_k\beta_{keb}$ (linearised) |
| $y_{ke}$ | $\{0,1\}$ | $e$ is in slot $k$'s support |
| $\lambda^{\mathrm{lin}}_k$ | $\{0,1\}$ | slot $k$ is a linear walk |
| $r_e$ | $\ge 0$ | residual on edge $e$ |

Constraints, for every slot $k$:

```math
u_k \le w_k \le W u_k,\qquad
\zeta_{keb}\le W\beta_{keb},\quad \zeta_{keb}\le w_k,\quad \zeta_{keb}\ge w_k-W(1-\beta_{keb})
```

```math
m_{ke}\le \bar m\, y_{ke},\qquad y_{ke}\le m_{ke},\qquad y_{ke}\le u_k,\qquad \sum_e y_{ke}\ge 2u_k
```

```math
m_{k,\mathrm{seg}(x)} = \sum_{e\ \text{adjacency-like at}\ x} m_{ke}\quad\text{for every phased extremity } x,
\qquad \sum_{e\ \text{linear end}} m_{ke}=2\lambda^{\mathrm{lin}}_k,\quad \lambda^{\mathrm{lin}}_k\le u_k
```

```math
w_k\ge w_{k+1},\qquad u_k\ge u_{k+1}
```

and the residual and objective:

```math
r_e \ge \left|\sum_k\sum_b 2^b\zeta_{keb} - F_e\right| = \left|\sum_k w_k m_{ke}-F_e\right|,\qquad
\min\ \sum_e \omega_e r_e + \rho\sum_k u_k + \mu\sum_{k,e} m_{ke} + P_{\mathrm{reads}}
```

where $\mathrm{seg}(x)$ is the segment haplotype at phased extremity $x$; "adjacency-like" edges are
adjacency labellings and linear ends, a self-loop counted twice; $\rho$ = `--species-penalty`
(1.0); $\mu=10^{-4}$ (`multiplicity_penalty`), a tie-break towards the walk with fewer
traversals; $\omega_e$ is the residual weight of §5.3; and $P_{\mathrm{reads}}$ is §5.4. The balance
row makes each slot an alternating multiset (each phased extremity has as many adjacency traversals as
segment traversals). A slot has no linear end (a circle) or exactly two (one linear walk). The
product $w_k m_{ke}$ is linearised exactly through the binary expansion of $m_{ke}$
[Glover 1975](pipeline/references.md#glover1975). The ordering constraints break the symmetry between slots.

**Differences from CoRAL's `minimize_cycles`** (`coral/cycle_decomposition.py::minimize_cycles`,
`coral/models/concrete.py::get_minimize_objective`, `coral/models/constraints.py`, CoRAL 3.1.1):

| | CoRAL 3.1.1 | `rck-ecdna` |
|---|---|---|
| graph | AA breakpoint graph, total copy number | phased edges (segment haplotypes, adjacency labellings) |
| weights $w$ | continuous ($\ge 0$) | integer in the MILP; fractional after the refit (§6) |
| multiplicity $x$ | integer, 0–10 | integer, 0–`--max-multiplicity` (3) |
| fit | $\sum_i w_i x_{ei}\le \mathrm{CN}_e$ (quadratic); cover ≥ `p_total_weight` (0.9) of length-weighted copy number and ≥ `p_bp_cn` (0.9) of every discordant edge's | $\lvert\sum_k w_k m_{ke}-F_e\rvert$ penalised; exact linearisation |
| objective | number of cycles minus the covered length-weighted fraction | residual + $\rho$ × species + tie-break (+ reads) |
| connectivity | spanning-tree / node-order constraints | lazy cuts (§5.2) |
| ends | source/sink edges $s$, $t$ | linear ends at $\mathcal L$, through one virtual end node |

### 5.2 Connectivity

> **Provenance:** (a) subtour-elimination cuts [Dantzig et al. 1954](pipeline/references.md#dantzig1954), added lazily [Gurobi](pipeline/references.md#gurobi).
> **Code:** `rck/ecdna/decompose.py::CycleDecomposer._callback`, `::_components`.

Balance alone allows a slot's support to fall apart into several closed walks. At every integer
solution (`MIPSOL`), each slot's support is split into connected components. For each pair
$(C_i, C_j)$ of components of slot $k$, a lazy cut

```math
\sum_{f\in\delta(C_i)} y_{kf} \ \ge\ y_{k,e_i}+y_{k,e_j}-1
```

is added, where $\delta(C_i)$ is the set of edges with exactly one endpoint in $C_i$ and $e_i\in C_i$,
$e_j\in C_j$. A connected, balanced, alternating multiset is one alternating closed walk (§7).

### 5.3 Residual weights, including the noise-tolerant modes

> **Provenance:** default weights (c) new, `3c82929`; noise-tolerant modes (c) new, `2fdf73d`.
> **Code:** `rck/ecdna/decompose.py::CycleDecomposer._omega`, `::ranking_key`.

```math
\omega_e =
\begin{cases}
\max\!\left(\min\!\left(\tfrac{\ell_s}{\ell_{\mathrm{cap}}},1\right), 0.05\right)\cdot \sigma_{\mathrm{seg}} & e=(s,h)\ \text{segment}\\
\omega_{\mathrm{adj}} & \text{otherwise}
\end{cases}
\qquad\text{and, with }\texttt{--residual-mode relative:}\quad
\omega_e \leftarrow \frac{\omega_e}{\sigma\,\max(F_e,1)}
```

where $\ell_{\mathrm{cap}}$ = 50 kb, $\sigma_{\mathrm{seg}}$ = `--seg-residual-scale` (1),
$\omega_{\mathrm{adj}}$ = `--adj-residual-weight` (1), and $\sigma$ = `--residual-sigma` (0.1). The
relative mode measures the residual in units of the expected noise. It includes adjacencies because
the flow is balanced, so a dip on a segment comes with dips on its adjacencies. Below
$1/\sigma$ = 10 copies the relative weight exceeds the absolute one.

`--rank-residual sum` ranks on residual + structure cost (§5.7) instead of residual first. Lower
weights alone only change what reaches the pool; with the default ranking an exact fit in the pool
still ranks first.

*Why (bench_nested v2, step 2 diagnostic).* By default one copy of residual on a segment of 50 kb or
more costs as much as one species. A 1–3 copy dip from depth bias is then cheaper to cover exactly
with a leftover or a merged species. *Effect:* on the four decomposition failures (`s301`, `s703`,
`s505`, `s403`), no setting tried changed rank 1 or added a better pool answer
(`results/bench_nested/diag/rck_noise/SUMMARY.md`). The options stay opt-in, and defaults and their
outputs are unchanged.

### 5.4 Long-read constraints

> **Provenance:** (c) modified from CoRAL's subpath constraints [Zhu et al. 2024](pipeline/references.md#zhu2024)
> (`coral/models/constraints.py::get_shared_subpath_edge_constraints`,
> `::get_addtl_subpath_edge_constraints`), `c5e0b9c`. **Code:**
> `rck/ecdna/decompose.py::read_constraints`, `::CycleDecomposer._read_terms`, `::explained_reads`.

With `--paths`, a long-read path $p$ becomes a read constraint when it spans at least
`--min-path-junctions` (2) novel adjacencies, has at least `--min-path-support` (1) reads, and every
one of its adjacencies carries ecDNA flow. Identical paths, in either direction, are merged and their
support $s_p$ summed. For each adjacency $g$ of $p$, $G_g\subseteq E$ is the set of its labellings
with flow. Then

```math
o_{kg}=\bigvee_{e\in G_g} y_{ke},\qquad z_{kp}=\bigwedge_{g\in p} o_{kg},\qquad x_p=\bigvee_k z_{kp}
```

each written as linear constraints on binaries, and

```math
P_{\mathrm{reads}} = \pi_{\mathrm{path}}\left(1-\frac{\sum_p s_p x_p}{\sum_p s_p}\right),\qquad
\text{optionally}\quad \sum_p s_p x_p \ \ge\ \phi_{\mathrm{path}}\sum_p s_p
```

with $\pi_{\mathrm{path}}$ = `--path-weight` (0.5) and $\phi_{\mathrm{path}}$ = `--path-fraction`
(0, off). The reward is capped at 0.5, below one species, so reads choose between decompositions
that fit the flow equally well but never add a species on their own. If no decomposition satisfies
`--path-fraction`, the run logs a warning and solves without it. $z$ and $x$ are fixed by the slot
supports, so they do not multiply the solution pool.

Containing every adjacency of a read is necessary, not sufficient, for the read to fit a walk.
After solving, each read is checked in order against each species' top ordering
(`rck/ecdna/walk.py::path_matches`). `reads_explained` is the share of support found.

**Differences from CoRAL.** CoRAL constrains every edge of a path with its count
($x_{ei}\ge r_{pi}\,\mathrm{count}_e$). In its exact model every path constraint must be satisfied.
In its greedy and post-processing models a fraction `p_subpaths` (0.9) of the *number* of path
constraints must be. `rck-ecdna` requires only the path's novel adjacencies to be in one species,
weighs paths by read support, and checks order afterwards. `--path-fraction 0.9` is therefore a
share of read support, not of path count.

### 5.5 Alternatives

> **Provenance:** (a) Gurobi's solution pool [Gurobi](pipeline/references.md#gurobi) and a standard no-good cut, `3c82929`.
> **Code:** `rck/ecdna/decompose.py::CycleDecomposer.solve`.

- `--alternatives-mode pool` (default): Gurobi `PoolSearchMode=2`, `PoolSolutions =
  max(20, 4 × --num-alternatives)`, `PoolGap` = `--pool-gap` (0.10). With `--step-penalty` > 0 the
  pool is widened to an absolute window: `PoolGap` = ∞, `PoolGapAbs` = `--ambiguity-window` (or
  2ρ if that is 0), `PoolSolutions = max(100, …)` (§5.7).
- `--alternatives-mode nogood`: up to `--num-alternatives` re-solves, each adding the cut
  $\sum_{(k,e):\,y_{ke}=1}(1-y_{ke})+\sum_{(k,e):\,y_{ke}=0}y_{ke}\ge 1$, stopped when the objective
  leaves the window.

Every solution is re-checked for balance, connectivity and 0 or 2 ends, and deduplicated by its
multiset of (weight, edge multiplicities).

### 5.6 Post-processing: repeated circles and `low_weight`

> **Provenance:** (c) new, `41f5804`. **Code:** `rck/ecdna/decompose.py::CycleDecomposer._read`,
> `::_flag_low_weight`.

- **Repeated circle.** A circular slot whose multiplicities share a divisor $g>1$ goes round the same
  circle $g$ times. It is reported as that circle at weight $g\,w_k$.
- **`low_weight`.** With $n>1$ species, species $k$ is flagged `low_weight` when

  ```math
  w_k < \tau,\qquad \tau =
  \begin{cases}
  \texttt{--min-species-weight} & \text{if set}\\
  \max\!\left(1,\ \phi_{\mathrm{lw}}\max_{k'} w_{k'}\right) & \text{otherwise}
  \end{cases}
  ```

  with $\phi_{\mathrm{lw}}$ = `--min-species-weight-frac` (0.10). A flagged species is kept in the
  output but is usually flow residue, not a molecule. The ranking (§5.7) and the evaluation ignore it
  (`rck-ecdna-eval --include-low-weight` scores it).

### 5.7 Ranking and the step penalty

> **Provenance:** ranking by significant species (c) new, `f84faff`; step penalty (c) new,
> `90ac98a`. **Code:** `rck/ecdna/decompose.py::ranking_key`, `::structure_cost`, `::pool_keeps`,
> `::prune_pool`, `::decompose_with_pool`.

Each pool solution $D$ (integer weights) gets

```math
R(D)=\sum_e \omega_e\left|\sum_{k\in D} w_k m_{ke}-F_e\right|,\qquad
n_{\mathrm{sig}}(D)=\#\{k\in D: k\ \text{not flagged low-weight}\},\qquad
S(D)=\sum_{k\ \mathrm{sig}}\ \sum_{(s,h)} m_{k,(s,h)}
```

```math
C(D)=\rho\, n_{\mathrm{sig}}(D)+\lambda\, S(D)
```

where $S$ counts the **steps** of the significant species, one per segment traversal (a walk's
length in segments, not weighted by copies), and $\lambda$ = `--step-penalty` (default 0).
Decompositions are sorted lexicographically by

```math
\bigl(R,\ C,\ n_{\mathrm{sig}},\ n,\ \mathrm{objective}\bigr)\quad\text{(default)},\qquad
\bigl(R+C,\ C,\ n_{\mathrm{sig}},\ n,\ \mathrm{objective}\bigr)\quad(\texttt{--rank-residual sum})
```

and cut to `--num-alternatives` (10). Rank 1 is the answer.

With $\lambda>0$ the widened pool is first pruned (`pool_keeps`). Solution $D$ is kept if

```math
\mathrm{obj}(D)\le \mathrm{obj}^\ast + \lvert\mathrm{obj}^\ast\rvert\,g_{\mathrm{pool}}
\quad\text{or}\quad
\lambda S(D)+\mathrm{obj}(D)-\mathrm{obj}^\ast < \lambda S^\ast
```

where $\ast$ is the best pool solution and $g_{\mathrm{pool}}$ is `--pool-gap` (0.10). That is, it is kept if it is within the usual gap, or if
the prior can prefer it.

*Why ranking by significant species* (`f84faff`): "one molecule plus weight-1 leftovers" should beat
"the molecule split in two" when the solver scores them equally (demo3).

*Why the step penalty.* Species that share reference sequence, such as nested ecDNA, can often be
explained in several ways that reproduce the flow exactly. demo8 has an outer cycle at 19.5 copies
and a nested one at 6.5 on the same haplotype. "Outer 18 + nested 6" and "outer 12 + one walk
through outer and nested at 6" predict the same copy number everywhere (outer-only segments
$18=12+6$, shared segments $18+6=12+2\times 6$). Neither a better fit nor more samples can separate
them; only reads that span the difference, or a prior, can. In demo7 with total-only observations,
one walk round the outer cycle twice and the nested one once needed one species fewer, and the
species penalty chose it. The step penalty is a description-length prior on walks:

- splitting a species into pieces keeps the total steps about the same, so it does not push towards
  splits;
- a merged walk re-traverses the shared sequence and pays for the extra steps (demo8: 47 against 28
  for the separate circles);
- merging saves one species ($\rho$) and costs the extra steps, so it is undone once $\lambda\times$
  extra steps $>\rho$;
- a true internal duplication has the same steps as its split into two species, so $\rho$ keeps it
  whole.

*Why only in the ranking.* The first design put the penalty in the MILP. There it also charged the
light leftover species that absorb unevenness in the integer flow. Their walks are often long (45
steps on demo3), so the solver reshaped real species to shorten them: it split demo3's duplication
(composite 0.96 → 0.72 at λ ≥ 0.02) and demo7's tandem duplication with block flips (`ecdna2` 1.00 →
0.77 at λ ≥ 0.1; `results/step_sweep/sweep_v1_milp_penalty.tsv`). A weight threshold inside the MILP
was gamed by splitting a merged walk into copies too light to pay. With λ > 0 the pool is widened
(§5.5) so that alternatives with more species reach the ranking.

*Effect* (`results/step_sweep/sweep.tsv`, rank-1 composite per truth species):

| run | λ = 0 | λ = 0.02 | λ ≥ 0.05 |
|---|---|---|---|
| demo3 | 0.960 | 0.959 | 0.959 (weight 12.77 → 12.70, an equivalent solution) |
| demo4 amplicon_0 / amplicon_1 | 0.983 / 0.997 | same | same |
| demo5 | 0.990 | same | same |
| demo6 union | 0.984, 0.994 | same | same |
| demo7 union, `--block-flip unoriented` | 0.994, 1.000 | same | same |
| demo7 union, total-only | 0.510, 0.335 (merged walk, 67 steps) | same | **0.982, 0.971** (two species, 31 steps) |
| demo8 union (3:1, one haplotype) | 0.985, 0.989 (tie-break, gap 0.004) | same | same, decided by the step penalty (gap 0.95 at 0.05, 1.90 at 0.1) |

demo7's total-only merge needs λ > 1/36 (merging saves one species and adds 36 steps). The wider
pool costs time: demo8 went from 8 s to 31 s. On the 40-case `results/bench_nested/` benchmark,
λ = 0.2 is the smallest value on the plateau of best mean composite, with no control or demo case
losing against λ = 0, and is the recommended setting for nested amplicons. RCK's default stays 0, so
earlier results are unchanged.

### 5.8 Ambiguity report

> **Provenance:** (c) new, `90ac98a`. **Code:** `rck/ecdna/decompose.py::CycleDecomposer.same_fit_alternatives`,
> `::same_fit_structures`, `::ambiguity_summary`, `::ambiguity_report`.

After ranking, one more MILP (§5.1 with two additions) looks for alternatives the data cannot tell
apart from the top decomposition $D^\ast$:

```math
\sum_k w_k m_{ke} = \sum_{k\in D^\ast} w^\ast_k m^\ast_{ke}\quad\forall e,\qquad
\sum_e \left|y_{ke}-y_{le}\right| \ge u_l\quad\forall k<l
```

The first fixes the flow each edge covers, so the fit is identical. The second forbids two species
on the same support, which excludes one walk split into copies. Solutions within
`--ambiguity-window` (2.0; 0 turns the report off) of the best are enumerated (up to 50). Each
distinct set of significant structures other than the top one's is an alternative $A$ with

```math
\mathrm{gap}(A) = \bigl(\mathrm{obj}(A)+\lambda S(A)\bigr)-\bigl(\mathrm{obj}(D^\ast)+\lambda S(D^\ast)\bigr)
```

The report (`summary.json` `ambiguity`) gives `ambiguous`, `n_alternatives`, the closest
alternative's `gap`, and `decided_by`:

- `species_penalty` if the closest alternative has another number of significant species;
- else `step_penalty` if its steps differ and λ > 0;
- else `tie-break` if its steps differ;
- else `reads or tie-break`.

It also lists the first 10 alternatives. A negative gap means the ranking preferred fewer species
over a lower score. An ambiguous top decomposition logs a warning. The extra MILP adds
$K(K-1)/2$ binaries per edge.

*Why:* demo8's correct answer won by 0.004 through the $10^{-4}$ edge-traversal tie-break, and
nothing in the output said so.

### 5.9 Offline re-ranking over λ (`--keep-pool`, `rerank`)

> **Provenance:** (c) new, `949f796`. **Code:** `rck/ecdna/decompose.py::rank_pool`, `::_pool_record`,
> `rck/ecdna/rerank.py`.

λ acts only in the ranking, and every λ > 0 widens the pool to the same window, so the
decomposition MILP and its pool do not depend on λ. The pruning of the widened pool does depend on
λ. `--keep-pool` (it needs λ > 0) writes `pool/`: every pool solution before pruning, annotated,
refitted and evaluated like a ranked decomposition, plus, for each solution that is rank 1 at some λ
in `--pool-lambdas`, its same-fit alternatives. `python -m rck.ecdna.rerank RUN --lambdas …` applies
the run's own `pool_keeps` and `ranking_key` per λ. λ = 0 is only emulated (solutions within
`--pool-gap` of the best): a λ = 0 run solves a narrower pool whose ties can come out differently
(demo3: 12.77 copies at λ = 0, an equivalent 12.70 at λ > 0). *Effect:* on the bench, re-ranked
composites and ambiguity fields matched direct runs at 32/32 (case, λ), so each case needs two runs
(λ = 0 and one pool run) instead of one per λ.

### 5.10 Amplicon call

> **Provenance:** (c) new, `10416aa`. AmpliconArchitect and CoRAL also classify amplicons
> [Deshpande et al. 2019](pipeline/references.md#deshpande2019) [Zhu et al. 2024](pipeline/references.md#zhu2024); the 80% rule here is this pipeline's own.
> **Code:** `rck/ecdna/decompose.py::classify`.

From the significant species of a decomposition, with $L_k$ the walk length in bp:

```math
f_{\mathrm{circ}}=\frac{\sum_{k\ \mathrm{circular},\ w_k\ge c_{\min}} w_k L_k}{\sum_{k\ \mathrm{circular},\ w_k\ge c_{\min}} w_k L_k+\sum_{k\ \mathrm{linear}} w_k L_k}
```

The call is `ecDNA` if $f_{\mathrm{circ}}\ge\phi_{\mathrm{ec}}$, `linear` if
$f_{\mathrm{circ}}\le 1-\phi_{\mathrm{ec}}$, `mixed` in between, and `none` when the denominator is
0. Here $c_{\min}$ = `--min-ecdna-copies` (2) and $\phi_{\mathrm{ec}}$ = `--ecdna-call-fraction`
(0.8). It is recomputed after the refit.

---

## 6. Refit of species copy numbers

> **Provenance:** (c) new, `c5e0b9c`. Fractional cycle weights as in CoRAL (`w` is continuous in
> `coral/models/concrete.py::CycleLPModel`) [CoRAL 3.1.1 code](pipeline/references.md#coral_code); the LP is solved by HiGHS
> [Huangfu & Hall 2018](pipeline/references.md#huangfu2018) through SciPy [Virtanen et al. 2020](pipeline/references.md#virtanen2020). **Code:** `rck/ecdna/refit.py::observation_rows`,
> `::_solve_lp`, `::refit_weights`, `::refit_decompositions`.

The decomposition reproduces the integer flow, which is fitted edge by edge to noisy observations
and can be uneven (15 copies on most edges of a molecule, 28 on one segment). Matching it exactly
needs leftover species, and an integer leftover can take copies from the real molecule (demo4
tumor2: 6 + 1 against a true 7.44). Each decomposition therefore keeps its structures $m_{ke}$ and
gets real weights $w_k\ge 0$ fitted directly to the observations:

```math
\min_{w\ge 0}\ \eta\sum_k w_k + \sum_{r} c_r\left|\sum_k a_{rk}\, w_k - t_r\right|,\qquad
a_{rk}=\sum_{e\in E_r} m_{ke}
```

The rows $r$ are the flow MILP's own data terms with the background held at the flow solution:

| Row | $c_r$ | $t_r$ | $E_r$ |
|---|---|---|---|
| segment haplotype $(s,h)$, allele data | $w_s n_s$ | $\hat c_{s,h}-x^b_{s,h}$ | $\{(s,h)\}$ |
| segment, total-only | $w_s n_s$ | $\hat c_s-\sum_h x^b_{s,h}$ | $\{(s,A),(s,B)\}$ |
| piece (with `--block-flip`, flips applied) | $w_s\sigma_\pi n_s$ | $\tau_{\pi,h}-x^b_{s,h}$ or $\hat c_\pi-\sum_h x^b_{s,h}$ | as above |
| adjacency with reads | $\alpha_j/\sqrt{\max(\hat c_j,1)}$ | $\hat c_j-\sum_\varphi y^b_{j,\varphi}$ | its labellings |

with $\eta=10^{-3}$ (`weight_penalty`). If the lightest species is refitted below
`--refit-min-copies` (0.5), it is dropped and the LP solved again (the last species is never
dropped). The ranking is not changed. Weights are rounded to 2 digits, `low_weight` flags and the
amplicon call are recomputed from them, and `data_objective` / `data_objective_int` report the misfit
with the refitted and integer weights. `--no-refit` reports the integer weights. The LP needs no
Gurobi licence.

---

## 7. Walks and orderings

> **Provenance:** alternating circuit (a) Hierholzer's algorithm [Hierholzer & Wiener 1873](pipeline/references.md#hierholzer1873), adapted to
> alternating (segment / adjacency) walks `[citation needed: Euler walks without forbidden transitions; Kotzig 1968 proposed, not verified]`; enumeration and read ranking (c) new,
> `3c82929`. **Code:** `rck/ecdna/walk.py::alternating_circuit`, `::enumerate_circuits`,
> `::canonical`, `::path_support`, `::path_matches`.

A species is a balanced, connected multiset of phased edges. A walk must alternate: after arriving
at an extremity on an adjacency it leaves on that extremity's segment, and vice versa. Every phased
extremity has exactly one segment and as many adjacency traversals as segment traversals, so such
a closed walk exists.

- `alternating_circuit` builds a greedy closed trail, then splices sub-trails in before any segment
  traversal whose segment still has copies left. Splicing before a segment step keeps alternation.
- A linear species is closed with a pseudo segment between its two ends, extracted as a circuit
  starting on it, and opened again.
- `enumerate_circuits` runs a bounded depth-first search (at most `--max-orderings`, 20, and
  200,000 expansions) for distinct orderings of the same multiset, up to rotation and reversal (only
  reversal for a linear walk).
- Orderings of one multiset have identical copy numbers, so only long reads can tell them apart.
  They are ranked by the summed support of `--paths` reads with ≥ 2 novel adjacencies found in that
  order (cyclically for a circle), then by the number of matched paths.

---

## 8. Changes to original RCK: patches and fixes

> **Provenance:** see each item. **Code:** `rck/core/ilp_gurobi.py`, `rck/core/graph.py`,
> `rck/utils/adj/*.py`, `rck/rck_run.py`, `setup.py`.

**`rck/version-fix` (0b41976, Runpeng Luo), (b) patch.** RCK as published does not import on a
current scientific-Python stack. The commit changes no model, and reproduces upstream `master`'s
solution byte for byte on the same input:

- Gurobi 9 renamed the module `gurobi` to `gurobipy` and moved the four-argument
  `addConstr(lhs, sense, rhs, name)` to `addLConstr`. All 29 constraints in `ilp_gurobi.py` are
  linear, so `addLConstr` is an exact substitution;
- networkx 2.4 removed `connected_component_subgraphs` (five uses); a local helper restores it
  (floor raised to networkx ≥ 2.4);
- PyVCF is unmaintained; `setup.py` asks for `pyvcf3`, which provides the same `vcf` module;
- two invalid escape sequences in `rck_run.py` made raw strings.

**`b26bb96`, (c) fix to original RCK; changes results when the option is used.**

- `rck_run.py` read the `--scnb` copy-number boundaries from the `--scnt` path, so `--scnb` was
  silently ignored. It now reads `--scnb`.
- `MIPGapAbs` was set to the relative `--run-g-mip-gap` value (0.015 copies), which forced
  near-exact optimality. It is now set only from the new `--run-g-mip-gap-abs` (default: Gurobi's
  own).

**Packaging, part of (c).** `setup.py` adds the `rck.ecdna` package and the `rck-ecdna`,
`rck-ecdna-eval`, `rck-ecdna-plot` and `rck-ecdna-rerank` console scripts. `README.md` links this
document.

**Running on ionic, (b) patch** (`rck/ecdna/gurobi_env.py`, `41f5804`, `8f0f810`):
`default_threads` uses `--threads`, else `SLURM_CPUS_PER_TASK`, else the CPUs the process may run
on. `quiet_gurobi_logger` keeps gurobipy's Python logger off the console. `close_model` closes the
log file and disposes each model as soon as its values are read, because an open log on NFS leaves
a `.nfsXXXX` file. Test folders fall back when `/tmp` is full. `evaluation.json` writes missing
scores as `null` instead of `NaN`, which strict JSON forbids.

---

## 9. Evaluation

Scoring against simulator truth is documented in
[docs/pipeline/evaluation.md](pipeline/evaluation.md), with every metric defined:

- `rck-ecdna-eval` and `--truth-dir` (`evaluate.py`): truth reading, scope, per-pair scores, the
  recovery matrix and composite;
- the method-agnostic AA-format scorer `rck.ecdna.aa_eval`, with `aa_export`;
- the graph route (`aa_export subgraph`, `aa_import`, `coral_flank`).

Passing `--truth-dir` (and optionally `--truth-sample`, `--eval-tol-bp`, `--min-truth-fraction`,
`--all-truth-species`) to `rck-ecdna` writes `evaluation.tsv`, `evaluation.json`,
`truth_cycles.tsv` and `recovery_matrix.tsv`.

---

## 10. Known limitations

- **Nested species on one haplotype.** When allele data cannot separate species that share
  reference segments, copy numbers alone cannot either (§5.7). The step penalty decides by a prior,
  and the ambiguity report says so. Long reads that span the difference are the only evidence-based
  remedy. demo7's and demo8's three multi-junction reads each lie inside one species.
- **Allele observations are trusted as phased** in blocks stage 2 oriented, unless
  `--flip-oriented-penalty` is set. `--block-flip unoriented` remains opt-in (default `off`) and
  needs a subgraph with `segment_blocks.tsv`.
- **Species multiplicity is capped** at `--max-multiplicity` (3). bench_nested `ctrl_dup_s708` and
  `s710` traverse a segment 4 times and cannot be expressed as one species.
- **Light nested species** (bench `nest_samehap_6to1`) are split into several sub-threshold
  `low_weight` walks, so rank 1 has one significant species.
- **`low_weight` is flagged twice:** on integer weights for the ranking, and on refitted weights for
  the output. The two disagree in 8 of the 40 bench rank-1 answers, both ways.
- **Exact decomposition of a noisy integer flow** trades one copy of residual against one species.
  bench `s301` and `s703` have a near-correct answer in the pool that loses on the ranking. The
  noise-tolerant options did not change that (§5.3).
- **One sample per run.** A subgraph describes one sample. The simulator writes adjacency and path
  observations for the first tumour sample only (`scale.json` `junction_sample`).
- **Junction reads per copy fallback.** Without `scale.json`, `--extremities` or an override, θ is
  the median `#REF` of the subgraph's own reference adjacencies over $c_0$. Those adjacencies lie
  mostly inside the amplicon, so this can overestimate θ. The run logs a warning.
- **Scale errors upstream propagate.** A wrong ploidy class in `rck_inputs.py` doubles every segment
  copy number against the junction ones (bench `nest_cross_3sp_s403` at class 6). Rebuild the
  subgraph with `--rck-ploidy`.
- Long decompositions: three-species nested cases can take tens of minutes (`--time-limit` per MILP).

---

## Appendix A. Install

`rck-ecdna` needs Python 3.8 or newer, `numpy`, `pandas`, `networkx`, `scipy` (the refit) and
`gurobipy`.

```bash
cd RCK
pip install -e .          # rck, rck-ecdna, rck-ecdna-eval, rck-ecdna-plot, rck-ecdna-rerank
pip install gurobipy      # or: conda install -c gurobi gurobi
```

The figures also need `matplotlib` and `cnplot` [cnplot](pipeline/references.md#cnplot), in amplicon-allele-analysis
`third_party/cnplot` (Python 3.10 or newer). Both are optional: without them the run finishes and
logs that the plots were skipped.

```bash
pip install -e /path/to/amplicon-allele-analysis/third_party/cnplot   # also installs matplotlib and seaborn
```

Instead of installing cnplot you can point `--cnplot-path` (or `$CNPLOT_PATH`) at that checkout.

`pip install gurobipy` ships a size-limited licence (2,000 variables and 2,000 constraints). It is
enough for the unit tests but not for the demo decompositions, which need an academic licence
(Appendix F).

## Appendix B. How to run

From the root of amplicon-allele-analysis, with RCK installed, on the demo5 subgraph (one simple
cycle), which carries everything RCK needs:

```bash
rck-ecdna \
    --subgraph-dir simulations/runs/demo5_cycle1_ecdna/subgraph/amplicon_0 \
    --truth-dir simulations/runs/demo5_cycle1_ecdna \
    --outdir rck_out/demo5 --threads 4
```

Leave out `--truth-dir` on real data. Nested amplicons (recommended settings from
`results/bench_nested/`):

```bash
rck-ecdna --subgraph-dir SUB --outdir OUT --block-flip unoriented --step-penalty 0.2
# sweep λ offline from one pool run
rck-ecdna --subgraph-dir SUB --truth-dir RUN --outdir OUT --step-penalty 0.1 --keep-pool
python -m rck.ecdna.rerank OUT --lambdas 0.02 0.05 0.1 0.2 0.5    # one row per λ
```

The demo3 and demo4 fixtures live in RCK (`tests/ecdna/data/`); they have no windows, so give the
scales (the values `results/step_sweep/sweep.py` uses):

```bash
rck-ecdna --subgraph-dir tests/ecdna/data/amplicon_0 --truth-dir tests/ecdna/data/demo3_truth \
    --depth-per-copy 15 --junction-reads-per-copy 6 --outdir rck_out/demo3
```

**Windows mode.** With a simulator run's observations, `--obs-dir RUN/obs` sets `--windows-dir
RUN/obs/windows`, `--extremities RUN/obs/extremity.tsv` and `--paths RUN/obs/path.tsv`, and
`--inputs windows` forces the scales to be estimated from them (§3).

**Subgraph-only mode.** Subgraphs built by amplicon-allele-analysis `simulations/subgraph_construct.py`
(`yubi_dev`) carry, computed from the observations without truth (`simulations/rck_inputs.py`):

| File in the subgraph folder | Content |
|---|---|
| `scale.json` | depth per copy and junction reads per copy, the purity / tumour ploidy solution they come from, every ploidy class considered, `background_cn` (cell ploidy), the sample, and the sample the junction reads come from |
| `edges.tsv` columns `CN_TOTAL`, `CN_A`, `CN_B` | cell-averaged segment copy number and its allele split |
| `edges.tsv` column `CN` | observed adjacency copy number (`#REF` or `#NOV` / junction reads per copy) |
| `cn_profile.tsv` | copy number in 100 kb bins along each chromosome of the subgraph, for the region flanks |
| `paths.tsv` | long-read paths whose extremities are all in the subgraph |
| `phase_blocks.tsv` | one row per phase block of the subgraph's segments and junction flanks: `votes_same`, `votes_flip`, `orientation` (1 as observed, −1 flipped, 0 not oriented), `n_snps`, `status` (`voted` / `no_votes`). Older folders list only the voted blocks |
| `segment_blocks.tsv` | one row per (segment edge, phase block it overlaps): the parent edge's `NODE_1 NODE_2 CHROM_1 POS_1 POS_2` and `edge_index`, `phase_block`, `piece_start`, `piece_end` (half-open; the pieces tile the edge), `piece_length`, `n_snps`, `A_count_obs` / `B_count_obs` (observed phase), `orientation`, `A_count` / `B_count` / `BAF` and `CN_TOTAL` / `CN_A` / `CN_B` (oriented; `CN_A` / `CN_B` empty below 5 SNPs). Read with `--block-flip` |

`--inputs auto` (default) uses these files whenever `scale.json` exists. `--inputs windows` goes
back to estimating from `--windows-dir` / `--obs-dir`. `--depth-per-copy`,
`--junction-reads-per-copy` and `--baseline-cn` still override. A subgraph describes one sample
(`subgraph_construct.py --sample`, default the first tumour), so build one subgraph folder per sample.
`--sample` is ignored with a warning when it differs from `scale.json`. The depth fit in
`rck_inputs.py` (purity × ploidy grid, ploidy class chosen jointly over tumour samples, `--rck-ploidy`
and `--rck-depth-per-copy` overrides) is documented with stage 2.

## Appendix C. Inputs and options

| Option | Content |
|---|---|
| `--subgraph-dir` | `nodes.tsv`, `edges.tsv` (and optionally `phase_blocks.tsv`, `segment_blocks.tsv`) of one amplicon; with `scale.json` etc. it is the only input needed |
| `--inputs` | `auto` (subgraph files when `scale.json` exists), `subgraph`, or `windows` |
| `--obs-dir` | sets `--windows-dir OBS/windows`, `--extremities OBS/extremity.tsv`, `--paths OBS/path.tsv` |
| `--windows-dir` | UGP windows: `window.tsv.gz`, `window.depth.npz`, `sample_ids.tsv`, optional `snp.tsv.gz` + `snp.[AB]allele.npz` |
| `--extremities` | `extremity.tsv` with genome-wide `#REF`, for junction reads per copy |
| `--paths` | `path.tsv` long-read paths: read constraints (§5.4) and ordering ranking (§7) |
| `--hatchet-seg` | HATCHet `best.seg.ucn`, for depth per copy without a matched normal (cell lines such as COLO320) |
| `--sample`, `--outdir`, `--name` | sample name; output folder; name in output files (default: subgraph folder name) |
| `--truth-dir` / `--truth-cycles` | ground truth for the optional evaluation |

Options by group (defaults in brackets; section in parentheses):

- **Scale (§3):** `--depth-per-copy`, `--junction-reads-per-copy`, `--depth-mode tumor|ratio`
  [tumor], `--baseline-cn` [scale.json cell ploidy, else 2], `--block-bp` [1,000,000],
  `--min-snps-allele` [5], `--block-flip off|unoriented` [off], `--flip-oriented-penalty` [off],
  `--flip-tol` [0.01].
- **Graph (§2):** `--phantom-max-bp` [2], `--soft-junction-hap`.
- **Flow MILP (§4):** `--imputed-penalty` [2.0], `--novel-penalty` [0.05], `--alpha-novel` [0.5],
  `--alpha-ref` [0.25], `--bg-max` [4], `--bg-prior-weight` [0.1], `--soft-junction-penalty` [1.0],
  `--enforce-is`, `--no-background-novel`, `--linear-ends boundaries|all|none` [boundaries],
  `--linear-end-penalty` [0.25], `--boundary-min-copies` [2.0], `--boundary-frac` [0.25].
- **Decomposition (§5):** `--max-species` [4], `--max-multiplicity` [3], `--species-penalty` [1.0],
  `--step-penalty` [0], `--ambiguity-window` [2.0], `--keep-pool`, `--pool-lambdas` [0 0.02 0.05
  0.1 0.2 0.5], `--min-species-weight`, `--min-species-weight-frac` [0.10], `--min-ecdna-copies`
  [2], `--ecdna-call-fraction` [0.8], `--num-alternatives` [10], `--pool-gap` [0.10],
  `--alternatives-mode pool|nogood` [pool], `--adj-residual-weight` [1.0], `--seg-residual-scale`
  [1.0], `--residual-mode absolute|relative` [absolute], `--residual-sigma` [0.1],
  `--rank-residual first|sum` [first], `--max-orderings` [20].
- **Long reads (§5.4):** `--path-weight` [0.5], `--path-fraction` [0], `--min-path-support` [1],
  `--min-path-junctions` [2].
- **Refit (§6):** `--no-refit`, `--refit-min-copies` [0.5].
- **Evaluation:** `--truth-dir`, `--truth-cycles`, `--truth-sample`, `--eval-tol-bp` [10,000],
  `--min-truth-fraction` [0.5], `--all-truth-species`.
- **Figures:** `--no-plots`, `--plot-ranks` [1], `--plot-low-weight`, `--cnplot-path`.
- **Gurobi:** `--threads`, `--time-limit` (seconds per MILP), `--mip-gap` [1e-4], `-v`.

## Appendix D. Outputs

| File | Content |
|---|---|
| `decompositions.tsv` | one row per decomposition, species and ordering: `amplicon_call`, `circular_fraction`, `reads_explained`, `data_objective` (observation misfit after the refit), `topology`, `status` (`species` or `low_weight`), `weight` (copies per cell, bulk average, fractional after the refit), `weight_int` (integer weight of the structure search), `length_bp`, `path_support`, `uses_imputed`, `junctions`, `walk` (`chrom:start-end:hap:strand`, comma-separated) and `joins` (adjacency after each step, `kind:aid:phasing`; a linear walk's last join is `end:<extremity>:<hap>`) |
| `flow_segments.tsv` | per segment: observed total and allele copy numbers, inferred ecDNA and background copy numbers for A and B. With `--block-flip`, a segment with pieces reports its total split by the BAF of its piece counts summed after the chosen flips |
| `flow_adjacencies.tsv` | per adjacency labelling: read support, observed copy number, ecDNA and background flow; `linear_end` rows for copy-ends |
| `aa/<name>_graph.txt`, `aa/<name>_decomp<k>_cycles.txt` | AmpliconArchitect-format graph and cycles [Deshpande et al. 2019](pipeline/references.md#deshpande2019) (total copy number; cycle lines also carry `Haplotypes=`; linear walks start and end with AA's `0+`/`0-` source) |
| `observations.json` | scales, how each was estimated, and every observation; with `--block-flip`, also `pieces`: each piece's segment, block, span, totals, A / B in its base frame, `fixed`, `flipped` and the A / B after the flip (`a_final`, `b_final`) |
| `block_flips.tsv` | with `--block-flip` only: one row per phase block with a piece: `phase_block`, `orientation` (from the subgraph), `flip`, `status` (`fixed`, `flipped`, `kept` or `unresolved`), `n_pieces`, `n_snps`, `delta_if_reversed` |
| `summary.json` | graph normalisation log, flow status (with `linear_end_copies`), top-level `amplicon_call` and `circular_fraction`, decompositions (with `reads_total`, `reads_explained`, `data_objective`, `data_objective_int`, `pruned`: species the refit dropped, and `steps`), top-level `ambiguity` (§5.8); with `--block-flip`, `block_flip` (mode, number of pieces and blocks, counts per status, the flipped blocks) |
| `pool/pool.tsv`, `pool/pool.json`, `pool/decompositions.tsv`, `pool/evaluation.tsv` | with `--keep-pool` only: every pool solution before the λ-dependent pruning (`pool_index`, objective, residual, species, significant species and their steps, integer weights), the same solutions as decompositions (`decomposition_rank` = pool index + 1) and their evaluation; `pool.json` adds the settings, the rank-1 solution per `--pool-lambdas` value and those solutions' same-fit alternatives |
| `gurobi_flow.log`, `gurobi_decompose.log` | solver logs |
| `evaluation.tsv`, `evaluation.json`, `truth_cycles.tsv`, `recovery_matrix.tsv` | only with `--truth-dir` or `--truth-cycles` ([evaluation](pipeline/evaluation.md)) |
| (`python -m rck.ecdna.aa_export`) | not written by a run: the answer, its same-fit alternatives, the graph and allele copy numbers in the AA format the method-agnostic scorer reads ([evaluation](pipeline/evaluation.md)) |
| `plots/<name>_ecdna.pdf`, `plots/<name>_rank<r>_<species>.png`, `plots/<name>_recovery_matrix.png` | figures (Appendix E) |

## Appendix E. Figures

> **Provenance:** (c) new, `4e53747`; the topology track follows CoRAL's cycle figure
> [Zhu et al. 2024](pipeline/references.md#zhu2024); axes from cnplot [cnplot](pipeline/references.md#cnplot). **Code:** `rck/ecdna/plot.py`.

At the end of a run, `rck-ecdna` draws one page for each species of the top decomposition. Species
flagged `low_weight` are left out unless `--plot-low-weight` is passed, and `--plot-ranks N` also
draws the next N − 1 alternatives. Every track shares one cnplot `GenomeAxis` that covers only the
segments the species and its matched truth species use, padded by 5% (at least 5 kb), with the gaps
between them collapsed. The tracks are:

- **topology**: one row per step of the top ordering, numbered. Each step is a block at its reference
  position, blue for haplotype A and orange for B, with an arrow for the strand. Connectors follow
  the traversal order and are dotted for an imputed join and for the join that closes the circle.
  amplicon-allele-analysis draws the same track for simulated truth (`simulations/src/plots.py`), so
  the two can be compared side by side;
- **truth**: the matched truth species drawn the same way, when the run was evaluated;
- **junctions**: one arc per novel adjacency of the walk, coloured by the haplotype it joins. Grey
  marks a haplotype crossing (AB / BA) and a dotted arc an imputed one;
- **CN-A / CN-B**: allele-specific copy number, B mirrored under A. Points are the observations,
  drawn with cnplot `plot_scatter_1d` along each segment (with `--block-flip`, along each piece,
  after the chosen flips). The step line is the inferred ecDNA + background copy number. Shaded
  blocks are this species' own copies (weight × the number of times its walk uses the segment).

Beside the tracks, the walk is drawn as a ring: clockwise from the top in traversal order, one arc
per step coloured by haplotype, with a tick at every novel adjacency (dotted for an imputed one). A
linear species is drawn as an open ring. Arc lengths are to scale, except that every step gets at
least 1.5% of the ring, which the page notes when it applies. The title gives the decomposition's
call, the species' weight (refitted and integer), its length and, when the run was evaluated, the
matched truth species and its scores. After an evaluation, the last page is a heatmap of the
recovery matrix.

`rck-ecdna-plot --outdir <run>` redraws an existing output folder (`--plot-ranks`,
`--plot-low-weight`, `--cnplot-path`, `--no-png`), and `--no-plots` turns the figures off. A plotting
error is logged as a warning and never fails the run.

## Appendix F. Running on Princeton ionic

ionic schedules jobs with Slurm. `docs/examples/ecdna_ionic.sbatch` is a template;
`docs/examples/ecdna_demo3_tests.sbatch` runs the tests; `docs/examples/ecdna_bench.sbatch` and
`collect_bench.py` run and summarise a benchmark array.

- **Gurobi licence:** a named-user academic licence is tied to the machine where you ran `grbgetkey`,
  so it does not work on compute nodes. For batch jobs, request a Gurobi academic WLS licence and put
  its `gurobi.lic` (which holds `WLSACCESSID`, `WLSSECRET`, `LICENSEID`) in your home directory, or
  point `GRB_LICENSE_FILE` at it. If ionic provides a Gurobi module (`module avail gurobi`), its
  floating licence works as well. The WLS licence allows 2 concurrent sessions: keep array jobs to 2
  concurrent tasks (`--array=0-N%2`).
- **Environment:** the project environment is `/n/fs/yubi-iw/conda/envs/gurobi-env`. To create one:

  ```bash
  conda create -n rck-ecdna python=3.11 numpy pandas networkx scipy -y
  conda activate rck-ecdna
  pip install gurobipy
  pip install -e /path/to/RCK
  pip install -e /path/to/amplicon-allele-analysis/third_party/cnplot   # optional: figures
  ```

- **Batch jobs:** put the environment's `bin/` on `PATH`
  (`export PATH=/n/fs/yubi-iw/conda/envs/<env>/bin:$PATH`), as `results/step_sweep/sweep.sbatch`
  does. `source ~/.bashrc; conda activate` can fail on a compute node, where conda's shell setup is not
  loaded. Gurobi uses `--threads`, else `SLURM_CPUS_PER_TASK`, else the CPUs the process may run on.
- **Nested-ecDNA benchmark** (`results/bench_nested/`): simulations and subgraphs run as one array
  without Gurobi (amplicon-allele-analysis `simulations/bench_nested.sbatch`), RCK as a second
  (`results/bench_nested/bench_nested.sbatch`), and `collect.py` summarises; the order is in
  `results/bench_nested/README.md`.

## Appendix G. Tests

The tests are plain `unittest` modules under `tests/ecdna/`. Run them from the RCK root, with the
standard library or with pytest:

```bash
python -m unittest discover -s tests/ecdna -t . -v
# or
python -m pytest tests/ecdna -v
```

`discover` also runs the tests that solve MILPs (Gurobi licence). Without Gurobi, name the modules
marked "no" below, e.g. `python -m unittest tests.ecdna.test_aa_eval tests.ecdna.test_evaluate`.

Gurobi's log goes to `gurobi_*.log`; pass `-v` to also see it on the console. Test outputs go to
`$RCK_ECDNA_TEST_TMP` if set, else the system temp directory if it has at least 200 MB free, else
`~/.cache/rck-ecdna-tests`, and temp-directory cleanup errors (such as `.nfsXXXX` placeholders on NFS
homes) never fail a test. Tests that solve a MILP are skipped when `gurobipy` cannot be imported.
`test_demo3.py` also needs a full licence; setting `RCK_ECDNA_DEMO3_OBS` to a demo3 `obs/` folder
enables its second test, which estimates the scales from the windows.

| File | Needs Gurobi | Covers |
|---|---|---|
| `test_subgraph_io.py` | no | synthetic chains; phantom contraction and edge relocation; parallel imputed merging; junction labellings; the `amplicon_0` fixture (18 segments; 20 reference, 14 novel and 12 imputed adjacencies; 4 terminals; crossings AB/BA) |
| `test_scale.py` | no | depth mode from balanced blocks (rejecting LOH and CN 3 blocks); allele split; flank background; override and fallback; ratio mode cancelling shared bias; θ from extremities; sample selection |
| `test_walk.py` | no | alternating Hierholzer; canonical form under rotation and reversal; a tandem-duplicated segment; a hub with three loops (2 orderings); balance checks; path parsing and ranking; linear walks from end to end |
| `test_evaluate.py` | no | per-pair scores, truth reading and scoping ([evaluation](pipeline/evaluation.md)) |
| `test_ilp_circular.py` | yes (small) | single circle (10 copies on A plus a 1+1 background); haplotype crossings; imputed closing edge used when cheap, linear flow when expensive, no flow with `--linear-ends none`; linear ends only at boundaries and flanks; a linear amplicon; `--enforce-is` |
| `test_decompose.py` | yes (small) | repeated segment (one species, multiplicity 2); two species sharing a segment, with every optimal answer reported; disconnected circles; pool and no-good alternatives; orderings; a linear species; a circle plus a linear piece (called ecDNA); end-to-end runs for a circle (`ecDNA`) and an amplicon with no closing junction (`linear`, and `none` with `--linear-ends none`); read constraints (merging a read and its reverse, dropping single-junction, weak and off-flow reads); a read through a>x>b choosing the merged circle (reward and `--path-fraction`); a read never adding a species, and the order post-check; noise tolerance (weights and ranking key without Gurobi; a 2-copy dip on one loop) |
| `test_refit.py` | LP: no; CLI: yes (small) | a fractional weight (10.4); an unsupported leftover dropped; a supported one kept; the `--refit-min-copies` threshold; the last species never dropped; observation rows use the flow background; `--no-refit` |
| `test_subgraph_inputs.py` | observations: no; run: yes (small) | scales, segment and junction copy numbers and flank background from `scale.json`, `CN_*` columns and `cn_profile.tsv`; CN columns win over the scale; `--depth-per-copy` recomputes from `WEIGHT`; `--inputs auto/subgraph/windows`; a full run from the subgraph folder alone |
| `test_recovery.py`, `test_aa_eval.py`, `test_aa_import.py`, `test_coral_flank.py` | no | the scorers and the graph route ([evaluation](pipeline/evaluation.md)) |
| `test_plot.py` | no (figures: matplotlib and cnplot) | axis regions padded and merged; figures skipped when cnplot cannot be imported; a circular species with truth and recovery-matrix pages; a linear species without truth |
| `test_step_penalty.py` | yes (small) | a hub segment with two 3-segment loops: at 2:1 the merged walk wins without `--step-penalty` (ambiguous, decided by the species penalty) and the separate circles with it; at 3:1 tie-break against steps, gaps = λ × extra steps; default ranking unchanged; `--ambiguity-window 0`; pool ranking without Gurobi; `--keep-pool` re-ranks to the direct runs at λ = 0.02–0.5; the CLI refuses it at λ = 0; `rerank` reproduces rank 1 |
| `test_block_flip.py` | inputs: no; flow and run: yes (small) | pieces mapped through phantom contraction; orientation −1 swapped and fixed; a block missing from `phase_blocks.tsv` is unoriented; low-SNP pieces total-only; missing `segment_blocks.tsv` falls back with a warning; a swapped unoriented block is flipped (all ecDNA on A, no imputed edge); `--block-flip off` unchanged; an oriented block never flipped (unless `--flip-oriented-penalty`); a block without evidence `unresolved`; the refit uses flipped observations; the outputs |
| `test_gurobi_env.py` | no | thread default (explicit, `SLURM_CPUS_PER_TASK`, CPU affinity); gurobipy logger kept off the console |
| `test_demo4.py` | yes (full licence) | both demo4 subgraphs: evaluation never errors, each subgraph scored against its own species, per-sample truth weights (`RCK_ECDNA_DEMO4_OBS` enables estimated scales) |
| `test_demo7.py` | yes (full licence) | `--block-flip unoriented` on the rebuilt demo7 union subgraph (`RCK_ECDNA_DEMO7_SUBGRAPH`): chr4:185 and chr4:186 flipped, chr4:171, chr4:172 and chr4:187 fixed, no imputed edge |
| `test_demo3.py` | yes (full licence) | `amplicon_0` against the demo3 truth, rank 1: cycle order matches, haplotype Jaccard ≥ 0.95, junction recall ≥ 0.9, copy-number error ≤ 30%; the same with long-read paths (≥ 90% of read support explained) |

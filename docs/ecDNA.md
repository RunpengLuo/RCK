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
unchanged apart from the patches and fixes in §8; the multi-clone route (§11) adds a thin subclass
of it, for off-by-default additions, in a new file, `rck/core/ilp_gurobi_ecdna.py`.

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
11. [Multi-clone model (`--model multiclone`)](#11-multi-clone-model---model-multiclone)
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
| Alternating Euler walk | §7 | (a) Hierholzer's algorithm [Hierholzer & Wiener 1873](pipeline/references.md#hierholzer1873), adapted to alternation [Kotzig 1968] | `3c82929` |
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

In the below table, some properties were renamed to match the terminology in RCK (Aganezov & Raphael, 2020), but they are referred to as different names in the code, which is why there is a (code: ...) section.
| `edges.tsv` `EDGE_TYPE` | DIAG element | Observation carried |
|---|---|---|
| `reference` | segment (code: segment edge, reference edge) | `WEIGHT` = mean depth; `BAF`, `n_snps`; `CN_TOTAL`, `CN_A`, `CN_B` in subgraph mode |
| `extremity` | reference adjacency between the two extremities at a breakpoint (code: reference join) | `WEIGHT` = `#REF`, reads spanning the reference join |
| `variant` | observed novel adjacency (code: junction, variant edge) | `WEIGHT` = `#NOV`; `NOV_A`, `NOV_B`; `junction_hap` |
| `imputed_variant` | candidate novel adjacency proposed by stage 2's imputation (code: imputed edge) | none (its `WEIGHT` is kept as `imputed_support`, not as evidence) |

A reference adjacency touching an `imputed_ext*` extremity carries no read count. During subgraph construction, the algorithm finds amplified segments from read depth (section 2.4 of subgraph_construct.md). Each segment end is connected to an observed extremity within --tolerance bp. When no observed extremity is close enough, it creates an imputed extremity, named imputed_ext_<chr>_<n>, which is a breakpoint inferred from a change in depth, with no reads to back it. Its #REF is written as 0 only because the field needs a value. RCK should not use that 0, though, because if the flow model read the 0 as data, it would try to give the reference adjacency 0 copies, which cuts the genome incorrectly. So subgraph_io.py sets reads=None for any reference adjacency with an imputed_ext end. The model then has no read-based target there. Its copy number is set only by copy-number balance with the neighbouring segments and adjacencies.

### 2.2 Phantom segments

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::_contract_phantoms`.

The simulator names both sides of a breakpoint as separate extremities (`x+` and `x+1-`). Between
them sits a segment of 1 bp (`seg_hap = no_span`) whose two ends each have exactly one reference
adjacency, so it only carries reference flow. Every segment of at most `--phantom-max-bp` (default
2) bp with that shape is contracted. Its two reference adjacencies become one, whose read count is
their mean; its base joins the right-hand segment; and any other edge on its tail (head) moves to
the true tail (head) on the same side. `--phantom-max-bp 0` turns this off.

*Why:* This matches the graph input with the expected input from the original RCK algorithm. If the pantom segments aren't contracted, then every breakpoint adds a segment, so there becomes two extremities and a reference adjacency that the MILPs must balance which increases the size, leading to slower optimization. Also, a novel adjacency may attach to the phantom’s end rather than to the real segment it joins. The flow balance and the decomposition’s walks then route copies through a 1 bp piece, and segments don’t line up with the truth, which causes errors in deconvolution step.

### 2.3 Parallel adjacencies

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::_merge_parallel`.

Sometimes, there are edges between the same pair of extremities. This comes from subgraph construction’s imputation: its rules for proposing candidate novel adjacencies sometimes propose an edge that already exist.

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
This is borrowed from RCK telomeres. In RCK, telomeres are used to end paths, which is similar to terminals in this version.

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
   mode) and $\mathop{\mathrm{mode}}_{\log}$ is the mode of a histogram of $\log \bar d_B$ (bin width
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
*Why the mode:* Original RCK doesn't have a depth per copy (γ). Its input is segment copy numbers that another tool has already computed. RCK-ecDNA needs γ because its pipeline begins further upstream at raw observations. It has to convert mean read depth per segment, allele counts per segment, and spanning-read counts per novel adjacency into copy numbers.


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

*Why the mode:* Original RCK doesn't have a junction reads per copy (θ). Its input is segment copy numbers that another tool has already computed. RCK-ecDNA needs θ because its pipeline begins further upstream at raw observations. It has to convert these into copy numbers.

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
(c) new to the RCK model; read-based adjacency copy numbers follow AmpliconArchitect’s read model for breakpoint-graph edges [Deshpande et al. 2019], also used by CoRAL [Zhu et al. 2024]; θ estimation and the haplotype split are new.

$\hat c_j = r_j/\theta$ for every reference or novel adjacency with a read count $r_j$; imputed
candidates have none. This is the *observed adjacency copy number*. The `CN` column of
`edges.tsv` holds the same quantity. The flow's integer value (§4) is the *inferred* adjacency
copy number in RCK's sense.

### 3.5 Flank background

(c) new to the RCK model

For each region $[a,b]$ on chromosome $c$, the flank copy number is the bin-weighted mean of the
copy number over $[a-1\,\mathrm{Mb},\, a-10\,\mathrm{kb}]$ and $[b+10\,\mathrm{kb},\, b+1\,\mathrm{Mb}]$
(windows, or `cn_profile.tsv` in subgraph mode), else $c_0$. Every segment in the region gets it as
$\hat c^{\mathrm{bg}}_s$, the background prior of §4.5.

*Why it was added:* RCK models the whole genome. Every chromosome’s copies are in the model end to end, so the copy number under an amplicon is fixed by the same chromosome’s flow elsewhere. RCK-ecDNA works on a subgraph cut out of the genome, so it loses that context. The flanks stand in for the rest of the chromosome.

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

This is build off the original RCK algorithm, which has a a free flip per segment. RCK’s input is allele-specific, not haplotype-specific. So RCK decides, per segment, which allele belongs to which haplotype, as part of its MILP. In the RCK-ecDNA algorith, our observations are partly phased. Within a phase block, read-based phasing ties the segments’ A and B together, and only the block as a whole can be flipped relative to its neighbours. So, this step is RCK’s allele-to-haplotype assignment, constrained by the phasing information we have. The constraint makes the decision much smaller than RCK’s, with only one choice per unoriented block instead of one per segment.

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

Section 4.1 is for reference and does not describe the modified RCK-ecDNA pipeline.

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

### 4.2 Variables in RCK-ecDNA

| Symbol | Domain | Meaning | Code |
|---|---|---|---|
| $x^{e}_{s,h}$ | $\mathbb Z\cap[0,U_s]$ | ecDNA-layer copy number of segment $s$ on haplotype $h\in\{A,B\}$. This is similar to $c_{i,s,h}$ in the original RCK (integer copy number of segment $s$ on haplotype $h$ in clone $i$). The clone index is replaced by a layer index. Although RCK has one karyotype per clone, RCK-ecDNA has one bulk sample split into two components that add up: $x_{s,h} = x^e_{s,h} + x^b_{s,h}$ | `xe` |
| $x^{b}_{s,h}$ | $\mathbb Z\cap[0,\beta_{\max}]$ | background copy number. See above. | `xb` |
| $y^{e}_{j,\varphi}$ | $\mathbb Z\cap[0,M_j]$ | ecDNA-layer copy number of adjacency $j$ with labelling $\varphi\in\Phi(j)$. This is similar to $y_{i,j,\varphi}$: copy number of adjacency $j$ with labelling $\varphi$ in clone $i$. The clone index is replaced by a layer index: $y_{j,\varphi} = y^e_{j,\varphi} + y^b_{j,\varphi}$. | `ye` |
| $y^{b}_{j,\varphi}$ | $\mathbb Z\cap[0,\beta_{\max}]$ | background copy number; only for reference adjacencies and (unless `--no-background-novel`) observed novel ones, never imputed | `yb` |
| $p_{j,\varphi}$ | $\{0,1\}$ | labelling $\varphi$ of novel or imputed $j$ is used. This is the same variable and use case as RCK’s presence indicator per novel-adjacency labelling. New uses: it carries the --novel-penalty and --imputed-penalty costs, and it no longer feeds RCK’s at-least-$1-\mathrm{fp}$ usage constraint. | `pv` |
| $t_{v,h}$ | $\mathbb Z\cap[0,4\beta_{\max}]$ | background copies ending at terminal $v$. $t$ is similar to the telomeres used in RCK and the source and sink nodes in CoRAL and AmpliconArchitect. In RCK-ecDNA, t is allowed at subgraph terminals. | `t` |
| $e_{v,h}$ | $\mathbb Z\cap[0,U_{s(v)}]$ | ecDNA-layer copies ending at linear-end extremity $v\in\mathcal L$ (§4.4). This is based on AA’s source edges and CoRAL’s s–t paths, which allow linear amplicons to end instead of becoming circles. RCK-ecDNA adopts this same method here to let the ecDNA layer explain a linear amplicon, at a cost per copy-end (--linear-end-penalty) so circles stay the default. | `te` |
| $f_b$ | $\{0,1\}$ | phase block $b$'s A and B are swapped (§4.6). This is similar to $b_f$ in the original RCK, which was the allele flip of fragment $f$. In RCK-ecDNA, it is the same variable and has the same mechanism, ther eis one binary flip shared by every segment in a span whose A/B labelling is arbitrary. The span is a phase block instead of a copy-number fragment, and only blocks subgraph construction couldn’t orient get a flip. | `flip` |

$U_s=\lceil 1.5\max(\hat c_s,1)\rceil+3$ (`ub_factor`, `ub_add`). This is the same as RCK’s copy-number bounds per segment haplotype (--scnb). It has the same role in RCK-ecDNA, but this algorithm derives the bound automatically from the observed copy number instead of needing the --scnb flag. \\
$\beta_{\max}$ is `--bg-max` (4). This is a new flag that was added. Since the background is chromosomal, its copy number should stay small and is capped at 4 copies per haplotype. My reasoning for this was that the background layer holds the chromosomal copies under an amplicon, which is normally 1 per haplotype, but clones could gain copies, so the cap leaves room up to 4 per haplotype. I wanted to keep it small to stop the background from absorbing amplified copies that should be part of the ecDNA layer. I'm open to other suggestions for the value, but so far, tests run with this cap have successfully deconvolved the ecDNA species. \\
$M_j=\max$ of $U_s$ over the segments at $j$'s ends (or $\beta_{\max}$ if none), halved for a
self-loop (`_big_m`). This is based on the same big-M constraint from the original RCK algorithm.  \\
$x_{s,h}=x^e_{s,h}+x^b_{s,h}$ and $y_j=\sum_{\varphi}(y^e_{j,\varphi}+y^b_{j,\varphi})$ are the inferred segment and adjacency copy numbers. \\

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
end through their source edges [Deshpande et al. 2019](pipeline/references.md#deshpande2019) [Zhu et al. 2024](pipeline/references.md#zhu2024). Therefore, the RCK-ecDNA algorithm adds an addition to RCK that allows the ecDNA layer to start or stop at a set $\mathcal L$ of extremities chosen by `--linear-ends`:

- `boundaries` (default): extremity $v$ of segment $s$ is in $\mathcal L$ when it has no reference
  adjacency, or its reference neighbour has no segment (a region flank), or, for a neighbouring
  segment $s'$, $\lvert \hat c_s-\hat c_{s'}\rvert \ge \max(m, \phi\max(\hat c_s,\hat c_{s'}))$ with
  $m$ = `--boundary-min-copies` (2) and $\phi$ = `--boundary-frac` (0.25). That is where the end of a
  linear amplicon shows in the depth;
- `all`: every extremity of a segment;
- `none`: the strictly circular model.

Each copy-end costs $\kappa_{\mathrm{end}}$ = `--linear-end-penalty` (0.25), so a linear walk costs
0.5 per copy. These values are heuristic and I tested them on a bench of simulated data, but I'm open to tweaking these values. An imputed closing edge costs `--imputed-penalty` (2.0) once, so with 4 or more copies
an affordable imputed closure wins, while a long gap with no other support is explained as linear.

*Why:* without ends, an amplicon with no closing junction can only be explained by forcing a
circle, an imputed join, or deviations, which results in poor results. I first introduced this idea after poor results in an initial run on demo4, and after adding linear ends, was able to deconvolve demo4 correctly, it is in amplicon-allele-analysis `results/rck_ecdna_plan1/` ("Plan 1").

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

The segment data term is adopted from the original RCK algorithm. It was modified in the following ways: \\

1. RCK sums deviations over clones, but RCK-ecDNA compares the sum of its two layers (background and ecDNA), $x_{s,h} = x^e_{s,h} + x^b_{s,h}$, with one bulk observation. \\
2. RCK weights by $L_s = \lceil \ell_s/10\rceil$, which is linear in length with no cap. RCK-ecDNA uses $w_s = \min(\ell_s/50\text{ kb}, 1)$ with a floor of 0.02, so long background segments can’t outweigh short amplified ones. This was created in response to initial poor results on demo3 and demo4. \\
3. RCK-ecDNA is now scaling by observation size, $n_s = 1/\sqrt{\max(\hat c_s,1)}$. A deviation counts relative to the observed copy number, as with a Poisson standard deviation. On the other hand, RCK took copy numbers from a caller and didn’t model read noise. \\
4. A total-only case for segments without enough heterozygous SNPs. RCK didn't need this because it always received allele-specific copy numbers. \\

The adjacency read term, $\dfrac{\alpha_j}{\sqrt{\max(\hat c_j,1)}}\lvert y_j - \hat c_j\rvert$	is based on AmpliconArchitect’s read model for breakpoint-graph edges, which fits copy numbers to read counts on discordant and concordant edges (CoRAL keeps the idea, with long reads). AA maximises a Poisson likelihood, but RCK-ecDNA uses an absolute deviation scaled by $1/\sqrt{\hat c}$, a linear stand-in for Poisson noise that keeps the model a MILP. Weights $\alpha$ = 0.5 (novel) and 0.25 (reference) are hand-set and can be tweaked. \\
The linear copy-end cost, $\kappa_{\mathrm{end}}\sum e_{v,h}$, is based on AA’s source edges and CoRAL’s source/sink nodes, which let linear structures end. The per-copy cost and its value 0.25 are new (section 4.4). \\
The novel-adjacency cost, $\pi_j p_{j,\varphi}$ with $\pi_j$ = 0.05, replaces RCK’s constraint that at least $1-\mathrm{fp}$ (90%) of novel adjacencies be used. This is based on the parsimony principle and the value was hand-set. \\
The background prior, $\beta\sum_s w_s \lvert\sum_h x^b_{s,h} - \hat c^{\mathrm{bg}}_s\rvert$, pulls the background layer towards the flank copy number, so chromosomal copies stay out of the ecDNA layer. It stands in for the whole-genome context RCK has and that my subgraph does not have (section 3.5). $\beta = 0.1$, hand-set. \\
The balance tie-break, $\varepsilon_{\mathrm{bal}}\sum\lvert x^b_{s,A}-x^b_{s,B}\rvert$, stops background copies from swapping haplotype at no cost. Added together with the linear ends. \\
The imputed-candidate cost, $\pi_j$ = 2.0, makes a candidate adjacency without read evidence win only when it explains several copies (4 or more, against linear ends). encourages parsimony. value hand-set. \\
The Soft-label cost, $\pi_{\mathrm{soft}}$ (only with --soft-junction-hap), allows a labelling that the reads don’t support, at a cost, instead of forbidding it. Labels are defined in (subgraph_construct.md, section 3.3). This flag is off by default and honestly not used much, I only used it to try exploring other options for one of the demos, but I don't use it for solving for ecDNA and wouldn't recommend using it. \\
The flip cost, $\kappa_b f_b$ ($10^{-6}$), is a very small tie-break for a block with no evidence either way. When flipping and not flipping fit the data equally well, that cost makes pushes the labeling to not flip so that the block keeps its observed A/B labelling. I added this because for a lot of the unoriented blocks, the data gives no reason either way. The segments may have little allele imbalance, or no novel adjacency with a haplotype label nearby. Then $f_b = 0$ and $f_b = 1$ have exactly the same objective. Without the cost, the solver’s choice would be arbitrary and the block_flips.tsv report would be misleading because it would list blocks as “flipped” when nothing supported flipping them. \\



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


For --novel-penalty: The cost is added per label used: an indicator $p_{j,\varphi}$ for each novel adjacency $j$ and labelling $\varphi$, weighted by $\pi_j$ in the objective. There are two reasons for changing this strategy from the original RCK. First, not every novel adjacency belongs to the ecDNA. A union subgraph includes SVs on background chromosomes and in other clones. Forcing 90% of them into the flow like original RCK would push copies through joins that may not be on the amplicon. Secondly, RCK has no adjacency read counts, but RCK-ecDNA has an observed copy number for each novel adjacency ($\hat c_j = r_j/\theta$), fitted by the read term with weight $\alpha = 0.5$. So the data already says how many copies each one should carry. The penalty only has to stop the model using a novel adjacency the data doesn’t support. The value 0.05 for this flag is a heuristic and was hand-set to be a small number.

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

This table summarises §4.1–§4.6. The last column says why each change was made, and which other
method, if any, the change is based on.

| | RCK `OptModelMultiClone` | `rck-ecdna` `OptModelCircular` | Why, and what it is based on |
|---|---|---|---|
| copy numbers | per clone; integer | one bulk sample; integer, two layers (ecDNA + background) | A subgraph is one bulk sample, and the ecDNA layer must close into cycles while the chromosomal copies under the amplicon must not. The split into a circular part and a linear part follows AA's and CoRAL's separation of cycles from s–t paths (§4.2) |
| telomeres | allowed at listed telomere extremities ($\ge$ balance) | ecDNA layer: only at linear ends, each copy-end costs `--linear-end-penalty`; background: at terminals | RCK's telomere relaxation, moved to the edges of a cut-out subgraph (terminals, §2.4) for the background, and made a per-copy cost for the ecDNA layer after AA's source edges and CoRAL's source and sink (§4.4) |
| novel adjacency use | at least $1-\mathrm{fp}$ of them used | no requirement; each used costs `--novel-penalty` | Not every novel adjacency in a union subgraph is on the ecDNA, and read-based adjacency copy numbers now decide use (§4.5) |
| missing adjacencies | n/a | imputed candidates, each costing `--imputed-penalty`; background may not use them | Stage 2 proposes candidate joins where no novel adjacency was observed; each must explain several copies to be used (§4.5) |
| allele flips | per input fragment, free | per unoriented phase block, $10^{-6}$; oriented blocks fixed unless `--flip-oriented-penalty` | RCK's fragment flip $b_f$, applied to phase blocks the subgraph could not orient; the tiny cost keeps a block with no evidence as observed (§4.6) |
| objective | length-weighted L1 to input allele copy numbers | length-capped, $1/\sqrt{c}$-normalised L1 to allele-specific (or total) segment observations, plus adjacency read term, background prior and penalties | RCK's length-weighted L1, adapted to raw observations; the adjacency read term follows AA's read model for breakpoint-graph edges; the other terms are new and heuristic (§4.5) |
| bounds | input `--scnb`; big-M from them | $U_s$ from the observations; same big-M rule; background capped at $\beta_{\max}$ | RCK's bounds and big-M, derived automatically; the background cap is new (§4.2) |
| infinite sites | on | optional (`--enforce-is`) | RCK's extremity exclusivity; off by default because chromothripsis-derived ecDNA violates it (§4.3) |

---

### 4.8 Depth bias (`--depth-bias-window`, opt-in)

> **Provenance:** (c) new, 2026-11-12. **Code:** `rck/ecdna/ilp_circular.py::depth_bias_windows`,
> `::OptModelCircular._beta`, `::_biased`, `::_flip_beta`, `::_bias_terms`; `rck/ecdna/refit.py::observation_rows`.
> **Evidence:** `results/bench_nested/NEXT_STEPS.md`, "Flow-stage noise, step 0".

Read depth carries a smooth multiplicative bias (GC content, mappability) that nearby segments
share. The integer flow follows it, and a few-copy dip that no junction supports then has to be
covered by an extra species in the decomposition. On bench_nested, removing the simulator's
correlated bias (σ = 0) fixed 5 of 8 failing cases, and the segment depth error was correlated
between neighbouring segments (lag-1 r = 0.60) while the junction-read error was not affected.

With `--depth-bias-window W` > 0 (default 0, off):

- segments are grouped, per chromosome in coordinate order, into windows of at most about W bp
  that start at a segment; a segment longer than W, and a gap over W, start a new window.
  Consecutive windows not separated by such a gap are neighbours $\mathcal N$;
- each window $k$ gets a continuous $\beta_k\in[-\beta_{\max},\beta_{\max}]$ (`--depth-bias-max`, 0.3);
- every segment and piece target $\tau$ of §4.5–§4.6 in window $k$ becomes $\tau\,(1-\beta_k)$.
  This is linear because $\tau$ is data. A flippable piece's target
  $(a_\pi+f_b(b_\pi-a_\pi))(1-\beta_k)$ uses $z_{b,k}=f_b\beta_k$, linearised exactly
  (McCormick, binary $f_b$, bounded $\beta_k$);
- the objective adds

  ```math
  \kappa_\beta\sum_k|\beta_k|+\kappa_{\mathrm{sm}}\sum_{(k,j)\in\mathcal N}|\beta_k-\beta_j|
  ```

  with $\kappa_\beta$ = `--depth-bias-cost` (1.0; a 10% shift of one window costs 0.1) and
  $\kappa_{\mathrm{sm}}$ = `--depth-bias-smooth` (1.0).

Junction reads, the background prior and the linear-end boundaries are unchanged. The refit (§6) and
the phase-block report use the corrected targets. `flow_segments.tsv` gains `depth_bias` and
`obs_corrected`, and `summary.json` gains `depth_bias`, only when the option is on. With the option
off the model and every output are as before.

*Choice of costs.* The offline pricing in step 0 found that the bias removes most of the noisy
flow's advantage only when it is cheap (κ_β ≈ 1–2). At 5 or more it does almost nothing. The risk
of a cheap bias is that a true copy-number change without a junction inside the window is absorbed.
The bound, the smoothness term and the bench guard cases are what check this. It is not yet
validated by a solve on the bench.

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
| $\zeta_{keb}$ | $[0,W]$ | $=w_k\beta_{keb}$ (linearised), this is new to RCK-ecDNA and added for linearization |
| $y_{ke}$ | $\{0,1\}$ | $e$ is in slot $k$'s support |
| $\lambda^{\mathrm{lin}}_k$ | $\{0,1\}$ | slot $k$ is a linear walk. this is the same as CoRAL’s s–t paths. |
| $r_e$ | $\ge 0$ | residual on edge $e$. this replaces CoRAL’s coverage constraints. |

Constraints, for every slot $k$:

```math
u_k \le w_k \le W u_k,\qquad
\zeta_{keb}\le W\beta_{keb},\quad \zeta_{keb}\le w_k,\quad \zeta_{keb}\ge w_k-W(1-\beta_{keb})
```
Left most constraint means that a used slot has a positive weight and an unused one has zero. This is from CoRAL.
The other 3 constraints write the multiplicity in binary so each product $w_k\beta_{keb}$ of an integer and a binary can be expressed exactly with linear constraints, technique published in https://pubsonline.informs.org/doi/10.1287/opre.23.3.434. CoRAL keeps $\sum_i w_i x_{ei}$ as a quadratic (bilinear) term, which Gurobi has to solve as a nonconvex problem. I changed it in RCK-ecDNA to use linearising to keep the decomposition a plain MILP, which is usually faster and more reliable to solve to optimality. The cost is integer weights.

```math
m_{ke}\le \bar m\, y_{ke},\qquad y_{ke}\le m_{ke},\qquad y_{ke}\le u_k,\qquad \sum_e y_{ke}\ge 2u_k
```
These constraints mark which edges a species uses, and require a used species to have at least two (one segment and one adjacency). CoRAL has the same inside its connectivity constraints.

```math
m_{k,\mathrm{seg}(x)} = \sum_{e\ \text{adjacency-like at}\ x} m_{ke}\quad\text{for every phased extremity } x,
\qquad \sum_{e\ \text{linear end}} m_{ke}=2\lambda^{\mathrm{lin}}_k,\quad \lambda^{\mathrm{lin}}_k\le u_k
```
At each phased extremity $x$, the multiplicity of its segment is the summed multiplicity of its adjacencies and linear ends (a self-loop counted twice). This is the same as CoRAL’s per-cycle Eulerian constraint, which makes each species an alternating walk (segment, adjacency, segment, …). But in RCK-ecDNA, I apply it at phased extremities, so a walk stays on one haplotype through a segment, and can switch only through an AB or BA labelling. \\
A species has either no ends (a circle) or two (one linear walk). This is like the source and sink edges ($s$, $t$) and an s–t path for linear structure in CoRAL. RCK-ecDNA routes all linear ends through one virtual node and counts ends instead. This is the decomposition side of section 4.4’s linear ends, so the flow and the decomposition use the same mechanism.

```math
w_k\ge w_{k+1},\qquad u_k\ge u_{k+1}
```
Symmetry breaking for MILP \\

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
[Glover 1975](pipeline/references.md#glover1975). The ordering constraints break the symmetry between slots. \\

**Explanation for Residual** \\
For CoRAL, the residual is a one-sided fit, $\sum_i w_i x_{ei} \le \mathrm{CN}_e$ where cycles can not exceed an edge’s copy number. It requires the cycles to cover at least 90% of the graph’s length-weighted copy number (p_total_weight), and at least 90% of every discordant edge’s (p_bp_cn). In RCK-ecDNA, it penalises the deviation in both directions, weighted per edge ($\omega_e$, section 5.3). This is because the input is a fitted, balanced integer flow, not noisy raw copy numbers, so it can in principle be reproduced exactly. A two-sided residual aims for that, and measures how far an answer is from it. The 90% thresholds were fixed choices by CoRAL and a hard cap can make a decomposition infeasible when copy numbers are uneven. I chose to use a penalty so that we can always get an answer and lets the objective weigh misfit against parsimony. Also, in RCK-ecDNA, the decomposition MILP gives several alternative decompositions, and RCK-ecDNA then sorts them to choose rank 1 using the residual. \\
**Explanation for Objective** \\
The weighted misfit is replaces CoRAL’s covered-fraction term and its 90% constraints. The species penalty is from CoRAL. CoRAL minimises the number of cycles minus the covered fraction, here, I just put both into one weighted sum. The tie-break is new and was added to prefer walks with fewer traversals among equally good answers. The long-read paths are also from CoRAL. \\


**Differences from CoRAL's `minimize_cycles`** (`coral/cycle_decomposition.py::minimize_cycles`,
`coral/models/concrete.py::get_minimize_objective`, `coral/models/constraints.py`, CoRAL 3.1.1):

| | CoRAL 3.1.1 | `rck-ecdna` |
|---|---|---|
| graph | decomposes the copy numbers on AA’s breakpoint graph, with total copy number per edge | decomposes a haplotype-resolved flow, so each species comes out with a haplotype on every segment and adjacency. decomposes the fitted integer flow rather than raw copy numbers because the flow MILP has already reconciled the observations, so the decomposition only has to split a balanced flow |
| weights $w$ | continuous ($\ge 0$) | integer in the MILP in $[0,W]$; fractional after the refit (§6) |
| multiplicity $x$ | integer, 0–10 | integer, times species $k$ traverses edge $e$, 0–`--max-multiplicity` (3). 3 was set at the start of setting up the simulations because it is the largest cap at 2 bits, the smallest that allows a segment traversed more than twice, like in a tandem duplication inside a circle. CoRAL allows 10 (4 bits), but CoRAL solves on a smaller, unphased graph. |
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

**Where it comes from.** Requiring a species to be one connected walk is needed in every cycle
decomposition, and CoRAL has it too. The method here is different from CoRAL's:

- CoRAL writes connectivity as spanning-tree / node-order constraints, added to the model before
  it is solved (`coral/models/constraints.py`).
- `rck-ecdna` uses **subtour-elimination cuts**, the standard way to forbid disconnected pieces
  in the travelling-salesman problem [Dantzig et al. 1954](pipeline/references.md#dantzig1954). The cut above is
  the cut-set form for walks that need not visit every node: if a slot uses an edge in $C_i$ and
  one in $C_j$, at least one used edge must leave $C_i$.
- The cuts are **lazy**: Gurobi calls back at each integer solution, and a cut is added only for
  the components that solution actually has [Gurobi](pipeline/references.md#gurobi). This keeps
  the model small, because most of the exponentially many possible cuts are never needed. Writing
  connectivity up front, as CoRAL does, adds variables and constraints for every slot and node
  whether or not they bind.

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
outputs are unchanged. \\

Weighting misfit by segment length idea comes from RCK and CoRAL. RCK weights each segment’s deviation by $L_s = \lceil \ell_s/10\rceil$ (section 4.1). CoRAL measures fit by length-weighted copy number: its objective and its 90% coverage target count each segment’s copies times its length. \\
The 50 kb cap and the 0.05 floor are added for RCK-ecDNA. This is the same cap the flow MILP uses ($w_s$, section 4.5), so long background segments don’t outweigh short amplified ones. \\
The adjacency weight is similar to CoRAL’s coverage target of at least 90% of every discordant edge’s copy number (p_bp_cn). In RCK-ecDNA, the fits on breakpoints is judged per edge, not by length, because a novel adjacency has no length. Also, I use a soft weight of 1 instead of a hard 90% threshold, and the option --adj-residual-weight. \\
--seg-residual-scale and --adj-residual-weight are used to scale the segment or adjacency weights down. These were added for tuning. \\
--residual-mode relative measures misfit in units of the expected noise. A deviation of 10% of an edge’s copy number costs about one unit, at any copy number. This come from the idea that weighting residuals by their expected standard deviation, like weighted least squares or a χ²-type fit. It is based on an assumption of proportional noise (10% of the copy number), which is usually what is produced by multiplicative biases like GC content. \\
--rank-residual sum ranks on residual plus structure cost.

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

**Where each part comes from.**

- *Using multi-junction long reads as constraints on the decomposition* is CoRAL's idea: a read
  that spans several novel adjacencies shows they lie on one molecule. Its optional threshold,
  `--path-fraction`, plays the role of CoRAL's `p_subpaths` (0.9).
- *Writing "slot $k$ contains read $p$" with binaries* ($o$, $z$, $x$ as ORs and ANDs of other
  binaries) is a standard way to express logical conditions as linear constraints.
- *What a read must match* was loosened. CoRAL ties each path edge to the multiplicity of the
  cycle that explains it. `rck-ecdna` only requires the read's novel adjacencies to be in one
  species, and checks their order afterwards against the walks (§7). The decomposition MILP has no
  order, only multiplicities, so order cannot be checked inside it. Adjacency containment is the
  part that can, and order is checked where it exists.
- *Weighting reads by support* (summed over identical or reversed paths) is new. A path seen in 10
  reads is stronger evidence than one seen once, and CoRAL counts paths, not reads.
- *A soft reward, capped below one species* ($\pi_{\mathrm{path}}$ = 0.5 against $\rho$ = 1) is
  new. Reads then choose between decompositions that fit the flow equally well, which is where
  they are most useful (§5.7), but cannot add a species on their own. The cap is a design choice,
  not a tuned value. CoRAL instead makes path constraints hard in its exact model.
- *Fixed by the slot supports:* $z$ and $x$ are fully determined by the $y_{ke}$, so they do not
  create extra pool solutions that differ only in these binaries (`c5e0b9c`).

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

**Where each part comes from.**

- *The solution pool* is a Gurobi feature: with `PoolSearchMode=2` Gurobi searches for the best
  `PoolSolutions` solutions within `PoolGap` of the optimum [Gurobi](pipeline/references.md#gurobi).
- *The no-good cut* is the standard cut that removes one binary solution and nothing else: at least
  one support indicator must change. It is the fallback when the pool is not wanted.
- *Reporting ranked alternatives, rather than one answer,* is this pipeline's choice. Nested and
  repeated structures often have several decompositions that fit equally well (§5.7, §5.8), and a
  single answer would hide that.
- *Re-checking and deduplicating* is new. Pool solutions are integer-feasible but can differ only
  in slot order or in a repeated circle written two ways, so they are validated and merged before
  ranking (`3c82929`, `f84faff`).

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

**Where each part comes from.** Both rules are new, and both came from the first runs on ionic
(`41f5804`).

- *Repeated circle.* A slot that goes round one circle $g$ times reproduces the same flow as one
  circle at $g$ times the weight, and it is the same molecule. Before this rule, a species at
  weight 5 traversing its circle twice was reported as its own species, instead of one at weight
  10. Dividing the multiplicities by their greatest common divisor gives the simplest description.
  This is plain arithmetic, not taken from another method.
- *`low_weight`.* The integer flow is uneven (§6), and reproducing it exactly needs small leftover
  species. On demo3 the true species came with two weight-1 multi-megabase circles that were only
  residue. Flagging, rather than deleting, keeps them visible in the output while keeping them out
  of the ranking and the scores.
- *The 10% threshold* is heuristic, chosen so that a weight-1 leftover next to a 10-copy molecule is
  flagged. The Task 1 scorer uses the same 10% rule as its default *strict* cut for AA's and CoRAL's
  cycles (§9), so every method is judged by the same rule. It is a known weak point for light
  nested species (§10).

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

**Where each part comes from.** Everything in this section is new to this pipeline.

- *Residual first.* A decomposition that fits the flow better beats a simpler one that fits worse.
  This needs the two-sided residual of §5.1, which gives every pool solution a number. CoRAL's
  one-sided fit with 90% coverage thresholds is pass or fail, so it cannot rank a pool.
- *Then structure cost, then species count.* Among equally good fits, the simpler explanation wins.
  The species count $\rho\, n_{\mathrm{sig}}$ is the same parsimony as CoRAL's cycle count (§5.1),
  but only significant species count, so weight-1 leftovers (§5.6) do not make the true answer lose
  (demo3, `f84faff`).
- *The step penalty $\lambda S$* is a minimum-description-length prior on walks: among explanations
  that fit equally well, prefer the one that is shorter to write down [MacKay 2003](pipeline/references.md#mackay2003).
  A walk is described by its steps, so a merged walk that re-traverses shared sequence costs more.
  CoRAL's objective has no term like it.
- *The pruning rule* keeps only solutions the prior could still prefer, so the widened pool does not
  slow the ranking.
- *λ = 0.2 recommended* comes from the 40-case bench (`results/bench_nested/`), as the smallest λ
  on the plateau of best mean composite with no control or demo case losing. The default stays 0, so
  earlier results do not change.

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

**Where each part comes from.** The report is new; no other method here reports it.

- *Fixing the covered flow* ($\sum_k w_k m_{ke}$ equal to the top decomposition's on every edge)
  makes the fit identical by construction. Any alternative found can only differ in structure, which
  is exactly the question: is there another set of molecules the data cannot tell apart?
- *The distinctness constraint* forbids two slots with the same support, so a walk split into
  identical copies is not counted as a different answer.
- *`decided_by` and the gap* say which part of the ranking (§5.7) chose rank 1, and by how much.
  This makes visible when rank 1 rests on a prior or a tie-break rather than on the data, which is
  what most of the bench's remaining failures turned out to be.
- *The window* (`--ambiguity-window`, 2.0, about two species) is heuristic.

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

**Where it comes from.** New, and an engineering step rather than a modelling one. It relies on one
property of the design in §5.7: because λ acts only in the ranking, the MILP and its pool are the
same for every λ > 0, so a saved pool can be re-ranked for any λ without solving again. Putting the
penalty in the MILP, as the first design did, would have made this impossible. It made the λ sweep
on the bench affordable under a 2-session Gurobi licence.

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

**Where each part comes from.**

- *Classifying an amplicon from its reconstruction* follows AmpliconArchitect
  [Deshpande et al. 2019](pipeline/references.md#deshpande2019), whose amplicons are classified as ecDNA,
  breakage–fusion–bridge, complex non-cyclic or linear by its companion tool AmpliconClassifier
  `[citation needed: AmpliconClassifier paper and repository, not verified]`, and CoRAL, which
  reports cycles and linear paths [Zhu et al. 2024](pipeline/references.md#zhu2024). AmpliconClassifier's rules are not
  reproduced here.
- *The rule here is new:* the share of copy-weighted length ($w_k L_k$) that lies on circles. Weighting
  by copies and length measures how much of the amplified DNA is circular, so a long, high-copy
  circle outweighs a short linear leftover.
- *Circles below $c_{\min}$ = 2 copies do not count* towards the circular share, so a weight-1
  leftover circle cannot turn a linear amplicon into an ecDNA call.
- *The 80% threshold* and the 2-copy minimum are heuristic. They were added with the linear ends
  (`10416aa`), when the model first had to choose between circular and linear explanations.

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

**Where each part comes from.**

- *Fractional species weights* follow CoRAL, whose cycle weights are continuous
  [CoRAL 3.1.1 code](pipeline/references.md#coral_code). Real ecDNA copy numbers are cell averages, not integers. The
  decomposition MILP uses integer weights only because it splits an integer flow and because of the
  exact linearisation (§5.1); the refit restores continuous weights once the structures are chosen.
- *Fitting to the observations, not to the flow,* is new. The integer flow is already a rounded fit
  of the observations, so fitting weights to it would carry its unevenness into the copy numbers.
  The rows are the flow MILP's own data terms (§4.5), so "fit" means the same thing in both
  stages.
- *Least absolute deviations* (the weighted L1 misfit) is the same misfit as the flow MILP, and is
  less pulled by a single outlying segment than least squares. Written as an LP, it is solved by
  HiGHS [Huangfu & Hall 2018](pipeline/references.md#huangfu2018) through SciPy [Virtanen et al. 2020](pipeline/references.md#virtanen2020).
- *The small weight penalty* $\eta = 10^{-3}$ per copy prefers fewer copies among equal fits
  (`weight_penalty`), the same tie-break device as elsewhere in the model.
- *Dropping species refitted below 0.5 copies* is new. It removes leftovers that only absorbed flow
  unevenness, so they no longer take copies from the real molecule: on demo4 tumor2 the integer
  answer was 6 + 1 copies against a true 7.44 (`c5e0b9c`). The 0.5 threshold is heuristic: half a
  copy per cell.
- *The structures and the ranking are not changed,* so the refit only adjusts copy numbers; it
  cannot turn a wrong structure into a right one.

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

**Where each part comes from.**

- *Turning a balanced, connected multiset into a closed walk* is Hierholzer's algorithm for Euler
  circuits [Hierholzer & Wiener 1873](pipeline/references.md#hierholzer1873): build a closed trail, then splice in
  sub-trails until every edge is used.
- *Alternation* is the adaptation needed for a breakpoint graph: segments and adjacencies must take
  turns, so sub-trails are spliced in only before a segment step. That such an alternating circuit
  always exists when every extremity is balanced is a known result for Euler trails in graphs with
  two edge colours `[citation needed: Kotzig 1968 proposed, not verified]`. AA's and CoRAL's cycles
  are alternating walks in the same sense.
- *Closing a linear species with a pseudo segment* reduces a linear walk to the circular case, so
  one algorithm serves both. It was added with the linear ends (`10416aa`).
- *Enumerating several orderings* is new. One multiset can often be walked in more than one order
  (for example when a segment is visited twice), and the orders imply different molecules. A bounded
  depth-first search lists them, up to rotation and reversal, which describe the same molecule.
- *Ranking orderings by long reads* is new, and uses the same reads as §5.4. Copy numbers cannot
  separate orderings, so reads that span several novel adjacencies in a specific order are the only
  evidence. CoRAL uses the same kind of evidence inside its decomposition, with path constraints on
  cycle edges; `rck-ecdna` checks order here, where a walk with an order exists.

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

**`yubi_CNP_ionic` (2026-11-13), (c) fixes to original RCK.**

- `60ca24d`: after `b26bb96`, `rck_run.py` still failed on any `--scnb` file:
  `SegmentCopyNumberBoundaries.update()` iterates its argument as a dict, but was given the
  `SegmentCopyNumberBoundaries` object (`TypeError`). It now passes the object's records. Checked on
  the demo5 DIAG with the bounds only in `--scnb`: they are applied, and the solution is unchanged.
- `18dd969`: `rck/utils/adj/convert.py` imported PyVCF at module level, so `rck_run.py`, which
  never reads a VCF, failed in an environment without it (the ionic gurobi env). The import now
  happens inside the two VCF readers; the VCF converters still need PyVCF.

**Packaging, part of (c).** `setup.py` adds the `rck.ecdna` package and the `rck-ecdna`,
`rck-ecdna-eval`, `rck-ecdna-plot` and `rck-ecdna-rerank` console scripts. `README.md` links this
document.

**Running on ionic, (b) patch** (`rck/ecdna/gurobi_env.py`, `41f5804`, `8f0f810`):
`default_threads` uses `--threads`, else `SLURM_CPUS_PER_TASK`, else the CPUs the process may run
on. `quiet_gurobi_logger` keeps gurobipy's Python logger off the console. `close_model` closes the
log file and disposes each model as soon as its values are read, because an open log on NFS leaves
a `.nfsXXXX` file. Test folders fall back when `/tmp` is full. `evaluation.json` writes missing
scores as `null` instead of `NaN`, which strict JSON forbids.

**Why these matter for the comparison.** The `rck/version-fix` patches and the `b26bb96` fixes touch
only original RCK (`rck`, `rck/core/*`), not `rck-ecdna`. They matter because original RCK is also
a baseline in the benchmark, run on the same graphs as an allele-specific karyotype method with no
ecDNA species (results/graph_route/NEXT_STEPS.md). The patches let it run on current software
without changing its answers; the `--scnb` fix changes answers only when that option is used, and
should be stated whenever original RCK is run with it.

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

**Which metrics come from which method.** The definitions, equations and sources are in
[evaluation.md](pipeline/evaluation.md); in short:

| Metric | Level | From |
|---|---|---|
| Segment and haplotype Jaccard, phasing accuracy, junction F1, order match, copy-number error, composite | species | (c) new to this pipeline, built on the Jaccard index [Jaccard 1912](pipeline/references.md#jaccard1912) |
| composite-H0 (no phasing) | species | (c) new, for comparisons with methods that do not phase |
| Error-free reconstruction | graph | (a) AmpliconArchitect's simulation benchmark [Deshpande et al. 2019](pipeline/references.md#deshpande2019); `error_free_merged` is (c) |
| Breakpoint-graph accuracy, cycle-interval overlap, cyclic LCS, length error | graph, species | (a) CoRAL [Zhu et al. 2024](pipeline/references.md#zhu2024); both the code and the supplementary-methods versions are reported |
| Novel-adjacency FNR and FDR; allele-specific copy-number error; infinite-sites violations | graph | (a) RCK [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020); FDR instead of FPR is (c) |
| Strict and best-of modes | all | best-of follows CoRAL's paper; strict uses the 10% rule of §5.6 for every method |

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
- **The answer can depend on which equally good flow the solver returns.** The flow MILP often has
  several integer flows with the same objective, and Gurobi's choice depends on model size and
  variable order. In the graph-route round trip, removing candidate adjacencies the flow never used
  left the flow objective unchanged but changed the decomposition in 29 of 36 cases
  (`results/graph_route/`). Results are therefore single-seed; repeats with several Gurobi seeds are
  planned, and preferring the most decomposable flow among equal ones is a candidate fix.
- **Most weights are heuristic and untuned.** The objective weights of §4.5 and §5.1 (for example
  `--novel-penalty` 0.05, `--linear-end-penalty` 0.25, `--imputed-penalty` 2.0, `--bg-max` 4,
  `--species-penalty` 1.0, `--max-multiplicity` 3) were set by hand, most of them on demo3 and
  demo4, and kept fixed. Only
  λ (`--step-penalty`) has been swept on the bench. Linear ends (`--linear-ends boundaries`) are on in
  every run and have not been compared against `none`.

---

## 11. Multi-clone model (`--model multiclone`)

*Added on `yubi_CNP_ionic` (2026-11-13), step 3 of [`multiclone_plan.md`](multiclone_plan.md);
status and next steps in [`multiclone_NEXT_STEPS.md`](multiclone_NEXT_STEPS.md). Off unless
`--model multiclone` is given; the bulk route (§2–§7) is unchanged.*

The bulk route reinvents parts of RCK because a subgraph is one bulk sample. The multi-clone route
instead runs RCK's own model [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020). Its
only data term is per-clone allele-specific segment copy numbers from an upstream caller (HATCHet;
on the simulator, the truth as a stand-in), and it infers per-clone integer karyotypes. Which
copies form ecDNA circles and which form linear chromosomes is decided afterwards, by decomposing
the karyotypes heaviest walk first, as AmpliconArchitect [Deshpande et al. 2019](pipeline/references.md#deshpande2019)
and CoRAL [Zhu et al. 2024](pipeline/references.md#zhu2024) do (§11.2, §11.3).

**Input: the annotated DIAG.** `--rck-input <subgraph>/rck_input/` is written by
amplicon-allele-analysis `subgraph_construct.py --clone-profiles`
([its §5.6](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_CNP/simulations/subgraph_construct.md)),
in RCK's native formats:
- `rck.scnt.tsv`: per-clone copy numbers in the observed phase frame, with the bounds as `cnb=`;
- `rck.adj.tsv`: observed novel adjacencies;
- `rck.fragments.tsv`: one fragment per phase block;
- `rck.telomeres.tsv`: the outer flank ends;
- `rck.adjgroups.tsv`: molecule groups from long reads;
- `rck.scnb.tsv`: the bounds again, for `--scnb`;
- two files RCK does not read: `rck.adj.imputed.tsv` and `clone_proportions.tsv`.

`rck/ecdna/multiclone.py` reads them exactly as `rck/rck_run.py --no-pre` does, with RCK's readers,
alignment checks and run defaults (`--run-nas-fp 0.1`, `--run-group-m-default-fp 0.1`,
`--run-segment-length-attr length_10`).

**Model** (`rck/core/ilp_gurobi_ecdna.py::OptModelMultiCloneEcDNA`, a thin subclass of
`rck/core/ilp_gurobi.py::OptModelMultiClone`).
- With its options off, it builds RCK's model: the same variables, constraints and objective
  (`test_defaults_build_rcks_model`). On demo5, the route's `rck.scnt.tsv` and `rck.acnt.tsv` are
  byte-identical to `rck_run.py --no-pre`'s.
- The only always-on difference is the empty-novel-adjacency case in §11.1.

**Outputs** (`--outdir`):
- `rck.scnt.tsv` and `rck.acnt.tsv`, as `rck_run.py` writes them;
- `summary.json`: `model: multiclone`, status, objective, flipped fragments, and the novel
  adjacencies used per clone;
- `gurobi_multiclone.log`, `gurobi_greedy.log`;
- `decompositions.tsv` and, with `--truth-dir`, the evaluation (§11.3).

```bash
rck-ecdna --model multiclone --rck-input subgraphs/demo5/amplicon_0/rck_input \
    --outdir rck_mc_demo5 --time-limit 60 --threads 2
```

### 11.1 Deviations from `OptModelMultiClone`

| | RCK `OptModelMultiClone` | `OptModelMultiCloneEcDNA` | Default | Why |
|---|---|---|---|---|
| 1 − fp presence rule | over all novel adjacencies; $0 \ge 1-\mathrm{fp}$ (infeasible) when there are none | over observed novel adjacencies only; skipped when there are none | always (identical to RCK whenever there are observed novel adjacencies and imputed ones are off) | imputed candidates are not SV calls; the empty case is an RCK edge case (an amplicon whose junctions are all imputed, e.g. demo5 `amplicon_2`) |
| imputed adjacencies | n/a | `rck.adj.imputed.tsv` joins the graph; cost `--mc-imputed-penalty` per adjacency used (any labelling, any clone) | off | recovers a circularising junction the SV caller missed (plan 4.2) |
| novel adjacency cost | none | `--mc-novel-penalty` per observed novel adjacency used | 0 (off) | kept as an option (plan decision 3) |
| junction read term | none (RCK has no read data term) | `--mc-junction-weight` $\cdot\,\lvert v_a - \sum_k u_k \sum_\varphi N_{k,a,\varphi}\rvert / \sqrt{\max(v_a,1)}$ for observed novel adjacencies only (the DIAG has no read count per reference adjacency). $v_a$ = `#NOV`/θ; θ comes from `--junction-reads-per-copy`, else `scale.json` next to the DIAG; $u_k$ from `clone_proportions.tsv` | 0 (off) | the bulk route's `--alpha-novel` term, per clone (plan 4.2). Its units are copies, while RCK's objective is copies × `length_10`, so at equal weight it acts as a tie-break |
| extra ends | excess only at telomere positions, free | also at `--mc-linear-ends boundaries` (a total copy-number step of at least max(`--boundary-min-copies`, `--boundary-frac` × the larger) in any clone) or at `all` extremities, at `--mc-linear-end-penalty` per copy | `none` (off) | explains an amplicon without a closing junction as linear (HSR-like) instead of charging the deltas (plan 4.2) |
| solver settings | `--run-g-mip-gap 0.015`, `--run-g-threads 4` | `--mip-gap 1e-4`, `--threads` | — | the `rck-ecdna` Gurobi options; not a model change |

Unchanged from RCK:
- the variables, the balance at every extremity (excess at telomeres) and the objective
  $\sum_k \sum_s L_s(\delta_{k,s,A} + \delta_{k,s,B})$;
- the fragment flips;
- the bounds (`--scnb` is now readable, see §8);
- the reference labelling (AA/BB) and one labelling per novel adjacency shared by all clones;
- the reciprocal-location rules;
- the molecule, general and labelling groups.

The clone proportions $u_k$ enter only through the junction term, and later the bulk proportions.

### 11.2 Why there is no chromosome/ecDNA layer in the model (risk 1)

The first step-3 version split every clone's copy numbers into a chromosome layer, with RCK's
telomere rule, and an ecDNA layer, a strict circulation (plan, original decision 2). Probing every
optimal solution showed:
- **The flanks pinned linear copies.** The chromosome layer's telomere excess never exceeded the
  flank copy number.
- **But a circulation balances in either layer**, and RCK's objective sees only the totals. So
  ecDNA copies could sit in the chromosome layer as a closed cycle, at the same objective.
  - On the synthetic DIAG, the ecDNA layer ranged from 0 to all $e$ copies.
  - On demo5, Gurobi put all 20 copies in the chromosome layer.

Original RCK cannot separate the two. It infers copy numbers, not chromosomes, and its only rule
about ends (telomere positions) does not apply to a cycle. AmpliconArchitect and CoRAL do not
split copy numbers into layers either: they decompose the amplicon into walks, a closed walk is a
cycle, and a linear walk must end at a source/sink. They take the heaviest walk first. The route
therefore keeps RCK's model as it is, and step 4 (plan 4.3) decomposes the per-clone karyotypes
heaviest first:
- cycles are ecDNA species;
- linear walks may end only at RCK's telomere positions and are chromosomal.

The commits with the layer version are `8e6fb74` and `d3b0448`; `ea99498` removed it.

### 11.3 Decomposition: heaviest walk first (`rck/ecdna/greedy.py`)

*Step 4 (2026-11-13).* The decomposition takes every clone's RCK karyotype and splits it into
walks shared by the clones, heaviest first. This is AmpliconArchitect's rule
[Deshpande et al. 2019](pipeline/references.md#deshpande2019), run with the one-walk MILP of
CoRAL's greedy mode [Zhu et al. 2024](pipeline/references.md#zhu2024)
(`coral/cycle_decomposition.py::maximize_weights_greedy`, CoRAL 3.1.1).

**Input.** RCK's solution for every clone $k$, as copy numbers on phased edges:
- segment haplotypes, $c_{k,s,h}$;
- adjacency labellings, $N_{k,a,\varphi}$;
- *linear ends*: the copy-number excess at an extremity, which RCK allows only at telomere
  positions (and at `--mc-linear-ends` extremities). These use the phased-edge vocabulary of
  §5.

The remaining copy number $R_k(e)$ starts at these values.

**One round** solves one MILP for a single walk $\sigma$:
- integer multiplicities $m_e \le$ `--max-multiplicity` (binary expansion);
- alternation and balance at every phased extremity;
- connectivity through lazy cuts;
- no ends (a **cycle**) or exactly two (a **path**, through one virtual end node, as in §5.1);
- integer copies $\psi_k \ge 0$ per clone, with $\psi_k m_e \le R_k(e)$ on every edge (products
  linearised exactly).

**Which walk is heaviest** (`--mc-heaviest`):
- `copies` (default), AmpliconArchitect's rule. First maximise $\sum_k \psi_k$. With that fixed,
  maximise CoRAL's length-weighted copy number $\sum_k \psi_k \sum_{s} m_s \ell_s$, plus CoRAL's
  long-read bonus. Ties then go to the walk with fewer edge traversals (a tie-break below one bp
  per traversal).
- `length`, CoRAL's rule: the second objective alone.

The walk is subtracted from every clone, and the next round starts. The loop stops in any of
these cases:
- nothing is left (RCK's karyotype is an exact balanced integer flow, so this is the normal end);
- the remaining length-weighted copy number is at most `--mc-stop-fraction` of the total
  (CoRAL's `p_total_weight`; off by default);
- `--mc-max-walks` walks were taken.

**Calls.**
- Cycles are ecDNA candidates. A cycle lighter than `--min-species-weight-frac` of the heaviest
  cycle is `low_weight`, and the amplicon call uses the cycles only (`--min-ecdna-copies`), as in
  §5.10.
- Paths are chromosomes (`status` `chromosomal`). The evaluation skips them.
- A cycle's `weight` is its bulk copies per cell, $\sum_k u_k \psi_k$.

**Why `copies` is the default.** The DIAG has 1 Mb flank segments, because it needs them as
chromosome anchors (plan 3.2). CoRAL's graphs hold only the amplified intervals. With CoRAL's
length weight, a chromosome path through both flanks can therefore outweigh a short ecDNA.
- It is then extracted first, and it loops through the amplicon up to `--max-multiplicity`
  times (an HSR-like walk), taking copies from the cycle.
- On the synthetic DIAG (100 kb flanks, a 20 kb ecDNA of 5 copies), `length` takes such a path
  first, while `copies` takes the 5-copy cycle (`test_greedy.py`).
- On demo7, `length` scores a mean composite of 0.256 and 0/4 exact ψ, against 0.747 and 2/4
  for `copies`.
- With the labelling below held to the truth, `copies` is exact on demo7 (composite 1.0, 4/4 ψ),
  while `length` stays at 0.256.

| | CoRAL greedy (`maximize_weights_greedy`) | `rck/ecdna/greedy.py` |
|---|---|---|
| input | one breakpoint graph with fractional copy numbers | every clone's integer RCK karyotype |
| heaviest | length-weighted copy number + long-read bonus | `copies`: most copies first, then CoRAL's weight; `length`: CoRAL's weight |
| weights | one fractional weight | integer $\psi_k$ per clone, shared walk |
| ends | source edges | RCK's telomere positions (and `--mc-linear-ends`) |
| multiplicity | free | at most `--max-multiplicity` |
| long reads | ordered subpaths | RCK molecule groups (`rck.adjgroups.tsv`): satisfied when the walk uses all their novel adjacencies; bonus $\max(\alpha R / n_{\text{unsatisfied}}, 1)$, $\alpha$ = `--mc-read-alpha` (0.01) |
| stop | remaining $\le$ 10 %, or weight below resolution | nothing left (or `--mc-stop-fraction`, `--mc-max-walks`) |
| ties | solver | fewer edge traversals |

**Outputs.**
- `decompositions.tsv` has the bulk route's columns (Appendix D), one species per walk in
  extraction order, named `ec1, ec2, …` for cycles and `path1, …` for paths, plus:
  - `kind` (`cycle`/`path`) and `order`;
  - `psi_<clone>`;
  - `proportion_<clone>`: $\psi_{\sigma,k}$ over the clone's cycles;
  - `proportion_bulk`: $\sum_k u_k \psi_{\sigma,k}$ over the cycles.
- `summary.json` has a `decomposition` block.
- With `--truth-dir`: `evaluation.*` (structure, cycles only) and `evaluation_psi.tsv`, which
  scores ψ and the proportions against the `ecdna_psi` rows (`docs/pipeline/evaluation.md`).

**Results (2026-11-13; 60 s limit, 2 threads, about 1 s each).**
- demo5: the 20-copy cycle first, then one telomere-to-telomere path per haplotype and region.
  Composite 1.0, and ψ is exact.
- demo7 (2 clones, nested species with ψ 26/13 and 14/7): composite 0.747, 2/4 exact ψ.
  - The inner species (14/7) comes out merged with a pass round part of the outer one, and the
    outer species keeps 12/6 of its 26/13 copies.
  - The cause is RCK's labelling, not the greedy. RCK never reads `junction_hap`, so it labelled
    the inner species' junctions sim00005–8 AB/BA, while the truth (and `junction_hap`) say AA.
    It reaches objective 0 either way.
  - The karyotype therefore has no inner cycle on one haplotype.
  - Using `junction_hap` would be a new deviation from RCK, not made
    (`multiclone_NEXT_STEPS.md`).

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
| `--model` | `bulk` (default; this document §2–§7) or `multiclone` (§11) |
| `--subgraph-dir` | `nodes.tsv`, `edges.tsv` (and optionally `phase_blocks.tsv`, `segment_blocks.tsv`) of one amplicon; with `scale.json` etc. it is the only input needed. Required with `--model bulk` |
| `--rck-input` | `--model multiclone`: the annotated DIAG folder `<subgraph>/rck_input/` (§11) |
| `--mc-nas-fp`, `--mc-group-m-fp`, `--mc-segment-length-attr` | `--model multiclone`: RCK's `--run-nas-fp` (0.1), `--run-group-m-default-fp` (0.1), `--run-segment-length-attr` (`length_10`) |
| `--mc-imputed-penalty`, `--mc-novel-penalty`, `--mc-junction-weight`, `--mc-linear-ends`, `--mc-linear-end-penalty` | `--model multiclone`: the ecDNA additions of §11.1, all off by default |
| `--mc-heaviest`, `--mc-max-walks`, `--mc-read-alpha`, `--mc-stop-fraction`, `--mc-no-decompose` | `--model multiclone`: the heaviest-first decomposition (§11.3): `copies` (default) or `length`; at most 50 walks; CoRAL's α 0.01; stop share 0 (use everything up); skip it. `--max-multiplicity`, `--min-ecdna-copies`, `--min-species-weight-frac` and `--max-orderings` apply as in the bulk route |
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
| `test_multiclone.py` | yes (small) | `--model multiclone` (§11) on tiny synthetic DIAGs: the default model equals `OptModelMultiClone` (variables, constraints, objective, solution); 1 clone and 2 clones (one labelling and one flip shared by the clones); a flipped phase block; no novel adjacencies (feasible); each off-by-default addition; the CLI outputs and `--model` checks |
| `test_greedy.py` | yes (small) | the heaviest-first decomposition (§11.3): one cycle in two clones (ψ per clone, then one telomere-to-telomere path per haplotype), two species with different ψ per clone, `copies` vs `length` (the looping path), `decompositions.tsv` columns and the evaluation reader skipping paths |
| `test_demo3.py` | yes (full licence) | `amplicon_0` against the demo3 truth, rank 1: cycle order matches, haplotype Jaccard ≥ 0.95, junction recall ≥ 0.9, copy-number error ≤ 30%; the same with long-read paths (≥ 90% of read support explained) |

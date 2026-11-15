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
| Allele split $o_{j,A} = o_j(1-\mathrm{BAF}_j)$ | §3.3 | (a) | the BAF–allele copy-number relation used by ASCAT [Van Loo et al. 2010](pipeline/references.md#vanloo2010) and HATCHet [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020) |
| Subgraph-only inputs (`scale.json`, `CN_*`, `cn_profile.tsv`) | §3 | (c) new | `a987ffa` |
| Two-layer flow MILP | §4 | (c) modified from RCK `OptModelMultiClone` | `7d8dff9` |
| Background slack at terminals | §4.3 | (a) RCK's telomere relaxation, moved to the background layer | `7d8dff9` |
| Imputed candidate adjacencies with a usage penalty | §4.5 | (c) new | `7d8dff9` |
| Length-capped, $1/\sqrt{o_j}$-normalised data terms; junction read term | §4.5 | (c) new | `7d8dff9` |
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
extremity has exactly one segment. Nodes left without any segment or adjacency are dropped.
`ECGraph.validate` logs any node that still breaks this rule as a warning; it does not repair it.

**Notation.** RCK's DIAG is $D(R,\tilde A_N)=(V,E)$. For segments $j\in\{1,\dots,m\}$ and
haplotypes $H\in\{A,B\}$, the haplotype-specific segment is $j^H=[j^H_t,j^H_h]$, with
extremities $j^H_\sigma$, $\sigma\in\{t,h\}$ (tail = left end, head = right end). The vertices
$v\in V$ are these haplotype-specific extremities. $E=E_S\cup E_A$: $E_S$ holds the segment
edges $j^H$; $E_R\subseteq E_A$ holds the reference adjacency edges; and $E_N=E_A\setminus E_R$
holds the novel adjacency edges. The adjacency set is $\mathcal A=A(R)\cup H(\tilde A_N)$:
the reference adjacencies $A(R)$ and the labelled versions of the measured, unlabelled novel
adjacencies $\tilde A_N$. An unlabelled adjacency $a=\{j_\sigma,k_{\sigma'}\}$ has labelled
versions $h(a)=\{\{j^H_\sigma,k^{H'}_{\sigma'}\}: H,H'\in\{A,B\}\}$. $T(R)$ is the set of
telomeres. In this document a labelling is written $\varphi\in\{AA,AB,BA,BB\}$ (the names come
from RCK's code, `rck/core/structures.py::Phasing`), and the labelled edge is $a^\varphi$.
In the code, a node of `ECGraph` is an unlabelled extremity $j_\sigma$ (subgraph side `_-` = $t$,
`_+` = $h$). A phased node `(nid, H)` is the vertex $j^H_\sigma$; a phased edge `("S", sid, H)` is
$j^H$; and `("J", aid, φ)` is $a^\varphi$ (`rck/ecdna/phased.py`). The terminals $\mathcal T$
(§2.4) play the role of $T(R)$. The linear ends $\mathcal L$ (§4.4) are new.

### 2.1 From subgraph edges to DIAG elements

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::build_graph`.

In the table below, some properties were renamed to match the terminology in RCK (Aganezov & Raphael, 2020). The code calls them by other names, which is why there is a (code: ...) section.

| `edges.tsv` `EDGE_TYPE` | DIAG element | Observation carried |
|---|---|---|
| `reference` | segment (code: segment edge, reference edge) | `WEIGHT` = mean depth; `BAF`, `n_snps` (the BAF is ignored when `n_snps` = 0); `CN_TOTAL`, `CN_A`, `CN_B` in subgraph mode |
| `extremity` | reference adjacency between the two extremities at a breakpoint (code: reference join) | `WEIGHT` = `#REF`, reads spanning the reference join |
| `variant` | observed novel adjacency (code: junction, variant edge) | `WEIGHT` = `#NOV`; `junction_hap`. `NOV_A`, `NOV_B` are read but no model step uses them |
| `imputed_variant` | candidate novel adjacency proposed by stage 2's imputation (code: imputed edge) | none (its `WEIGHT` is kept as `imputed_support`, not as evidence; only the AA-format export writes it) |

Rows of any other `EDGE_TYPE` are skipped, with a note in the normalisation log. A node that appears only in `edges.tsv` is added, with its side taken from its `_+` / `_-` suffix.

A reference adjacency touching an `imputed_ext*` extremity carries no read count. During subgraph construction, the algorithm finds amplified segments from read depth (section 2.4 of subgraph_construct.md). Each segment end is connected to an observed extremity within --tolerance bp. When no observed extremity is close enough, it creates an imputed extremity, named imputed_ext_<chr>_<n>, which is a breakpoint inferred from a change in depth, with no reads to back it. Its #REF is written as 0 only because the field needs a value. RCK should not use that 0, though, because if the flow model read the 0 as data, it would try to give the reference adjacency 0 copies, which cuts the genome incorrectly. So subgraph_io.py sets reads=None for any reference adjacency with an end whose name starts with imputed_ext. The model then has no read-based target there. Its copy number is set only by copy-number balance with the neighbouring segments and adjacencies.

### 2.2 Phantom segments

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::_contract_phantoms`.

The simulator names both sides of a breakpoint as separate extremities (`x+` and `x+1-`). Between
them sits a segment of 1 bp (`seg_hap = no_span` in the simulator's output) whose two ends each
have exactly one reference adjacency, so it only carries reference flow. A segment is contracted
when it is at most `--phantom-max-bp` (default 2) bp long, each of its ends has exactly one
reference adjacency, and both of those adjacencies lead to the end of another segment. The code
tests only this shape, not `seg_hap`. Segments are visited in genome order. On contraction:

- its two reference adjacencies become one, joining the head of the left-hand segment to the
  tail of the right-hand segment. The merged read count is the mean of the read counts the two
  have (`None` if neither has one);
- its base joins the right-hand segment;
- any other edge on its tail moves to the tail of the right-hand segment, and any other edge on
  its head moves to the head of the left-hand segment, so each edge keeps its orientation. A moved
  novel adjacency has its labellings (§2.5) recomputed, because it may have become a self-loop.

`--phantom-max-bp 0` turns this off.

*Why:* This makes the graph look like the input the original RCK algorithm expects. If the phantom segments aren't contracted, then every breakpoint adds a segment, which brings two more extremities and a reference adjacency that the MILPs must balance. That makes the model larger and the optimisation slower. Also, a novel adjacency may attach to the phantom’s end rather than to the real segment it joins. The flow balance and the decomposition’s walks then route copies through a 1 bp piece, and segments don’t line up with the truth, which causes errors in the deconvolution step.

### 2.3 Parallel adjacencies

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::_merge_parallel`.

Sometimes there are several edges between the same pair of extremities. This comes from subgraph construction’s imputation: its rules for proposing candidate novel adjacencies sometimes propose an edge that already exists.

Edges between the same unordered pair of extremities are merged (after phantom contraction):

- several reference adjacencies become one, with the mean of their read counts (counting only
  those that have one);
- several observed novel adjacencies become one, with the summed `#NOV`. The first one's id,
  `junction_hap`, labellings and `NOV_A`/`NOV_B` are kept;
- an imputed candidate parallel to a reference or observed novel adjacency is merged into it
  (it is that adjacency). If both kinds are present, it goes to the reference adjacency. Its
  weight is kept only as `imputed_support` (the largest, if several are merged). Stage 2 produces
  such edges as `inter` edges at imputed extremities and as anchors that loop back onto their own
  extremity;
- several parallel imputed candidates become one.

A reference and an observed novel adjacency on the same pair are not merged with each other.

### 2.4 Terminals

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::build_graph`.

An extremity is *terminal* when it has no segment (the outer flanks of each region), is a
chromosome-boundary node (`SOURCE = chromosome_boundary`), has no side (`ORIENT` is neither `+`
nor `-`), or has no reference adjacency. Terminals form the set $\mathcal T$. Only the background
layer of the flow MILP may start or end at a terminal, through the slack $s^{\mathrm{bg}}_v$
(§4.3). The ecDNA layer may stop only at linear ends $\mathcal L$ (§4.4). These are chosen
separately, among extremities that have a segment, so they can include terminals that have a
segment.
This is borrowed from RCK's telomeres $T(R)$. In RCK, telomeres are where paths may end, which is what terminals do for the background layer here.

### 2.5 Allowed haplotype labellings

> **Provenance:** labelling rules (a) RCK [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020); labels from read evidence (c) new,
> `c257e54`; soft mode (c) new, `7d8dff9`. **Code:** `rck/ecdna/subgraph_io.py::_junction_phasings`.

Each adjacency $a=\{j_\sigma,k_{\sigma'}\}$ gets a set $\Phi(a)$ of allowed labellings, a subset
of $h(a)$. Labelling $\varphi=(H,H')$ joins $j^H_\sigma$ to $k^{H'}_{\sigma'}$, where $j_\sigma$
is the adjacency's first end (`NODE_1`, code `u`), as in RCK (`rck/core/structures.py::Phasing`).

- Reference adjacencies: $\Phi=\{AA, BB\}$ (RCK: a reference adjacency never joins homologs).
- Observed novel adjacencies follow stage 2's read-based label `junction_hap`, read
  case-insensitively: `A` → $\{AA\}$, `B` → $\{BB\}$, `mixed` → $\{AB, BA\}$ (a haplotype
  crossing), missing or any other value → all four.
- Imputed candidates: $\{AA, BB\}$.
- A novel self-loop ($j_\sigma=k_{\sigma'}$, a fold-back) can never be AB or BA, as in RCK
  (`OptModelMultiClone.define_constraints_for_novel_labeling`); a `mixed` self-loop falls back to
  $\{AA,BB\}$.

With `--soft-junction-hap`, $\Phi(a)$ is all four labellings for observed novel adjacencies (only
AA and BB for a self-loop), and the read-based ones become the *preferred* set $\Phi^*(a)$: any
other labelling used costs `--soft-junction-penalty` (default 1.0; §4.5). An adjacency with no
label prefers all four, so it pays nothing.

### 2.6 Regions

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/subgraph_io.py::ECGraph.regions`.

The *regions* of a graph are the spans its segments cover, split wherever two segments on one
chromosome are more than 1 Mb apart. They define the flank background (§3.5) and which truth
species a run is scored against (evaluation). They also set the background exclusion zones of
§3.1–§3.2.

---

## 3. Copy-number observations

Two scales turn read counts into copy numbers (`rck/ecdna/scale.py`):

- $\gamma$, **depth per copy**: tumour read depth contributed by one copy of a locus, averaged over
  the bulk sample. Observed segment copy number $o_j = \bar d_j/\gamma$, with $\bar d_j$ the mean
  depth of segment $j$.
- $\theta$, **junction reads per copy**: reads spanning one copy of an adjacency. Observed
  adjacency copy number $o_a = \mathrm{NOV}(a)/\theta$ (novel, `#NOV`) or
  $\mathrm{REF}(a)/\theta$ (reference, `#REF`).

All copy numbers here are cell-averaged over the sample, so they need not be integers.

The baseline copy number $c_0$ is `--baseline-cn`. By default it is `scale.json`'s
`background_cn` in subgraph mode (2 if that key is missing), else 2. It is used in the windows
route for γ (§3.1), in the θ fallback (§3.2), and as the flank background when no flank data
exist (§3.5).

The input route is chosen by `--inputs` (default `auto`). `subgraph` reads only the subgraph
folder: `scale.json`, the `CN_*` edge columns, `cn_profile.tsv` and `paths.tsv`. It ignores
`--windows-dir`, `--obs-dir`, `--extremities`, `--hatchet-seg` and `--sample`. `auto` means
`subgraph` when `scale.json` exists, and `windows` otherwise.

### 3.1 Depth per copy γ

> **Provenance:** (c) new, `c257e54` (`_log_mode` fixed in `41f5804`); subgraph route `a987ffa`.
> The HATCHet route reads HATCHet's `best.seg.ucn` [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020).
> **Code:** `rck/ecdna/scale.py::estimate_depth_scale`, `::_gamma_from_hatchet`, `::_log_mode`,
> `::build_observations`.

The override `--depth-per-copy` always wins. Otherwise γ comes from the first available source:

0. **Subgraph folder** (`--inputs subgraph`, or `auto` when `scale.json` exists): γ =
   `scale.json` `depth_per_copy`, fitted by stage 2's `rck_inputs.py` over purity and ploidy from
   genome-wide depth and SNPs (documented in
   [subgraph_construct.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/subgraph_construct.md)).
1. **Override** `--depth-per-copy`. This applies in either route and needs no windows.
2. **HATCHet** (`--hatchet-seg best.seg.ucn`, windows route; it needs `--windows-dir`). The rows
   are those of `--sample` in the `SAMPLE` column, else of the first sample. Over HATCHet
   segments $\xi$ that are at least 1 Mb long, have positive copy number, lie more than 1 Mb
   from every region, and have window coverage, γ is the length-weighted median of

   ```math
   \frac{\bar d_\xi}{\sum_{i} u_i \left(\hat c_{i,\xi}+\check c_{i,\xi}\right)}
   ```

   where $\bar d_\xi$ is the mean tumour depth over $\xi$, $u_i$ is clone $i$'s proportion
   (`u_*`), and $\hat c_{i,\xi}$, $\check c_{i,\xi}$ are its two allele-specific copy numbers in
   `best.seg.ucn` (`cn_*` = `a|b`). This is HATCHet's output, which is RCK's usual input
   $\tilde C$. The denominator is the cell-averaged total copy number. This is the route without
   a matched normal (cell lines). `--depth-mode ratio` is ignored here, with a warning.
3. **Windows mode** (default without the above): 1 kb windows of the tumour sample are grouped
   into blocks $\xi$ of `--block-bp` (default 1 Mb) by start position. A block's depth
   $\bar d_\xi$ is the mean of its windows' depths. A block is *background* when it has at least
   half its windows (`block_bp`/2 divided by the median window length) and lies more than 1 Mb
   from every region. It is *balanced* when it also has at least 20 SNPs and mirrored BAF
   $\sum\min(A,B)/\sum(A+B)\ge 0.40$, summed over its SNPs with $A+B>0$. Then

   ```math
   \gamma = \frac{\operatorname{mode}_{\log}\{\bar d_\xi : \xi \text{ balanced}\}}{c_0}
   ```

   where $\mathop{\mathrm{mode}}_{\log}$ is the mode of a histogram of $\log \bar d_\xi$ (bin
   width about 0.02, range padded by 3 bin widths on each side). The histogram is smoothed by a
   Gaussian kernel of σ = 2 bins, truncated at ±6 bins, and the mode is reported at the centre of
   the highest bin. With fewer than 10 balanced blocks every background block is used. Only
   positive, finite values count. With fewer than 3 values, or values within one bin width, the
   median is used.

`--depth-mode ratio` (windows mode only) instead divides tumour by matched-normal depth window by
window. Segment and flank copy numbers are then $c_0\, \rho_j/\rho_{\mathrm{mode}}$. Here $\rho_j$
is the ratio of summed tumour to summed normal depth over the segment, and $\rho_{\mathrm{mode}}$ is
the log-mode of the block ratios over the same blocks used for γ. γ itself is still the tumour mode
above. It is still used for θ (§3.2) and wherever the windows have no coverage. This removes
coverage bias that both libraries share (GC, mappability). The amplicon-allele-analysis
simulator draws the bias independently per sample, so on simulated data the ratio only adds noise,
and `tumor` is the default.

*Why the mode:* an amplicon sample is mostly unamplified background, and its commonest copy number
on allelically balanced blocks is the cell ploidy. The 2026-10-19 fix replaced `np.convolve(...,
mode="same")`, which returned the wrong length whenever the kernel was longer than the histogram.

*Why γ:* Original RCK doesn't have a depth per copy (γ). Its input is segment copy numbers that another tool has already computed. RCK-ecDNA needs γ because its pipeline begins further upstream at raw observations. It has to convert mean read depth per segment, allele counts per segment, and spanning-read counts per novel adjacency into copy numbers.


### 3.2 Junction reads per copy θ

> **Provenance:** (c) new, `c257e54`. **Code:** `rck/ecdna/scale.py::estimate_junction_scale`.

1. The override `--junction-reads-per-copy`; else, in subgraph mode, `scale.json`
   `junction_reads_per_copy` (when present and non-zero).
2. **Background extremities** (windows route only; needs `--extremities` and windows): this uses
   the genome-wide extremities $v$ of `extremity.tsv` (position column `POS0`, else `START`).
   Take those more than 1 Mb from every region whose local copy number
   $o_v = \bar d(v\pm 25\,\mathrm{kb})/\gamma$ lies in $[0.5, 8]$. Then

   ```math
   \theta = \frac{\sum_v \mathrm{REF}(v)}{\sum_v o_v}
   ```

   It needs at least 5 such extremities, or it falls through to step 3 with a warning.
3. **Fallback** (logged as a warning; also used in subgraph mode when `scale.json` has no θ): the
   median of the positive `#REF` counts of the subgraph's own reference adjacencies (after §2.2–§2.3),
   divided by $c_0$. With none, the run stops and asks for `--junction-reads-per-copy`.

*Why θ:* Original RCK doesn't have a junction reads per copy (θ). Its input is segment copy numbers that another tool has already computed. RCK-ecDNA needs θ because its pipeline begins further upstream at raw observations. It has to convert these into copy numbers.

### 3.3 Segment observations and the allele split

> **Provenance:** allele split (a), the BAF relation used by ASCAT [Van Loo et al. 2010](pipeline/references.md#vanloo2010) and HATCHet
> [Zaccaria & Raphael 2020](pipeline/references.md#zaccaria2020); subgraph-mode columns (c) new, `a987ffa`.
> **Code:** `rck/ecdna/scale.py::build_observations`.

For segment $j$:

```math
o_j = \frac{\bar d_j}{\gamma}, \qquad o_{j,A} = o_j\,(1-\mathrm{BAF}_j), \qquad o_{j,B} = o_j\,\mathrm{BAF}_j
```

where $\bar d_j$ is the overlap-weighted mean tumour depth of the windows over $j$. When there are
no windows (subgraph mode, or an override without `--windows-dir`), or the windows have no bins
there, it is the subgraph's `WEIGHT`. $\mathrm{BAF}_j$ is the phase-oriented B-allele fraction
from stage 2, clipped to $[0,1]$. In subgraph mode, $o_j$, $o_{j,A}$, $o_{j,B}$ are the
`CN_TOTAL`, `CN_A`, `CN_B` columns. If `CN_A`/`CN_B` are missing, the BAF split of `CN_TOTAL` is
used. If `CN_TOTAL` is missing, `WEIGHT`/γ is used. With `--depth-per-copy` in subgraph mode, the
`CN_*` columns are not used: $o_j$ = `WEIGHT`/γ, split by the BAF. A segment with fewer than
`--min-snps-allele` (default 5) SNPs gets a total-only observation.
`--min-snps-allele 100000000` makes every observation total-only, a diagnostic used below.

### 3.4 Adjacency observations
(c) new to the RCK model; read-based adjacency copy numbers follow AmpliconArchitect’s read model for breakpoint-graph edges [Deshpande et al. 2019], also used by CoRAL [Zhu et al. 2024]; θ estimation and the haplotype split are new.

$o_a = \mathrm{NOV}(a)/\theta$ or $\mathrm{REF}(a)/\theta$ for every reference or novel adjacency
with a read count. Imputed candidates and reference adjacencies at an `imputed_ext` extremity
(§2.1) have none. A read count of 0 on any other adjacency is an observation of 0. This is the
*observed adjacency copy number*. The `CN` column of `edges.tsv` holds the same quantity (it is
written by stage 2; `rck-ecdna` does not read it). The flow's integer value (§4) is the *inferred*
adjacency copy number in RCK's sense.

### 3.5 Flank background

(c) new to the RCK model

For each region $[\mathrm{start}_r,\mathrm{end}_r]$ on one chromosome, the flank copy number is
the mean copy number over $[\mathrm{start}_r-1\,\mathrm{Mb},\, \mathrm{start}_r-10\,\mathrm{kb}]$
(clipped at 0) and $[\mathrm{end}_r+10\,\mathrm{kb},\, \mathrm{end}_r+1\,\mathrm{Mb}]$. The copy
number comes from the windows (depth/γ, or the ratio of §3.1), or from `cn_profile.tsv` in subgraph
mode. The two flanks are weighted by their covered bins (windows) or covered kb (`cn_profile.tsv`).
If neither flank has data, the value is $c_0$. Every segment in the region gets this value as
$o^{\mathrm{bg}}_j$, the background prior of §4.5.

*Why it was added:* RCK models the whole genome. Every chromosome’s copies are in the model end to end, so the copy number under an amplicon is fixed by the same chromosome’s flow elsewhere. RCK-ecDNA works on a subgraph cut out of the genome, so it loses that context. The flanks stand in for the rest of the chromosome.

### 3.6 Per-phase-block observations (`--block-flip unoriented`)

> **Provenance:** (c) new, `a422c39`. **Code:** `rck/ecdna/scale.py::piece_observations`,
> `rck/ecdna/subgraph_io.py::read_segment_blocks`.

The segment BAF is on the true haplotype scale only where stage 2 could orient every phase block the
segment touches. With `--block-flip unoriented`, each segment with rows in `segment_blocks.tsv` gets
one observation per *piece* $\pi$, its part inside one phase block $b(\pi)$:

- pieces are matched to segments by their parent edge key (`NODE_1 NODE_2 CHROM_1 POS_1 POS_2`),
  through phantom contraction (a contracted phantom's pieces go to the segment that absorbed its
  base). Rows that match no segment are ignored, with a warning;
- the piece total $o_\pi$ is its `CN_TOTAL`, rescaled by (`scale.json` γ)/γ under
  `--depth-per-copy`. In the windows route it is the windows' copy number over the piece. If
  neither is available, it is the segment's $o_j$. $(o_{\pi,A},o_{\pi,B})$ are its `CN_A`/`CN_B`.
  In windows mode, with `--depth-per-copy`, or when those columns are missing, the total is split
  by the piece BAF `B_count_obs`/(`A_count_obs`+`B_count_obs`). A piece with fewer than
  `--min-snps-allele` SNPs, or with no allele counts, is total-only;
- its share $\sigma_\pi = L(\pi)/\sum_{\pi'\in j}L(\pi')$ is its fraction of the summed length
  of the segment's pieces (`piece_end` − `piece_start`), not of the segment's own length;
- `CN_A`/`CN_B` are first returned to the observed frame, using `segment_blocks.tsv`'s own
  `orientation` column. A piece in a block with orientation 1 or −1 in `phase_blocks.tsv` is
  then *fixed* and taken in the oriented frame (A/B swapped for −1). A piece in a block with
  orientation 0, or one missing from `phase_blocks.tsv`, is taken in the observed frame. Its A/B
  may be swapped by the block's flip $\chi_b$ (§4.6). If `phase_blocks.tsv` is absent, the
  orientations come from `segment_blocks.tsv` (with a warning). If it is empty, every block counts
  as unoriented.

Without `segment_blocks.tsv` the run logs a warning and uses the segment observations, as with
`--block-flip off` (the default).

This is built on the original RCK algorithm, which has a free A/B flip per fragment $f$. Without `--fragments`, each segment is its own fragment, so the flip is per segment. RCK’s input is allele-specific, not haplotype-specific. So RCK decides, per fragment, which allele belongs to which haplotype, as part of its MILP. In the RCK-ecDNA algorithm, our observations are partly phased. Within a phase block, read-based phasing ties the segments’ A and B together, and only the block as a whole can be flipped relative to its neighbours. So, this step is RCK’s allele-to-haplotype assignment, constrained by the phasing information we have. The phase block $b$ plays the role of RCK's fragment $f$. The constraint makes the decision much smaller than RCK’s, with only one choice per unoriented block instead of one per segment.

---

## 4. Flow MILP: inferred copy numbers

> **Provenance:** (c) modified from RCK's `OptModelMultiClone` [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020), `7d8dff9`; linear
> ends and `bg_balance_weight` `10416aa`; block flips `a422c39`.
> **Code:** `rck/ecdna/ilp_circular.py::OptModelCircular` (`build`, `_big_m`, `_balance`, `_abs_dev`, `_objective`,
> `_infinite_sites`, `_piece_terms`, `block_report`, `linear_end_nodes`). Solved with Gurobi
> [Gurobi](pipeline/references.md#gurobi).

The flow MILP infers integer copy numbers per segment haplotype and per adjacency labelling, split
into two layers: an **ecDNA layer**, strictly balanced so that it decomposes into closed walks, and
a **background layer** for the chromosomal copies under the amplicon.

### 4.1 What RCK solves

Section 4.1 is for reference and does not describe the modified RCK-ecDNA pipeline.

For reference, RCK's model (`rck/core/ilp_gurobi.py::OptModelMultiClone` at upstream `4c6c288`;
`rck/version-fix` (`0b41976`) only moves its constraints to Gurobi's `addLConstr` and leaves the
model unchanged) infers, for each clone $i$, integer copy numbers $a_{i,j}=\mu_i(j^A)$ and
$b_{i,j}=\mu_i(j^B)$ of every segment haplotype, and $\mu_i(a^\varphi)$ of every labelled adjacency
$a^\varphi$. The input is the allele-specific (unlabelled) segment copy numbers
$\hat c_{i,j},\check c_{i,j}$ (from HATCHet, TitanCNA, Battenberg or ReMixT).

**Paper.** RCK minimises the distance (Eq. 11)

```math
\lVert C-\tilde C\rVert_{\mathcal F}=\sum_{f\in\mathcal F}\ \min_{\{d,d'\}=\{\hat c,\check c\}}\ \sum_{j\in f}\sum_i \left(\lvert a_{i,j}-d_{i,j}\rvert+\lvert b_{i,j}-d'_{i,j}\rvert\right)L(j)
```

subject to balance at every vertex $v$ (Eq. 5),
$\mu_i(e_S(v))=\mu_i(e_R(v))+\sum_{e\in E_N(v)} l(e)\,\mu_i(e)$, with excess allowed only at
telomeres (Eq. 6), exactly one labelled version of every measured novel adjacency
($\lvert h^+_E(a)\rvert=1$), adjacency utilisation $P$, molecule groups, and the generalised
infinite-sites constraints.

**Code.** `OptModelMultiClone` minimises (`define_constraints_on_deltas`, `define_objective`)

```math
\min \sum_{i}\sum_{j} L_{10}(j)\left(\delta_{i,j,A}+\delta_{i,j,B}\right),\qquad
\delta_{i,j,A} \ge \left\lvert \mu_i(j^A) - \left(b_{f(j)}\,\hat c_{i,j} + (1-b_{f(j)})\,\check c_{i,j}\right)\right\rvert,\quad
\delta_{i,j,B} \ge \left\lvert \mu_i(j^B) - \left((1-b_{f(j)})\,\hat c_{i,j} + b_{f(j)}\,\check c_{i,j}\right)\right\rvert
```

with $\delta_{i,j,H}\in\mathbb Z_{\ge 0}$ (each $\lvert\cdot\rvert$ written as two inequalities) and
$b_f\in\{0,1\}$ the allele flip of fragment $f$ ($b_f=1$ keeps the input order, $b_f=0$ swaps it).
The binary $b_f$ is the paper's $\min$ over $\{d,d'\}$: the minimisation picks the better
assignment for each fragment. The weight $L_{10}(j)$ is the segment attribute chosen by
`--run-segment-length-attr`, default `length_10` $=\lceil (L(j)-1)/10\rceil$, i.e. about
$\lceil L(j)/10\rceil$, in place of the paper's $L(j)$. The constraints are:

- **balance** at every vertex $v$ (`define_constraints_on_nodes`, `get_lin_expr_for_node_haplotype`):
  $\mu_i(e_S(v))-\mu_i(e_R(v))-\sum_{e\in E_N(v)} l(e)\,\mu_i(e)=0$, or $\ge 0$ when $v$ is a
  telomere (Eq. 5 and 6). The telomeres are the chromosome ends of the input segments, plus
  `--telomere-positions` / `--telomere-segments`, plus (unless `--telomere-positions-no-reciprocal`)
  the reference partner of each telomere;
- **adjacency copy numbers** (`define_big_m_constraints_for_adjacencies_copy_numbers`):
  $\mu_i(a^\varphi)=p_{i,a,\varphi}\,\mu'_i(a)$, where $p_{i,a,\varphi}\in\{0,1\}$ says that labelling
  $\varphi$ of $a$ is present in clone $i$, and $\mu'_i(a)\ge 1$ is one integer per novel adjacency and
  clone (shared by its labellings) and one per labelling for a reference adjacency (code `N'`). The
  product is linearised with $M_a$ = the largest upper copy-number bound of the two segments at
  $a$'s ends: $\mu_i(a^\varphi)\le M_a\,p_{i,a,\varphi}$, $\mu_i(a^\varphi)\le\mu'_i(a)$,
  $\mu_i(a^\varphi)\ge\mu'_i(a)-M_a(1-p_{i,a,\varphi})$, $\mu_i(a^\varphi)\ge 0$;
- **labelling** (`define_constraints_for_ref_labeling`, `define_constraints_for_novel_labeling`):
  $p_{i,a,AB}=p_{i,a,BA}=0$ for a reference adjacency. For a novel adjacency,
  $p_{a,\varphi}=\bigvee_i p_{i,a,\varphi}$ (Gurobi OR constraint), $\sum_\varphi p_{a,\varphi}\le 1$
  (at most one labelling over all clones), and a novel self-loop has $p_{a,AB}=p_{a,BA}=0$. The
  paper's "exactly one" becomes "at most one" plus the utilisation constraint;
- **utilisation** (`define_constraints_nov_adjacency_overall_presence`):
  $\frac{1}{\lvert\tilde A_N\rvert}\sum_{a\in\tilde A_N}\sum_\varphi p_{a,\varphi}\ge 1-\mathrm{fp}$, with
  $\mathrm{fp}=1-P$ = `--run-nas-fp` (default 0.1, so $P=0.9$);
- **reciprocal locations** (`define_constraints_on_each_location`, always on): for every reference
  adjacency $\{u,v\}$ with novel adjacencies at both $u$ and $v$, and every pair $a$ at $u$, $a'$ at
  $v$, the labellings that put $a$ on A at $u$ and $a'$ on B at $v$ sum to at most 1, and likewise
  with A and B swapped. This is the paper's homologous-reciprocal-extremity-exclusivity. Extremity
  exclusivity is not a constraint of this model; `rck_run` only checks it after solving
  (`--post-check-labeling`);
- **copy-number bounds** (`define_segment_copy_number_boundary_constraints`): lower and upper bounds
  per segment haplotype from the input or `--scnb`, missing ones filled by `--pre-scnb-strategy`
  (default `uniform-min-max`: $[0,10]$, `--pre-scnb-uniform-min` / `--pre-scnb-uniform-max`). They
  follow the flip: $\mu_i(j^A)\le b_f\,\mathrm{ub}_A+(1-b_f)\,\mathrm{ub}_B$, and so on;
- **adjacency groups** (`define_constraints_for_adjacency_groups_*`): a molecule group $u$ has, in
  each clone where its indicator is 1, at least $(1-\mathrm{fp}_u)\lvert u\rvert$ labelled versions
  present, and its indicator is 1 in at least one clone ($\mathrm{fp}_u$ default
  `--run-group-m-default-fp` 0.1; the paper requires all of $u$); general groups use
  `--run-group-n-default-fp` (0.1); a labelling group puts all novel adjacencies at its extremities
  on one haplotype (preprocessing adds inferred groups for short novel adjacencies).

Gurobi's relative gap is `--run-g-mip-gap` (0.015). Above `--run-presolve-nas-threshold` (5000)
novel adjacencies, `rck_run` first solves each chromosome on its own as a warm start.

A fragment is a span on which the input copy-number tool inferred allele-specific copy numbers, so
its A/B labelling is arbitrary and shared by its segments (`docs/Segments.md`, "Segments vs
fragments"). In `rck_run`, with preprocessing on (the default), the fragments are the segments before
they are split at novel-adjacency and telomere positions, and a `--fragments` file is overwritten by
them; `--fragments` is used only with `--no-pre`.

### 4.2 Variables in RCK-ecDNA

Notation: a vertex $v$ is a haplotype-specific extremity $j^H_\sigma$, so every slack below exists
once per haplotype. $\Phi(a)$ is the set of allowed labellings of adjacency $a$ and $\Phi^*(a)$ the
preferred (read-based) ones (§2.5). Reference adjacencies and imputed candidates have
$\Phi=\{AA,BB\}$. $\mathcal T$ is the set of terminals (§2.4) and $\mathcal L$ the linear ends
(§4.4).

| Symbol | Domain | Meaning | RCK or new | Code |
|---|---|---|---|---|
| $\mu^{ec}(j^H)$ | $\mathbb Z\cap[0,U_j]$ | ecDNA-layer copy number of segment $j$ on haplotype $H\in\{A,B\}$. This is similar to $\mu_i(j^H)$ in the original RCK (integer copy number of segment $j$ on haplotype $H$ in clone $i$). The clone index is replaced by a layer index. Although RCK has one karyotype per clone, RCK-ecDNA has one bulk sample split into two components that add up: $\mu(j^H) = \mu^{ec}(j^H) + \mu^{bg}(j^H)$ | RCK's $\mu_i(j^H)$, per layer | `xe` |
| $\mu^{bg}(j^H)$ | $\mathbb Z\cap[0,U^{bg}]$ | background copy number. See above. | new layer | `xb` |
| $\mu^{ec}(a^\varphi)$ | $\mathbb Z\cap[0,M_a]$ | ecDNA-layer copy number of adjacency $a$ with labelling $\varphi\in\Phi(a)$, for every adjacency (reference, novel or imputed). This is similar to $\mu_i(a^\varphi)$: copy number of adjacency $a$ with labelling $\varphi$ in clone $i$. The clone index is replaced by a layer index: $\mu(a^\varphi) = \mu^{ec}(a^\varphi) + \mu^{bg}(a^\varphi)$. | RCK's $\mu_i(a^\varphi)$, per layer | `ye` |
| $\mu^{bg}(a^\varphi)$ | $\mathbb Z\cap[0,U^{bg}]$ | background copy number; only for reference adjacencies and (unless `--no-background-novel`) observed novel ones, never imputed | new layer | `yb` |
| $p_{a,\varphi}$ | $\{0,1\}$ | labelling $\varphi$ of novel or imputed $a$ is used. This is the same variable and use case as RCK’s presence indicator per novel-adjacency labelling. New uses: it carries the `--novel-penalty`, `--imputed-penalty` and `--soft-junction-penalty` costs, and it no longer feeds RCK’s at-least-$1-\mathrm{fp}$ usage constraint. Reference adjacencies have no presence indicator. | RCK | `pv` |
| $s^{bg}_v$ | $\mathbb Z\cap[0,4U^{bg}]$ | background copies ending at terminal $v\in\mathcal T$. $s^{bg}$ is similar to the telomeres used in RCK and the source and sink nodes in CoRAL and AmpliconArchitect. In RCK-ecDNA, $s^{bg}$ is allowed at subgraph terminals. | new (RCK's telomere excess) | `t` |
| $s^{ec}_v$ | $\mathbb Z\cap[0,U_{j(v)}]$ | ecDNA-layer copies ending at linear-end extremity $v\in\mathcal L$ (§4.4). This is based on AA’s source edges and CoRAL’s s–t paths, which allow linear amplicons to end instead of becoming circles. RCK-ecDNA adopts this same method here to let the ecDNA layer explain a linear amplicon, at a cost per copy-end (`--linear-end-penalty`) so circles stay the default. | new | `te` |
| $\chi_b$ | $\{0,1\}$ | phase block $b$'s A and B are swapped (§4.6). This is similar to $b_f$ in the original RCK, which was the allele flip of fragment $f$. In RCK-ecDNA, it is the same variable and has the same mechanism: there is one binary flip shared by every segment in a span whose A/B labelling is arbitrary. The span is a phase block instead of a copy-number fragment, and only blocks with allele data that subgraph construction couldn’t orient get a flip (also oriented ones with `--flip-oriented-penalty`). | RCK's $b_f$, per phase block | `flip` |

$U_j=\lceil 1.5\max(o_j,1)\rceil+3$ (`ub_factor`, `ub_add`), with $o_j$ the total observed copy
number. This is the same as RCK’s copy-number bounds per segment haplotype (`--scnb`). It has the same role in RCK-ecDNA, but this algorithm derives the bound automatically from the observed copy number instead of needing the `--scnb` flag. \\
$U^{bg}$ is `--bg-max` (4). This is a new flag that was added. Since the background is chromosomal, its copy number should stay small and is capped at 4 copies per haplotype. My reasoning for this was that the background layer holds the chromosomal copies under an amplicon, which is normally 1 per haplotype, but clones could gain copies, so the cap leaves room up to 4 per haplotype. I wanted to keep it small to stop the background from absorbing amplified copies that should be part of the ecDNA layer. I'm open to other suggestions for the value, but so far, tests run with this cap have successfully deconvolved the ecDNA species. \\
$M_a=\max$ of $U_j$ over the segments at $a$'s ends (or $U^{bg}$ if neither end has a segment);
for a self-loop $M_a=\max(\lfloor M_a/2\rfloor,1)$ (`_big_m`). This is based on the same big-M constraint from the original RCK algorithm, where $M_a$ is the largest upper bound of the two end segments; the self-loop halving and the $U^{bg}$ fallback are new.  \\
$\mu(j^H)=\mu^{ec}(j^H)+\mu^{bg}(j^H)$ and $\mu(a)=\sum_{\varphi\in\Phi(a)}\left(\mu^{ec}(a^\varphi)+\mu^{bg}(a^\varphi)\right)$ are the inferred segment and adjacency copy numbers. \\

### 4.3 Constraints

Let $l(e)\in\{1,2\}$ be the number of times the labelled adjacency $e=a^\varphi$ occupies vertex
$v$: 2 for a self-loop labelled AA (or BB) at $v=j^A_\sigma$ (or $j^B_\sigma$), else 1
(`Adjacency.hap_at`). $e_R(v)$ is the reference adjacency labelled $HH$ at $v=j^H_\sigma$ and $E_N(v)$
the labelled novel and imputed adjacencies at $v$. The code sums over every labelled adjacency at
$v$, so it does not assume a single reference adjacency per extremity.

**Balance** (`_balance`). For every vertex $v$, if $v$ has a segment edge $e_S(v)=j^H$:

```math
\mu^{ec}(e_S(v)) = \mu^{ec}(e_R(v)) + \sum_{e\in E_N(v)} l(e)\,\mu^{ec}(e) + s^{ec}_v,
\qquad
\mu^{bg}(e_S(v)) = \mu^{bg}(e_R(v)) + \sum_{e\in E_N(v)} l(e)\,\mu^{bg}(e) + s^{bg}_v
```

and if $v$ has no segment (a region flank, always a terminal):
$\mu^{ec}(e_R(v)) + \sum l(e)\,\mu^{ec}(e) = 0$ and
$\mu^{bg}(e_R(v)) + \sum l(e)\,\mu^{bg}(e) = s^{bg}_v$. Here $s^{ec}_v=0$ unless $v\in\mathcal L$,
$s^{bg}_v=0$ unless $v\in\mathcal T$, and $\mu^{bg}(e)=0$ for an edge without a background variable.
This is RCK's haplotype-aware balance (Eq. 5), applied per layer. Because the ecDNA layer has no
telomeres (apart from linear ends), it decomposes into closed alternating walks, plus walks between
linear ends (§5, §7). The terminal slack $s^{bg}$ plays the role of RCK's telomere excess (Eq. 6),
for the background layer only. The terminals are the flanks, chromosome-boundary nodes and segment
ends without a reference adjacency (§2.4).

**Labelling and presence** (`build`). For every novel or imputed $a$ and $\varphi\in\Phi(a)$:

```math
\mu^{ec}(a^\varphi)+\mu^{bg}(a^\varphi) \le (M_a+U^{bg})\,p_{a,\varphi},\qquad
\mu^{ec}(a^\varphi)+\mu^{bg}(a^\varphi) \ge p_{a,\varphi},\qquad
\sum_{\varphi\in\Phi(a)} p_{a,\varphi}\le 1
```

(an imputed $a$ has no $\mu^{bg}$ term). At most one labelling per novel adjacency, as in RCK.
Reference adjacencies have no presence indicator: AA and BB may both carry copies, and the bounds
are the variable bounds of §4.2. RCK's minimum share of novel adjacencies used is dropped:
chromosomal SVs that touch the amplicon may stay unused.

**Infinite sites** (`--enforce-is`, off by default; `_infinite_sites`). For every vertex
$v=j^H_\sigma$ with more than one term:
$\sum_{a\ni v,\ a\ \text{novel or imputed}}\ \sum_{\varphi\in\Phi(a):\ a^\varphi\ \text{occupies}\ v} p_{a,\varphi}\le 1$.
This is the paper's extremity exclusivity [Aganezov & Raphael 2020](pipeline/references.md#aganezov2020) at a haplotype-specific extremity. It is off because chromothripsis-derived
ecDNA violates it. RCK's code does not put extremity exclusivity in its MILP (it only checks it after
solving, §4.1), and RCK's reciprocal-location constraint is not used here.

### 4.4 Linear ends

> **Provenance:** (c) new, `10416aa`. **Code:** `rck/ecdna/ilp_circular.py::OptModelCircular.linear_end_nodes`.

Not every focal amplification is circular. AmpliconArchitect and CoRAL also report linear
amplicons (breakage–fusion–bridge, heavily rearranged or tandem-duplicated), and reach a segment's
end through their source edges [Deshpande et al. 2019](pipeline/references.md#deshpande2019) [Zhu et al. 2024](pipeline/references.md#zhu2024). Therefore, the RCK-ecDNA algorithm adds an addition to RCK that allows the ecDNA layer to start or stop at a set $\mathcal L$ of extremities chosen by `--linear-ends` (both haplotypes of each chosen extremity; only extremities with a segment):

- `boundaries` (default): extremity $v$ of segment $j$ is in $\mathcal L$ when it has no reference
  adjacency, or a reference neighbour has no segment (a region flank), or, for a neighbouring
  segment $j'$, $\lvert o_j-o_{j'}\rvert \ge \max(m_{bd}, \varphi_{bd}\max(o_j,o_{j'}))$ with
  $m_{bd}$ = `--boundary-min-copies` (2) and $\varphi_{bd}$ = `--boundary-frac` (0.25), on the total
  observed copy numbers. That is where the end of a linear amplicon shows in the depth;
- `all`: every extremity of a segment;
- `none`: the strictly circular model.

Each copy-end costs $\kappa_{\mathrm{end}}$ = `--linear-end-penalty` (0.25), so a linear walk costs
0.5 per copy. These values are heuristic and I tested them on a bench of simulated data, but I'm open to tweaking these values. An imputed closing edge costs $\kappa_{\mathrm{imp}}$ = `--imputed-penalty` (2.0) once, so with more than 4 copies
an affordable imputed closure wins (at exactly 4 copies the two cost the same), while a long gap with no other support is explained as linear.

*Why:* without ends, an amplicon with no closing junction can only be explained by forcing a
circle, an imputed join, or deviations, which results in poor results. I first introduced this idea after poor results in an initial run on demo4, and after adding linear ends, was able to deconvolve demo4 correctly, it is in amplicon-allele-analysis `results/rck_ecdna_plan1/` ("Plan 1").

### 4.5 Objective

> **Provenance:** (c) new, `7d8dff9`; `bg_balance_weight` `10416aa`. **Code:**
> `rck/ecdna/ilp_circular.py::OptModelCircular._objective`, `::_seg_weight`.

```math
\begin{aligned}
\min\ & \sum_{j} D_j
+ \omega_{\mathrm{bg}} \sum_j w_j \left\lvert\textstyle\sum_H \mu^{bg}(j^H) - o^{\mathrm{bg}}_j\right\rvert
+ \varepsilon_{\mathrm{bal}} \sum_{j\in\mathcal J_{AB}} \left\lvert\mu^{bg}(j^A)-\mu^{bg}(j^B)\right\rvert
+ \sum_{a:\,o_a\ \text{defined}} \frac{\alpha_a}{\sqrt{\max(o_a,1)}}\left\lvert\mu(a)-o_a\right\rvert\\
&+ \sum_{a\ \text{novel or imputed}}\ \sum_{\varphi\in\Phi(a)} \left(\kappa_a+\kappa_{\mathrm{soft}}\,[a\ \text{novel},\ \varphi\notin\Phi^{\ast}(a)]\right) p_{a,\varphi}
+ \kappa_{\mathrm{end}} \sum_{v\in\mathcal L} s^{ec}_v
+ \sum_b \kappa_{\mathrm{flip},b}\, \chi_b
\ \bigl[+\ \text{depth-bias terms, §4.8}\bigr]
\end{aligned}
```

with $\kappa_a=\kappa_{\mathrm{imp}}$ for an imputed and $\kappa_{\mathrm{nov}}$ for an observed novel
adjacency, $\alpha_a=\alpha_{\mathrm{ref}}$ for a reference and $\alpha_{\mathrm{nov}}$ for a novel
adjacency. $o_a$ is defined for every reference or novel adjacency with a read count (§3.4); imputed
candidates have none. $\mathcal J_{AB}$ is the set of segments with an allele-specific observation (of
the segment, or of at least one of its pieces); the term is dropped when $\varepsilon_{\mathrm{bal}}=0$.
The segment data term is, in this order of precedence,

```math
D_j =
\begin{cases}
\sum_{\pi\in j} w_j\,\sigma_\pi\, n_j\, D_\pi & \text{with pieces (§4.6)}\\[2pt]
w_j\, n_j \sum_{H} \left\lvert\mu(j^H)-o_{j,H}\right\rvert & \text{allele-specific observation}\\[2pt]
w_j\, n_j \left\lvert\mu(j^A)+\mu(j^B)-o_j\right\rvert & \text{total-only observation}
\end{cases}
\qquad
w_j=\max\!\left(\min\!\left(\tfrac{L(j)}{L_{\mathrm{cap}}},1\right), w_{\min}\right),\quad
n_j=\frac{1}{\sqrt{\max(o_j,1)}}
```

where $L(j)$ = end − start + 1 and $n_j$ always uses the total $o_j$. Only $D_j$ is normalised by
$n_j$; the background prior is weighted by $w_j$ alone and the balance tie-break by neither.

The segment data term is adopted from the original RCK algorithm. It was modified in the following ways: \\

1. RCK sums deviations over clones, but RCK-ecDNA compares the sum of its two layers (background and ecDNA), $\mu(j^H) = \mu^{ec}(j^H) + \mu^{bg}(j^H)$, with one bulk observation. \\
2. RCK weights by $L_{10}(j) \approx \lceil L(j)/10\rceil$, which is linear in length with no cap. RCK-ecDNA uses $w_j = \min(L(j)/50\text{ kb}, 1)$ with a floor of 0.02, so long background segments can’t outweigh short amplified ones. This was created in response to initial poor results on demo3 and demo4. \\
3. RCK-ecDNA is now scaling by observation size, $n_j = 1/\sqrt{\max(o_j,1)}$. A deviation counts relative to the observed copy number, as with a Poisson standard deviation. On the other hand, RCK took copy numbers from a caller and didn’t model read noise. \\
4. A total-only case for segments without enough heterozygous SNPs. RCK didn't need this because it always received allele-specific copy numbers. \\

The adjacency read term, $\dfrac{\alpha_a}{\sqrt{\max(o_a,1)}}\lvert \mu(a) - o_a\rvert$	is based on AmpliconArchitect’s read model for breakpoint-graph edges, which fits copy numbers to read counts on discordant and concordant edges (CoRAL keeps the idea, with long reads). AA maximises a Poisson likelihood, but RCK-ecDNA uses an absolute deviation scaled by $1/\sqrt{o_a}$, a linear stand-in for Poisson noise that keeps the model a MILP. Weights $\alpha_{\mathrm{nov}}$ = 0.5 (novel) and $\alpha_{\mathrm{ref}}$ = 0.25 (reference) are hand-set and can be tweaked. \\
The linear copy-end cost, $\kappa_{\mathrm{end}}\sum s^{ec}_v$, is based on AA’s source edges and CoRAL’s source/sink nodes, which let linear structures end. The per-copy cost and its value 0.25 are new (section 4.4). \\
The novel-adjacency cost, $\kappa_{\mathrm{nov}} p_{a,\varphi}$ with $\kappa_{\mathrm{nov}}$ = 0.05, replaces RCK’s constraint that at least $1-\mathrm{fp}$ (90%) of novel adjacencies be used. This is based on the parsimony principle and the value was hand-set. \\
The background prior, $\omega_{\mathrm{bg}}\sum_j w_j \lvert\sum_H \mu^{bg}(j^H) - o^{\mathrm{bg}}_j\rvert$, pulls the background layer towards the flank copy number, so chromosomal copies stay out of the ecDNA layer. It stands in for the whole-genome context RCK has and that my subgraph does not have (section 3.5). $\omega_{\mathrm{bg}} = 0.1$, hand-set. \\
The balance tie-break, $\varepsilon_{\mathrm{bal}}\sum\lvert \mu^{bg}(j^A)-\mu^{bg}(j^B)\rvert$, stops background copies from swapping haplotype at no cost. Added together with the linear ends. \\
The imputed-candidate cost, $\kappa_{\mathrm{imp}}$ = 2.0, makes a candidate adjacency without read evidence win only when it explains several copies (more than 4, against linear ends; §4.4). It encourages parsimony. The value is hand-set. \\
The soft-label cost, $\kappa_{\mathrm{soft}}$ (only with `--soft-junction-hap`), allows a labelling that the reads don’t support, at a cost, instead of forbidding it. Labels are defined in (subgraph_construct.md, section 3.3) and §2.5. This flag is off by default and honestly not used much, I only used it to try exploring other options for one of the demos, but I don't use it for solving for ecDNA and wouldn't recommend using it. Without the flag $\Phi(a)=\Phi^*(a)$, so the term is always 0. \\
The flip cost, $\kappa_{\mathrm{flip}} \chi_b$ ($10^{-6}$), is a very small tie-break for a block with no evidence either way. When flipping and not flipping fit the data equally well, that cost pushes the labelling to not flip so that the block keeps its observed A/B labelling. I added this because for a lot of the unoriented blocks, the data gives no reason either way. The segments may have little allele imbalance, or no novel adjacency with a haplotype label nearby. Then $\chi_b = 0$ and $\chi_b = 1$ have exactly the same objective. Without the cost, the solver’s choice would be arbitrary and the block_flips.tsv report would be misleading because it would list blocks as “flipped” when nothing supported flipping them. \\



| Symbol | Meaning | RCK or new | Default (option) | Code |
|---|---|---|---|---|
| $L(j)$ | segment length (bp), end − start + 1 | RCK ($L(j)$) | | `Segment.length` |
| $L_{\mathrm{cap}}$, $w_{\min}$ | length cap and floor of the segment weight $w_j$ | new | 50,000; 0.02 (no option) | `seg_length_cap`, `min_seg_weight` |
| $w_j$, $n_j$, $D_j$ | segment length weight, observation normaliser, segment data term | new ($D_j$ adapted from RCK's $\lVert C-\tilde C\rVert$) | | `_seg_weight`, `norm` |
| $o_j$, $o_{j,H}$, $o_a$, $o^{\mathrm{bg}}_j$ | observations of §3: segment total, segment allele-specific, adjacency ($\mathrm{NOV}(a)/\theta$ or $\mathrm{REF}(a)/\theta$), flank background | new | | `so.total`, `so.a`/`so.b`, `obs.adjacencies`, `obs.background` |
| $o_\pi$, $o_{\pi,H}$, $\sigma_\pi$ | piece total, piece allele-specific, piece share of the segment weight (§3.6, §4.6) | new | | `pc.total`, `pc.a`/`pc.b`, `pc.share` |
| $\omega_{\mathrm{bg}}$ | pull of the background towards the flank copy number | new | 0.1 (`--bg-prior-weight`) | `bg_prior_weight` |
| $\varepsilon_{\mathrm{bal}}$, $\mathcal J_{AB}$ | tie-break towards a heterozygous background, over segments with allele data | new | 0.01 (no option) | `bg_balance_weight`, `obs.has_alleles` |
| $\alpha_a$: $\alpha_{\mathrm{nov}}$ / $\alpha_{\mathrm{ref}}$ | weight of the read term: novel / reference adjacency | new | 0.5 / 0.25 (`--alpha-novel`, `--alpha-ref`) | `alpha_novel`, `alpha_ref` |
| $\kappa_a$: $\kappa_{\mathrm{imp}}$ / $\kappa_{\mathrm{nov}}$ | cost per labelling used: imputed / observed novel adjacency | new (replaces RCK's utilisation $P$) | 2.0 / 0.05 (`--imputed-penalty`, `--novel-penalty`) | `imputed_penalty`, `novel_penalty` |
| $\kappa_{\mathrm{soft}}$, $\Phi^{\ast}(a)$ | cost of a labelling outside the read-based ones (`--soft-junction-hap` only); the preferred labellings | new | 1.0 (`--soft-junction-penalty`) | `soft_junction_penalty`, `Adjacency.preferred` |
| $\kappa_{\mathrm{end}}$ | cost per linear copy-end | new | 0.25 (`--linear-end-penalty`) | `linear_end_penalty` |
| $\kappa_{\mathrm{flip},b}$ | flip cost (§4.6): unoriented / oriented block | new (RCK's $b_f$ is free) | $10^{-6}$ (no option) / off (`--flip-oriented-penalty`) | `flip_penalty`, `flip_oriented_penalty`, `_flip_cost` |


For `--novel-penalty`: The cost is added per label used: an indicator $p_{a,\varphi}$ for each novel adjacency $a$ and labelling $\varphi$, weighted by $\kappa_{\mathrm{nov}}$ in the objective. There are two reasons for changing this strategy from the original RCK. First, not every novel adjacency belongs to the ecDNA. A union subgraph includes SVs on background chromosomes and in other clones. Forcing 90% of them into the flow like original RCK would push copies through joins that may not be on the amplicon. Secondly, RCK has no adjacency read counts, but RCK-ecDNA has an observed copy number for each novel adjacency ($o_a = \mathrm{NOV}(a)/\theta$), fitted by the read term with weight $\alpha_{\mathrm{nov}} = 0.5$. So the data already says how many copies each one should carry. The penalty only has to stop the model using a novel adjacency the data doesn’t support. The value 0.05 for this flag is a heuristic and was hand-set to be a small number.

Absolute values are linearised with one continuous auxiliary variable each (`_abs_dev`):
$d\ge 0$, $d\ge \text{expr}-\text{target}$, $d\ge\text{target}-\text{expr}$, and $d$ enters the
objective with the term's weight. Gurobi's relative MIP gap is `--mip-gap` ($10^{-4}$).

*Why these terms.* The background layer and its flank prior let the chromosomal copies under an
amplicon (for example 1A + 1B) stay out of the cycles. Without them those copies would be forced
into the ecDNA layer or charged to the deviation term. Capping the length weight at 50 kb stops a
long background segment from dominating short amplified ones. The $1/\sqrt{\max(o,1)}$
normaliser makes a deviation count relative to the size of the observation, as a Poisson standard
deviation would. The read term lets junction reads, not only depth, decide adjacency copy numbers.
The novel-adjacency cost is small parsimony. The imputed cost makes a candidate edge win only when
it explains several copies. `bg_balance_weight` stops background copies from swapping haplotype at no
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

- one binary $\chi_b$ per phase block with orientation 0 (or missing) that has at least one
  allele-specific piece, shared by its pieces across segments; with `--flip-oriented-penalty`, also
  per oriented block with an allele-specific piece;
- an allele-specific piece with observed $(o_{\pi,A},o_{\pi,B})$ fits the targets

  ```math
  \tau_{\pi,A}=o_{\pi,A}+\chi_{b(\pi)}\,(o_{\pi,B}-o_{\pi,A}),\qquad \tau_{\pi,B}=o_{\pi,B}+\chi_{b(\pi)}\,(o_{\pi,A}-o_{\pi,B}),\qquad
  D_\pi=\sum_H \left\lvert\mu(j^H)-\tau_{\pi,H}\right\rvert
  ```

  which are linear in $\chi_b$, so the absolute-deviation constraints stay linear. A piece of a block
  without a flip variable uses $\chi=0$. A total-only piece has
  $D_\pi=\lvert \mu(j^A)+\mu(j^B)-o_\pi\rvert$. Every piece of a segment shares the segment's
  inferred $\mu(j^H)$;
- the segment's weight is split over its pieces by length ($\sigma_\pi$ = piece length / summed
  length of the segment's pieces), so each segment's total weight and the objective scale are
  unchanged. $n_j$ uses the segment total, so it does not depend on $\chi_b$;
- $\kappa_{\mathrm{flip}}=10^{-6}$ for a block none of whose pieces is oriented (`flip_penalty`),
  `--flip-oriented-penalty` for an oriented one. With no evidence either way a block stays as
  observed, and the adjacency labellings decide its orientation through the balance constraints;
- after solving, with every copy number (and every $\beta_g$, §4.8) fixed, $\Delta_b$
  (`delta_if_reversed`) is the change in the objective from reversing $\chi_b$: the change in the
  data terms of the block's allele-specific pieces plus the flip cost. The statuses are:
  `fixed` (oriented by the subgraph and not flipped); `unresolved` (orientation 0, and either no flip
  variable or $\Delta_b\le$ `--flip-tol` (0.01); a flipped block is set back to $\chi_b=0$, kept as
  observed); otherwise `flipped` or `kept` by the final $\chi_b$ (an oriented block flipped under
  `--flip-oriented-penalty` is `flipped`).

The flips are fixed after the flow. The refit (§6), `flow_segments.tsv` and the figure points use
the observations with the final flips applied.

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
| telomeres | excess ($\ge$ balance) at reference telomeres and listed telomere positions (and their reference partners) | ecDNA layer: only at linear ends, each copy-end costs `--linear-end-penalty`; background: at terminals, free | RCK's telomere relaxation, moved to the edges of a cut-out subgraph (terminals, §2.4) for the background, and made a per-copy cost for the ecDNA layer after AA's source edges and CoRAL's source and sink (§4.4) |
| novel adjacency use | at most one labelling each; at least $1-\mathrm{fp}$ of them used ($\mathrm{fp}=1-P$, `--run-nas-fp` 0.1) | at most one labelling each; no requirement; each labelling used costs `--novel-penalty` | Not every novel adjacency in a union subgraph is on the ecDNA, and read-based adjacency copy numbers now decide use (§4.5) |
| missing adjacencies | n/a | imputed candidates (AA/BB), each labelling used costing `--imputed-penalty`; background may not use them | Stage 2 proposes candidate joins where no novel adjacency was observed; each must explain several copies to be used (§4.5) |
| allele flips | per input fragment, free | per phase block with allele data that the subgraph left unoriented, $10^{-6}$; oriented blocks fixed unless `--flip-oriented-penalty` | RCK's fragment flip $b_f$, applied to phase blocks the subgraph could not orient; the tiny cost keeps a block with no evidence as observed (§4.6) |
| objective | length-weighted ($L_{10}(j)$) L1 to input allele copy numbers | length-capped, $1/\sqrt{o_j}$-normalised L1 to allele-specific (or total) segment observations, plus adjacency read term, background prior, balance tie-break and penalties | RCK's length-weighted L1, adapted to raw observations; the adjacency read term follows AA's read model for breakpoint-graph edges; the other terms are new and heuristic (§4.5) |
| bounds | input / `--scnb` (default fill $[0,10]$), following the fragment flip; big-M $M_a$ = largest end bound | $U_j$ from the observations; same big-M rule, halved for a self-loop; background capped at $U^{bg}$ | RCK's bounds and big-M, derived automatically; the background cap is new (§4.2) |
| infinite sites | reciprocal-location constraint (homologous-reciprocal-extremity-exclusivity) always in the MILP; extremity exclusivity only checked after solving | per-haplotype extremity exclusivity, optional (`--enforce-is`); no reciprocal-location constraint | Off by default because chromothripsis-derived ecDNA violates it (§4.3) |

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

- segments are grouped, per chromosome in coordinate order, into windows $g$. A window starts at a
  segment, and the next segment joins it only if it ends within W bp of the window's start
  (end − window start + 1 ≤ W) and the gap after the previous segments is at most W bp; otherwise
  it starts a new window. So a window spans at most W bp, except that a segment longer than W is a
  window of its own. Consecutive windows on the same chromosome are neighbours $\mathcal N$ unless a
  gap over W separates them. Windows are named `<chrom>:<n>`;
- each window $g$ gets a continuous $\beta_g\in[-\beta_{\max},\beta_{\max}]$ ($\beta_{\max}$ =
  `--depth-bias-max`, 0.3);
- every segment and piece target of §4.5–§4.6 in window $g$ ($o_j$, $o_{j,H}$, $o_\pi$,
  $o_{\pi,H}$) becomes target $\times(1-\beta_g)$. This is linear because the target is data. A
  flippable piece's target

  ```math
  \big(o_{\pi,H}+\chi_b(o_{\pi,H'}-o_{\pi,H})\big)(1-\beta_g)
  = o_{\pi,H}(1-\beta_g)+(o_{\pi,H'}-o_{\pi,H})(\chi_b-\varpi_{b,g})
  ```

  uses one continuous
  $\varpi_{b,g}=\chi_b\beta_g\in[-\beta_{\max},\beta_{\max}]$ per block and window, linearised exactly
  (McCormick, binary $\chi_b$, bounded $\beta_g$):
  $-\beta_{\max}\chi_b\le\varpi_{b,g}\le\beta_{\max}\chi_b$,
  $\beta_g-\beta_{\max}(1-\chi_b)\le\varpi_{b,g}\le\beta_g+\beta_{\max}(1-\chi_b)$;
- the objective adds

  ```math
  \kappa_\beta\sum_g|\beta_g|+\kappa_{\mathrm{sm}}\sum_{(g,g')\in\mathcal N}|\beta_g-\beta_{g'}|
  ```

  with $\kappa_\beta$ = `--depth-bias-cost` (1.0; a 10% shift of one window costs 0.1) and
  $\kappa_{\mathrm{sm}}$ = `--depth-bias-smooth` (1.0).

Junction reads, the background prior, the balance tie-break and the linear-end boundaries are
unchanged. The refit (§6, `observation_rows`) and the phase-block report use the corrected targets.
`flow_segments.tsv` gains `depth_bias` and `obs_corrected` ($o_j(1-\beta_g)$), and `summary.json`
gains `depth_bias`, only when the option is on. With the option off the model and every output are
as before.

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
(`FlowSolution.ecdna_edges`, `rck/ecdna/phased.py`). There are three kinds of phased edge:

- segment haplotypes $j^H$ (key `("S", sid, H)`), with $F_{j^H}=\mu^{\mathrm{ec}}(j^H)$ (`xe`); $j^H$ joins
  vertex $j^H_t$ to vertex $j^H_h$;
- labelled adjacencies $a^\varphi$ (key `("J", aid, φ)`), with $F_{a^\varphi}=\mu^{\mathrm{ec}}(a^\varphi)$ (`ye`);
- linear ends at vertex $v\in\mathcal L$ (key `("T", nid, H)`), with $F=s^{\mathrm{ec}}_v$ (`te`). A linear end
  joins $v$ to one **virtual end node** shared by all linear ends (`END_NODE`).

A vertex $v$ is a haplotype-specific extremity, so it already carries its haplotype. Only edges with
$F_e>0$ enter the model: $E$ is the set of edges with positive ecDNA flow, and an edge with no ecDNA
flow has no variables, so no species can use it. There are $K$ = `--max-species` (4) slots;
$W=\max_{e\in E}F_e$; $\bar x$ = `--max-multiplicity` (3); $n_{\mathrm{bit}}=\max(1,\lceil\log_2(\bar x+1)\rceil)$
(2 at the default).

| Symbol | Domain | Meaning | Origin | Code |
|---|---|---|---|---|
| $k$, $K$ | $k=1..K$ | slot index; number of slots | CoRAL ($i$, $k$) | `k`, `K` = `max_species` |
| $F_e$ | $\mathbb Z_{>0}$ | ecDNA flow on phased edge $e$ (input) | CoRAL $\mathrm{CN}_e$, here a fitted flow | `flow[e]` |
| $W$ | | $\max_e F_e$, bound on weights | CoRAL `max_cn` | `W` |
| $z_k$ | $\{0,1\}$ | slot $k$ is used | CoRAL $z_i$ | `u[k]` |
| $w_k$ | $\mathbb Z\cap[0,W]$ | integer copy weight of slot $k$ | CoRAL $w_i$ (continuous there) | `w[k]` |
| $\nu_{keb}$ | $\{0,1\}$ | bit $b$ of the multiplicity of edge $e$ in slot $k$ | new | `beta[k,e,b]` |
| $x_{ke}=\sum_b 2^b\nu_{keb}$ | $\{0,\dots,\bar x\}$ | multiplicity: times species $k$ traverses $e$ | CoRAL $x_{ei}$ (integer 0–10 there) | `mult[k,e]` (expression) |
| $\zeta_{keb}$ | $[0,W]$, continuous | $=w_k\nu_{keb}$ (linearised) | new, for the linearisation | `zeta[k,e,b]` |
| $y_{ke}$ | $\{0,1\}$ | $e$ is in slot $k$'s support | new (CoRAL's $y_1,y_2$ are spanning-tree arcs, a different meaning) | `y[k,e]` |
| $\lambda^{\mathrm{lin}}_k$ | $\{0,1\}$ | slot $k$ is a linear walk | new; plays the role of CoRAL's s–t paths | `linear[k]` |
| $r_e$ | $\ge 0$ | residual on edge $e$ | new; replaces CoRAL's coverage constraints (CoRAL's $r_{pi}$ is a path indicator, not this) | `r[e]`, `res[e]` |
| $\bar x$, $n_{\mathrm{bit}}$ | | multiplicity cap; bits per multiplicity | CoRAL caps at 10 | `max_multiplicity`, `nbits` |
| $\omega_e$ | $>0$ | residual weight per copy (§5.3) | new | `_omega(e)` |
| $\kappa_{\mathrm{sp}}$ | | species penalty | CoRAL (coefficient 1 on $z_i$) | `species_penalty` |
| $\kappa_{\mathrm{mult}}$ | | multiplicity tie-break | new | `multiplicity_penalty` |
| $\xi_{ka}$, $\Xi_{k\upsilon}$, $\Xi_\upsilon$, $N_\upsilon$ | $\{0,1\}$; $N_\upsilon>0$ | read terms (§5.4) | CoRAL's $r_{pi}$, $R_p$, modified | `o`, `z`, `x`, `support` |

Constraints, for every slot $k$:

```math
z_k \le w_k \le W z_k,\qquad
\zeta_{keb}\le W\nu_{keb},\quad \zeta_{keb}\le w_k,\quad \zeta_{keb}\ge w_k-W(1-\nu_{keb})
\qquad\forall e\in E,\ b<n_{\mathrm{bit}}
```
Left most constraint means that a used slot has a positive weight and an unused one has zero. The upper half, $w_k\le W z_k$, is from CoRAL (`ConstraintWZ`, with CoRAL's `max_cn` for $W$). The lower half, $z_k\le w_k$, is added in RCK-ecDNA; CoRAL's exact model has no lower bound on a used cycle's weight (its greedy and post-processing models use $w_i\ge 0.1\,z_i$).
The other 3 constraints write the multiplicity in binary so each product $w_k\nu_{keb}$ of an integer and a binary can be expressed exactly with linear constraints, technique published in https://pubsonline.informs.org/doi/10.1287/opre.23.3.434 [Glover 1975](pipeline/references.md#glover1975). CoRAL keeps $\sum_i w_i x_{ei}$ as a quadratic (bilinear) term, which Gurobi has to solve as a nonconvex problem. I changed it in RCK-ecDNA to use linearising to keep the decomposition a plain MILP, which is usually faster and more reliable to solve to optimality. The cost is integer weights.

```math
x_{ke}\le \bar x,\qquad x_{ke}\le \bar x\, y_{ke},\qquad y_{ke}\le x_{ke},\qquad y_{ke}\le z_k\qquad\forall e\in E,
\qquad\qquad \sum_{e\in E} y_{ke}\ge 2z_k
```
These constraints mark which edges a species uses, and require a used species to use at least two edges. With balance (below), the smallest possible species is a one-segment circle closed by one adjacency. The first row is needed because $n_{\mathrm{bit}}$ bits can encode more than $\bar x$ when $\bar x+1$ is not a power of two. CoRAL has no explicit counterpart: in CoRAL an empty cycle with $z_i=1$ is only discouraged by its objective. Here $z_k\le w_k$ and $\sum_e y_{ke}\ge 2z_k$ make $z_k$ count only real species.

```math
x_{k,e_S(v)} = \sum_{e\ \text{adjacency-like at}\ v} l(e)\, x_{ke}\quad\text{for every vertex } v\ne\text{end node},
\qquad \sum_{e\ \text{linear end}} x_{ke}=2\lambda^{\mathrm{lin}}_k,\quad \lambda^{\mathrm{lin}}_k\le z_k
```
At each vertex $v$, the multiplicity of its segment haplotype is the summed multiplicity of its labelled adjacencies and linear ends (a self-loop counted twice, $l(e)=2$). This is the same as CoRAL’s per-cycle Eulerian constraint, which makes each species an alternating walk (segment, adjacency, segment, …). But in RCK-ecDNA, I apply it at haplotype-specific vertices, so a walk stays on one haplotype through a segment, and can switch only through an AB or BA labelling. \\
A species has either no ends (a circle) or two (one linear walk). This is like the source and sink edges ($s$, $t$) and an s–t path for linear structure in CoRAL. RCK-ecDNA routes all linear ends through one virtual node and counts ends instead. This is the decomposition side of section 4.4’s linear ends, so the flow and the decomposition use the same mechanism.

How the code builds the balance rows (`build`, `bal[k,j]`):

- The rows are written for every vertex that is an endpoint of an edge in $E$, except the virtual end
  node, which has no segment to balance against.
- If the vertex's segment haplotype $e_S(v)$ carries no ecDNA flow, the left side is 0, so every
  adjacency-like edge at $v$ is forced to 0 in every slot.
- If a vertex has no adjacency-like edge in $E$, the right side is 0, so the segment is forced to 0.
- The end-count row and $\lambda^{\mathrm{lin}}_k$ exist only when $E$ contains at least one linear end. With
  none, every species is a circle.
- Two ends of one species may sit on the same vertex (one linear-end edge with $x_{ke}=2$).

```math
w_k\ge w_{k+1},\qquad z_k\ge z_{k+1}\qquad (k=1..K-1)
```
Symmetry breaking for MILP. Used slots come first, in non-increasing weight. Slots of equal weight are not ordered further. \\

and the residual and objective:

```math
r_e \ge \sum_k\sum_b 2^b\zeta_{keb} - F_e,\qquad r_e \ge F_e-\sum_k\sum_b 2^b\zeta_{keb}
\qquad\Bigl(\text{so } r_e\ge\bigl|\textstyle\sum_k w_k x_{ke}-F_e\bigr|\Bigr)
```

```math
\min\ \sum_{e\in E} \omega_e r_e + \kappa_{\mathrm{sp}}\sum_k z_k + \kappa_{\mathrm{mult}}\sum_{k,e} x_{ke} + \Omega_{\mathrm{reads}}
```

where $e_S(v)$ is the segment haplotype at vertex $v$; "adjacency-like" edges are
labelled adjacencies and linear ends, a self-loop counted twice; $\kappa_{\mathrm{sp}}$ = `--species-penalty`
(1.0); $\kappa_{\mathrm{mult}}=10^{-4}$ (`multiplicity_penalty`, not a command-line option), a tie-break towards the walk with fewer
traversals; $\omega_e$ is the residual weight of §5.3; and $\Omega_{\mathrm{reads}}$ is §5.4 (0 without `--paths`). The balance
row makes each slot an alternating multiset (each vertex has as many adjacency traversals as
segment traversals). A slot has no linear end (a circle) or exactly two (one linear walk). The
product $w_k x_{ke}$ is linearised exactly through the binary expansion of $x_{ke}$
[Glover 1975](pipeline/references.md#glover1975). The ordering constraints break the symmetry between slots.
The MILP is always feasible (all slots empty, every $r_e=F_e$); only `--path-fraction` (§5.4) can make it infeasible.
Gurobi runs with `LazyConstraints=1`, `MIPGap` = `--mip-gap` ($10^{-4}$), `TimeLimit` = `--time-limit` (none) and
`--threads`. The model also has an option `distinct_supports` (no two used slots with the same support), used
only by the ambiguity check (§5.8). \\

**Explanation for Residual** \\
For CoRAL, the residual is a one-sided fit, $\sum_i w_i x_{ei} \le \mathrm{CN}_e$ where cycles can not exceed an edge’s copy number. In its exact model it requires the cycles to cover at least 90% of the graph’s length-weighted copy number (p_total_weight), and at least 90% of every discordant edge’s (p_bp_cn). In RCK-ecDNA, it penalises the deviation in both directions, weighted per edge ($\omega_e$, section 5.3). This is because the input is a fitted, balanced integer flow, not noisy raw copy numbers, so it can in principle be reproduced exactly. A two-sided residual aims for that, and measures how far an answer is from it. The 90% thresholds were fixed choices by CoRAL and a hard cap can make a decomposition infeasible when copy numbers are uneven. I chose to use a penalty so that we can always get an answer and lets the objective weigh misfit against parsimony. Also, in RCK-ecDNA, the decomposition MILP gives several alternative decompositions, and RCK-ecDNA then sorts them to choose rank 1 using the residual. \\
**Explanation for Objective** \\
The weighted misfit replaces CoRAL’s covered-fraction term and its 90% constraints. The species penalty is from CoRAL. CoRAL minimises the number of cycles minus the covered fraction, here, I just put both into one weighted sum. The tie-break is new and was added to prefer walks with fewer traversals among equally good answers. The long-read paths are also from CoRAL. \\


**Differences from CoRAL's `minimize_cycles`** (`coral/cycle_decomposition.py::minimize_cycles`,
`coral/models/concrete.py::get_minimize_objective`, `coral/models/constraints.py`, CoRAL 3.1.1):

| | CoRAL 3.1.1 | `rck-ecdna` |
|---|---|---|
| graph | decomposes the copy numbers on AA’s breakpoint graph, with total copy number per edge | decomposes a haplotype-resolved flow, so each species comes out with a haplotype on every segment and adjacency. decomposes the fitted integer flow rather than raw copy numbers because the flow MILP has already reconciled the observations, so the decomposition only has to split a balanced flow. Only edges with positive ecDNA flow are in the model |
| weights $w$ | continuous, in $[0,\texttt{max\_cn}]$ | integer in the MILP in $[0,W]$; fractional after the refit (§6) |
| multiplicity $x$ | integer variable, 0–10; discordant edges also capped by their inferred multiplicity | integer, times species $k$ traverses edge $e$, 0–`--max-multiplicity` (3), written with $n_{\mathrm{bit}}$ binary bits. 3 was set at the start of setting up the simulations because it is the largest cap at 2 bits, the smallest that allows a segment traversed more than twice, like in a tandem duplication inside a circle. CoRAL allows 10 (which would need 4 bits here), but CoRAL solves on a smaller, unphased graph. |
| fit | $\sum_i w_i x_{ei}\le \mathrm{CN}_e$ (quadratic); cover ≥ `p_total_weight` (0.9) of length-weighted copy number and ≥ `p_bp_cn` (0.9) of every discordant edge's | $\lvert\sum_k w_k x_{ke}-F_e\rvert$ penalised; exact linearisation |
| objective | number of cycles minus the covered length-weighted fraction | residual + $\kappa_{\mathrm{sp}}$ × species + tie-break (+ reads) |
| connectivity | spanning-tree / node-order constraints | lazy cuts (§5.2) |
| ends | source/sink edges $s$, $t$ | linear ends at $\mathcal L$, through one virtual end node |

### 5.2 Connectivity

> **Provenance:** (a) subtour-elimination cuts [Dantzig et al. 1954](pipeline/references.md#dantzig1954), added lazily [Gurobi](pipeline/references.md#gurobi).
> **Code:** `rck/ecdna/decompose.py::CycleDecomposer._callback`, `::_components`.

Balance alone allows a slot's support to fall apart into several closed walks. The callback runs at
every new integer solution (`MIPSOL`) and does the following:

1. For each slot $k$, read the support $\{e: y_{ke}>0.5\}$. Skip the slot if it is empty.
2. Split the support into connected components by union-find over the edges' endpoints
   (`_components`). The virtual end node is an ordinary node here, so the two ends of a linear
   walk join its component.
3. If there is one component, the slot is fine. Otherwise, for every **ordered** pair
   $(C_i, C_j)$, $i\ne j$, of components of slot $k$, add the lazy cut

```math
\sum_{f\in\delta(C_i)} y_{kf} \ \ge\ y_{k,e_i}+y_{k,e_j}-1
```

where $C_i$ also denotes the set of vertices touched by its edges, $\delta(C_i)$ is the set of edges of $E$
(all edges with flow, not only the support) with exactly one endpoint in $C_i$, and $e_i\in C_i$,
$e_j\in C_j$ are the first edge of each component in the code's edge order. A self-loop is never in
$\delta(C_i)$. A slot with $c$ components gets $c(c-1)$ cuts. A connected, balanced, alternating
multiset is one alternating closed walk (§7); for a linear slot, the walk through the virtual end
node is one open walk between its two ends.

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
\max\!\left(\min\!\left(\tfrac{L(j)}{L_{\mathrm{cap}}},1\right), w^{\mathrm{dec}}_{\min}\right)\cdot \sigma_{\mathrm{seg}} & e=j^H\ \text{segment haplotype}\\
\omega_{\mathrm{adj}} & \text{otherwise (labelled adjacency or linear end)}
\end{cases}
\qquad\text{and, with }\texttt{--residual-mode relative:}\quad
\omega_e \leftarrow \frac{\omega_e}{\sigma_{\mathrm{res}}\,\max(F_e,1)}
```

where $L_{\mathrm{cap}}$ = 50 kb (`seg_length_cap`), $w^{\mathrm{dec}}_{\min}$ = 0.05 (`min_seg_weight` of
`DecompositionParams`; neither is a command-line option), $\sigma_{\mathrm{seg}}$ = `--seg-residual-scale` (1),
$\omega_{\mathrm{adj}}$ = `--adj-residual-weight` (1), and $\sigma_{\mathrm{res}}$ = `--residual-sigma` (0.1). Linear ends
take the adjacency weight $\omega_{\mathrm{adj}}$. Every $F_e$ in the model is at least 1, so
$\max(F_e,1)=F_e$. The relative mode measures the residual in units of the expected noise. It includes adjacencies because
the flow is balanced, so a dip on a segment comes with dips on its adjacencies. Below
$1/\sigma_{\mathrm{res}}$ = 10 copies the relative weight exceeds the absolute one.

`--rank-residual sum` ranks on residual + structure cost (§5.7) instead of residual first. Lower
weights alone only change what reaches the pool; with the default ranking an exact fit in the pool
still ranks first.

*Why (bench_nested v2, step 2 diagnostic).* By default one copy of residual on a segment of 50 kb or
more costs as much as one species. A 1–3 copy dip from depth bias is then cheaper to cover exactly
with a leftover or a merged species. *Effect:* on the four decomposition failures (`s301`, `s703`,
`s505`, `s403`), no setting tried changed rank 1 or added a better pool answer
(`results/bench_nested/diag/rck_noise/SUMMARY.md`). The options stay opt-in, and defaults and their
outputs are unchanged. \\

Weighting misfit by segment length idea comes from RCK and CoRAL. RCK's paper weights each segment’s deviation by its length $L(j)$; RCK's code weights it by the segment attribute chosen by `--run-segment-length-attr` (default `length_10` = $\lceil (L(j)-1)/10\rceil$) (section 4.1). CoRAL measures fit by length-weighted copy number: its objective and its 90% coverage target count each segment’s copies times its length. \\
The 50 kb cap and the 0.05 floor are added for RCK-ecDNA. The cap is the same one the flow MILP uses ($w_j$, section 4.5), so long background segments don’t outweigh short amplified ones. The floor is not: the flow MILP's floor is 0.02. \\
The adjacency weight is similar to CoRAL’s coverage target of at least 90% of every discordant edge’s copy number (p_bp_cn). In RCK-ecDNA, the fits on breakpoints is judged per edge, not by length, because a novel adjacency has no length. Also, I use a soft weight of 1 instead of a hard 90% threshold, and the option --adj-residual-weight. The same weight applies to reference adjacencies and linear ends. \\
--seg-residual-scale and --adj-residual-weight are used to scale the segment or adjacency weights down. These were added for tuning. \\
--residual-mode relative measures misfit in units of the expected noise. A deviation of 10% of an edge’s copy number costs about one unit, at any copy number. This comes from the idea of weighting residuals by their expected standard deviation, like weighted least squares or a χ²-type fit. It is based on an assumption of proportional noise (10% of the copy number), which is usually what is produced by multiplicative biases like GC content. \\
--rank-residual sum ranks on residual plus structure cost.

### 5.4 Long-read constraints

> **Provenance:** (c) modified from CoRAL's subpath constraints [Zhu et al. 2024](pipeline/references.md#zhu2024)
> (`coral/models/constraints.py::get_shared_subpath_edge_constraints`,
> `::get_addtl_subpath_edge_constraints`), `c5e0b9c`. **Code:**
> `rck/ecdna/decompose.py::read_constraints`, `::CycleDecomposer._read_terms`, `::explained_reads`.

With `--paths`, each long-read path is first reduced to its **junctions**: the joins between
consecutive path steps that are not reference adjacencies of the graph (`rck/ecdna/walk.py::read_paths`;
its support is the `#NOV` column). A path $\upsilon$ becomes a read constraint when

- it has at least `--min-path-junctions` (2) junctions (counted with repeats, per path row);
- its own support is at least `--min-path-support` (1) reads (checked per row, before merging);
- every junction matches at least one labelled non-reference (novel or imputed) adjacency $a^\varphi$ with
  ecDNA flow, joining the same two extremities. Otherwise the path is skipped and counted as off the flow.

Paths with the same junction sequence, in either direction, are merged and their support $N_\upsilon$
summed; the first one seen is kept as the representative. For each distinct junction $a$ of $\upsilon$,
$G_a\subseteq E$ is the set of its labellings with flow. Read constraints are built only when
$\kappa_{\mathrm{path}}>0$ or $\phi_{\mathrm{path}}>0$. Then

```math
\xi_{ka}=\bigvee_{e\in G_a} y_{ke},\qquad \Xi_{k\upsilon}=\bigwedge_{a\in \upsilon} \xi_{ka},\qquad \Xi_\upsilon=\bigvee_k \Xi_{k\upsilon}
```

each written as linear constraints on binaries:

```math
\xi_{ka}\ge y_{ke}\ (e\in G_a),\quad \xi_{ka}\le\sum_{e\in G_a}y_{ke};\qquad
\Xi_{k\upsilon}\le\xi_{ka}\ (a\in \upsilon),\quad \Xi_{k\upsilon}\ge\sum_{a\in \upsilon}\xi_{ka}-(|\upsilon|-1);\qquad
\Xi_\upsilon\ge\Xi_{k\upsilon},\quad \Xi_\upsilon\le\sum_k\Xi_{k\upsilon}
```

where $|\upsilon|$ is the number of distinct junctions of $\upsilon$. When $G_a$ has one edge, $\xi_{ka}$ is that
$y_{ke}$ itself (no new variable). $\xi_{ka}$ is shared by all reads that contain $a$. The objective term and the optional constraint are

```math
\Omega_{\mathrm{reads}} = \kappa_{\mathrm{path}}\left(1-\frac{\sum_\upsilon N_\upsilon \Xi_\upsilon}{\sum_\upsilon N_\upsilon}\right),\qquad
\text{optionally}\quad \sum_\upsilon N_\upsilon \Xi_\upsilon \ \ge\ \phi_{\mathrm{path}}\sum_\upsilon N_\upsilon
```

with $\kappa_{\mathrm{path}}$ = `--path-weight` (0.5) and $\phi_{\mathrm{path}}$ = `--path-fraction`
(0, off). The term is never negative, so it does not change the pool's relative gap. The cost of unexplained reads is capped at 0.5, below one species, so reads choose between decompositions
that fit the flow equally well but never add a species on their own. If no decomposition satisfies
`--path-fraction`, the run logs a warning and solves without it (only in `--alternatives-mode pool`: the retry is triggered by
the pool solve returning no solution; in `nogood` mode it does not happen, see the findings note). $\Xi$ is fixed by the slot
supports, so it does not multiply the solution pool.

Containing every adjacency of a read is necessary, not sufficient, for the read to fit a walk.
After solving, each read is checked against each species' top ordering
(`rck/ecdna/walk.py::path_matches`): its junctions must appear consecutively and in order along the walk, in either
direction (cyclically for a circle). `reads_explained` is the share of the read-constraint support found this way.

**Differences from CoRAL.** CoRAL constrains every edge of a path with its count
($x_{ei}\ge r_{pi}\,\mathrm{count}_e$). In its exact model every path constraint must be satisfied.
Its greedy model adds cycles until a fraction `p_subpaths` (0.9) of the *number* of path
constraints is satisfied (together with the weight target), and its post-processing model must
keep at least the share the greedy solution reached. `rck-ecdna` requires only the path's
junctions (non-reference adjacencies) to be in one species,
weighs paths by read support, and checks order afterwards. `--path-fraction 0.9` is therefore a
share of read support, not of path count.

**Where each part comes from.**

- *Using multi-junction long reads as constraints on the decomposition* is CoRAL's idea: a read
  that spans several novel adjacencies shows they lie on one molecule. Its optional threshold,
  `--path-fraction`, plays the role of CoRAL's `p_subpaths` (0.9).
- *Writing "slot $k$ contains read $\upsilon$" with binaries* ($\xi$, $\Xi$ as ORs and ANDs of other
  binaries) is a standard way to express logical conditions as linear constraints.
- *What a read must match* was loosened. CoRAL ties each path edge to the multiplicity of the
  cycle that explains it. `rck-ecdna` only requires the read's junctions to be in one
  species, and checks their order afterwards against the walks (§7). The decomposition MILP has no
  order, only multiplicities, so order cannot be checked inside it. Adjacency containment is the
  part that can, and order is checked where it exists.
- *Weighting reads by support* (summed over identical or reversed paths) is new. A path seen in 10
  reads is stronger evidence than one seen once, and CoRAL counts paths, not reads.
- *A soft cost, capped below one species* ($\kappa_{\mathrm{path}}$ = 0.5 against $\kappa_{\mathrm{sp}}$ = 1) is
  new. Reads then choose between decompositions that fit the flow equally well, which is where
  they are most useful (§5.7), but cannot add a species on their own. The cap is a design choice,
  not a tuned value. CoRAL instead makes path constraints hard in its exact model.
- *Fixed by the slot supports:* $\Xi_{k\upsilon}$ and $\Xi_\upsilon$ are fully determined by the $y_{ke}$, so they do not
  create extra pool solutions that differ only in these binaries (`c5e0b9c`).

### 5.5 Alternatives

> **Provenance:** (a) Gurobi's solution pool [Gurobi](pipeline/references.md#gurobi) and a standard no-good cut, `3c82929`.
> **Code:** `rck/ecdna/decompose.py::CycleDecomposer.solve`.

- `--alternatives-mode pool` (default): Gurobi `PoolSearchMode=2`, `PoolSolutions =
  max(20, 4 × --num-alternatives)` (40 at the default 10), `PoolGap` = `--pool-gap` (0.10). With `--step-penalty` > 0 the
  pool is widened to an absolute window: `PoolGap` = ∞ ($10^{100}$), `PoolGapAbs` = `--ambiguity-window` (2; or
  $2\kappa_{\mathrm{sp}}$ if that is 0), `PoolSolutions = max(100, 4 × --num-alternatives)` (§5.7). One solve; if it returns no solution the run
  stops with an error (or retries without `--path-fraction`, §5.4).
- `--alternatives-mode nogood`: up to `--num-alternatives` solves of the same model. After each solve, the cut
  $\sum_{(k,e):\,\hat y_{ke}=1}(1-y_{ke})+\sum_{(k,e):\,\hat y_{ke}=0}y_{ke}\ge 1$ over all $K\cdot|E|$ support
  binaries is added, where $\hat y$ is that solve's support. The loop stops when a solve finds no solution or when its
  objective exceeds the first optimum by more than the window $|\mathrm{best}|\cdot$`--pool-gap` (with
  `--step-penalty` > 0, the larger of that and `--ambiguity-window`, or $2\kappa_{\mathrm{sp}}$ if that is 0). A solution that fails the
  re-check below is not kept, but it still gets its cut and uses up one solve.

Every solution is read slot by slot (`_read`): weights and multiplicities are rounded to integers, a
repeated circle is reduced (§5.6), and the solution is discarded if any species is unbalanced, has more
than one connected component, or has a number of ends other than 0 or 2. Solutions are then
deduplicated by the set of their species, each species being (integer weight, edge multiplicities)
after the reduction, so slot order does not matter. All valid solutions are kept (no truncation
here) and sorted by (objective, number of species); pruning, ranking and the cut to
`--num-alternatives` follow in §5.7.

**Where each part comes from.**

- *The solution pool* is a Gurobi feature: with `PoolSearchMode=2` Gurobi searches for the best
  `PoolSolutions` solutions within `PoolGap` of the optimum [Gurobi](pipeline/references.md#gurobi).
- *The no-good cut* is the standard cut that removes one binary solution and nothing else: at least
  one support indicator must change. It is the fallback when the pool is not wanted. It acts on the
  slot supports, so a solution with two equal-weight slots swapped is not cut off and can come back;
  deduplication then drops it.
- *Reporting ranked alternatives, rather than one answer,* is this pipeline's choice. Nested and
  repeated structures often have several decompositions that fit equally well (§5.7, §5.8), and a
  single answer would hide that.
- *Re-checking and deduplicating* is new. Pool solutions are integer-feasible but can differ only
  in slot order or in a repeated circle written two ways, so they are validated and merged before
  ranking (`3c82929`, `f84faff`).

### 5.6 Post-processing: repeated circles and `low_weight`

> **Provenance:** (c) new, `41f5804`. **Code:** `rck/ecdna/decompose.py::CycleDecomposer._read`,
> `::_flag_low_weight`.

- **Repeated circle.** A circular slot (no linear end) whose multiplicities $x_{ke}$ share a greatest
  common divisor $g>1$ goes round the same circle $g$ times. It is reported as that circle, with
  multiplicities $x_{ke}/g$, at weight $g\,w_k$. Linear slots are not reduced.
- **`low_weight`.** With $n>1$ species in a decomposition, species $k$ is flagged `low_weight` when

  ```math
  w_k < \tau_{\mathrm{lw}},\qquad \tau_{\mathrm{lw}} =
  \begin{cases}
  \texttt{--min-species-weight} & \text{if set}\\
  \max\!\left(1,\ \phi_{\mathrm{lw}}\max_{k'} w_{k'}\right) & \text{otherwise}
  \end{cases}
  ```

  with $\phi_{\mathrm{lw}}$ = `--min-species-weight-frac` (0.10) and the maximum over the species of the
  same decomposition. A lone species is never flagged. The rule is applied to the integer weights (after the
  repeated-circle reduction) for the ranking, and applied again to the fractional weights after the refit (§6,
  `rck/ecdna/refit.py::refit_decompositions`); the `status` written out is the post-refit one. A flagged species is kept in the
  output but is usually flow residue, not a molecule. The ranking (§5.7) does not count it as a species
  and does not charge its steps, but its flow still enters the residual, and the total number of species is a
  later tie-break key. The evaluation leaves it out (`rck-ecdna-eval --include-low-weight` scores it).

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
  of the species count of the ranking and out of the scores.
- *The 10% threshold* is heuristic, chosen so that a weight-1 leftover next to a heavy molecule is
  flagged. The comparison is strict and the threshold is at least 1, so a weight-1 species is flagged
  only when the heaviest species has more than 10 copies (next to exactly 10 copies, $\tau_{\mathrm{lw}}=1$
  and it is not flagged). The Task 1 scorer uses the same 10% rule as its default *strict* cut for AA's and CoRAL's
  cycles (§9), so every method is judged by the same rule. It is a known weak point for light
  nested species (§10).

### 5.7 Ranking and the step penalty

> **Provenance:** ranking by significant species (c) new, `f84faff`; step penalty (c) new,
> `90ac98a`. **Code:** `rck/ecdna/decompose.py::ranking_key`, `::structure_cost`, `::pool_keeps`,
> `::prune_pool`, `::decompose_with_pool`.

Write $\Gamma$ for one pool solution: a set of species (slots) $k$ with integer weights $w_k$
(after the repeated-circle rule of §5.6) and multiplicities $x_{ke}$. $\Gamma_{\mathrm{sig}}\subseteq\Gamma$
is its set of species not flagged `low_weight` (§5.6), with the flag decided on $\Gamma$'s own integer
weights. Each pool solution gets

```math
\mathrm{res}(\Gamma)=\sum_e \omega_e\left|\sum_{k\in \Gamma} w_k x_{ke}-F_e\right|,\qquad
\lvert\Gamma_{\mathrm{sig}}\rvert,\qquad
\mathrm{steps}(\Gamma)=\sum_{k\in\Gamma_{\mathrm{sig}}}\ \sum_{j,H} x_{k,j^H}
```

```math
\mathrm{SC}(\Gamma)=\kappa_{\mathrm{sp}}\, \lvert\Gamma_{\mathrm{sig}}\rvert+\lambda\, \mathrm{steps}(\Gamma)
```

where $\omega_e$ is the residual weight of §5.3, $\mathrm{steps}$ counts the **steps** of the
significant species, one per segment traversal (a walk's length in segments, not weighted by copies),
$\kappa_{\mathrm{sp}}$ = `--species-penalty` (1.0; the species penalty of §5.1) and $\lambda$ =
`--step-penalty` (default 0). $\mathrm{SC}$ is the structure cost. With $\mathrm{obj}(\Gamma)$ the
decomposition MILP's objective (§5.1, Gurobi's pool objective value), decompositions are sorted
lexicographically by

```math
\bigl(\mathrm{res},\ \mathrm{SC},\ \lvert\Gamma_{\mathrm{sig}}\rvert,\ \lvert\Gamma\rvert,\ \mathrm{obj}\bigr)\quad\text{(default)},\qquad
\bigl(\mathrm{res}+\mathrm{SC},\ \mathrm{SC},\ \lvert\Gamma_{\mathrm{sig}}\rvert,\ \lvert\Gamma\rvert,\ \mathrm{obj}\bigr)\quad(\texttt{--rank-residual sum})
```

with $\mathrm{res}$ and $\mathrm{SC}$ (and their sum) rounded to 6 decimals and $\mathrm{obj}$ to 9. The
sort is stable, so full ties keep the pool order (pool solutions sorted by objective, then by fewer
species). The list is cut to `--num-alternatives` (10). Rank 1 is the answer.

With $\lambda>0$ the widened pool (§5.5) is first pruned (`prune_pool`, `pool_keeps`). Solution
$\Gamma$ is kept if

```math
\mathrm{obj}(\Gamma)\le \mathrm{obj}^\ast + \lvert\mathrm{obj}^\ast\rvert\,\delta_{\mathrm{pool}}
\quad\text{or}\quad
\lambda\, \mathrm{steps}(\Gamma)+\mathrm{obj}(\Gamma)-\mathrm{obj}^\ast < \lambda\, \mathrm{steps}^\ast
```

where $\ast$ is the best pool solution (the first in pool order: lowest objective, then fewest
species), $\delta_{\mathrm{pool}}$ is `--pool-gap` (0.10), and the first test has a $10^{-6}$
tolerance. That is, it is kept if it is within the usual gap, or if the prior can prefer it. With
$\lambda=0$, or a single pool solution, nothing is pruned.

In order, `decompose_with_pool`:

1. solves the pool (§5.5), validated, deduplicated and sorted by objective, then by fewer species;
2. with `--keep-pool`, records the whole unpruned pool (§5.9);
3. prunes it (only with $\lambda>0$);
4. for each remaining solution, sorts its species by weight, then by total multiplicity (both
   descending), flags `low_weight` (§5.6) and computes the sort key above;
5. sorts, keeps the first `--num-alternatives` (at least 1), and annotates each kept decomposition
   with its orderings (§7) and amplicon call (§5.10);
6. builds the ambiguity report for rank 1 (§5.8), if rank 1 has species and `--ambiguity-window` > 0.

The refit (§6) runs afterwards and does not change the ranking.

**Where each part comes from.** Everything in this section is new to this pipeline.

- *Residual first.* A decomposition that fits the flow better beats a simpler one that fits worse.
  This needs the two-sided residual of §5.1, which gives every pool solution a number. CoRAL's
  one-sided fit with 90% coverage thresholds is pass or fail, so it cannot rank a pool.
- *Then structure cost, then species count.* Among equally good fits, the simpler explanation wins.
  The species count $\kappa_{\mathrm{sp}}\, \lvert\Gamma_{\mathrm{sig}}\rvert$ is the same parsimony as
  CoRAL's cycle count (§5.1), but only significant species count, so weight-1 leftovers (§5.6) do not
  make the true answer lose (demo3, `f84faff`).
- *The step penalty $\lambda\,\mathrm{steps}$* is a minimum-description-length prior on walks: among
  explanations that fit equally well, prefer the one that is shorter to write down [MacKay 2003](pipeline/references.md#mackay2003).
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
- merging saves one species ($\kappa_{\mathrm{sp}}$) and costs the extra steps, so it is undone once
  $\lambda\times$ extra steps $>\kappa_{\mathrm{sp}}$;
- a true internal duplication has the same steps as its split into two species, so
  $\kappa_{\mathrm{sp}}$ keeps it whole.

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

After ranking, one more MILP (§5.1, including its read terms and connectivity cuts, with two
additions) looks for alternatives the data cannot tell apart from the top decomposition $\Gamma^\ast$:

```math
\sum_k w_k x_{ke} = \sum_{k\in \Gamma^\ast} w^\ast_k x^\ast_{ke}\quad\forall e,\qquad
\sum_e \left|y_{ke}-y_{le}\right| \ge z_l\quad\forall k<l
```

Here $w^\ast_k$ are the top decomposition's integer weights (the report is built before the refit
of §6). The first fixes the flow each edge covers, so the fit is identical. The second forbids two
species on the same support, which excludes one walk split into copies. Each $|y_{ke}-y_{le}|$ is a
continuous variable $d_{kle}\in[0,1]$ with $d_{kle}\ge y_{ke}-y_{le}$, $d_{kle}\ge y_{le}-y_{ke}$,
$d_{kle}\le y_{ke}+y_{le}$, $d_{kle}\le 2-y_{ke}-y_{le}$, which is exactly the XOR of the two binaries.
Gurobi's pool (`PoolSearchMode=2`) enumerates up to 50 solutions (`ambiguity_pool`, not a CLI
option) within `--ambiguity-window` (2.0, `PoolGapAbs`; 0 turns the report off) of this MILP's best
objective. Each distinct set of significant structures (the multiplicities of the species not flagged
`low_weight` on the solution's own integer weights, weights ignored) other than the top one's is an
alternative $\Gamma'$, kept at its best objective, with

```math
\mathrm{gap}(\Gamma') = \bigl(\mathrm{obj}(\Gamma')+\lambda\, \mathrm{steps}(\Gamma')\bigr)-\bigl(\mathrm{obj}(\Gamma^\ast)+\lambda\, \mathrm{steps}(\Gamma^\ast)\bigr)
```

where $\mathrm{obj}(\Gamma')$ is its objective in this MILP (the §5.1 objective; the residual term is
the same for every solution because the cover is fixed) and $\mathrm{obj}(\Gamma^\ast)$ the top
decomposition's pool objective. The report (`summary.json` `ambiguity`, also per decomposition
in `decompositions`) gives `ambiguous`, `window`, `n_alternatives`, the top decomposition's
`n_species` (significant) and `steps`, the closest alternative's `gap` (closest = smallest
$|\mathrm{gap}|$), and `decided_by`:

- `species_penalty` if the closest alternative has another number of significant species;
- else `step_penalty` if its steps differ and λ > 0;
- else `tie-break` if its steps differ;
- else `reads or tie-break`.

It also lists the first 10 alternatives (best objective first), each with `objective`, `gap`,
`n_species`, `steps`, `rank` (its rank in the ranked list, if it is listed there) and its integer
`weights`. A negative gap means the ranking preferred fewer species over a lower score. An ambiguous
top decomposition logs a warning. The extra MILP adds $K(K-1)/2$ continuous variables and four rows
per edge, plus one row per slot pair. The report is not recomputed after the refit.

*Why:* demo8's correct answer won by 0.004 through the $10^{-4}$ edge-traversal tie-break, and
nothing in the output said so.

**Where each part comes from.** The report is new; no other method here reports it.

- *Fixing the covered flow* ($\sum_k w_k x_{ke}$ equal to the top decomposition's on every edge)
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

λ acts only in the ranking, and every λ > 0 widens the pool to the same window
(`--ambiguity-window`, or $2\kappa_{\mathrm{sp}}$ if that is 0; §5.5), so the decomposition MILP and
its pool do not depend on λ. The pruning of the widened pool does depend on λ. `--keep-pool` (it
needs λ > 0; the CLI refuses it otherwise) writes `pool/`: every pool solution before pruning
(`pool.tsv`, with what the ranking reads), annotated and refitted like a ranked decomposition
(`decompositions.tsv`, rank = pool index + 1) and, with `--truth-dir`, evaluated (`evaluation.tsv`);
plus, in `pool.json`, the settings, the rank-1 pool index per λ and, for each solution that is rank 1
at some λ in `--pool-lambdas` (default 0 0.02 0.05 0.1 0.2 0.5), its same-fit alternatives (§5.8;
only with `--ambiguity-window` > 0). `python -m rck.ecdna.rerank RUN --lambdas …` (default: the
pool's `--pool-lambdas`) applies the run's own `pool_keeps` and `ranking_key` per λ (`rank_pool`,
relative to the first pool solution) and recomputes `ambiguous`, `decided_by` and `gap` with
`ambiguity_summary`; these are missing for a λ whose rank-1 solution was not rank 1 at any λ of
`--pool-lambdas`. λ = 0 is only emulated (solutions within `--pool-gap` of the best): a λ = 0 run
solves a narrower pool whose ties can come out differently (demo3: 12.77 copies at λ = 0, an
equivalent 12.70 at λ > 0). *Effect:* on the bench, re-ranked composites and ambiguity fields
matched direct runs at 32/32 (case, λ), so each case needs two runs (λ = 0 and one pool run) instead
of one per λ.

**Where it comes from.** New, and an engineering step rather than a modelling one. It relies on one
property of the design in §5.7: because λ acts only in the ranking, the MILP and its pool are the
same for every λ > 0, so a saved pool can be re-ranked for any λ without solving again. Putting the
penalty in the MILP, as the first design did, would have made this impossible. It made the λ sweep
on the bench affordable under a 2-session Gurobi licence.

### 5.10 Amplicon call

> **Provenance:** (c) new, `10416aa`. AmpliconArchitect and CoRAL also classify amplicons
> [Deshpande et al. 2019](pipeline/references.md#deshpande2019) [Zhu et al. 2024](pipeline/references.md#zhu2024); the 80% rule here is this pipeline's own.
> **Code:** `rck/ecdna/decompose.py::classify`.

From the significant species of a decomposition (those not flagged `low_weight`), with
$L_k=\sum_{j,H} x_{k,j^H}\,L(j)$ the walk length in bp (segment lengths summed over the walk's
segment traversals):

```math
\xi_{\mathrm{circ}}=\frac{\sum_{k\ \mathrm{circular},\ w_k\ge c_{\min}} w_k L_k}{\sum_{k\ \mathrm{circular},\ w_k\ge c_{\min}} w_k L_k+\sum_{k\ \mathrm{linear}} w_k L_k}
```

The call is `ecDNA` if $\xi_{\mathrm{circ}}\ge\phi_{\mathrm{ec}}$, `linear` if
$\xi_{\mathrm{circ}}\le 1-\phi_{\mathrm{ec}}$, `mixed` in between, and `none` (no fraction) when the
denominator is 0. Here $c_{\min}$ = `--min-ecdna-copies` (2) and $\phi_{\mathrm{ec}}$ =
`--ecdna-call-fraction` (0.8). It is recomputed after the refit, from the refitted weights and the
recomputed `low_weight` flags. It is reported as `amplicon_call` and `circular_fraction`.

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
tumor2: 6 + 1 against a true 7.44). Each decomposition (every ranked one, and every pool solution
with `--keep-pool`) therefore keeps its structures $x_{ke}$ and gets real weights $w_k\ge 0$ fitted
directly to the observations:

```math
\min_{w\ge 0}\ \eta\sum_k w_k + \sum_{\iota} \kappa_\iota\left|\sum_k X_{\iota k}\, w_k - \tilde o_\iota\right|,\qquad
X_{\iota k}=\sum_{e\in E_\iota} x_{ke}
```

The rows $\iota$ are the flow MILP's own data terms (§4.5) with the background held at the flow
solution; $\kappa_\iota$ is the row's cost, $\tilde o_\iota$ its observation minus the background,
and $E_\iota$ the phased edges whose ecDNA copies add up to that observation. `observation_rows`
builds them as follows. For each segment $j$ (in sorted order), exactly one of the first three
cases applies:

| Row | $\kappa_\iota$ | $\tilde o_\iota$ | $E_\iota$ |
|---|---|---|---|
| each piece $\pi$ of $j$ (with `--block-flip`), allele data: one row per $H$ | $w_j\sigma_\pi n_j$ | $(1-\beta_{g(j)})\,o_{\pi,H}-\mu^{\mathrm{bg}}(j^H)$, with $o_{\pi,A}$ and $o_{\pi,B}$ swapped when $\chi_b=1$ for $\pi$'s block $b$ | $\{j^H\}$ |
| each piece $\pi$ of $j$, total only | $w_j\sigma_\pi n_j$ | $(1-\beta_{g(j)})\,o_\pi-\sum_H\mu^{\mathrm{bg}}(j^H)$ | $\{j^A,j^B\}$ |
| segment haplotype $j^H$ (no pieces), allele data: one row per $H$ | $w_j n_j$ | $(1-\beta_{g(j)})\,o_{j,H}-\mu^{\mathrm{bg}}(j^H)$ | $\{j^H\}$ |
| segment $j$ (no pieces), total only | $w_j n_j$ | $(1-\beta_{g(j)})\,o_j-\sum_H\mu^{\mathrm{bg}}(j^H)$ | $\{j^A,j^B\}$ |
| each adjacency $a$ with an observation $o_a$ (novel or reference with a read count; imputed adjacencies have none) | $\alpha_a/\sqrt{\max(o_a,1)}$ | $o_a-\sum_{\varphi\in\Phi(a)}\mu^{\mathrm{bg}}(a^\varphi)$ | $\{a^\varphi:\varphi\in\Phi(a)\}$ |

Here $w_j=\max(\min(L(j)/L_{\mathrm{cap}},1),w_{\min})$ and $n_j=1/\sqrt{\max(o_j,1)}$ (from the
segment total $o_j$, in every segment and piece row) are as in §4.5; $\mu^{\mathrm{bg}}$ are the flow
solution's background copy numbers; $\chi_b$ are the flow solution's block flips; $\beta_{g(j)}$ is
the flow solution's depth bias for the window $g(j)$ containing $j$ (§4.8; 0 without
`--depth-bias-window`), so segment targets are the bias-corrected observations the flow fitted;
adjacency targets are not bias-corrected, as in the flow. $o_a=\mathrm{NOV}(a)/\theta$ or $\mathrm{REF}(a)/\theta$; $\alpha_a$ is `--alpha-novel` for novel
adjacencies and `--alpha-ref` for reference ones. $\eta=10^{-3}$ (`weight_penalty`, not a CLI
option).

The absolute values are linearised in the usual way: variables $w_k\ge 0$ and $d_\iota\ge 0$, rows
$d_\iota\ge \sum_k X_{\iota k}w_k-\tilde o_\iota$ and $d_\iota\ge \tilde o_\iota-\sum_k X_{\iota k}w_k$,
objective $\eta\sum_k w_k+\sum_\iota\kappa_\iota d_\iota$, solved by `scipy.optimize.linprog` with
`method="highs"`. If the lightest species is refitted below `--refit-min-copies` (0.5), it is dropped
and the LP solved again (the last species is never dropped). If the LP fails, the integer weights are
kept with a warning. The ranking is not changed. Weights are rounded to 2 digits, the species are
re-sorted by weight, `low_weight` flags and the amplicon call are recomputed from them, dropped
species are listed under `pruned`, and `data_objective` / `data_objective_int` report the misfit
(the sum over rows, without the $\eta$ term) with the refitted (unrounded) and the integer weights.
`--no-refit` reports the integer weights. The LP needs no Gurobi licence.

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
extremity has exactly one segment and as many adjacency traversals as segment traversals, so a
connected multiset has such a closed walk. A walk is written as steps; a step is one segment
traversal (segment, haplotype, strand) and the adjacency taken right after it.

- `alternating_circuit`:
  1. start at the first endpoint of the first segment edge in sorted order (or of the pseudo segment
     for a linear species, below);
  2. build a greedy closed trail: traverse the current extremity's segment, then take the first
     adjacency (in sorted order) at the far extremity that still has copies left, and repeat until
     the walk is back at an extremity whose segment has no copies left (this can only be the start;
     anything else is reported as unbalanced);
  3. scan the trail; before any segment traversal whose segment still has copies left, splice in a
     new trail started at the extremity that traversal starts from, and rescan from there. A spliced
     trail starts with a segment step and returns by an adjacency to the same extremity, so splicing
     before a segment step keeps alternation;
  4. if copies are still left, the multiset is not connected (an error).
- A linear species is closed with a pseudo segment between its two ends, extracted as a circuit
  starting on it, and opened again: the pseudo step is dropped, and the walk starts at one end and
  stops at the other.
- `enumerate_circuits` runs a depth-first search from the same start. At each extremity the segment
  step is forced; the search branches over the distinct adjacencies (sorted) at the far extremity
  that still have copies left. A complete walk that returns to the start is stored under its
  canonical form (`canonical`: the smallest token sequence over all rotations of the walk and its
  reverse for a circle, the smaller of the walk and its reverse for a linear walk), so orderings are
  distinct up to rotation and reversal. The search stops at `--max-orderings` (20) distinct orderings
  or 200,000 expansions (adjacency choices tried; `max_expansions`, not a CLI option), and returns
  the orderings sorted by canonical form. If it finds none, the single walk of `alternating_circuit`
  is used.
- Orderings of one multiset have identical copy numbers, so only long reads can tell them apart.
  Each `--paths` path is reduced to its non-reference junctions (novel or imputed adjacencies) and
  compared with the walk's own non-reference junction sequence, at the level of extremities without
haplotype (a junction is the pair of extremities it joins). A path is *applicable* if all its
  junction extremities occur in the walk's junctions, and *matched* if its junctions occur
  consecutively and in order along the walk, in either direction (cyclically for a circle). Orderings
  are ranked by the summed support (`#NOV` of `path.tsv`) of matched paths with ≥ 2 junctions, then by
  the number of matched paths (of any number of junctions); remaining ties keep the canonical order.
  The ≥ 2 here is fixed in the code; `--min-path-junctions` and `--min-path-support` apply only to the
  decomposition's read terms (§5.4). The top ordering gives the species' `length_bp` and is used for
  `reads_explained`.

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
- networkx 2.4 removed `connected_component_subgraphs` (five uses, in `rck/core/graph.py` and
  `rck/utils/adj/*.py`); a local helper restores it (floor raised to networkx ≥ 2.4);
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

**Why these matter for the comparison.** The `rck/version-fix` patches and the `b26bb96`, `60ca24d`
and `18dd969` fixes touch only original RCK (`rck/rck_run.py`, `rck/core/*`, `rck/utils/adj/*`,
`setup.py`), not the bulk route of `rck-ecdna`. The multi-clone route (§11) does run on them: it
builds the patched `OptModelMultiClone` and reads its inputs with RCK's readers (it does not call
`rck_run.py`, and passes the `--scnb` records to `update()` itself, as `60ca24d` does). They matter
because original RCK is also a baseline in the benchmark, run on the same graphs as an
allele-specific karyotype method with no ecDNA species (results/graph_route/NEXT_STEPS.md). The
patches let it run on current software without changing its answers; the `--scnb` fix changes
answers only when that option is used, and should be stated whenever original RCK is run with it.

---

## 9. Evaluation

Scoring against simulator truth is documented in
[docs/pipeline/evaluation.md](pipeline/evaluation.md), with every metric defined:

- `rck-ecdna-eval` and `--truth-dir` (`evaluate.py`): truth reading, scope, per-pair scores, the
  recovery matrix and composite;
- the method-agnostic AA-format scorer `rck.ecdna.aa_eval`, with `aa_export`;
- the graph route (`aa_export subgraph`, `aa_import`, `coral_flank`).

Passing `--truth-dir` or `--truth-cycles` (and optionally `--truth-sample`, `--eval-tol-bp`,
`--min-truth-fraction`, `--all-truth-species`) to `rck-ecdna` writes `evaluation.tsv`,
`evaluation.json`, `truth_cycles.tsv` and `recovery_matrix.tsv`; the multi-clone route also writes
`evaluation_psi.tsv` with `--truth-dir` (§11.3; [evaluation.md §3.9](pipeline/evaluation.md)).

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
  the median `REF(a)` (`#REF`) of the subgraph's own reference adjacencies over $c_0$. Those adjacencies lie
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
only data term is the per-clone allele-specific segment copy numbers $\{\hat c_{i,j}, \check c_{i,j}\}$
of clones $i \in \{1..n\}$ and segments $j$, from an upstream caller (HATCHet; on the simulator, the
truth as a stand-in), and it infers per-clone integer haplotype-specific karyotypes: segment copy
numbers $a_{i,j}$, $b_{i,j}$ and adjacency copy numbers $\mu_i(e)$. Which
copies form ecDNA circles and which form linear chromosomes is decided afterwards, by decomposing
the karyotypes heaviest walk first, as AmpliconArchitect [Deshpande et al. 2019](pipeline/references.md#deshpande2019)
and CoRAL [Zhu et al. 2024](pipeline/references.md#zhu2024) do (§11.2, §11.3).

**Input: the annotated DIAG.** `--rck-input <subgraph>/rck_input/` is written by
amplicon-allele-analysis `subgraph_construct.py --clone-profiles`
([its §5.6](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_CNP/simulations/subgraph_construct.md)),
in RCK's native formats:
- `rck.scnt.tsv`: per-clone copy numbers $\{\hat c_{i,j}, \check c_{i,j}\}$ in the observed phase
  frame, with the bounds as `cnb=`;
- `rck.adj.tsv`: observed novel adjacencies $\tilde A_N$;
- `rck.fragments.tsv`: one fragment $f \in \mathcal F$ per phase block;
- `rck.telomeres.tsv`: the outer flank ends. With the first and last segment extremity of every
  chromosome they form the telomere positions $T(R)$;
- `rck.adjgroups.tsv`: molecule groups $u \in \mathcal U$ from long reads and, with
  `--diag-labeling-groups block`, labelling groups from junction reads (§11.4);
- `rck.scnb.tsv`: the bounds again, for `--scnb` (the route always reads it when present);
- two files RCK does not read: `rck.adj.imputed.tsv` and `clone_proportions.tsv` (the cell
  fractions $u_i$).

`rck/ecdna/multiclone.py` reads them exactly as `rck/rck_run.py --no-pre` does, with RCK's readers,
alignment checks and run defaults (`--run-nas-fp 0.1`, i.e. utilisation $P = 0.9$;
`--run-group-m-default-fp 0.1`; `--run-group-n-default-fp 0.1`, which has no `rck-ecdna` option;
`--run-segment-length-attr length_10`). Bounds missing from both files are filled as `rck_run.py`
fills them (uniform 0 to 10).

**Model** (`rck/core/ilp_gurobi_ecdna.py::OptModelMultiCloneEcDNA`, a thin subclass of
`rck/core/ilp_gurobi.py::OptModelMultiClone`).
- With its options off, it builds RCK's model: the same variables, constraints and objective
  (`test_defaults_build_rcks_model`). On demo5, the route's `rck.scnt.tsv` and `rck.acnt.tsv` are
  byte-identical to `rck_run.py --no-pre`'s.
- The only always-on difference is the empty-novel-adjacency case in §11.1.

**Outputs** (`--outdir`):
- `rck.scnt.tsv` and `rck.acnt.tsv`, as `rck_run.py` writes them;
- `summary.json`: `model: multiclone`, status, objective, flipped fragments, the novel
  adjacencies used per clone (`novel_used`), the copies ending at `--mc-linear-ends` extremities
  (`linear_ends`), θ and the settings;
- `gurobi_multiclone.log`, `gurobi_greedy.log`;
- `decompositions.tsv` (unless `--mc-no-decompose`) and, with `--truth-dir` or `--truth-cycles`,
  the evaluation (§11.3).

The multi-clone route writes no `flow_*.tsv`, `aa/` folder or figures.

```bash
rck-ecdna --model multiclone --rck-input subgraphs/demo5/amplicon_0/rck_input \
    --outdir rck_mc_demo5 --time-limit 60 --threads 2
```

### 11.1 Deviations from `OptModelMultiClone`

Below, $p_{a,\varphi} \in \{0,1\}$ says that labelling $\varphi$ of novel adjacency $a$ is present in
some clone (RCK's `PP` variable, the OR of the per-clone presence binaries), and $\mu_i(a^\varphi)$
is clone $i$'s copy number of the labelled edge $a^\varphi$.

| | RCK `OptModelMultiClone` | `OptModelMultiCloneEcDNA` | Default | Why |
|---|---|---|---|---|
| 1 − fp presence rule | $\frac{1}{\lvert\tilde A_N\rvert}\sum_{a \in \tilde A_N}\sum_\varphi p_{a,\varphi} \ge 1-\mathrm{fp}$ over all novel adjacencies of the graph (the code's form of the paper's utilisation $P = 1-\mathrm{fp}$); $0 \ge 1-\mathrm{fp}$ (infeasible) when there are none | the same mean over observed novel adjacencies only (imputed ones are left out of the sum and the count); skipped when there are none | always (identical to RCK whenever there are observed novel adjacencies and imputed ones are off) | imputed candidates are not SV calls; the empty case is an RCK edge case (an amplicon whose junctions are all imputed, e.g. demo5 `amplicon_2`) |
| imputed adjacencies | n/a | `rck.adj.imputed.tsv` joins the graph, and its adjacencies take part in RCK's other constraints like any novel adjacency; cost `--mc-imputed-penalty` $\cdot \sum_\varphi p_{a,\varphi}$ per imputed adjacency (any labelling, any clone). Unset: the file is not read; 0: read, and free | off (unset) | recovers a circularising junction the SV caller missed (plan 4.2) |
| novel adjacency cost | none | `--mc-novel-penalty` $\cdot \sum_\varphi p_{a,\varphi}$ per observed novel adjacency | 0 (off) | kept as an option (plan decision 3) |
| junction read term | none (RCK has no read data term) | `--mc-junction-weight` $\cdot\,\lvert o_a - \sum_i u_i \sum_\varphi \mu_i(a^\varphi)\rvert / \sqrt{\max(o_a,1)}$ for observed novel adjacencies with a `#NOV` value only (the DIAG has no read count per reference adjacency). $o_a = \mathrm{NOV}(a)/\theta$; θ comes from `--junction-reads-per-copy`, else `junction_reads_per_copy` in the subgraph's `scale.json` (the parent folder of `--rck-input`); $u_i$ from `clone_proportions.tsv` (a clone missing there counts as 0). The run stops with an error when θ or `clone_proportions.tsv` is missing | 0 (off) | the bulk route's `--alpha-novel` term, per clone (plan 4.2). Its units are copies, while RCK's objective is copies × `length_10`, so at equal weight it acts as a tie-break |
| extra ends | excess only at telomere positions $T(R)$, free | also at `--mc-linear-ends boundaries` (both extremities of a reference junction between contiguous segments whose input total copy number $\hat c_{i,j}+\check c_{i,j}$ steps by at least max(`--boundary-min-copies`, `--boundary-frac` × the larger) in any clone) or at `all` segment extremities; telomere positions stay free. The excess at such an extremity, per clone and haplotype, is an integer charged `--mc-linear-end-penalty` (0.25) per copy | `none` (off) | explains an amplicon without a closing junction as linear (HSR-like) instead of charging the deltas (plan 4.2) |
| solver settings | `--run-g-mip-gap 0.015`, `--run-g-threads 4`, `--run-g-time-limit 28800` | `--mip-gap 1e-4`, `--threads` (default `SLURM_CPUS_PER_TASK`, else the CPUs the process may use), `--time-limit` (default none) | — | the `rck-ecdna` Gurobi options; not a model change |

Unchanged from RCK:
- the variables, the balance at every extremity (the paper's Eq. 5,
  $\mu_i(e_S(v)) = \mu_i(e_R(v)) + \sum_{e \in E_N(v)} l(e)\,\mu_i(e)$, with excess allowed only at
  $T(R)$) and the objective
  $\sum_i \sum_j w(j)\,(\lvert a_{i,j} - d_{i,j}\rvert + \lvert b_{i,j} - d'_{i,j}\rvert)$, where
  $\{d, d'\} = \{\hat c, \check c\}$ is chosen per fragment $f \ni j$ by its flip. The paper weights by
  $L(j)$; RCK's code weights by the segment attribute `--mc-segment-length-attr`, default
  `length_10` $= \lceil (L(j)-1)/10 \rceil$ (as `rck_run.py`; `OptModelMultiClone` itself falls back to
  `length_100` when no attribute is passed);
- the fragment flips (one per fragment, shared by the clones);
- the bounds (`--scnb` is now readable, see §8);
- the reference labelling (AA/BB); at most one labelling per novel adjacency, shared by all clones
  ($\sum_\varphi p_{a,\varphi} \le 1$: the paper's "exactly one labelled version", $\lvert h^+_E(a)\rvert = 1$,
  relaxed to at most one, with the 1 − fp rule above); no AB/BA labelling for a self-loop;
- the reciprocal-location rules;
- the molecule, general and labelling groups.

The clone proportions $u_i$ enter only through the junction term, and later the bulk weights and
proportions of the decomposition (§11.3).

### 11.2 Why there is no chromosome/ecDNA layer in the model (risk 1)

The first step-3 version split every clone's copy numbers into a chromosome layer, with RCK's
telomere rule, and an ecDNA layer, a strict circulation (plan, original decision 2). Probing every
optimal solution showed:
- **The flanks pinned linear copies.** The chromosome layer's telomere excess never exceeded the
  flank copy number.
- **But a circulation balances in either layer**, and RCK's objective sees only the totals. So
  ecDNA copies could sit in the chromosome layer as a closed cycle, at the same objective.
  - On the synthetic DIAG, the ecDNA layer ranged from 0 to all of the ecDNA copies.
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
walks $q \in \mathcal W$ shared by the clones, heaviest first. This is AmpliconArchitect's rule
[Deshpande et al. 2019](pipeline/references.md#deshpande2019), run with the one-walk MILP of
CoRAL's greedy mode [Zhu et al. 2024](pipeline/references.md#zhu2024)
(`coral/cycle_decomposition.py::maximize_weights_greedy`, CoRAL 3.1.1).

**Input.** RCK's solution $\mu_i$ for every clone $i$, as copy numbers on phased edges:
- segment haplotypes, $\mu_i(j^A) = a_{i,j}$ and $\mu_i(j^B) = b_{i,j}$;
- adjacency labellings, $\mu_i(a^\varphi)$;
- *linear ends*: the copy-number excess at an extremity $v$, which RCK allows only at telomere
  positions $T(R)$ (and at `--mc-linear-ends` extremities). These use the phased-edge vocabulary of
  §5.

The remaining copy number $R_i(e)$ starts at these values. A round may use only edges with
$R_i(e) > 0$ in some clone, so a path can end only where a karyotype still has excess.

**One round** solves one MILP for a single walk $q$:
- integer multiplicities $x_e \in \{0, \dots, \bar x\}$, $\bar x$ = `--max-multiplicity` (binary
  expansion), with a support binary $y_e$; at least two edges are used;
- alternation and balance at every phased extremity;
- connectivity through lazy cuts;
- no ends (a **cycle**) or exactly two (a **path**, through one virtual end node, as in §5.1);
- integer copies $\psi_{q,i} \ge 0$ per clone, at most $\max_e R_i(e)$, with $\sum_i \psi_{q,i} \ge 1$
  and $\psi_{q,i}\,x_e \le R_i(e)$ on every edge (products linearised exactly).

**Which walk is heaviest** (`--mc-heaviest`):
- `copies` (default), AmpliconArchitect's rule. First maximise $\sum_i \psi_{q,i}$. With that fixed
  (a second solve with $\sum_i \psi_{q,i}$ at least the optimum), maximise CoRAL's length-weighted
  copy number $\sum_i \psi_{q,i} \sum_{j,H} x_{j^H} L(j)$ (in bp, $L(j)$ = end − start + 1; not RCK's
  `length_10`), plus CoRAL's long-read bonus, minus $\epsilon \sum_e x_e$ with
  $\epsilon = 1/(\bar x\,\lvert E_q\rvert + 1)$, where $E_q$ are the edges the round may use. All
  traversals together therefore cost less than one bp-copy, so ties go to the walk with fewer edge
  traversals without overriding a difference in weight.
- `length`, CoRAL's rule: the second objective alone.

The walk is subtracted from every clone, $R_i(e) \leftarrow R_i(e) - \psi_{q,i}\,x_e$, and the next
round starts. A cycle whose multiplicities share a factor $g > 1$ is stored as the cycle walked
once, at $g\,\psi_{q,i}$ copies. The loop stops (`summary.json` `decomposition.stopped`) in any of
these cases:
- nothing is left (`exhausted`; RCK's karyotype is an exact balanced integer flow, so this is the
  normal end);
- the remaining length-weighted copy number
  $\lVert R \rVert_L = \sum_i \sum_{j,H} L(j)\,R_i(j^H)$ is at most `--mc-stop-fraction` of the total
  and every molecule group is satisfied (`stop_fraction`; CoRAL's `p_total_weight`; off by default);
- `--mc-max-walks` walks were taken (`max_walks`);
- the MILP finds no walk (`no_walk`).

**Calls.**
- Cycles are ecDNA candidates. When there is more than one cycle, a cycle lighter than
  `--min-species-weight-frac` of the heaviest cycle is `low_weight`. The amplicon call uses the
  significant cycles only, as in §5.10: `ecDNA` when one of them has a bulk weight of at least
  `--min-ecdna-copies`, else `none` (with paths left out, the call is never `linear` or `mixed`,
  and `--ecdna-call-fraction` and `--min-species-weight` have no effect).
- Paths are chromosomes (`status` `chromosomal`). The evaluation skips them.
- A walk's `weight` is its bulk copies per cell, $\sum_i u_i \psi_{q,i}$ (`weight_int` is empty).

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
| weights | one fractional weight | integer $\psi_{q,i}$ per clone, shared walk |
| ends | source edges | RCK's telomere positions (and `--mc-linear-ends`) |
| multiplicity | free | at most `--max-multiplicity` |
| long reads | ordered subpaths | RCK molecule groups (`rck.adjgroups.tsv`) with at least two novel adjacencies: satisfied when the walk uses all their novel adjacencies (any labelling); bonus $\max(\alpha \lVert R \rVert_L / n_{\text{unsatisfied}}, 1)$ per newly satisfied group, $\alpha$ = `--mc-read-alpha` (0.01); a group with an adjacency that has no copies left cannot be satisfied in that round |
| stop | remaining $\le$ 10 %, or weight below resolution | nothing left (or `--mc-stop-fraction`, `--mc-max-walks`, no walk found) |
| ties | solver | fewer edge traversals |

**Outputs.**
- `decompositions.tsv` has the bulk route's columns (Appendix D), one species per walk in
  extraction order, named `ec1, ec2, …` for cycles and `path1, …` for paths, plus:
  - `kind` (`cycle`/`path`) and `order`;
  - `psi_<clone>`: $\psi_{q,i}$;
  - `proportion_<clone>`: $\psi_{q,i} / \sum_{q' \text{ cycle}} \psi_{q',i}$ (empty for paths);
  - `proportion_bulk`: $\sum_i u_i \psi_{q,i}$ over the same sum for all cycles (empty for paths).
- `summary.json` has a `decomposition` block.
- With `--truth-dir` or `--truth-cycles`: `evaluation.*`, `truth_cycles.tsv` and
  `recovery_matrix.tsv` (structure, cycles only; phase blocks from `rck.fragments.tsv` for the
  within-block scores). With `--truth-dir`, also `evaluation_psi.tsv`, which scores ψ and the
  proportions against the `ecdna_psi` rows (`docs/pipeline/evaluation.md` §3.9).

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
  - Fixed by passing the junction read labels as RCK labelling groups within phase blocks
    (§11.4): demo7 then gives ψ exact 4/4 and 1.0 within phase blocks.

### 11.4 Haplotype labelling from junction reads, within phase blocks (DIAG labelling groups)

*Added 2026-11-13. Input-side only: RCK's model is unchanged.*

**What RCK infers, and what it leaves open.** RCK's haplotype-specific reconstruction infers each
novel adjacency's labelling (AA/AB/BA/BB) from the allele-specific segment copy numbers.
- Every fragment $f \in \mathcal F$ (here, a phase block) has a free A/B flip $b_f$ (RCK's
  `FRAGMENT_ALLELE` variable), because the input caller's allele labels are arbitrary per fragment.
- A labelling is therefore fixed only where the copy numbers and balance force it.
- When a novel adjacency joins two fragments that nothing else ties together, several
  combinations of flips and labellings explain the copy numbers equally well. This is common for
  ecDNA, whose junctions join distant amplified pieces in different phase blocks.

**Why this matters for ecDNA but not for RCK's own goal.**
- RCK reports a copy-number karyotype graph and never claims to resolve molecules, so for its
  purpose these alternatives are equivalent.
- RCK-ecDNA decomposes the karyotype into molecules (§11.3). When species overlap (nested
  ecDNA), the alternatives route copies differently and give different species.
- On demo7, RCK labelled the inner species' junctions sim00005–8 AB/BA, while the truth (and the
  junction reads) say AA, at the same objective, 0. The decomposition then merged the inner species
  into the outer one: ψ exact for 2 of 4 (species, clone) pairs.

This is a limit of the information in allele-specific copy numbers, not a defect in RCK's
algorithm.

**RCK's own remedy: labelling groups.** RCK accepts adjacency groups of type `L` (`agt=L`, with
`agl` naming the end of each listed adjacency) and enforces them in
`OptModelMultiClone.define_constraints_for_adjacency_groups_labeling`: all listed ends lie on one
shared haplotype, in every clone, without fixing which. In the code, one indicator per haplotype
$H$ is the OR of the presence binaries $p_{a,\varphi}$ of every labelling that puts a listed
extremity on $H$, over every novel adjacency at that extremity (not only the listed one), and the
two indicators sum to at most 1. Since $p_{a,\varphi}$ is shared by the clones, so is the choice. A
group whose `agl` list does not match its adjacencies in length is skipped. Plain RCK's default
preprocessing infers such groups itself for short-range adjacencies
(`infer_short_nas_labeling_groups`), which the DIAG turns off by needing `--no-pre`.

amplicon-allele-analysis `subgraph_construct.py --diag-labeling-groups block` (its §5.6) supplies
them from **per-end** junction read labels (`junction_hap_obs_end1/2`): the phased reads at each
junction end, in the observed phase frame of that end's block (simulator file
`obs/adjacency_allele_ends.tsv`; on real data, a long-read haplotagger run per phase block):
- **one group per phase block and label (A or B)**, over the junction ends with that label. A
  junction that crosses haplotypes puts its two ends in different groups. The first version gave
  one label per junction to both ends, which is wrong for such junctions;
- **never across blocks**, because read haplotags are comparable only within a phase block;
- groups that RCK's refinement (`refined_labeling_groups`) would merge across blocks or labels are
  dropped.

The multi-clone route reads them like any adjacency group: RCK's readers, `projected_groups`,
`refined_labeling_groups`, unchanged.

**Deviation status.** None in the model, so there is no new row in §11.1. It is an input
RCK already defines. On the DIAG side it is one more annotation, off by default, documented there.

**Results** (2026-11-13; 60 s limit, 2 threads, under 1 s each; within-block scores from
`docs/pipeline/evaluation.md` §3.10):

| case | labelling groups | cycles ψ (clone1, clone2) | ψ exact | composite | composite within blocks |
|---|---|---|---|---|---|
| demo5 | none | 20 | 1/1 | 1.000 | 1.000 |
| demo5 | 1 (per end; chr17:5) | 20 | 1/1 | 1.000 | 1.000 |
| demo7 (truth 26/13, 14/7) | none | 14/7 (merged walk), 12/6 | 2/4 | 0.747 | 0.747 |
| demo7 | 1 (first version: one label per junction; chr4:171, 8 ends) | 26/13, 14/7 | 4/4 | 0.693 | 1.000 |
| demo7 | **5 (per end: chr4:171 A and B, 172, 186, 187)** | **26/13, 14/7** | **4/4** | **1.000** | **1.000** |

With the groups, demo7's two species are separate cycles with exact per-clone copies, all
junctions, and the right order.
- With the first, one-label version, the absolute composite was 0.693 (phasing accuracy 0.58 and
  0.35): its single group left RCK free to choose the A/B of most blocks, so the species'
  haplotypes differed from the truth's *between* blocks. The within-block score (1.0) showed that
  this was the only error.
- With per-end labels, demo7 has groups in five blocks. Joined by the junctions themselves and by
  the copy numbers, they also fix the relative phase of those blocks, and the absolute score is
  1.0 too. No phase between blocks is given to RCK: every label is in its own block's observed
  frame.
- Without the groups, the within-block score equals the absolute one, 0.747: the merge is a real
  structural error, not a frame choice.
- **Label accuracy** (simulations, exact haplotags): on demo5, demo7 and 10 regenerated
  bench_nested cases, 234 of 234 labelled junction ends are right, including the 58 ends of
  crossing junctions.

**Limits.**
- A group says only "same haplotype". That the A ends and the B ends of a block lie on opposite
  haplotypes is not passed on.
- The constraint is hard, so a wrong per-end label would force a wrong labelling. Real haplotags
  carry switch errors, which the simulator leaves out.
- Junction ends alone in their block, like the far ends of sim00006–8 on demo7, are left free.


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
| `--mc-nas-fp`, `--mc-group-m-fp`, `--mc-segment-length-attr` | `--model multiclone`: RCK's `--run-nas-fp` (0.1), `--run-group-m-default-fp` (0.1), `--run-segment-length-attr` (`length_10`; also `length`, `length_100`, `length_1000`). RCK's `--run-group-n-default-fp` is fixed at 0.1 |
| `--mc-imputed-penalty`, `--mc-novel-penalty`, `--mc-junction-weight`, `--mc-linear-ends`, `--mc-linear-end-penalty` | `--model multiclone`: the ecDNA additions of §11.1, all off by default: unset (imputed adjacencies not read), 0, 0, `none` (`none`, `boundaries` or `all`), and 0.25 per copy once `--mc-linear-ends` is on. They also read `--junction-reads-per-copy`, `--boundary-min-copies` and `--boundary-frac` |
| `--mc-heaviest`, `--mc-max-walks`, `--mc-read-alpha`, `--mc-stop-fraction`, `--mc-no-decompose` | `--model multiclone`: the heaviest-first decomposition (§11.3): `copies` (default) or `length`; at most 50 walks; CoRAL's α 0.01; stop share 0 (use everything up); skip it. `--max-multiplicity`, `--min-ecdna-copies`, `--min-species-weight-frac` and `--max-orderings` apply as in the bulk route; `--min-species-weight` and `--ecdna-call-fraction` are not used |
| `--inputs` | `auto` (subgraph files when `scale.json` exists), `subgraph`, or `windows` |
| `--obs-dir` | sets `--windows-dir OBS/windows`, `--extremities OBS/extremity.tsv`, `--paths OBS/path.tsv` |
| `--windows-dir` | UGP windows: `window.tsv.gz`, `window.depth.npz`, `sample_ids.tsv`, optional `snp.tsv.gz` + `snp.[AB]allele.npz` |
| `--extremities` | `extremity.tsv` with genome-wide `#REF`, for junction reads per copy |
| `--paths` | `path.tsv` long-read paths: read constraints (§5.4) and ordering ranking (§7) |
| `--hatchet-seg` | HATCHet `best.seg.ucn`, for depth per copy without a matched normal (cell lines such as COLO320) |
| `--sample`, `--outdir`, `--name` | sample name; output folder; name in output files (default: subgraph folder name; with `--model multiclone`, the folder that holds `--rck-input`) |
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
- **Depth bias (§4.8, opt-in):** `--depth-bias-window` [0 = off; bp], `--depth-bias-cost` [1.0],
  `--depth-bias-smooth` [1.0], `--depth-bias-max` [0.3].
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
- **Gurobi:** `--threads` [`SLURM_CPUS_PER_TASK`, else the CPUs the process may use],
  `--time-limit` (seconds per MILP) [none], `--mip-gap` [1e-4], `-v` / `--verbose`.

## Appendix D. Outputs

| File | Content |
|---|---|
| `decompositions.tsv` | one row per decomposition, species and ordering: `decomposition_rank`, `objective`, `residual`, `n_species`, `species`, `ordering_rank`, `n_orderings`, `paths_matched`, `paths_applicable`, `n_steps`, `amplicon_call`, `circular_fraction`, `reads_explained`, `data_objective` (observation misfit after the refit), `topology`, `status` (`species` or `low_weight`), `weight` (copies per cell, bulk average, fractional after the refit), `weight_int` (integer weight of the structure search), `length_bp`, `path_support`, `uses_imputed`, `junctions`, `walk` (`chrom:start-end:hap:strand`, comma-separated) and `joins` (adjacency after each step, `kind:aid:phasing`; a linear walk's last join is `end:<extremity>:<hap>`) |
| `flow_segments.tsv` | per segment: observed total and allele copy numbers, inferred ecDNA and background copy numbers for A and B. With `--block-flip`, a segment with pieces reports its total split by the BAF of its piece counts summed after the chosen flips. With `--depth-bias-window`, also `depth_bias` and `obs_corrected` (the observed total × (1 − bias)) |
| `flow_adjacencies.tsv` | per adjacency labelling: read support, observed copy number, ecDNA and background flow; `linear_end` rows for copy-ends |
| `aa/<name>_graph.txt`, `aa/<name>_decomp<k>_cycles.txt` | AmpliconArchitect-format graph and cycles [Deshpande et al. 2019](pipeline/references.md#deshpande2019) (total copy number; cycle lines also carry `Haplotypes=`; linear walks start and end with AA's `0+`/`0-` source) |
| `observations.json` | scales, how each was estimated, and every observation; with `--block-flip`, also `pieces`: each piece's segment, block, span, totals, A / B in its base frame, `fixed`, `flipped` and the A / B after the flip (`a_final`, `b_final`) |
| `block_flips.tsv` | with `--block-flip` only: one row per phase block with a piece: `phase_block`, `orientation` (from the subgraph), `flip`, `status` (`fixed`, `flipped`, `kept` or `unresolved`), `n_pieces`, `n_snps`, `delta_if_reversed` |
| `summary.json` | graph normalisation log, flow status (with `linear_end_copies`), top-level `amplicon_call` and `circular_fraction`, decompositions (with `reads_total`, `reads_explained`, `data_objective`, `data_objective_int`, `pruned`: species the refit dropped, and `steps`), top-level `ambiguity` (§5.8); with `--depth-bias-window`, `depth_bias` (number of windows, min, max, number nonzero); with `--block-flip`, `block_flip` (mode, number of pieces and blocks, counts per status, the flipped blocks) |
| `pool/pool.tsv`, `pool/pool.json`, `pool/decompositions.tsv`, `pool/evaluation.tsv` | with `--keep-pool` only: every pool solution before the λ-dependent pruning (`pool_index`, objective, residual, species, significant species and their steps, integer weights), the same solutions as decompositions (`decomposition_rank` = pool index + 1) and their evaluation; `pool.json` adds the settings, the rank-1 solution per `--pool-lambdas` value and those solutions' same-fit alternatives |
| `gurobi_flow.log`, `gurobi_decompose.log` | solver logs |
| `evaluation.tsv`, `evaluation.json`, `truth_cycles.tsv`, `recovery_matrix.tsv` | only with `--truth-dir` or `--truth-cycles` ([evaluation](pipeline/evaluation.md)) |
| (`python -m rck.ecdna.aa_export`) | not written by a run: the answer, its same-fit alternatives, the graph and allele copy numbers in the AA format the method-agnostic scorer reads ([evaluation](pipeline/evaluation.md)) |
| `--model multiclone` | `rck.scnt.tsv`, `rck.acnt.tsv`, `summary.json`, `decompositions.tsv` (with `kind`, `order`, `psi_<clone>`, `proportion_<clone>`, `proportion_bulk`), `gurobi_multiclone.log`, `gurobi_greedy.log`, and with truth `evaluation_psi.tsv` (§11, §11.3); no `flow_*.tsv`, `aa/`, `observations.json`, `pool/` or figures |
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
| `test_evaluate.py` | no | per-pair scores, truth reading and scoping, within-phase-block phasing ([evaluation](pipeline/evaluation.md)) |
| `test_ilp_circular.py` | yes (small) | single circle (10 copies on A plus a 1+1 background); haplotype crossings; imputed closing edge used when cheap, linear flow when expensive, no flow with `--linear-ends none`; linear ends only at boundaries and flanks; a linear amplicon; `--enforce-is` |
| `test_decompose.py` | yes (small) | repeated segment (one species, multiplicity 2); two species sharing a segment, with every optimal answer reported; disconnected circles; pool and no-good alternatives; orderings; a linear species; a circle plus a linear piece (called ecDNA); end-to-end runs for a circle (`ecDNA`) and an amplicon with no closing junction (`linear`, and `none` with `--linear-ends none`); read constraints (merging a read and its reverse, dropping single-junction, weak and off-flow reads); a read through a>x>b choosing the merged circle (reward and `--path-fraction`); a read never adding a species, and the order post-check; noise tolerance (weights and ranking key without Gurobi; a 2-copy dip on one loop) |
| `test_refit.py` | LP: no; CLI: yes (small) | a fractional weight (10.4); an unsupported leftover dropped; a supported one kept; the `--refit-min-copies` threshold; the last species never dropped; observation rows use the flow background; `--no-refit` |
| `test_subgraph_inputs.py` | observations: no; run: yes (small) | scales, segment and junction copy numbers and flank background from `scale.json`, `CN_*` columns and `cn_profile.tsv`; CN columns win over the scale; `--depth-per-copy` recomputes from `WEIGHT`; `--inputs auto/subgraph/windows`; a full run from the subgraph folder alone |
| `test_recovery.py`, `test_aa_eval.py`, `test_aa_import.py`, `test_coral_flank.py` | no | the scorers and the graph route ([evaluation](pipeline/evaluation.md)) |
| `test_plot.py` | no (figures: matplotlib and cnplot) | axis regions padded and merged; figures skipped when cnplot cannot be imported; a circular species with truth and recovery-matrix pages; a linear species without truth |
| `test_step_penalty.py` | yes (small) | a hub segment with two 3-segment loops: at 2:1 the merged walk wins without `--step-penalty` (ambiguous, decided by the species penalty) and the separate circles with it; at 3:1 tie-break against steps, gaps = λ × extra steps; default ranking unchanged; `--ambiguity-window 0`; pool ranking without Gurobi; `--keep-pool` re-ranks to the direct runs at λ = 0.02–0.5; the CLI refuses it at λ = 0; `rerank` reproduces rank 1 |
| `test_block_flip.py` | inputs: no; flow and run: yes (small) | pieces mapped through phantom contraction; orientation −1 swapped and fixed; a block missing from `phase_blocks.tsv` is unoriented; low-SNP pieces total-only; missing `segment_blocks.tsv` falls back with a warning; a swapped unoriented block is flipped (all ecDNA on A, no imputed edge); `--block-flip off` unchanged; an oriented block never flipped (unless `--flip-oriented-penalty`); a block without evidence `unresolved`; the refit uses flipped observations; the outputs |
| `test_depth_bias.py` | grouping, defaults, outputs, refit: no; solves: yes (small) | `--depth-bias-window` (§4.8): off by default; window grouping and neighbours (also from unsorted input); `flow_segments.tsv` unchanged by default and with `depth_bias` when on; refit targets; the flags; unbiased data fit no bias; off leaves the flow unchanged; a low middle segment absorbed by its window's bias |
| `test_aa_eval_blocks.py`, `test_aa_partition.py`, `test_species_matrix.py`, `test_bench_metrics.py`, `test_aa_export_baseline.py` | no | within-block phased scores of `aa_eval`; the confusion-matrix metrics M1 and M2; the per-cell species matrix; the bench_IW metrics (`bench_metrics`); `aa_export.baseline_graph` for the decomposition-only baselines ([evaluation](pipeline/evaluation.md)) |
| `test_gurobi_env.py` | no | thread default (explicit, `SLURM_CPUS_PER_TASK`, CPU affinity); gurobipy logger kept off the console |
| `test_demo4.py` | yes (full licence) | both demo4 subgraphs: evaluation never errors, each subgraph scored against its own species, per-sample truth weights (`RCK_ECDNA_DEMO4_OBS` enables estimated scales) |
| `test_demo7.py` | yes (full licence) | `--block-flip unoriented` on the rebuilt demo7 union subgraph (`RCK_ECDNA_DEMO7_SUBGRAPH`): chr4:185 and chr4:186 flipped, chr4:171, chr4:172 and chr4:187 fixed, no imputed edge |
| `test_multiclone.py` | yes (small) | `--model multiclone` (§11) on tiny synthetic DIAGs: the default model equals `OptModelMultiClone` (variables, constraints, objective, solution); 1 clone and 2 clones (one labelling and one flip shared by the clones); a flipped phase block; no novel adjacencies (feasible); each off-by-default addition; a DIAG labelling group (`agt=L`) becoming RCK's labelling constraint; the CLI outputs and `--model` checks |
| `test_greedy.py` | yes (small) | the heaviest-first decomposition (§11.3): one cycle in two clones (ψ per clone, then one telomere-to-telomere path per haplotype), two species with different ψ per clone, `copies` vs `length` (the looping path), `decompositions.tsv` columns and the evaluation reader skipping paths |
| `test_demo3.py` | yes (full licence) | `amplicon_0` against the demo3 truth, rank 1: cycle order matches, haplotype Jaccard ≥ 0.95, junction recall ≥ 0.9, copy-number error ≤ 30%; the same with long-read paths (≥ 90% of read support explained) |

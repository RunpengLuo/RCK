# Glossary: code names and published terms

The documentation uses the published term for each concept, mostly RCK's
([Aganezov & Raphael 2020](references.md#aganezov2020)), and keeps the code's name wherever a file,
column or option uses it. This table maps one to the other. "This pipeline" means the term has no
published equivalent and is defined in the doc given.

Abbreviations: **RCK** = [Aganezov & Raphael 2020](references.md#aganezov2020) and its code;
**AA** = AmpliconArchitect, [Deshpande et al. 2019](references.md#deshpande2019);
**CoRAL** = [Zhu et al. 2024](references.md#zhu2024) and its 3.1.1 code;
**HATCHet** = [Zaccaria & Raphael 2020](references.md#zaccaria2020);
**ASCAT** = [Van Loo et al. 2010](references.md#vanloo2010).
Stage docs: **S1** = [simulations/README.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/README.md),
**S2** = [simulations/subgraph_construct.md](https://github.com/RunpengLuo/amplicon-allele-analysis/blob/yubi_dev/simulations/subgraph_construct.md),
**S3** = [docs/ecDNA.md](../ecDNA.md), **E** = [evaluation.md](evaluation.md).

## Genome graph

| Code name | Where it appears | Published term | Source | Note |
|---|---|---|---|---|
| node, `NODE_ID` (`<EID>_+` / `<EID>_-`), breakend `(chrom, pos, ±)`, EID | `nodes.tsv`, `obs/extremity.tsv`, AA graph nodes | **extremity** | RCK | `+` ends the segment to its left, `-` starts the segment to its right |
| `reference` edge, segment edge; AA *sequence edge* | `edges.tsv` `EDGE_TYPE`, AA graph | **segment** | RCK; AA | A stretch of the reference genome between two extremities |
| fragment | RCK `docs/Segments.md` | **fragment** | RCK | A coarser span that segments refine; S1 §S4 |
| `extremity` edge, reference join; AA *concordant edge* | `edges.tsv`, AA graph | **reference adjacency** | RCK; AA | Joins the two extremities of one breakpoint in reference order |
| junction, `variant` edge, `adjacency.tsv` row, `AID`; AA *discordant edge* | `edges.tsv`, `obs/adjacency.tsv`, AA graph | **novel adjacency** | RCK; AA | Code and AA files say "junction" or "discordant edge" |
| self-loop discordant edge | AA graph | **fold-back novel adjacency** | RCK; AA | Joins an extremity to itself; only AA or BB labellings are allowed |
| AA *source edge*; `0+` … `0-` in a cycle | AA graph and cycles files | end of a **path** (linear amplicon) | AA | |
| `imputed_variant`, imputed edge, imputed candidate | `edges.tsv`, `uses_imputed` | candidate novel adjacency | this pipeline (S2 §2.4) | No read support; added where depth changes |
| `imputed_ext_*` | `nodes.tsv` | imputed extremity | this pipeline (S2 §2.4) | Depth breakpoint on a 10 kb grid |
| phantom segment (`seg_hap = no_span`) | `edges.tsv` | — | this pipeline (S2, S3 §2.2) | Segment of ≤ 2 bp (`--phantom-max-bp`) between the two sides of a junction; contracted |
| terminal | `summary.json` | analogous to a **telomere** | RCK | Where the background layer may end; S3 §2.4 |
| linear end | `flow_adjacencies.tsv` `linear_end` | end of a linear amplicon; AA/CoRAL source edge | AA; CoRAL | S3 §4.4 |
| DIAG | S3 | **diploid interval adjacency graph** | RCK | |
| breakpoint (sub)graph | S2, AA files | **breakpoint graph** | AA; CoRAL | |
| region | `summary.json` | — | this pipeline (S3 §2.6) | Segments on one chromosome with gaps ≤ 1 Mb |

## Haplotypes

| Code name | Where it appears | Published term | Source | Note |
|---|---|---|---|---|
| A / B | everywhere | haplotype A / B | RCK | |
| `Phasing`, `PHASING`, `PHASINGS`, `joins` | `adjacency_clone.tsv`, `rck.acnt.tsv`, outputs | **haplotype labelling** (phasing) AA, AB, BA, BB | RCK | `PHASING` keeps the most common labelling; `PHASINGS` keeps all of them |
| `junction_hap` A / B / `mixed` | `edges.tsv` | read-based labelling of a novel adjacency | this pipeline (S2 §3.3) | `mixed` means AB, BA, or several labellings |
| `seg_hap` A / B / AB | `edges.tsv` | haplotype composition of a segment | this pipeline (S2 §3.2) | |
| `B_SHARE` | `verify_subgraph.py` | expected B share of a novel adjacency's reads | this pipeline (S2 §6.2) | |
| phase block, `orientation` (1, −1, 0) | `phase_blocks.tsv` | phase block and its A/B frame against the truth | — (S2 §3.4) | |
| phase-block piece | `segment_blocks.tsv` | — | this pipeline (S2 §3.5, S3 §3.6) | The part of a segment inside one phase block; compare RCK's per-fragment allele flip |
| `XI` | `truth/phasing.tsv` | phase orientation of a SNP (latent) | this pipeline (S1 §S0) | |
| `hap_labels: per_step` / `inherit` | simulator config | — | this pipeline (S1 §S3) | |

## Copy numbers and scales

| Code name | Where it appears | Published term | Source | Note |
|---|---|---|---|---|
| `CN` column of junction rows | `edges.tsv` | **observed adjacency copy number** | RCK concept | reads / θ; an observation, not an inference |
| `observed`, `obs_cn` | `flow_adjacencies.tsv` | observed adjacency copy number | RCK | |
| flow value `ye + yb` | `flow_adjacencies.tsv` | **inferred adjacency copy number** | RCK | Integer solution of the flow MILP (S3 §4) |
| `CN_TOTAL`, `CN_A`, `CN_B` | `edges.tsv`, `segment_blocks.tsv` | observed (allele-specific) **segment copy number**, averaged over cells | RCK; HATCHet | |
| flow segment copy numbers | `flow_segments.tsv` | inferred (allele-specific) segment copy number | RCK | |
| `cn={...}` | `rck.scnt.tsv`, `rck.acnt.tsv` | clone-specific segment / adjacency copy number | RCK | Truth |
| `depth_per_copy` γ | `scale.json` | read depth per copy | — (S2 §5.1, S3 §3.1) | |
| `junction_reads_per_copy` θ, `theta_used`, `theta_fitted` | `scale.json`, `pipeline.tsv` | reads per adjacency copy | — (S2 §5.2) | |
| `#NOV`, `#REF` | `obs/*.tsv`, `WEIGHT` | read support of a novel / reference adjacency | RCK | |
| purity, `tumor_ploidy`, `cell_ploidy`, ploidy class | `scale.json` | tumour purity and ploidy | ASCAT; HATCHet | S2 §5.1 |
| `U`, `V` | `proportions.tsv` | clone proportions / read proportions | HATCHet / MASCoTE | S1 §S5–S6 |
| `PSI` ψ | `proportions.tsv` | ecDNA copies per ecDNA-positive cell of a clone | this pipeline (S1 §S5) | |
| `PHI` φ, `ecdna_burden` | `proportions.tsv`, mixture `Copy_count` | bulk ecDNA copies per cell | this pipeline; AA/CoRAL `Copy_count` | |

## Structures and answers

| Code name | Where it appears | Published term | Source | Note |
|---|---|---|---|---|
| species | `truth/ecdna.tsv`, `decompositions.tsv` | ecDNA species; AA/CoRAL **cycle** (or **path** if linear) | AA; CoRAL; Decoil | One distinct molecule type |
| walk, ordering | `walks.walk`, `walk` column | ordered cycle (AA cycles file) | AA | |
| slot | decomposition MILP | — | this pipeline (S3 §5.1) | One candidate species in the MILP |
| `weight`, `Copy_count` | `decompositions.tsv`, cycles files | species copy number (copies per cell) | AA | |
| `weight_int` | `decompositions.tsv` | integer species weight before the refit | this pipeline (S3 §6) | |
| `low_weight` | `status`, `Status=` | — (species that is not significant) | this pipeline (S3 §5.6) | |
| `topology` circular / linear; amplicon call `ecDNA` / `linear` / `mixed` / `none` | outputs | ecDNA vs linear amplicon | AA; CoRAL concept | The threshold is this pipeline's own (S3 §5.10) |
| path, `PATH`, `PID`, `paths.tsv` | `obs/path.tsv`, subgraph | long-read path; RCK molecule adjacency group (`agt=M`); CoRAL **subpath/path constraint** | RCK; CoRAL | |
| mixture amplicon, truth mixture graph | `truth/aa/mixture_amplicon<N>_*` | **amplicon** (ecSimulator/AA) | ecSimulator; AA | Built by this pipeline (S1 §4.3, E §4.2) |
| `graph_reads_paths`, `graph_plain` | graph route | truth mixture graph with / without synthetic reads and path constraints | this pipeline | |
| step | `steps` column | one segment traversal of a walk | this pipeline (S3 §5.7) | |
| λ, `--step-penalty` | options, bench tables | — | this pipeline (S3 §5.7) | |

## Metrics

| Code name | Where it appears | Published term | Source | Defined in |
|---|---|---|---|---|
| `jaccard_total` / `segment_jaccard` | `evaluation.tsv`, recovery matrix | base-pair Jaccard index | [Jaccard 1912](references.md#jaccard1912) | E §3.1 |
| `jaccard_hap` / `haplotype_jaccard` | same | haplotype-aware Jaccard | this pipeline | E §3.1 |
| `composite` | recovery matrix | — | this pipeline | E §3.7 |
| `composite_h0` | `aa_species_metrics.tsv` | — (haplotype-blind composite) | this pipeline | E §4.4 |
| `bp_accuracy` | `aa_graph_metrics.tsv` | breakpoint-graph accuracy | CoRAL | E §4.3 |
| `error_free` | same | error-free reconstruction | AA | E §4.3 |
| `error_free_merged` | same | — | this pipeline, modified from AA | E §4.3 |
| novel-adjacency FNR / FDR | same | novel-adjacency FNR / FDR | RCK | E §4.3 |
| `cycle_overlap_coral_code`, `cycle_overlap_supp` | `aa_species_metrics.tsv` | CoRAL fragment overlap; cycle-interval overlap | CoRAL code; CoRAL supplement | E §4.4 |
| `cyclic_lcs_coral_code`, `cyclic_lcs_supp` | same | CoRAL normalised longest cycle subsequence | CoRAL code; CoRAL supplement | E §4.4 |
| `length_error` | same | reconstruction length ratio (log2 form) | CoRAL | E §4.4 |
| `cn_explained_top1/3` | `aa_sample_metrics.tsv` | best / top-three copy-number ratio | CoRAL | E §4.5 |
| `subwalks_satisfied` | same | subwalk (path) constraints satisfied | CoRAL | E §4.5 |
| strict / best_of | all `aa_*` tables | — (best_of follows CoRAL's per-truth best match) | this pipeline; CoRAL | E §4.6 |
| `merged`, "carries" | `bench_sweep.tsv`, h2h tables | — | this pipeline | bench_nested README §3.2 |
| set metric | `grow_rules.py` | — | this pipeline | E §6.2 |

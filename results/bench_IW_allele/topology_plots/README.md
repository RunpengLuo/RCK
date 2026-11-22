# Topology figures (make_topology_plots.py)

Simple → hard; the script's docstring says what each file is. rck_bulk_* / rck_mc_* pages (rck-ecdna-plot, RCK bulk and multi-clone routes): haplotype A blue, B orange. coral_* figures: CoRAL 3.1.1's own style, which does not show alleles.

| folder | case |
|---|---|
| `1_simple_S1_s11001` | 1 species, one allele |
| `2_hard_S1_s14003` | 1 species with internal duplications / fold-back, one allele |
| `3_simple_S2_s12003` | 2 distinct species (no shared bases), one allele each |
| `4_simple_S4_s13001` | 4 distinct species, one allele each |
| `5_hard_S2_s15001` | 2 nested species, one allele |
| `6_hard_S2_s15003` | 2 nested species, both alleles amplified (62% of ecDNA bases) |
| `7_hard_S4_s16005` | 4 nested species, both alleles amplified (56%) |
| `8_hard_S4_s16003` | 4 nested species, both alleles amplified (70%); the hardest case, multi-clone's weakest |

## topology_compare.png in the numbered folders (make_compare_plots.py --topology-index)

Added 2026-11-21. `topology_compare.png` as in `compare/`: one genome axis; the truth species, RCK MC's strict species (haplotypes after the scorer's one A/B swap per phase block) and CoRAL's strict cycles (`coral:walks_kept`, decomposition-only); "→ ecdnaN" is the benchmark scorer's pairing (`score/aa_eval/aa_species_metrics.tsv`, strict). CoRAL steps over 5 Mb (background in its linear walks) do not set the axis and are drawn clipped to it (said in the CoRAL heading when it happens: s15001). Each folder's README.txt has the full score rows. RCK MC exact = phased exact cycles (within phase blocks); CoRAL exact = blind. Scores: `bench_IW/bench_iw_intpsi_table.tsv`, `bench_iw_allelic_table.tsv`.

| folder | RCK MC | CoRAL |
|---|---|---|
| `1_simple_S1_s11001` | 1 / 1 exact, J_seg 1.00, J_adj 1.00, 0 missed, 0 false, mean composite 1.000 | 1 / 1 exact, J_seg 1.00, J_adj 1.00, 0 missed, 1 false, CN precision 0.19 |
| `2_hard_S1_s14003` | 1 / 1 exact, J_seg 1.00, J_adj 1.00, 0 missed, 0 false, mean composite 1.000 | 0 / 1 exact, J_seg 1.00, J_adj 1.00, 0 missed, 1 false, CN precision 0.27 |
| `3_simple_S2_s12003` | 2 / 2 exact, J_seg 1.00, J_adj 1.00, 0 missed, 0 false, mean composite 1.000 | 2 / 2 exact, J_seg 1.00, J_adj 1.00, 0 missed, 0 false, CN precision 1.00 |
| `4_simple_S4_s13001` | 3 / 4 exact, J_seg 1.00, J_adj 0.88, 0 missed, 0 false, mean composite 0.933 | 2 / 4 exact, J_seg 0.50, J_adj 0.50, 2 missed, 0 false, CN precision 1.00 |
| `5_hard_S2_s15001` | 2 / 2 exact, J_seg 1.00, J_adj 1.00, 0 missed, 0 false, mean composite 1.000 | 0 / 2 exact, J_seg 0.59, J_adj 0.50, 0 missed, 1 false, CN precision 0.36 |
| `6_hard_S2_s15003` | 1 / 2 exact, J_seg 0.82, J_adj 0.72, 0 missed, 0 false, mean composite 0.790 | 0 / 2 exact, J_seg 0.21, J_adj 0.53, 0 missed, 0 false, CN precision 0.29 |
| `7_hard_S4_s16005` | 3 / 4 exact, J_seg 1.00, J_adj 0.90, 0 missed, 0 false, mean composite 0.937 | 1 / 4 exact, J_seg 0.50, J_adj 0.25, 3 missed, 0 false, CN precision 1.00 |
| `8_hard_S4_s16003` | 2 / 4 exact, J_seg 0.76, J_adj 0.66, 0 missed, 1 false, mean composite 0.787 | 0 / 4 exact, J_seg 0.00, J_adj 0.24, 3 missed, 0 false, CN precision 0.03 |
| `item1_dryrun/sib_opp_shared_s18002` | 2 / 2 exact, J_seg 1.00, J_adj 1.00, 0 missed, 1 false, mean composite 0.851 | 0 / 2 exact, J_seg 0.40, J_adj 0.50, 1 missed, 0 false, CN precision 1.00 |
| `item1_dryrun/sib_mixed_mirror_s19006` | 1 / 2 exact, J_seg 1.00, J_adj 0.83, 0 missed, 0 false, mean composite 0.876 | 1 / 2 exact, J_seg 0.69, J_adj 0.50, 1 missed, 0 false, CN precision 1.00 |

`item1_dryrun/`: the three dry-run candidates' truth figures; all three are finished iw_allelic cases (the same simulations), so s18002 and s19006 got a `topology_compare.png` there; s17004's is `compare/2_sib_opp_mirror_s17004/`. `species_matrix/`: matrix figures only (no topology figure needed).

**RCK MC pages here are within phase blocks too (2026-11-21, `make_compare_plots.py --pages-case iw_intpsi <case> <folder>`):** `rck_mc_rank1_ec<k>.png` pairs and colours as the benchmark table (checked on all 8 against the block TSV, the scorer's rows and `mean_composite`); the run's own pages are `*_absolute.png`. `rck_mc_pages_frame_check.tsv` (audit_page_pairing.py) lists old vs new pair per page. Pair changed on `4_simple_S4_s13001` ec3 (none → ecdna3, 0.733) and `7_hard_S4_s16005` ec3 (ecdna1 → ecdna4) and ec4 (ecdna4 → ecdna1). Same pair, higher composite: s12003 ec2, s13001 ec1/ec2/ec4, s16005 ec1/ec2. Unchanged: s11001, s14003, s15001, s15003, s16003. `rck_mc_recovery_matrix.png` here is still the run's (absolute); the within-block matrix is in `species_matrix/iw_intpsi/<case>/`.

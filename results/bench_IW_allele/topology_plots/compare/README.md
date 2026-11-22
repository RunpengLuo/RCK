# RCK MC vs CoRAL topology comparison (make_compare_plots.py)

Cases chosen 2026-11-21 (Yubi: the six plus the two extras) to show RCK MC on nested ecDNA with both alleles amplified. make_compare_plots.py's docstring says what each file is; each folder's README.txt has the case's scores. Haplotype A blue, B orange; CoRAL has no alleles (grey / CoRAL's own style).

**RCK MC pages are within phase blocks (2026-11-21).** `rck_mc_rank1_ec<k>.png` pairs each inferred species with the truth species the benchmark table's scorer pairs it with, haplotypes relabelled with the scorer's one A/B swap per phase block (walk colours and A/B copy-number tracks); the header says "within phase blocks" and gives the scorer's composite and the case mean (= the table's `mean_composite`). Checked on all 8 cases against `rck_mc_recovery_matrix_block.tsv`, the scorer's per-species rows and the table. The run's own absolute-frame pages are kept as `rck_mc_rank1_ec<k>_absolute.png`. Pairs that changed: s18004 ec1 ecdna1 → ecdna2, ec2 ecdna2 → ecdna1; s18006 the same; hard_S4_s16005 ec3 ecdna1 → ecdna4, ec4 ecdna4 → ecdna1; hard_S4_s16002 ec4 (none → ecdna4); hard_S4_s16006 ec1 ecdna4 → ecdna3, ec2 (none → ecdna1), ec3 (none → ecdna2), ec4 ecdna3 → ecdna4 (`../rck_mc_pages_frame_check.tsv`). Redo: `python make_compare_plots.py --index K --pages`.

| folder | set | why |
|---|---|---|
| `1_sib_opp_shared_s18004` | iw_allelic | outer + nested sibling on the opposite allele, equal clone weights (no clone-wise allele signal in the input), clean. RCK MC: both species exact. CoRAL: one cycle over both, sibling missed. |
| `2_sib_opp_mirror_s17004` | iw_allelic | outer + nested sibling on the opposite allele, mirrored clones (outer in clone 1, sibling in clone 2), clean; the pilot case. RCK MC: both exact, psi per clone exact. CoRAL: 1 of 2 exact, a false species, CN precision 0.35. |
| `3_sib_mixed_mirror_s19012` | iw_allelic | sibling on the opposite allele on some shared intervals and the same allele on others, mirrored clones, clean. RCK MC: both exact. CoRAL: sibling missed. |
| `4_sib_opp_shared_s18006` | iw_allelic | opposite-allele sibling, equal clone weights, reuse stratum (the outer walk reuses an extremity). RCK MC: both exact. CoRAL: none exact, CN precision 0.15. |
| `5_hard_S4_s16005` | iw_intpsi | 4 nested species, both alleles amplified on 56% of ecDNA bases. RCK MC: 3 of 4 exact, J_seg 1.00. CoRAL: one cycle, 3 species missed. |
| `6_hard_S4_s16002` | iw_intpsi | 4 nested species, both alleles amplified on 39%. RCK MC: 3 of 4 exact, no false species. CoRAL: 1 missed, none exact, CN precision 0.30. |
| `7_hard_S4_s16006` | iw_intpsi | (extra) 4 nested species, both alleles on 39%. RCK MC: 2 of 4 exact, J_seg 1.00. CoRAL: 2 missed, none exact. |
| `8_hard_S2_s15003` | iw_intpsi | (extra) 2 nested species, both alleles on 62%. RCK MC: 1 of 2 exact, J_seg 0.82. CoRAL: J_seg 0.21, none exact. |

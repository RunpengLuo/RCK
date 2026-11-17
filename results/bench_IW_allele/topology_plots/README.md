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

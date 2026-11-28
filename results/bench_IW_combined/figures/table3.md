# Table 3: RCK-style evaluation (combined 36-case benchmark, within phase blocks)

Scores of Aganezov & Raphael (2020): haplotype-specific novel-adjacency rates (their Fig. 2A), the Eq. 11 allele-specific copy-number distance (Fig. 2B), and generalized infinite-sites violations per case (Fig. 4; the last three columns).

**single-allele cases** (n = 18 cases)

| method | phasing | adjacency FNR (unlabelled) | labelled adjacency FNR | labelled adjacency FDR | Eq. 11 CN distance (rel., median) | extremity reuse | homologous | homologous-reciprocal |
|---|---|---|---|---|---|---|---|---|
| truth | – | – | – | – | – | 0.56 | 0.00 | 0.00 |
| RCK MC | own | 0.12 | 0.15 | 0.04 | 0.000 | 0.44 | 0.00 | 0.00 |
| RCK MC | BAF | 0.12 | 0.15 | 0.11 | 0.000 | 0.44 | 0.06 | 0.67 |
| RCK MC | oracle | 0.12 | 0.12 | 0.00 | 0.000 | 0.44 | 0.00 | 0.00 |
| CoRAL | BAF | 0.02 | 0.08 | 0.21 | 2.383 | 1.17 | 0.06 | 1.56 |
| CoRAL | oracle | 0.02 | 0.02 | 0.01 | 2.383 | 0.61 | 0.00 | 0.11 |
| CE | BAF | 0.01 | 0.07 | 0.22 | 2.390 | 1.06 | 0.06 | 1.56 |
| CE | oracle | 0.01 | 0.01 | 0.01 | 2.389 | 0.61 | 0.00 | 0.11 |

**both alleles amplified** (n = 18 cases)

| method | phasing | adjacency FNR (unlabelled) | labelled adjacency FNR | labelled adjacency FDR | Eq. 11 CN distance (rel., median) | extremity reuse | homologous | homologous-reciprocal |
|---|---|---|---|---|---|---|---|---|
| truth | – | – | – | – | – | 1.22 | 0.17 | 0.00 |
| RCK MC | own | 0.11 | 0.14 | 0.03 | 0.000 | 0.67 | 0.06 | 0.00 |
| RCK MC | BAF | 0.11 | 0.38 | 0.47 | 0.347 | 2.89 | 0.00 | 3.33 |
| RCK MC | oracle | 0.11 | 0.34 | 0.43 | 0.342 | 3.39 | 0.00 | 3.00 |
| CoRAL | BAF | 0.10 | 0.36 | 0.44 | 0.500 | 3.67 | 0.17 | 3.72 |
| CoRAL | oracle | 0.10 | 0.32 | 0.37 | 0.498 | 3.67 | 0.00 | 2.78 |
| CE | BAF | 0.06 | 0.37 | 0.48 | 0.553 | 3.50 | 0.00 | 3.56 |
| CE | oracle | 0.06 | 0.32 | 0.41 | 0.550 | 3.67 | 0.06 | 2.72 |

Means over cases (Eq. 11: median, as CoRAL / CE's over-counted copies make the mean uninformative). own: the method's haplotypes; BAF / oracle: post-hoc phasing (Fig. 7). The input subgraphs contain no spurious candidate junctions (0 of 9 median per case), so RCK's input-level FPR is undefined here (PLAN.md Part B adds spurious candidates). n = 18 mixed cases (12 sibling, 6 nested): descriptive only.

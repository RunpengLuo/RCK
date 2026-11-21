# Table 3: RCK-style evaluation (bench_IW_intpsi, within phase blocks)

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

**both alleles amplified** (n = 6 cases)

| method | phasing | adjacency FNR (unlabelled) | labelled adjacency FNR | labelled adjacency FDR | Eq. 11 CN distance (rel., median) | extremity reuse | homologous | homologous-reciprocal |
|---|---|---|---|---|---|---|---|---|
| truth | – | – | – | – | – | 1.67 | 0.50 | 0.00 |
| RCK MC | own | 0.17 | 0.21 | 0.05 | 0.000 | 0.67 | 0.17 | 0.00 |
| RCK MC | BAF | 0.17 | 0.43 | 0.46 | 0.306 | 3.67 | 0.00 | 4.17 |
| RCK MC | oracle | 0.17 | 0.44 | 0.48 | 0.305 | 3.67 | 0.00 | 4.00 |
| CoRAL | BAF | 0.12 | 0.40 | 0.41 | 2.651 | 4.67 | 0.50 | 4.33 |
| CoRAL | oracle | 0.12 | 0.41 | 0.39 | 2.651 | 4.00 | 0.00 | 3.33 |
| CE | BAF | 0.14 | 0.41 | 0.40 | 1.587 | 4.00 | 0.00 | 3.50 |
| CE | oracle | 0.14 | 0.41 | 0.39 | 1.585 | 4.00 | 0.17 | 3.17 |

Means over cases (Eq. 11: median, as CoRAL / CE's over-counted copies make the mean uninformative). own: the method's haplotypes; BAF / oracle: post-hoc phasing (Fig. 7). The input subgraphs contain no spurious candidate junctions (0 of 9 median per case), so RCK's input-level FPR is undefined here (PLAN.md Part B adds spurious candidates). n = 6 mixed cases: descriptive only.

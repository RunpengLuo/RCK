# Table 2: allele-specific ecDNA reconstruction (bench_IW_intpsi, within phase blocks)

**single-allele cases** (n = 18 cases)

| method | phasing | allele CN recall | allele CN precision | allele CN error (copies) | ecDNA BAF error | major allele right | phased J_seg | phased J_adj | phased exact-cycle rate |
|---|---|---|---|---|---|---|---|---|---|
| RCK bulk | own | 0.97 ± 0.04 | 0.98 ± 0.03 | 0.7 ± 0.6 | 0.000 ± 0.000 | 1.00 ± 0.00 | 0.98 ± 0.07 | 0.95 ± 0.14 | 0.92 |
| RCK bulk | BAF | 0.97 ± 0.04 | 0.98 ± 0.03 | 0.7 ± 0.6 | 0.000 ± 0.001 | 1.00 ± 0.00 | 0.98 ± 0.07 | 0.84 ± 0.21 | 0.69 |
| RCK bulk | oracle | 0.97 ± 0.04 | 0.98 ± 0.03 | 0.7 ± 0.6 | 0.000 ± 0.000 | 1.00 ± 0.00 | 0.98 ± 0.07 | 0.95 ± 0.14 | 0.92 |
| RCK multi-clone* | own | 1.00 ± 0.01 | 1.00 ± 0.00 | 0.1 ± 0.2 | 0.000 ± 0.000 | 1.00 ± 0.00 | 0.99 ± 0.04 | 0.89 ± 0.14 | 0.75 |
| RCK multi-clone* | BAF | 1.00 ± 0.01 | 1.00 ± 0.00 | 0.1 ± 0.2 | 0.000 ± 0.001 | 1.00 ± 0.00 | 0.99 ± 0.04 | 0.80 ± 0.22 | 0.64 |
| RCK multi-clone* | oracle | 1.00 ± 0.01 | 1.00 ± 0.00 | 0.1 ± 0.2 | 0.000 ± 0.000 | 1.00 ± 0.00 | 0.99 ± 0.04 | 0.89 ± 0.14 | 0.75 |
| CoRAL† | BAF | 0.89 ± 0.15 | 0.48 ± 0.41 | 2.2 ± 1.2 | 0.000 ± 0.001 | 0.98 ± 0.08 | 0.74 ± 0.36 | 0.79 ± 0.23 | 0.44 |
| CoRAL† | oracle | 0.89 ± 0.15 | 0.48 ± 0.41 | 2.2 ± 1.2 | 0.000 ± 0.000 | 0.98 ± 0.08 | 0.74 ± 0.36 | 0.91 ± 0.18 | 0.58 |
| CE† | BAF | 0.99 ± 0.02 | 0.33 ± 0.26 | 9.7 ± 13.5 | 0.000 ± 0.001 | 0.99 ± 0.06 | 0.86 ± 0.19 | 0.82 ± 0.20 | 0.47 |
| CE† | oracle | 0.99 ± 0.02 | 0.33 ± 0.26 | 9.7 ± 13.5 | 0.000 ± 0.000 | 0.99 ± 0.06 | 0.86 ± 0.19 | 0.92 ± 0.15 | 0.61 |

**both alleles amplified** (n = 6 cases)

| method | phasing | allele CN recall | allele CN precision | allele CN error (copies) | ecDNA BAF error | major allele right | phased J_seg | phased J_adj | phased exact-cycle rate |
|---|---|---|---|---|---|---|---|---|---|
| RCK bulk | own | 0.95 ± 0.06 | 0.96 ± 0.05 | 3.0 ± 2.3 | 0.012 ± 0.017 | 0.97 ± 0.08 | 0.71 ± 0.23 | 0.65 ± 0.33 | 0.30 |
| RCK bulk | BAF | 0.82 ± 0.12 | 0.83 ± 0.11 | 12.0 ± 7.1 | 0.131 ± 0.086 | 0.98 ± 0.03 | 0.50 ± 0.25 | 0.38 ± 0.29 | 0.10 |
| RCK bulk | oracle | 0.82 ± 0.12 | 0.83 ± 0.11 | 11.9 ± 7.1 | 0.130 ± 0.084 | 1.00 ± 0.00 | 0.51 ± 0.24 | 0.40 ± 0.28 | 0.10 |
| RCK multi-clone* | own | 0.95 ± 0.12 | 1.00 ± 0.00 | 1.6 ± 3.6 | 0.028 ± 0.065 | 0.94 ± 0.16 | 0.92 ± 0.12 | 0.80 ± 0.15 | 0.65 |
| RCK multi-clone* | BAF | 0.81 ± 0.14 | 0.86 ± 0.08 | 10.8 ± 6.3 | 0.132 ± 0.086 | 0.98 ± 0.03 | 0.68 ± 0.15 | 0.45 ± 0.19 | 0.20 |
| RCK multi-clone* | oracle | 0.81 ± 0.14 | 0.86 ± 0.08 | 10.7 ± 6.4 | 0.130 ± 0.084 | 1.00 ± 0.00 | 0.66 ± 0.16 | 0.47 ± 0.28 | 0.25 |
| CoRAL† | BAF | 0.69 ± 0.22 | 0.45 ± 0.43 | 13.9 ± 6.7 | 0.129 ± 0.087 | 0.90 ± 0.20 | 0.32 ± 0.33 | 0.26 ± 0.30 | 0.10 |
| CoRAL† | oracle | 0.69 ± 0.22 | 0.45 ± 0.43 | 13.9 ± 6.6 | 0.128 ± 0.085 | 0.92 ± 0.20 | 0.32 ± 0.33 | 0.26 ± 0.30 | 0.10 |
| CE† | BAF | 0.76 ± 0.11 | 0.59 ± 0.39 | 13.4 ± 5.5 | 0.130 ± 0.089 | 0.90 ± 0.20 | 0.44 ± 0.28 | 0.30 ± 0.28 | 0.10 |
| CE† | oracle | 0.76 ± 0.11 | 0.59 ± 0.39 | 13.3 ± 5.5 | 0.129 ± 0.086 | 0.91 ± 0.20 | 0.45 ± 0.26 | 0.30 ± 0.28 | 0.10 |

Mean ± SD over cases. own: the method's haplotypes; BAF / oracle: post-hoc phasing (fig. 7 caption). Exact-cycle rate = phased exact cycles / truth species, summed over cases. n = 6 mixed cases is too few for tests; report it descriptively.

* given the true clone profiles (upper bound for the multi-clone route).  † decomposition only, on RCK bulk's graph.

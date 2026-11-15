# Table 1: bench_IW_intpsi, per method variant

**all cases** (n = 24 per method)

| method | ecDNA CN recall | ecDNA CN precision | proportion TV | species F1 (blind) | J_seg (blind) | J_adj (blind) | exact-cycle rate |
|---|---|---|---|---|---|---|---|
| RCK-ecDNA bulk | 0.96 ± 0.04 | 0.98 ± 0.03 | 0.07 ± 0.13 | **0.76 ± 0.38** | 0.92 ± 0.16 | **0.89 ± 0.20** | 0.71 |
| RCK-ecDNA multi-clone* | **0.99 ± 0.06** | **1.00 ± 0.00** | **0.06 ± 0.13** | 0.70 ± 0.28 | **0.97 ± 0.07** | 0.87 ± 0.13 | 0.71 |
| CoRAL† (linear walks kept) | 0.85 ± 0.19 | 0.48 ± 0.40 | 0.23 ± 0.22 | 0.43 ± 0.43 | 0.66 ± 0.38 | 0.80 ± 0.30 | 0.43 |
| CoRAL† (linear walks dropped) | 0.69 ± 0.31 | 0.91 ± 0.26 | 0.26 ± 0.27 | 0.47 ± 0.47 | 0.65 ± 0.38 | 0.64 ± 0.38 | 0.43 |
| CycleExtractor† (linear walks kept) | 0.96 ± 0.10 | 0.40 ± 0.32 | 0.33 ± 0.24 | 0.38 ± 0.30 | 0.79 ± 0.24 | 0.80 ± 0.29 | 0.45 |
| CycleExtractor† (linear walks dropped) | 0.95 ± 0.09 | 0.42 ± 0.33 | 0.26 ± 0.22 | 0.46 ± 0.35 | 0.79 ± 0.24 | 0.80 ± 0.29 | 0.45 |

**simple cases** (n = 12 per method)

| method | ecDNA CN recall | ecDNA CN precision | proportion TV | species F1 (blind) | J_seg (blind) | J_adj (blind) | exact-cycle rate |
|---|---|---|---|---|---|---|---|
| RCK-ecDNA bulk | 0.96 ± 0.04 | 0.98 ± 0.04 | **0.07 ± 0.15** | **0.87 ± 0.29** | **0.99 ± 0.03** | **0.98 ± 0.08** | 0.96 |
| RCK-ecDNA multi-clone* | **1.00 ± 0.01** | **1.00 ± 0.00** | 0.10 ± 0.17 | 0.70 ± 0.25 | 0.98 ± 0.05 | 0.87 ± 0.14 | 0.75 |
| CoRAL† (linear walks kept) | 0.88 ± 0.15 | 0.55 ± 0.47 | 0.15 ± 0.15 | 0.62 ± 0.42 | 0.77 ± 0.35 | 0.96 ± 0.14 | 0.68 |
| CoRAL† (linear walks dropped) | 0.76 ± 0.29 | 0.91 ± 0.29 | 0.20 ± 0.26 | 0.67 ± 0.43 | 0.76 ± 0.36 | 0.75 ± 0.35 | 0.68 |
| CycleExtractor† (linear walks kept) | 0.99 ± 0.03 | 0.32 ± 0.29 | 0.24 ± 0.25 | 0.48 ± 0.24 | 0.85 ± 0.19 | 0.96 ± 0.10 | 0.64 |
| CycleExtractor† (linear walks dropped) | 0.98 ± 0.03 | 0.33 ± 0.28 | 0.16 ± 0.19 | 0.56 ± 0.25 | 0.85 ± 0.19 | 0.96 ± 0.10 | 0.64 |

**nested cases** (n = 12 per method)

| method | ecDNA CN recall | ecDNA CN precision | proportion TV | species F1 (blind) | J_seg (blind) | J_adj (blind) | exact-cycle rate |
|---|---|---|---|---|---|---|---|
| RCK-ecDNA bulk | 0.97 ± 0.05 | 0.97 ± 0.03 | 0.07 ± 0.12 | 0.65 ± 0.44 | 0.84 ± 0.20 | 0.81 ± 0.25 | 0.46 |
| RCK-ecDNA multi-clone* | **0.97 ± 0.08** | **1.00 ± 0.00** | **0.01 ± 0.03** | **0.70 ± 0.32** | **0.96 ± 0.08** | **0.87 ± 0.14** | 0.68 |
| CoRAL† (linear walks kept) | 0.83 ± 0.22 | 0.40 ± 0.33 | 0.31 ± 0.26 | 0.23 ± 0.36 | 0.55 ± 0.38 | 0.64 ± 0.33 | 0.18 |
| CoRAL† (linear walks dropped) | 0.62 ± 0.32 | 0.92 ± 0.24 | 0.31 ± 0.28 | 0.28 ± 0.45 | 0.54 ± 0.38 | 0.52 ± 0.39 | 0.18 |
| CycleExtractor† (linear walks kept) | 0.93 ± 0.13 | 0.48 ± 0.34 | 0.43 ± 0.21 | 0.29 ± 0.34 | 0.74 ± 0.27 | 0.65 ± 0.33 | 0.25 |
| CycleExtractor† (linear walks dropped) | 0.92 ± 0.13 | 0.50 ± 0.36 | 0.36 ± 0.20 | 0.35 ± 0.41 | 0.73 ± 0.27 | 0.64 ± 0.33 | 0.25 |

Mean ± SD over cases; best mean per column in bold. * given the true clone profiles (upper bound for the multi-clone route).  † decomposition only, on RCK bulk's graph.
Strict cycles: RCK `status`, CoRAL / CE `cn-frac:0.1`. Exact-cycle rate = exact cycles / true species, summed over cases.

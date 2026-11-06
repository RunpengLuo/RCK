# References

Every reference cited in the pipeline documentation. Each entry was checked on 2026-11-06
against Crossref (`api.crossref.org/works/<DOI>`) or Europe PMC (which also gives the PMID), and
each software entry against its repository. Entries that could only be checked in part say so.
References that could not be verified are not listed here. Where the docs need one, they show
`[citation needed: …]`; the [last section](#not-verified) lists those.

Code is cited as `path/file.py::function`, at these commits:

| Code | Repository | Commit |
|---|---|---|
| RCK-ecDNA, scorer, benchmarks | [RunpengLuo/RCK](https://github.com/RunpengLuo/RCK/tree/yubi_dev), branch `yubi_dev` | b4ff12d |
| Simulator, subgraph construction | [RunpengLuo/amplicon-allele-analysis](https://github.com/RunpengLuo/amplicon-allele-analysis/tree/yubi_dev), branch `yubi_dev` | eea71aa |
| Original RCK | [aganezov/RCK](https://github.com/aganezov/RCK/tree/4c6c288bfe1e20905069d842bc58609b10568c7e) `master` | 4c6c288 (2020-03-11) |
| Original RCK, ported to current Python/Gurobi | [RunpengLuo/RCK](https://github.com/RunpengLuo/RCK/tree/rck/version-fix) `rck/version-fix` | 0b41976 |
| CoRAL | [AmpliconSuite/CoRAL](https://github.com/AmpliconSuite/CoRAL/tree/v3.1.1) tag `v3.1.1`, as installed in `/n/fs/yubi-iw/conda/envs/coral-3.1.1` | 7f3bf31 |

## Methods and tools this pipeline builds on or compares against

<a id="aganezov2020"></a>**[Aganezov & Raphael 2020] RCK.**
Aganezov S, Raphael BJ. Reconstruction of clone- and haplotype-specific cancer genome karyotypes
from bulk tumor samples. *Genome Research* 30(9):1274–1290 (2020).
doi:[10.1101/gr.256701.119](https://doi.org/10.1101/gr.256701.119), PMID 32887685.
Preprint: bioRxiv (2019), doi:[10.1101/560839](https://doi.org/10.1101/560839).

<a id="rck_upstream"></a>**[RCK code]** https://github.com/aganezov/RCK at 4c6c288. Docs used
for terminology: `README.md`, `docs/Segments.md`, `docs/Adjacencies.md`, `docs/AdjacencyGroups.md`.

<a id="deshpande2019"></a>**[Deshpande et al. 2019] AmpliconArchitect (AA).**
Deshpande V, Luebeck J, Nguyen ND, et al. Exploring the landscape of focal amplifications in
cancer using AmpliconArchitect. *Nature Communications* 10(1):392 (2019).
doi:[10.1038/s41467-018-08200-y](https://doi.org/10.1038/s41467-018-08200-y), PMID 30674876.
Code: https://github.com/AmpliconSuite/AmpliconArchitect.

<a id="zhu2024"></a>**[Zhu et al. 2024] CoRAL.**
Zhu K, Jones MG, Luebeck J, et al. CoRAL accurately resolves extrachromosomal DNA genome
structures with long-read sequencing. *Genome Research* 34(9):1344–1354 (2024).
doi:[10.1101/gr.279131.124](https://doi.org/10.1101/gr.279131.124), PMID 38981681.

<a id="coral_code"></a>**[CoRAL 3.1.1 code]** https://github.com/AmpliconSuite/CoRAL, tag
`v3.1.1` (commit 7f3bf31). The docs cite functions in the installed package
(`coral/<module>.py::<function>`).

<a id="cycleextractor"></a>**[Faizrahnemoon et al. 2026] CycleExtractor.**
Faizrahnemoon M, Luebeck J, Hung KL, et al. Fast and accurate resolution of ecDNA sequence using
Cycle-Extractor. bioRxiv preprint, posted 2026-03-13.
doi:[10.64898/2026.03.10.710955](https://doi.org/10.64898/2026.03.10.710955).
Code: https://github.com/AmpliconSuite/CycleExtractor (local copy
`/n/fs/yubi-iw/tools/CycleExtractor`, commit 1b29b33).

<a id="giurgiu2024"></a>**[Giurgiu et al. 2024] Decoil.**
Giurgiu M, Wittstruck N, Rodriguez-Fos E, et al. Reconstructing extrachromosomal DNA structural
heterogeneity from long-read sequencing data using Decoil. *Genome Research* 34(9):1355–1364
(2024). doi:[10.1101/gr.279123.124](https://doi.org/10.1101/gr.279123.124), PMID 39111816.

<a id="zaccaria2020"></a>**[Zaccaria & Raphael 2020] HATCHet, and the MASCoTE simulation framework.**
Zaccaria S, Raphael BJ. Accurate quantification of copy-number aberrations and whole-genome
duplications in multi-sample tumor sequencing data. *Nature Communications* 11(1):4301 (2020).
doi:[10.1038/s41467-020-17967-y](https://doi.org/10.1038/s41467-020-17967-y), PMID 32879317.

<a id="mascote"></a>**[MASCoTE code]** https://github.com/raphael-group/mascote (default branch
`master`). Not vendored in this project. The function names the simulator docs attribute to it
(`RandomTree.buildRandom`, `Mutation.mutate`, `MixBAMs.computeLengths`, …) were not checked
against the code.

<a id="ecsimulator"></a>**[ecSimulator]** https://github.com/AmpliconSuite/ecSimulator (README: version
0.7.1). Vendored at amplicon-allele-analysis `third_party/ecSimulator` (pinned commit 07a0e3c). No
separate paper was found; ecSimulator is the simulator used in the CoRAL benchmarks
([Zhu et al. 2024](#zhu2024)).

<a id="vanloo2010"></a>**[Van Loo et al. 2010] ASCAT.**
Van Loo P, Nordgard SH, Lingjærde OC, et al. Allele-specific copy number analysis of tumors.
*PNAS* 107(39):16910–16915 (2010).
doi:[10.1073/pnas.1009843107](https://doi.org/10.1073/pnas.1009843107), PMID 20837533.

<a id="gridss2"></a>**[Cameron et al. 2021] GRIDSS2.**
Cameron DL, Baber J, Shale C, et al. GRIDSS2: comprehensive characterisation of somatic structural
variation using single breakend variants and structural variant phasing. *Genome Biology* 22(1):202
(2021). doi:[10.1186/s13059-021-02423-x](https://doi.org/10.1186/s13059-021-02423-x), PMID 34253237.

<a id="lange2022"></a>**[Lange et al. 2022] ecDNA segregation and evolution.**
Lange JT, Rose JC, Chen CY, et al. The evolutionary dynamics of extrachromosomal DNA in human
cancers. *Nature Genetics* 54(10):1527–1533 (2022).
doi:[10.1038/s41588-022-01177-x](https://doi.org/10.1038/s41588-022-01177-x), PMID 36123406.

<a id="yang2017"></a>**[Yang et al. 2017] NanoSim.**
Yang C, Chu J, Warren RL, Birol I. NanoSim: nanopore sequence read simulator based on statistical
characterization. *GigaScience* 6(4):1–6 (2017).
doi:[10.1093/gigascience/gix010](https://doi.org/10.1093/gigascience/gix010), PMID 28327957.

<a id="byrskabishop2022"></a>**[Byrska-Bishop et al. 2022] 1000 Genomes, high coverage.**
Byrska-Bishop M, Evani US, Zhao X, et al. High-coverage whole-genome sequencing of the expanded
1000 Genomes Project cohort including 602 trios. *Cell* 185(18):3426–3440.e19 (2022).
doi:[10.1016/j.cell.2022.08.004](https://doi.org/10.1016/j.cell.2022.08.004), PMID 36055201.

<a id="cnplot"></a>**[cnplot]** https://github.com/raphael-group/cnplot (default branch `main`).
Vendored at amplicon-allele-analysis `third_party/cnplot` (pinned commit 1e8bfef).

## Solvers and libraries

<a id="gurobi"></a>**[Gurobi]** Gurobi Optimizer Reference Manual, version 13.0.
https://docs.gurobi.com/projects/optimizer/en/current/. The pipeline used gurobipy 13.0.3 with an
academic WLS licence. The publisher's legal name was not shown on the page that was checked.

<a id="huangfu2018"></a>**[Huangfu & Hall 2018] HiGHS.**
Huangfu Q, Hall JAJ. Parallelizing the dual revised simplex method. *Mathematical Programming
Computation* 10(1):119–142 (2018).
doi:[10.1007/s12532-017-0130-5](https://doi.org/10.1007/s12532-017-0130-5).

<a id="virtanen2020"></a>**[Virtanen et al. 2020] SciPy.**
Virtanen P, Gommers R, Oliphant TE, et al. SciPy 1.0: fundamental algorithms for scientific
computing in Python. *Nature Methods* 17(3):261–272 (2020).
doi:[10.1038/s41592-019-0686-2](https://doi.org/10.1038/s41592-019-0686-2).

<a id="truong2020"></a>**[Truong et al. 2020] ruptures (change-point detection).**
Truong C, Oudre L, Vayatis N. Selective review of offline change point detection methods.
*Signal Processing* 167:107299 (2020).
doi:[10.1016/j.sigpro.2019.107299](https://doi.org/10.1016/j.sigpro.2019.107299).

## Algorithms and statistics

<a id="hierholzer1873"></a>**[Hierholzer & Wiener 1873] Euler circuits.**
Hierholzer C, Wiener C. Ueber die Möglichkeit, einen Linienzug ohne Wiederholung und ohne
Unterbrechung zu umfahren. *Mathematische Annalen* 6(1):30–32 (1873).
doi:[10.1007/BF01442866](https://doi.org/10.1007/BF01442866).

<a id="vatinlen2008"></a>**[Vatinlen et al. 2008] Minimum flow decomposition.**
Vatinlen B, Chauvet F, Chrétienne P, Mahey P. Simple bounds and greedy algorithms for decomposing a
flow into a minimal set of paths. *European Journal of Operational Research* 185(3):1390–1401
(2008). doi:[10.1016/j.ejor.2006.05.043](https://doi.org/10.1016/j.ejor.2006.05.043).

<a id="glover1975"></a>**[Glover 1975] Linearisation of products with binaries.**
Glover F. Improved linear integer programming formulations of nonlinear integer problems.
*Management Science* 22(4):455–460 (1975).
doi:[10.1287/mnsc.22.4.455](https://doi.org/10.1287/mnsc.22.4.455).

<a id="dantzig1954"></a>**[Dantzig et al. 1954] Subtour elimination cuts.**
Dantzig G, Fulkerson R, Johnson S. Solution of a large-scale traveling-salesman problem.
*Journal of the Operations Research Society of America* 2(4):393–410 (1954).
doi:[10.1287/opre.2.4.393](https://doi.org/10.1287/opre.2.4.393).

<a id="scott1974"></a>**[Scott & Knott 1974] Binary segmentation.**
Scott AJ, Knott M. A cluster analysis method for grouping means in the analysis of variance.
*Biometrics* 30(3):507 (1974). doi:[10.2307/2529204](https://doi.org/10.2307/2529204).

<a id="yao1988"></a>**[Yao 1988] BIC penalty for change points.**
Yao Y-C. Estimating the number of change-points via Schwarz' criterion. *Statistics & Probability
Letters* 6(3):181–189 (1988).
doi:[10.1016/0167-7152(88)90118-6](https://doi.org/10.1016/0167-7152(88)90118-6).

<a id="killick2012"></a>**[Killick et al. 2012] PELT.**
Killick R, Fearnhead P, Eckley IA. Optimal detection of changepoints with a linear computational
cost. *Journal of the American Statistical Association* 107(500):1590–1598 (2012).
doi:[10.1080/01621459.2012.737745](https://doi.org/10.1080/01621459.2012.737745).

<a id="rousseeuw1993"></a>**[Rousseeuw & Croux 1993] MAD as a robust scale.**
Rousseeuw PJ, Croux C. Alternatives to the median absolute deviation. *Journal of the American
Statistical Association* 88(424):1273–1283 (1993).
doi:[10.1080/01621459.1993.10476408](https://doi.org/10.1080/01621459.1993.10476408).

<a id="mahalanobis1936"></a>**[Mahalanobis 1936] Mahalanobis distance.**
Mahalanobis PC. On the generalised distance in statistics. *Proceedings of the National Institute
of Sciences of India* 2:49 (1936). Only the reprint was verified: *Sankhya A* 80(S1):1–7 (2018),
doi:[10.1007/s13171-019-00164-5](https://doi.org/10.1007/s13171-019-00164-5). The original's
volume and page come from secondary records.

<a id="brown2001"></a>**[Brown et al. 2001] Jeffreys interval for a binomial proportion.**
Brown LD, Cai TT, DasGupta A. Interval estimation for a binomial proportion. *Statistical Science*
16(2):101–133 (2001). doi:[10.1214/ss/1009213286](https://doi.org/10.1214/ss/1009213286).

<a id="jaccard1912"></a>**[Jaccard 1912] Jaccard index.**
Jaccard P. The distribution of the flora in the alpine zone. *New Phytologist* 11(2):37–50 (1912).
doi:[10.1111/j.1469-8137.1912.tb05611.x](https://doi.org/10.1111/j.1469-8137.1912.tb05611.x).

<a id="wagnerfischer1974"></a>**[Wagner & Fischer 1974] Dynamic programming over strings.**
Wagner RA, Fischer MJ. The string-to-string correction problem. *Journal of the ACM*
21(1):168–173 (1974). doi:[10.1145/321796.321811](https://doi.org/10.1145/321796.321811).

<a id="mackay2003"></a>**[MacKay 2003]** MacKay DJC. *Information Theory, Inference, and Learning
Algorithms*. Cambridge University Press (2003). Title, publisher and year were confirmed only
through published reviews (e.g. doi:[10.1145/1189056.1189063](https://doi.org/10.1145/1189056.1189063)).
No book DOI or ISBN was found.

## Not verified

These items appear in the docs as `[citation needed: …]`:

- **Kotzig A., "Moves without forbidden transitions in a graph", Matematický časopis (1968).**
  Proposed for Euler walks under transition constraints (the alternation in
  `rck/ecdna/walk.py::alternating_circuit`). Not found in Crossref.
- **A published method that orients phase blocks using reads across structural-variant junctions**
  (`simulations/subgraph_construct.md`, phase-block orientation), if one exists.
- **The source of the 1.15 over-dispersion factor** (`simulations/subgraph_construct.md`).
- **A published rule for choosing λ from a plateau of benchmark scores**
  (`results/bench_nested/README.md`, λ selection).
- **CoRAL supplementary methods.** The sections that define the cycle-interval overlap, the cyclic
  LCS normalised by the true length, and the log2 length error were not checked against the PDF
  (`Supplemental_Methods.pdf` is not on disk). The same applies to the section of the AA paper
  that defines error-free reconstruction and the section of the RCK paper that defines the
  novel-adjacency FNR and FDR. The papers are verified; the section pointers are not.

Also checked: the DOI the previous `simulations/README.md` gave for NanoSim
(10.1093/bioinformatics/btw753) is a different paper (Šošić & Šikić, Edlib, *Bioinformatics*
2017). The correct DOI is above.

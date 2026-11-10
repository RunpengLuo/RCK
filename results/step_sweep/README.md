# `--step-penalty` sweep (rck-ecdna decomposition)

λ = 0, 0.02, 0.05, 0.1 and 0.2 on every demo with a truth, run on Princeton ionic (Slurm job
31981121, 4 CPUs, 4.5 min) with `sweep.py` / `sweep.sbatch`. Paths in both point at the
original working folders.

- `sweep.tsv`: one row per case and λ. Columns: the rank-1 significant species, their weights
  and steps, the composite score per truth species, and the ambiguity report (`decided_by`,
  `gap`, `n_alternatives`).
- `runs/<case>/lambda_<λ>/`: the full rck-ecdna output of each run (`--no-plots`).
- `sweep.log`: the job log.
- `sweep_v1_milp_penalty.tsv` / `.log`: the first design, with the penalty inside the
  decomposition MILP (Slurm job 31981081; that code was replaced, so only the summary is kept).
  It split true duplications on demo3 (composite 0.96 → 0.72 at λ ≥ 0.02) and demo7 with
  `--block-flip` (ecdna2 1.00 → 0.77 at λ ≥ 0.1). This is why the penalty now acts in the
  ranking only (see `docs/ecDNA.md`, section 4).

Cases (subgraphs from amplicon-allele-analysis `simulations/runs/` unless noted):

| case | subgraph | extra options |
|---|---|---|
| demo3 | `tests/ecdna/data/amplicon_0` | `--depth-per-copy 15 --junction-reads-per-copy 6` |
| demo4_amp0, demo4_amp1 | `tests/ecdna/data/demo4/amplicon_{0,1}` | `--depth-per-copy 13.5 --truth-sample tumor1 --min-species-weight 2` |
| demo5 | `demo5_cycle1_ecdna/subgraph/amplicon_0` | |
| demo6_union | `demo6_cycle2_ecdna/subgraph_union/union` | |
| demo7_blockflip | `demo7_nestedcycles_ecdna/subgraph_union_blocks/union` | `--block-flip unoriented` |
| demo7_totalonly | `demo7_nestedcycles_ecdna/subgraph_union/union` | `--min-snps-allele 100000000` |
| demo8 | `demo8_nestedcycles_ecdna/subgraph_union/union` | |

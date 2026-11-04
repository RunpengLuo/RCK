# Environments of the graph route (Task 2)

One conda environment, in home: `~/.conda/envs/coral-3.1.1`.

| what | how it is run | pinned to |
|---|---|---|
| CoRAL 3.1.1 | `~/.conda/envs/coral-3.1.1/bin/coral cycle ...` (Gurobi: Slurm only) | CoRAL commit `7f3bf317`, runtime packages at its `poetry.lock` versions |
| CycleExtractor 1.0.0 | `~/.conda/envs/coral-3.1.1/bin/python -I ~/tools/CycleExtractor/CE.py --graph G --output O [--solver highs]` | code: git checkout of commit `1b29b33` in `~/tools/CycleExtractor` |

CycleExtractor needs numpy, networkx, gurobipy and (for `--solver highs`) highspy. The
coral-3.1.1 environment has the same versions the separate `cycle-extractor` environment had
(numpy 2.4.3, networkx 3.6.1, gurobipy 12.0.3, highspy 1.13.1; their `conda list --explicit`
files were identical), so that environment was removed on 2026-11-03 to free home quota
(5 GB; 287 MB -> 642 MB free) and CycleExtractor runs with the coral-3.1.1 python. Check, done
2026-11-03: CycleExtractor's README sample (`GBM39EC_EGFR_graph.txt`, HiGHS) gives the README's
cycle, `Segments=6+,2+,4+` (a rotation of `2+,4+,6+`), `Copy_count=85.18981770205714`
(README: 85.18981770204846, the graph's discordant CN is ...05714).

Lock files:

* `coral-3.1.1.conda-explicit.txt`: `conda list --explicit -p ~/.conda/envs/coral-3.1.1`
* `coral-3.1.1.pip-freeze.txt`: `~/.conda/envs/coral-3.1.1/bin/pip freeze --all`

Rebuild: `conda create -p ENV --file coral-3.1.1.conda-explicit.txt --override-channels -c conda-forge`,
then `ENV/bin/pip install -r coral-3.1.1.pip-freeze.txt` (CoRAL's line pins the commit; install it
with `--no-deps`), then `git clone https://github.com/AmpliconSuite/CycleExtractor ~/tools/CycleExtractor`
and `git -C ~/tools/CycleExtractor checkout 1b29b33`. (`--override-channels -c conda-forge`: the
host's conda refuses the default Anaconda channels until their terms of service are accepted.)

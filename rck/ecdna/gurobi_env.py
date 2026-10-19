"""Gurobi settings that depend on where the code runs."""

from __future__ import annotations

import logging
import os
from typing import Optional


def default_threads(requested: Optional[int] = None) -> int:
    """Threads Gurobi should use.

    An explicit request wins. Otherwise use the Slurm allocation (``SLURM_CPUS_PER_TASK``), then
    the CPUs this process may run on, so a 4-CPU job on an 80-core node does not start 32
    solver threads.
    """
    if requested:
        return int(requested)
    slurm = os.environ.get("SLURM_CPUS_PER_TASK")
    if slurm:
        try:
            value = int(slurm.split("(")[0])
            if value > 0:
                return value
        except ValueError:
            pass
    try:
        return max(len(os.sched_getaffinity(0)), 1)
    except AttributeError:  # macOS has no sched_getaffinity
        return max(os.cpu_count() or 1, 1)


def quiet_gurobi_logger(verbose: bool = False) -> None:
    """Keep gurobipy's Python logger off the console unless verbose.

    gurobipy mirrors its whole solver log to the ``gurobipy`` logger, so a root logger at INFO
    prints it even with ``LogToConsole=0``. The full log still goes to each model's LogFile.
    """
    logging.getLogger("gurobipy").setLevel(logging.DEBUG if verbose else logging.WARNING)

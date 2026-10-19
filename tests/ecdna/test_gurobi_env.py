import logging
import os
import unittest
from unittest import mock

from rck.ecdna.gurobi_env import close_model, default_threads, quiet_gurobi_logger


class TestGurobiEnv(unittest.TestCase):
    def test_explicit_request_wins(self):
        with mock.patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "4"}):
            self.assertEqual(default_threads(2), 2)

    def test_slurm_allocation(self):
        with mock.patch.dict(os.environ, {"SLURM_CPUS_PER_TASK": "4"}):
            self.assertEqual(default_threads(None), 4)

    def test_without_slurm(self):
        env = {k: v for k, v in os.environ.items() if k != "SLURM_CPUS_PER_TASK"}
        with mock.patch.dict(os.environ, env, clear=True):
            n = default_threads(None)
        self.assertGreaterEqual(n, 1)
        self.assertLessEqual(n, os.cpu_count() or n)

    def test_close_model_tolerates_none(self):
        close_model(None)

    def test_quiet_logger(self):
        quiet_gurobi_logger(verbose=False)
        self.assertEqual(logging.getLogger("gurobipy").level, logging.WARNING)
        quiet_gurobi_logger(verbose=True)
        self.assertEqual(logging.getLogger("gurobipy").level, logging.DEBUG)
        quiet_gurobi_logger(verbose=False)


if __name__ == "__main__":
    unittest.main()

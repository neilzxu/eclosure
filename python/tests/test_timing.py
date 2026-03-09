from __future__ import annotations

import os
import time

import numpy as np
import pytest

from eclosure import closedBY, closedeBH


def _skip_unless_enabled() -> None:
  if os.getenv("ECLOSURE_RUN_TIMING_TESTS") != "1":
    pytest.skip("Set ECLOSURE_RUN_TIMING_TESTS=1 to run timing benchmarks")


def test_timing_cebh_on_1e6_evalues():
  _skip_unless_enabled()

  rng = np.random.default_rng(1)
  evalues = rng.exponential(scale=5.0, size=1_000_000)

  start = time.perf_counter()
  res = closedeBH(evalues, alpha=0.05)
  elapsed = time.perf_counter() - start

  assert isinstance(res, (int, np.integer))
  assert res >= 0
  print(f"\nclosedeBH timing (1e6 e-values): {elapsed:.3f} seconds")


def test_timing_cby_on_1e5_pvalues():
  _skip_unless_enabled()

  rng = np.random.default_rng(1)
  pvalues = rng.random(100_000) ** 4

  start = time.perf_counter()
  res = closedBY(pvalues, alpha=0.05)
  elapsed = time.perf_counter() - start

  assert isinstance(res, (int, np.integer))
  assert res >= 0
  print(f"\nclosedBY timing (1e5 p-values): {elapsed:.3f} seconds")

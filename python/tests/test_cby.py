from __future__ import annotations

import numpy as np

from eclosure import closedBY


def test_cby_empty_subset_is_always_significant():
  pvalues = [0.01, 0.2, 0.8]
  assert closedBY(pvalues, subset=[False, False, False], alpha=0.05) is True
  assert closedBY(pvalues, subset=[], alpha=0.05) is True


def test_cby_empty_input_returns_zero():
  assert closedBY([], alpha=0.05) == 0


def test_cby_alpha_zero_matches_r_behavior():
  pvalues = [0.001, 0.01, 0.5]
  assert closedBY(pvalues, alpha=0.0) == 0
  assert closedBY(pvalues, subset=[0, 1], alpha=0.0) is False


def test_cby_exact_zeros_are_handled():
  pvalues = [0.0, 0.01, 0.02, 0.5]
  assert closedBY(pvalues, alpha=0.05) >= 0
  assert isinstance(closedBY(pvalues, subset=[0, 1], alpha=0.05), bool)


def test_cby_discovery_mode_top_r_set_is_significant():
  pvalues = np.array([0.8] * 15 + [1e-8] * 5)
  discoveries = closedBY(pvalues, alpha=0.05)
  assert discoveries > 0

  top_r = np.argsort(pvalues)[:discoveries]
  assert closedBY(pvalues, subset=top_r, alpha=0.05) is True


def test_cby_approximate_never_exceeds_exact():
  rng = np.random.default_rng(11)
  for _ in range(20):
    pvalues = rng.uniform(size=40)
    exact = closedBY(pvalues, alpha=0.05, approximate=False)
    approximate = closedBY(pvalues, alpha=0.05, approximate=True)
    assert approximate <= exact


def test_cby_discovery_result_is_invariant_to_input_order():
  rng = np.random.default_rng(23)
  pvalues = rng.uniform(size=30)
  shuffled = rng.permutation(pvalues)
  assert closedBY(pvalues, alpha=0.05) == closedBY(shuffled, alpha=0.05)

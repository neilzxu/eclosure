from __future__ import annotations

import numpy as np

from eclosure import closedeBH


def test_cebh_empty_subset_is_always_mean_consistent():
  evalues = [0.5, 1.0, 2.0]
  assert closedeBH(evalues, subset=[False, False, False], alpha=0.05) is True
  assert closedeBH(evalues, subset=[], alpha=0.05) is True


def test_cebh_empty_input_returns_zero():
  assert closedeBH([], alpha=0.05) == 0


def test_cebh_single_hypothesis_boundary_behavior():
  alpha = 0.05
  assert closedeBH([1 / alpha], alpha=alpha) == 1
  assert closedeBH([1 / alpha - 1e-6], alpha=alpha) == 0
  assert closedeBH([100.0], alpha=alpha) == 1
  assert closedeBH([0.0], alpha=alpha) == 0


def test_cebh_all_ones_gives_no_rejections():
  assert closedeBH([1.0] * 20, alpha=0.05) == 0


def test_cebh_all_boundary_values_rejects_everything():
  alpha = 0.05
  assert closedeBH([1 / alpha] * 10, alpha=alpha) == 10


def test_cebh_discovery_mode_top_r_set_is_consistent():
  evalues = np.array([1.0] * 15 + [500.0] * 5)
  discoveries = closedeBH(evalues, alpha=0.05)
  assert discoveries > 0

  threshold = np.sort(evalues)[::-1][discoveries - 1]
  top_r = evalues >= threshold
  assert closedeBH(evalues, subset=top_r, alpha=0.05) is True


def test_cebh_approximate_never_exceeds_exact():
  rng = np.random.default_rng(7)
  for _ in range(30):
    rate = rng.uniform(0.1, 2.0)
    evalues = rng.exponential(scale=1.0 / rate, size=50)
    exact = closedeBH(evalues, alpha=0.05, approximate=False)
    approximate = closedeBH(evalues, alpha=0.05, approximate=True)
    assert approximate <= exact


def test_cebh_discovery_result_is_invariant_to_input_order():
  rng = np.random.default_rng(9)
  evalues = rng.exponential(scale=1.0 / 0.4, size=30)
  shuffled = rng.permutation(evalues)
  assert closedeBH(evalues, alpha=0.05) == closedeBH(shuffled, alpha=0.05)

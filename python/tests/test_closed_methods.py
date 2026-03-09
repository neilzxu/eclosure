from __future__ import annotations

import numpy as np
import pytest

pd = pytest.importorskip("pandas")
pytest.importorskip("scipy")
pytest.importorskip("tqdm")

from eclosure import closedeBH
from eclosure.experiments import (compute_by, compute_cby,
                                  compute_cebh_discovery_set, compute_ebh,
                                  simulate_statistics)
from eclosure.resources import data_dir


def _load_evalue_frames() -> list[pd.DataFrame]:
  base = data_dir("evalues")
  return [pd.read_csv(path) for path in sorted(base.glob("*_evalues.csv"))]


def _load_pvalue_arrays() -> list[np.ndarray]:
  base = data_dir("pvalues")
  arrays = []
  for path in sorted(base.glob("*.csv")):
    df = pd.read_csv(path)
    arrays.append(df["pvalue"].to_numpy())
  return arrays


def _cebh_discoveries(evalues: np.ndarray | list[float],
                      alpha: float) -> set[int]:
  arr = np.asarray(evalues, dtype=float)
  k = int(closedeBH(arr, alpha=alpha))
  if k <= 0:
    return set()
  sorted_idx = np.argsort(arr)[::-1]
  return set(map(int, sorted_idx[:k]))


def test_closed_ebh_python_and_cpp_agree_on_simulations():
  np.random.seed(2024)
  configs = [(25, 0.5, 1.5), (40, 0.6, 2.0)]
  for K, pi0, mu in configs:
    evalues, _ = simulate_statistics(K=K,
                                     null_prop=pi0,
                                     signal_strength=mu,
                                     mode='simple')
    cpp = _cebh_discoveries(evalues, alpha=0.1)
    py = compute_cebh_discovery_set(evalues, alpha=0.1, mode='simple')
    ebh = compute_ebh(evalues, alpha=0.1)

    assert cpp == py
    assert ebh.issubset(cpp)


def test_closed_ebh_superset_on_simulations():
  np.random.seed(123)
  evalues, _ = simulate_statistics(K=40,
                                   null_prop=0.6,
                                   signal_strength=1.5,
                                   mode='simple')
  closed = _cebh_discoveries(evalues, alpha=0.1)
  ebh = compute_ebh(evalues, alpha=0.1)
  assert ebh.issubset(closed)


def test_closed_ebh_superset_on_real_data():
  datasets = _load_evalue_frames()
  assert datasets, "No packaged e-value datasets found"
  for alpha in (0.05, 0.1):
    for df in datasets:
      evalues = df["evalue"].dropna().to_numpy()
      closed = _cebh_discoveries(evalues, alpha=alpha)
      ebh = compute_ebh(evalues, alpha=alpha)
      assert ebh.issubset(
          closed), f"Closed eBH failed superset check at alpha={alpha}"


def test_closed_ebh_python_and_cpp_agree_on_real_data():
  datasets = _load_evalue_frames()
  assert datasets, "No packaged e-value datasets found"
  for alpha in (0.05, 0.1):
    for df in datasets:
      evalues = df["evalue"].dropna().to_numpy()
      cpp = _cebh_discoveries(evalues, alpha=alpha)
      py = compute_cebh_discovery_set(evalues, alpha=alpha, mode='simple')
      ebh = compute_ebh(evalues, alpha=alpha)

      assert cpp == py
      assert ebh.issubset(cpp)


def test_closed_by_superset_on_simulations():
  np.random.seed(321)
  pvalues, _ = simulate_statistics(K=25,
                                   null_prop=0.5,
                                   signal_strength=3.0,
                                   mode='simple-p',
                                   alpha=0.1,
                                   rho=0.05)
  pvalues = np.asarray(pvalues)
  closed = compute_cby(pvalues, alpha=0.1)
  by = compute_by(pvalues, alpha=0.1)
  assert by.issubset(closed)


def test_closed_by_superset_on_real_data():
  datasets = _load_pvalue_arrays()
  assert datasets, "No packaged p-value datasets found"
  for alpha in (0.05, 0.1):
    for pvalues in datasets:
      closed = compute_cby(pvalues, alpha=alpha)
      by = compute_by(pvalues, alpha=alpha)
      assert by.issubset(
          closed), f"Closed BY failed superset check at alpha={alpha}"

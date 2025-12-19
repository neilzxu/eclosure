from __future__ import annotations

from dataclasses import dataclass
from typing import Iterable, Optional, Sequence, Set

import numpy as np

try:
  from . import _cebh
except Exception as exc:  # pragma: no cover - import errors should be explicit
  raise ImportError(
      "The eclosure C++ bindings could not be imported. "
      "Run `pip install .` from the repository root to build the extension."
  ) from exc


class ClosedeBHError(RuntimeError):
  """Raised when the closed e-BH backend fails."""


def _ensure_array(evalues: Sequence[float]) -> np.ndarray:
  arr = np.asarray(evalues, dtype=float)
  if arr.ndim != 1:
    raise ValueError("evalues must be a 1-D sequence")
  return arr


def largest_mean_consistent_size(
    evalues: Sequence[float],
    alpha: float = 0.05,
    approximate: bool = False,
) -> int:
  arr = _ensure_array(evalues)
  if not (0 < alpha <= 1):
    raise ValueError("alpha must lie in (0, 1]")
  try:
    return int(
        _cebh.largest_mean_consistent_k(arr.tolist(), float(alpha), bool(approximate)))
  except Exception as exc:  # pragma: no cover - surfaces binding failures
    raise ClosedeBHError("C++ backend failed to compute a consistent size") from exc


def is_mean_consistent_subset(
    evalues: Sequence[float],
    subset: Iterable[int],
    alpha: float = 0.05,
) -> bool:
  arr = _ensure_array(evalues)
  subset_arr = np.asarray(list(subset), dtype=int)
  if subset_arr.size == 0:
    return True
  if subset_arr.min(initial=0) < 0 or subset_arr.max(initial=0) >= arr.size:
    raise ValueError("subset indices must lie in [0, len(evalues))")
  try:
    return bool(
        _cebh.is_subset_mean_consistent(arr.tolist(), subset_arr.tolist(), float(alpha)))
  except Exception as exc:  # pragma: no cover
    raise ClosedeBHError(
        "C++ backend failed to check subset mean consistency") from exc


def closede_bh_discoveries(
    evalues: Sequence[float],
    alpha: float = 0.05,
    approximate: bool = False,
) -> Set[int]:
  arr = _ensure_array(evalues)
  k = largest_mean_consistent_size(arr, alpha=alpha, approximate=approximate)
  if k <= 0:
    return set()
  sorted_idx = np.argsort(arr)[::-1]
  return set(map(int, sorted_idx[:k]))


@dataclass
class ClosedeBHRunner:
  """State-free runner mirroring the old subprocess-based API."""

  alpha: float = 0.05
  approximate: bool = False

  def discovery_set(self, evalues: Sequence[float]) -> Set[int]:
    return closede_bh_discoveries(
        evalues, alpha=self.alpha, approximate=self.approximate)

  def largest_mean_consistent_size(self, evalues: Sequence[float]) -> int:
    return largest_mean_consistent_size(
        evalues, alpha=self.alpha, approximate=self.approximate)

  def check_subset(self, evalues: Sequence[float], subset: Iterable[int]) -> bool:
    return is_mean_consistent_subset(evalues, subset, alpha=self.alpha)


_DEFAULT_RUNNER: Optional[ClosedeBHRunner] = None


def get_default_runner() -> ClosedeBHRunner:
  global _DEFAULT_RUNNER
  if _DEFAULT_RUNNER is None:
    _DEFAULT_RUNNER = ClosedeBHRunner()
  return _DEFAULT_RUNNER

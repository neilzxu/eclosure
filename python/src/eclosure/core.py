"""Core Python wrappers around the compiled closed e-BH and closed BY routines."""

from __future__ import annotations

from typing import Iterable, Sequence

import numpy as np

try:
  from . import _cebh
except Exception as exc:  # pragma: no cover - import errors should be explicit
  raise ImportError(
      "The eclosure C++ bindings could not be imported. "
      "Run `pip install .` from the repository root to build the extension."
  ) from exc


def _ensure_array(evalues: Sequence[float]) -> np.ndarray:
  arr = np.asarray(evalues, dtype=float)
  if arr.ndim != 1:
    raise ValueError("evalues must be a 1-D sequence")
  return arr


def _validate_alpha(alpha: float) -> None:
  if not (0 < alpha <= 1):
    raise ValueError("alpha must lie in (0, 1]")


def _normalize_subset(subset: Iterable[int], size: int) -> np.ndarray:
  raw = np.asarray(list(subset))
  if raw.ndim != 1:
    raise ValueError("subset must be a 1-D sequence")
  if np.issubdtype(raw.dtype, np.bool_):
    if raw.size != size:
      raise ValueError("boolean subset must have length len(evalues)")
    return np.flatnonzero(raw)

  subset_arr = np.asarray(raw, dtype=int)
  if subset_arr.size == 0:
    return subset_arr
  if subset_arr.min(initial=0) < 0 or subset_arr.max(initial=0) >= size:
    raise ValueError("subset indices must lie in [0, len(evalues))")
  return subset_arr


def _harmonic_numbers(size: int) -> np.ndarray:
  if size <= 0:
    return np.array([0.0], dtype=float)
  return np.concatenate(([0.0], np.cumsum(1.0 / np.arange(1, size + 1))))


def closedeBH(
    evalues: Sequence[float],
    subset: Iterable[int] | None = None,
    alpha: float = 0.05,
    approximate: bool = False,
) -> int | bool:
  """Run closed e-BH in discovery mode or set-checking mode.

  Args:
    evalues: One-dimensional sequence of e-values.
    subset: If `None`, return the number of discoveries. Otherwise check
      whether the specified subset is mean consistent. A subset may be given
      either as a boolean mask or as 0-based integer indices.
    alpha: Target FDR level in `(0, 1]`.
    approximate: If `True`, use the approximate discovery-mode routine.

  Returns:
    An integer discovery count when `subset is None`, otherwise a single
    boolean.

  Examples:
    >>> closedeBH([1.0, 100.0], alpha=0.05)
    1
    >>> closedeBH([1.0, 100.0], subset=[1], alpha=0.05)
    True
  """
  arr = _ensure_array(evalues)
  _validate_alpha(alpha)

  if subset is None:
    sorted_arr = np.sort(arr)
    try:
      if approximate:
        return int(
            _cebh.largestmeanconsistent_approximate_cpp(
                sorted_arr.tolist(), float(alpha)))
      return int(_cebh.largestmeanconsistent_cpp(sorted_arr.tolist(), float(alpha)))
    except Exception as exc:  # pragma: no cover - surfaces binding failures
      raise RuntimeError(
          "C++ backend failed to compute a consistent size") from exc

  subset_arr = _normalize_subset(subset, arr.size)
  if subset_arr.size == 0:
    return True

  membership = np.ones(arr.size, dtype=bool)
  membership[subset_arr] = False
  e_out = np.sort(arr[membership])
  e_in = np.sort(arr[subset_arr])
  cum_e = np.cumsum(np.concatenate(([0.0], e_out, e_in)))

  try:
    return bool(
        _cebh.meanconsistent_wrapper_cpp(
            cum_e.tolist(), int(subset_arr.size), float(alpha)))
  except Exception as exc:  # pragma: no cover
    raise RuntimeError(
        "C++ backend failed to check subset mean consistency") from exc


def closedBY(
    pvalues: Sequence[float],
    subset: Iterable[int] | None = None,
    alpha: float = 0.05,
    approximate: bool = False,
) -> int | bool:
  """Run closed BY in discovery mode or set-checking mode.

  Args:
    pvalues: One-dimensional sequence of p-values.
    subset: If `None`, return the number of discoveries. Otherwise check
      whether the specified subset is significant. A subset may be given
      either as a boolean mask or as 0-based integer indices.
    alpha: Target FDR level in `[0, 1]`.
    approximate: If `True`, use the approximate discovery-mode routine.

  Returns:
    An integer discovery count when `subset is None`, otherwise a single
    boolean.

  Examples:
    >>> closedBY([0.01, 0.5], alpha=0.05)
    1
    >>> closedBY([0.01, 0.5], subset=[0], alpha=0.05)
    True
  """
  arr = _ensure_array(pvalues).copy()
  if not (0 <= alpha <= 1):
    raise ValueError("alpha must lie in [0, 1]")

  m = arr.size
  arr[arr == 0] = 1e-12

  if subset is not None:
    subset_arr = _normalize_subset(subset, m)
    if subset_arr.size == 0:
      return True
    if alpha == 0:
      return False

    harmonic = _harmonic_numbers(m).tolist()
    outside_mask = np.ones(m, dtype=bool)
    outside_mask[subset_arr] = False
    ordered = np.concatenate(
        (np.sort(arr[outside_mask])[::-1], np.sort(arr[subset_arr])[::-1]))
    try:
      result = _cebh.cBY_check_cpp(
          ordered.tolist(), int(subset_arr.size), harmonic, 0, float(alpha))
      return bool(result["res"])
    except Exception as exc:  # pragma: no cover
      raise RuntimeError("C++ backend failed to check closedBY significance") from exc

  if alpha == 0:
    return 0

  sorted_arr = np.sort(arr)[::-1]
  try:
    if approximate:
      return int(
          _cebh.largestcBYsignificant_approximate_cpp(
              sorted_arr.tolist(), float(alpha)))
    return int(_cebh.largestcBYsignificant_cpp(sorted_arr.tolist(), float(alpha)))
  except Exception as exc:  # pragma: no cover
    raise RuntimeError("C++ backend failed to compute closedBY size") from exc

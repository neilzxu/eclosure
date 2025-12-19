from __future__ import annotations

from pathlib import Path
from typing import Literal

_THIS_DIR = Path(__file__).resolve().parent


def _repo_root() -> Path:
  # When installed as a package, fall back to the site-packages parent.
  parts = list(_THIS_DIR.parents)
  if len(parts) >= 3:
    return parts[2]
  return parts[-1]


def data_dir(kind: Literal["evalues", "pvalues"]) -> Path:
  local = _THIS_DIR / "data" / kind
  if local.exists():
    return local
  repo_candidate = _repo_root() / "r" / "eclosure" / "inst" / "extdata" / kind
  if repo_candidate.exists():
    return repo_candidate
  raise FileNotFoundError(f"Unable to locate packaged {kind} data.")


def repo_root() -> Path:
  return _repo_root()

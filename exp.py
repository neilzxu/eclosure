"""
Backward-compatible shim that delegates to the packaged experiments CLI.

This keeps `python exp.py ...` working while relying on the canonical
implementation in `eclosure.experiments`.
"""

from eclosure.experiments import *  # noqa: F401,F403
from eclosure.experiments import main as _main


if __name__ == "__main__":
    raise SystemExit(_main())

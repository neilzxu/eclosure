# Repository Scripts

The repository ships a few research-oriented Python scripts outside the
published package. They are useful for reproducibility and exploratory work,
but they are not included in wheels or treated as stable package interfaces.

## `eclosure.experiments`

This module powers the packaged CLI documented in [Scripts](scripts.md). It
contains simulation grids, real-data utilities, and optional plotting helpers.
Use it for reproducing paper experiments and generating result tables.

Source: <https://github.com/neilzxu/eclosure/blob/master/python/src/eclosure/experiments.py>

## `scripts/jelle.py`

This script contains small-scale exploratory utilities for comparing BH-style
procedures and e-value constructions on toy examples. Run it from the repo
root with:

```bash
python scripts/jelle.py
```

Source: <https://github.com/neilzxu/eclosure/blob/master/scripts/jelle.py>

## `scripts/sweep.py`

This script contains exploratory sweep code for simulation studies and plotting
under optional dependencies. It remains repository support code and is not
packaged for installation. Run it from the repo root with:

```bash
python scripts/sweep.py --help
```

Source: <https://github.com/neilzxu/eclosure/blob/master/scripts/sweep.py>

# eclosure

`eclosure` provides Python bindings and research utilities for e-Closure based
multiple testing procedures. The package wraps a compiled C++ backend for the
core closed e-BH and closed BY routines and ships example datasets used in the
paper.

- Documentation: <https://neilzxu.github.io/eclosure/>
- Source: <https://github.com/neilzxu/eclosure>
- Paper: <https://arxiv.org/abs/2509.02517>

## Installation

Install the core package from PyPI:

```bash
python -m pip install eclosure
```

Optional extras:

```bash
# experiment and real-data helpers
python -m pip install "eclosure[experiments]"

# plotting helpers used by the experiment code
python -m pip install "eclosure[plotting]"
```

Supported releases ship binary wheels for common CPython platforms. If you
build from source, you will need a working C++17 toolchain.

## Quick start

```python
import eclosure

eclosure.closedeBH([1.0, 100.0], alpha=0.05)
eclosure.closedBY([0.01, 0.5], alpha=0.05)
```

The top-level API focuses on the stable testing procedures. The experiments
workflow is available through the packaged module:

```bash
python -m eclosure.experiments --mode real --alpha 0.1 --out-dir figures
```

## Learn more

The full package documentation includes:

- API reference for `closedBY`, `closedeBH`, and packaged data helpers
- Script documentation for `python -m eclosure.experiments`
- Testing guidance for unit, real-data, and timing suites
- Notes on repository-only exploratory scripts that are not part of the package

See <https://neilzxu.github.io/eclosure/> for the full docs.

## Citation

```bibtex
@article{xu2025bringing,
  title={Bringing Closure to False Discovery Rate Control: A General Principle for Multiple Testing},
  author={Xu, Ziyu and Solari, Aldo and Fischer, Lasse and de Heide, Rianne and Ramdas, Aaditya and Goeman, Jelle},
  journal={arXiv preprint arXiv:2509.02517},
  year={2025}
}
```

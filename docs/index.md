# eclosure

`eclosure` packages the Python interface for the e-Closure procedures together
with research helpers and packaged example datasets. The stable public surface
is intentionally small: the top-level API exposes `closedBY` and `closedeBH`,
while the repository also bundles experiment-oriented modules used in the
paper.

## Installation

Install the core package:

```bash
python -m pip install eclosure
```

Install optional extras for experiments, plotting, or docs work:

```bash
python -m pip install "eclosure[experiments]"
python -m pip install "eclosure[plotting]"
python -m pip install -e ".[docs]"
```

Source builds require a working C++17 toolchain. Tagged releases are intended
to ship wheels for Linux, macOS, and Windows on supported CPython versions.

## Quick start

```python
import eclosure

discoveries = eclosure.closedeBH([1.0, 100.0], alpha=0.05)
significant = eclosure.closedBY([0.01, 0.5], subset=[0], alpha=0.05)
```

## Documentation map

- [API](api.md): stable Python API reference for the closed testing procedures.
- [Scripts](scripts.md): CLI usage for `python -m eclosure.experiments`.
- [Tests](tests.md): how the Python test suites are organized and run.

## Related links

- GitHub repository: <https://github.com/neilzxu/eclosure>
- Paper: <https://arxiv.org/abs/2509.02517>

## Citation

```bibtex
@article{xu2025bringing,
  title={Bringing Closure to False Discovery Rate Control: A General Principle for Multiple Testing},
  author={Xu, Ziyu and Solari, Aldo and Fischer, Lasse and de Heide, Rianne and Ramdas, Aaditya and Goeman, Jelle},
  journal={arXiv preprint arXiv:2509.02517},
  year={2025}
}
```

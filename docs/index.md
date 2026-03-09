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

- [API](api.md): stable Python API reference and packaged data helpers.
- [Scripts](scripts.md): CLI usage for `python -m eclosure.experiments`.
- [Tests](tests.md): how the Python test suites are organized and run.
- [Advanced Modules](advanced.md): research-oriented modules included with the
  repository but not treated as part of the stable surface.

## Related links

- GitHub repository: <https://github.com/neilzxu/eclosure>
- Release artifacts: <https://github.com/neilzxu/eclosure/releases>
- Paper: <https://arxiv.org/abs/2509.02517>

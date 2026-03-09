# Procedures based on the e-Closure Principle

This repository hosts R and Python interfaces for the e-Closure procedures in
this [paper](https://arxiv.org/abs/2509.02517), together with simulations and
reproducible real-data examples. Both frontends call the same underlying C++
implementations in `r/eclosure/src/`.

## Layout

- `r/eclosure/` – the R package, published as `eClosure`
- `python/` – the Python package source tree, with the package under
  `python/src/eclosure/`
- `r/eclosure/src/` – the C++ implementations used by both the R package and
  the Python extension module
- `python/tests/` – pytest-based Python tests
- `r/eclosure/tests/testthat/` – `testthat`-based R tests
- `data/` – shared datasets used in tests and examples
- `raw_data/` – one-off scripts for regenerating packaged datasets

## Python package

Package documentation for the Python interface lives at
<https://neilzxu.github.io/eclosure/>. Tagged release artifacts are published
at <https://github.com/neilzxu/eclosure/releases>.
Repository-only exploratory scripts live under `scripts/` and are not shipped
in the Python package or wheels.

Install only the core Python package from the repository root:

```bash
python -m pip install -e .
```

Optional extras:

```bash
# experiment/data helpers
python -m pip install -e .[experiments]

# plotting helpers used by the figure-generation code
python -m pip install -e .[plotting]

# Python tests that rely on the experiment helpers
python -m pip install -e .[test]

# full development environment
python -m pip install -e .[dev]
```

Build the docs locally with:

```bash
python -m pip install -e .[docs]
mkdocs build
```

Run the Python unit tests with:

```bash
python -m pip install -e .[test]
python -m pytest python/tests
```

Run the Python timing benchmarks with:

```bash
ECLOSURE_RUN_TIMING_TESTS=1 python -m pytest -s python/tests/test_timing.py
```

Build a wheel and source distribution with:

```bash
python -m pip install build twine
python -m build .
python -m twine check dist/*
```

The GitHub Actions release workflow in the public `neilzxu/eclosure` repo
builds the same artifacts on tag pushes matching `v*` or on manual
`workflow_dispatch`. It produces:

- one source distribution
- CPython 3.9-3.13 wheels for Linux `x86_64`
- CPython 3.9-3.13 wheels for Windows `AMD64`
- CPython 3.9-3.13 wheels for macOS `universal2`

On `v*` tag pushes, the workflow also creates or updates the corresponding
GitHub Release and attaches the built `sdist` and wheel files as release
assets.

Each wheel job installs the built wheel and runs the non-timing smoke tests:

```bash
python -m pytest python/tests/test_cebh.py python/tests/test_cby.py python/tests/test_closed_methods.py
```

The docs workflow in the same public repo builds the MkDocs site on pull
requests and on pushes to the default branch, then deploys GitHub Pages from
the default branch build.

After downloading verified artifacts from the release workflow or from a GitHub
Release page, the final PyPI upload remains a manual step:

```bash
python -m twine upload dist/*
```

### Release from GitHub Actions

To build release artifacts in GitHub Actions instead of locally:

1. Push your branch to the public `neilzxu/eclosure` repository.
2. Open `Actions` in GitHub.
3. Select the `release-artifacts` workflow.
4. Click `Run workflow` and choose the branch to build.

You can also trigger the same workflow by pushing a version tag:

```bash
git tag v0.1.1
git push origin v0.1.1
```

If you triggered the workflow with a pushed `v*` tag, open the corresponding
GitHub Release page to download the attached assets. You can also download the
same files from the workflow run page. For a manual `Run workflow` build, use
the workflow run page artifacts:

- `sdist`
- `wheels-ubuntu-latest`
- `wheels-windows-latest`
- `wheels-macos-latest`

Unpack all downloaded artifacts into a local `dist/` directory, then verify and
upload them:

```bash
python -m pip install --upgrade twine
python -m twine check dist/*
python -m twine upload --repository-url https://test.pypi.org/legacy/ dist/*
```

For production PyPI, use:

```bash
python -m twine upload dist/*
```

The `release-artifacts` workflow is not triggered by pull requests, so use
`Run workflow` or a pushed `v*` tag when you want downloadable wheel artifacts.

Example:

```python
import eclosure

eclosure.closedeBH([1.0, 100.0], alpha=0.05)
eclosure.closedBY([0.01, 0.5], alpha=0.05)
```

Repository-only exploratory scripts can be run from the checkout root:

```bash
python scripts/jelle.py
python scripts/sweep.py --help
```

## R package

Install the R package from the repository root with:

```bash
R CMD INSTALL r/eclosure
```

Run the R unit tests with:

```bash
Rscript -e 'devtools::test("r/eclosure")'
```

Run the R timing benchmarks with:

```bash
ECLOSURE_RUN_TIMING_TESTS=1 Rscript -e 'devtools::test("r/eclosure", filter = "timing")'
```

Build the R source package with:

```bash
R CMD build r/eclosure
```

For interactive development, load the package in-place with:

```r
devtools::load_all("r/eclosure")
```

Example:

```r
library(eClosure)

closedeBH(c(1, 100), alpha = 0.05)
closedBY(c(0.01, 0.5), alpha = 0.05)
```

## Testing

The repository has separate unit tests and timing benchmarks for both
frontends. The timing benchmarks are guarded by the
`ECLOSURE_RUN_TIMING_TESTS=1` environment variable so they stay out of the
default test run.

- Python unit tests: `python -m pip install -e .[test] && python -m pytest python/tests`
- Python timing benchmarks: `python -m pip install -e .[test] && ECLOSURE_RUN_TIMING_TESTS=1 python -m pytest -s python/tests/test_timing.py`
- R unit tests: `Rscript -e 'devtools::test("r/eclosure")'`
- R timing benchmarks: `ECLOSURE_RUN_TIMING_TESTS=1 Rscript -e 'devtools::test("r/eclosure", filter = "timing")'`

The packaged Python test guidance is mirrored in the published docs under
<https://neilzxu.github.io/eclosure/tests/>.

## Real-data result tables

To regenerate the tables for all real datasets used in the paper:

- E-value experiments (i.e., crypto)

```bash
python -m pip install -e .[experiments]
python -m eclosure.experiments --mode real --alpha 0.1 --out-dir figures
```

This scans `data/evalues` by default (override with `--data-dir` or
`--col-name`), prints the table to stdout, and writes a CSV to
`figures/real/real_results_alpha_0.1.csv`. Repeat with another `--alpha` to
produce additional tables.

- p-value datasets for closed BY (APSAC, NAEP, PADJUST, PVALUES, VANDEVIJVER,
GOLUB): install the R package first, then run:

```bash
R CMD INSTALL r/eclosure
Rscript r/eclosure/inst/scripts/examples_BYbar.R > figures/by_real_data_table.tex
```

The script prints the summary data frame and the LaTeX table; drop the output
redirection if you just want to inspect it in the console.

# Citation

```latex
@article{xu2025bringing,
  title={Bringing Closure to False Discovery Rate Control: A General Principle for Multiple Testing},
  author={Xu, Ziyu and Solari, Aldo and Fischer, Lasse and de Heide, Rianne and Ramdas, Aaditya and Goeman, Jelle},
  journal={arXiv preprint arXiv:2509.02517},
  year={2025}
}
```

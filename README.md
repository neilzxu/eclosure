# Procedures based on the e-Closure Principle

This repository hosts a shared C++ implementation of multiple testing procedures for this [paper](https://arxiv.org/abs/2509.02517) together with
packaged R and Python interfaces, simulations, and reproducible real-data
examples.

## Layout

- `r/eclosure/` – R package exposing closed e-BH and closed BY utilities,
including `inst/extdata` with e-value and p-value datasets plus `testthat`
regression tests. These compare the cpp implementation of closed eBH with existing implementations of closed eBH and eBH in R.
Note that the cpp source file is also contained directly in the R package at the moment.
- `r/eclosure/inst/scripts/examples_BYbar.R` – additional dataset suite that generates results of closed BY procedures on real datasets.
- `python/src/eclosure/` – Python package with the pybind11 extension,
simulation helpers (former `exp.py`), and data copies under `data/`.
- `python/tests/` – pytest suite covering simulation and real-data checks. These compare the cpp implementation of closed eBH with existing implementations of closed eBH and eBH in python.
- `raw_data/` – one-off scripts for regenerating packaged datasets:
- `raw_data/crypto.py`: generates crypto e-values into `data/evalues/`
  from the CSVs in `raw_data/`.
- `raw_data/run_real_data_analysis.R`: downloads/cleans NYC taxi and NASDAQ
  data and writes e-values to `data/evalues/` (runs from repo root).
- `data/` - data generated from raw data that is used by python and R packages for unit tests.

## Python package

Install in editable mode from the repository root:

```bash
pip install -e .
```

Run the shared tests with:

```bash
PYTHONPATH=python/src pytest python/tests
```

## R package

From `r/eclosure`, use:

```r
devtools::load_all()
devtools::test()
```

The `tests/testthat` directory exercises the same real and simulated datasets as
the Python tests, confirming that the Rcpp bindings respect the original R
implementations.

## Real-data result tables

To regenerate the tables for all real datasets used in the paper:

- E-value experiments (i.e., crypto)

```bash
python -m eclosure.experiments --mode real --alpha 0.1 --out-dir figures
```

This scans `data/evalues` by default (override with `--data-dir` or
`--col-name`), prints the table to stdout, and writes a CSV to
`figures/real/real_results_alpha_0.1.csv`. Repeat with another `--alpha` to
produce additional tables.

- p-value datasets for closed BY (APSAC, NAEP, PADJUST, PVALUES, VANDEVIJVER,
GOLUB): run from `r/eclosure/R` so the helper functions are found, e.g.

```bash
cd r/eclosure/R
Rscript ../inst/scripts/examples_BYbar.R > ../../figures/by_real_data_table.tex
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

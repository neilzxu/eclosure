# Scripts

The supported command-line workflow for the Python package is the packaged
module entrypoint:

```bash
python -m eclosure.experiments --help
```

Install the experiment dependencies first:

```bash
python -m pip install "eclosure[experiments]"
```

## Simulation mode

Run the simulation grid used by the paper:

```bash
python -m eclosure.experiments --mode simulation --alpha 0.1 --out-dir figures
```

Useful options:

- `--sim-config`: choose named simulation configurations, or use `all`.
- `--n-trials`: set the number of trials per grid point.
- `--n-cores`: control parallel worker count.
- `--out-dir`: choose the directory where generated outputs are written.

## Real-data mode

Analyze packaged or custom e-value CSV files:

```bash
python -m eclosure.experiments --mode real --alpha 0.1 --out-dir figures
```

Useful options:

- `--data-dir`: point to a directory containing `*_evalues.csv` files.
- `--col-name`: choose the e-value column name in those CSV files.
- `--alpha`: set the target FDR level.

The real-data workflow writes results under `OUT_DIR/real/` and also prints a
summary table to stdout.

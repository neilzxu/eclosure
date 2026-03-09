# Tests

The Python test suite is organized under `python/tests/` and split into fast
unit coverage, higher-level behavioral checks, and opt-in timing benchmarks.

## Install test dependencies

```bash
python -m pip install -e ".[test]"
```

## Standard test run

Run the default Python tests with:

```bash
python -m pytest python/tests
```

This includes:

- `test_cebh.py`: closed e-BH boundary cases and invariants.
- `test_cby.py`: closed BY boundary cases and invariants.
- `test_closed_methods.py`: higher-level checks that compare C++ bindings with
  Python helpers and packaged real-data inputs.

## Timing benchmarks

Timing tests are intentionally excluded from the default run. Enable them with
the `ECLOSURE_RUN_TIMING_TESTS=1` guard:

```bash
ECLOSURE_RUN_TIMING_TESTS=1 python -m pytest -s python/tests/test_timing.py
```

These tests print timing information for large synthetic inputs and are meant
for performance checks rather than routine CI.

## Wheel smoke tests

Release CI installs each built wheel and runs:

- `python/tests/test_cebh.py`
- `python/tests/test_cby.py`
- `python/tests/test_closed_methods.py`

That keeps wheel validation fast while still checking the compiled extension,
optional test dependencies, and packaged datasets.

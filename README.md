# Pergamon

Pergamon is a population-analysis and visualization package for astrophysical samples. Its scientific role is to take a population dictionary of sample properties, summarize the key features, visualize the one- and two-dimensional marginal distributions, and diagnose how those populations compare to one another.

## Scientific purpose

The package is intended for population-level analysis of astrophysical systems, including comparisons across samples, feature extraction, and diagnostic visualization of multi-dimensional parameter distributions. It is especially suited to checking whether samples differ in physically meaningful ways and produces summary plots that make the underlying population structure visible.

## Current ecosystem role

Within the active scientific stack, Pergamon sits in the population-analysis and visualization layer. It depends on the shared numerical and plotting utilities in `tdpy`, uses the same scientific-library conventions as the broader workflow stack, and provides a workflow layer for summarizing and comparing astrophysical populations rather than implementing low-level numerical primitives.

## Installation

```bash
cd pergamon
python -m pip install -e .
export PERGAMON_PATH=/path/to/pergamon
```

`PERGAMON_PATH` identifies the repository root. Keep runtime inputs in `data/` and generated pipeline outputs in `visuals/`; both directories are ignored by Git.

## Minimal workflow

The package is designed around a dictionary of population features and a workflow entry point, `pergamon.init(...)`, which performs the analysis and optional plotting. A lightweight smoke test is sufficient to confirm the package loads and exposes the expected functionality:

```python
import pergamon
print(hasattr(pergamon, '__file__'))
print(hasattr(pergamon, 'init'))
```

## Main modules

- `pergamon/main.py`: the core population-analysis workflow and plotting logic.
- `examples/`: notebooks or example-driving files for population visualization.
- `tests/`: regression and import checks.

## Dependencies

The package uses the standard scientific stack and ecosystem utilities, including:

- `numpy`
- `pandas`
- `matplotlib`
- `tdpy`
- `ephesos`
- `miletos`
- `nicomedia`

## Output behavior

The package produces visual summaries of population properties and feature distributions. The important scientific objective is to make the population-level patterns visible to the researcher without forcing them to read the underlying implementation.

## Development status

Pergamon is maintained as a focused population-analysis workflow rather than a monolithic general-purpose repository. It is useful when used with the appropriate population dictionaries and feature conventions, and it should continue to rely on shared utility layers rather than duplicating low-level numerical functionality.


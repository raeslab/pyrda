# PyRDA: Python Distance-Based Redundancy Analysis

A Python implementation of Distance-Based Redundancy Analysis (db-RDA) as practiced in the RaesLab ([source code](https://github.com/raeslab/pyrda)). This package provides a streamlined workflow for conducting db-RDA analysis on microbial datasets to identify and visualize the impact of metadata variables on community composition.

## Features

- **Same results as RLdbRDA**: R², adjusted R², F, AIC and the forward selection match
  `vegan::capscale`/`anova.cca`/`ordiR2step` to ~1e-14; with `test="permutation"` and the
  same seed the p-values are identical to R's
- **All vegan distances**: every `vegan::vegdist` index (except `robust.aitchison`), same names
- **Statistical Testing**: by default a permutation-free analytic test that gives R's
  permutation p-values without permuting (fast and deterministic); R's permutation tests
  are available too. Benjamini-Hochberg correction as in R
- **Forward Selection**: Stepwise model building to identify non-redundant effects
- **Visualization**: Create publication-ready plots showing individual and cumulative effect sizes
- **Validation**: Comprehensive input validation and error checking

## Development environment (dev container)

Stack: Python 3.13 (uv), R 4.5 with vegan + [RLdbRDA](https://github.com/raeslab/RLdbRDA) preinstalled, Claude Code, VS Code dev container.

Requirements: Docker Desktop, and VS Code with the **Dev Containers** extension (`ms-vscode-remote.remote-containers`).

1. Clone the repository and open the folder in VS Code.
2. Run `Dev Containers: Reopen in Container` from the Command Palette (`Ctrl+Shift+P`).
3. The first build takes a few minutes; `uv sync` runs automatically.

Key commands inside the container:

```bash
uv run python example.py                          # Python implementation
uv run pytest                                     # tests
uv run ruff check . && uv run ruff format .       # lint / format
Rscript -e 'library(RLdbRDA); ?rldbrda'           # R reference implementation
claude                                            # Claude Code
```

## Installation

### From GitHub

```bash
uv add git+https://github.com/raeslab/pyrda        # or: pip install git+https://github.com/raeslab/pyrda
```

### From Source

```bash
git clone https://github.com/raeslab/pyrda
cd pyrda
uv sync          # or: pip install -e .
```

### Dependencies

- Python ≥ 3.10
- numpy ≥ 1.20.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0
- matplotlib ≥ 3.5.0

## Quick Start

```python
import pandas as pd
from pyrda import rldbrda, prepare_plot_data, plot_dbrda

# Load your data
species_data = pd.read_csv("species_abundance.csv", index_col=0)
environmental_data = pd.read_csv("environmental_variables.csv", index_col=0)

# Perform db-RDA analysis (analytic p-values, no permutations)
results = rldbrda(species_data, environmental_data, method="bray", p_cutoff=0.05)

# Exactly the same output as set.seed(1); RLdbRDA::rldbrda(...) in R
results_r = rldbrda(species_data, environmental_data, test="permutation", seed=1)

# Prepare and create visualization
plot_data = prepare_plot_data(results)
fig = plot_dbrda(plot_data)
fig.show()
```

## Usage

### Basic Analysis

The main function `rldbrda()` performs the complete analysis workflow:

```python
from pyrda import rldbrda

results = rldbrda(
    df=species_data,  # Species abundance matrix
    meta=environmental_data,  # Environmental metadata
    method="bray",  # Distance metric (vegan::vegdist name)
    p_cutoff=0.05,  # Significance threshold
    test="analytic",  # default; "permutation" runs R's permutation tests
    seed=None,  # with test="permutation": seed=1 reproduces set.seed(1) in R
)
```

### Parameters

- **df**: pandas.DataFrame with samples as rows and species as columns
- **meta**: pandas.DataFrame with samples as rows and environmental variables as columns
- **method**: any `vegan::vegdist` method ("bray", "jaccard", "euclidean", "hellinger",
  "gower", ...; scipy's "braycurtis" is accepted too). Other names go to
  `scipy.spatial.distance.pdist` and have no R equivalent.
- **p_cutoff**: Significance threshold for variable selection (default: 0.05)
- **test**: how p-values are computed
  - `"analytic"` (default): no permutations. Uses the exact mean, variance and skewness
    of the permutation distribution (Pearson type III approximation); p-values agree with
    R's permutation test within its own Monte Carlo error and are deterministic. When a
    test's permutation distribution is strongly skewed (a few samples dominate the
    variation, e.g. Euclidean distance on raw counts) the approximation can be off; a
    warning then names the affected variables so you can check them with
    `test="permutation"`
  - `"permutation"`: vegan's permutation tests, 999 permutations per variable and 9999
    per forward-selection step, exactly as RLdbRDA
- **permutations** / **cumul_permutations**: number of permutations (default 999 / 9999);
  the analytic test reports its p-values on the same scale (smallest value 1/1000)
- **seed**: seed of the R-compatible random generator used by `test="permutation"`;
  `seed=s` gives the same permutations as `set.seed(s)` in R

### Results Interpretation

The results DataFrame contains:

- **r2**: Individual R² values for each variable
- **r2adj**: Adjusted R² values
- **pval**: p-values of the (analytic or permutation) pseudo-F test
- **padj**: Benjamini-Hochberg adjusted p-values
- **RDAcumul_R2.adj**: Cumulative R² from forward selection (for significant variables);
  the `<All variables>` row holds the adjusted R² of the model with all significant variables
- **RDAcumul_Df / RDAcumul_AIC / RDAcumul_F / RDAcumul_pval**: test of each forward-selection
  step, as reported by `ordiR2step`

### Visualization

```python
from pyrda import prepare_plot_data, plot_dbrda

# Prepare data for plotting
plot_data = prepare_plot_data(results)

# Create horizontal bar plot
fig = plot_dbrda(plot_data, figsize=(10, 8))

# Customize colors if desired
fig = plot_dbrda(
    plot_data,
    colors={
        "r2adj": "#60A68B",  # Individual R²
        "RDAcumul_R2.adj": "#1F4068",  # Cumulative R²
    },
)
```

## Example

See `example.py` for a complete working example using the varespec/varechem datasets.

```bash
uv run python example.py
```

## Data Requirements

### Input Data Format

1. **Species/Community Data**: Rows = samples, columns = species/taxa
2. **Environmental Metadata**: Rows = samples, columns = environmental variables

### Requirements

- Both datasets must have identical row indices (sample names)
- No missing values in environmental metadata
- Column names must be valid Python identifiers

### Data Validation

The package includes comprehensive validation:

```python
from pyrda.validation import validate_input

# This will raise informative errors if data doesn't meet requirements
validate_input(species_data, environmental_data)
```

## Methodology

PyRDA implements a two-step db-RDA procedure:

1. **Individual Effects**: Calculate R² for each environmental variable using constrained ordination
2. **Cumulative Effects**: Forward stepwise selection to identify non-redundant variable contributions

### Statistical Approach

- **Distance Matrix**: port of `vegan::vegdist`
- **Constrained Ordination**: as `vegan::capscale`: principal coordinates on all axes with
  positive eigenvalues, then redundancy analysis; negative eigenvalues are left out of R²
- **Significance Testing**: pseudo-F as `anova.cca`; p-values from the exact moments of the
  permutation distribution (default) or from the permutations themselves
- **Multiple Testing**: Benjamini-Hochberg false discovery rate correction
- **Forward Selection**: as `vegan::ordiR2step(Pin = 1, R2scope = TRUE)`: add the variable
  with the largest adjusted R² while it increases and stays at or below the adjusted R² of
  the model with all significant variables

## Comparison with R Implementation

PyRDA reproduces RLdbRDA. Statistics agree to ~1e-14. With `test="permutation", seed=s`
the p-values are identical to `set.seed(s); rldbrda(...)` in R; with the default analytic
test they agree within R's own run-to-run variation, and the selected variables match R's
most frequent result.

### Known differences

- **Ties in forward selection**: when adding a variable gives exactly the adjusted R² of
  the full model, R decides `R2 <= R2.all` by floating-point rounding, so its selection can
  change with the column order of `meta`. PyRDA treats equal values as equal and adds the
  variable.
- **Crash cases**: where R stops with an error (a constant metadata column, or no variable
  improving the model) PyRDA warns and continues.
- **Row order**: R sorts on R² as text, which misplaces variables with R² < 1e-4; PyRDA
  sorts numerically.
- **`robust.aitchison`** is not available (vegan imputes zeros with OptSpace).

## Testing

Run the test suite:

```bash
uv run pytest tests/ -v
```

The tests include:

- Unit tests for individual functions
- Integration tests with example data
- Exact comparison with RLdbRDA output (`tests/data/r_reference`, regenerate with
  `Rscript generate.R` in that folder)
- Edge case validation

## Development

### Code Quality

The package uses ruff for code formatting and linting:

```bash
uv run ruff check .
uv run ruff format .
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Run ruff for code formatting
6. Submit a pull request

## Support

- **Issues**: Report bugs and request features on [GitHub Issues](https://github.com/raeslab/pyrda/issues)
- **Documentation**: See docstrings in the source code
- **Examples**: Check the `example.py` script for usage patterns

## Acknowledgments


This Python implementation is based on the original R package RLdbRDA developed by Sebastian Proost for the RaesLab. The statistical methodology follows established practices in microbial ecology and community analysis.

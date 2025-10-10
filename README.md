# PyRDA: Python Distance-Based Redundancy Analysis

A Python implementation of Distance-Based Redundancy Analysis (db-RDA) as practiced in the RaesLab. This package provides a streamlined workflow for conducting db-RDA analysis on microbial datasets to identify and visualize the impact of metadata variables on community composition.

## Features

- **Distance-Based Redundancy Analysis**: Perform constrained ordination to identify environmental drivers of community composition
- **Multiple Distance Metrics**: Support for Bray-Curtis, Euclidean, Manhattan, and other distance measures
- **Statistical Testing**: Permutation-based significance testing with multiple testing correction
- **Forward Selection**: Stepwise model building to identify non-redundant effects
- **Visualization**: Create publication-ready plots showing individual and cumulative effect sizes
- **Validation**: Comprehensive input validation and error checking

## Installation

### From Source

```bash
git clone https://github.com/raeslab/RLdbRDA.git
cd RLdbRDA/pyrda
pip install -e .
```

### Dependencies

- Python ≥ 3.8
- numpy ≥ 1.20.0
- pandas ≥ 1.3.0
- scipy ≥ 1.7.0
- scikit-learn ≥ 1.0.0
- matplotlib ≥ 3.5.0

## Quick Start

```python
import pandas as pd
from pyrda import rldbrda, prepare_plot_data, plot_dbrda

# Load your data
species_data = pd.read_csv("species_abundance.csv", index_col=0)
environmental_data = pd.read_csv("environmental_variables.csv", index_col=0)

# Perform db-RDA analysis
results = rldbrda(species_data, environmental_data, method="braycurtis", p_cutoff=0.05)

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
    df=species_data,           # Species abundance matrix
    meta=environmental_data,   # Environmental metadata
    method="braycurtis",       # Distance metric
    p_cutoff=0.05             # Significance threshold
)
```

### Parameters

- **df**: pandas.DataFrame with samples as rows and species as columns
- **meta**: pandas.DataFrame with samples as rows and environmental variables as columns
- **method**: Distance metric ("braycurtis", "euclidean", "manhattan", etc.)
- **p_cutoff**: Significance threshold for variable selection (default: 0.05)

### Results Interpretation

The results DataFrame contains:

- **r2**: Individual R² values for each variable
- **r2adj**: Adjusted R² values
- **pval**: Raw p-values from permutation tests
- **padj**: Benjamini-Hochberg adjusted p-values
- **RDAcumul_R2.adj**: Cumulative R² from forward selection (for significant variables)

### Visualization

```python
from pyrda import prepare_plot_data, plot_dbrda

# Prepare data for plotting
plot_data = prepare_plot_data(results)

# Create horizontal bar plot
fig = plot_dbrda(plot_data, figsize=(10, 8))

# Customize colors if desired
fig = plot_dbrda(plot_data, colors={
    "r2adj": "#60A68B",           # Individual R²
    "RDAcumul_R2.adj": "#1F4068"  # Cumulative R²
})
```

## Example

See `example.py` for a complete working example using the varespec/varechem datasets.

```bash
cd pyrda
python example.py
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

- **Distance Matrix**: Computed using scipy.spatial.distance
- **Constrained Ordination**: Approximated using Principal Coordinates Analysis + linear regression
- **Significance Testing**: Permutation-based F-tests
- **Multiple Testing**: Benjamini-Hochberg false discovery rate correction
- **Forward Selection**: Greedy stepwise addition based on R² improvement

## Comparison with R Implementation

PyRDA is designed to produce results comparable to the original R package RLdbRDA:

- **Same statistical methodology**: db-RDA with forward selection
- **Compatible output format**: Similar result structure and column names
- **Equivalent visualizations**: Matching plot aesthetics and interpretation

### Key Differences

- **Implementation**: Python/scikit-learn vs R/vegan
- **Ordination**: Classical MDS approximation vs direct capscale
- **Dependencies**: Pure Python ecosystem vs R/Bioconductor

## Testing

Run the test suite:

```bash
cd pyrda
pytest tests/ -v
```

The tests include:

- Unit tests for individual functions
- Integration tests with example data
- Comparison with expected R outputs (when available)
- Edge case validation

## Development

### Code Quality

The package uses ruff for code formatting and linting:

```bash
pip install ruff
ruff check .
ruff format .
```

### Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Run ruff for code formatting
6. Submit a pull request

## License

MIT License - see LICENSE file for details.

## Citation

If you use PyRDA in your research, please cite:

```bibtex
@software{pyrda2024,
  title={PyRDA: Python Distance-Based Redundancy Analysis},
  author={Claude Code},
  year={2024},
  url={https://github.com/raeslab/RLdbRDA}
}
```

## Support

- **Issues**: Report bugs and request features on GitHub Issues
- **Documentation**: See docstrings in the source code
- **Examples**: Check the `example.py` script for usage patterns

## Acknowledgments

This Python implementation is based on the original R package RLdbRDA developed by Sebastian Proost for the RaesLab. The statistical methodology follows established practices in microbial ecology and community analysis.
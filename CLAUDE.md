# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Project Overview

RLdbRDA is an R package for Distance-Based Redundancy Analysis (db-RDA) as practiced in the RaesLab. The package provides a streamlined workflow for conducting db-RDA analysis on microbial datasets to identify and visualize the impact of metadata variables on community composition.

## Development Commands

### Package Development
- `devtools::load_all()` - Load the package for development/testing
- `devtools::document()` - Generate documentation from roxygen2 comments
- `devtools::check()` - Run R CMD check on the package
- `devtools::test()` - Run all tests using testthat
- `R CMD build .` - Build the package tarball
- `R CMD check --as-cran .` - Full CRAN-level package check

### Testing
- `devtools::test()` - Run all tests in tests/testthat/
- `testthat::test_file("tests/testthat/test-rldbrda.R")` - Run specific test file

### Installation
- `devtools::install()` - Install the package locally
- `renv::install("raeslab/RLdbRDA")` - Install from GitHub using renv

## Package Architecture

The package implements a two-step db-RDA procedure:

### Core Functions (in R/rldbrda.R):

1. **Input Validation**: `validate_input()` - Validates data and metadata compatibility
2. **Individual Effects**: `get_r2_single()` and `get_r2()` - Calculate R² for each metadata variable using `vegan::capscale`
3. **Cumulative Effects**: `get_cumul()` - Forward stepwise db-RDA using `vegan::ordiR2step`
4. **Main Function**: `rldbrda()` - Orchestrates the full analysis workflow
5. **Visualization**: `prepare_plot_data()` and `plot_dbrda()` - Prepare and create effect size plots

### Workflow:
1. Validate input data and metadata (matching rownames, no missing values, valid column names)
2. Compute distance matrix using `vegan::vegdist` (default: Bray-Curtis)
3. Calculate individual R² values for each metadata variable
4. Filter significant variables (default p < 0.05, Benjamini-Hochberg adjusted)
5. Perform forward stepwise db-RDA on significant variables
6. Combine individual and cumulative results for visualization

### Dependencies
- Core: vegan, dplyr, ggplot2, tibble, tidyr, magrittr
- Suggests: testthat (≥ 3.0.0)

### Key Constraints
- Metadata cannot contain missing values
- Column names must be valid R names (use `make.names()` if needed)
- Data and metadata must have identical rownames
- Requires at least one significant variable to proceed with cumulative analysis

## Testing Strategy

Tests use the `varespec` and `varechem` datasets from the vegan package as test data. All major functions are tested for basic functionality and error conditions.
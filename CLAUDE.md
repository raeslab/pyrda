# CLAUDE.md

Guidance for Claude Code when working in this repository.

## Project Overview

pyrda is a Python port of the RaesLab R package [RLdbRDA](https://github.com/raeslab/RLdbRDA):
Distance-Based Redundancy Analysis (db-RDA) to quantify how metadata variables explain
microbial community composition. A main goal is to **match RLdbRDA's results**; the R
package is installed in the dev container as the reference implementation.

## Environment

Everything runs in the VS Code dev container (root `Dockerfile`, `.devcontainer/`):

- Python 3.13 managed by **uv** (never pip/poetry/conda). Venv: `/home/vscode/.venv`.
- R 4.5 (rocker/r-ver) with vegan and RLdbRDA pinned to a commit (`RLDBRDA_REF` in the Dockerfile).
- Dependencies live in `pyproject.toml` (dev tools in `[dependency-groups] dev`); `uv.lock` is committed.

## Commands

- `uv sync` - install/update the environment
- `uv add <pkg>` / `uv add --dev <pkg>` - add a dependency
- `uv run pytest` - run all tests (`uv run pytest tests/test_core.py -k name` for one)
- `uv run ruff check .` and `uv run ruff format .` - lint and format
- `uv run python example.py` - run the varespec/varechem example
- `Rscript -e 'library(RLdbRDA); ...'` - run the R reference implementation

## Package Architecture (`src/pyrda/`)

- `validation.py` - `validate_input()` (matching indices, no missing metadata, Python
  identifiers as column names), `make_valid_names()`. Deliberately *not* R's `make.names()`.
- `distance.py` - `vegdist()`: port of `vegan::vegdist` (all indices except
  `robust.aitchison`, same names, R's partial matching; scipy names as fallback)
- `rrandom.py` - `RRandom`: R's Mersenne-Twister + `sample.int` + `permute::shuffleSet`,
  bit-for-bit (`RRandom(s)` == `set.seed(s)`)
- `moments.py` - exact permutation moments of `tr(A P G P')` + Pearson type III
  (`test="analytic"`); graph contractions in closed form (O(n^2 r), never einsum on Y).
  Skewness > `_MAX_SKEW` only *warns* (names the variables); never auto-switch to
  permutations - the user decides
- `core.py`
  - `_pcoa()` - capscale's PCoA (all positive axes; negative eigenvalues kept apart,
    distances / sqrt(n-1) when max >= 4)
  - `_fit()` / `_orthonormal_basis()` - RDA of the axes, R's qr() aliasing rule (tol 1e-7)
  - `_test_model()` / `_test_nested()` - `anova.cca` / `anovaCCAlist` (analytic = default,
    permutation). The parametric F test was removed on purpose (3-30x too conservative).
  - `get_r2_single()` / `get_r2()` - per-variable R², adjusted R², F, p, BH (`p.adjust`)
  - `get_cumulative()` - `ordiR2step(Pin = 1, R2scope = TRUE)`, incl. `<All variables>` row
  - `rldbrda()` - full workflow; `prepare_plot_data()` - long format for plotting
- `plotting.py` - `plot_dbrda()`, `save_plot()`, `create_summary_plot()` (matplotlib)

### Workflow

1. Validate data and metadata
2. Compute distance matrix
3. Individual R² per metadata variable
4. Keep significant variables (default padj < 0.05, Benjamini-Hochberg)
5. Forward stepwise selection on those variables (cumulative R²)
6. Merge individual and cumulative results for plotting

### Key Constraints

- Metadata cannot contain missing values
- Column names must be valid Python identifiers (use `make_valid_names()`)
- Data and metadata must share identical row indices
- Cumulative analysis needs at least one significant variable

## Comparing with R

- Only the **statistics** must match R (R², F, p, AIC, selection); R's naming/validation
  conventions are not ported.
- Same seed => identical p-values: `rldbrda(..., test="permutation", seed=s)` ==
  `set.seed(s); rldbrda(...)`. The default `test="analytic"` matches within Monte Carlo error.
  R's RNG is consumed in the same order (999 perms per variable, 9999 per selection step).
- Exact regression tests: `tests/test_r_reference.py` against `tests/data/r_reference/`
  (regenerate with `Rscript generate.R` in that folder). Tolerance 1e-12 (R prints 15 digits).
- `data/example_output/` is an *unseeded* R run: p-values and the selection there are one
  random draw (R gives 5 different selections over 200 seeds for varespec).
- Known, intended differences: forward-selection ties with `<All variables>` (R decides by
  rounding, pyrda uses tolerance `_TIE_TOL`), crash cases in R, numeric vs text sorting of r2.
  The README's "Known differences" section lists them for users.

## Conventions

- Ruff, line length 88, `requires-python >= 3.10`
- Tests in `tests/` with pytest; add tests for new functionality

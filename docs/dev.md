# Developing pyrda

## Development environment (dev container)

Stack: Python 3.13 (uv), R 4.5 with vegan + [RLdbRDA](https://github.com/raeslab/rldbrda)
preinstalled, Claude Code, VS Code dev container. R and RLdbRDA are there as the reference
implementation that pyrda must match; RLdbRDA is pinned to a commit (`RLDBRDA_REF` in the
[Dockerfile](../Dockerfile)).

Requirements: Docker Desktop, and VS Code with the **Dev Containers** extension
(`ms-vscode-remote.remote-containers`).

1. Clone the repository and open the folder in VS Code.
2. Run `Dev Containers: Reopen in Container` from the Command Palette (`Ctrl+Shift+P`).
3. The first build takes a few minutes; `uv sync` runs automatically.

Key commands inside the container:

```bash
uv run python docs/example.py                     # worked example (rewrites docs/output/)
uv run pytest                                     # tests
uv run ruff check . && uv run ruff format .       # lint / format
Rscript -e 'library(RLdbRDA); ?rldbrda'           # R reference implementation
claude                                            # Claude Code
```

## Testing

Run the test suite:

```bash
uv run pytest tests/ -v
```

The tests include:

- Unit tests for individual functions
- Integration tests with the example data in [`docs/data`](data)
- Exact comparison with RLdbRDA output
  ([`tests/test_r_reference.py`](../tests/test_r_reference.py), relative tolerance 1e-12)
- Edge case validation

### R reference data

The files in [`tests/data/r_reference`](../tests/data/r_reference) are frozen RLdbRDA
output; the repository holds no R code to regenerate them. They were made in the dev
container (R 4.5.3, vegan 2.7.3, RLdbRDA at `RLDBRDA_REF`), with `set.seed()` before every
`rldbrda()` call so that pyrda's `test="permutation", seed=...` must reproduce them
exactly:

| Files | Input | R call |
|---|---|---|
| `varespec_seed1.tsv`, `varespec_seed42.tsv` | [`docs/data`](data) | `rldbrda()` with defaults, seeds 1 and 42 |
| `mixed_{bray,jaccard,euclidean}_seed7.tsv`, `collinear_{bray,jaccard}_seed7.tsv` | `mixed_*.tsv`, `collinear_*.tsv` | `rldbrda(method = ..., p_cutoff = 0.5)`, seed 7 (no collinear/euclidean file: R stops with an error there) |
| `small_vegdist.tsv` | `small_counts.tsv` | `vegdist()` for 20 indices, one column per index |

The synthetic inputs (`mixed_*`, `collinear_*`, `small_counts.tsv`) are made by
[`make_inputs.py`](../tests/data/r_reference/make_inputs.py). Data files the tests read
(`docs/data`, `tests/data`) must not be edited.

## Docs

- [`docs/example.py`](example.py) writes [`docs/output`](output); the README's example
  table and figure come from it. Rerun it after changes that affect the results or the
  plot, and update the table in the README.
- [`docs/comparison.md`](comparison.md) is the validation report against RLdbRDA, with its
  figures in [`docs/img`](img). It is a snapshot of pyrda 0.1.0.

## Code Quality

The package uses ruff for code formatting and linting:

```bash
uv run ruff check .
uv run ruff format .
```

## Contributing

1. Fork the repository
2. Create a feature branch
3. Add tests for new functionality
4. Ensure all tests pass
5. Run ruff for code formatting
6. Submit a pull request

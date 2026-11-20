"""Compare pyrda with RLdbRDA output stored in tests/data/r_reference.

The reference files are frozen RLdbRDA output, made with ``set.seed()`` (see
docs/dev.md for how); with ``test="permutation"`` and the same seed the
p-values must match exactly, not just approximately. The default analytic test
is checked against the same files with a tolerance.
"""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest

from pyrda import rldbrda

ROOT = Path(__file__).parent.parent
REF = ROOT / "tests" / "data" / "r_reference"
STAT_COLUMNS = [
    "Fa",
    "r2",
    "r2adj",
    "N",
    "pval",
    "padj",
    "RDAcumul_R2.adj",
    "RDAcumul_Df",
    "RDAcumul_AIC",
    "RDAcumul_F",
    "RDAcumul_pval",
    "RDAcumul_N",
]


def read(path):
    return pd.read_csv(path, sep="\t", index_col=0)


def run_case(name):
    if name.startswith("varespec"):
        df, meta = (
            read(ROOT / "docs/data/varespec.tsv"),
            read(ROOT / "docs/data/varechem.tsv"),
        )
        seed = int(name.removeprefix("varespec_seed"))
        return rldbrda(df, meta, test="permutation", seed=seed)
    case, method, _ = name.split("_")
    df, meta = read(REF / f"{case}_data.tsv"), read(REF / f"{case}_meta.tsv")
    return rldbrda(df, meta, method=method, p_cutoff=0.5, test="permutation", seed=7)


def assert_same(py, r, rows):
    for col in STAT_COLUMNS:
        # R writes 15 significant digits
        np.testing.assert_allclose(
            py.loc[rows, col].astype(float),
            r.loc[rows, col].astype(float),
            rtol=1e-12,
            atol=1e-14,
            equal_nan=True,
            err_msg=col,
        )


@pytest.mark.parametrize(
    "name",
    [
        "varespec_seed1",
        "varespec_seed42",
        "mixed_euclidean_seed7",
        "mixed_jaccard_seed7",
        "collinear_bray_seed7",
        "collinear_jaccard_seed7",
    ],
)
def test_identical_to_r(name):
    r = read(REF / f"{name}.tsv")
    py = run_case(name)
    assert list(py.index) == list(r.index)
    assert_same(py, r, r.index)


def test_known_difference_scope_tie():
    """R drops the last variable when its R2 ties with <All variables>.

    In ``mixed_bray`` adding ``moist`` gives the full model, so its adjusted
    R-squared equals that of ``<All variables>``. R's ``<=`` comparison is then
    decided by floating-point rounding (in R the outcome changes with the column
    order of ``meta``); pyrda treats equal as equal.
    Everything else is identical.
    """
    r = read(REF / "mixed_bray_seed7.tsv")
    py = run_case("mixed_bray_seed7")
    assert np.isnan(r.loc["moist", "RDAcumul_R2.adj"])
    assert py.loc["moist", "RDAcumul_R2.adj"] == pytest.approx(
        r.loc["<All variables>", "RDAcumul_R2.adj"], abs=1e-12
    )
    common = [i for i in r.index if i != "moist"]
    assert_same(py, r, common)
    cols = ["Fa", "r2", "r2adj", "pval", "padj"]
    np.testing.assert_allclose(
        py.loc["moist", cols].astype(float),
        r.loc["moist", cols].astype(float),
        rtol=1e-12,
    )


@pytest.mark.filterwarnings("ignore:.*analytic p-value")
@pytest.mark.parametrize("name", ["varespec_seed1", "mixed_jaccard_seed7"])
def test_analytic_default_close_to_r(name):
    """Default test: same statistics, p-values within R's Monte Carlo error."""
    r = read(REF / f"{name}.tsv")
    if name.startswith("varespec"):
        df, meta = (
            read(ROOT / "docs/data/varespec.tsv"),
            read(ROOT / "docs/data/varechem.tsv"),
        )
        py = rldbrda(df, meta)
    else:
        case, method, _ = name.split("_")
        df, meta = read(REF / f"{case}_data.tsv"), read(REF / f"{case}_meta.tsv")
        py = rldbrda(df, meta, method=method, p_cutoff=0.5)
    feats = [i for i in r.index if i != "<All variables>"]
    for col in ["Fa", "r2", "r2adj"]:
        np.testing.assert_allclose(py.loc[feats, col], r.loc[feats, col], rtol=1e-12)
    # 999 permutations: standard error at most ~0.016; allow 4 SE
    np.testing.assert_allclose(py.loc[feats, "pval"], r.loc[feats, "pval"], atol=0.06)

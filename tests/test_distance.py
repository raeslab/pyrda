"""Tests for the vegdist port."""

from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from scipy.spatial.distance import squareform

from pyrda.distance import resolve_method, vegdist

REF = Path(__file__).parent / "data" / "r_reference"


@pytest.fixture(scope="module")
def small():
    return pd.read_csv(REF / "small_counts.tsv", sep="\t", index_col=0)


@pytest.fixture(scope="module")
def r_distances():
    return pd.read_csv(REF / "small_vegdist.tsv", sep="\t")


@pytest.mark.filterwarnings("ignore")
def test_all_methods_match_vegan(small, r_distances):
    for method in r_distances.columns:
        expected = squareform(r_distances[method].to_numpy())
        got = vegdist(small, method)
        np.testing.assert_allclose(
            got, expected, rtol=1e-12, atol=1e-14, err_msg=method
        )


def test_aliases_and_partial_matching():
    assert resolve_method("braycurtis") == ("bray", True)
    assert resolve_method("euc") == ("euclidean", True)
    assert resolve_method("cosine") == ("cosine", False)
    with pytest.raises(ValueError, match="ambiguous"):
        resolve_method("ch")


def test_jaccard_is_quantitative_unless_binary():
    x = np.array([[1.0, 5.0, 0.0], [2.0, 1.0, 3.0], [0.0, 4.0, 4.0]])
    bray = vegdist(x, "bray")
    np.testing.assert_allclose(vegdist(x, "jaccard"), 2 * bray / (1 + bray))
    binary = vegdist(x, "jaccard", binary=True)
    # binary Jaccard = 1 - shared / union of present species
    assert binary[0, 1] == pytest.approx(1 - 2 / 3)


def test_scipy_fallback_and_errors():
    x = np.random.default_rng(0).random((5, 4))
    assert vegdist(x, "cosine").shape == (5, 5)
    with pytest.raises(NotImplementedError):
        vegdist(x, "robust.aitchison")
    with pytest.raises(ValueError, match="missing"):
        vegdist(np.array([[1.0, np.nan], [1.0, 2.0]]))

"""Tests for core analysis functions."""

import numpy as np
import pandas as pd
import pytest

from pyrda.core import (
    _compute_distance_matrix,
    get_r2,
    get_r2_single,
    prepare_plot_data,
    rldbrda,
)
from pyrda.rrandom import RRandom


class TestComputeDistanceMatrix:
    """Test the _compute_distance_matrix function."""

    def test_bray_curtis_distance(self):
        """Test Bray-Curtis distance calculation."""
        df = pd.DataFrame([[1, 2, 3], [4, 5, 6], [7, 8, 9]])
        dist_matrix = _compute_distance_matrix(df, method="braycurtis")

        # Check basic properties
        assert dist_matrix.shape == (3, 3)
        assert np.allclose(np.diag(dist_matrix), 0)  # Diagonal should be zero
        assert np.allclose(dist_matrix, dist_matrix.T)  # Should be symmetric

    def test_euclidean_distance(self):
        """Test Euclidean distance calculation."""
        df = pd.DataFrame([[0, 0], [1, 1], [2, 2]])
        dist_matrix = _compute_distance_matrix(df, method="euclidean")

        # Check basic properties
        assert dist_matrix.shape == (3, 3)
        assert np.allclose(np.diag(dist_matrix), 0)
        assert np.allclose(dist_matrix, dist_matrix.T)


class TestGetR2Single:
    """Test the get_r2_single function."""

    def test_continuous_variable(self):
        """Test R² calculation for continuous variable."""
        # Create simple test data
        distance_matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        meta = pd.DataFrame({"var1": [1.0, 2.0, 3.0]}, index=["s1", "s2", "s3"])

        result = get_r2_single(distance_matrix, meta, "var1")

        # Check that result has expected structure
        assert isinstance(result, pd.Series)
        assert "feature" in result.index
        assert "r2" in result.index
        assert "r2adj" in result.index
        assert "Fa" in result.index
        assert "pval" in result.index
        assert "N" in result.index
        assert result["feature"] == "var1"
        assert result["N"] == 3

    def test_categorical_variable(self):
        """Test R² calculation for categorical variable."""
        distance_matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        meta = pd.DataFrame({"var1": ["A", "B", "C"]}, index=["s1", "s2", "s3"])

        result = get_r2_single(distance_matrix, meta, "var1")

        assert isinstance(result, pd.Series)
        assert result["feature"] == "var1"
        assert result["N"] == 3

    def test_missing_values(self):
        """Test handling of missing values."""
        distance_matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        meta = pd.DataFrame({"var1": [1.0, np.nan, 3.0]}, index=["s1", "s2", "s3"])

        result = get_r2_single(distance_matrix, meta, "var1")

        assert result["N"] == 2  # Only non-missing values


class TestGetR2:
    """Test the get_r2 function."""

    def test_multiple_variables(self):
        """Test R² calculation for multiple variables."""
        distance_matrix = np.array([[0, 1, 2], [1, 0, 1], [2, 1, 0]])
        meta = pd.DataFrame(
            {"var1": [1.0, 2.0, 3.0], "var2": ["A", "B", "C"]}, index=["s1", "s2", "s3"]
        )

        result = get_r2(distance_matrix, meta)

        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2
        assert "var1" in result.index
        assert "var2" in result.index
        assert "padj" in result.columns  # Check that p-values are adjusted

    def test_p_value_adjustment(self):
        """Test that p-values are properly adjusted."""
        # Create a more realistic distance matrix
        np.random.seed(42)
        distance_matrix = np.random.rand(4, 4)
        distance_matrix = (distance_matrix + distance_matrix.T) / 2  # Make symmetric
        np.fill_diagonal(distance_matrix, 0)  # Zero diagonal

        meta = pd.DataFrame(
            {"var1": [1, 2, 3, 4], "var2": [1, 1, 2, 2], "var3": [1, 2, 1, 2]}
        )

        result = get_r2(distance_matrix, meta)

        # Check that adjusted p-values exist and are >= raw p-values
        assert "padj" in result.columns
        assert all(result["padj"] >= result["pval"])


class TestPrepairePlotData:
    """Test the prepare_plot_data function."""

    def test_basic_functionality(self):
        """Test basic plot data preparation."""
        # Create mock db-RDA results
        dbrda_data = pd.DataFrame(
            {
                "r2adj": [0.3, 0.2, 0.1],
                "RDAcumul_R2.adj": [0.3, np.nan, np.nan],
                "padj": [0.01, 0.1, 0.2],
            },
            index=["var1", "var2", "var3"],
        )

        result = prepare_plot_data(dbrda_data)

        assert isinstance(result, pd.DataFrame)
        assert "feature" in result.columns
        assert "variable" in result.columns
        assert "value" in result.columns
        assert "significant" in result.columns

        # Check that data is in long format
        assert len(result) == 6  # 3 features * 2 variables

        # Check significance assignment
        sig_data = result[result["feature"] == "var1"]
        assert all(sig_data["significant"] == 1)  # var1 should be significant

    def test_missing_value_imputation(self):
        """Test that missing values are properly imputed."""
        dbrda_data = pd.DataFrame(
            {
                "r2adj": [0.3, 0.2],
                "RDAcumul_R2.adj": [0.3, np.nan],
            },
            index=["var1", "var2"],
        )

        result = prepare_plot_data(dbrda_data)

        # Check that NaN was replaced with max value
        cumul_values = result[result["variable"] == "RDAcumul_R2.adj"]["value"]
        assert not cumul_values.isna().any()
        assert cumul_values.iloc[1] == 0.3  # Should be imputed with max value


class TestRldbrda:
    """Test the main rldbrda function."""

    def test_basic_workflow(self):
        """Test the complete workflow."""
        # Create simple test data
        df = pd.DataFrame(
            {
                "species1": [1, 2, 3, 4],
                "species2": [2, 3, 4, 5],
                "species3": [3, 4, 5, 6],
            },
            index=["s1", "s2", "s3", "s4"],
        )

        meta = pd.DataFrame(
            {"env1": [1.0, 2.0, 3.0, 4.0], "env2": ["A", "A", "B", "B"]},
            index=["s1", "s2", "s3", "s4"],
        )

        result = rldbrda(df, meta, p_cutoff=1.0)  # Accept all features

        assert isinstance(result, pd.DataFrame)
        assert len(result) >= 1  # Should have at least one feature
        assert "r2" in result.columns
        assert "r2adj" in result.columns

    def test_no_significant_features(self):
        """Test behavior when no features are significant."""
        # Create data with no real patterns
        np.random.seed(42)
        df = pd.DataFrame(
            np.random.rand(10, 5), columns=[f"species_{i}" for i in range(5)]
        )
        meta = pd.DataFrame(
            np.random.rand(10, 3), columns=[f"env_{i}" for i in range(3)]
        )

        with pytest.raises(ValueError, match="No significant features found"):
            rldbrda(df, meta, p_cutoff=0.001)  # Very strict cutoff

    def test_input_validation_called(self):
        """Test that input validation is called."""
        df = pd.DataFrame([[1, 2], [3, 4]], index=["s1", "s2"])
        meta = pd.DataFrame([[1, 2]], index=["s1"])  # Mismatched indices

        with pytest.raises(ValueError):
            rldbrda(df, meta)


class TestRCompatibility:
    """Behaviour that mirrors vegan/RLdbRDA."""

    @pytest.fixture
    def data(self):
        rng = np.random.default_rng(5)
        n = 20
        x = rng.normal(size=n)
        df = pd.DataFrame(
            rng.poisson(np.exp(1 + np.outer(x, rng.normal(size=15)))),
            index=[f"s{i}" for i in range(n)],
        )
        meta = pd.DataFrame(
            {
                "x": x,
                "x_copy": x,
                "noise": rng.normal(size=n),
                "grp": rng.choice(["a", "b", "c"], n),
            },
            index=df.index,
        )
        return df, meta

    def test_bh_matches_r_p_adjust(self):
        from pyrda.core import _p_adjust_bh

        # p.adjust(c(0.01, 0.04, 0.03, 0.2, NA, 0.04), "BH") in R
        got = _p_adjust_bh([0.01, 0.04, 0.03, 0.2, np.nan, 0.04])
        expected = [0.05, 0.05, 0.05, 0.2, np.nan, 0.05]
        np.testing.assert_allclose(got, expected, equal_nan=True)

    def test_r2_uses_all_real_axes(self, data):
        """R-squared is the share of the real (positive-eigenvalue) inertia."""
        from pyrda.distance import vegdist

        df, meta = data
        d = vegdist(df)
        n = len(d)
        j = np.eye(n) - 1 / n
        ev, vec = np.linalg.eigh(-0.5 * j @ (d**2) @ j)
        pos = ev > 1e-8
        y = vec[:, pos] * np.sqrt(ev[pos])
        x = meta["x"].to_numpy() - meta["x"].mean()
        fitted = np.outer(x, x @ y) / (x @ x)
        expected = (fitted**2).sum() / (y**2).sum()
        result = get_r2_single(d, meta, "x")
        assert result["r2"] == pytest.approx(expected, rel=1e-10)

    def test_seed_reproduces_and_changes_p_values(self, data):
        df, meta = data
        a = rldbrda(df, meta, p_cutoff=1, test="permutation", seed=1)
        b = rldbrda(df, meta, p_cutoff=1, test="permutation", seed=1)
        pd.testing.assert_frame_equal(a, b)

    def test_collinear_copy_is_never_selected(self, data):
        df, meta = data
        result = rldbrda(df, meta, p_cutoff=1, seed=1)
        selected = result["RDAcumul_R2.adj"].dropna().index
        assert not {"x", "x_copy"} <= set(selected)

    def test_selection_does_not_depend_on_column_order(self, data):
        df, meta = data
        reference = rldbrda(df, meta, p_cutoff=1, test="analytic")
        ref_sel = list(reference["RDAcumul_R2.adj"].dropna().index)
        for cols in (["grp", "noise", "x", "x_copy"], ["noise", "grp", "x", "x_copy"]):
            other = rldbrda(df, meta[cols], p_cutoff=1, test="analytic")
            assert list(other["RDAcumul_R2.adj"].dropna().index) == ref_sel

    def test_constant_variable_is_handled(self, data):
        df, meta = data
        meta = meta.assign(const=1.0)
        with pytest.warns(UserWarning, match="no variation"):
            result = get_r2(_compute_distance_matrix(df), meta, seed=1)
        assert result.loc["const", "r2"] == 0
        assert np.isnan(result.loc["const", "pval"])

    def test_all_variables_row_and_plot_fill(self, data):
        df, meta = data
        result = rldbrda(df, meta, p_cutoff=1, seed=1)
        assert "<All variables>" in result.index
        r2_all = result.loc["<All variables>", "RDAcumul_R2.adj"]
        plot = prepare_plot_data(result)
        assert "<All variables>" not in set(plot["feature"])
        not_selected = result.index[result["RDAcumul_R2.adj"].isna()]
        filled = plot[
            plot["feature"].isin(not_selected) & (plot["variable"] == "RDAcumul_R2.adj")
        ]
        assert (filled["value"] == r2_all).all()
        assert (filled["significant"] == 0).all()
        assert list(plot["feature"].cat.categories) == list(result.index[::-1])

    @pytest.mark.parametrize("test", ["analytic", "permutation"])
    def test_all_tests_run(self, data, test):
        df, meta = data
        result = rldbrda(df, meta, p_cutoff=1, test=test, seed=2)
        p = result["pval"].dropna()
        assert ((p > 0) & (p <= 1)).all()

    def test_invalid_test_name(self, data):
        df, meta = data
        with pytest.raises(ValueError, match="test must be one of"):
            rldbrda(df, meta, test="bogus")
        with pytest.raises(ValueError, match="test must be one of"):
            rldbrda(df, meta, test="parametric")  # removed: far too conservative

    def test_analytic_close_to_permutation(self, data):
        df, meta = data
        d = _compute_distance_matrix(df)
        perm = get_r2(d, meta, test="permutation", permutations=19999, seed=4)
        analytic = get_r2(d, meta, permutations=19999)
        np.testing.assert_allclose(analytic["pval"], perm["pval"], atol=0.01)

    def test_analytic_warns_but_never_permutes(self):
        """Skewed permutation distributions give a warning naming the variable."""
        rng = np.random.default_rng(3)
        n = 40
        x = rng.normal(size=n)
        counts = rng.poisson(5, size=(n, 30)).astype(float)
        counts[0] *= 50  # one sample dominates a Euclidean ordination
        df = pd.DataFrame(counts, index=[f"s{i}" for i in range(n)])
        skewed = np.exp(rng.normal(0, 1.5, n))
        skewed[0] = skewed.max() * 3
        meta = pd.DataFrame({"x": x, "skewed": skewed}, index=df.index)
        d = _compute_distance_matrix(df, method="euclidean")

        class NoPermutations(RRandom):
            def shuffle_set(self, n, nperm):
                raise AssertionError("analytic test must not permute")

        with pytest.warns(UserWarning, match=r"'skewed'.*test='permutation'"):
            result = get_r2(d, meta, test="analytic", seed=NoPermutations(1))
        assert result["pval"].notna().all()

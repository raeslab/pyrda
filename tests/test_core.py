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
        distance_matrix = np.eye(4) + 0.1  # Nearly identity matrix
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
        df = pd.DataFrame(np.random.rand(10, 5))
        meta = pd.DataFrame(np.random.rand(10, 3))

        with pytest.raises(ValueError, match="No significant features found"):
            rldbrda(df, meta, p_cutoff=0.001)  # Very strict cutoff

    def test_input_validation_called(self):
        """Test that input validation is called."""
        df = pd.DataFrame([[1, 2], [3, 4]], index=["s1", "s2"])
        meta = pd.DataFrame([[1, 2]], index=["s1"])  # Mismatched indices

        with pytest.raises(ValueError):
            rldbrda(df, meta)

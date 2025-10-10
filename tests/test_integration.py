"""Integration tests using the example varespec/varechem data."""

import os

import numpy as np
import pandas as pd
import pytest

from pyrda import plot_dbrda, prepare_plot_data, rldbrda


class TestIntegrationWithExampleData:
    """Test the full workflow with the example data."""

    @pytest.fixture
    def varespec_data(self):
        """Load varespec test data."""
        data_path = os.path.join(
            os.path.dirname(__file__), "../../../data/varespec.tsv"
        )
        if os.path.exists(data_path):
            df = pd.read_csv(data_path, sep="\t", index_col=0)
            return df
        else:
            # Create mock data for CI environments
            np.random.seed(42)
            return pd.DataFrame(
                np.random.rand(24, 44),
                index=[str(i) for i in range(18, 42)],
                columns=[f"species_{i}" for i in range(44)],
            )

    @pytest.fixture
    def varechem_data(self):
        """Load varechem test data."""
        data_path = os.path.join(
            os.path.dirname(__file__), "../../../data/varechem.tsv"
        )
        if os.path.exists(data_path):
            df = pd.read_csv(data_path, sep="\t", index_col=0)
            return df
        else:
            # Create mock data for CI environments
            np.random.seed(42)
            return pd.DataFrame(
                np.random.rand(24, 14),
                index=[str(i) for i in range(18, 42)],
                columns=[
                    "N",
                    "P",
                    "K",
                    "Ca",
                    "Mg",
                    "S",
                    "Al",
                    "Fe",
                    "Mn",
                    "Zn",
                    "Mo",
                    "Baresoil",
                    "Humdepth",
                    "pH",
                ],
            )

    @pytest.fixture
    def expected_output(self):
        """Load expected output for comparison."""
        data_path = os.path.join(
            os.path.dirname(__file__), "../../../data/example_output/out.tsv"
        )
        if os.path.exists(data_path):
            df = pd.read_csv(data_path, sep="\t", index_col=0)
            return df
        else:
            return None

    def test_full_workflow(self, varespec_data, varechem_data):
        """Test the complete workflow with example data."""
        # Run the analysis
        result = rldbrda(
            varespec_data, varechem_data, method="braycurtis", p_cutoff=0.05
        )

        # Basic checks
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "r2" in result.columns
        assert "r2adj" in result.columns
        assert "pval" in result.columns
        assert "padj" in result.columns

        # Check that results are ordered properly
        assert result["r2"].is_monotonic_decreasing

        # Check for significant features
        significant_features = result[result["padj"] < 0.05]
        assert len(significant_features) > 0

    def test_plot_data_preparation(self, varespec_data, varechem_data):
        """Test plot data preparation."""
        result = rldbrda(varespec_data, varechem_data)
        plot_data = prepare_plot_data(result)

        assert isinstance(plot_data, pd.DataFrame)
        assert "feature" in plot_data.columns
        assert "variable" in plot_data.columns
        assert "value" in plot_data.columns
        assert "significant" in plot_data.columns

        # Check that we have both r2adj and cumulative data
        variables = plot_data["variable"].unique()
        assert "r2adj" in variables
        assert "RDAcumul_R2.adj" in variables

    def test_plotting_function(self, varespec_data, varechem_data):
        """Test that plotting function works without errors."""
        result = rldbrda(varespec_data, varechem_data)
        plot_data = prepare_plot_data(result)

        # This should not raise any errors
        fig = plot_dbrda(plot_data)
        assert fig is not None

    def test_comparison_with_expected_output(
        self, varespec_data, varechem_data, expected_output
    ):
        """Test comparison with R output if available."""
        if expected_output is None:
            pytest.skip("Expected output file not available")

        result = rldbrda(varespec_data, varechem_data)

        # Check that we get similar features as significant
        our_significant = set(result[result["padj"] < 0.05].index)
        r_significant = set(expected_output[expected_output["padj"] < 0.05].index)

        # Allow some differences due to implementation differences
        # But there should be substantial overlap
        overlap = len(our_significant.intersection(r_significant))
        total_unique = len(our_significant.union(r_significant))

        if total_unique > 0:
            overlap_ratio = overlap / total_unique
            assert overlap_ratio > 0.5, (
                f"Significant features overlap too low: {overlap_ratio}"
            )

    def test_reproducibility(self, varespec_data, varechem_data):
        """Test that results are reproducible."""
        # Run analysis twice
        result1 = rldbrda(varespec_data, varechem_data)
        result2 = rldbrda(varespec_data, varechem_data)

        # Results should be identical
        pd.testing.assert_frame_equal(result1, result2)

    def test_different_distance_methods(self, varespec_data, varechem_data):
        """Test with different distance methods."""
        methods = ["braycurtis", "euclidean", "manhattan"]

        for method in methods:
            result = rldbrda(varespec_data, varechem_data, method=method)
            assert isinstance(result, pd.DataFrame)
            assert len(result) > 0

    def test_different_p_cutoffs(self, varespec_data, varechem_data):
        """Test with different p-value cutoffs."""
        cutoffs = [0.01, 0.05, 0.1]

        results = {}
        for cutoff in cutoffs:
            try:
                result = rldbrda(varespec_data, varechem_data, p_cutoff=cutoff)
                results[cutoff] = result
            except ValueError:
                # No significant features - this is okay for very strict cutoffs
                continue

        # More lenient cutoffs should generally give more features
        if 0.01 in results and 0.1 in results:
            assert len(results[0.1]) >= len(results[0.01])

    def test_edge_cases(self, varechem_data):
        """Test edge cases."""
        # Test with minimal data
        small_varespec = pd.DataFrame(
            [[1, 0, 0], [0, 1, 0], [0, 0, 1]],
            index=varechem_data.index[:3],
            columns=["sp1", "sp2", "sp3"],
        )
        small_varechem = varechem_data.iloc[:3, :3]

        result = rldbrda(small_varespec, small_varechem, p_cutoff=1.0)
        assert isinstance(result, pd.DataFrame)

    def test_categorical_variables(self, varespec_data):
        """Test with categorical metadata variables."""
        # Create categorical metadata
        categorical_meta = pd.DataFrame(
            {
                "group": ["A"] * 12 + ["B"] * 12,
                "treatment": (["T1", "T2", "T3", "T4"] * 6)[:24],
            },
            index=varespec_data.index,
        )

        result = rldbrda(varespec_data, categorical_meta, p_cutoff=1.0)
        assert isinstance(result, pd.DataFrame)
        assert len(result) == 2

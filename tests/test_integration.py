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
        data_path = os.path.join(os.path.dirname(__file__), "../docs/data/varespec.tsv")
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
        data_path = os.path.join(os.path.dirname(__file__), "../docs/data/varechem.tsv")
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

    def test_full_workflow(self, varespec_data, varechem_data):
        """Test the complete workflow with example data."""
        # Run the analysis
        result = rldbrda(
            varespec_data, varechem_data, method="braycurtis", p_cutoff=0.2
        )

        # Basic checks
        assert isinstance(result, pd.DataFrame)
        assert len(result) > 0
        assert "r2" in result.columns
        assert "r2adj" in result.columns
        assert "pval" in result.columns
        assert "padj" in result.columns

        # Check for significant features
        significant_features = result[result["padj"] < 0.2]
        assert len(significant_features) > 0

    def test_plot_data_preparation(self, varespec_data, varechem_data):
        """Test plot data preparation."""
        result = rldbrda(varespec_data, varechem_data, p_cutoff=0.2)
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
        result = rldbrda(varespec_data, varechem_data, p_cutoff=0.2)
        plot_data = prepare_plot_data(result)

        # This should not raise any errors
        fig = plot_dbrda(plot_data)
        assert fig is not None

    def test_reproducibility(self, varespec_data, varechem_data):
        """Analytic results are deterministic; permutations need a seed."""
        result1 = rldbrda(varespec_data, varechem_data, p_cutoff=0.2)
        result2 = rldbrda(varespec_data, varechem_data, p_cutoff=0.2)
        pd.testing.assert_frame_equal(result1, result2)
        kw = {"p_cutoff": 0.2, "test": "permutation", "seed": 3}
        result1 = rldbrda(varespec_data, varechem_data, **kw)
        result2 = rldbrda(varespec_data, varechem_data, **kw)
        pd.testing.assert_frame_equal(result1, result2)

    def test_tests_agree_on_statistics(self, varespec_data, varechem_data):
        """The choice of test only changes p-values, not the effect sizes."""
        perm = rldbrda(varespec_data, varechem_data, test="permutation", seed=1)
        analytic = rldbrda(varespec_data, varechem_data)
        cols = ["Fa", "r2", "r2adj"]
        pd.testing.assert_frame_equal(perm.loc[analytic.index, cols], analytic[cols])

    def test_different_distance_methods(self, varespec_data, varechem_data):
        """Test with different distance methods."""
        methods = ["braycurtis", "euclidean"]

        for method in methods:
            result = rldbrda(varespec_data, varechem_data, method=method, p_cutoff=0.2)
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

        # with 3 samples every permutation p-value is large: nothing passes,
        # as in R
        with pytest.raises(ValueError, match="No significant features"):
            rldbrda(small_varespec, small_varechem, p_cutoff=0.5, seed=1)

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

        result = rldbrda(varespec_data, categorical_meta, p_cutoff=1.0, seed=1)
        assert isinstance(result, pd.DataFrame)
        assert set(result.index) == {"group", "treatment", "<All variables>"}
        # factors use k - 1 degrees of freedom, as in R
        n = len(varespec_data)
        r2 = result.loc["treatment", "r2"]
        expected = 1 - (1 - r2) * (n - 1) / (n - 3 - 1)
        assert result.loc["treatment", "r2adj"] == pytest.approx(expected)

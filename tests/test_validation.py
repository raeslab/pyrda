"""Tests for validation module."""

import numpy as np
import pandas as pd
import pytest

from pyrda.validation import make_valid_names, validate_input


class TestValidateInput:
    """Test the validate_input function."""

    def test_valid_input(self):
        """Test with valid input data."""
        df = pd.DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6]}, index=["sample1", "sample2", "sample3"]
        )
        meta = pd.DataFrame(
            {"var1": [1.0, 2.0, 3.0], "var2": ["a", "b", "c"]},
            index=["sample1", "sample2", "sample3"],
        )
        # Should not raise any exception
        validate_input(df, meta)

    def test_mismatched_indices(self):
        """Test with mismatched indices."""
        df = pd.DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6]}, index=["sample1", "sample2", "sample3"]
        )
        meta = pd.DataFrame(
            {"var1": [1.0, 2.0, 3.0]}, index=["sample1", "sample2", "sample4"]
        )
        with pytest.raises(ValueError, match="Data and metadata don't match"):
            validate_input(df, meta)

    def test_missing_values_in_metadata(self):
        """Test with missing values in metadata."""
        df = pd.DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6]}, index=["sample1", "sample2", "sample3"]
        )
        meta = pd.DataFrame(
            {"var1": [1.0, np.nan, 3.0]}, index=["sample1", "sample2", "sample3"]
        )
        with pytest.raises(ValueError, match="Missing values found in metadata"):
            validate_input(df, meta)

    def test_invalid_column_names(self):
        """Test with invalid column names."""
        df = pd.DataFrame(
            {"A": [1, 2, 3], "B": [4, 5, 6]}, index=["sample1", "sample2", "sample3"]
        )
        meta = pd.DataFrame(
            {"var-1": [1.0, 2.0, 3.0], "2var": [4.0, 5.0, 6.0]},
            index=["sample1", "sample2", "sample3"],
        )
        with pytest.raises(
            ValueError, match="Column names in metadata contain invalid characters"
        ):
            validate_input(df, meta)


class TestMakeValidNames:
    """Test the make_valid_names function."""

    def test_valid_names_unchanged(self):
        """Test that valid names remain unchanged."""
        names = ["var1", "var2", "valid_name"]
        result = make_valid_names(names)
        assert result == names

    def test_invalid_characters_replaced(self):
        """Test that invalid characters are replaced with underscores."""
        names = ["var-1", "var.2", "var@3"]
        expected = ["var_1", "var_2", "var_3"]
        result = make_valid_names(names)
        assert result == expected

    def test_names_starting_with_numbers(self):
        """Test that names starting with numbers get 'X' prefix."""
        names = ["1var", "2test"]
        expected = ["X1var", "X2test"]
        result = make_valid_names(names)
        assert result == expected

    def test_empty_names(self):
        """Test that empty names become 'X'."""
        names = ["", " ", "valid"]
        expected = ["X", "X", "valid"]
        result = make_valid_names(names)
        assert result == expected

"""Input validation functions for pyrda package."""


def validate_input(df, meta):
    """
    Validate input data and metadata for compatibility.

    This function checks if the input data frame matches the corresponding metadata
    in terms of row indices and also verifies whether the metadata contains any
    missing values. It ensures that the metadata column names are valid Python names.

    Parameters
    ----------
    df : pandas.DataFrame
        A data frame containing the data to be validated.
    meta : pandas.DataFrame
        A data frame containing metadata associated with the data. Both `df` and
        `meta` must have identical row indices.

    Raises
    ------
    ValueError
        If data and metadata don't match, if metadata contains missing values,
        or if column names contain invalid characters.
    """
    # Check if indices match
    if not df.index.equals(meta.index):
        raise ValueError(
            "Data and metadata don't match, should have the same row indices."
        )

    # Check for missing values in metadata
    if meta.isnull().any().any():
        raise ValueError(
            "Missing values found in metadata. Metadata cannot contain missing values."
        )

    # Check for valid column names (Python identifiers)
    invalid_cols = []
    for col in meta.columns:
        if not str(col).isidentifier():
            invalid_cols.append(col)

    if invalid_cols:
        raise ValueError(
            f"Column names in metadata contain invalid characters: {invalid_cols}. "
            "Use valid Python identifiers for column names."
        )


def make_valid_names(names):
    """
    Convert a list of names to valid Python identifiers.

    Parameters
    ----------
    names : list
        List of strings to convert to valid identifiers.

    Returns
    -------
    list
        List of valid Python identifiers.
    """
    valid_names = []
    for name in names:
        # Convert to string first
        name_str = str(name).strip()

        # If empty or only whitespace, use "X"
        if not name_str or name_str.isspace():
            valid_names.append("X")
            continue

        # Replace invalid characters with underscores
        valid_name = "".join(c if c.isalnum() else "_" for c in name_str)
        # Ensure it doesn't start with a number
        if valid_name and valid_name[0].isdigit():
            valid_name = "X" + valid_name
        # Ensure it's not empty after processing
        if not valid_name:
            valid_name = "X"
        valid_names.append(valid_name)
    return valid_names

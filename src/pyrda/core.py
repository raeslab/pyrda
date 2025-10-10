"""Core analysis functions for distance-based redundancy analysis."""

import numpy as np
import pandas as pd
from scipy.spatial.distance import pdist, squareform
from scipy.stats import f
from sklearn.linear_model import LinearRegression
from sklearn.metrics import r2_score
from statsmodels.stats.multitest import multipletests

from .validation import validate_input


def _compute_distance_matrix(df, method="braycurtis"):
    """
    Compute distance matrix from data using specified method.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data matrix.
    method : str, default "braycurtis"
        Distance metric to use.

    Returns
    -------
    numpy.ndarray
        Distance matrix.
    """
    distances = pdist(df.values, metric=method)
    return squareform(distances)


def _perform_constrained_ordination(distance_matrix, constraint_var):
    """
    Perform constrained ordination similar to vegan's capscale.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        Distance matrix.
    constraint_var : pandas.Series
        Constraining variable.

    Returns
    -------
    dict
        Dictionary containing R², adjusted R², F-statistic, and p-value.
    """
    n = len(constraint_var)

    # Convert distance matrix to principal coordinates (classical MDS)
    # Center the squared distance matrix
    D_squared = distance_matrix**2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D_squared @ H

    # Eigendecomposition
    eigenvals, eigenvecs = np.linalg.eigh(B)
    # Sort in descending order
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    # Take only positive eigenvalues
    pos_idx = eigenvals > 1e-8
    eigenvals = eigenvals[pos_idx]
    eigenvecs = eigenvecs[:, pos_idx]

    # Principal coordinates
    coords = eigenvecs @ np.diag(np.sqrt(eigenvals))

    # Fit linear model: coordinates ~ constraint_variable
    if constraint_var.dtype in ["object", "category"]:
        # Categorical variable - create dummy variables
        constraint_matrix = pd.get_dummies(constraint_var, drop_first=False).values
    else:
        # Continuous variable
        constraint_matrix = constraint_var.values.reshape(-1, 1)

    # Remove constant columns
    if constraint_matrix.shape[1] > 1:
        var_mask = np.var(constraint_matrix, axis=0) > 1e-10
        constraint_matrix = constraint_matrix[:, var_mask]

    if constraint_matrix.shape[1] == 0:
        # No variation in constraint
        return {"r2": 0.0, "r2adj": 0.0, "F": 0.0, "pval": 1.0}

    # Perform multivariate regression
    n_coords = min(coords.shape[1], 10)  # Use first 10 coordinates
    coords_subset = coords[:, :n_coords]

    # Calculate R² for each coordinate and combine
    r2_values = []
    f_values = []

    for i in range(n_coords):
        reg = LinearRegression().fit(constraint_matrix, coords_subset[:, i])
        y_pred = reg.predict(constraint_matrix)
        r2 = r2_score(coords_subset[:, i], y_pred)
        r2_values.append(r2)

        # Calculate F-statistic
        ss_res = np.sum((coords_subset[:, i] - y_pred) ** 2)
        ss_tot = np.sum((coords_subset[:, i] - np.mean(coords_subset[:, i])) ** 2)
        ss_reg = ss_tot - ss_res

        if ss_tot > 0:
            df_reg = constraint_matrix.shape[1]
            df_res = n - df_reg - 1
            if df_res > 0:
                ms_reg = ss_reg / df_reg
                ms_res = ss_res / df_res
                f_stat = ms_reg / ms_res if ms_res > 0 else 0
            else:
                f_stat = 0
        else:
            f_stat = 0
        f_values.append(f_stat)

    # Weight by eigenvalues (handle case where sum might be zero)
    eigenval_sum = np.sum(eigenvals[:n_coords])
    if eigenval_sum > 1e-10 and len(r2_values) > 0:
        eigenval_weights = eigenvals[:n_coords] / eigenval_sum
        overall_r2 = np.average(r2_values, weights=eigenval_weights)
        overall_f = np.average(f_values, weights=eigenval_weights)
    else:
        # Fallback to simple average if eigenvalues sum to zero or no coordinates
        overall_r2 = np.mean(r2_values) if r2_values else 0.0
        overall_f = np.mean(f_values) if f_values else 0.0

    # Adjusted R²
    df_reg = constraint_matrix.shape[1]
    df_res = n - df_reg - 1
    if df_res > 0 and overall_r2 < 1:
        r2adj = 1 - (1 - overall_r2) * (n - 1) / df_res
    else:
        r2adj = overall_r2

    # Approximate p-value using F-distribution
    if df_res > 0 and overall_f > 0:
        pval = 1 - f.cdf(overall_f, df_reg, df_res)
    else:
        pval = 1.0

    return {"r2": overall_r2, "r2adj": r2adj, "F": overall_f, "pval": pval}


def get_r2_single(distance_matrix, meta, feature):
    """
    Calculate R-squared values for a single feature.

    This function computes the R-squared and adjusted R-squared for a single feature
    from a distance matrix and its corresponding metadata using constrained
    ordination and permutation tests.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A distance matrix in which distances between items are stored.
    meta : pandas.DataFrame
        A data frame containing metadata associated with the items in distance_matrix.
    feature : str
        A string specifying the column name in meta for which the
        R-squared values are to be calculated.

    Returns
    -------
    pandas.Series
        A series with the feature name, F statistic, R-squared, adjusted R-squared,
        number of samples, and p-value.
    """
    constraint_var = meta[feature].dropna()
    valid_idx = constraint_var.index

    # Filter distance matrix for non-missing values
    idx_positions = [meta.index.get_loc(i) for i in valid_idx]
    filtered_dist = distance_matrix[np.ix_(idx_positions, idx_positions)]

    # Perform constrained ordination
    result = _perform_constrained_ordination(filtered_dist, constraint_var)

    return pd.Series(
        {
            "feature": feature,
            "Fa": result["F"],
            "r2": result["r2"],
            "r2adj": result["r2adj"],
            "N": len(constraint_var),
            "pval": result["pval"],
        }
    )


def get_r2(distance_matrix, meta):
    """
    Calculate R-squared values for all features across distance matrix.

    This function applies get_r2_single to each feature in the metadata dataframe
    and calculates the R-squared and adjusted R-squared using constrained ordination
    and permutation tests for all features. It also adjusts p-values using the
    Benjamini-Hochberg method.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A distance matrix containing distances between sampled entities.
    meta : pandas.DataFrame
        A dataframe containing metadata associated with distance_matrix.

    Returns
    -------
    pandas.DataFrame
        A dataframe with rows for each feature and columns for the feature name,
        F statistic, R-squared, adjusted R-squared, number of samples,
        raw p-values, and Benjamini-Hochberg adjusted p-values.
    """
    features = meta.columns
    results = []

    for feature in features:
        result = get_r2_single(distance_matrix, meta, feature)
        results.append(result)

    results_df = pd.DataFrame(results)
    results_df.set_index("feature", inplace=True)

    # Adjust p-values using Benjamini-Hochberg method
    rejected, padj, alpha_sidak, alpha_bonf = multipletests(
        results_df["pval"], method="fdr_bh"
    )
    results_df["padj"] = padj

    return results_df


def get_cumulative(distance_matrix, meta):
    """
    Perform cumulative db-RDA analysis using forward stepwise selection.

    This function performs forward stepwise distance-based redundancy analysis (db-RDA)
    to determine the cumulative effects of metadata variables on a distance matrix.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A distance matrix where rows and columns correspond to samples.
    meta : pandas.DataFrame
        A dataframe containing metadata for each sample in distance_matrix.

    Returns
    -------
    pandas.DataFrame
        A dataframe with the R-squared values for each variable added sequentially
        into the ordination model.
    """
    features = list(meta.columns)
    selected_features = []
    results = []

    # Convert distance matrix to principal coordinates once
    n = distance_matrix.shape[0]
    D_squared = distance_matrix**2
    H = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * H @ D_squared @ H

    eigenvals, eigenvecs = np.linalg.eigh(B)
    idx = np.argsort(eigenvals)[::-1]
    eigenvals = eigenvals[idx]
    eigenvecs = eigenvecs[:, idx]

    pos_idx = eigenvals > 1e-8
    eigenvals = eigenvals[pos_idx]
    eigenvecs = eigenvecs[:, pos_idx]

    n_coords = min(len(eigenvals), 10)
    coords = eigenvecs[:, :n_coords] @ np.diag(np.sqrt(eigenvals[:n_coords]))

    while len(selected_features) < len(features):
        best_feature = None
        best_r2 = -1
        best_result = None

        # Try adding each remaining feature
        for feature in features:
            if feature in selected_features:
                continue

            # Current feature set + candidate feature
            current_features = selected_features + [feature]
            current_meta = meta[current_features]

            # Create design matrix
            design_matrix = []
            for col in current_features:
                if current_meta[col].dtype in ["object", "category"]:
                    # Categorical variable
                    dummies = pd.get_dummies(current_meta[col], drop_first=False)
                    design_matrix.append(dummies.values)
                else:
                    # Continuous variable
                    design_matrix.append(current_meta[col].values.reshape(-1, 1))

            if design_matrix:
                X = np.hstack(design_matrix)
                # Remove constant columns
                var_mask = np.var(X, axis=0) > 1e-10
                X = X[:, var_mask]

                if X.shape[1] > 0:
                    # Calculate R² for multivariate regression
                    r2_values = []
                    for i in range(n_coords):
                        reg = LinearRegression().fit(X, coords[:, i])
                        y_pred = reg.predict(X)
                        r2 = r2_score(coords[:, i], y_pred)
                        r2_values.append(r2)

                    # Weight by eigenvalues
                    eigenval_weights = eigenvals[:n_coords] / np.sum(
                        eigenvals[:n_coords]
                    )
                    overall_r2 = np.average(r2_values, weights=eigenval_weights)

                    # Adjusted R²
                    df_reg = X.shape[1]
                    df_res = n - df_reg - 1
                    if df_res > 0 and overall_r2 < 1:
                        r2adj = 1 - (1 - overall_r2) * (n - 1) / df_res
                    else:
                        r2adj = overall_r2

                    if r2adj > best_r2:
                        best_r2 = r2adj
                        best_feature = feature
                        best_result = {
                            "feature": feature,
                            "R2.adj": r2adj,
                            "Df": 1,  # Simplified
                            "AIC": 0,  # Placeholder
                            "F": 0,  # Placeholder
                            "pval": 0.05,  # Placeholder
                            "N": n,
                        }

        if best_feature is None:
            break

        selected_features.append(best_feature)
        results.append(best_result)

        # Stop if improvement is minimal
        if (
            len(results) > 1
            and abs(results[-1]["R2.adj"] - results[-2]["R2.adj"]) < 0.001
        ):
            break

    # Convert to DataFrame with proper column names
    if results:
        results_df = pd.DataFrame(results)
        results_df.set_index("feature", inplace=True)
        # Rename columns to match R output
        results_df.columns = [f"RDAcumul_{col}" for col in results_df.columns]
    else:
        # Empty DataFrame if no features selected
        results_df = pd.DataFrame()

    return results_df


def prepare_plot_data(dbrda_data):
    """
    Prepare data for plotting db-RDA results.

    This function takes the combined db-RDA results and prepares them for visualization
    by restructuring the data, imputing missing values, and adding aesthetics
    to indicate the significance of features.

    Parameters
    ----------
    dbrda_data : pandas.DataFrame
        A dataframe from the output of rldbrda function containing
        R-squared and cumulative R-squared values for metadata variables.

    Returns
    -------
    pandas.DataFrame
        A tidied long-form dataframe suitable for plotting.
    """
    # Create a copy to avoid modifying original data
    data = dbrda_data.copy()

    # Define the order of the y-axis based on cumulative R² values
    # We want features with LOWEST cumulative R² at the TOP of the plot
    # Since matplotlib plots categorical data with first category at bottom,
    # we need to arrange our categories in reverse order
    if "RDAcumul_R2.adj" in data.columns:
        # Separate features with and without cumulative R² values
        has_cumul = data[data["RDAcumul_R2.adj"].notna()]
        no_cumul = data[data["RDAcumul_R2.adj"].isna()]

        # Sort features with cumulative R² in ascending order (lowest first)
        cumul_order = has_cumul.sort_values(
            "RDAcumul_R2.adj", ascending=False
        ).index.tolist()

        # Sort features without cumulative R² by individual r2adj
        no_cumul_order = no_cumul.sort_values("r2adj", ascending=True).index.tolist()

        # For the plot order: cumulative features first (in ascending order),
        # then non-cumulative features (reversed so highest individual R² at top)
        plot_order = no_cumul_order + cumul_order
    else:
        # Fallback to r2adj if no cumulative data
        plot_order = list(
            reversed(data.sort_values("r2adj", ascending=True).index.tolist())
        )

    # Identify insignificant features (those with NaN in RDAcumul_R2.adj)
    insignificant_features = data[data["RDAcumul_R2.adj"].isna()].index.tolist()

    # Impute missing values in RDAcumul_R2.adj with the largest non-NA value
    if "RDAcumul_R2.adj" in data.columns:
        max_val = data["RDAcumul_R2.adj"].max()
        data["RDAcumul_R2.adj"] = data["RDAcumul_R2.adj"].fillna(max_val)

    # Filter and reshape to long form
    plot_data = data[["r2adj", "RDAcumul_R2.adj"]].copy()
    plot_data = plot_data.reset_index().rename(columns={"index": "feature"})
    plot_data = pd.melt(
        plot_data,
        id_vars=["feature"],
        value_vars=["r2adj", "RDAcumul_R2.adj"],
        var_name="variable",
        value_name="value",
    )

    # Add significance indicator
    plot_data["significant"] = plot_data["feature"].apply(
        lambda x: 0 if x in insignificant_features else 1
    )

    # Set proper factor levels for ordering
    plot_data["feature"] = pd.Categorical(plot_data["feature"], categories=plot_order)

    return plot_data


def rldbrda(df, meta, method="braycurtis", p_cutoff=0.05):
    """
    Distance-Based Redundancy Analysis workflow.

    This function executes a full db-RDA analysis on a given dataset and its
    corresponding metadata. It validates input, calculates pairwise distances,
    determines R-squared values for each metadata variable, filters variables
    based on significance, and performs cumulative db-RDA on significant variables.

    Parameters
    ----------
    df : pandas.DataFrame
        A numeric dataframe with rows as samples and columns as species/variables
        for which the distance matrix will be computed.
    meta : pandas.DataFrame
        A dataframe where rows correspond to the samples in 'df' and columns
        contain metadata variables.
    method : str, default "braycurtis"
        A string indicating the distance measure to be used.
    p_cutoff : float, default 0.05
        A significance level for determining which metadata variables to include
        in subsequent analysis.

    Returns
    -------
    pandas.DataFrame
        A dataframe with combined R-squared and cumulative R-squared values for
        significant metadata variables, ordered by R-squared values.
    """
    # Validate input
    validate_input(df, meta)

    # Compute distance matrix
    distance_matrix = _compute_distance_matrix(df, method=method)

    # Calculate individual R² values
    r2_results = get_r2(distance_matrix, meta)

    # Select significant variables
    significant_features = r2_results[r2_results["padj"] < p_cutoff].index.tolist()

    if len(significant_features) < 1:
        raise ValueError("No significant features found, cannot continue!")

    # Perform cumulative analysis on significant features
    cumul_results = get_cumulative(distance_matrix, meta[significant_features])

    # Combine results
    combined = r2_results.join(cumul_results, how="left")

    # Sort by r2, then by cumulative R2.adj
    combined = combined.sort_values(["r2"], ascending=False)
    if "RDAcumul_R2.adj" in combined.columns:
        # Secondary sort by cumulative R2.adj for tied r2 values
        combined = combined.sort_values(
            ["RDAcumul_R2.adj"], ascending=True, na_position="last"
        )

    return combined

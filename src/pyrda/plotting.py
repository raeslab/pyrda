"""Visualization functions for db-RDA results."""

import matplotlib.pyplot as plt
import numpy as np


def plot_dbrda(plot_data, figsize=(10, 8), colors=None, alpha_range=(0.5, 1.0)):
    """
    Visualize db-RDA effect sizes.

    This function creates a horizontal bar plot to visualize the effect sizes from
    db-RDA analysis. Features (variables) are displayed on the y-axis and their
    effect sizes on the x-axis. The plot differentiates between individual adjusted
    R-squared values and cumulative R-squared values, and indicates the significance
    of features.

    Parameters
    ----------
    plot_data : pandas.DataFrame
        A long-form dataframe prepared by prepare_plot_data containing
        the variables with their effect sizes and an indicator for significance.
        Expected columns are 'feature' for features, 'variable' for effect type,
        'value' for effect size, and 'significant' for significance indication.
    figsize : tuple, default (10, 8)
        Figure size as (width, height).
    colors : dict, optional
        Colors for different variable types. Default uses green and blue.
    alpha_range : tuple, default (0.5, 1.0)
        Range of alpha values for non-significant and significant features.

    Returns
    -------
    matplotlib.figure.Figure
        The matplotlib figure object.
    """
    if colors is None:
        colors = {
            "r2adj": "#60A68B",  # Green for individual R²
            "RDAcumul_R2.adj": "#1F4068",  # Blue for cumulative R²
        }

    # Create figure and axis
    fig, ax = plt.subplots(figsize=figsize)

    # Get unique features and their order
    features = plot_data["feature"].cat.categories
    n_features = len(features)

    # Set up y positions
    y_pos = np.arange(n_features)
    bar_height = 0.35

    # Plot bars for each variable type
    for i, var_type in enumerate(["r2adj", "RDAcumul_R2.adj"]):
        data_subset = plot_data[plot_data["variable"] == var_type]

        # Get values and significance for each feature
        values = []
        alphas = []
        for feature in features:
            feature_data = data_subset[data_subset["feature"] == feature]
            if not feature_data.empty:
                values.append(feature_data["value"].iloc[0])
                significant = feature_data["significant"].iloc[0]
                alphas.append(alpha_range[1] if significant else alpha_range[0])
            else:
                values.append(0)
                alphas.append(alpha_range[0])

        # Create bars with individual alpha values
        bars = ax.barh(
            y_pos + i * bar_height - bar_height / 2,
            values,
            bar_height,
            color=colors[var_type],
            alpha=1.0,  # We'll set alpha individually
            label="R²" if var_type == "r2adj" else "Cumulative R²",
        )

        # Set individual alpha for each bar
        for bar, alpha in zip(bars, alphas):
            bar.set_alpha(alpha)

    # Customize the plot
    ax.set_yticks(y_pos)
    ax.set_yticklabels(features)
    ax.set_xlabel("Effect size")
    ax.set_ylabel("Feature")
    ax.legend(loc="lower right")

    # Apply minimal theme
    ax.spines["top"].set_visible(False)
    ax.spines["right"].set_visible(False)
    ax.grid(True, alpha=0.3, axis="x")
    ax.set_axisbelow(True)

    # Adjust layout
    plt.tight_layout()

    return fig


def save_plot(fig, filename, dpi=300, bbox_inches="tight"):
    """
    Save plot to file.

    Parameters
    ----------
    fig : matplotlib.figure.Figure
        The figure to save.
    filename : str
        Output filename.
    dpi : int, default 300
        Resolution in dots per inch.
    bbox_inches : str, default "tight"
        Bounding box setting.
    """
    fig.savefig(filename, dpi=dpi, bbox_inches=bbox_inches)


def create_summary_plot(dbrda_results, save_path=None, show_plot=True):
    """
    Create a complete summary plot from db-RDA results.

    Parameters
    ----------
    dbrda_results : pandas.DataFrame
        Results from rldbrda function.
    save_path : str, optional
        Path to save the plot. If None, plot is not saved.
    show_plot : bool, default True
        Whether to display the plot.

    Returns
    -------
    matplotlib.figure.Figure
        The matplotlib figure object.
    """
    from .core import prepare_plot_data

    # Prepare plot data
    plot_data = prepare_plot_data(dbrda_results)

    # Create plot
    fig = plot_dbrda(plot_data)

    # Save if requested
    if save_path is not None:
        save_plot(fig, save_path)

    # Show if requested
    if show_plot:
        plt.show()

    return fig

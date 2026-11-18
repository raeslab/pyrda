#!/usr/bin/env python3
"""
Example usage of the pyrda package for distance-based redundancy analysis.

This script demonstrates how to use pyrda to perform distance-based redundancy
analysis on microbial community data, similar to the R implementation in RLdbRDA.
"""

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from pyrda import plot_dbrda, prepare_plot_data, rldbrda


def load_example_data():
    """
    Load the example varespec and varechem datasets.

    Returns
    -------
    tuple
        (varespec, varechem) DataFrames
    """
    try:
        # Try to load the actual data files
        varespec = pd.read_csv("./data/varespec.tsv", sep="\t", index_col=0)
        varechem = pd.read_csv("./data/varechem.tsv", sep="\t", index_col=0)
        print("Loaded example datasets from ./data/")
    except FileNotFoundError:
        # Create simulated data if files not found
        print("Example data files not found, creating simulated data...")
        np.random.seed(42)

        # Create simulated species abundance data
        n_samples = 24
        n_species = 44
        sample_names = [str(i) for i in range(18, 42)]
        species_names = [f"Species_{i:02d}" for i in range(1, n_species + 1)]

        varespec = pd.DataFrame(
            np.random.lognormal(mean=1, sigma=1, size=(n_samples, n_species)),
            index=sample_names,
            columns=species_names,
        )

        # Create simulated environmental data
        env_vars = [
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
        ]
        varechem = pd.DataFrame(
            np.random.normal(loc=50, scale=20, size=(n_samples, len(env_vars))),
            index=sample_names,
            columns=env_vars,
        )
        # Make sure pH is realistic
        varechem["pH"] = np.random.normal(loc=6, scale=1, size=n_samples)
        varechem = varechem.abs()  # Ensure all values are positive

    return varespec, varechem


def main():
    """Main example workflow."""
    print("PyRDA Example: Distance-Based Redundancy Analysis")
    print("=" * 50)

    # Load example data
    print("\n1. Loading example data...")
    varespec, varechem = load_example_data()

    print(f"   Species data shape: {varespec.shape}")
    print(f"   Environmental data shape: {varechem.shape}")
    print("   Species data preview:")
    print(varespec.iloc[:5, :5])
    print("\n   Environmental data preview:")
    print(varechem.head())

    # Perform db-RDA analysis
    print("\n2. Performing distance-based redundancy analysis...")
    try:
        results = rldbrda(
            varespec,
            varechem,
            method="bray",  # Bray-Curtis dissimilarity (vegan::vegdist names)
            p_cutoff=0.05,  # Significance threshold
            # test="permutation", seed=1 would reproduce set.seed(1) +
            # RLdbRDA::rldbrda() in R exactly; the default analytic test needs
            # no permutations
        )
        print("   Analysis completed successfully!")
    except ValueError as e:
        print(f"   Error in analysis: {e}")
        return

    # Display results
    print("\n3. Analysis results:")
    print(f"   Total variables analyzed: {len(results)}")

    significant_vars = results[results["padj"] < 0.05]
    print(f"   Significant variables (p < 0.05): {len(significant_vars)}")

    if len(significant_vars) > 0:
        print("\n   Top significant variables:")
        print(significant_vars[["r2adj", "padj", "RDAcumul_R2.adj"]].head())

    print("\n   All results (top 10):")
    print(results[["r2adj", "padj"]].head(10))

    # Prepare data for plotting
    print("\n4. Preparing visualization...")
    plot_data = prepare_plot_data(results)
    print(f"   Plot data shape: {plot_data.shape}")

    # Create plot
    print("\n5. Creating effect size plot...")
    plot_dbrda(plot_data, figsize=(10, 8))
    plt.title("Distance-Based Redundancy Analysis Results")

    # Save plot
    output_file = "pyrda_example_plot.png"
    plt.savefig(output_file, dpi=300, bbox_inches="tight")
    print(f"   Plot saved as: {output_file}")

    # Show plot (comment out if running in headless environment)
    try:
        plt.show()
    except Exception:
        print("   (Plot display not available in this environment)")

    # Save results to file
    results_file = "pyrda_example_results.csv"
    results.to_csv(results_file)
    print(f"\n6. Results saved to: {results_file}")

    # Summary statistics
    print("\nSummary Statistics:")
    print(f"   Mean R² (adjusted): {results['r2adj'].mean():.4f}")
    print(f"   Max R² (adjusted): {results['r2adj'].max():.4f}")
    print(f"   Variables with R² > 0.1: {sum(results['r2adj'] > 0.1)}")

    if "RDAcumul_R2.adj" in results.columns:
        cumul_data = results["RDAcumul_R2.adj"].dropna()
        if len(cumul_data) > 0:
            print(f"   Cumulative R² (final): {cumul_data.iloc[-1]:.4f}")

    print("\nAnalysis complete!")


if __name__ == "__main__":
    main()

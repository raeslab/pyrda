"""Worked example of the README: db-RDA of the varespec/varechem data.

Reads the example data in ``docs/data/`` and rewrites ``docs/output/`` (the
results table and the figure shown in the README). Run it from the repository
root with ``uv run python docs/example.py``.

Data: vegan's ``varespec`` (cover of 44 plant species at 24 sites) and
``varechem`` (14 soil variables at the same sites), from Väre, Ohtonen &
Oksanen (1995), J. Veg. Sci. 6: 523-530.
"""

from pathlib import Path

import pandas as pd

from pyrda import plot_dbrda, prepare_plot_data, rldbrda

DOCS = Path(__file__).parent

species = pd.read_csv(DOCS / "data" / "varespec.tsv", sep="\t", index_col=0)
soil = pd.read_csv(DOCS / "data" / "varechem.tsv", sep="\t", index_col=0)

# db-RDA on Bray-Curtis distances, analytic p-values (deterministic)
results = rldbrda(species, soil, method="bray", p_cutoff=0.05)
print(results[["r2", "r2adj", "pval", "padj", "RDAcumul_R2.adj"]].round(3))

out = DOCS / "output"
out.mkdir(exist_ok=True)
results.to_csv(out / "dbrda_results.tsv", sep="\t")

fig = plot_dbrda(prepare_plot_data(results), figsize=(8, 6))
fig.savefig(out / "dbrda_plot.png", dpi=150, bbox_inches="tight")
print(f"\nWrote {out / 'dbrda_results.tsv'} and {out / 'dbrda_plot.png'}")

"""Create the synthetic inputs for the R reference tests (run from this folder)."""

import numpy as np
import pandas as pd


def community(n, n_taxa, drivers, rng, noise=1.0, depth=5000):
    """Negative-binomial counts whose log-abundance depends on `drivers`."""
    base = rng.normal(0, 1.5, n_taxa)
    load = rng.normal(0, 1, (drivers.shape[1], n_taxa))
    eta = base + drivers @ load + rng.normal(0, noise, (len(drivers), n_taxa))
    p = np.exp(eta)
    p /= p.sum(1, keepdims=True)
    mu = p * depth
    counts = rng.negative_binomial(0.8, 0.8 / (0.8 + mu))
    index = [f"s{i + 1}" for i in range(len(drivers))]
    return pd.DataFrame(
        counts, index=index, columns=[f"t{j + 1}" for j in range(n_taxa)]
    )


rng = np.random.default_rng(2024)

# mixed: continuous, factor and logical metadata
n = 30
g1, g2 = rng.normal(size=n), rng.normal(size=n)
site = rng.choice(["A", "B", "C"], n)
mixed = pd.DataFrame(
    {
        "ph": g1 + rng.normal(0, 0.5, n),
        "temp": g2 + rng.normal(0, 1, n),
        "moist": rng.normal(size=n),
        "site": site,
        "treated": rng.random(n) > 0.5,
        "noise": rng.normal(size=n),
    }
)
spec = community(n, 60, np.column_stack([g1, 0.6 * g2, (site == "A") * 1.0]), rng)
mixed.index = spec.index
spec.to_csv("mixed_data.tsv", sep="\t")
mixed.to_csv("mixed_meta.tsv", sep="\t")

# collinear: duplicated and linearly dependent variables
n = 20
x, z = rng.normal(size=n), rng.normal(size=n)
coll = pd.DataFrame(
    {"x": x, "x_copy": x, "z": z, "xz": x + 2 * z, "w": rng.normal(size=n)}
)
spec = community(n, 40, np.column_stack([x, z]), rng, noise=0.6)
coll.index = spec.index
spec.to_csv("collinear_data.tsv", sep="\t")
coll.to_csv("collinear_meta.tsv", sep="\t")

# small data set for the distance indices
small = pd.DataFrame(
    rng.negative_binomial(0.7, 0.2, size=(8, 12)), columns=[f"t{j}" for j in range(12)]
)
small.index = [f"s{i}" for i in range(8)]
small.to_csv("small_counts.tsv", sep="\t")

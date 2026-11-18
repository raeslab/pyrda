"""Dissimilarity indices matching ``vegan::vegdist()``.

The formulas follow vegan's C implementation (``src/vegdist.c``) and the
transformations done in R before the C call (``decostand()``), so that
``vegdist(x, "bray")`` here gives the same matrix as ``vegdist(x, "bray")`` in R.
"""

import warnings

import numpy as np
from scipy.spatial.distance import pdist, squareform
from scipy.stats import hypergeom

VEGAN_METHODS = (
    "manhattan",
    "euclidean",
    "canberra",
    "bray",
    "kulczynski",
    "gower",
    "morisita",
    "horn",
    "mountford",
    "jaccard",
    "raup",
    "binomial",
    "chao",
    "altGower",
    "cao",
    "mahalanobis",
    "clark",
    "chisq",
    "chord",
    "hellinger",
    "aitchison",
    "robust.aitchison",
)

# scipy names that refer to the same index as a vegan method
_ALIASES = {"braycurtis": "bray", "cityblock": "manhattan", "euclidian": "euclidean"}

# vegan methods that do not warn about negative entries / empty rows
_NEG_OK = {"manhattan", "euclidean", "canberra", "gower", "mahalanobis", "chord"}
_NEG_OK |= {"hellinger"}
_EMPTY_OK = {"manhattan", "euclidean", "gower", "mahalanobis", "chisq"}

_ZAP = 1e-15
_LN2 = np.log(2.0)


def resolve_method(method):
    """Map a method name to a vegan method (R's ``pmatch``) or a scipy metric.

    Returns ``(name, is_vegan)``.
    """
    name = _ALIASES.get(method, method)
    if name in VEGAN_METHODS:
        return name, True
    prefixed = [m for m in VEGAN_METHODS if m.startswith(name)]
    if len(prefixed) == 1:
        return prefixed[0], True
    if len(prefixed) > 1:
        raise ValueError(f"ambiguous distance method {method!r}: {prefixed}")
    return method, False


def vegdist(x, method="bray", binary=False):
    """Compute a dissimilarity matrix like ``vegan::vegdist()``.

    Parameters
    ----------
    x : array-like or pandas.DataFrame
        Samples in rows, variables (species) in columns.
    method : str, default "bray"
        Any ``vegdist`` method name (partial matching as in R), or the scipy
        aliases ``"braycurtis"``/``"cityblock"``. Names that are not vegan
        methods are passed to :func:`scipy.spatial.distance.pdist` (these have
        no R equivalent). ``"robust.aitchison"`` is not supported because vegan
        imputes zeros with OptSpace matrix completion.
    binary : bool, default False
        Convert data to presence/absence first (``decostand(x, "pa")``).

    Returns
    -------
    numpy.ndarray
        Square, symmetric dissimilarity matrix.
    """
    x = np.asarray(x, dtype=float)
    if x.ndim != 2:
        raise ValueError("input data must be a 2-D matrix")
    if np.isnan(x).any():
        raise ValueError("missing values are not allowed in the data")

    name, is_vegan = resolve_method(method)
    if binary:
        x = (x > 0).astype(float)
    if not is_vegan:
        return squareform(pdist(x, metric=name))
    if name == "robust.aitchison":
        raise NotImplementedError(
            "'robust.aitchison' is not supported: vegan imputes zeros with "
            "OptSpace. Compute the distances in R and pass the matrix to "
            "get_r2()/get_cumulative() instead."
        )

    if name not in _EMPTY_OK and (x.sum(axis=1) == 0).any():
        warnings.warn(
            f"you have empty rows: their dissimilarities may be meaningless in "
            f"method {name!r}",
            stacklevel=2,
        )
    if name not in _NEG_OK | {"aitchison"} and (x < 0).any():
        warnings.warn(
            f"results may be meaningless because data have negative entries in "
            f"method {name!r}",
            stacklevel=2,
        )

    d = _DISPATCH[name](x)
    d = np.asarray(d, dtype=float)
    with np.errstate(invalid="ignore"):
        d[d < _ZAP] = 0.0
    if np.isnan(d).any():
        warnings.warn("missing values in results", stacklevel=2)
    np.fill_diagonal(d, 0.0)
    return d


# --------------------------------------------------------------------------
# helpers


def _pairwise(x, func):
    """Fill a distance matrix row by row, ``func(xi, xj_block) -> distances``.

    Mirrors vegan's loop where the first argument is the later row.
    """
    n, p = x.shape
    block = max(1, 2_000_000 // max(1, p))  # bound temporaries to ~16 MB each
    d = np.zeros((n, n))
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        for i in range(1, n):
            for start in range(0, i, block):
                stop = min(i, start + block)
                row = func(x[i], x[start:stop])
                d[i, start:stop] = row
                d[start:stop, i] = row
    return d


def _square(cond):
    return squareform(cond, checks=False)


# --------------------------------------------------------------------------
# individual indices (see vegan/src/vegdist.c)


def _manhattan(x):
    return _square(pdist(x, "cityblock"))


def _euclidean(x):
    return _square(pdist(x, "euclidean"))


def _bray(x):
    d = _square(pdist(x, "cityblock"))
    s = x.sum(axis=1)
    with np.errstate(divide="ignore", invalid="ignore"):
        for i in range(len(d)):  # row by row: no extra n x n temporaries
            d[i] /= s[i] + s
    return d


def _jaccard(x):
    b = _bray(x)
    return 2 * b / (1 + b)


def _kulczynski(x):
    def f(xi, xj):
        sim = np.minimum(xi, xj).sum(axis=1)
        dist = 1 - sim / xi.sum() / 2 - sim / xj.sum(axis=1) / 2
        return np.where(dist < 0, 0.0, dist)

    return _pairwise(x, f)


def _canberra(x):
    def f(xi, xj):
        nz = (xi != 0) | (xj != 0)
        denom = np.abs(xi) + np.abs(xj)
        term = np.where(denom > 0, np.abs(xi - xj) / denom, np.inf)
        count = nz.sum(axis=1)
        dist = np.where(nz, term, 0.0).sum(axis=1) / count
        return np.where(count == 0, np.nan, dist)

    return _pairwise(x, f)


def _clark(x):
    def f(xi, xj):
        nz = (xi != 0) | (xj != 0)
        denom = xi + xj
        term = np.where(denom > 0, ((xi - xj) / denom) ** 2, np.inf)
        count = nz.sum(axis=1)
        dist = np.sqrt(np.where(nz, term, 0.0).sum(axis=1) / count)
        return np.where(count == 0, np.nan, dist)

    return _pairwise(x, f)


def _gower(x):
    # decostand(x, "range", MARGIN = 2); constant columns become NaN and are
    # skipped by the C code (they do not count towards the number of columns)
    with np.errstate(divide="ignore", invalid="ignore"):
        lo = x.min(axis=0)
        x = (x - lo) / (x.max(axis=0) - lo)
    x = x[:, ~np.isnan(x).all(axis=0)]
    if x.shape[1] == 0:
        return np.full((len(x), len(x)), np.nan)
    return _manhattan(x) / x.shape[1]


def _alt_gower(x):
    def f(xi, xj):
        valid = (xi > 0) | (xj > 0)
        count = valid.sum(axis=1)
        dist = np.where(valid, np.abs(xi - xj), 0.0).sum(axis=1) / count
        return np.where(count == 0, np.nan, dist)

    return _pairwise(x, f)


def _morisita(x):
    def f(xi, xj):
        sim = (xi * xj).sum(axis=1)
        t1, t2 = xi.sum(), xj.sum(axis=1)
        tlam1, tlam2 = (xi * (xi - 1)).sum(), (xj * (xj - 1)).sum(axis=1)
        dist = 1 - 2 * sim / (tlam1 / t1 / (t1 - 1) + tlam2 / t2 / (t2 - 1)) / t1 / t2
        dist = np.where(sim <= 0, 1.0, dist)
        return np.where(dist < 0, 0.0, dist)

    return _pairwise(x, f)


def _horn(x):
    def f(xi, xj):
        sim = (xi * xj).sum(axis=1)
        t1, t2 = xi.sum(), xj.sum(axis=1)
        sq1, sq2 = (xi * xi).sum(), (xj * xj).sum(axis=1)
        dist = 1 - 2 * sim / (sq1 / t1 / t1 + sq2 / t2 / t2) / t1 / t2
        return np.where(dist < 0, 0.0, dist)

    return _pairwise(x, f)


def _counts(xi, xj):
    sim = ((xi > 0) & (xj > 0)).sum(axis=1).astype(float)
    t1 = np.full(len(xj), float((xi > 0).sum()))
    t2 = (xj > 0).sum(axis=1).astype(float)
    return sim, t1, t2


def _mountford(x):
    def f(xi, xj):
        sim, t1, t2 = _counts(xi, xj)
        dist = np.full(len(xj), np.nan)
        dist[sim == 0] = 0.0
        full = (sim > 0) & ((sim == t1) | (sim == t2))
        dist[full] = _LN2
        newton = (sim > 0) & ~full
        dist[(t1 == 0) | (t2 == 0)] = np.nan
        newton &= (t1 > 0) & (t2 > 0)
        if newton.any():
            j, a, b = sim[newton], t1[newton], t2[newton]
            th = 2 * j / (2 * a * b - (a + b) * j)
            active = np.ones(len(th), dtype=bool)
            for _ in range(20):
                old = th.copy()
                fun = np.exp(th * a) + np.exp(th * b) - np.exp(th * (a + b - j)) - 1
                der = a * np.exp(th * a) + b * np.exp(th * b)
                der -= (a + b - j) * np.exp(th * (a + b - j))
                th = np.where(active, th - fun / der, th)
                done = (np.abs(old - th) / old < 1e-5) | (np.abs(old - th) < 1e-12)
                active &= ~done
                if not active.any():
                    break
            dist[newton] = th
        return _LN2 - dist

    return _pairwise(x, f)


def _raup(x):
    count = x.shape[1]

    def f(xi, xj):
        sim, t1, t2 = _counts(xi, xj)
        a = np.minimum(t1, t2)
        b = np.maximum(t1, t2)
        return 1 - hypergeom.cdf(sim - 1, count, a, b)

    return _pairwise(x, f)


def _binomial(x):
    def f(xi, xj):
        nk = xi + xj
        used = nk != 0
        lognk = np.log(np.where(used, nk, 1.0))
        t1 = np.where(xi > 0, xi * (np.log(np.where(xi > 0, xi, 1.0)) - lognk), 0.0)
        xj_safe = np.where(xj > 0, xj, 1.0)
        t2 = np.where(xj > 0, xj * (np.log(xj_safe) - lognk), 0.0)
        terms = np.where(used, (t1 + t2 + nk * _LN2) / np.where(used, nk, 1.0), 0.0)
        dist = terms.sum(axis=1)
        dist = np.where(used.sum(axis=1) == 0, np.nan, dist)
        return np.where(dist < 0, 0.0, dist)

    return _pairwise(x, f)


def _chao(x):
    def f(xi, xj):
        # xi is vegan's "i1" (later row), xj its "i2"
        shared = (xi > 0) & (xj > 0)
        itot, jtot = xi.sum(), xj.sum(axis=1)
        ishare = np.where(shared, xi, 0.0).sum(axis=1)
        jshare = np.where(shared, xj, 0.0).sum(axis=1)
        j1 = shared & (np.abs(xj - 1) < 0.01)
        j2 = shared & ~j1 & (np.abs(xj - 2) < 0.01)
        i1 = shared & (np.abs(xi - 1) < 0.01)
        i2 = shared & ~i1 & (np.abs(xi - 2) < 0.01)
        ishar1 = np.where(j1, xi, 0.0).sum(axis=1)
        jonce, jtwice = j1.sum(axis=1).astype(float), j2.sum(axis=1).astype(float)
        jshar1 = np.where(i1, xj, 0.0).sum(axis=1)
        ionce, itwice = i1.sum(axis=1).astype(float), i2.sum(axis=1).astype(float)

        u = ishare / itot
        jonce, jtwice = np.maximum(jonce, 1), np.maximum(jtwice, 1)
        u_add = (jtot - 1) / jtot * jonce / jtwice / 2.0 * ishar1 / itot
        u = np.where(ishar1 > 0, u + u_add, u)
        u = np.where(u > 1, 1.0, u)
        v = jshare / jtot
        ionce, itwice = np.maximum(ionce, 1), np.maximum(itwice, 1)
        v_add = (itot - 1) / itot * ionce / itwice / 2.0 * jshar1 / jtot
        v = np.where(jshar1 > 0, v + v_add, v)
        v = np.where(v > 1, 1.0, v)
        dist = np.where((u <= 0) | (v <= 0), 1.0, 1 - u * v / (u + v - u * v))
        return np.where(dist < 0, 0.0, dist)

    return _pairwise(x, f)


def _cao(x):
    def f(xi, xj):
        used = ~((xi == 0) & (xj == 0))
        x1 = np.where(xi < 0.1, 0.1, xi)
        x2 = np.where(xj < 0.1, 0.1, xj)
        t1 = x1 + x2
        t2 = x1 * np.log(x2) + x2 * np.log(x1)
        dist = np.where(used, np.log(t1) - _LN2 - t2 / t1, 0.0).sum(axis=1)
        count = used.sum(axis=1)
        dist = np.where(dist < 0, 0.0, dist) / count
        return np.where(count == 0, np.nan, dist)

    return _pairwise(x, f)


def _chord(x):
    def f(xi, xj):
        cp = (xi * xj).sum(axis=1)
        ss1, ss2 = (xi * xi).sum(), (xj * xj).sum(axis=1)
        dist = np.sqrt(np.maximum(2.0 * (1.0 - cp / np.sqrt(ss1 * ss2)), 0.0))
        return np.where(ss1 * ss2 <= 0, np.nan, dist)

    return _pairwise(x, f)


def _hellinger(x):
    def f(xi, xj):
        cp = np.sqrt(xi * xj).sum(axis=1)
        ss1, ss2 = xi.sum(), xj.sum(axis=1)
        dist = np.sqrt(np.maximum(2.0 * (1.0 - cp / np.sqrt(ss1 * ss2)), 0.0))
        return np.where(ss1 * ss2 <= 0, np.nan, dist)

    return _pairwise(x, f)


def _euclidean_nan(x):
    """Euclidean distance skipping NaN columns pairwise (veg_euclidean)."""
    if not np.isnan(x).any():
        return _euclidean(x)

    def f(xi, xj):
        valid = ~np.isnan(xi) & ~np.isnan(xj)
        dist = np.sqrt(np.where(valid, (xi - xj) ** 2, 0.0).sum(axis=1))
        return np.where(valid.sum(axis=1) == 0, np.nan, dist)

    return _pairwise(x, f)


def _chisq(x):
    # decostand(x, "chi.square"): empty species give NaN columns, skipped later
    if (x.sum(axis=0) == 0).any():
        warnings.warn(
            "data have empty species which influence the results in method 'chisq'",
            stacklevel=3,
        )
    rs, cs, tot = x.sum(axis=1), x.sum(axis=0), x.sum()
    with np.errstate(divide="ignore", invalid="ignore"):
        x = np.sqrt(tot) * x / np.outer(rs, np.sqrt(cs))
    return _euclidean_nan(x)


def _mahalanobis(x):
    x = x - x.mean(axis=0)
    s2 = np.cov(x, rowvar=False, ddof=1).reshape(x.shape[1], x.shape[1])
    ev, vec = np.linalg.eigh(s2)
    ev, vec = ev[::-1], vec[:, ::-1]
    tol = np.sqrt(np.finfo(float).eps)
    keep = ev > max(tol, tol * ev[0])
    sisqr = vec[:, keep] @ (np.sqrt(1 / ev[keep])[:, None] * vec[:, keep].T)
    return _euclidean(x @ sisqr)


def _aitchison(x):
    if (x <= 0).any():
        raise ValueError(
            "method 'aitchison' (clr) cannot be used with non-positive data: "
            "add a pseudocount first"
        )
    clog = np.log(x)
    return _euclidean(clog - clog.mean(axis=1, keepdims=True))


_DISPATCH = {
    "manhattan": _manhattan,
    "euclidean": _euclidean,
    "canberra": _canberra,
    "bray": _bray,
    "kulczynski": _kulczynski,
    "gower": _gower,
    "morisita": _morisita,
    "horn": _horn,
    "mountford": _mountford,
    "jaccard": _jaccard,
    "raup": _raup,
    "binomial": _binomial,
    "chao": _chao,
    "altGower": _alt_gower,
    "cao": _cao,
    "mahalanobis": _mahalanobis,
    "clark": _clark,
    "chisq": _chisq,
    "chord": _chord,
    "hellinger": _hellinger,
    "aitchison": _aitchison,
}

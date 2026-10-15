"""Core analysis functions for distance-based redundancy analysis.

The numerics reproduce what RLdbRDA does in R:

* ``vegan::capscale()`` - principal coordinates of the distance matrix (keeping
  all axes with positive eigenvalues), followed by a redundancy analysis of
  those coordinates. Negative eigenvalues ("imaginary" inertia) are excluded
  from R-squared and from the residual of the F statistic.
* ``vegan::anova.cca()`` - pseudo-F with a permutation test (999 permutations
  for the individual tests, 9999 in the forward selection).
* ``vegan::ordiR2step()`` - forward selection on adjusted R-squared, bounded by
  the adjusted R-squared of the model with all candidate variables.

By default (``test="analytic"``) the permutation tests are not run: their
p-values are computed from the exact moments of the permutation distribution
(see :mod:`pyrda.moments`). With ``test="permutation"`` the permutations are
drawn with an exact copy of R's random number generator, so ``seed=s`` gives
the same p-values as ``set.seed(s)`` before ``RLdbRDA::rldbrda()`` in R.
"""

import warnings
from dataclasses import dataclass

import numpy as np
import pandas as pd
from scipy.linalg import eigh

from .distance import vegdist
from .moments import GowerMoments
from .rrandom import RRandom
from .validation import validate_input

_ZERO = np.sqrt(np.finfo(float).eps)  # vegan's tolerance for zero eigenvalues
_EPS = _ZERO  # vegan's tolerance when counting permuted F >= observed F
_QR_TOL = 1e-7  # tolerance of R's qr() for aliased (collinear) columns
# Adjusted R-squared values closer than this are equal. ordiR2step compares
# R-squared of models spanning the same space (e.g. the last candidate and the
# full scope); in R the outcome then depends on rounding, here on exact equality.
_TIE_TOL = 1e-10
# test="analytic" warns when the permutation distribution's skewness exceeds
# this: the Pearson type III curve may no longer follow it. Calibrated against
# 99,999-permutation p-values; it flags the outlier-dominated cases where the
# approximation fails and not the naturally skewed small-n ones.
_MAX_SKEW = 2.0
_TESTS = ("analytic", "permutation")
_ALL_VARIABLES = "<All variables>"
_CUMUL_COLUMNS = ["R2.adj", "Df", "AIC", "F", "pval", "N"]


def _compute_distance_matrix(df, method="bray"):
    """
    Compute distance matrix from data using specified method.

    Parameters
    ----------
    df : pandas.DataFrame
        Input data matrix.
    method : str, default "bray"
        Distance metric, see :func:`pyrda.distance.vegdist`.

    Returns
    -------
    numpy.ndarray
        Distance matrix.
    """
    return vegdist(df, method=method)


@dataclass
class _Ordination:
    """Principal coordinates as used by capscale."""

    Y: np.ndarray  # real axes, column-centred (n x k)
    total: float  # real inertia, sum(Y**2)
    imaginary: float  # sum of negative eigenvalues
    adjust: float  # capscale divides distances by sqrt(n - 1) if max(d) >= 4

    @property
    def n(self):
        return self.Y.shape[0]

    def moments(self):
        """Permutation moments of the Gower matrix (for ``test="analytic"``)."""
        if not hasattr(self, "_moments"):
            self._moments = GowerMoments(self.Y)
        return self._moments


@dataclass
class _Fit:
    """A capscale model: constrained (CCA) and residual (CA) inertia."""

    Q: np.ndarray  # orthonormal basis of the centred constraints (n x rank)
    rank: int  # rank of the constraints (vegan's qrank)
    cca: float  # constrained inertia
    ca: float  # residual (real) inertia

    def df_residual(self, n):
        return n - self.rank - 1


def _pcoa(distance_matrix):
    """Principal coordinates exactly as ``capscale`` -> ``wcmdscale`` does."""
    d = np.asarray(distance_matrix, dtype=float)
    n = d.shape[0]
    adjust = 1.0
    if n > 1 and d.max() >= 4 + np.finfo(float).eps:
        adjust = np.sqrt(n - 1)
    # Gower centring in place: one n x n copy plus the eigenvectors
    m = np.square(d)
    if adjust != 1.0:
        m /= adjust**2
    m -= m.mean(axis=0)
    m -= m.mean(axis=1, keepdims=True)
    m *= -0.5
    ev, vec = eigh(m, overwrite_a=True, check_finite=False)
    del m
    ev, vec = ev[::-1], vec[:, ::-1]
    keep = np.abs(ev) > max(_ZERO, _ZERO * ev[0])
    pos = keep & (ev > _ZERO)
    neg = keep & (ev < 0)
    if pos.any():
        y = vec[:, pos] * np.sqrt(ev[pos])
    else:
        y = np.zeros((n, 1))
    y = y - y.mean(axis=0)
    return _Ordination(
        Y=y, total=float((y**2).sum()), imaginary=float(ev[neg].sum()), adjust=adjust
    )


def _is_categorical(series):
    return not (
        pd.api.types.is_numeric_dtype(series) or pd.api.types.is_bool_dtype(series)
    )


def _design_columns(series):
    """Model-matrix columns for one variable (no intercept), like R's model.matrix.

    Numeric and logical variables enter as one column; factors (strings,
    categories) as treatment-coded dummies of the observed levels.
    """
    if not _is_categorical(series):
        return np.asarray(series, dtype=float).reshape(-1, 1)
    values = (
        series.astype(str)
        if not isinstance(series.dtype, pd.CategoricalDtype)
        else series
    )
    levels = pd.unique(values.dropna())
    levels = sorted(levels, key=str)
    if len(levels) <= 1:
        # vegan adds a dummy level, giving an all-zero (aliased) column
        return np.zeros((len(series), 1))
    return np.column_stack([(values == lev).to_numpy(float) for lev in levels[1:]])


def _design_matrix(meta, columns):
    return np.hstack([_design_columns(meta[c]) for c in columns])


def _orthonormal_basis(x):
    """Basis of the centred columns of ``x``, dropping aliased ones like R's qr().

    A column is aliased when its norm after removing the previous columns falls
    below ``1e-7`` times its original norm (LINPACK ``dqrdc2``).
    """
    x = x - x.mean(axis=0)
    basis = []
    for col in x.T:
        norm0 = np.linalg.norm(col)
        ref = norm0 if norm0 > 0 else 1.0
        resid = col.copy()
        if basis:
            q = np.column_stack(basis)
            resid -= q @ (q.T @ resid)
            resid -= q @ (q.T @ resid)  # re-orthogonalise for stability
        norm = np.linalg.norm(resid)
        if norm >= _QR_TOL * ref and norm > 0:
            basis.append(resid / norm)
    if not basis:
        return np.zeros((x.shape[0], 0))
    return np.column_stack(basis)


def _fit(ordination, x):
    """Fit constraints ``x`` to the ordination (vegan's ordConstrain)."""
    q = _orthonormal_basis(x)
    rank = q.shape[1]
    cca = 0.0
    if rank:
        sv = np.linalg.svd(q.T @ ordination.Y, compute_uv=False)
        lam = sv**2
        lam = lam[np.abs(lam) >= max(_ZERO, _ZERO * lam[0])]
        cca = float(lam.sum())
    return _Fit(Q=q, rank=rank, cca=cca, ca=ordination.total - cca)


def _r2(ordination, fit):
    """R-squared and adjusted R-squared (vegan's RsquareAdj for capscale)."""
    n = ordination.n
    r2 = fit.cca / ordination.total
    m = fit.rank
    if m >= n - 1:
        return r2, np.nan
    return r2, 1 - (1 - r2) * (n - 1) / (n - m - 1)


def _aic(ordination, fit):
    """extractAIC() of a capscale model (as reported by add1/ordiR2step)."""
    n = ordination.n
    rss = fit.ca * (n - 1)
    edf = fit.rank + 1
    return n * np.log(rss / n) + 2 * edf


def _permuted_fitted_ss(y, q, perms, max_elements=20_000_000):
    """Constrained inertia ``||Q' Y[perm]||^2`` for every permutation."""
    inv = np.argsort(perms, axis=1)  # Q' Y[perm] == (Q[inv])' Y
    n, k = y.shape
    rank = q.shape[1]
    chunk = max(1, max_elements // max(1, (n + k) * rank))
    out = np.empty(len(perms))
    for start in range(0, len(perms), chunk):
        idx = inv[start : start + chunk]
        # one large matrix product instead of many small ones
        z = q[idx].transpose(1, 0, 2).reshape(n, -1)  # (n, c * rank)
        ss = ((y.T @ z) ** 2).sum(axis=0).reshape(len(idx), rank)
        out[start : start + chunk] = ss.sum(axis=1)
    return out


class _Tester:
    """How p-values are computed: the test, permutations and random stream."""

    def __init__(self, test, permutations, seed):
        if test not in _TESTS:
            raise ValueError(f"test must be one of {_TESTS}, got {test!r}")
        self.test = test
        self.permutations = permutations
        self.rng = seed if isinstance(seed, RRandom) else RRandom(seed)
        self.skewed = []  # labels of analytic tests that may be inaccurate

    def draw(self, n):
        return self.rng.shuffle_set(n, self.permutations)

    def analytic_pvalue(self, ordination, factor, weights, threshold, label):
        """Expected permutation p-value without permuting.

        Pearson type III approximation of the share of permutations reaching
        the observed statistic, reported as ``(N p + 1) / (N + 1)``: the value
        R's test with ``N`` permutations gives on average. Tests whose
        permutation distribution is too skewed for the approximation (a few
        samples dominate) are remembered for :meth:`warn_skewed`.
        """
        p, skew = ordination.moments().tail_and_skewness(factor, weights, threshold)
        if abs(skew) > _MAX_SKEW:
            self.skewed.append(label)
        return (self.permutations * p + 1) / (self.permutations + 1)

    def warn_skewed(self, where):
        if self.skewed:
            names = ", ".join(repr(str(s)) for s in self.skewed)
            its, them = ("their", "them") if len(self.skewed) > 1 else ("its", "it")
            warnings.warn(
                f"{where}: the analytic p-value may be inaccurate for {names}: "
                f"{its} permutation distribution is strongly skewed (a few samples "
                f"dominate the variation). Use test='permutation' to check {them}.",
                stacklevel=3,
            )


def _test_model(ordination, fit, tester, label=None):
    """F statistic and p-value of a single model (vegan's anova.cca)."""
    n = ordination.n
    perms = None
    if tester.test == "permutation":
        # vegan draws the permutations before checking the model
        perms = tester.draw(n)
    if fit.rank == 0 or fit.cca == 0:
        return 0.0, np.nan
    r = fit.df_residual(n)
    if r <= 0 or fit.ca <= 0:
        return np.nan, np.nan
    f0 = (fit.cca / fit.rank) / (fit.ca / r)
    if tester.test == "analytic":
        # F(perm) >= F0  <=>  tr(H P G P') >= observed constrained inertia
        weights = np.ones(fit.rank)
        return f0, tester.analytic_pvalue(ordination, fit.Q, weights, fit.cca, label)
    num = _permuted_fitted_ss(ordination.Y, fit.Q, perms)
    den = fit.cca + fit.ca - num
    f_perm = (num / fit.rank) / (den / r)
    pval = (np.sum(f_perm >= f0 - _EPS) + 1) / (len(perms) + 1)
    return f0, float(pval)


def _test_nested(ordination, small, big, tester, label=None):
    """F test of ``big`` against nested ``small`` (vegan's anovaCCAlist)."""
    n = ordination.n
    perms = tester.draw(n) if tester.test == "permutation" else None
    df = small.df_residual(n) - big.df_residual(n)
    resdf = big.df_residual(n)
    if df < 1 or resdf <= 0:
        return np.nan, np.nan
    fval = (small.ca - big.ca) / df / (big.ca / resdf)
    if tester.test == "analytic":
        # F(perm) >= fval  <=>  tr(((1 + c) H_big - H_small) P G P') >= c * total
        c = fval * df / resdf
        factor = np.hstack([big.Q, small.Q])
        weights = np.r_[np.full(big.rank, 1 + c), np.full(small.rank, -1.0)]
        total = big.cca + big.ca
        p = tester.analytic_pvalue(ordination, factor, weights, c * total, label)
        return fval, p
    num_small = _permuted_fitted_ss(ordination.Y, small.Q, perms)
    num_big = _permuted_fitted_ss(ordination.Y, big.Q, perms)
    pscale = (big.cca + big.ca - num_big) / resdf
    f_perm = (num_big - num_small) / df / pscale
    pval = (np.sum(f_perm >= fval - _EPS) + 1) / (len(perms) + 1)
    return fval, float(pval)


def _r2_single(ordination_cache, distance_matrix, meta, feature, tester):
    values = meta[feature]
    keep = values.notna().to_numpy()
    if keep.all():
        ordination = ordination_cache()
    else:
        idx = np.flatnonzero(keep)
        ordination = _pcoa(np.asarray(distance_matrix)[np.ix_(idx, idx)])
    fit = _fit(ordination, _design_columns(values[keep]))
    if fit.rank == 0:
        warnings.warn(
            f"variable {feature!r} has no variation: it cannot explain anything",
            stacklevel=3,
        )
    r2, r2adj = _r2(ordination, fit)
    f0, pval = _test_model(ordination, fit, tester, label=feature)
    return pd.Series(
        {
            "feature": feature,
            "Fa": f0,
            "r2": r2,
            "r2adj": r2adj,
            "N": int(keep.sum()),
            "pval": pval,
        }
    )


def get_r2_single(
    distance_matrix, meta, feature, test="analytic", permutations=999, seed=None
):
    """
    Calculate R-squared values for a single feature.

    Equivalent to ``capscale(distmat ~ meta[, feature], na.action = na.omit)``
    followed by ``anova.cca()`` and ``RsquareAdj()`` in R. Samples with a
    missing value for ``feature`` are dropped.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A square distance matrix.
    meta : pandas.DataFrame
        Metadata, rows matching the rows of ``distance_matrix``.
    feature : str
        The column of ``meta`` to test.
    test : {"analytic", "permutation"}, default "analytic"
        "analytic" gives the p-value of vegan's permutation test without
        permuting, from the exact mean, variance and skewness of the
        permutation distribution (Pearson type III). It warns, naming the
        variables, when that distribution is too skewed for the approximation;
        check those with "permutation". "permutation" runs vegan's permutation
        test itself, with R's random number generator.
    permutations : int, default 999
        Number of permutations (as ``anova.cca`` in R). "analytic" reports the
        p-value on the same scale: ``(N p + 1) / (N + 1)``, the value a test
        with ``N`` permutations gives on average.
    seed : int or RRandom, optional
        Seed for R's random number generator (``set.seed()`` equivalent); only
        used by ``test="permutation"``.

    Returns
    -------
    pandas.Series
        The feature name, F statistic (``Fa``), R-squared, adjusted R-squared,
        number of samples and p-value.
    """
    tester = _Tester(test, permutations, seed)
    full = _cached_pcoa(distance_matrix)
    result = _r2_single(full, distance_matrix, meta, feature, tester)
    tester.warn_skewed("get_r2_single")
    return result


def _cached_pcoa(distance_matrix):
    cache = {}

    def get():
        if "pcoa" not in cache:
            cache["pcoa"] = _pcoa(distance_matrix)
        return cache["pcoa"]

    return get


def _p_adjust_bh(p):
    """R's ``p.adjust(p, "BH")``; missing p-values are left out of the count."""
    p = np.asarray(p, dtype=float)
    out = np.full_like(p, np.nan)
    ok = ~np.isnan(p)
    pv = p[ok]
    lp = len(pv)
    if lp == 0:
        return out
    if lp == 1:
        out[ok] = pv
        return out
    o = np.argsort(-pv, kind="stable")
    i = np.arange(lp, 0, -1)
    adj = np.minimum(1, np.minimum.accumulate(lp / i * pv[o]))
    res = np.empty(lp)
    res[o] = adj
    out[ok] = res
    return out


def get_r2(distance_matrix, meta, test="analytic", permutations=999, seed=None):
    """
    Calculate R-squared values for all features across distance matrix.

    Applies :func:`get_r2_single` to each column of ``meta`` (in order, sharing
    one random stream as R does) and adds Benjamini-Hochberg adjusted
    p-values.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A square distance matrix.
    meta : pandas.DataFrame
        Metadata, rows matching the rows of ``distance_matrix``.
    test : {"analytic", "permutation"}, default "analytic"
        See :func:`get_r2_single`.
    permutations : int, default 999
        Number of permutations per feature.
    seed : int or RRandom, optional
        Seed for R's random number generator (``set.seed()`` equivalent).

    Returns
    -------
    pandas.DataFrame
        One row per feature with columns feature, Fa, r2, r2adj, N, pval and
        padj (Benjamini-Hochberg).
    """
    tester = _Tester(test, permutations, seed)
    full = _cached_pcoa(distance_matrix)
    rows = [
        _r2_single(full, distance_matrix, meta, feature, tester)
        for feature in meta.columns
    ]
    tester.warn_skewed("get_r2")
    results = pd.DataFrame(rows)
    results.index = results["feature"].tolist()
    results["N"] = results["N"].astype(int)
    results["padj"] = _p_adjust_bh(results["pval"])
    return results


def get_cumulative(
    distance_matrix, meta, test="analytic", permutations=9999, seed=None
):
    """
    Perform cumulative db-RDA analysis using forward stepwise selection.

    Mirrors RLdbRDA's ``get_cumul()``: ``ordiR2step(capscale(d ~ 1),
    scope = all variables, Pin = 1, R2scope = TRUE, permutations = 9999)``.
    At each step the variable giving the largest adjusted R-squared is added,
    as long as it increases the adjusted R-squared and does not exceed that of
    the model with all variables. Because ``Pin = 1`` the p-values are
    reported but never used to stop the selection.

    Parameters
    ----------
    distance_matrix : numpy.ndarray
        A square distance matrix.
    meta : pandas.DataFrame
        Candidate variables (no missing values).
    test : {"analytic", "permutation"}, default "analytic"
        Test used for the reported F and p-value of each step (see
        :func:`get_r2_single`); it never changes which variables are selected.
    permutations : int, default 9999
        Number of permutations per step.
    seed : int or RRandom, optional
        Seed for R's random number generator (``set.seed()`` equivalent).

    Returns
    -------
    pandas.DataFrame
        One row per selected variable (in order of selection) plus a final
        ``"<All variables>"`` row, with columns ``RDAcumul_R2.adj``,
        ``RDAcumul_Df``, ``RDAcumul_AIC``, ``RDAcumul_F``, ``RDAcumul_pval``
        and ``RDAcumul_N``.
    """
    tester = _Tester(test, permutations, seed)
    if meta.isnull().any().any():
        raise ValueError("Missing values found in metadata.")
    ordination = _pcoa(distance_matrix)
    n = ordination.n
    scope = list(meta.columns)
    columns = {c: _design_columns(meta[c]) for c in scope}

    def fit_of(terms):
        if not terms:
            return _Fit(Q=np.zeros((n, 0)), rank=0, cca=0.0, ca=ordination.total)
        return _fit(ordination, np.hstack([columns[t] for t in terms]))

    r2_all = _r2(ordination, fit_of(scope))[1]
    if np.isnan(r2_all):
        raise ValueError("R2 cannot be adjusted in upper scope (too many terms?)")

    selected = []
    current = fit_of(selected)
    r2_previous = 0.0
    rows = []
    while len(selected) < len(scope):
        candidates = [t for t in scope if t not in selected]
        fits = [fit_of(selected + [t]) for t in candidates]
        r2_adds = np.array([_r2(ordination, ft)[1] for ft in fits])
        r2_adds = np.where(np.isnan(r2_adds), 0.0, r2_adds)
        # which.max(): first of the (numerically) tied maxima
        best = int(np.flatnonzero(r2_adds >= r2_adds.max() - _TIE_TOL)[0])
        improves = r2_adds[best] > r2_previous + _TIE_TOL
        if not (improves and r2_adds[best] <= r2_all + _TIE_TOL):
            break
        new = fits[best]
        if current.rank == 0:
            f_val, p_val = _test_model(ordination, new, tester, candidates[best])
        else:
            f_val, p_val = _test_nested(
                ordination, current, new, tester, candidates[best]
            )
        df_added = new.rank - current.rank
        if df_added < 1:
            f_val, p_val = np.nan, np.nan
        selected.append(candidates[best])
        current = new
        r2_previous = min(r2_adds[best], r2_all)  # a tie with the scope is equal
        rows.append(
            [r2_previous, df_added, _aic(ordination, new), f_val, p_val, len(meta)]
        )

    tester.warn_skewed("get_cumulative")
    if not rows:
        warnings.warn(
            "no variable improved the adjusted R-squared: nothing was selected "
            "(RLdbRDA stops with an error in this case)",
            stacklevel=2,
        )
    rows.append([r2_all, np.nan, np.nan, np.nan, np.nan, len(meta)])
    result = pd.DataFrame(
        rows,
        index=selected + [_ALL_VARIABLES],
        columns=[f"RDAcumul_{c}" for c in _CUMUL_COLUMNS],
    )
    result["RDAcumul_N"] = result["RDAcumul_N"].astype(int)
    return result


def _combine(r2, cumul):
    """R's combine_data(): merge, then order by r2 and cumulative R2.adj."""
    combined = r2.join(cumul, how="outer")
    combined = combined.loc[sorted(combined.index, key=str)]
    combined = combined.sort_values(
        "r2", ascending=False, kind="stable", na_position="last"
    )
    combined = combined.sort_values(
        "RDAcumul_R2.adj", ascending=True, kind="stable", na_position="last"
    )
    return combined


def prepare_plot_data(dbrda_data):
    """
    Prepare data for plotting db-RDA results.

    Same as RLdbRDA's ``prepare_plot_data()``: the plot order is the reversed
    row order of the results, features that were not selected in the
    cumulative model are flagged as not significant, and their cumulative
    R-squared is set to the largest one (that of ``"<All variables>"``).

    Parameters
    ----------
    dbrda_data : pandas.DataFrame
        Output of :func:`rldbrda`.

    Returns
    -------
    pandas.DataFrame
        Long-form data with columns ``feature``, ``variable``, ``value`` and
        ``significant``; ``feature`` is categorical in plot order.
    """
    data = dbrda_data.copy()
    plot_order = list(reversed(data.index.tolist()))
    insignificant = set(data.index[data["RDAcumul_R2.adj"].isna()])
    max_val = data["RDAcumul_R2.adj"].max()
    data["RDAcumul_R2.adj"] = data["RDAcumul_R2.adj"].fillna(max_val)
    data = data[data.index != _ALL_VARIABLES]

    rows = []
    for name, rec in data.iterrows():
        sig = 0 if name in insignificant else 1
        for variable in ("r2adj", "RDAcumul_R2.adj"):
            rows.append((name, variable, float(rec[variable]), sig))
    plot_data = pd.DataFrame(
        rows, columns=["feature", "variable", "value", "significant"]
    )
    plot_data["feature"] = pd.Categorical(plot_data["feature"], categories=plot_order)
    return plot_data


def rldbrda(
    df,
    meta,
    method="bray",
    p_cutoff=0.05,
    test="analytic",
    permutations=999,
    cumul_permutations=9999,
    seed=None,
):
    """
    Distance-Based Redundancy Analysis workflow.

    Python equivalent of ``RLdbRDA::rldbrda()``: validates the input,
    computes the distance matrix, tests every metadata variable on its own,
    and runs forward selection on the variables with a BH-adjusted p-value
    below ``p_cutoff``.

    Parameters
    ----------
    df : pandas.DataFrame
        A numeric dataframe with rows as samples and columns as species/variables
        for which the distance matrix will be computed.
    meta : pandas.DataFrame
        A dataframe where rows correspond to the samples in 'df' and columns
        contain metadata variables.
    method : str, default "bray"
        Distance measure, any ``vegan::vegdist`` method (see
        :func:`pyrda.distance.vegdist`).
    p_cutoff : float, default 0.05
        Variables with ``padj < p_cutoff`` enter the forward selection.
    test : {"analytic", "permutation"}, default "analytic"
        "analytic" gives R's permutation p-values within their own Monte Carlo
        error without permuting: much faster, and deterministic.
        "permutation" runs the permutation tests and, with ``seed``,
        reproduces R exactly. See :func:`get_r2_single`.
    permutations : int, default 999
        Permutations for the per-variable tests.
    cumul_permutations : int, default 9999
        Permutations for the tests reported in the forward selection.
    seed : int, optional
        Seed for R's random number generator, used by ``test="permutation"``:
        ``seed=s`` reproduces ``set.seed(s); rldbrda(df, meta)`` in R exactly.

    Returns
    -------
    pandas.DataFrame
        The same table as RLdbRDA: one row per feature plus an
        ``"<All variables>"`` row, ordered by cumulative adjusted R-squared
        (selected features) and then by R-squared.
    """
    rng = _Tester(test, permutations, seed).rng  # also validates `test`
    validate_input(df, meta)

    distance_matrix = _compute_distance_matrix(df, method=method)
    r2_results = get_r2(distance_matrix, meta, test, permutations, rng)

    significant_features = r2_results.index[r2_results["padj"] < p_cutoff].tolist()
    if len(significant_features) < 1:
        raise ValueError("No significant features found, cannot continue!")

    cumul_results = get_cumulative(
        distance_matrix, meta[significant_features], test, cumul_permutations, rng
    )
    return _combine(r2_results, cumul_results)

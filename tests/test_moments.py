"""Exact permutation moments, checked against full enumeration."""

import itertools

import numpy as np
import pytest

from pyrda.moments import GowerMoments


def centred(rng, n, k):
    x = rng.normal(size=(n, k))
    return x - x.mean(axis=0)


@pytest.mark.parametrize("n", [5, 7])
def test_moments_match_enumeration(n):
    rng = np.random.default_rng(n)
    y = centred(rng, n, 4)
    q, _ = np.linalg.qr(centred(rng, n, 2))
    factor = np.hstack([q, q[:, :1]])
    weights = np.array([1.7, 1.7, -1.0])  # the nested-test matrix
    a = factor @ np.diag(weights) @ factor.T
    g = y @ y.T
    stats = np.array(
        [np.trace(a @ g[np.ix_(p, p)]) for p in itertools.permutations(range(n))]
    )
    expected = [np.mean(stats**k) for k in (1, 2, 3)]
    got = GowerMoments(y).raw_moments(factor, weights)
    np.testing.assert_allclose(got, expected, rtol=1e-10)


def test_upper_tail_close_to_enumeration():
    rng = np.random.default_rng(1)
    n = 8
    y = centred(rng, n, 5)
    q, _ = np.linalg.qr(centred(rng, n, 1))
    g = y @ y.T
    h = q @ q.T
    stats = np.array(
        [np.trace(h @ g[np.ix_(p, p)]) for p in itertools.permutations(range(n))]
    )
    moments = GowerMoments(y)
    for quantile in (0.5, 0.9, 0.95):
        t = np.quantile(stats, quantile)
        exact = np.mean(stats >= t)
        approx = moments.upper_tail(q, np.ones(1), t)
        assert approx == pytest.approx(exact, abs=0.02)

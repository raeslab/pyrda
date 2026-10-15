"""Permutation-free approximation of vegan's permutation F tests.

Every F test in the workflow (a single model, or a model against a nested one)
rejects when a statistic ``W(P) = tr(A P G P')`` exceeds a threshold, where
``G = Y Y'`` is the (double-centred) Gower matrix of the principal coordinates,
``P`` a permutation of the samples and ``A`` a double-centred symmetric matrix
built from the projections onto the constraints. The first three moments of
``W`` over *all* ``n!`` permutations can be computed exactly (Kazi-Aoual et al.
1995, Comput. Stat. Data Anal. 20:643-656); the p-value is then taken from a
Pearson type III distribution with those moments.

The moments are obtained generically: ``E[W^k]`` is a sum over set partitions
of the ``2k`` row/column indices, and each term is a sum over distinct indices
recovered from ordinary (unrestricted) graph contractions by Moebius inversion
on the partition lattice. Contractions containing a vertex of degree one vanish
because both matrices are double-centred.
"""

import itertools
import math
from functools import cache

import numpy as np
from scipy import stats


def _set_partitions(m):
    """All set partitions of ``range(m)`` as restricted growth strings."""

    def rec(prefix, top):
        if len(prefix) == m:
            yield tuple(prefix)
            return
        for v in range(top + 2):
            yield from rec(prefix + [v], max(top, v))

    yield from rec([0], 0)


def _canonical_rgs(labels):
    mapping = {}
    return tuple(mapping.setdefault(lab, len(mapping)) for lab in labels)


def _graph_key(rgs, k):
    """Canonical form of the multigraph with edges (rgs[2e], rgs[2e+1])."""
    n_vertices = max(rgs) + 1
    edges = [(rgs[2 * e], rgs[2 * e + 1]) for e in range(k)]
    best = None
    for perm in itertools.permutations(range(n_vertices)):
        relabelled = tuple(sorted(tuple(sorted((perm[u], perm[v]))) for u, v in edges))
        if best is None or relabelled < best:
            best = relabelled
    return n_vertices, best


@cache
def _structure(k):
    """Data-independent part of E[W^k].

    Returns the surviving graphs, a matrix ``C`` mapping graph contractions to
    distinct-index sums for every set partition, and the partition sizes.
    """
    partitions = list(_set_partitions(2 * k))
    graphs = {}
    coefficients = []
    for part in partitions:
        n_blocks = max(part) + 1
        coef = {}
        for merge in _set_partitions(n_blocks):
            coarser = _canonical_rgs([merge[b] for b in part])
            mobius = 1
            for size in np.bincount(merge):
                mobius *= (-1) ** (size - 1) * math.factorial(size - 1)
            key = _graph_key(coarser, k)
            degree = np.zeros(key[0], dtype=int)
            for u, v in key[1]:
                degree[u] += 1
                degree[v] += 1
            if (degree == 1).any():
                continue  # zero for double-centred matrices
            idx = graphs.setdefault(key, len(graphs))
            coef[idx] = coef.get(idx, 0) + mobius
        coefficients.append(coef)
    graph_list = [None] * len(graphs)
    for key, idx in graphs.items():
        graph_list[idx] = key
    c = np.zeros((len(partitions), len(graph_list)))
    for i, coef in enumerate(coefficients):
        for j, value in coef.items():
            c[i, j] = value
    sizes = np.array([max(p) + 1 for p in partitions])
    return graph_list, c, sizes


def _component_type(vertices, edges):
    """Name of a connected multigraph with at most three edges."""
    loops = sorted(u for u, v in edges if u == v)
    links = [(u, v) for u, v in edges if u != v]
    if len(vertices) == 1:
        return f"loops{len(loops)}"
    if len(vertices) == 2:
        if not loops:
            return f"links{len(links)}"  # double or triple edge
        if len(loops) == 1:
            return "loop+double"
        return "loop+loop+link"
    return "triangle"


def _components(n_vertices, edges):
    """Split a multigraph into connected components (vertices, edges)."""
    parent = list(range(n_vertices))

    def find(v):
        while parent[v] != v:
            v = parent[v]
        return v

    for u, v in edges:
        parent[find(u)] = find(v)
    groups = {}
    for v in range(n_vertices):
        groups.setdefault(find(v), ([], []))[0].append(v)
    for u, v in edges:
        groups[find(u)][1].append((u, v))
    return list(groups.values())


class _Contractions:
    """Sums over the connected graphs that occur, for ``M = F diag(w) F'``.

    Every surviving multigraph with up to three edges is a product of these
    components, each computed in closed form. A low-rank ``M`` (the
    constraints) is handled in the rank space, a full-rank one (the Gower
    matrix) through the ``n x n`` matrix; nothing costs more than ``O(n^2 r)``.
    """

    def __init__(self, factor, weights):
        n, r = factor.shape
        fw = factor * weights
        wg = weights[:, None] * (factor.T @ factor)  # W F'F, r x r
        d = (factor * fw).sum(axis=1)  # diagonal of M
        self.values = {
            "loops1": d.sum(),
            "loops2": (d**2).sum(),
            "loops3": (d**3).sum(),
            "triangle": np.trace(wg @ wg @ wg),
            "links2": np.trace(wg @ wg),
        }
        fd = factor.T @ d
        self.values["loop+loop+link"] = fd @ (weights * fd)
        # (M^2)_aa = f_a W F'F W f_a'
        m2_diag = ((fw @ (factor.T @ fw)) * factor).sum(axis=1)
        self.values["loop+double"] = d @ m2_diag
        if r * r <= n:
            # sum_ab M_ab^3 = sum_{rst} w_r w_s w_t (sum_a F_ar F_as F_at)^2
            t = np.einsum("ar,as,at->rst", factor, factor, factor)
            ww = np.einsum("r,s,t->rst", weights, weights, weights)
            self.values["links3"] = (ww * t**2).sum()
        else:
            # build M in row blocks to keep memory at O(n) rows at a time
            step = max(1, 4_000_000 // n)
            self.values["links3"] = sum(
                float(((fw[i : i + step] @ factor.T) ** 3).sum())
                for i in range(0, n, step)
            )

    def graph_sum(self, graph):
        n_vertices, edges = graph
        total = 1.0
        for vertices, comp_edges in _components(n_vertices, edges):
            total *= self.values[_component_type(vertices, comp_edges)]
        return total


def _graph_sums(factor, weights, graph_list):
    """Contractions of ``M = factor @ diag(weights) @ factor.T`` over each graph."""
    contractions = _Contractions(np.asarray(factor, float), np.asarray(weights, float))
    return np.array([contractions.graph_sum(g) for g in graph_list])


class GowerMoments:
    """Precomputed contractions of the Gower matrix ``G = Y Y'``."""

    def __init__(self, y):
        self.n = y.shape[0]
        w = np.ones(y.shape[1])
        self._sums = {k: _graph_sums(y, w, _structure(k)[0]) for k in (1, 2, 3)}

    def raw_moments(self, factor, weights):
        """``E[W^k]`` for ``k = 1, 2, 3`` with ``A = factor diag(weights) factor'``."""
        n = self.n
        moments = []
        for k in (1, 2, 3):
            graph_list, c, sizes = _structure(k)
            s_a = c @ _graph_sums(factor, weights, graph_list)
            s_g = c @ self._sums[k]
            valid = sizes <= n
            falling = np.array([math.perm(n, b) for b in sizes[valid]], dtype=float)
            moments.append(float(np.sum(s_a[valid] * s_g[valid] / falling)))
        return moments

    def upper_tail(self, factor, weights, threshold):
        """Pearson type III approximation of ``P(W >= threshold)``."""
        return self.tail_and_skewness(factor, weights, threshold)[0]

    def tail_and_skewness(self, factor, weights, threshold):
        """``P(W >= threshold)`` (Pearson type III) and the exact skewness of W."""
        m1, m2, m3 = self.raw_moments(factor, weights)
        var = m2 - m1**2
        if not var > 0:
            return float(threshold <= m1), 0.0
        sd = np.sqrt(var)
        skew = (m3 - 3 * m1 * m2 + 2 * m1**3) / sd**3
        z = (threshold - m1) / sd
        if abs(skew) < 1e-8:
            return float(stats.norm.sf(z)), skew
        shape = 4 / skew**2
        if skew > 0:
            return float(stats.gamma.sf(shape + z * np.sqrt(shape), shape)), skew
        return float(stats.gamma.cdf(shape - z * np.sqrt(shape), shape)), skew

"""The R-compatible random stream must reproduce R exactly.

Expected values were produced in R 4.5 (Mersenne-Twister, Rejection sampling).
"""

import numpy as np

from pyrda.rrandom import RRandom


def test_runif_matches_r():
    # set.seed(42); runif(5)
    expected = [0.914806043496355, 0.937075413297862, 0.286139534786344]
    expected += [0.830447626067325, 0.641745518893003]
    np.testing.assert_allclose(RRandom(42).runif(5), expected, rtol=1e-14)


def test_sample_int_matches_r():
    rng = RRandom(42)
    # set.seed(42); sample.int(10, 10); sample.int(24, 24)
    assert (rng.sample_int(10) + 1).tolist() == [1, 5, 10, 8, 2, 4, 6, 9, 7, 3]
    expected = [15, 7, 4, 5, 14, 18, 24, 3, 9, 16, 11, 22, 21, 12, 13, 2, 10, 17]
    expected += [19, 1, 6, 8, 20, 23]
    assert (rng.sample_int(24) + 1).tolist() == expected


def test_runif_continues_after_sampling():
    # set.seed(42); sample.int(10, 10); runif(1)  -> stream continues
    rng = RRandom(42)
    rng.sample_int(10)
    reference = RRandom(42).runif(11)[-1]
    assert rng.runif(1)[0] == reference


def test_shuffle_set_shapes_follow_permute():
    rng = RRandom(1)
    assert rng.shuffle_set(24, 999).shape == (999, 24)
    # n = 5: all 5! - 1 permutations (complete enumeration)
    perms = rng.shuffle_set(5, 999)
    assert perms.shape == (119, 5)
    assert len({tuple(p) for p in perms}) == 119
    # n = 7: 999 drawn from the 5039 enumerated permutations
    assert rng.shuffle_set(7, 999).shape == (999, 7)
    assert rng.shuffle_set(7, 9999).shape == (5039, 7)


def test_permutations_are_valid():
    perms = RRandom(3).shuffle_set(30, 50)
    assert all(sorted(p) == list(range(30)) for p in perms)

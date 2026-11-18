"""Bit-for-bit reproduction of the random numbers R uses for permutation tests.

R's default generator is Mersenne-Twister with ``sample.kind = "Rejection"``
(R >= 3.6). vegan's ``anova.cca()`` draws its permutations through
``permute::shuffleSet()``, which calls ``sample.int(n, n)`` once per permutation
(or enumerates all permutations when there are few of them). This module
replicates that chain, so ``RRandom(seed)`` produces the same permutation
matrices as ``set.seed(seed)`` followed by the same calls in R.
"""

import itertools
import math

import numpy as np

# permute::how() defaults used by vegan
_MINPERM = 5040
_CHECK_EPS = 0.5


class RRandom:
    """R-compatible random stream, equivalent to ``set.seed(seed)`` in R.

    Parameters
    ----------
    seed : int, optional
        Seed as passed to R's ``set.seed()``. When omitted a random seed is
        drawn, so results are random but the stream is still R's.
    """

    def __init__(self, seed=None):
        if seed is None:
            seed = int(np.random.SeedSequence().generate_state(1)[0])
        self.seed = int(seed)
        # R's set.seed(): scramble the seed with an LCG, then fill the
        # 625-word Mersenne-Twister state (mti followed by mt[624]).
        s = self.seed & 0xFFFFFFFF
        for _ in range(50):
            s = (69069 * s + 1) & 0xFFFFFFFF
        key = np.empty(625, dtype=np.uint32)
        for j in range(625):
            s = (69069 * s + 1) & 0xFFFFFFFF
            key[j] = s
        self._bitgen = np.random.MT19937()
        self._bitgen.state = {
            "bit_generator": "MT19937",
            "state": {"key": key[1:], "pos": 624},
        }
        self._buf = []
        self._pos = 0

    def _words16(self, size):
        # floor(unif_rand() * 65536) == (32-bit MT output) >> 16 exactly
        return (self._bitgen.random_raw(size) >> 16).tolist()

    def runif(self, size):
        """Uniform draws identical to R's ``runif(size)``."""
        # Flush the 16-bit word buffer back to the generator position first.
        self._rewind_buffer()
        u = self._bitgen.random_raw(size).astype(np.float64) * 2.3283064365386963e-10
        u[u <= 0.0] = 0.5 * 2.328306437080797e-10
        return u

    def _rewind_buffer(self):
        unused = len(self._buf) - self._pos
        if unused:
            self._bitgen.state = self._state_before_buffer
            self._bitgen.random_raw(self._pos)
        self._buf, self._pos = [], 0

    def _refill(self):
        self._state_before_buffer = self._bitgen.state
        self._buf = self._words16(8192)
        self._pos = 0

    def _sample(self, n, size, out, row):
        """R's do_sample() without replacement; writes 0-based values."""
        x = list(range(n))
        buf, pos = self._buf, self._pos
        m = n
        for i in range(size):
            bits = (m - 1).bit_length()  # == ceil(log2(m)) for m >= 1
            if bits > 15:
                # rbits() needs more than one 16-bit word: slow generic path
                self._buf, self._pos = buf, pos
                v = self._unif_index_wide(m, bits)
                buf, pos = self._buf, self._pos
            else:
                mask = (1 << bits) - 1
                while True:
                    if pos >= len(buf):
                        self._refill()
                        buf, pos = self._buf, 0
                    v = buf[pos] & mask
                    pos += 1
                    if v < m:
                        break
            out[row, i] = x[v]
            m -= 1
            x[v] = x[m]
        self._buf, self._pos = buf, pos

    def _unif_index_wide(self, m, bits):
        while True:
            v = 0
            for _ in range(bits // 16 + 1):
                if self._pos >= len(self._buf):
                    self._refill()
                v = 65536 * v + self._buf[self._pos]
                self._pos += 1
            v &= (1 << bits) - 1
            if v < m:
                return v

    def sample_int(self, n, size=None):
        """Equivalent of R's ``sample.int(n, size)`` (0-based result)."""
        size = n if size is None else size
        if size > n:
            raise ValueError("cannot take a sample larger than the population")
        out = np.empty((1, size), dtype=np.int64)
        self._sample(n, size, out, 0)
        return out[0]

    def shuffle_set(self, n, nperm):
        """Equivalent of ``permute::shuffleSet(n, control = how(nperm = nperm))``.

        Returns an ``(nperm_used, n)`` array of 0-based row permutations
        (row ``b`` means the data are reordered as ``Y[perm[b]]``). As in R,
        small designs switch to complete enumeration of all permutations.
        """
        num_pos = math.factorial(n)
        if nperm + _CHECK_EPS > num_pos - 1 or num_pos - 1 < _MINPERM + _CHECK_EPS:
            # complete enumeration (lexicographic order, observed one removed)
            perms = np.array(list(itertools.permutations(range(n)))[1:], dtype=np.int64)
            if len(perms) > nperm:
                perms = perms[self.sample_int(len(perms), nperm)]
            return perms
        out = np.empty((nperm, n), dtype=np.int64)
        for b in range(nperm):
            self._sample(n, n, out, b)
        return out

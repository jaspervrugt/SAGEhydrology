# ==========================================================================
# SAGEhydrology | Python source
# File: backend/matlab_random.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from matlab_random."""
import numpy as np

def randperm(n, seed):
    """Randperm."""
    if isinstance(n, (bool, np.bool_)) or int(n) != n or n < 0:
        raise ValueError('Permutation length must be a nonnegative integer.')
    if not np.isfinite(seed) or int(seed) != seed or (not 0 <= seed <= 2 ** 32 - 1):
        raise ValueError('MATLAB split seed must be an integer from 0 to 4294967295.')
    rng = np.random.RandomState(5489 if seed == 0 else int(seed))
    return np.argsort(rng.random_sample(int(n)), kind='stable')

def split_permutation(n, c):
    """Split permutation."""
    method = c.get('temporal_rng', 'matlab')
    if method == 'matlab':
        return randperm(n, c['split_seed'])
    if method == 'numpy':
        return np.random.RandomState(int(c['split_seed'])).permutation(n)
    raise ValueError('Unknown temporal permutation policy: ' + str(method))

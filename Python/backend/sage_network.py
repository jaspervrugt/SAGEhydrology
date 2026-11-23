# ==========================================================================
# SAGEhydrology | Python source
# File: backend/sage_network.py
#
# Headless SAGE numerical backend.
#
# Copyright (c) 2026 Jasper A. Vrugt
# University of California, Irvine | jasper@uci.edu
# ==========================================================================

"""Headless numerical routines extracted from sage_network."""
import numpy as np
from scipy.special import ndtr

class Network:

    def __init__(self, inputs, outputs, widths=(32,), activations=('tanh',), seed=0):
        """Initialize the evaluator or network state."""
        if not 1 <= len(widths) <= 5 or len(widths) != len(activations):
            raise ValueError('Choose one to five hidden layers and one activation per layer.')
        if any((int(n) != n or not 1 <= n <= 4096 for n in widths)):
            raise ValueError('Hidden widths must be integers from 1 to 4096.')
        if any((a not in ('tanh', 'relu', 'gelu') for a in activations)):
            raise ValueError('Hidden activations must be tanh, relu, or gelu.')
        rng = np.random.default_rng(seed)
        sizes = [inputs, *widths, outputs]
        self.activations = list(activations)
        self.weights = []
        self.biases = []
        for i, (ni, no) in enumerate(zip(sizes[:-1], sizes[1:])):
            limit = np.sqrt(6 / (ni + no))
            w = rng.uniform(-limit, limit, (ni, no))
            if i < len(widths) and activations[i] == 'relu':
                w = rng.normal(size=(ni, no)) * np.sqrt(2 / ni)
            self.weights.append(w)
            self.biases.append(np.zeros(no))

    @property
    def arrays(self):
        """Return the network parameter arrays."""
        return [v for pair in zip(self.weights, self.biases) for v in pair]

    def forward(self, attributes):
        """Evaluate the feedforward network and retain its layer cache."""
        a = np.asarray(attributes, float)
        cache = []
        for i, (w, b) in enumerate(zip(self.weights, self.biases)):
            z = a @ w + b
            cache.append((a, z))
            if i == len(self.weights) - 1:
                with np.errstate(over='ignore'):
                    a = 1 / (1 + np.exp(-z))
            elif self.activations[i] == 'tanh':
                a = np.tanh(z)
            elif self.activations[i] == 'gelu':
                a = z * ndtr(z)
            else:
                a = np.maximum(z, 0)
        return (a, (cache, a))

    def backward(self, cache, gradient):
        """Propagate the supplied sensitivities through the network."""
        layers, y = cache
        dz = np.asarray(gradient, float) * y * (1 - y)
        grads = []
        for i in range(len(self.weights) - 1, -1, -1):
            a, z = layers[i]
            grads[0:0] = [a.T @ dz, dz.sum(axis=0)]
            if i:
                dz = dz @ self.weights[i].T
                prev = layers[i - 1][1]
                activation = self.activations[i - 1]
                if activation == 'tanh':
                    dz *= 1 - np.tanh(prev) ** 2
                elif activation == 'gelu':
                    density = np.exp(-0.5 * prev ** 2) / np.sqrt(2 * np.pi)
                    dz *= ndtr(prev) + prev * density
                else:
                    dz *= prev > 0
        return grads

def standardize(raw, training_rows, mode='training'):
    """
    Fit imputations/scaling on training basins only; retain for
    inference.
    """
    if isinstance(mode, dict):
        from sage_attribute_preprocessing import fit
        return fit(raw, training_rows, mode)
    x = np.asarray(raw, float)
    if mode == 'matlab':
        if not np.isfinite(x).all():
            raise ValueError('Selected attributes contain missing/nonfinite values. Choose complete attributes or remove affected basins, as required by the MATLAB attribute workflow.')
        mean = x.mean(axis=0)
        scale = x.std(axis=0, ddof=1) if len(x) > 1 else np.ones(x.shape[1])
        scale[~np.isfinite(scale) | (scale == 0)] = 1
        return ((x - mean) / scale, dict(median=mean.tolist(), mean=mean.tolist(), scale=scale.tolist(), mode='matlab'))
    train = x[np.asarray(training_rows, int)]
    if not len(train) or np.any(~np.isfinite(train).any(axis=0)):
        raise ValueError('Each selected attribute needs a finite value in a training basin.')
    median = np.nanmedian(np.where(np.isfinite(train), train, np.nan), axis=0)
    filled = np.where(np.isfinite(x), x, median)
    mean = filled[training_rows].mean(axis=0)
    scale = filled[training_rows].std(axis=0)
    scale[scale == 0] = 1
    return ((filled - mean) / scale, dict(median=median.tolist(), mean=mean.tolist(), scale=scale.tolist()))

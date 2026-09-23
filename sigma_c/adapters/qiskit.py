# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Qiskit adapter — a parameter sweep to (sigma, O), plus the jitter channel from
per-shot memory. Opt-in: `pip install "sigma-c-framework[qiskit]"`.

The dial is a circuit `Parameter`; sweeping it and reading an observable's
expectation gives O(sigma), ready for `analyze()`.
"""
from __future__ import annotations
from typing import Optional, Sequence

import numpy as np


def observable_sweep(circuit, parameter, values, observable, *, estimator=None):
    """(sigma, O) from a parameterized Qiskit circuit swept over `values`.

    - `circuit`: a QuantumCircuit whose ONLY free parameter is `parameter`.
    - `parameter`: the swept `Parameter` (this is the dial sigma).
    - `values`: the dial values (1-D, ascending).
    - `observable`: a `SparsePauliOp` (or BaseOperator) whose expectation is O.
    - `estimator`: any BaseEstimatorV2; defaults to the exact `StatevectorEstimator`.

    Returns `(sigma, O)` as numpy arrays; feed them to `analyze(sigma, O)`.
    """
    from qiskit.primitives import StatevectorEstimator
    est = estimator or StatevectorEstimator()
    sigma = np.asarray(values, dtype=float)
    res = est.run([(circuit, observable, sigma.reshape(-1, 1))]).result()
    O = np.asarray(res[0].data.evs, dtype=float).ravel()
    return sigma, O


def analyze_sweep(circuit, parameter, values, observable, *, estimator=None,
                  window=None, **analyze_kwargs):
    """Convenience: `observable_sweep` then `analyze`. Returns a `Result`."""
    from sigma_c import analyze, bare
    sigma, O = observable_sweep(circuit, parameter, values, observable, estimator=estimator)
    return analyze(sigma, O, window=window or bare(), **analyze_kwargs)


def pair_parities_from_memory(memory: Sequence[str], pairs: Sequence[tuple]) -> np.ndarray:
    """Per-shot pair parities for the jitter channel, from a Sampler's per-shot
    bitstrings (`result.get_memory()` / a list of measurement strings, MSB-left as
    Qiskit returns them). `pairs` is a list of (i, j) qubit index pairs.

    Returns an (n_shots, n_pairs) array of +/-1 parities to feed
    `sigma_c.experimental.jitter.pair_parity_rank`.
    """
    shots = [s.replace(" ", "") for s in memory]
    n = len(shots[0]) if shots else 0
    out = np.empty((len(shots), len(pairs)), dtype=float)
    for r, bits in enumerate(shots):
        # Qiskit strings are MSB-left: qubit k is bits[n-1-k].
        b = [1 - 2 * int(bits[n - 1 - k]) for k in range(n)]   # 0->+1, 1->-1
        for c, (i, j) in enumerate(pairs):
            out[r, c] = b[i] * b[j]
    return out

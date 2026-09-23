# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Amazon Braket adapter — a FreeParameter sweep to (sigma, O). Opt-in:
`pip install "sigma-c-framework[braket]"`.

The dial is a Braket `FreeParameter`; sweeping it and reading a result type's value
(an expectation, exact at shots=0) gives O(sigma), ready for `analyze()`. This is the
path to the framework's real hardware validation (Rigetti/IQM run through Braket).
"""
from __future__ import annotations

import numpy as np


def observable_sweep(circuit, parameter_name, values, *, device=None):
    """(sigma, O) from a parameterized Braket `Circuit` swept over `values`.

    - `circuit`: a `Circuit` with a `FreeParameter` named `parameter_name` and a
      single result type attached (e.g. `.expectation(Observable.Z(), target=0)`).
    - `parameter_name`: the swept FreeParameter's name (the dial sigma).
    - `values`: the dial values (1-D, ascending).
    - `device`: any Braket device; defaults to the exact `LocalSimulator("braket_sv")`
      at shots=0.

    Returns `(sigma, O)` as numpy arrays; feed them to `analyze(sigma, O)`.
    """
    from braket.devices import LocalSimulator
    dev = device or LocalSimulator("braket_sv")
    sigma = np.asarray(values, dtype=float)
    O = np.array([
        float(dev.run(circuit, inputs={parameter_name: float(v)}, shots=0).result().values[0])
        for v in sigma
    ])
    return sigma, O


def analyze_sweep(circuit, parameter_name, values, *, device=None, window=None,
                  **analyze_kwargs):
    """Convenience: `observable_sweep` then `analyze`. Returns a `Result`."""
    from sigma_c import analyze, bare
    sigma, O = observable_sweep(circuit, parameter_name, values, device=device)
    return analyze(sigma, O, window=window or bare(), **analyze_kwargs)

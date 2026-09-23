# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
PennyLane adapter — a QNode sweep to (sigma, O), plus an analytic susceptibility.
Opt-in: `pip install "sigma-c-framework[pennylane]"`.

PennyLane's autodiff hands you dO/dsigma directly, so the susceptibility
chi = |sigma * dO/dsigma| can be computed EXACTLY (no finite difference) and passed
to `analyze(..., chi=chi)`.
"""
from __future__ import annotations

import numpy as np


def observable_sweep(qnode, values):
    """(sigma, O) from a scalar PennyLane QNode `qnode(theta)` swept over `values`.
    `qnode` must return a single expectation. Returns numpy arrays for `analyze()`."""
    sigma = np.asarray(values, dtype=float)
    O = np.array([float(qnode(float(v))) for v in sigma])
    return sigma, O


def susceptibility(qnode, values):
    """(sigma, O, chi) where chi = |sigma * dO/dsigma| is built from PennyLane's
    ANALYTIC gradient of the QNode (no finite difference). Feed via
    `analyze(sigma, O, chi=chi)`. The QNode must be differentiable (e.g. the autograd
    or a diff-enabled interface)."""
    import pennylane as qml
    sigma = np.asarray(values, dtype=float)
    O = np.array([float(qnode(float(v))) for v in sigma])
    grad = qml.grad(qnode)
    # The dial must be a trainable PennyLane array for the analytic gradient;
    # a plain float carries no requires_grad flag and qml.grad returns ().
    dO = np.array([float(grad(qml.numpy.array(v, requires_grad=True))) for v in sigma])
    chi = np.abs(sigma * dO)
    return sigma, O, chi


def analyze_sweep(qnode, values, *, use_autodiff_chi=False, window=None, **analyze_kwargs):
    """Convenience: sweep the QNode then `analyze`. With `use_autodiff_chi=True` the
    susceptibility is the analytic one (no finite difference). Returns a `Result`."""
    from sigma_c import analyze, bare
    win = window or bare()
    if use_autodiff_chi:
        sigma, O, chi = susceptibility(qnode, values)
        return analyze(sigma, O, chi=chi, window=win, **analyze_kwargs)
    sigma, O = observable_sweep(qnode, values)
    return analyze(sigma, O, window=win, **analyze_kwargs)

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Known-answer tests for the opt-in SDK adapters, each run END-TO-END against the SDK's
own local simulator (no hardware). The dial is a circuit parameter; an RY(theta)
sweep gives <Z> = cos(theta), so chi = |theta*dO/dtheta| = |theta*sin(theta)|, whose
peak is the solution of tan(theta) = -theta, sigma_c ~ 2.029. Each adapter must
recover that from real SDK output.

Every test skips cleanly if its SDK is not installed (the adapters are optional).
"""
import numpy as np
import pytest

from sigma_c import analyze, bare

SIGMA = np.linspace(0.1, 3.0, 200)
KNOWN_SIGMA_C = 2.029          # argmax |theta * sin(theta)| on (0, pi)


def _check(O, sigma):
    assert np.allclose(O, np.cos(SIGMA), atol=1e-6)           # the SDK gave <Z> = cos
    r = analyze(sigma, O, window=bare())
    assert r.sigma_c is not None
    assert abs(float(r.sigma_c) - KNOWN_SIGMA_C) < 0.05       # recovered from real data
    return r


def test_qiskit_adapter_recovers_known_sigma_c():
    pytest.importorskip("qiskit")
    from qiskit.circuit import QuantumCircuit, Parameter
    from qiskit.quantum_info import SparsePauliOp
    from sigma_c.adapters.qiskit import observable_sweep, analyze_sweep

    th = Parameter("theta"); qc = QuantumCircuit(1); qc.ry(th, 0)
    sigma, O = observable_sweep(qc, th, SIGMA, SparsePauliOp("Z"))
    _check(O, sigma)
    r = analyze_sweep(qc, th, SIGMA, SparsePauliOp("Z"))
    assert abs(float(r.sigma_c) - KNOWN_SIGMA_C) < 0.05


def test_pennylane_adapter_and_autodiff_chi():
    pytest.importorskip("pennylane")
    import pennylane as qml
    from sigma_c.adapters.pennylane import observable_sweep, susceptibility, analyze_sweep

    dev = qml.device("default.qubit", wires=1)

    @qml.qnode(dev)
    def qn(theta):
        qml.RY(theta, wires=0)
        return qml.expval(qml.PauliZ(0))

    sigma, O = observable_sweep(qn, SIGMA)
    _check(O, sigma)
    # PennyLane's ANALYTIC gradient gives chi with no finite difference.
    s2, O2, chi = susceptibility(qn, SIGMA)
    assert np.allclose(chi, np.abs(SIGMA * np.sin(SIGMA)), atol=1e-5)
    r = analyze_sweep(qn, SIGMA, use_autodiff_chi=True)
    assert abs(float(r.sigma_c) - KNOWN_SIGMA_C) < 0.05


def test_braket_adapter_recovers_known_sigma_c():
    pytest.importorskip("braket")
    from braket.circuits import Circuit, FreeParameter, Observable
    from sigma_c.adapters.braket import observable_sweep

    c = Circuit().ry(0, FreeParameter("theta")).expectation(Observable.Z(), target=0)
    sigma, O = observable_sweep(c, "theta", SIGMA)
    _check(O, sigma)


def test_qiskit_pair_parities_for_jitter():
    pytest.importorskip("qiskit")
    from sigma_c.adapters.qiskit import pair_parities_from_memory
    from sigma_c.experimental.jitter import pair_parity_covariance

    mem = ["000", "111", "000", "111"]          # MSB-left, 3 qubits
    p = pair_parities_from_memory(mem, [(0, 1), (0, 2)])
    assert p.shape == (4, 2)
    assert np.all(p == 1.0)                       # fully correlated -> parities all +1
    assert pair_parity_covariance(p).shape == (2, 2)

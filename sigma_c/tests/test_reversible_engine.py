# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Unit tests for the reversible/self-adjoint spectral engine (core/reversible).

This is a NEW pure module, not yet wired into analyze(); it is unit-tested here
directly. Once wired, the same behaviours are additionally exercised through the
public analyze() (the half-fix lesson). It returns FACTS, not verdicts.

Covers the separate-checker test cases at the engine level:
  - reversible chain given as a euclidean-NON-symmetric stochastic P -> verified
    (in the pi-weighted space), NOT rejected;
  - drift walk (circulant, euclidean-NORMAL but NOT reversible) -> none;
  - non-reversible chain -> none;
  - fast-Jordan / non-normal operator -> none (euclidean, not self-adjoint);
  - degenerate lambda_2 -> degenerate with a distinct next level;
plus negative-dominant (not embeddable), reducible leading (pi not unique),
and the spectrum-only "declared" tier.
"""
import math

import numpy as np
import pytest

from sigma_c.core.reversible import (
    analyze_reversible,
    stationary_by_detailed_balance,
)


# --- a reversible birth-death chain, euclidean-NON-symmetric -------------------
# pi = [1,2,1]/4 ; detailed balance holds ; P is not a symmetric matrix.
P_REV = np.array([
    [0.5,  0.5,  0.0],
    [0.25, 0.5,  0.25],
    [0.0,  0.5,  0.5],
])


def test_pi_from_detailed_balance_matches_and_is_consistent():
    pi, reversible, irreducible, reason = stationary_by_detailed_balance(P_REV, 1e-9)
    assert reversible and irreducible, reason
    assert np.allclose(pi, [0.25, 0.5, 0.25])
    # it is genuinely stationary too
    assert np.allclose(pi @ P_REV, pi)


def test_reversible_euclidean_nonsymmetric_is_verified_in_pi_space():
    # euclidean check would FAIL (P != P^T); auto/pi check must SUCCEED.
    assert not np.allclose(P_REV, P_REV.T)  # euclidean-non-symmetric, by construction
    r = analyze_reversible(operator=P_REV, T_star=1.0, inner_product="auto")
    assert r.selfadjoint == "verified", r.reason
    assert r.irreducible
    # eigenvalues are real and lambda_1 = 1
    assert r.lambda1 == pytest.approx(1.0, abs=1e-9)
    assert all(abs(v.imag) < 1e-12 for v in np.asarray(r.eigenvalues, dtype=complex))


def test_declared_pi_path_matches_auto():
    r = analyze_reversible(operator=P_REV, T_star=1.0, inner_product=[0.25, 0.5, 0.25])
    assert r.selfadjoint == "verified", r.reason
    assert r.lambda2 is not None and 0 < r.lambda2 < 1


def test_euclidean_check_rejects_the_reversible_nonsymmetric_chain():
    # The OLD (euclidean) notion of self-adjointness wrongly rejects it -- this test
    # pins WHY the inner_product declaration is necessary, not optional.
    r = analyze_reversible(operator=P_REV, T_star=1.0, inner_product="euclidean")
    assert r.selfadjoint == "none"
    assert "self-adjoint" in r.reason


# --- drift walk: circulant (euclidean-normal) but NOT reversible ---------------
def _drift_walk(n=3, p=0.7):
    P = np.zeros((n, n))
    for i in range(n):
        P[i, (i + 1) % n] = p
        P[i, (i - 1) % n] = 1 - p
    return P


def test_drift_walk_is_normal_but_not_reversible_none():
    P = _drift_walk()
    # euclidean-normal (circulant): P P^T == P^T P
    assert np.allclose(P @ P.T, P.T @ P)
    r = analyze_reversible(operator=P, T_star=1.0, inner_product="auto")
    assert r.selfadjoint == "none", r.reason
    assert "detailed balance" in r.reason or "reversible" in r.reason


def test_one_directional_edge_is_structural_none():
    P = np.array([[0.5, 0.5, 0.0],
                  [0.0, 0.5, 0.5],
                  [0.5, 0.0, 0.5]])  # a directed cycle: P_ij>0 xor P_ji>0
    r = analyze_reversible(operator=P, T_star=1.0, inner_product="auto")
    assert r.selfadjoint == "none"
    assert "asymmetr" in r.reason.lower() or "one-directional" in r.reason


# --- fast-Jordan / non-normal operator: not self-adjoint euclidean -------------
def test_fast_jordan_operator_is_none():
    A = np.array([[1.0, 0.0, 0.0],
                  [0.0, 0.5, 1.0],
                  [0.0, 0.0, 0.5 + 1e-6]])  # near-defective, non-normal
    r = analyze_reversible(operator=A, T_star=1.0, inner_product="euclidean")
    assert r.selfadjoint == "none", r.reason


# --- structured spectra via the symmetric-operator path ------------------------
def _sym_with_eigs(eigs):
    """A symmetric matrix with the given (real) eigenvalues via a random orthogonal Q."""
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((len(eigs), len(eigs))))
    return Q @ np.diag(eigs) @ Q.T


def test_degenerate_lambda2_gap_to_next_distinct_level():
    A = _sym_with_eigs([1.0, 0.5, 0.5, 0.1])
    r = analyze_reversible(operator=A, T_star=1.0, inner_product="euclidean")
    assert r.selfadjoint == "verified"
    assert r.lambda2 == pytest.approx(0.5, abs=1e-9)
    assert r.degenerate and r.lambda2_multiplicity == 2
    assert r.lambda_next == pytest.approx(0.1, abs=1e-9)
    # delta uses the next DISTINCT level, not the partner
    assert r.delta == pytest.approx(math.log(0.5) - math.log(0.1), abs=1e-9)


def test_negative_dominant_eigenvalue_not_embeddable():
    A = _sym_with_eigs([1.0, -0.9, 0.3])
    r = analyze_reversible(operator=A, T_star=1.0, inner_product="euclidean")
    assert r.selfadjoint == "verified"
    assert r.lambda2 == pytest.approx(0.3, abs=1e-9)   # largest BELOW 1, not |-0.9|
    assert r.embeddable is False                       # -0.9 dominates in modulus


def test_reducible_leading_eigenvalue_not_irreducible():
    A = _sym_with_eigs([1.0, 1.0, 0.3])
    r = analyze_reversible(operator=A, T_star=1.0, inner_product="euclidean")
    assert r.irreducible is False
    assert "degenerate" in r.reason or "not unique" in r.reason or "reducible" in r.reason


def test_clean_gap_gives_tau_and_band():
    A = _sym_with_eigs([1.0, 0.5, 0.1])
    r = analyze_reversible(operator=A, T_star=1.0, inner_product="euclidean")
    assert r.tau_spectral == pytest.approx(-1.0 / math.log(0.5), abs=1e-9)
    assert r.band_rel is not None and r.band_rel > 0    # eps/(1-lambda2), small but present
    assert r.embeddable and r.irreducible and not r.degenerate


# --- spectrum-only tier: self-adjointness DECLARED, not verified ---------------
def test_spectrum_only_is_declared():
    r = analyze_reversible(spectrum=[1.0, 0.5, 0.1], T_star=1.0, selfadjoint_declared=True)
    assert r.selfadjoint == "declared"
    assert r.lambda2 == pytest.approx(0.5, abs=1e-9)
    assert any("DECLARED" in n for n in r.notes)


def test_spectrum_only_complex_is_none():
    r = analyze_reversible(spectrum=[1.0, 0.5 + 0.3j, 0.1], T_star=1.0)
    assert r.selfadjoint == "none"
    assert "non-real" in r.reason or "self-adjoint" in r.reason


# --- hardening (separate-checker requests) -------------------------------------

def test_metastable_chain_is_verified_not_tripped_by_symmetry_tolerance():
    # 1 - lambda_2 ~ 4e-9 (deeply metastable). The eps*||S||*n symmetry tolerance
    # must NOT reject it, and the tiny gap must NOT read as below-band.
    a, b = 1e-9, 3e-9  # 2-state reversible; lambda_2 = 1-(a+b), pi=[b,a]/(a+b)
    P = np.array([[1 - a, a], [b, 1 - b]])
    r = analyze_reversible(operator=P, T_star=1.0, inner_product="auto")
    assert r.selfadjoint == "verified", r.reason
    assert r.lambda2 == pytest.approx(1 - (a + b), abs=1e-12)
    assert r.gap_below_band is False
    assert np.allclose(r.pi, [0.75, 0.25])


def test_weak_edge_pi_is_path_independent():
    # asymmetric reversible triangle with a WEAK edge 0-2 (conductance 1e-15). The
    # spanning-tree walk may route pi through the weak edge; reversibility makes pi
    # path-independent, so the result must be exact regardless (max-spanning-tree
    # over min(P_ij,P_ji) is therefore unnecessary for CORRECTNESS).
    pi_true = np.array([0.1, 0.6, 0.3])
    c01 = c12 = 0.05
    c02 = 1e-15
    P = np.zeros((3, 3))
    P[0, 1], P[1, 0] = c01 / pi_true[0], c01 / pi_true[1]
    P[1, 2], P[2, 1] = c12 / pi_true[1], c12 / pi_true[2]
    P[0, 2], P[2, 0] = c02 / pi_true[0], c02 / pi_true[2]
    for i in range(3):
        P[i, i] = 1 - P[i].sum()
    pi, reversible, irreducible, _ = stationary_by_detailed_balance(P, 1e-9)
    assert reversible and irreducible
    assert np.max(np.abs(pi - pi_true)) < 1e-12   # exact despite the weak-edge route
    r = analyze_reversible(operator=P, T_star=1.0, inner_product="auto")
    assert r.selfadjoint == "verified"

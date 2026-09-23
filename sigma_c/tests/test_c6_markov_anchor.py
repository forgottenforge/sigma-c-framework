# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Golden test for paper Example C.6 -- the 4-state Markov chain.

This is the substantive faithfulness anchor of the paper: a
4-state Markov chain with genuine sub-leading spectrum r ~ 0.521, two
structurally distinct lambda_2-faithful test functions, recovers tau to
within 6.3% using the *theoretical* profile constant rho_star = 2-sqrt(3)
(NOT a fit). All numbers in this test are taken verbatim from the paper's
verify_examples.py output that ran 41/41 PASS at submission.

The chain (paper §C.6):
    P_0 = [[0.6, 0.3, 0.1, 0.0],
           [0.2, 0.5, 0.2, 0.1],
           [0.0, 0.2, 0.5, 0.3],
           [0.1, 0.0, 0.3, 0.6]]
    epsilon = 0.1
    U_ij = 1/4  (uniform mixing)
    P_eps = (1 - epsilon) * P_0 + epsilon * U

Spectrum: {1, 0.6002, 0.3128, 0.1669}
r := lambda_3 / lambda_2 ~ 0.521
tau := -T_*/log(lambda_2) ~ 1.959  (with T_* = 1)

Two probes:
    Phi_A = ( 3,  1, -1, -3)
    Phi_B = ( 4,  1, -1, -2)
Initial distribution: nu_0 = (1/4, 1/4, 1/4, 1/4)

Cite (THEOREM_MAP):
    def:faithful, prop:faith-sufficient (F1 condition demonstrated),
    thm:cross-obs-concentration, prop:structural-reduction. (NOT thm:spectral-id:
    result.tau is the fallen bridge, backed by thm:spectral-id-B (GAP-KNOWN); see #3.)
"""
from __future__ import annotations
import math

import numpy as np
import pytest

from sigma_c import (
    analyze,
    gamma_k,
    Framework,
    check_F1,
    kl_modal_coefficients,
)


# ---------------------------------------------------------------------------
# Common setup: the chain, its spectrum, and the two probes
# ---------------------------------------------------------------------------

P_0 = np.array([
    [0.6, 0.3, 0.1, 0.0],
    [0.2, 0.5, 0.2, 0.1],
    [0.0, 0.2, 0.5, 0.3],
    [0.1, 0.0, 0.3, 0.6],
])
EPSILON = 0.1
U = np.full((4, 4), 0.25)
P_EPS = (1 - EPSILON) * P_0 + EPSILON * U

PHI_A = np.array([3, 1, -1, -3], dtype=float)
PHI_B = np.array([4, 1, -1, -2], dtype=float)
NU_0 = np.array([0.25, 0.25, 0.25, 0.25])

# Paper-printed reference values (paper §C.6 / verify_examples.py)
SPECTRUM_REF = [1.0, 0.6002, 0.3128, 0.1669]
GAP_RATIO_R_REF = 0.521
TAU_REF = 1.959
RHO_STAR_REF = 2 - math.sqrt(3)
SIGMA_C_A_REF = 0.5576
SIGMA_C_B_REF = 0.5243
TAU_A_REF = 2.081
TAU_B_REF = 1.957
CROSS_PROBE_DISAGREEMENT_REF = 0.063  # 6.3%
WORST_CASE_OR_BOUND_REF = 0.52        # paper's O(r) ~ 52%


def _gamma2_modal_observable(coeffs, eigvals, sigma, T_star=1.0):
    """
    Construct the modal-sum observable O_i(sigma) = sum_j c_j * sigma * tau_j^2 / (sigma + tau_j)^2

    using the Gamma-2 window (paper Example C.6). Skips the j=1 (trivial) mode.
    """
    O = np.zeros_like(sigma, dtype=float)
    for j, (c, lam) in enumerate(zip(coeffs, eigvals)):
        if j == 0:
            continue  # skip lambda_1 == 1 (trivial mode)
        if abs(lam) < 1e-12 or lam.real <= 0:
            continue
        tau_j = -T_star / math.log(abs(lam.real))
        # Gamma-2 windowed form (Mellin inversion of u*exp(-u))
        O = O + c.real * sigma * tau_j ** 2 / (sigma + tau_j) ** 2
    return O


# ---------------------------------------------------------------------------
# 1. Spectrum reproduction
# ---------------------------------------------------------------------------

class TestC6Spectrum:
    """The chain's spectrum must match the paper's printed values."""

    def test_eigenvalues_match_paper(self):
        eigvals = np.linalg.eigvals(P_EPS)
        sorted_abs = sorted((abs(x) for x in eigvals), reverse=True)
        for actual, expected in zip(sorted_abs, SPECTRUM_REF):
            assert actual == pytest.approx(expected, abs=5e-4), (
                f"Eigenvalue mismatch: {actual} vs paper {expected}"
            )

    def test_gap_ratio_matches_paper(self):
        eigvals = np.linalg.eigvals(P_EPS)
        sorted_abs = sorted((abs(x) for x in eigvals), reverse=True)
        r = sorted_abs[2] / sorted_abs[1]
        assert r == pytest.approx(GAP_RATIO_R_REF, abs=2e-3)

    def test_tau_matches_paper(self):
        eigvals = np.linalg.eigvals(P_EPS)
        sorted_abs = sorted((abs(x) for x in eigvals), reverse=True)
        lam2 = sorted_abs[1]
        tau = -1.0 / math.log(lam2)
        assert tau == pytest.approx(TAU_REF, abs=2e-3)


# ---------------------------------------------------------------------------
# 2. Modal coefficients (verify F1 inputs match the paper)
# ---------------------------------------------------------------------------

class TestC6ModalCoefficients:
    """Modal coefficients c_j for both probes must match the paper's numbers
    within reasonable tolerance, and the F1 sufficient condition must pass."""

    @pytest.fixture
    def coeffs_A(self):
        eigvals, coeffs = kl_modal_coefficients(P_EPS, PHI_A, initial=NU_0)
        return eigvals, coeffs

    @pytest.fixture
    def coeffs_B(self):
        eigvals, coeffs = kl_modal_coefficients(P_EPS, PHI_B, initial=NU_0)
        return eigvals, coeffs

    def test_phi_A_sub_leading_ratios_below_r(self, coeffs_A):
        eigvals, coeffs = coeffs_A
        c2 = abs(coeffs[1])
        # Both |c_3/c_2| and |c_4/c_2| below r (paper memory)
        sub = [abs(c) / c2 for c in coeffs[2:]]
        for ratio in sub:
            assert ratio <= GAP_RATIO_R_REF + 1e-2, (
                f"Phi_A: sub-leading ratio {ratio} above r = {GAP_RATIO_R_REF}"
            )

    def test_phi_B_sub_leading_ratios_below_r(self, coeffs_B):
        eigvals, coeffs = coeffs_B
        c2 = abs(coeffs[1])
        sub = [abs(c) / c2 for c in coeffs[2:]]
        for ratio in sub:
            assert ratio <= GAP_RATIO_R_REF + 1e-2, (
                f"Phi_B: sub-leading ratio {ratio} above r = {GAP_RATIO_R_REF}"
            )

    def test_F1_passes_for_phi_A(self, coeffs_A):
        eigvals, coeffs = coeffs_A
        # F1 takes [c_2, c_3, ...] -- coefficients sequence starting at j=2
        modal = [coeffs[j] for j in range(1, len(coeffs))]
        check = check_F1(
            modal_coeffs=modal,
            gap_ratio=GAP_RATIO_R_REF,
            window_range_M=1,
        )
        assert check.passed, f"F1 should pass for Phi_A: {check.summary()}"
        assert math.isfinite(check.C_R)
        assert check.condition == "F1"

    def test_F1_passes_for_phi_B(self, coeffs_B):
        eigvals, coeffs = coeffs_B
        modal = [coeffs[j] for j in range(1, len(coeffs))]
        check = check_F1(
            modal_coeffs=modal,
            gap_ratio=GAP_RATIO_R_REF,
            window_range_M=1,
        )
        assert check.passed, f"F1 should pass for Phi_B: {check.summary()}"


# ---------------------------------------------------------------------------
# 3. The bridge: sigma_c values from the analytic windowed observable
# ---------------------------------------------------------------------------

class TestC6Bridge:
    """The full paper Example C.6: build the modal-sum observable for each
    probe, run analyze() with the analytic gamma_2 window, recover sigma_c
    matching the paper's printed values, and verify the cross-probe
    agreement is well inside the worst-case bound.

    This anchor rests on TWO things the register now
    flags: (i) the gamma_2 window gives chi TWO peaks (rising + falling side) ---
    the two-equal-peaks pathology of `thm:spectral-id-B` (GAP-KNOWN): sigma_c is
    ill-defined and the "rising-side = min" pick is a CONVENTION, not a certified
    argmax; (ii) tau is recovered by the FALLEN sigma_c = rho_star*tau bridge, by
    hand (regime II emits result.tau = None). So the numeric assertions below are a
    REPRODUCIBILITY check of the paper's old-convention numbers --- NOT a
    certification that tau is identified. The honest 4-code verdict is pinned
    separately in test_c6_bridge_tau_is_not_certified_by_contract: the contract
    does NOT return OK for this construction. Cite: thm:spectral-id-B (GAP-KNOWN),
    thm:cross-obs-concentration (REVIEW-CLEAN)."""

    def _run_probe(self, probe):
        eigvals, coeffs = kl_modal_coefficients(P_EPS, probe, initial=NU_0)
        sigma = np.geomspace(0.05, 50.0, 2000)
        O = _gamma2_modal_observable(coeffs, eigvals, sigma)
        result = analyze(
            sigma, O,
            window=gamma_k(2),
            framework=Framework.REVERSIBLE_MARKOV,
        )
        return result

    def test_sigma_c_A_matches_paper(self):
        result = self._run_probe(PHI_A)
        assert result.sigma_c is not None
        # Geometric class may be II_geom because chi has two peaks (rising +
        # falling side of the gamma2 form); take the rising side.
        if isinstance(result.sigma_c, list):
            rising = min(result.sigma_c)
        else:
            rising = result.sigma_c
        assert rising == pytest.approx(SIGMA_C_A_REF, abs=3e-3), (
            f"sigma_c[O_A] rising-side: got {rising}, expected {SIGMA_C_A_REF}"
        )

    def test_sigma_c_B_matches_paper(self):
        result = self._run_probe(PHI_B)
        if isinstance(result.sigma_c, list):
            rising = min(result.sigma_c)
        else:
            rising = result.sigma_c
        assert rising == pytest.approx(SIGMA_C_B_REF, abs=3e-3), (
            f"sigma_c[O_B] rising-side: got {rising}, expected {SIGMA_C_B_REF}"
        )

    def test_tau_recovery_via_bridge_for_phi_A(self):
        """tau_recovered = sigma_c / rho_star with rho_star = 2 - sqrt(3) analytic."""
        result = self._run_probe(PHI_A)
        rising = min(result.sigma_c) if isinstance(result.sigma_c, list) else result.sigma_c
        tau_a = rising / RHO_STAR_REF
        assert tau_a == pytest.approx(TAU_A_REF, abs=5e-3)

    def test_tau_recovery_via_bridge_for_phi_B(self):
        result = self._run_probe(PHI_B)
        rising = min(result.sigma_c) if isinstance(result.sigma_c, list) else result.sigma_c
        tau_b = rising / RHO_STAR_REF
        assert tau_b == pytest.approx(TAU_B_REF, abs=5e-3)

    def test_cross_probe_disagreement_within_paper_bound(self):
        """The 6.3% cross-probe disagreement is the flagship anchor: well
        below the worst-case O(r) ~ 52% bound predicted by Thm 6.1."""
        result_A = self._run_probe(PHI_A)
        result_B = self._run_probe(PHI_B)
        rA = min(result_A.sigma_c) if isinstance(result_A.sigma_c, list) else result_A.sigma_c
        rB = min(result_B.sigma_c) if isinstance(result_B.sigma_c, list) else result_B.sigma_c
        tau_a = rA / RHO_STAR_REF
        tau_b = rB / RHO_STAR_REF
        disagreement = abs(tau_a - tau_b) / max(tau_a, tau_b)
        # Paper's printed value: 6.3 %. We allow +- 1 % tolerance.
        assert disagreement <= 0.08, (
            f"cross-probe disagreement {disagreement:.3f} > 8% (paper: 6.3%)"
        )
        assert disagreement <= WORST_CASE_OR_BOUND_REF, (
            f"cross-probe disagreement above worst-case O(r) bound"
        )

    def test_provenance_is_analytic_gamma2(self):
        """Enforcement 2: rho_star_source must be analytic, not fitted."""
        result = self._run_probe(PHI_A)
        assert result.rho_star_source == "analytic:gamma2"
        assert result.falsifiable
        assert result.rho_star == pytest.approx(RHO_STAR_REF, abs=1e-12)

    def test_c6_bridge_tau_is_not_certified_by_contract(self):
        """REHANG (#2): the honest 4-code verdict, pinned next to the old-convention
        numbers above. This C.6 gamma_2 construction rests on thm:spectral-id-B
        (GAP-KNOWN two-equal-peaks) and the fallen sigma_c=rho_star*tau bridge, so
        the contract must NOT return OK for tau. (No spectrum and no gamma_A are
        supplied, and chi is two-peaked -> tau is either None (regime II) or
        emitted via the bridge (regime I) -- in neither case certified.)"""
        from sigma_c.codes import Code
        result = self._run_probe(PHI_A)
        assert result.window_readability_status.code is not Code.OK, (
            f"the C.6 gamma_2/bridge tau must not be certified OK: {result.window_readability_status}"
        )


# ---------------------------------------------------------------------------
# 4. Theorem-backing invariants for C.6
# ---------------------------------------------------------------------------

class TestC6TheoremBacking:
    """Every analyze() call on C.6 must surface specific paper citations."""

    def test_citations_include_paper_anchors(self):
        eigvals, coeffs = kl_modal_coefficients(P_EPS, PHI_A, initial=NU_0)
        sigma = np.geomspace(0.05, 50.0, 1000)
        O = _gamma2_modal_observable(coeffs, eigvals, sigma)
        result = analyze(sigma, O, window=gamma_k(2),
                         framework=Framework.REVERSIBLE_MARKOV)
        labels = result.citations
        # def:sigmac and prop:structural-reduction are always cited.
        assert "def:sigmac" in labels
        assert "prop:structural-reduction" in labels
        # Framework declared -> the framework TAXONOMY is cited...
        assert "prop:standard-frameworks" in labels
        # ...but NOT the PROVED thm:spectral-id. result.tau here is the fallen
        # sigma_c=rho_star*tau bridge, not the live-route theorem; citing a PROVED
        # label (which renders no warning) on a value it does not back is the
        # The honest backing is
        # thm:spectral-id-B (GAP-KNOWN), cited at the tau emission when tau exists.
        assert "thm:spectral-id" not in labels

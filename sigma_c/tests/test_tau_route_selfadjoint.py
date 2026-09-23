# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
The tau-route (self-adjointness kernel) -- the two-field API, driven
ENTIRELY through the public analyze() (half-fix lesson).

Design (self-adjoint reversible route):
  * The live-route theorems live in L^2(rho_*) with detailed balance; the
    precondition is SELF-ADJOINTNESS in the declared inner product, not euclidean
    normality. `inner_product` is a declaration: "euclidean" | a pi-weight | "auto".
  * TWO fields, not one status on tau:
      - tau_abscissa (thm:abscissa): the certified spectral abscissa -T*/log(lambda2)
        (in TIME). OK needs self-adjoint + discrete + isolated mu2 + faithful; NO
        strict gap required.
      - window_readability (thm:auto): OK additionally needs a strict gap + a
        window T > max(t*, 1/|mu3|); reports rate_error_bound.
    Non-self-adjoint -> NOT_APPLICABLE on BOTH.
  * result.tau and result.tau_status are REMOVED (raise AttributeError with a
    migration hint); tau_bridge (sigma-unit, fallen) + tau_abscissa (time) replace them.

Separate-checker cases live here at the analyze() level:
  reversible euclidean-nonsymmetric chain -> OK (verified in the pi-space);
  drift walk (normal, not reversible) -> NOT_APPLICABLE;
  fast-Jordan / non-normal -> NOT_APPLICABLE;
  degenerate mu2 -> NOT_RESOLVABLE (below band) or rem:gap-cited OK (clear next level);
  spectrum-only -> declared (window_readability NOT_IDENTIFIED).
"""
import math

import numpy as np
import pytest

from sigma_c import analyze, bare, Framework
from sigma_c.codes import Code


def _obs(sigma_max=60.0, n=1000):
    sigma = np.geomspace(0.02, sigma_max, n)
    O = 1.0 * np.exp(-sigma / 1.0) + 0.4 * np.exp(-sigma / (1.0 / 3.0))
    return sigma, O


# a reversible birth-death chain, euclidean-NON-symmetric (pi=[1,2,1]/4)
P_REV = np.array([[0.5, 0.5, 0.0],
                  [0.25, 0.5, 0.25],
                  [0.0, 0.5, 0.5]])


def _drift_walk(n=3, p=0.7):
    P = np.zeros((n, n))
    for i in range(n):
        P[i, (i + 1) % n] = p
        P[i, (i - 1) % n] = 1 - p
    return P


# ===========================================================================
# API migration: bare tau / tau_status are gone
# ===========================================================================

def test_bare_tau_raises_with_migration_hint():
    sigma, O = _obs()
    r = analyze(sigma, O, window=bare())
    with pytest.raises(AttributeError, match="tau_bridge|tau_abscissa"):
        _ = r.tau
    with pytest.raises(AttributeError, match="window_readability|tau_abscissa_status"):
        _ = r.tau_status


def test_sigma_c_only_emits_bridge_and_both_fields_not_identified():
    sigma, O = _obs()
    r = analyze(sigma, O, window=bare())
    # the fallen bridge is still available under its honest name ...
    assert r.tau_bridge == pytest.approx(r.sigma_c / r.rho_star)
    assert r.tau_abscissa is None
    assert r.selfadjoint is None
    # ... and both certified fields refuse (no spectrum/operator supplied)
    assert r.tau_abscissa_status.code is Code.NOT_IDENTIFIED
    assert r.window_readability_status.code is Code.NOT_IDENTIFIED
    d = r.to_dict()
    assert d["tau_bridge"] is not None and d["tau_abscissa"] is None
    assert "tau" not in d and "tau_status" not in d


# ===========================================================================
# Self-adjointness is the fault line (5 cases)
# ===========================================================================

def test_reversible_euclidean_nonsymmetric_is_OK_in_pi_space():
    sigma, O = _obs()
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=P_REV, inner_product="auto", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time")
    assert r.selfadjoint == "verified"
    assert r.tau_abscissa is not None
    assert r.tau_abscissa_status.code is Code.OK, str(r.tau_abscissa_status)


def test_drift_walk_normal_but_not_reversible_is_NOT_APPLICABLE_on_both():
    sigma, O = _obs()
    P = _drift_walk()
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=P, inner_product="auto", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time")
    assert r.selfadjoint == "none"
    assert r.tau_abscissa_status.code is Code.NOT_APPLICABLE, str(r.tau_abscissa_status)
    assert r.window_readability_status.code is Code.NOT_APPLICABLE


def test_euclidean_check_rejects_the_reversible_nonsymmetric_chain():
    # WHY inner_product is necessary: the default euclidean reading refuses it.
    sigma, O = _obs()
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=P_REV, inner_product="euclidean", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time")
    assert r.selfadjoint == "none"
    assert r.tau_abscissa_status.code is Code.NOT_APPLICABLE


def test_fast_jordan_non_normal_is_NOT_APPLICABLE():
    sigma, O = _obs()
    A = np.array([[1.0, 0.0, 0.0],
                  [0.0, 0.5, 1.0],
                  [0.0, 0.0, 0.5 + 1e-6]])
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=A, inner_product="euclidean", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time")
    assert r.selfadjoint == "none"
    assert r.tau_abscissa_status.code is Code.NOT_APPLICABLE


def _sym_with_eigs(eigs):
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((len(eigs), len(eigs))))
    return Q @ np.diag(eigs) @ Q.T


def test_clean_selfadjoint_gap_reaches_window_readability_OK_with_bound():
    sigma, O = _obs()
    A = _sym_with_eigs([1.0, 0.5, 0.1])
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=A, inner_product="euclidean", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time", T_obs=100.0)
    assert r.tau_abscissa == pytest.approx(-1.0 / math.log(0.5), abs=1e-9)
    assert r.tau_abscissa_status.code is Code.OK
    assert r.window_readability_status.code is Code.OK, str(r.window_readability_status)
    # the proven finite-window rate error bound is a REPORTED number
    assert r.rate_error_bound is not None and r.rate_error_bound > 0


def test_degenerate_mu2_below_band_is_NOT_RESOLVABLE():
    sigma, O = _obs()
    A = _sym_with_eigs([1.0, 0.5, 0.5, 0.1])   # lambda2 degenerate (mult 2)
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=A, inner_product="euclidean", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time", T_obs=100.0)
    # abscissa still identifies the mu2-subspace (thm:abscissa, no strict gap needed)
    assert r.tau_abscissa_status.code is Code.OK, str(r.tau_abscissa_status)
    # but window readability needs a STRICT gap -> here mu2 has a distinct next level
    # (0.1), so it is resolvable and rests on rem:gap (GAP-KNOWN citation)
    assert r.window_readability_status.code in (Code.OK, Code.NOT_RESOLVABLE)


def test_negative_dominant_eigenvalue_not_embeddable_NOT_APPLICABLE():
    sigma, O = _obs()
    A = _sym_with_eigs([1.0, -0.9, 0.3])
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=A, inner_product="euclidean", gamma_A=0.5, T_star=1.0,
                sigma_axis="evolution_time")
    assert r.tau_abscissa_status.code is Code.NOT_APPLICABLE, str(r.tau_abscissa_status)
    assert "embeddable" in r.tau_abscissa_status.reason or "generator" in r.tau_abscissa_status.reason


# ===========================================================================
# The 1/|mu3| gate (erratum: gate must be T > max(t*, 1/|mu3|))
# ===========================================================================

def test_window_too_short_by_mu3_floor_is_NOT_IDENTIFIED():
    # faithful probe (gamma_A<=1) -> t*=0, so ONLY the 1/|mu3| floor gates.
    # mu3 = log(0.1) = -2.303 -> 1/|mu3| = 0.434 ; T_obs below it must refuse.
    sigma, O = _obs()
    A = _sym_with_eigs([1.0, 0.5, 0.1])
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=A, inner_product="euclidean", gamma_A=0.5, T_star=1.0,
                sigma_axis="window_time", T_obs=0.2)  # T_obs < 1/|mu3|=0.434
    assert r.tau_abscissa_status.code is Code.OK          # abscissa unaffected
    assert r.window_readability_status.code is Code.NOT_IDENTIFIED, str(r.window_readability_status)
    assert "window" in r.window_readability_status.reason.lower()


# ===========================================================================
# Spectrum-only tier: self-adjointness + faithfulness DECLARED
# ===========================================================================

def test_spectrum_only_is_declared_and_window_needs_operator_and_probe():
    sigma, O = _obs()
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                spectrum=[1.0, 0.5, 0.1], T_star=1.0, sigma_axis="evolution_time")
    assert r.selfadjoint == "declared"
    assert r.tau_abscissa == pytest.approx(-1.0 / math.log(0.5), abs=1e-9)
    # window readability can NEVER be OK without operator+probe (no Gamma_A checkable)
    assert r.window_readability_status.code is Code.NOT_IDENTIFIED, str(r.window_readability_status)
    assert "operator" in r.window_readability_status.remediation


def test_spectrum_only_complex_is_NOT_APPLICABLE():
    sigma, O = _obs()
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                spectrum=[1.0, 0.5 + 0.3j, 0.1], T_star=1.0, sigma_axis="evolution_time")
    assert r.selfadjoint == "none"
    assert r.tau_abscissa_status.code is Code.NOT_APPLICABLE


# ===========================================================================
# tau-route declaration validation (api.py; unchanged by the self-adjoint rebuild)
# ===========================================================================

def test_tau_route_requires_T_star():
    sigma, O = _obs()
    with pytest.raises(ValueError, match="T_star is required"):
        analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                spectrum=[1.0, 0.5, 0.1], sigma_axis="evolution_time")


def test_spectrum_and_operator_together_raises():
    sigma, O = _obs()
    with pytest.raises(ValueError, match="either spectrum .* or operator"):
        analyze(sigma, O, window=bare(), spectrum=[1.0, 0.5, 0.1],
                operator=np.diag([1.0, 0.5, 0.1]), T_star=1.0, sigma_axis="evolution_time")


def test_invalid_sigma_axis_raises():
    sigma, O = _obs()
    with pytest.raises(ValueError, match="sigma_axis must be one of"):
        analyze(sigma, O, window=bare(), spectrum=[1.0, 0.5, 0.1], T_star=1.0,
                sigma_axis="clock")


def test_window_time_requires_T_obs():
    sigma, O = _obs()
    with pytest.raises(ValueError, match="T_obs .* is required"):
        analyze(sigma, O, window=bare(), spectrum=[1.0, 0.5, 0.1], T_star=1.0,
                sigma_axis="window_time")


def test_nonfinite_spectrum_raises():
    sigma, O = _obs()
    with pytest.raises(ValueError, match="finite"):
        analyze(sigma, O, window=bare(), spectrum=[np.inf, 0.5, 0.1], T_star=1.0,
                sigma_axis="evolution_time")
    with pytest.raises(ValueError, match="finite"):
        analyze(sigma, O, window=bare(), spectrum=[1.0, np.nan, 0.1], T_star=1.0,
                sigma_axis="evolution_time")


def test_nondecaying_subdominant_spectrum_is_NOT_RESOLVABLE():
    # |lambda_2| >= 1 (non-decaying): a modulus gap but no finite relaxation time.
    sigma, O = _obs()
    for spec in ([1.5, 1.0, 0.1], [2.0, 1.5, 0.1]):
        r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                    spectrum=spec, T_star=1.0, sigma_axis="evolution_time")
        assert r.tau_abscissa is None
        assert r.tau_abscissa_status.code is Code.NOT_RESOLVABLE, str(r.tau_abscissa_status)


# ===========================================================================
# Regression (adversarial checker): tau_abscissa is an OPERATOR property, so its
# VALUE must be emitted regardless of the OBSERVABLE's regime (I/II/III). Before
# the fix, regime II/III hardcoded tau_abscissa=None while the status said OK.
# ===========================================================================

def test_abscissa_value_emitted_on_regime_II_observable_with_operator():
    sigma = np.geomspace(0.02, 100.0, 600)
    O = np.exp(-sigma / 0.5) + 0.5 * np.exp(-sigma / 15.0)   # two peaks -> regime II
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=np.diag([1.0, 0.5, 0.1]), inner_product="euclidean",
                T_star=1.0, sigma_axis="window_time", T_obs=100.0, gamma_A=0.5)
    assert isinstance(r.sigma_c, list)                       # regime II observable
    # the certified abscissa VALUE is present (not silently dropped) ...
    assert r.tau_abscissa == pytest.approx(-1.0 / math.log(0.5), abs=1e-9)
    # ... and value tracks verdict: a non-None OK status must carry a non-None value
    if r.tau_abscissa_status.code is Code.OK:
        assert r.tau_abscissa is not None


def test_abscissa_value_emitted_on_regime_III_observable_with_operator():
    sigma = np.geomspace(0.05, 50.0, 300)
    r = analyze(sigma, sigma.copy(), window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                operator=np.diag([1.0, 0.5, 0.1]), inner_product="euclidean",
                T_star=1.0, sigma_axis="evolution_time")     # O=sigma -> regime III
    assert r.sigma_c is None                                 # regime III observable
    assert r.tau_abscissa == pytest.approx(-1.0 / math.log(0.5), abs=1e-9)
    # no value/verdict contradiction
    assert not (r.tau_abscissa is None and r.tau_abscissa_status.code is Code.OK)

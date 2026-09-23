# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
FIRST test of the 4-code contract: INVARIANCE on the golden anchors.

The sigma_c code is a FORMALISATION of the verdict analyze()/classify() already
produces --- not better, not stricter, the SAME. This suite pins that: on the
golden anchors the code must reproduce today's regime verdict, and in general the
code must be a pure relabel of (A1 flag, geometric regime, sigma_c defined?).

Per the discipline: the question is whether the CLAIM
shifted or only its FORM. Any deviation of the code from today's verdict is a
finding about classify(), to be recorded in THEOREM_MAP BEFORE it is coded ---
not a silent improvement. This test runs before any *new* case is added.
"""
import math

import numpy as np

from sigma_c import analyze, bare, two_probe_test, Framework
from sigma_c.codes import Code


# --- golden anchors, reconstructed exactly as in test_golden_paper_anchors.py ---

def _ising_bare_correlator():
    """Regime I anchor: 1D Ising bare correlator (paper Example C.7)."""
    beta_j = 0.5
    sigma = np.geomspace(0.02, 100.0, 1000)
    C = math.tanh(beta_j) ** sigma
    return sigma, C


def _power_law():
    """Regime III anchor: pure power law, chi_O monotone (paper Thm 9.1)."""
    r_arr = np.geomspace(1e-3, 1.0, 400)
    return r_arr, r_arr ** 0.6


# ---------------------------------------------------------------------------
# 1. The code reproduces today's verdict on the golden anchors
# ---------------------------------------------------------------------------

def test_regime_I_anchor_sigma_c_code_is_OK():
    sigma, C = _ising_bare_correlator()
    r = analyze(sigma, C, window=bare(), framework=Framework.TRANSFER_MATRIX_1D)
    # today's verdict (unchanged):
    assert r.regime.geometric == "I_geom"
    assert r.sigma_c is not None
    # the code is its formalisation:
    assert r.sigma_c_status.code is Code.OK, str(r.sigma_c_status)


def test_regime_III_anchor_sigma_c_code_is_NOT_RESOLVABLE():
    r_arr, W = _power_law()
    r = analyze(r_arr, W, window=bare())
    # today's verdict (unchanged):
    assert r.regime.geometric == "III_geom"
    assert r.sigma_c is None
    # the code is its formalisation:
    assert r.sigma_c_status.code is Code.NOT_RESOLVABLE, str(r.sigma_c_status)


def test_A1_violation_is_NOT_APPLICABLE_not_a_silent_value():
    # Same regime-I anchor, but the caller declares preprocessing NOT
    # scale-equivariant. Today: .falsifiable == False, "exploratory". The code
    # formalises that as NOT_APPLICABLE (a precondition of the method is violated).
    sigma, C = _ising_bare_correlator()
    r = analyze(sigma, C, window=bare(), preprocessing_scale_equivariant=False)
    assert r.falsifiable is False           # today's verdict
    assert r.sigma_c_status.code is Code.NOT_APPLICABLE, str(r.sigma_c_status)


# ---------------------------------------------------------------------------
# 2. The code is a PURE RELABEL of the existing verdict (no new threshold)
# ---------------------------------------------------------------------------

def _all_anchor_results():
    sigma, C = _ising_bare_correlator()
    r_arr, W = _power_law()
    return [
        analyze(sigma, C, window=bare(), framework=Framework.TRANSFER_MATRIX_1D),
        analyze(r_arr, W, window=bare()),
        analyze(sigma, C, window=bare(), preprocessing_scale_equivariant=False),
        analyze(sigma, C, window=bare(), preprocessing_scale_equivariant=True),
    ]


def test_sigma_c_code_is_a_pure_function_of_the_existing_verdict():
    """The biconditionals that make the code a relabel, not a new decision:

      NOT_APPLICABLE  <=>  A1 declared violated (preprocessing_scale_equivariant is False)
      NOT_RESOLVABLE  <=>  (A1 ok) and geometric regime III
      OK              <=>  (A1 ok) and regime I/II with sigma_c defined
    """
    for r in _all_anchor_results():
        code = r.sigma_c_status.code
        a1_violated = (r.preprocessing_scale_equivariant is False)
        geom = r.regime.geometric

        if a1_violated:
            assert code is Code.NOT_APPLICABLE, str(r.sigma_c_status)
            continue
        if geom == "III_geom":
            assert code is Code.NOT_RESOLVABLE, str(r.sigma_c_status)
            assert r.sigma_c is None
        else:  # I_geom or II_geom
            assert code is Code.OK, str(r.sigma_c_status)
            assert r.sigma_c is not None
        # sigma_c NEVER carries NOT_IDENTIFIED (that gate is tau's):
        assert code is not Code.NOT_IDENTIFIED


def test_ok_verdict_carries_a_reason_and_no_remediation():
    sigma, C = _ising_bare_correlator()
    v = analyze(sigma, C, window=bare()).sigma_c_status
    assert v.reason  # non-empty human reason
    assert v.remediation == ""  # remediation is a NOT_IDENTIFIED-only field


# ---------------------------------------------------------------------------
# 3. tau: the baseline is NOT today's (broken) unconditional tau. It is
#    NOT_IDENTIFIED, lifted only by a proven condition. The guard has the
#    OPPOSITE sign to sigma_c: an unexpected OK on a REAL anchor is the alarm.
# ---------------------------------------------------------------------------

def _synthetic_gapped_chain(gamma_A=0.63, sigma_max=60.0):
    """The ONE anchor allowed to reach tau OK: the live-route 3-state chain,
    spectrum {0,-1,-3} known BY CONSTRUCTION (mu_2=-1 isolated, gap Delta=2).
    Autocorrelation observable + the transfer eigenvalues as a supplied spectrum
    + a reversible non-experimental framework + the contamination-to-signal gamma_A
    (known by construction) so the tail-window T>t* can be MACHINE-CHECKED.
    gamma_A default 0.63 (<1 -> t*=0); sigma_max sets the observation window T."""
    sigma = np.geomspace(0.02, sigma_max, 1000)
    O = 1.0 * np.exp(-sigma / 1.0) + 0.4 * np.exp(-sigma / (1.0 / 3.0))
    lam = [1.0, math.exp(-1.0), math.exp(-3.0)]  # |lambda| decreasing, isolated gap
    # self-adjoint route: an OPERATOR (self-adjoint diag) is VERIFIED, so it can
    # reach window_readability OK; a bare spectrum would only be 'declared' (no probe
    # -> no Gamma_A -> window can never be OK). evolution_time gives T=max(sigma-grid).
    A = np.diag(lam)
    return analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                   operator=A, inner_product="euclidean", gamma_A=gamma_A,
                   T_star=1.0, sigma_axis="evolution_time")


def test_real_anchor_emitting_tau_today_is_NOT_IDENTIFIED():
    # The Ising bare correlator emits result.tau_bridge (the fallen load-dominant
    # read-out). The contract refuses to certify a relaxation time: tau_abscissa is
    # NOT_IDENTIFIED (no spectrum/operator), with a remediation naming what to supply.
    sigma, C = _ising_bare_correlator()
    r = analyze(sigma, C, window=bare(), framework=Framework.TRANSFER_MATRIX_1D)
    assert r.tau_bridge is not None, "the bridge read-out is emitted here"
    assert r.tau_abscissa is None, "no spectrum/operator -> no certified abscissa"
    v = r.tau_abscissa_status
    assert v.code is Code.NOT_IDENTIFIED, str(v)
    assert v.remediation, "NOT_IDENTIFIED must name what to supply"
    assert "two_probe" in v.remediation and "spectrum" in v.remediation


def test_regime_III_tau_abscissa_is_NOT_IDENTIFIED():
    # sigma_c-only regime III: no spectrum/operator -> tau_abscissa NOT_IDENTIFIED.
    r_arr, W = _power_law()
    r = analyze(r_arr, W, window=bare())
    assert r.tau_abscissa is None
    assert r.tau_abscissa_status.code is Code.NOT_IDENTIFIED, str(r.tau_abscissa_status)


def test_A1_violation_tau_is_NOT_APPLICABLE():
    sigma, C = _ising_bare_correlator()
    r = analyze(sigma, C, window=bare(), preprocessing_scale_equivariant=False)
    assert r.tau_abscissa_status.code is Code.NOT_APPLICABLE, str(r.tau_abscissa_status)


def test_only_the_synthetic_gapped_anchor_reaches_window_OK():
    # The ONE window-readability OK: a self-adjoint OPERATOR + isolated gap + gamma_A,
    # with the tail-window T > max(t*, 1/|mu_3|) MACHINE-CHECKED.
    r = _synthetic_gapped_chain()
    assert r.tau_abscissa is not None
    assert r.regime.spectral == "gap", r.regime.spectral
    assert r.tau_abscissa_status.code is Code.OK, str(r.tau_abscissa_status)
    v = r.window_readability_status
    assert v.code is Code.OK, str(v)
    assert "T=" in v.reason and "t*=" in v.reason, "OK must show the machine-checked gate"


def test_bare_spectrum_gap_is_declared_window_NOT_IDENTIFIED():
    # A bare SPECTRUM (no operator) is only 'declared': self-adjointness + faithfulness
    # unverified -> window readability can never be OK (no probe -> no Gamma_A).
    sigma = np.geomspace(0.02, 60.0, 1000)
    O = 1.0 * np.exp(-sigma) + 0.4 * np.exp(-3.0 * sigma)
    lam = [1.0, math.exp(-1.0), math.exp(-3.0)]
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                spectrum=lam, T_star=1.0, sigma_axis="evolution_time")  # no gamma_A, spectrum-only
    assert r.selfadjoint == "declared"
    assert r.window_readability_status.code is Code.NOT_IDENTIFIED, str(r.window_readability_status)


def test_verified_gap_but_window_too_short_is_NOT_IDENTIFIED():
    # gamma_A large -> t* = (2/Delta) log gamma_A large; a short window T < t* fails
    # the feasibility gate even with a verified operator gap. Delta=2, gamma_A=100 -> t*~4.6.
    r = _synthetic_gapped_chain(gamma_A=100.0, sigma_max=2.0)  # T=2 < t*~4.6
    assert r.regime.spectral == "gap"
    v = r.window_readability_status
    assert v.code is Code.NOT_IDENTIFIED, str(v)
    assert "too short" in v.reason, str(v)


def test_no_REAL_anchor_reaches_tau_OK_the_guard():
    # THE GUARD (opposite sign to sigma_c): any real anchor (no supplied spectrum/
    # operator) that reaches OK is the alarm, meaning a condition is too lax. All real
    # anchors must be NOT_IDENTIFIED / NOT_RESOLVABLE / NOT_APPLICABLE -- never OK, on
    # BOTH tau fields.
    sigma, C = _ising_bare_correlator()
    r_arr, W = _power_law()
    real_anchors = [
        analyze(sigma, C, window=bare(), framework=Framework.TRANSFER_MATRIX_1D),
        analyze(sigma, C, window=bare()),                       # no framework at all
        analyze(r_arr, W, window=bare()),                       # regime III
        analyze(sigma, C, window=bare(), preprocessing_scale_equivariant=False),
    ]
    for r in real_anchors:
        assert r.tau_abscissa_status.code is not Code.OK, (
            f"real anchor unexpectedly abscissa-OK -- a condition is too lax: {r.tau_abscissa_status}"
        )
        assert r.window_readability_status.code is not Code.OK, (
            f"real anchor unexpectedly window-OK: {r.window_readability_status}"
        )


# ---------------------------------------------------------------------------
# 4. Two-probe tau verdict: a PASS is OK-with-one-sided-caveat, NEVER bare OK.
# ---------------------------------------------------------------------------

def test_two_probe_pass_is_OK_with_one_sided_caveat():
    # A passing two-probe result from two ANALYTIC probes. (Constructed directly:
    # a clean pass from the registry windows is hard to build --- single-peaked
    # windows all have rho_star=1 so sigma_c collapses [probes_not_distinct], and
    # the distinct-rho_star windows are regime II. That is why the golden test
    # computes tau by hand. Here we unit-test the verdict logic on the outcome.)
    from sigma_c.result import TwoProbeResult
    tp = TwoProbeResult(
        passed=True, delta=0.01, delta_threshold=0.20, cause=None,
        tau_1=1.30, tau_2=1.30,
        rho_star_source_1="analytic:bare", rho_star_source_2="analytic:exponential",
    )
    v = tp.tau_status
    assert v.code is Code.OK, str(v)
    # ...but OK WITH the one-sided caveat spelled out --- never a bare OK:
    assert "ONE-SIDED" in v.reason and "does NOT certify" in v.reason, str(v)


def test_two_probe_fit_based_pass_cannot_certify():
    # A pass where a probe is fit-based (not analytic) cannot certify tau.
    from sigma_c.result import TwoProbeResult
    tp = TwoProbeResult(
        passed=True, delta=0.01, delta_threshold=0.20, cause=None,
        tau_1=1.30, tau_2=1.30,
        rho_star_source_1="analytic:bare", rho_star_source_2="fitted",
    )
    assert tp.tau_status.code is Code.NOT_IDENTIFIED, str(tp.tau_status)


def test_two_probe_deviation_is_NOT_IDENTIFIED():
    # Two observables with genuinely different tau -> deviation -> trouble.
    sigma = np.geomspace(0.02, 100.0, 1000)
    C_fast = math.tanh(0.5) ** sigma          # tau ~ 1.30
    C_slow = math.tanh(0.75) ** sigma         # different tau
    r1 = analyze(sigma, C_fast, window=bare())
    r2 = analyze(sigma, C_slow, window=bare())
    tp = two_probe_test(r1, r2)
    assert not tp.passed, "probes with different tau should deviate"
    assert tp.tau_status.code is Code.NOT_IDENTIFIED, str(tp.tau_status)


# ---------------------------------------------------------------------------
# 5. The tail-onset numerics themselves (the only genuinely computed condition).
# ---------------------------------------------------------------------------

def test_tail_onset_t_star_formula():
    from sigma_c.core.faithfulness import tail_onset_t_star, tail_window_feasible
    # Gamma_A <= 1 (well-loaded slow mode) -> t* = 0 (tail dominates from the start)
    assert tail_onset_t_star(0.5, 2.0) == 0.0
    assert tail_onset_t_star(1.0, 2.0) == 0.0
    # Gamma_A > 1 -> t* = (2/Delta) log Gamma_A
    assert tail_onset_t_star(math.e, 2.0) == pytest_approx(1.0)   # (2/2)*ln(e)=1
    assert tail_onset_t_star(100.0, 2.0) == pytest_approx(math.log(100.0))  # (2/2)*ln100
    # no resolved faster mode (Delta inf / <=0) -> no contamination onset
    assert tail_onset_t_star(100.0, float("inf")) == 0.0
    assert tail_onset_t_star(100.0, 0.0) == 0.0
    # feasibility gate is strict T > t*
    assert tail_window_feasible(5.0, math.e, 2.0) is True    # 5 > 1
    assert tail_window_feasible(0.5, math.e, 2.0) is False   # 0.5 < 1


def pytest_approx(x):
    import pytest
    return pytest.approx(x, rel=1e-9, abs=1e-12)

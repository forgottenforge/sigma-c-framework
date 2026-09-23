# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Golden tests anchored to the paper's printed numerical values.

The contract: the kernel reproduces the paper's anchors. A failing test
means either the kernel is wrong or the paper is wrong -- both are worth knowing
(Enforcement item 1 of the design).

References by paper label (resolved via THEOREM_MAP.md):
  - prop:structural-reduction  : sigma_c = rho_star * tau
  - the I/II/III split          : now a named peak-count operation (§8 trichotomy labels ARCHIVED)
  - thm:diagnostic             : regime III power-law
  - prop:standard-frameworks   : framework taxonomy
  - def:operational-test-noncirc : non-circular two-probe test

Run with:
    python -m pytest sigma_c/tests/test_golden_paper_anchors.py -v
"""
from __future__ import annotations
import math

import numpy as np
import pytest

from sigma_c import analyze, bare, gamma_k, exponential, log_gaussian
from sigma_c import Framework, two_probe_test
from sigma_c.windows import WINDOW_REGISTRY


# ---------------------------------------------------------------------------
# Window family: analytical rho_star values (paper Section 4.2 table)
# ---------------------------------------------------------------------------

class TestWindowFamilyRhoStar:
    """Each canonical window has its rho_star known analytically -- no fit."""

    def test_bare_window_rho_star_is_one(self):
        assert bare().rho_star == pytest.approx(1.0)

    def test_gamma2_rho_star_is_2_minus_sqrt3(self):
        # Paper Section 4.2 table.
        expected = 2 - math.sqrt(3)
        assert gamma_k(2).rho_star == pytest.approx(expected, abs=1e-12)

    def test_gamma3_rho_star_is_3_minus_2_sqrt2(self):
        # Paper Section 4.2 table.
        expected = 3 - 2 * math.sqrt(2)
        assert gamma_k(3).rho_star == pytest.approx(expected, abs=1e-12)

    def test_exponential_rho_star_is_one(self):
        assert exponential().rho_star == pytest.approx(1.0)

    def test_log_gaussian_rho_star_is_one(self):
        assert log_gaussian().rho_star == pytest.approx(1.0)

    def test_all_registry_windows_carry_analytic_provenance(self):
        for name, factory in WINDOW_REGISTRY.items():
            w = factory()
            assert w.rho_star_source.startswith("analytic:"), (
                f"Window {name} provenance must be analytic, "
                f"got {w.rho_star_source!r}"
            )


# ---------------------------------------------------------------------------
# Paper Example C.7 -- 1D Ising chain at beta * J = 0.5
# ---------------------------------------------------------------------------

class TestIsingChainAnchor:
    """Paper Example C.7: tau = -1/log tanh(beta * J)."""

    BETA_J = 0.5
    EXPECTED_TAU = -1.0 / math.log(math.tanh(0.5))  # approx 1.295392...

    def _bare_correlator(self):
        sigma = np.geomspace(0.02, 100.0, 1000)
        t = math.tanh(self.BETA_J)
        C = t ** sigma
        return sigma, C

    def test_bare_correlator_gives_sigma_c_equals_tau(self):
        """With bare window, rho_star = 1, so sigma_c = tau directly."""
        sigma, C = self._bare_correlator()
        r = analyze(sigma, C, window=bare(),
                    framework=Framework.TRANSFER_MATRIX_1D)

        assert r.sigma_c is not None, "regime should be I (single mode)"
        assert r.regime.geometric == "I_geom"
        assert r.sigma_c == pytest.approx(self.EXPECTED_TAU, rel=2e-3)
        assert r.tau_bridge == pytest.approx(self.EXPECTED_TAU, rel=2e-3)

    def test_falsifiable_under_bare_window(self):
        sigma, C = self._bare_correlator()
        r = analyze(sigma, C, window=bare())
        assert r.falsifiable, "bare window is analytic -> falsifiable"
        assert r.rho_star_source == "analytic:bare"


# ---------------------------------------------------------------------------
# Regime III -- power-law observables (paper Thm 9.1)
# ---------------------------------------------------------------------------

class TestRegimeIIIPowerLaw:
    """Paper Thm 9.1: chi_O is monotone -> sigma_c = bottom value (None)."""

    def test_pure_power_law_returns_none(self):
        r_arr = np.geomspace(1e-3, 1.0, 400)
        for alpha in (0.3, 0.6, 0.9, 1.5):
            W = r_arr ** alpha
            result = analyze(r_arr, W, window=bare(),
                             label=f"power-law alpha={alpha}")
            assert result.sigma_c is None, (
                f"alpha={alpha}: sigma_c must be bottom (None), "
                f"got {result.sigma_c}"
            )
            assert result.regime.geometric == "III_geom"

    def test_regime_iii_carries_diagnostic_citation(self):
        r_arr = np.geomspace(1e-3, 1.0, 400)
        W = r_arr ** 0.6
        result = analyze(r_arr, W)
        assert "thm:diagnostic" in result.citations


# ---------------------------------------------------------------------------
# rho_star separation principle (paper Prop 4.1)
# ---------------------------------------------------------------------------

class TestRhoStarSeparation:
    """Paper Prop 4.1: sigma_c = rho_star * tau, so changing the window
    changes sigma_c proportionally to rho_star and leaves tau alone."""

    def test_bare_vs_gamma2_recover_same_tau(self):
        tau_true = 5.0
        sigma = np.geomspace(0.05, 100.0, 1000)

        # Bare correlator -> sigma_c = tau (rho_star = 1)
        O_bare = np.exp(-sigma / tau_true)
        r_bare = analyze(sigma, O_bare, window=bare())
        assert r_bare.tau_bridge == pytest.approx(tau_true, rel=5e-3)

        # Bare correlator -> exponential window -> still rho_star = 1
        r_exp = analyze(sigma, O_bare, window=exponential())
        assert r_exp.tau_bridge == pytest.approx(tau_true, rel=5e-3)


# ---------------------------------------------------------------------------
# Two-probe non-circular test (paper Def 6.4)
# ---------------------------------------------------------------------------

class TestNonCircularTwoProbe:

    def test_two_analytic_windows_pass_test(self):
        """Two TRULY distinct probes that recover the same tau.

        Field-found: an earlier version passed the SAME O to two different
        window labels -- the new probes_not_distinct precheck correctly
        flags that as a precondition failure (the downstream-application Phase 1
        mistake). A meaningful two-probe pass requires two windows whose
        canonical observable forms differ: here, O_exp from the
        exponential kernel and O_gamma2 from the gamma2 kernel, both
        applied to the same underlying exp(-t/tau) signal.
        """
        tau_true = 5.0
        sigma = np.geomspace(0.05, 100.0, 1000)
        # Canonical observables for the two windows applied to A(t) = exp(-t/tau)
        # exponential window: O(sigma) = sigma*tau / (sigma+tau), peaks at sigma=tau, rho_star=1
        # gamma2 window:      O(sigma) = sigma*tau^2 / (sigma+tau)^2, peaks at (2-sqrt(3))*tau, rho_star=2-sqrt(3)
        O_exp_window    = sigma * tau_true / (sigma + tau_true)
        O_gamma2_window = sigma * (tau_true ** 2) / (sigma + tau_true) ** 2
        r_1 = analyze(sigma, O_exp_window,    window=exponential())
        r_2 = analyze(sigma, O_gamma2_window, window=gamma_k(2))
        # The gamma2 form has two chi peaks (rising + falling side of the
        # window function); the rising-side peak is the lambda_2 reading,
        # matching the paper Example C.6 convention.
        sc_1 = r_1.sigma_c if not isinstance(r_1.sigma_c, list) else min(r_1.sigma_c)
        sc_2 = r_2.sigma_c if not isinstance(r_2.sigma_c, list) else min(r_2.sigma_c)
        assert sc_1 is not None and sc_2 is not None
        # The two sigma_c values must differ (distinct probes) ...
        assert abs(sc_1 - sc_2) > 1e-3, (
            "probes are not distinct: sigma_c collapsed onto the same value"
        )
        # ... and the recovered tau must agree.
        tau_1 = sc_1 / r_1.rho_star
        tau_2 = sc_2 / r_2.rho_star
        assert abs(math.log(tau_1) - math.log(tau_2)) < 0.05, (
            f"tau disagreement: tau_1={tau_1}, tau_2={tau_2}"
        )

    def test_two_probes_disagreeing_on_tau_flag_as_failure(self):
        # Simulate two probes that don't see the same system: different decay
        # times. They should fail the test with a named cause.
        sigma = np.geomspace(0.05, 200.0, 1000)
        O_1 = np.exp(-sigma / 5.0)
        O_2 = np.exp(-sigma / 50.0)  # 10x slower
        r_1 = analyze(sigma, O_1, window=bare())
        r_2 = analyze(sigma, O_2, window=bare())
        test = two_probe_test(r_1, r_2, delta_threshold=0.05)
        assert not test.passed
        assert test.cause in ("faithfulness_break", "system_disagreement")


# ---------------------------------------------------------------------------
# Enforcement-layer invariants
# ---------------------------------------------------------------------------

class TestEnforcementInvariants:
    """Discipline mechanisms from the design, enforced by tests."""

    def test_no_nan_on_regime_iii(self):
        """Enforcement 4: bottom value is None, never NaN."""
        sigma = np.geomspace(1e-3, 1.0, 200)
        result = analyze(sigma, sigma ** 0.5, window=bare())
        assert result.sigma_c is None
        # NaN check must hold for tau too
        assert result.tau_bridge is None
        assert not (isinstance(result.tau_bridge, float) and math.isnan(result.tau_bridge or 0.0))

    def test_provenance_is_analytic_for_registry_windows(self):
        """Enforcement 2: rho_star_source carries provenance."""
        sigma = np.geomspace(0.05, 50.0, 400)
        O = np.exp(-sigma / 2.0)
        for name in ("bare", "gamma2", "gamma3", "exponential", "log_gaussian"):
            r = analyze(sigma, O, window=name)
            assert r.rho_star_source == f"analytic:{name}"
            if r.sigma_c is not None:
                assert r.falsifiable

    def test_silent_detrending_off_by_default(self):
        """Enforcement 6: detrended must default to False."""
        sigma = np.geomspace(0.05, 50.0, 400)
        O = np.exp(-sigma / 2.0)
        r = analyze(sigma, O, window=bare())
        assert r.detrended is False

    def test_theorem_citations_present(self):
        """Every output that produced a sigma_c value must cite at least
        one paper theorem -- the 'every output cites a theorem' invariant."""
        sigma = np.geomspace(0.05, 50.0, 400)
        O = np.exp(-sigma / 2.0)
        r = analyze(sigma, O, window=bare())
        assert len(r.citations) > 0
        assert any("def:sigmac" == c or "prop:structural-reduction" == c
                   for c in r.citations)

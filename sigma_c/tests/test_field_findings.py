# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Tests for the three field-found patches.

Provenance: each of these was a real trap that a downstream-application Phase 1 walked
into; full narrative in `docs/V4_FRAMEWORK_FEEDBACK.md` of the downstream application
repo (and now also in this repo's CHANGELOG). The tests below pin the
fixes so any future refactor cannot quietly undo them.
"""
from __future__ import annotations

import math

import numpy as np
import pytest

from sigma_c import analyze, bare, exponential, gamma_k, two_probe_test
from sigma_c.core.susceptibility import quadratic_peak_in_log_sigma


# ---------------------------------------------------------------------------
# Finding 1 -- sub-grid peak localization
# ---------------------------------------------------------------------------

class TestSubGridPeakLocalization:
    """sigma_c is reported with sub-grid quadratic refinement; cross-tests
    on the same grid no longer see bit-identical sigma_c by coincidence."""

    def test_sub_grid_offset_is_in_minus_one_one(self):
        """Quadratic interpolation offset must clamp to [-1, 1] times the
        grid step; runaway sub-grid values are not allowed."""
        sigma = np.geomspace(0.1, 10.0, 32)
        # Construct chi with a peak between two grid samples
        i_peak = 16
        chi = np.zeros_like(sigma)
        # Skew: y_{-1} = 0.7, y_0 = 1.0, y_+1 = 0.9 -> true peak between i_peak and i_peak+1
        chi[i_peak - 1] = 0.7
        chi[i_peak] = 1.0
        chi[i_peak + 1] = 0.9
        sc_sub, at_boundary = quadratic_peak_in_log_sigma(sigma, chi, i_peak)
        assert not at_boundary
        # sub-grid sigma_c must be between sigma_grid[i_peak] and sigma_grid[i_peak+1]
        assert sigma[i_peak] <= sc_sub <= sigma[i_peak + 1]

    def test_boundary_flag_set_when_peak_at_grid_edge(self):
        """When chi peaks at index 0 or len-1, no sub-grid refinement
        is possible and the at_grid_boundary flag is set."""
        sigma = np.geomspace(0.1, 10.0, 32)
        chi = np.zeros_like(sigma)
        chi[0] = 1.0
        chi[1:] = 0.1
        sc_sub, at_boundary = quadratic_peak_in_log_sigma(sigma, chi, 0)
        assert at_boundary
        assert sc_sub == pytest.approx(sigma[0])

    def test_result_carries_at_grid_boundary_flag(self):
        """analyze() surfaces the boundary flag on the Result."""
        sigma = np.geomspace(0.1, 10.0, 50)
        # Observable whose chi peak lies INSIDE the grid -- no boundary
        O = sigma * 1.0 / (sigma + 1.0)
        r = analyze(sigma, O, window=exponential())
        assert r.sigma_c_at_grid_boundary is False

    def test_two_separate_analyses_no_longer_collide_on_grid(self):
        """Two analyses with chi peaks shifted by < 1 grid cell now report
        distinct sub-grid sigma_c values."""
        sigma = np.geomspace(0.1, 10.0, 32)
        # Two slightly-shifted analytical observables, both with peaks
        # inside the grid range but at slightly different positions.
        tau_1 = 1.0
        tau_2 = 1.05
        O_1 = sigma * tau_1 / (sigma + tau_1)
        O_2 = sigma * tau_2 / (sigma + tau_2)
        r_1 = analyze(sigma, O_1, window=exponential())
        r_2 = analyze(sigma, O_2, window=exponential())
        # Before sub-grid refinement, both would have snapped to the nearest grid sample
        # and matched. With sub-grid, the small difference is visible.
        assert r_1.sigma_c != r_2.sigma_c


# ---------------------------------------------------------------------------
# Finding 2 -- probes_not_distinct precheck on two_probe_test
# ---------------------------------------------------------------------------

class TestProbesNotDistinct:
    """A two-probe test where both probes collapse onto the same sigma_c
    is a precondition failure, not a system property. The kernel flags this
    instead of silently passing (or failing under another cause)."""

    def test_same_observable_two_windows_flagged_as_not_distinct(self):
        """The downstream-application Phase 1 mistake, in one assertion: pass the same
        observable through two windows. The kernel must flag this as
        probes_not_distinct, never as faithfulness_break."""
        tau_true = 5.0
        sigma = np.geomspace(0.05, 100.0, 1000)
        # ONE observable, used for both windows -- this is the trap.
        O = np.exp(-sigma / tau_true)
        r_1 = analyze(sigma, O, window=bare())
        r_2 = analyze(sigma, O, window=exponential())
        test = two_probe_test(r_1, r_2, delta_threshold=0.05)
        # sigma_c was the same (single chi peak shared between the two
        # Result objects, sub-grid refinement still puts both at the
        # same point because chi is identical) -- the precheck must fire.
        assert not test.passed
        assert test.cause == "probes_not_distinct"

    def test_truly_distinct_probes_pass_normally(self):
        """Sanity check: two windows applied to their OWN canonical
        observables produce distinct sigma_c values, and the test passes."""
        tau_true = 5.0
        sigma = np.geomspace(0.05, 100.0, 1000)
        # Different observables per window -- this is the correct usage.
        O_exp    = sigma * tau_true / (sigma + tau_true)
        O_gamma2 = sigma * (tau_true ** 2) / (sigma + tau_true) ** 2
        r_1 = analyze(sigma, O_exp,    window=exponential())
        r_2 = analyze(sigma, O_gamma2, window=gamma_k(2))
        # gamma2 may give regime II; collapse to rising side if so.
        sc_2 = r_2.sigma_c if not isinstance(r_2.sigma_c, list) else min(r_2.sigma_c)
        # The sigma_c values must differ by more than the grid step
        log_sigma_step = math.log(sigma[1]) - math.log(sigma[0])
        assert abs(math.log(r_1.sigma_c) - math.log(sc_2)) > log_sigma_step

        # With distinct probes, two-probe runs the actual delta check.
        test = two_probe_test(r_1, r_2, delta_threshold=0.05)
        assert test.cause != "probes_not_distinct"


# ---------------------------------------------------------------------------
# Finding 3 -- preprocessing_scale_equivariant flag as A1 guard
# ---------------------------------------------------------------------------

class TestA1PreprocessingGuard:
    """When the user declares the preprocessing as not scale-equivariant,
    the Result is exploratory (not falsifiable), even if the window is
    one of the analytic canonical five."""

    def test_falsifiable_unchanged_when_not_declared(self):
        """Default behaviour (preprocessing_scale_equivariant=None) must
        be identical to the default to preserve backward compatibility."""
        sigma = np.geomspace(0.1, 10.0, 200)
        O = sigma * 1.0 / (sigma + 1.0)
        r = analyze(sigma, O, window=exponential())  # no A1 declaration
        assert r.preprocessing_scale_equivariant is None
        assert r.falsifiable  # backward compat: analytic rho_star still falsifiable

    def test_falsifiable_unchanged_when_declared_true(self):
        """Explicit True declaration is the same as the default plus an
        affirmation note."""
        sigma = np.geomspace(0.1, 10.0, 200)
        O = sigma * 1.0 / (sigma + 1.0)
        r = analyze(sigma, O, window=exponential(),
                    preprocessing_scale_equivariant=True)
        assert r.preprocessing_scale_equivariant is True
        assert r.falsifiable
        # rho_star_source remains analytic
        assert r.rho_star_source.startswith("analytic:")

    def test_declared_false_downgrades_to_exploratory(self):
        """Explicit False declaration: even with an analytic window, the
        Result is exploratory and not falsifiable."""
        sigma = np.geomspace(0.1, 10.0, 200)
        O = sigma * 1.0 / (sigma + 1.0)
        r = analyze(sigma, O, window=exponential(),
                    preprocessing_scale_equivariant=False)
        assert r.preprocessing_scale_equivariant is False
        assert not r.falsifiable
        assert "exploratory:" in r.rho_star_source
        # The note must explain why
        joined_notes = " ".join(r.notes)
        assert "A1" in joined_notes or "absolute" in joined_notes.lower()

    def test_a1_violation_propagates_through_regime_iii(self):
        """If chi has no peak (regime III), the A1 flag should still
        propagate to the Result for provenance."""
        sigma = np.geomspace(0.1, 10.0, 200)
        # Power-law observable: no interior peak -> regime III
        O = sigma ** 0.5
        r = analyze(sigma, O, window=exponential(),
                    preprocessing_scale_equivariant=False)
        assert r.preprocessing_scale_equivariant is False
        assert r.sigma_c is None  # regime III

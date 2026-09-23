# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Tests for the F2/F3 faithfulness sufficient conditions (paper Prop 5.10).

These verdict-emitters were public-API exports with zero test coverage — i.e.
checkers nobody knew would run. This file exercises both across a pass case
and their failure/guard cases, so the exported checks are known to execute and
return sensible verdicts. (This tests that the CHECKERS run correctly; it does
NOT certify prop:faith-sufficient, whose register status cite() still renders.)
"""
from __future__ import annotations
import math

from sigma_c.core.faithfulness import check_F2, check_F3, FaithfulnessCheck


# ---------------------------------------------------------------------------
# F2 — spectrally filtered window: band must contain u_2 and exclude u_3
# ---------------------------------------------------------------------------

class TestCheckF2:
    # spectrum 1.0 / 0.5 / 0.1  ->  u_2 = -log(0.5) = 0.6931, u_3 = -log(0.1) = 2.3026
    SPECTRUM = [1.0, 0.5, 0.1]

    def test_passes_when_band_contains_u2_and_excludes_u3(self):
        r = check_F2(spectrum=self.SPECTRUM, window_mellin_support=(0.5, 1.0))
        assert isinstance(r, FaithfulnessCheck)
        assert r.condition == "F2"
        assert r.passed is True
        assert r.details["contains_u2"] and r.details["excludes_u3"]
        assert math.isfinite(r.C_R)            # achieved filter tolerance
        assert r.faithfulness_order == math.inf

    def test_fails_when_band_also_admits_u3(self):
        # b = 2.5 > u_3 = 2.303  ->  the sub-leading mode is not filtered out
        r = check_F2(spectrum=self.SPECTRUM, window_mellin_support=(0.5, 2.5))
        assert r.condition == "F2"
        assert r.passed is False
        assert r.details["contains_u2"] is True
        assert r.details["excludes_u3"] is False
        assert math.isinf(r.C_R)

    def test_fails_when_band_misses_u2(self):
        # band entirely above u_2 = 0.693
        r = check_F2(spectrum=self.SPECTRUM, window_mellin_support=(1.0, 2.0))
        assert r.passed is False
        assert r.details["contains_u2"] is False

    def test_refuses_with_fewer_than_three_eigenvalues(self):
        r = check_F2(spectrum=[1.0, 0.5], window_mellin_support=(0.5, 1.0))
        assert r.passed is False
        assert "reason" in r.details


# ---------------------------------------------------------------------------
# F3 — Doeblin-minorised system with bounded-variation probe
# ---------------------------------------------------------------------------

class TestCheckF3:
    def test_passes_for_valid_doeblin_probe(self):
        r = check_F3(doeblin_epsilon=0.5, n_0=1, probe_tv=1.0, abs_c_2=0.5)
        assert isinstance(r, FaithfulnessCheck)
        assert r.condition == "F3"
        assert r.passed is True
        # decay = (1-0.5)^1 = 0.5 ; series = 0.5/0.5 = 1.0 ; C_R = (1.0/0.5)*1.0 = 2.0
        assert r.details["decay_factor"] == 0.5
        assert r.C_R == 2.0
        assert r.faithfulness_order >= 1

    def test_refuses_epsilon_out_of_range(self):
        for bad in (0.0, 1.0, 1.5, -0.1):
            r = check_F3(doeblin_epsilon=bad, n_0=1, probe_tv=1.0, abs_c_2=0.5)
            assert r.passed is False
            assert math.isinf(r.C_R)

    def test_refuses_missing_or_nonpositive_c2(self):
        r_none = check_F3(doeblin_epsilon=0.5, n_0=1, probe_tv=1.0, abs_c_2=None)
        r_zero = check_F3(doeblin_epsilon=0.5, n_0=1, probe_tv=1.0, abs_c_2=0.0)
        assert r_none.passed is False and r_zero.passed is False

    def test_refuses_bad_n0(self):
        r = check_F3(doeblin_epsilon=0.5, n_0=0, probe_tv=1.0, abs_c_2=0.5)
        assert r.passed is False

    def test_smaller_epsilon_gives_larger_C_R(self):
        # weaker minorisation (smaller eps) -> slower decay -> larger constant
        strong = check_F3(doeblin_epsilon=0.7, n_0=1, probe_tv=1.0, abs_c_2=0.5)
        weak = check_F3(doeblin_epsilon=0.3, n_0=1, probe_tv=1.0, abs_c_2=0.5)
        assert strong.passed and weak.passed
        assert weak.C_R > strong.C_R

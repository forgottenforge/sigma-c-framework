# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Non-finite input must be REJECTED, not silently turned into a
spurious OK. NaN/Inf in sigma or O previously propagated through chi_O's gradient
and produced a clean-looking sigma_c with status OK and no warning — a false-green,
the same class of defect as a mis-anchored citation. Bad input is not a
measurement, so it is a ValueError (not one of the four applicability codes).
"""
from __future__ import annotations
import numpy as np
import pytest

from sigma_c import analyze, bare


class TestNonFiniteInputRejected:
    def _sig(self):
        s = np.geomspace(0.05, 50.0, 400)
        return s, np.exp(-s / 2.0)

    def test_nan_in_O_raises(self):
        s, O = self._sig()
        O[10] = np.nan
        with pytest.raises(ValueError):
            analyze(s, O, window=bare())

    def test_inf_in_O_raises(self):
        s, O = self._sig()
        O[10] = np.inf
        with pytest.raises(ValueError):
            analyze(s, O, window=bare())

    def test_nan_in_sigma_raises(self):
        s, O = self._sig()
        s[10] = np.nan
        with pytest.raises(ValueError):
            analyze(s, O, window=bare())

    def test_nan_from_callable_O_raises(self):
        # callable-O path (api.py chi_O_from_callable) -- bypassed chi_O's guard
        s, _ = self._sig()
        with pytest.raises(ValueError):
            analyze(s, lambda x: np.where(x > 5.0, np.nan, np.exp(-x / 2.0)),
                    window=bare())

    def test_nan_in_precomputed_chi_raises(self):
        # precomputed-chi path (api.py) -- never calls chi_O at all
        s, O = self._sig()
        chi = np.abs(np.gradient(O, np.log(s)))
        chi[10] = np.nan
        with pytest.raises(ValueError):
            analyze(s, O, chi=chi, window=bare())

    def test_clean_input_still_ok(self):
        # regression guard: the finite-check must not break the normal path
        s, O = self._sig()
        r = analyze(s, O, window=bare())
        assert r.sigma_c is not None

    def test_clean_callable_and_chi_paths_still_ok(self):
        # regression: the two non-array paths must still work on clean input
        s, O = self._sig()
        assert analyze(s, lambda x: np.exp(-x / 2.0), window=bare()).sigma_c is not None
        chi = np.abs(np.gradient(O, np.log(s)))
        assert analyze(s, O, chi=chi, window=bare()) is not None

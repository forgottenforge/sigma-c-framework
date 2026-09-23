# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
tau from a declared single-mode profile: the conditional yes-path (Prop 2) with a
checkable declaration. The residual check must (a) accept a genuine single mode
with the right tau, (b) reject a wrong declared profile, and (c) reject a
well-separated two-mode curve. Ratio-2 (barely separated) modes are accepted as an
effective single mode WITH a visible residual -- that is the feature by design.
"""
import numpy as np
import pytest

from sigma_c import tau_from_profile

S = np.geomspace(0.05, 200, 400)


def test_single_exponential_accepted_with_correct_tau():
    r = tau_from_profile(S, np.exp(-S / 7.0), profile="exponential")
    assert r.code == "OK_CONDITIONAL"
    assert abs(r.tau - 7.0) < 0.1
    assert r.residual < 0.01
    assert r.to_dict()["certified"] is False        # never claims certification


def test_debye_storage_accepted_with_correct_tau():
    r = tau_from_profile(S, 1.0 / (1.0 + (S / 5.0) ** 2), profile="debye")
    assert r.code == "OK_CONDITIONAL"
    assert abs(r.tau - 5.0) < 0.1
    assert r.residual < 0.02


def test_wrong_declared_profile_is_rejected():
    # Debye data, declared exponential: the shapes differ -> NOT_APPLICABLE.
    r = tau_from_profile(S, 1.0 / (1.0 + (S / 5.0) ** 2), profile="exponential")
    assert r.code == "NOT_APPLICABLE"
    assert r.tau is None
    assert r.residual > 0.1


def test_wide_two_mode_is_rejected():
    # Well-separated modes (ratio >= 3) must be rejected under the default residual.
    for ratio in (3, 5, 10, 20):
        for A in (0.3, 0.5, 1.0):
            O = np.exp(-S / 5.0) + A * np.exp(-S / (5.0 * ratio))
            r = tau_from_profile(S, O, profile="exponential")
            assert r.code != "OK_CONDITIONAL", f"leaked at ratio={ratio}, A={A}"


def test_barely_separated_two_mode_is_conditional_with_residual():
    # Ratio-2 modes look single-mode; accepted as an effective tau, but the
    # residual is non-trivial and reported (the honesty is the residual).
    O = np.exp(-S / 5.0) + 0.5 * np.exp(-S / 10.0)
    r = tau_from_profile(S, O, profile="exponential")
    if r.code == "OK_CONDITIONAL":
        assert r.residual > 0.02          # not a perfect fit, and it says so
        assert r.tau is not None


def test_regime_iii_and_ii_refuse_a_single_mode_tau():
    # monotone -> no single peak
    r = tau_from_profile(S, np.sqrt(S), profile="exponential")
    assert r.code == "NOT_RESOLVABLE" and r.tau is None


def test_summary_and_dict_never_bare_number():
    r = tau_from_profile(S, np.exp(-S / 7.0), profile="exponential")
    s = r.summary()
    assert "OK_CONDITIONAL" in s and "residual" in s and "conditional" in s.lower()
    d = r.to_dict()
    assert d["declared_profile"] == "exponential"
    assert d["relative_residual"] is not None
    assert d["library_version"] and d["schema_version"] == 1


def test_invalid_profile_name_raises():
    with pytest.raises(ValueError):
        tau_from_profile(S, np.exp(-S / 7.0), profile="not_a_profile")

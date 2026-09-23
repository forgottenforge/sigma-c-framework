# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Invariant guard for the two-field tau verdicts: the value/verdict
contradiction the adversarial checker found (regime II/III dropped tau_abscissa to
None while its status said OK) is a CLASS of bug, so it is secured with an
invariant over a diverse battery, not a single case.

The allowed cells (any result violating one fails):
  I1. window_readability == OK  =>  tau_abscissa_status == OK      (window rides on abscissa)
  I2. selfadjoint == "none"     =>  BOTH statuses NOT_APPLICABLE
  I3. selfadjoint == "declared" =>  window_readability != OK       (no probe -> no Gamma_A)
  I4. selfadjoint is None        =>  BOTH statuses NOT_IDENTIFIED   (sigma_c-only, no tau-route)
  I5. tau_abscissa_status == OK  =>  tau_abscissa is not None       (value tracks verdict)
  I6. window_readability == OK   =>  rate_error_bound is not None   (the reported number is present)

Everything runs through the public analyze().
"""
import math

import numpy as np
import pytest

from sigma_c import analyze, bare, Framework
from sigma_c.codes import Code


def _obs(sigma_max=60.0, n=1000):
    sigma = np.geomspace(0.02, sigma_max, n)
    return sigma, np.exp(-sigma / 1.0) + 0.4 * np.exp(-sigma / (1.0 / 3.0))


def _two_bump():
    sigma = np.geomspace(0.02, 100.0, 600)
    return sigma, np.exp(-sigma / 0.5) + 0.5 * np.exp(-sigma / 15.0)   # regime II


def _monotone():
    sigma = np.geomspace(0.05, 50.0, 300)
    return sigma, sigma.copy()                                        # regime III


def _sym(eigs):
    rng = np.random.default_rng(0)
    Q, _ = np.linalg.qr(rng.standard_normal((len(eigs), len(eigs))))
    return Q @ np.diag(eigs) @ Q.T


_P_REV = np.array([[0.5, 0.5, 0.0], [0.25, 0.5, 0.25], [0.0, 0.5, 0.5]])
_DRIFT = np.array([[0.0, 0.7, 0.3], [0.3, 0.0, 0.7], [0.7, 0.3, 0.0]])


def _battery():
    """Diverse (label, kwargs) covering every selfadjoint state x regime."""
    b = []
    # --- sigma_c-only (no tau-route): selfadjoint is None -----------------------
    for lbl, (sg, O) in [("only-I", _obs()), ("only-II", _two_bump()), ("only-III", _monotone())]:
        b.append((f"sigmac_{lbl}", dict(sigma=sg, O=O, window=bare())))
    # --- verified self-adjoint operators ---------------------------------------
    sg, O = _obs()
    common = dict(window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                  T_star=1.0, sigma_axis="evolution_time", gamma_A=0.5)
    b.append(("verified_gap_diag", dict(sigma=sg, O=O, operator=np.diag([1.0, 0.5, 0.1]),
                                        inner_product="euclidean", **common)))
    b.append(("verified_gap_revchain", dict(sigma=sg, O=O, operator=_P_REV,
                                             inner_product="auto", **common)))
    b.append(("verified_degenerate", dict(sigma=sg, O=O, operator=_sym([1.0, 0.5, 0.5, 0.1]),
                                           inner_product="euclidean", **common)))
    b.append(("verified_neg_dominant", dict(sigma=sg, O=O, operator=_sym([1.0, -0.9, 0.3]),
                                             inner_product="euclidean", **common)))
    b.append(("verified_nondecaying", dict(sigma=sg, O=O, operator=_sym([1.0, 1.0, 0.1]),
                                            inner_product="euclidean", **common)))
    # verified operator but the OBSERVABLE is regime II / III (the checker's bug class)
    sg2, O2 = _two_bump()
    b.append(("verified_regimeII_obs", dict(sigma=sg2, O=O2, operator=np.diag([1.0, 0.5, 0.1]),
                                            inner_product="euclidean", **common)))
    sg3, O3 = _monotone()
    b.append(("verified_regimeIII_obs", dict(sigma=sg3, O=O3, operator=np.diag([1.0, 0.5, 0.1]),
                                             inner_product="euclidean", **common)))
    # short window -> window not OK, abscissa still OK
    b.append(("verified_short_window", dict(sigma=sg, O=O, operator=np.diag([1.0, 0.5, 0.1]),
                                            inner_product="euclidean", window=bare(),
                                            framework=Framework.REVERSIBLE_MARKOV, T_star=1.0,
                                            sigma_axis="window_time", T_obs=0.2, gamma_A=0.5)))
    # --- selfadjoint = none (not self-adjoint in the declared inner product) ----
    b.append(("none_drift_auto", dict(sigma=sg, O=O, operator=_DRIFT, inner_product="auto", **common)))
    b.append(("none_revchain_euclid", dict(sigma=sg, O=O, operator=_P_REV,
                                            inner_product="euclidean", **common)))
    b.append(("none_nonnormal", dict(sigma=sg, O=O,
                                     operator=np.array([[1.0, 0.0, 0.0], [0.0, 0.5, 5.0], [0.0, 0.0, 0.3]]),
                                     inner_product="euclidean", **common)))
    # --- selfadjoint = declared (spectrum-only) --------------------------------
    b.append(("declared_gap", dict(sigma=sg, O=O, spectrum=[1.0, 0.5, 0.1], **common)))
    b.append(("declared_nondecaying", dict(sigma=sg, O=O, spectrum=[1.5, 1.0, 0.1], **common)))
    return b


@pytest.mark.parametrize("label,kwargs", _battery(), ids=[x[0] for x in _battery()])
def test_two_field_invariants_hold(label, kwargs):
    r = analyze(**kwargs)
    aa = r.tau_abscissa_status.code
    wr = r.window_readability_status.code
    sa = r.selfadjoint

    # I1: window OK rides on abscissa OK
    if wr is Code.OK:
        assert aa is Code.OK, f"{label}: window OK but abscissa {aa}"
    # I2: not self-adjoint -> NOT_APPLICABLE on both
    if sa == "none":
        assert aa is Code.NOT_APPLICABLE and wr is Code.NOT_APPLICABLE, f"{label}: {sa} {aa}/{wr}"
    # I3: declared (spectrum-only) -> window never OK
    if sa == "declared":
        assert wr is not Code.OK, f"{label}: declared but window OK"
    # I4: sigma_c-only (no tau-route) -> both NOT_IDENTIFIED
    if sa is None:
        assert aa is Code.NOT_IDENTIFIED and wr is Code.NOT_IDENTIFIED, f"{label}: {aa}/{wr}"
    # I5: abscissa OK => value present (tracks verdict; the checker's MAJOR)
    if aa is Code.OK:
        assert r.tau_abscissa is not None, f"{label}: abscissa OK but value None"
    # I6: window OK => the reported rate-error bound is present
    if wr is Code.OK:
        assert r.rate_error_bound is not None, f"{label}: window OK but no rate_error_bound"

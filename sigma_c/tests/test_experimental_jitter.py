# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
SMOKE tests for the EXPERIMENTAL jitter channel (sigma_c/experimental/jitter.py) —
the jitter theorem of The Parrot's Theorems, Section 5 + App B.3/B.4/C.

These are smoke + one paper-anchored known-answer (the five logbook rows of Sec. 5).
They are NOT the validation: the pre-registered validation is App C (peak agreement
within one grid cell, cos(v1,a) >= 0.9, shuffle-null, drift-off -> floor), which must
be run by an independent checker of different provenance. A green here does not mean
the channel is validated.
"""
import numpy as np

from sigma_c.experimental.jitter import (
    drift_covariance, jitter_rank, pair_parity_rank, pair_parity_covariance,
    tetrad, leading_coupling, cos_coupling,
    jitter_slope, two_channel_fdt, peak_agreement,
)


def test_not_reachable_from_the_kernel():
    import sigma_c
    for name in ("jitter_rank", "jitter_slope", "two_channel_fdt", "peak_agreement"):
        assert not hasattr(sigma_c, name), f"{name} leaked into the kernel namespace"


def test_nothing_experimental_claims_validated():
    bm = np.random.default_rng(0).standard_normal((300, 3))
    assert jitter_rank(bm).validated is False
    assert two_channel_fdt(0.14, 0.14, band_jitter=0.1, band_dial=0.1).validated is False
    assert peak_agreement(np.linspace(0, 1, 10), np.zeros(10), np.zeros(10)).validated is False


# --- (iv) known answer: the five logbook rows of Section 5 ------------------
def test_logbook_five_rows_slope_matches_paper():
    # Sec. 5 worked example: instrument 2 wanders, instrument 1 in step.
    M2 = np.array([0.42, 0.47, 0.50, 0.54, 0.58])
    M1 = np.array([0.9508, 0.9595, 0.9640, 0.9693, 0.9738])
    slope = jitter_slope(M1, M2)
    analytic = 2.0 / np.cosh(1 + 2 * M2.mean()) ** 2
    dial = 2.0 / np.cosh(1 + 2 * 0.5) ** 2
    assert abs(slope - 0.144) < 2e-3        # paper: jitter slope 0.144
    assert abs(analytic - 0.140) < 2e-3     # paper: analytic 0.140
    assert abs(dial - 0.141) < 2e-3         # paper: dial slope 0.141
    # (iv) cross-check: jitter slope == dial slope under detailed balance
    r = two_channel_fdt(slope, dial, band_jitter=0.03, band_dial=0.03)
    assert r.consistent_with_db is True


# --- (i)+(ii) detection + counting against the shuffle-null ------------------
def test_jitter_rank_detects_shared_drift_and_refuses_the_null():
    rng = np.random.default_rng(1)
    a = np.array([1.0, 0.8, 1.3])                 # dial susceptibilities
    delta = rng.standard_normal(600)              # ONE shared latent drift
    shared = delta[:, None] * a[None, :] + 0.05 * rng.standard_normal((600, 3))
    jr = jitter_rank(shared, n_shuffle=100)
    assert jr.detected is True and jr.rank >= 1
    indep = jitter_rank(rng.standard_normal((600, 3)), n_shuffle=100)
    assert indep.detected is False and indep.rank == 0


def test_jitter_rank_refuses_underdetermined():
    assert jitter_rank(np.random.default_rng(2).standard_normal((3, 5))).status == "insufficient_data"
    assert jitter_rank(np.zeros((10, 1))).status == "not_applicable"


def test_channel1_pair_parities_shot_level_cause():
    # App C Channel 1: a shared per-shot cause correlates the pair parities; the
    # per-shot covariance detects it, while independent pairs sit at the null.
    rng = np.random.default_rng(5)
    z = rng.choice([-1.0, 1.0], size=4000)              # one shared per-shot cause
    flips = rng.random((4000, 3)) < 0.2
    shared = z[:, None] * np.where(flips, -1.0, 1.0)    # correlated +/-1 pair parities
    assert pair_parity_rank(shared, n_shuffle=200).detected is True
    assert pair_parity_covariance(shared).shape == (3, 3)
    # detection is a ~alpha-level event under the null, so MOST independent draws
    # do not fire (a single draw is a coin flip against the 5% null; test the rate).
    fires = sum(pair_parity_rank(rng.choice([-1.0, 1.0], size=(2000, 3)),
                                 n_shuffle=100, seed=s).detected for s in range(12))
    assert fires <= 3


# --- (B.4) tetrad -----------------------------------------------------------
def test_tetrad_one_factor_holds_two_factors_break():
    a = np.array([1.0, 0.7, 1.2]); sy = 2.0
    one = sy * np.outer(a, a) + np.diag([0.1, 0.1, 0.1])
    t1 = tetrad(one)
    assert t1["status"] == "reported" and 0.7 <= t1["ratio"] <= 1.4
    # a second independent factor on instruments 2,3 breaks the identity
    b = np.array([0.0, 1.0, 1.0])
    two = sy * np.outer(a, a) + 1.5 * np.outer(b, b) + np.diag([0.1, 0.1, 0.1])
    t2 = tetrad(two)
    assert not (0.7 <= t2["ratio"] <= 1.4)


# --- (iii) identification ---------------------------------------------------
def test_leading_coupling_aligns_with_dial_susceptibilities():
    rng = np.random.default_rng(3)
    a = np.array([1.0, 0.8, 1.3])
    delta = rng.standard_normal(2000)
    shared = delta[:, None] * a[None, :] + 0.02 * rng.standard_normal((2000, 3))
    C = drift_covariance(shared)
    v1 = leading_coupling(C)
    assert cos_coupling(v1, a) >= 0.9              # App C criterion


# --- (iv) two_channel_fdt is reported, never pass/fail ----------------------
def test_two_channel_reported_not_verdict():
    assert two_channel_fdt(0.14, 0.145).consistent_with_db is None          # no bands
    assert two_channel_fdt(0.14, 0.142, band_jitter=0.03, band_dial=0.03).consistent_with_db is True
    assert two_channel_fdt(0.14, 0.40, band_jitter=0.03, band_dial=0.03).consistent_with_db is False
    assert two_channel_fdt(0.14, None).status == "needs_both_channels"


# --- App C peak agreement ---------------------------------------------------
def test_peak_agreement_same_grid_point():
    g = np.linspace(0, 1, 20)
    dial = np.exp(-((g - 0.674) ** 2) / 0.002)     # dial channel peaks at gamma_c
    jit = np.exp(-((g - 0.66) ** 2) / 0.002)       # jitter channel one cell away
    pa = peak_agreement(g, dial, jit, cos_v1_a=0.97)
    assert pa.within_one_cell is True and pa.validated is False

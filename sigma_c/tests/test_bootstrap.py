# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Statistical CI for sigma_c from replicates -- a new measurement statement, so it
carries its own coverage validation. The percentile bootstrap of a peak LOCATION
under-covers somewhat (documented); we validate that (a) it is clearly distinct
from resolution_band, and (b) its measured coverage on synthetic replicates whose
curves resolve is in a validated band (>= 0.85 at nominal 0.95), not a fabricated
95%.
"""
import numpy as np

from sigma_c import analyze
from sigma_c.bootstrap import bootstrap_sigma_c


def test_result_shape_and_honesty_flags():
    s = np.geomspace(0.5, 50, 80)
    base = 1.0 / (1.0 + np.exp(-(np.log(s) - np.log(5.0)) * 6.0))
    rng = np.random.default_rng(0)
    reps = np.array([base + rng.normal(0, 0.01, s.size) for _ in range(20)])
    r = bootstrap_sigma_c(s, reps, n_boot=200, seed=0)
    assert r["is_confidence_interval"] is True     # it IS a statistical CI ...
    assert r["approximate"] is True                # ... but honestly flagged approximate
    assert "resolution_band" in r["note"]          # points at the distinction
    assert r["ci_lo"] is not None and r["ci_lo"] < r["sigma_c"] < r["ci_hi"]
    assert r["n_resolved"] <= r["n_boot"]


def test_bad_shape_raises():
    import pytest
    from sigma_c.errors import InvalidInputError
    with pytest.raises(InvalidInputError):
        bootstrap_sigma_c(np.linspace(1, 5, 10), np.zeros((3, 7)))  # M mismatch


def test_measured_coverage_is_in_the_validated_band():
    # Coverage validation: many synthetic replicate-sets, true sigma_c known; the
    # nominal-95% interval must cover the truth in a validated band. We assert >=
    # 0.85 (argmax bootstrap under-covers; we do NOT claim a clean 95%).
    s = np.geomspace(0.5, 50, 80)
    s0, k, noise, R = 5.0, 6.0, 0.01, 20
    base = 1.0 / (1.0 + np.exp(-(np.log(s) - np.log(s0)) * k))
    true_sc = float(analyze(s, base).sigma_c)
    rng = np.random.default_rng(1)
    N = 120
    covered = resolved = 0
    for i in range(N):
        reps = np.array([base + rng.normal(0, noise, s.size) for _ in range(R)])
        r = bootstrap_sigma_c(s, reps, n_boot=200, seed=i)
        if r["ci_lo"] is not None:
            resolved += 1
            if r["ci_lo"] <= true_sc <= r["ci_hi"]:
                covered += 1
    assert resolved >= 0.9 * N, f"only {resolved}/{N} resolved"
    coverage = covered / resolved
    assert coverage >= 0.85, f"coverage {coverage:.2%} below validated floor 0.85"

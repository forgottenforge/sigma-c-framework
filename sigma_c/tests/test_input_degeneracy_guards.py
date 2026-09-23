# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Input-degeneracy guards (6.0): a constant/near-constant observable and a peak
count beyond the declared resolution must NOT read as a spurious OK. Siblings of
the NaN guard. Recorded in THEOREM_MAP as a deliberate verdict change from the
pre-6.0 behaviour (a constant curve used to return OK with phantom peaks).
"""
import numpy as np
import pytest

from sigma_c import analyze


SIGMA = np.linspace(1.0, 5.0, 50)


def _code(O, sig=None, **kw):
    s = SIGMA if sig is None else sig
    return analyze(s, O, **kw).to_dict()["sigma_c_status"]["code"]


def test_exact_constant_is_not_identified():
    r = analyze(SIGMA, np.ones(50)).to_dict()
    assert r["sigma_c_status"]["code"] == "NOT_IDENTIFIED"
    assert r["regime"]["degenerate_O"] is True
    assert r["sigma_c"] is None
    # the reason names the cause and the remediation points somewhere
    assert "no variation" in r["sigma_c_status"]["reason"]
    assert r["sigma_c_status"]["remediation"]


def test_zero_observable_is_not_identified():
    assert _code(np.zeros(50)) == "NOT_IDENTIFIED"


def test_constant_plus_tiny_noise_is_not_resolvable_via_peak_ceiling():
    rng = np.random.default_rng(0)
    code = _code(np.ones(50) + rng.normal(0, 1e-12, 50))
    assert code == "NOT_RESOLVABLE"  # many phantom peaks > ceiling


def test_straight_line_stays_not_resolvable_regime_iii():
    assert _code(np.linspace(1.0, 5.0, 50)) == "NOT_RESOLVABLE"


def test_real_small_bump_above_floor_is_still_found():
    # A genuine but small logistic step (range ~1e-3, orders of magnitude above
    # the eps*max|O| noise floor) MUST resolve to OK -- the guard must not eat it.
    O = 1.0 + 1e-3 / (1.0 + np.exp(-(SIGMA - 3.0) * 8.0))
    r = analyze(SIGMA, O).to_dict()
    assert r["sigma_c_status"]["code"] == "OK"
    assert r["regime"]["geometric"] == "I_geom"


def test_clean_logistic_unchanged_ok():
    O = 1.0 / (1.0 + np.exp(-(SIGMA - 3.0)))
    assert _code(O) == "OK"


def test_peak_ceiling_is_a_declared_convention():
    # A genuinely stable multimode (three clean, well-separated log-Gaussian bumps
    # -> 6 chi-lobes) is convention-window stable, so the ONLY thing refusing it is
    # the declared ceiling: NOT_RESOLVABLE at default 5 (peaks still listed),
    # OK when the caller raises the ceiling. This isolates the ceiling convention
    # from the (separate) noise / stability gate.
    x = np.geomspace(0.3, 30, 300)
    def lg(c):
        return np.exp(-((np.log(x) - np.log(c)) ** 2) / 0.03)
    O = lg(1.0) + lg(5.0) + lg(25.0)
    d_default = analyze(x, O).to_dict()
    assert d_default["regime"]["resolved_peak_count"] > 5
    assert d_default["sigma_c_status"]["code"] == "NOT_RESOLVABLE"
    assert isinstance(d_default["sigma_c"], list)  # peaks still listed
    d_loose = analyze(x, O, max_resolved_peaks=1000).to_dict()
    assert d_loose["sigma_c_status"]["code"] == "OK"  # stable + under ceiling


def test_white_noise_never_ok_across_convention_grid():
    # Decision: an OK that exists only under a narrow prominence convention is not
    # a measurement. Pure white noise must NEVER return OK, at ANY min_prominence_
    # ratio (a high r that keeps only the tallest noise spike is caught because the
    # peak count is not stable across the [r0/2, 2*r0] octave). Scale-free.
    s = np.linspace(1.0, 10.0, 80)
    r_grid = np.linspace(0.05, 0.99, 40)
    for seed in range(15):
        O = np.random.default_rng(seed).normal(0, 1, 80)
        for r in r_grid:
            assert _code(O, sig=s, min_prominence_ratio=float(r)) != "OK", (
                f"white noise returned OK at seed={seed}, r={r:.3f}"
            )


def test_convention_stable_signal_is_ok():
    # The mirror: a clean single bump is convention-window stable and stays OK
    # across a wide range of r (the gate must not eat real signal).
    s = np.linspace(0.5, 10.0, 60)
    O = 1.0 / (1.0 + np.exp(-(s - 4.0)))
    for r in (0.05, 0.1, 0.2, 0.4):
        assert _code(O, sig=s, min_prominence_ratio=float(r)) == "OK"


def test_tiny_linear_ramp_at_noise_floor_is_refused():
    # Adversarial-review regression: O = 1 + 1e-14*arange just above the ULP gate
    # used to slip an edge-differencing artifact through as a single-peak OK. It
    # must be refused (its only 'peak' is at the numerical floor, boundary-adjacent).
    O = 1.0 + 1e-14 * np.arange(50)
    assert _code(O) in ("NOT_RESOLVABLE", "NOT_IDENTIFIED")


def test_honest_amplitude_cutoff_is_consistent():
    # The cutoff sits near 1e-12 relative amplitude (~1000x machine eps): a real
    # bump above it resolves; below it, signal is indistinguishable from the
    # log-derivative's rounding noise and is refused -- CONSISTENTLY (same verdict
    # class for a sub-floor bump and a tiny ramp, no scale dependence).
    def bump(a):
        return 1.0 + a / (1.0 + np.exp(-(SIGMA - 3.0) * 8.0))
    assert _code(bump(1e-12)) == "OK"
    assert _code(bump(1e-13)) in ("NOT_RESOLVABLE", "NOT_IDENTIFIED")
    # scale invariance of the cutoff: multiply O by 1e9, same verdicts
    assert _code(1e9 * bump(1e-12)) == "OK"
    assert _code(1e9 * bump(1e-13)) in ("NOT_RESOLVABLE", "NOT_IDENTIFIED")


def test_many_peaks_still_lists_them():
    rng = np.random.default_rng(2)
    O = np.sin(np.linspace(0, 40, 80)) + rng.normal(0, 0.01, 80)
    d = analyze(np.linspace(1.0, 5.0, 80), O).to_dict()
    assert d["regime"]["resolved_peak_count"] > 5
    assert d["sigma_c_status"]["code"] == "NOT_RESOLVABLE"
    assert isinstance(d["sigma_c"], list) and len(d["sigma_c"]) > 5

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Unit tests for the TWO named convention-operations that replace the archived
§8 "trichotomy" formalism (paper Sätze 8.3/8.6/8.11/8.16; a private math
archive, not shipped; THEOREM_MAP.md marks the §8 labels ARCHIVED).

The live code runs, and this file pins, exactly two operations:

  (1) GEOMETRIC = a PEAK COUNT under the declared convention `min_prominence_ratio`
      (geometric_trichotomy). Not a theorem; the verdict rides on the convention,
      and the sweep test below pins the flip point at which it does.

  (2) SPECTRAL = an ISOLATED-LEADING-GAP boolean on a supplied spectrum
      (spectral_attribution): 'gap' iff |lambda_2| < |lambda_1| strictly,
      'no_gap' otherwise, None when no spectrum. Spectrum-only, NOT
      observable-faithful.

Plus a test for the ARCHIVED proof-status token rendered by theorem_map.cite().
"""
import numpy as np

from sigma_c.core.trichotomy import (
    geometric_trichotomy,
    spectral_attribution,
)


# ---------------------------------------------------------------------------
# (1) PEAK-COUNT operation on known I / II / III profiles
# ---------------------------------------------------------------------------

def test_peak_count_regime_I_single_peak():
    # one interior maximum -> I_geom
    chi = np.array([0.0, 0.5, 1.0, 0.5, 0.1])
    geom, peaks = geometric_trichotomy(chi)
    assert geom == "I_geom"
    assert len(peaks) == 1


def test_peak_count_regime_II_two_peaks():
    # two interior maxima both above the default 10% floor -> II_geom
    chi = np.array([0.0, 0.5, 1.0, 0.4, 0.1, 0.6, 0.9, 0.5, 0.05])
    geom, peaks = geometric_trichotomy(chi)
    assert geom == "II_geom"
    assert len(peaks) == 2


def test_peak_count_regime_III_no_interior_peak():
    # strictly monotone -> no interior maximum -> III_geom
    chi = np.array([0.1, 0.3, 0.5, 0.7, 0.9])
    geom, peaks = geometric_trichotomy(chi)
    assert geom == "III_geom"
    assert peaks == []


# ---------------------------------------------------------------------------
# PRECISION #1: the min_prominence_ratio sweep flip point.
# On ONE fixed two-bump observable, sweep the declared prominence convention
# and pin the ratio at which the peak count (hence I/II/III -> sigma_c type)
# flips. This IS the documented sensitivity of the geometric operation.
# ---------------------------------------------------------------------------

# Fixed multi-bump profile: main peak height 1.0 (index 2), secondary bump
# height 0.40 (index 6). Both are strict interior local maxima. The secondary
# bump sits at exactly 40% of max chi, so it is counted (II_geom) while the
# prominence floor r <= 0.40 and dropped (I_geom) once r > 0.40.
_SWEEP_CHI = np.array([0.0, 0.5, 1.0, 0.5, 0.2, 0.3, 0.40, 0.3, 0.1])


def _regime_at(ratio: float) -> str:
    return geometric_trichotomy(_SWEEP_CHI, min_prominence_ratio=ratio)[0]


def test_min_prominence_ratio_sweep_flip_point():
    # below the secondary bump's relative height -> two peaks (II_geom)
    assert _regime_at(0.35) == "II_geom"
    # above it -> the secondary bump drops out -> one peak (I_geom)
    assert _regime_at(0.45) == "I_geom"

    # Locate the flip by sweeping and pin it at the secondary bump's relative
    # height 0.40 (a bump at 40% of max chi survives r<=0.40, dies at r>0.40).
    ratios = np.linspace(0.30, 0.50, 2001)
    regimes = [_regime_at(float(r)) for r in ratios]
    flip_idx = next(
        i for i in range(1, len(regimes)) if regimes[i] != regimes[i - 1]
    )
    flip_ratio = float(ratios[flip_idx])
    # the count flips II -> I as the ratio crosses ~0.40
    assert regimes[flip_idx - 1] == "II_geom"
    assert regimes[flip_idx] == "I_geom"
    assert abs(flip_ratio - 0.40) < 1e-3, flip_ratio


# ---------------------------------------------------------------------------
# (2) ISOLATED-LEADING-GAP boolean on known gap / no-gap spectra
# ---------------------------------------------------------------------------

def test_gap_test_gap_present():
    # |lambda_2| = 0.5 < |lambda_1| = 1.0 -> gap
    assert spectral_attribution([1.0, 0.5, 0.1]) == "gap"


def test_gap_test_no_gap_leading_degenerate():
    # |lambda_2| reaches |lambda_1| -> no isolated leading gap
    assert spectral_attribution([1.0, 1.0, 0.5]) == "no_gap"


def test_gap_test_single_modulus_is_no_gap():
    # fewer than two moduli -> nothing to separate from lambda_1 -> no_gap
    assert spectral_attribution([1.0]) == "no_gap"


def test_gap_test_none_spectrum_is_none():
    assert spectral_attribution(None) is None


def test_gap_test_complex_moduli_sorted():
    # unsorted, complex: moduli {1.0, 0.6, 0.3} -> gap
    assert spectral_attribution([0.3j, 1.0, -0.6]) == "gap"


def test_gap_test_1e_12_tolerance():
    # just inside the tolerance band -> counts as no_gap (not a real gap)
    assert spectral_attribution([1.0, 1.0 - 1e-13]) == "no_gap"
    # clearly below the tolerance -> a real gap
    assert spectral_attribution([1.0, 1.0 - 1e-6]) == "gap"


# ---------------------------------------------------------------------------
# ARCHIVED proof-status token (theorem_map)
# ---------------------------------------------------------------------------

def test_archived_token_registered_and_rendered():
    from sigma_c import theorem_map
    from sigma_c.theorem_map import cite, reload_map

    assert "ARCHIVED" in theorem_map._PROOF_STATES
    reload_map()
    # the four §8 labels are ARCHIVED in THEOREM_MAP.md; cite() must render the
    # status like any other non-PROVED state (not silently as proved).
    rendered = cite("thm:trichotomy-geometric")
    assert "ARCHIVED" in rendered, rendered

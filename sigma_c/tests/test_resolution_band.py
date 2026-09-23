# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
P1-UX #1/#2/#4: the resolution/sensitivity band on sigma_c, its band-coupled
digits, and the plain-language + 4-code + tau-provenance surface in summary().

DISCIPLINE: everything runs through the public analyze()/Result, never against
core.resolution directly (the half-fix lesson).
"""
import numpy as np
import pytest

from sigma_c import analyze, gamma_k, bare, Framework
from sigma_c.core.susceptibility import find_interior_maxima


def _single_peak(n=400):
    # one susceptibility bump: a single exponential relaxation
    sigma = np.geomspace(0.05, 50.0, n)
    O = np.exp(-sigma / 3.0)
    return sigma, O


def _two_bumps(n=600):
    # two separate modes; the weaker bump sits BELOW the declared prominence (0.10)
    # so the declared regime is I (scalar sigma_c) but a lower r flips it to II --
    # a genuine per-input flip point (~0.079), not a hardcoded r.
    sigma = np.geomspace(0.02, 100.0, n)
    O = np.exp(-sigma / 0.5) + 0.08 * np.exp(-sigma / 15.0)
    return sigma, O


def _monotone(n=300):
    # O = sigma -> chi = |sigma * dO/dsigma| = sigma, strictly monotone -> no
    # interior chi peak -> regime III (sigma_c is None). (1/(1+sigma) does NOT
    # work: its chi = sigma/(1+sigma)^2 peaks at sigma=1.)
    sigma = np.geomspace(0.05, 50.0, n)
    return sigma, sigma.copy()


# ===========================================================================
# The band exists, is a band (not a CI), and brackets sigma_c
# ===========================================================================

def test_resolution_band_present_and_brackets_sigma_c():
    sigma, O = _single_peak()
    r = analyze(sigma, O, window=gamma_k(2))
    rb = r.resolution_band()
    assert rb is not None
    assert rb["is_confidence_interval"] is False
    assert rb["grid_cell_width"] > 0
    assert rb["band_lo"] < r.sigma_c < rb["band_hi"]
    assert rb["significant_digits"] >= 1
    # the value interval is exactly sigma_c +/- band_abs
    assert rb["band_hi"] - rb["band_lo"] == pytest.approx(2 * rb["band_abs"])


def test_resolution_band_in_to_dict_matches_method():
    sigma, O = _single_peak()
    r = analyze(sigma, O, window=gamma_k(2))
    d = r.to_dict()
    assert d["resolution_band"] == r.resolution_band()


def test_resolution_band_none_in_regime_III():
    sigma, O = _monotone()
    r = analyze(sigma, O, window=bare())
    assert r.sigma_c is None
    assert r.resolution_band() is None
    assert r.to_dict()["resolution_band"] is None


# ===========================================================================
# The prominence sweep is PER-INPUT (not a hardcoded r=0.40)
# ===========================================================================

def test_prominence_flip_is_computed_from_this_input():
    sigma, O = _two_bumps()
    r = analyze(sigma, O, window=bare())
    rb = r.resolution_band()
    assert rb is not None
    lo, hi = rb["prominence_stable_range"]
    # a real per-input flip exists somewhere in (0,1): the weaker bump crosses the
    # counting threshold. At least one side must be an interior flip point.
    assert (rb["prominence_flip_low"] is not None
            or rb["prominence_flip_high"] is not None)
    # and the reported flip is genuine: applying it changes the peak count vs the
    # declared convention -- i.e. it was read off THIS observable.
    declared = rb["min_prominence_ratio"]
    current = len(find_interior_maxima(sigma, r._profile_chi, min_prominence_ratio=declared))
    if rb["prominence_flip_high"] is not None:
        flipped = len(find_interior_maxima(sigma, r._profile_chi,
                                           min_prominence_ratio=rb["prominence_flip_high"]))
        assert flipped != current
    if rb["prominence_flip_low"] is not None:
        flipped = len(find_interior_maxima(sigma, r._profile_chi,
                                           min_prominence_ratio=rb["prominence_flip_low"]))
        assert flipped != current


def test_stable_range_is_full_when_single_dominant_peak():
    sigma, O = _single_peak()
    r = analyze(sigma, O, window=gamma_k(2))
    rb = r.resolution_band()
    lo, hi = rb["prominence_stable_range"]
    # one dominant peak: the count (1) never changes across the whole r range
    assert lo == 0.0 and hi == 1.0
    assert rb["prominence_flip_low"] is None and rb["prominence_flip_high"] is None


# ===========================================================================
# Significant digits are coupled to the band (no false precision)
# ===========================================================================

def test_significant_digits_coupled_to_band():
    sigma, O = _single_peak()
    r = analyze(sigma, O, window=gamma_k(2))
    rb = r.resolution_band()
    import math
    # the displayed value must not carry a place finer than the band's order
    order = math.floor(math.log10(rb["band_abs"]))
    # rounding sigma_c to the band's order reproduces the displayed rounded value
    assert rb["sigma_c_rounded"] == pytest.approx(round(rb["sigma_c"], -order))
    # and the display string has no more decimals than the band justifies
    if "." in rb["sigma_c_display"]:
        decimals = len(rb["sigma_c_display"].split(".")[1])
        assert decimals == max(0, -order)


# ===========================================================================
# summary() carries the plain language, the band, the legend, tau provenance
# ===========================================================================

def test_summary_has_plainlang_band_and_legend():
    sigma, O = _single_peak()
    r = analyze(sigma, O, window=gamma_k(2))
    s = r.summary()
    assert "what this reads" in s
    assert "resolution band, NOT a CI" in s
    assert "what each verdict means" in s          # plain-language 4-code legend
    assert "the number stands" in s and "wrong tool for this question" in s
    assert "NOT_IDENTIFIED" in s
    # tau is the fallen bridge on this path -> provenance must be visible
    assert "sigma_c/rho_star" in s


def test_summary_shows_certified_spectral_provenance_on_tau_route():
    import math
    sigma = np.geomspace(0.02, 60.0, 1000)
    O = 1.0 * np.exp(-sigma / 1.0) + 0.4 * np.exp(-sigma / (1.0 / 3.0))
    lam = [1.0, math.exp(-1.0), math.exp(-3.0)]
    r = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                spectrum=lam, gamma_A=0.5, T_star=1.0, sigma_axis="evolution_time")
    s = r.summary()
    assert "certified spectral" in s

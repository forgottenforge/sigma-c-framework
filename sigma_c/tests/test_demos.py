# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Demos are validation: each demo exposes run() returning known-answer facts, and
the suite asserts the instrument hits the known answer (a demo with
a known answer proves the instrument measures; the demo scripts also run end-to-end
via test_examples_run.py).
"""
import math

import pytest


def test_reversible_tau_demo_certifies_and_refuses():
    from sigma_c.examples import demo_reversible_tau as d
    r = d.run()
    # reversible chain: self-adjointness verified, certified tau matches -T*/log(lambda2)
    assert r["reversible"]["selfadjoint"] == "verified"
    assert r["reversible"]["abscissa_code"] == "OK"
    assert r["reversible"]["tau_abscissa"] == pytest.approx(r["tau_expected"], rel=1e-6)
    # drift walk: euclidean-normal but not reversible -> the honest refusal on BOTH
    assert r["drift"]["selfadjoint"] == "none"
    assert r["drift"]["abscissa_code"] == "NOT_APPLICABLE"
    assert r["drift"]["window_code"] == "NOT_APPLICABLE"


def test_raw_shot_demo_recovers_planted_tau_within_band():
    from sigma_c.examples import demo_raw_shot as d
    r = d.run()
    # sigma_c must equal the planted tau (chi peaks at sigma=tau, closed form) ...
    assert r["sigma_c"] == pytest.approx(r["tau_true"], rel=2e-3)
    # ... and the planted tau lies inside the (non-CI) resolution band
    assert r["band_lo"] <= r["tau_true"] <= r["band_hi"]
    assert r["sigma_c_status"] == "OK"


def test_percolation_demo_recovers_p_c_near_half():
    from sigma_c.examples import demo_percolation as d
    r = d.run()                          # deterministic (fixed seeds)
    # sigma_c lands near the exact bond threshold p_c = 1/2 (finite-grid, marches in)
    assert r["status"] == "OK"
    assert abs(r["sigma_c"] - r["p_c_exact"]) < 0.05, r["sigma_c"]


def test_no_scale_demo_refuses_honestly():
    from sigma_c.examples import demo_no_scale as d
    r = d.run()
    # a monotone curve has no interior peak -> no scale to name -> NOT_RESOLVABLE
    assert r["sigma_c"] is None
    assert r["status"] == "NOT_RESOLVABLE"


def test_seismic_demo_finds_Mc_not_physics():
    from sigma_c.examples import demo_seismic as d
    r = d.run()
    # generator is genuine Gutenberg-Richter (slope ~ -b) ...
    assert abs(r["gr_slope"] - (-r["b_true"])) < 0.15, r["gr_slope"]
    # ... standard tool recovers b ~ 1.0 above Mc ...
    assert abs(r["b_hat"] - r["b_true"]) < 0.15, r["b_hat"]
    # ... and sigma_c lands on the COMPLETENESS Mc (the detector's edge), not the earth
    assert not isinstance(r["sigma_c"], list)
    assert abs(r["sigma_c"] - r["mc_true"]) < 0.25, r["sigma_c"]

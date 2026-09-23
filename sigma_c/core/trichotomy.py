# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Regime classifier — two convention-operations (NOT theorems).

The paper's §8 "trichotomy" (Sätze 8.3/8.6/8.11/8.16) is a paper-only formalism
archived privately, not shipped (see THEOREM_MAP.md: those §8 labels are
ARCHIVED, not live theorems). The live code does NOT realise those theorems; it
runs two named, declared convention-operations:

  1. GEOMETRIC = a PEAK COUNT under the DECLARED convention `min_prominence_ratio`.
     Count strict interior maxima of chi_O whose height clears the declared
     fraction of max chi; 0 -> III_geom, 1 -> I_geom, >=2 -> II_geom. This is a
     COUNT under a declared convention, NOT a theorem: the verdict depends on
     `min_prominence_ratio` (default 0.10 = 10% of max chi; a bump below it is
     not counted), so that convention is recorded in the Trichotomy verdict (it
     governs the outcome and must be visible). See the min_prominence_ratio
     sensitivity sweep pinned in tests (test_trichotomy_ops.py).

  2. SPECTRAL = an ISOLATED-LEADING-GAP boolean on a supplied spectrum:
     "gap" iff |lambda_2| < |lambda_1| strictly. This is a SPECTRUM-ONLY test;
     it is NOT observable-faithful. Observable-faithfulness runs exclusively
     through Gamma_A (probe-relative), never through this gap. A "gap" is only a
     *candidate* for the tau OK-gate; OK still requires gamma_A AND T > t*.

Operational diagnostic (def:noise-floor-diagnostic, under eta_O):
  measurement-side check whether candidate peak amplitude is above the
  signal/noise floor.

The layers are reported as a single Trichotomy dict (Enforcement 1, Item A of
the prescription).
"""
from __future__ import annotations
from typing import Optional, Sequence, Tuple

import numpy as np

from sigma_c.core.susceptibility import find_interior_maxima
from sigma_c.result import Trichotomy


# Declared conventions for the two input-degeneracy guards (NOT theorems).
# The two guards are DECOUPLED (an adversarial-review finding): the degenerate
# gate only catches an observable that is constant to a few ULP, and the per-peak
# floor is the workhorse that rejects rounding-noise "peaks".
#
# DEGENERATE_ULP_FACTOR: O is "constant to numerical precision" iff
#   range(O) <= DEGENERATE_ULP_FACTOR * eps * max|O|.
# No factor of n: range = max-min does NOT accumulate over the sample count, and
# an n-factor wrongly refused real (small but resolvable) signal at large n.
DEGENERATE_ULP_FACTOR: float = 16.0
# NOISE_FLOOR_FACTOR: a peak in chi = |dO/dlog sigma| must clear
#   NOISE_FLOOR_FACTOR * eps * max|O| / min(dlog)
# to count. eps*max|O| is the rounding error in O; dividing by the SMALLEST
# log-step (strictest across the grid, not a global median) makes the floor at
# least as strict as any local point, so a boundary-differencing artifact at the
# fine end of a linear grid cannot sneak a peak through. 128 sits well below a
# real signal (ratio ~20+ at 1e-12 relative amplitude) and well above both pure
# rounding noise and edge artifacts (ratio ~2).
NOISE_FLOOR_FACTOR: float = 128.0
# Default ceiling on resolved peaks; above it the peak field is noise, not
# structure -> NOT_RESOLVABLE (peaks still listed). Tunable per call.
DEFAULT_MAX_RESOLVED_PEAKS: int = 5


def geometric_trichotomy(
    chi: np.ndarray,
    *,
    min_prominence_ratio: float = 0.10,
    chi_abs_floor: float = 0.0,
) -> Tuple[str, list]:
    """
    Return (geometric_regime, peak_indices).

    geometric_regime ∈ {"I_geom", "II_geom", "III_geom"}.

    This is a PEAK COUNT under the DECLARED convention `min_prominence_ratio`
    (see module docstring) — a named convention-operation, NOT a theorem. The
    threshold decides whether a bump counts as a peak, hence the regime.

    SENSITIVITY of `min_prominence_ratio`: on a fixed multi-bump observable the
    peak count (and hence I/II/III -> the sigma_c type) flips as this ratio
    sweeps past a bump's relative height. The flip point is pinned in
    test_trichotomy_ops.py::test_min_prominence_ratio_sweep_flip_point (a
    two-bump profile with a secondary bump at ~40% of max chi flips
    II_geom -> I_geom as the ratio crosses ~0.40). The verdict therefore rides
    on the declared value; it is surfaced in the Trichotomy verdict, not hidden.
    """
    peaks = find_interior_maxima(
        np.arange(len(chi)), chi,
        min_prominence_ratio=min_prominence_ratio,
        chi_abs_floor=chi_abs_floor,
    )
    if len(peaks) == 0:
        return "III_geom", []
    if len(peaks) == 1:
        return "I_geom", peaks
    return "II_geom", peaks


def operational_floor_check(
    chi: np.ndarray,
    O: np.ndarray,
    eta_O: float,
) -> bool:
    """
    True iff the candidate peak amplitude is below eta_O * ||O||_inf.

    Cite: def:noise-floor-diagnostic. Measurement-side diagnostic.
    """
    if eta_O <= 0:
        return False
    O_range = float(np.max(O) - np.min(O))
    if O_range == 0:
        return True
    candidate_peak = float(np.max(chi))
    return candidate_peak < eta_O * O_range


def spectral_attribution(
    spectrum: Optional[Sequence[complex]],
) -> Optional[str]:
    """
    Isolated-leading-gap test on a supplied spectrum. SPECTRUM-ONLY boolean.

    spectrum: full spectrum {lambda_1, lambda_2, ...} (any order; sorted here by
    decreasing |lambda|), OR None when no spectral data is available.

    Returns:
      "gap"    - isolated leading gap present: the leading non-trivial modulus is
                 strictly below |lambda_1|, i.e. |lambda_2| < |lambda_1| with the
                 existing 1e-12 tolerance.
      "no_gap" - no isolated leading gap / power-law: the leading non-trivial
                 modulus reaches |lambda_1| (or fewer than two moduli supplied).
      None     - no spectrum supplied.

    NOT observable-faithful. This test looks ONLY at the spectrum; it says
    nothing about whether the probe reads the slow mode. Observable-faithfulness
    is decided EXCLUSIVELY by Gamma_A (probe-relative, in codes.tau_abscissa_verdict),
    never by this gap. A "gap" verdict is at most a *candidate* for the tau
    OK-gate; OK additionally requires gamma_A and the tail window T > t*. This is
    a named convention-operation, NOT a theorem (paper §8.6 is ARCHIVED; see THEOREM_MAP.md).
    """
    if spectrum is None:
        return None
    abs_spec = sorted((abs(complex(x)) for x in spectrum), reverse=True)
    if len(abs_spec) < 2:
        return "no_gap"  # no non-trivial modulus to separate from lambda_1
    lam1 = abs_spec[0]
    lam2 = abs_spec[1]
    # Isolated leading gap: leading non-trivial modulus strictly under |lambda_1|.
    if lam2 < lam1 * (1 - 1e-12):
        return "gap"
    return "no_gap"


def _numerical_noise_floors(
    sigma: np.ndarray, O: np.ndarray, n: int
) -> Tuple[bool, float]:
    """(degenerate_O, chi_abs_floor) from the observable's variation and the grid.

    degenerate_O: O is constant to a few ULP,
        range(O) <= DEGENERATE_ULP_FACTOR * eps * max|O|  (or O == 0).
    chi_abs_floor: the per-peak floor NOISE_FLOOR_FACTOR * eps * max|O| / min(dlog),
        i.e. the rounding-error level of chi = |dO/dlog sigma| at the strictest
        (finest) log-step; peaks below it are indistinguishable from the
        log-derivative's own noise. `n` is unused (kept for call compatibility).
    """
    eps = float(np.finfo(float).eps)
    O_scale = float(np.max(np.abs(O))) if O.size else 0.0
    O_range = float(np.max(O) - np.min(O)) if O.size else 0.0
    # Guard 1 (constant to a few ULP) -- NO factor of n (max-min does not grow
    # with the sample count).
    degenerate = (O_scale == 0.0) or (
        O_range <= DEGENERATE_ULP_FACTOR * eps * O_scale
    )
    # Guard 2 floor uses the SMALLEST log-step (strictest across the grid), so a
    # boundary-differencing artifact at the fine end cannot pass.
    dlog_min = 0.0
    if sigma.size >= 2 and np.all(sigma > 0):
        dlog_min = float(np.min(np.abs(np.diff(np.log(sigma)))))
    chi_abs_floor = (
        NOISE_FLOOR_FACTOR * eps * O_scale / dlog_min if dlog_min > 0 else 0.0
    )
    return degenerate, chi_abs_floor


def classify(
    sigma: np.ndarray,
    chi: np.ndarray,
    O: np.ndarray,
    *,
    spectrum: Optional[Sequence[complex]] = None,
    eta_O: float = 0.0,
    min_prominence_ratio: float = 0.10,
    max_resolved_peaks: int = DEFAULT_MAX_RESOLVED_PEAKS,
) -> Tuple[Trichotomy, list]:
    """
    Full classification: peak-count operation + spectral gap-boolean + floor,
    guarded by two input-degeneracy checks (a constant/near-constant observable,
    and a peak count beyond the declared resolution).

    Returns (Trichotomy, peak_indices).
    """
    O = np.asarray(O, dtype=float)
    sigma = np.asarray(sigma, dtype=float)
    degenerate_O, chi_abs_floor = _numerical_noise_floors(sigma, O, len(chi))

    spec = spectral_attribution(spectrum)
    floor = operational_floor_check(chi, O, eta_O)

    if degenerate_O:
        # O carries no signal above the numerical noise; any peak is float noise
        # in the log-derivative. Force regime III (sigma_c = None) and flag it so
        # the verdict is NOT_IDENTIFIED (remediable), not a spurious OK.
        tri = Trichotomy(
            geometric="III_geom", spectral=spec,
            operational_floor_triggered=floor, degenerate_O=True,
            resolved_peak_count=0, max_resolved_peaks=max_resolved_peaks,
            peak_count_over_ceiling=False, chi_abs_floor=chi_abs_floor,
            eta_O=eta_O, min_prominence_ratio=min_prominence_ratio,
        )
        return tri, []

    geom, peaks = geometric_trichotomy(
        chi, min_prominence_ratio=min_prominence_ratio,
        chi_abs_floor=chi_abs_floor,
    )
    over_ceiling = len(peaks) > max_resolved_peaks

    tri = Trichotomy(
        geometric=geom,
        spectral=spec,
        operational_floor_triggered=floor,
        degenerate_O=False,
        resolved_peak_count=len(peaks),
        max_resolved_peaks=max_resolved_peaks,
        peak_count_over_ceiling=over_ceiling,
        chi_abs_floor=chi_abs_floor,
        eta_O=eta_O,
        min_prominence_ratio=min_prominence_ratio,
    )
    return tri, peaks

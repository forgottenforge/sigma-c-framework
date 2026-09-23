# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Resolution / sensitivity band on a scalar sigma_c -- the #1 adoption blocker
The blindness map says sigma_c is "the least trustworthy"
reading and "rides on min_prominence_ratio", but no number ever said HOW MUCH it
can move. This module supplies that number.

It is NOT a confidence interval. There is no sampling model here; a CI would be a
fabricated error bar (feedback_interval_must_carry_bias). Instead it reports the
two ways the REPORTED sigma_c can move with the data held fixed:

  (1) DISCRETIZATION -- half the local grid-cell width (in sigma) at the peak: the
      floor below which sub-grid quadratic refinement cannot locate the peak
      without a model. This sets the value interval [sigma_c - d, sigma_c + d] and,
      coupled to it, the number of significant digits it is honest to print
      (showing 2.000044 when the band is +-0.05 is false precision).

  (2) CONVENTION SENSITIVITY -- the range of min_prominence_ratio over which the
      GEOMETRIC peak-count (hence the I/II/III regime) is unchanged, computed FROM
      THIS observable (the per-input flip points, NOT a hardcoded r=0.40). Tells
      the user: your regime verdict holds for prominence r in [r_lo, r_hi]; outside
      it, a bump crosses the counting threshold and the regime flips.

DISCIPLINE: pure functions of the (immutable) profile + declared convention;
Result.resolution_band() derives them at call time from the retained profile, so
the band can never be a stale stored field (feedback_derived_not_stored_safeguard).
Exercised through the public Result.resolution_band() / to_dict() / summary(), not
tested against these internals directly (the half-fix lesson).
"""
from __future__ import annotations

import math
from typing import Optional

import numpy as np

# Declared convention (scale-free, NOT a theorem): OK requires the peak count to
# be stable across r in [r0/FACTOR, min(1, FACTOR*r0)]. Factor 2 -- an octave of
# the prominence convention each way.
CONVENTION_WINDOW_FACTOR: float = 2.0


def digits_for_band(value: float, band_abs: float) -> tuple:
    """(significant_digits, display_string, rounded_value) for `value` known to
    +-`band_abs`. Digits are coupled to the band: the last kept place is the
    order of magnitude of the band, so nothing below the resolution is printed.
    """
    if (not math.isfinite(value) or not math.isfinite(band_abs)
            or band_abs <= 0.0 or value == 0.0):
        return (4, f"{value:.4g}", float(value))
    order = math.floor(math.log10(band_abs))          # band 0.05 -> -2
    ndigits = -order                                   # decimals to keep
    rounded = round(float(value), ndigits)
    display = f"{rounded:.{max(0, ndigits)}f}" if ndigits > 0 else f"{int(round(rounded))}"
    lead = math.floor(math.log10(abs(value)))
    sig = max(1, lead - order + 1)
    return (sig, display, rounded)


def peak_count_stable_in_window(
    sigma_grid,
    chi,
    min_prominence_ratio: float,
    *,
    chi_abs_floor: float = 0.0,
    factor: float = CONVENTION_WINDOW_FACTOR,
    n_win: int = 33,
) -> bool:
    """Is the resolved peak COUNT unchanged across r in [r0/factor, min(1, factor*r0)]?

    The load-bearing OK gate for ANY regime (scalar or vector): a peak set that
    only holds under a narrow prominence-convention setting is not a measurement.
    Uses the SAME chi_abs_floor as the verdict, so the recount matches classify().
    White noise forced to a stable-looking count by a particular r fails this,
    because lowering/raising r by an octave changes how many noise bumps clear.
    """
    from sigma_c.core.susceptibility import find_interior_maxima
    sigma = np.asarray(sigma_grid, dtype=float)
    chi = np.asarray(chi, dtype=float)
    if sigma.size < 3 or chi.size != sigma.size:
        return True  # nothing to test

    def count(r: float) -> int:
        return len(find_interior_maxima(
            sigma, chi, min_prominence_ratio=r, chi_abs_floor=chi_abs_floor))

    r0 = float(min_prominence_ratio)
    c0 = count(r0)
    lo = r0 / factor
    hi = min(1.0, r0 * factor)
    return all(count(float(r)) == c0 for r in np.linspace(lo, hi, n_win))


def resolution_band(
    sigma_grid,
    chi,
    sigma_c: Optional[float],
    min_prominence_ratio: float,
    *,
    chi_abs_floor: float = 0.0,
    n_sweep: int = 400,
) -> Optional[dict]:
    """Resolution/sensitivity band for a scalar sigma_c. None unless sigma_c is a
    finite scalar and the profile has >= 3 samples."""
    if sigma_c is None or isinstance(sigma_c, (list, tuple)):
        return None
    try:
        sc = float(sigma_c)
    except (TypeError, ValueError):
        return None
    if not math.isfinite(sc):
        return None
    sigma = np.asarray(sigma_grid, dtype=float)
    chi = np.asarray(chi, dtype=float)
    n = sigma.size
    if n < 3 or chi.size != n:
        return None

    # locate the peak index nearest the reported sigma_c
    peak_idx = int(np.argmin(np.abs(sigma - sc)))

    # --- (1) discretization: half the LOCAL grid cell (in sigma) at the peak ---
    if 0 < peak_idx < n - 1:
        cell = max(sigma[peak_idx] - sigma[peak_idx - 1],
                   sigma[peak_idx + 1] - sigma[peak_idx])
    else:
        lo = max(peak_idx - 1, 0)
        hi = min(peak_idx + 1, n - 1)
        cell = abs(sigma[hi] - sigma[lo])
    band_abs = 0.5 * float(cell)
    band_lo = sc - band_abs
    band_hi = sc + band_abs

    # --- (2) convention sensitivity: per-input min_prominence_ratio sweep -------
    from sigma_c.core.susceptibility import find_interior_maxima

    def count(r: float) -> int:
        return len(find_interior_maxima(
            sigma, chi, min_prominence_ratio=r, chi_abs_floor=chi_abs_floor))

    current_count = count(min_prominence_ratio)
    grid = np.linspace(1e-3, 1.0, n_sweep)
    flip_low: Optional[float] = None    # largest r < declared where count differs
    flip_high: Optional[float] = None   # smallest r > declared where count differs
    for r in grid:
        c = count(float(r))
        if c == current_count:
            continue
        if r < min_prominence_ratio:
            flip_low = float(r)         # ascending grid -> ends at the largest
        elif r > min_prominence_ratio and flip_high is None:
            flip_high = float(r)
    stable_lo = flip_low if flip_low is not None else 0.0
    stable_hi = flip_high if flip_high is not None else 1.0

    # --- (2b) load-bearing convention stability for the OK gate -----------------
    # A peak that exists as a SINGLE peak only inside a narrow window of the
    # prominence convention is not a measurement (its Layer-3 location runs away
    # under the convention screw -- Parrot's own logic). OK requires the peak COUNT
    # to be unchanged across the DECLARED window r in [r0/FACTOR, min(1, FACTOR*r0)]
    # (FACTOR = CONVENTION_WINDOW_FACTOR, a declared, scale-free convention). Pure
    # white noise forced to a single peak by a high r fails this at every r0.
    r0 = float(min_prominence_ratio)
    win_lo = r0 / CONVENTION_WINDOW_FACTOR
    win_hi = min(1.0, r0 * CONVENTION_WINDOW_FACTOR)
    stable_in_window = peak_count_stable_in_window(
        sigma, chi, min_prominence_ratio, chi_abs_floor=chi_abs_floor)

    # --- significant digits coupled to the discretization band -----------------
    sig, display, rounded = digits_for_band(sc, band_abs)

    return {
        "sigma_c": sc,
        "is_confidence_interval": False,
        "grid_cell_width": float(cell),
        "band_abs": band_abs,
        "band_lo": band_lo,
        "band_hi": band_hi,
        "significant_digits": sig,
        "sigma_c_display": display,
        "sigma_c_rounded": rounded,
        "min_prominence_ratio": float(min_prominence_ratio),
        "regime_peak_count": current_count,
        "prominence_stable_range": (stable_lo, stable_hi),
        "prominence_flip_low": flip_low,
        "prominence_flip_high": flip_high,
        "convention_window": (win_lo, win_hi),
        "peak_count_stable_in_window": stable_in_window,
    }

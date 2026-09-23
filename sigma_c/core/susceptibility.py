# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
chi_O computation — paper Def 2.2 (def:sigmac), Def 2.1 (def:Onice).

Maps observable samples (sigma_k, O(sigma_k)) to chi_O(sigma) = |sigma * O'(sigma)|
on a regular log-spaced grid, with explicit smoothing provenance
(Enforcement 8 of the design).

The smoothing parameters are *logged*, not solved — the full numerical
robustness analysis is delegated to a companion paper (Rem 4.12 of the
paper).
"""
from __future__ import annotations
from dataclasses import dataclass
from typing import Callable, List, Optional, Tuple

import numpy as np


@dataclass(frozen=True)
class SmoothingSpec:
    """The reproduction recipe carried in every result (Enforcement 8)."""
    kernel: str = "log_linear"     # "log_linear" | "analytical" | "none"
    bandwidth: Optional[float] = None
    interpolant_order: int = 1

    def as_dict(self) -> dict:
        return {
            "kernel": self.kernel,
            "bandwidth": self.bandwidth,
            "interpolant_order": self.interpolant_order,
        }


def chi_O(
    sigma: np.ndarray,
    O: np.ndarray,
    *,
    smoothing: Optional[SmoothingSpec] = None,
) -> Tuple[np.ndarray, np.ndarray, SmoothingSpec]:
    """
    Compute chi_O(sigma) = |sigma * dO/dsigma|.

    Returns (sigma_grid, chi_values, smoothing_used) where smoothing_used
    records what was applied for downstream provenance.

    Default: numerical log-derivative on the provided grid; smoothing="none"
    if input is already analytical (callable).
    """
    from sigma_c.errors import InvalidInputError
    sigma = np.asarray(sigma, dtype=float)
    O = np.asarray(O, dtype=float)
    if sigma.shape != O.shape:
        raise InvalidInputError(
            f"sigma and O must have the same shape (got {sigma.shape} vs {O.shape})."
        )
    if not np.all(np.isfinite(sigma)) or not np.all(np.isfinite(O)):
        raise InvalidInputError(
            "sigma and O must be finite (no NaN/Inf): a non-finite value propagates "
            "through the log-derivative into a spurious OK verdict. Bad input is "
            "rejected, not read as a measurement."
        )
    if len(sigma) < 5:
        raise InvalidInputError("Need at least 5 sample points to compute chi_O.")
    if np.any(sigma <= 0):
        raise InvalidInputError(
            "sigma values must be strictly positive (log-domain): sigma_c lives on "
            "a log-spaced dial, so sigma <= 0 has no meaning. If sigma is a time "
            "axis starting at 0, drop the t=0 sample or shift the axis so all values "
            "are > 0. (And note: a single raw time series is NOT a (dial, observable) "
            "pairing -- see the README limitation on relaxation times.)"
        )

    if smoothing is None:
        smoothing = SmoothingSpec(kernel="log_linear", bandwidth=None,
                                  interpolant_order=1)

    # Sort by sigma if not already
    order = np.argsort(sigma)
    sigma_s = sigma[order]
    O_s = O[order]

    # Numerical d O / d log sigma using central differences on log sigma.
    log_sigma = np.log(sigma_s)
    dO_dlogsigma = np.gradient(O_s, log_sigma, edge_order=2)
    chi = np.abs(dO_dlogsigma)

    return sigma_s, chi, smoothing


def chi_O_from_callable(
    O_fn: Callable[[np.ndarray], np.ndarray],
    sigma_grid: np.ndarray,
) -> Tuple[np.ndarray, np.ndarray, SmoothingSpec]:
    """
    Compute chi_O from a callable observable O: sigma -> O(sigma).

    Uses analytical-precision finite differences. Marked smoothing="analytical".
    """
    sigma_grid = np.asarray(sigma_grid, dtype=float)
    h = sigma_grid * 1e-7  # relative step in log-domain
    O_plus = O_fn(sigma_grid + h)
    O_minus = O_fn(sigma_grid - h)
    dO = (O_plus - O_minus) / (2 * h)
    chi = np.abs(sigma_grid * dO)
    return sigma_grid, chi, SmoothingSpec(kernel="analytical", bandwidth=None,
                                          interpolant_order=2)


def find_interior_maxima(
    sigma: np.ndarray,
    chi: np.ndarray,
    *,
    min_prominence_ratio: float = 0.10,
    chi_abs_floor: float = 0.0,
    boundary_margin: int = 1,
) -> List[int]:
    """
    Indices of strict interior local maxima of chi.

    A peak must clear BOTH thresholds:
      - the RELATIVE prominence convention `min_prominence_ratio` * max(chi), and
      - the ABSOLUTE numerical-noise floor `chi_abs_floor`.

    The relative threshold alone normalises chi to its OWN maximum, so in pure
    rounding noise (a constant observable) it always finds "peaks" -- chi_max is
    then itself float noise. The absolute floor, set by the caller from the
    observable's variation and the log-step (~ eps * max|O| / dlog), rejects peaks
    that are indistinguishable from the log-derivative's rounding error. Cite:
    def:Onice, def:Onice-multi for the multi-peak case.
    """
    n = len(chi)
    if n < 3:
        return []
    peaks: List[int] = []
    chi_max = float(np.max(chi))
    threshold = min_prominence_ratio * chi_max if chi_max > 0 else 0.0
    threshold = max(threshold, float(chi_abs_floor))
    # Scan for interior maxima INCLUDING flat-topped (plateau) peaks: a smooth peak
    # whose true maximum lies between two grid points shows up as two (or more) tied
    # top samples, which a strict `chi[i] > chi[i+1]` test misses entirely (a real
    # peak then reads as regime III). A plateau that rises in and falls out counts
    # once, at its centre.
    i = boundary_margin
    while i <= n - 1 - boundary_margin:
        if chi[i] >= threshold and chi[i] > chi[i - 1]:
            j = i
            while j + 1 <= n - 1 - boundary_margin and chi[j + 1] == chi[i]:
                j += 1
            if j + 1 <= n - 1 and chi[j + 1] < chi[i]:
                peaks.append((i + j) // 2)   # plateau centre (== i for a lone peak)
                i = j + 1
                continue
        i += 1
    return peaks


def quadratic_peak_in_log_sigma(
    sigma_grid: np.ndarray,
    chi: np.ndarray,
    peak_idx: Optional[int] = None,
) -> Tuple[float, bool]:
    """
    Sub-grid refinement of a chi peak via quadratic interpolation in
    log(sigma). Returns (sigma_c_sub, at_grid_boundary).

    Without sub-grid refinement, two analyses on the same grid agree to
    within one grid cell trivially -- a cross-detector or cross-window
    `delta = 0.0` is then grid quantisation, not physics. This routine
    promotes the reported sigma_c to a sub-grid value when the peak is
    interior, and surfaces a boundary flag when it is at the grid edge
    (where no quadratic interpolation is defined). Cite: paper Rem 4.12
    (smoothing/numerical robustness).

    Parameters
    ----------
    sigma_grid : (N,) array
        Strictly positive, log-spaced grid.
    chi : (N,) array
        chi values on the grid.
    peak_idx : int, optional
        Index of the peak to refine. If None, taken to be argmax(chi).

    Returns
    -------
    sigma_c_sub : float
        The sub-grid peak location. Equals sigma_grid[peak_idx] when
        refinement is not possible (boundary or degenerate parabola).
    at_grid_boundary : bool
        True iff peak_idx is at index 0 or len(chi)-1 (cannot refine).
    """
    sigma_grid = np.asarray(sigma_grid, dtype=float)
    chi = np.asarray(chi, dtype=float)
    if peak_idx is None:
        peak_idx = int(np.argmax(chi))
    n = len(chi)
    if peak_idx <= 0 or peak_idx >= n - 1:
        return float(sigma_grid[peak_idx]), True

    log_sigma = np.log(sigma_grid)
    h = log_sigma[peak_idx + 1] - log_sigma[peak_idx]
    y0 = float(chi[peak_idx - 1])
    y1 = float(chi[peak_idx])
    y2 = float(chi[peak_idx + 1])
    denom = y0 - 2.0 * y1 + y2
    if abs(denom) < 1e-30:
        # Degenerate parabola; report the grid sample and don't fail.
        return float(sigma_grid[peak_idx]), False
    offset = 0.5 * (y0 - y2) / denom
    # Clamp to [-1, 1] -- a true peak should always have offset in (-1, 1)
    # and runaway values indicate a non-peak (e.g. a noisy plateau).
    offset = max(-1.0, min(1.0, offset))
    sigma_c_sub = float(np.exp(log_sigma[peak_idx] + offset * h))
    return sigma_c_sub, False

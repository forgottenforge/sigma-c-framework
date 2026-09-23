# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Stability indicator gamma_O — paper Prop 4.4 (prop:stability).

gamma_O = -d^2/dsigma^2 chi_O^2 at sigma = sigma_c is the strict-SOC
constant; it bounds |delta sigma_c| <= C * ||delta O||_{C^1} / gamma_O.

Low gamma_O ⟹ flat peak ⟹ noisy/transition-zone reading.
"""
from __future__ import annotations
from typing import Optional

import numpy as np


def compute_gamma_O(
    sigma: np.ndarray,
    chi: np.ndarray,
    peak_idx: int,
) -> Optional[float]:
    """
    Compute gamma_O = -d^2/dsigma^2 chi_O^2 at sigma[peak_idx].

    Returns None if the peak is too close to the boundary to compute SOC.
    Returns gamma_O in units of [chi^2] / [sigma]^2.

    Normalised: returns gamma_O * sigma_peak^2 / chi_peak^2, dimensionless,
    so values across problems are comparable.
    """
    n = len(sigma)
    if peak_idx <= 1 or peak_idx >= n - 2:
        return None
    sig0 = sigma[peak_idx]
    chi0 = chi[peak_idx]
    if chi0 <= 0:
        return None
    # Use chi^2 as paper uses (sigma O')^2 = chi^2
    chi2 = chi ** 2
    # Three-point second derivative on local grid
    s_l, s_c, s_r = sigma[peak_idx - 1], sigma[peak_idx], sigma[peak_idx + 1]
    f_l, f_c, f_r = chi2[peak_idx - 1], chi2[peak_idx], chi2[peak_idx + 1]
    # Non-uniform central difference for d^2 f / d sigma^2
    h_l = s_c - s_l
    h_r = s_r - s_c
    if h_l <= 0 or h_r <= 0:
        return None
    d2 = (
        2 * (h_r * f_l - (h_l + h_r) * f_c + h_l * f_r)
        / (h_l * h_r * (h_l + h_r))
    )
    gamma = -d2
    if gamma <= 0:
        return 0.0  # numerical noise at flat peak
    # Normalise by chi_peak^2 / sigma_peak^2 so the indicator is dimensionless
    return float(gamma * sig0 ** 2 / chi0 ** 2)

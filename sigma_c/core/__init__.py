# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""sigma_c core -- chi_O, trichotomy, stability, faithfulness."""
from sigma_c.core.susceptibility import chi_O, find_interior_maxima
from sigma_c.core.trichotomy import classify, geometric_trichotomy
from sigma_c.core.stability import compute_gamma_O
from sigma_c.core.faithfulness import (
    check_F1,
    check_F2,
    check_F3,
    kl_modal_coefficients,
    FaithfulnessCheck,
)

__all__ = [
    "chi_O",
    "find_interior_maxima",
    "classify",
    "geometric_trichotomy",
    "compute_gamma_O",
    "check_F1",
    "check_F2",
    "check_F3",
    "kl_modal_coefficients",
    "FaithfulnessCheck",
]

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Framework taxonomy — the six standard operator settings of paper Prop 5.4.

Cite: prop:standard-frameworks (THEOREM_MAP).

The declaration determines which kind of spectral identification the bridge
can deliver. None → restrict to dominant-scale-probe reading (Principle 2 of
the design).
"""
from __future__ import annotations
from enum import Enum

from sigma_c.theorem_map import cite


class Framework(str, Enum):
    """The six paper-verified transfer-operator settings (Prop 5.4)."""

    REVERSIBLE_MARKOV = "reversible_markov"
    DOEBLIN = "doeblin"
    COMPACT_HILBERT = "compact_hilbert"
    GNS_LINDBLAD = "gns_lindblad"
    TRANSFER_MATRIX_1D = "transfer_matrix_1d"
    # Prop 5.4 case (6): functional-space-dependent. Shipped experimental.
    ANISOTROPIC_BANACH = "anisotropic_banach"

    @property
    def is_experimental(self) -> bool:
        """Anisotropic-Banach requires user-supplied spectral hypotheses."""
        return self is Framework.ANISOTROPIC_BANACH

    @property
    def reading_kind(self) -> str:
        """
        Whether sigma_c/rho_star is a literal spectral-gap reading
        (cite: prop:probe, Caveat 2 of §7.2) or the operational
        dominant-scale reading.
        """
        return "spectral_gap" if not self.is_experimental else "spectral_gap_experimental"

    def cite_paper(self) -> str:
        return cite("prop:standard-frameworks", note=f"case for {self.value}")

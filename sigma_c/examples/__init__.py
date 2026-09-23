# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""Shared helpers for the shipped examples."""


def save_card(result, save_to):
    """Save a Result card, degrading gracefully when matplotlib (the ``[plot]``
    extra) is not installed — so every demo runs headless, printing its analysis
    and simply skipping the PNG. Returns True if a card was written."""
    try:
        result.card(save_to=save_to)
        print(f"card saved: {save_to}")
        return True
    except ImportError:
        print(f"(skipped card {save_to} — install \"sigma-c-framework[plot]\" for the PNG)")
        return False

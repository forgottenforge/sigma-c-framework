# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Typed exceptions, so a pipeline can branch on error *kind* without string-matching.

`InvalidInputError` subclasses the builtin `ValueError` on purpose: every guard that
used to raise a bare `ValueError` now raises this, so existing
`except ValueError:` code keeps working unchanged, while new code can catch the
specific type. Requested by a headless-integrator usability run.
"""
from __future__ import annotations


class SigmaCError(Exception):
    """Base class for all sigma_c-specific errors."""


class InvalidInputError(SigmaCError, ValueError):
    """Input to the kernel is malformed or degenerate: non-finite values, mismatched
    shapes, too few samples, non-positive dial values, or a contradictory tau-route
    declaration. Subclasses `ValueError` for backward compatibility."""

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Window family registry — paper §4.2. The ρ_⋆ identity follows from
prop:structural-reduction (σ_c = ρ_⋆·τ), not from thm:compat: thm:compat is the
four-axiom compatibility theorem (A1,A2,A4,A5) and has nothing to do with the
window ρ_⋆ values (THEOREM_MAP honesty note).

Each window carries its analytical rho_star, computed from the Mellin
transform, NOT from a fit. The rho_star_source field on the Result object
records exactly which window produced the value (Enforcement item 2).

Five canonical windows ship as first-class primitives:

| name         | w(u)                          | rho_star          |
|--------------|-------------------------------|-------------------|
| bare         | delta at u=1                  | 1                 |
| gamma2       | u * exp(-u)                   | 2 - sqrt(3)       |
| gamma3       | u**2 * exp(-u)                | 3 - 2*sqrt(2)     |
| exponential  | exp(-u)                       | 1                 |
| log_gaussian | (2*pi)^(-1/2) * exp(-(log u)^2/2)  | 1            |
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Callable, Dict, Optional, Union
import math
import numpy as np


@dataclass(frozen=True)
class Window:
    """
    A canonical analytical window with its proved rho_star.

    Attributes
    ----------
    name : str
        Identifier — also forms `rho_star_source = "analytic:<name>"`
        on the result object.
    rho_star : float
        Analytical profile constant sigma_c[g] for g(u) = u * w_mellin(u).
        Computed from the Mellin transform; never fitted.
    description : str
        Human-readable form of the window function (used on cards).
    weight : Callable[[np.ndarray], np.ndarray] | None
        Per-sample weight applied to the observable when smoothing.
        None means "no smoothing" (bare correlator path).
    """
    name: str
    rho_star: float
    description: str
    weight: Optional[Callable] = None

    def __repr__(self) -> str:
        return f"Window({self.name!r}, rho_star={self.rho_star:.4g})"

    @property
    def rho_star_source(self) -> str:
        """The provenance string carried in every result object."""
        return f"analytic:{self.name}"


# ---------------------------------------------------------------------------
# Canonical window constructors
# ---------------------------------------------------------------------------

def bare() -> Window:
    """
    Bare correlator: w = delta at u=1, rho_star = 1.

    The textbook reading: sigma_c equals tau directly with no profile factor.
    """
    return Window(
        name="bare",
        rho_star=1.0,
        description="delta(u-1) (no smoothing)",
        weight=None,
    )


def gamma_k(k: int = 2) -> Window:
    """
    Gamma-k window: w(u) = u^(k-1) * exp(-u).

    Analytical rho_star = k - sqrt(k**2 - 1).
    Cite: paper §4.2 table.
    """
    if k < 2:
        raise ValueError("gamma_k requires k >= 2 (k=1 collapses to exponential).")
    rho_star = float(k - math.sqrt(k * k - 1))

    def weight(u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        return (u ** (k - 1)) * np.exp(-u)

    return Window(
        name=f"gamma{k}",
        rho_star=rho_star,
        description=f"u^{k - 1} * exp(-u)  (Gamma-{k})",
        weight=weight,
    )


def exponential() -> Window:
    """w(u) = exp(-u), rho_star = 1."""
    def weight(u: np.ndarray) -> np.ndarray:
        return np.exp(-np.asarray(u, dtype=float))

    return Window(
        name="exponential",
        rho_star=1.0,
        description="exp(-u)",
        weight=weight,
    )


def log_gaussian() -> Window:
    """
    Log-Gaussian: w(u) = (2*pi)^(-1/2) * exp(-(log u)^2 / 2), rho_star = 1.

    Special property: log-scale-invariant shape → rho_star = 1 at all widths.
    """
    def weight(u: np.ndarray) -> np.ndarray:
        u = np.asarray(u, dtype=float)
        out = np.zeros_like(u)
        mask = u > 0
        out[mask] = math.exp(0) / math.sqrt(2 * math.pi) * np.exp(
            -np.log(u[mask]) ** 2 / 2
        )
        return out

    return Window(
        name="log_gaussian",
        rho_star=1.0,
        description="(2*pi)^(-1/2) * exp(-(log u)^2 / 2)",
        weight=weight,
    )


# ---------------------------------------------------------------------------
# Registry
# ---------------------------------------------------------------------------

WINDOW_REGISTRY: Dict[str, Callable[[], Window]] = {
    "bare": bare,
    "gamma2": lambda: gamma_k(2),
    "gamma3": lambda: gamma_k(3),
    "exponential": exponential,
    "log_gaussian": log_gaussian,
}


def resolve(window: Union[str, Window, None]) -> Window:
    """Look up a window by name, pass through if already a Window, default to bare."""
    if window is None:
        return bare()
    if isinstance(window, Window):
        return window
    if window in WINDOW_REGISTRY:
        return WINDOW_REGISTRY[window]()
    raise ValueError(
        f"Unknown window {window!r}. "
        f"Available: {sorted(WINDOW_REGISTRY)}; or pass a Window instance."
    )


__all__ = [
    "Window",
    "WINDOW_REGISTRY",
    "bare",
    "gamma_k",
    "exponential",
    "log_gaussian",
    "resolve",
]

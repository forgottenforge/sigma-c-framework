# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
tau from a DECLARED single-mode profile — the honest yes-path for the most common
user (a decay/relaxation curve who wants tau).

This is NOT a new theorem and NOT the certified spectral abscissa (that needs an
operator). It is the profile decomposition with a CHECKABLE DECLARATION:

  1. You DECLARE the profile you believe your observable follows (a single
     exponential relaxation, or a Debye loss peak). Declaring it is one of the
     five decision points; it is recorded in the result.
  2. The tool computes chi_O = |sigma dO/dsigma| and compares it to the analytic
     single-mode susceptibility TEMPLATE for your declared profile, centred at the
     data's own sigma_c. It reports the relative residual.
  3. If the residual is small (the single-mode declaration holds), Proposition 2
     makes the bridge exact, so tau = sigma_c / rho_star is returned CONDITIONAL
     ON the declared profile, with the residual attached.
  4. If the residual is large, the declared profile does NOT fit (your data is
     multi-mode or a different shape): tau is refused, not guessed.

The residual check is what catches a multi-mode profile; its reliability is
validated in test_profile_tau.py (single mode accepted with the right tau; a
two-mode curve rejected across a sweep of mode separations and amplitudes).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Sequence, Union

import numpy as np

# Declared conventions (NOT theorems), surfaced in the result.
DEFAULT_MAX_RESIDUAL: float = 0.10   # relative RMS residual to accept the profile

_PROFILES = ("exponential", "debye")


def _template(u: np.ndarray, profile: str) -> np.ndarray:
    """Analytic single-mode chi_O template as a function of u = sigma / sigma_c,
    normalised to peak value 1 at u = 1."""
    u = np.asarray(u, dtype=float)
    if profile == "exponential":
        # O = exp(-sigma/tau)  ->  chi = (sigma/tau) exp(-sigma/tau);  peak at u=1.
        return u * np.exp(1.0 - u)
    if profile == "debye":
        # Debye storage observable O = 1/(1+(sigma/tau)^2); its susceptibility
        # chi = |sigma dO/dsigma| = 4u^2/(1+u^2)^2, peak 1 at u=1 (sigma_c = tau).
        return 4.0 * u * u / (1.0 + u * u) ** 2
    raise ValueError(f"profile must be one of {_PROFILES}; got {profile!r}.")


def _relative_residual(chi: np.ndarray, template: np.ndarray) -> float:
    """Relative RMS residual between the observed chi and the template, each scaled
    to unit peak. Scale-free in the observable amplitude."""
    a = np.asarray(chi, dtype=float)
    b = np.asarray(template, dtype=float)
    amax = float(np.max(a))
    bmax = float(np.max(b))
    if amax <= 0 or bmax <= 0:
        return float("inf")
    a = a / amax
    b = b / bmax
    return float(np.sqrt(np.mean((a - b) ** 2)) / np.sqrt(np.mean(b ** 2)))


@dataclass
class ProfileTauResult:
    """tau conditional on a declared single-mode profile. Never a bare number:
    carries the code, the residual, and the declared profile."""
    tau: Optional[float]
    sigma_c: Optional[float]
    profile: str
    residual: float
    max_residual: float
    code: str
    reason: str
    window: str = "bare"
    rho_star: Optional[float] = None
    _extra: dict = field(default_factory=dict)

    @property
    def ok(self) -> bool:
        return self.code == "OK_CONDITIONAL"

    def summary(self) -> str:
        head = f"VERDICT : {self.code}"
        if self.tau is not None:
            head += f"  --  tau ~ {self.tau:.4g}  (conditional on declared '{self.profile}' profile)"
        else:
            head += f"  --  tau = _|_  (declared '{self.profile}' profile)"
        lines = [head,
                 f"fit     : relative residual {self.residual:.3g} "
                 f"(accept <= {self.max_residual:g})",
                 f"note    : {self.reason}",
                 f"(window={self.window}; this is CONDITIONAL on your declared profile, "
                 f"NOT the certified spectral abscissa)"]
        return "\n".join(lines)

    def to_dict(self) -> dict:
        from sigma_c import __version__
        return {
            "schema_version": 1,
            "library_version": __version__,
            "quantity": "tau_profile_conditional",
            "tau": self.tau,
            "sigma_c": self.sigma_c,
            "declared_profile": self.profile,
            "relative_residual": self.residual,
            "max_residual": self.max_residual,
            "window": self.window,
            "rho_star": self.rho_star,
            "code": self.code,
            "reason": self.reason,
            "certified": False,
            **self._extra,
        }


def tau_from_profile(
    sigma: Union[np.ndarray, Sequence[float]],
    O: Union[np.ndarray, Sequence[float]],
    *,
    profile: str = "exponential",
    window: Union[str, object, None] = "bare",
    max_residual: float = DEFAULT_MAX_RESIDUAL,
    min_prominence_ratio: float = 0.10,
) -> ProfileTauResult:
    """tau for a DECLARED single-mode profile (see module docstring).

    Parameters
    ----------
    sigma, O : arrays
        The dial and the observable, as for analyze().
    profile : {"exponential", "debye"}
        The single-mode shape you DECLARE your observable follows.
    window : window name or Window
        Sets rho_star (recorded). Default "bare" (rho_star = 1, tau = sigma_c).
    max_residual : float
        Declared acceptance threshold on the relative RMS residual.

    Returns a ProfileTauResult (never a bare number): OK_CONDITIONAL with tau and
    the residual, or a refusal (NOT_RESOLVABLE / NOT_APPLICABLE) with the residual.
    """
    if profile not in _PROFILES:
        raise ValueError(f"profile must be one of {_PROFILES}; got {profile!r}.")
    from sigma_c import analyze
    res = analyze(sigma, O, window=window, min_prominence_ratio=min_prominence_ratio)
    win_name = res.window or "bare"
    rho = res.rho_star

    # A single declared mode needs a single resolved peak.
    if res.sigma_c is None:
        return ProfileTauResult(
            tau=None, sigma_c=None, profile=profile, residual=float("inf"),
            max_residual=max_residual, code="NOT_RESOLVABLE", window=win_name,
            rho_star=rho,
            reason=("no interior susceptibility peak (regime III): there is no single "
                    "mode to attach a declared profile to."),
        )
    if isinstance(res.sigma_c, list):
        return ProfileTauResult(
            tau=None, sigma_c=None, profile=profile, residual=float("inf"),
            max_residual=max_residual, code="NOT_RESOLVABLE", window=win_name,
            rho_star=rho,
            reason=(f"multiple resolved peaks ({len(res.sigma_c)}): the observable is "
                    "not single-mode, so a single declared-profile tau does not apply."),
        )

    sc = float(res.sigma_c)
    sg = np.asarray(res._profile_sigma, dtype=float)
    chi = np.asarray(res._profile_chi, dtype=float)
    template = _template(sg / sc, profile)
    residual = _relative_residual(chi, template)

    if residual <= max_residual:
        tau = sc / rho if (rho is not None and rho > 0) else sc
        return ProfileTauResult(
            tau=tau, sigma_c=sc, profile=profile, residual=residual,
            max_residual=max_residual, code="OK_CONDITIONAL", window=win_name,
            rho_star=rho,
            reason=(f"the chi_O profile matches the declared single-mode "
                    f"'{profile}' template (residual {residual:.3g} <= {max_residual:g}); "
                    f"by Prop 2 the bridge is exact for a single mode, so "
                    f"tau = sigma_c/rho_star is valid CONDITIONAL on this declaration. "
                    f"It is NOT the certified spectral abscissa (which needs an operator)."),
        )
    return ProfileTauResult(
        tau=None, sigma_c=sc, profile=profile, residual=residual,
        max_residual=max_residual, code="NOT_APPLICABLE", window=win_name,
        rho_star=rho,
        reason=(f"the declared single-mode '{profile}' profile does NOT fit the chi_O "
                f"shape (residual {residual:.3g} > {max_residual:g}): your data is "
                f"multi-mode or a different shape, so a single-mode tau would be "
                f"misleading. Declare a different profile, or use the operator route "
                f"for a certified tau."),
    )

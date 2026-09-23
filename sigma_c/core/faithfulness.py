# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Faithfulness sufficient conditions F1/F2/F3 -- paper Prop 5.10
(prop:faith-sufficient).

These are the *explicit* checkers that close the design's Item B:
This module makes the three paper-specified
pathways into typed checks that each return an explicit residual constant
C_R, never inferred from a fit.

Cite (in docstrings):
    prop:faith-sufficient -- the three sufficient conditions
    def:faithful           -- the modal-sum form being verified
    thm:genericity         -- the F1 backbone (App D.1)
    prop:standard-frameworks (case 2) -- the Doeblin framework used by F3
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple, Union

import math
import numpy as np


# ---------------------------------------------------------------------------
# Common result type for all three checkers
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class FaithfulnessCheck:
    """
    Outcome of a single F1/F2/F3 check.

    Attributes
    ----------
    passed : bool
        True iff the sufficient condition holds with the supplied inputs.
    condition : str
        One of "F1", "F2", "F3" -- which paper sub-condition was checked.
    C_R : float
        Explicit residual constant. The paper's bound on the modal-sum
        residual: ||R||_inf <= C_R * A * r^n, where A = amplitude.
        C_R = +inf when the check fails (the modal-sum form is not
        certifiable from the given inputs).
    faithfulness_order : Union[int, float]
        The order n of the residual decay r^n. For F1 it is the window-range
        bound M; for F2 it is +inf (spectrally filtered); for F3 it is the
        derived geometric-decay exponent.
    details : dict
        Numerical evidence the user can inspect: modal ratios, decay
        factors, etc. Always populated, never opaque.
    citation : str
        The paper label this check is anchored to (always "prop:faith-sufficient").
    """
    passed: bool
    condition: str
    C_R: float
    faithfulness_order: Union[int, float]
    details: dict = field(default_factory=dict)
    citation: str = "prop:faith-sufficient"

    def summary(self) -> str:
        from sigma_c.theorem_map import cite

        status = "PASS" if self.passed else "FAIL"
        cr = "inf" if math.isinf(self.C_R) else f"{self.C_R:.4g}"
        n = ("inf" if math.isinf(self.faithfulness_order)
             else f"{self.faithfulness_order}")
        lines = [
            f"Faithfulness check {self.condition}: {status}",
            f"  C_R                  : {cr}",
            f"  faithfulness order n : {n}",
            f"  backed by            : {cite(self.citation)} (sub-condition {self.condition})",
        ]
        if self.details:
            lines.append("  details:")
            for k, v in self.details.items():
                if isinstance(v, float):
                    lines.append(f"    {k:<22}: {v:.4g}")
                else:
                    lines.append(f"    {k:<22}: {v}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# F1: Modal-coefficient bound with prescribed window range
# ---------------------------------------------------------------------------

def check_F1(
    modal_coeffs: Sequence[complex],
    *,
    gap_ratio: float,
    window_range_M: float = 1.0,
    window_mellin_bound: float = 1.0,
    shape_norm_bound: float = 1.0,
    c_min: Optional[float] = None,
    n_modes_total: Optional[int] = None,
) -> FaithfulnessCheck:
    """
    Check sufficient condition (F1) of paper Prop 5.10.

    The probe Phi has KL modal coefficients c_j(Phi) for j >= 2. (F1) holds
    when:
      |c_2(Phi)| >= c_min > 0                   (probe couples to lambda_2)
      sup_{j >= 3} |c_j(Phi)| <= C_Phi          (sub-leading mass bounded)
      window range sigma_max / T_* <= M         (bounded window)
      Schwartz Mellin: ||w_tilde||_inf bounded
    Then O = O[Phi, w] is in O_faith(x; lambda_2, M) with
      C_R <= kappa(X) * (C_Phi / c_min) * ||w_tilde||_inf / ||g||_inf

    Parameters
    ----------
    modal_coeffs : sequence of complex
        Modal coefficients [c_2, c_3, c_4, ...] of the probe in the KL
        expansion w.r.t. the transfer operator. Index 0 of the sequence
        is c_2 (the leading non-trivial mode).
    gap_ratio : float
        r = sup_{j >= 3} |lambda_j| / |lambda_2|, the spectral gap ratio
        of the underlying transfer operator. The condition is satisfied
        only if r < 1 AND the modal bound below holds.
    window_range_M : float, default 1.0
        Bound M >= sigma_max / T_* on the window range. Smaller M ->
        stronger faithfulness order (n grows with M in the bookkeeping;
        residual decays as r^M).
    window_mellin_bound : float, default 1.0
        Bound on ||w_tilde||_inf (the Mellin-transform amplitude of the
        window). For the canonical analytical windows of paper Section 4.2
        this is O(1).
    shape_norm_bound : float, default 1.0
        Bound on ||g||_inf, the L^inf norm of the profile shape
        g(u) = u * w_tilde(u). Used as the denominator in the C_R formula.
    c_min : float, optional
        Lower bound on |c_2(Phi)|. If None, taken to be |c_2| from the
        supplied modal_coeffs.
    n_modes_total : int, optional
        Bound kappa(X) on the number of sub-leading modes contributing. If
        None, taken to be len(modal_coeffs) - 1.

    Returns
    -------
    FaithfulnessCheck
        Cite: prop:faith-sufficient (F1), thm:genericity for the App D
        backbone.
    """
    coeffs = np.asarray(modal_coeffs, dtype=complex)
    if coeffs.size < 1:
        raise ValueError("modal_coeffs must contain at least c_2.")
    if not (0.0 <= gap_ratio < 1.0):
        return FaithfulnessCheck(
            passed=False, condition="F1",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "gap_ratio not in [0, 1)",
                     "gap_ratio_received": gap_ratio},
        )
    if shape_norm_bound <= 0:
        return FaithfulnessCheck(
            passed=False, condition="F1",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "shape_norm_bound must be > 0"},
        )
    if window_range_M < 1:
        return FaithfulnessCheck(
            passed=False, condition="F1",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "window_range_M must be >= 1"},
        )

    c2 = complex(coeffs[0])
    abs_c2 = abs(c2)
    if c_min is None:
        c_min = abs_c2
    if abs_c2 < c_min or abs_c2 == 0:
        return FaithfulnessCheck(
            passed=False, condition="F1",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "|c_2| below c_min (probe does not couple to lambda_2)",
                     "abs_c_2": abs_c2, "c_min": c_min},
        )

    sub = np.abs(coeffs[1:])
    C_phi = float(np.max(sub)) if sub.size > 0 else 0.0
    sub_ratios = (sub / abs_c2).tolist() if sub.size > 0 else []
    max_sub_ratio = max(sub_ratios) if sub_ratios else 0.0

    if n_modes_total is None:
        n_modes_total = max(1, coeffs.size - 1)
    kappa = float(n_modes_total)

    # paper's explicit formula
    C_R = kappa * (C_phi / c_min) * window_mellin_bound / shape_norm_bound
    # The faithfulness order in the bookkeeping is M.
    n_order = int(max(1, math.floor(window_range_M)))

    # The bookkeeping passes when the per-mode ratio is at or below r.
    # Stricter than the paper requires (which only uses the worst-case rho
    # bound), but matches what an honest check_F1 reports.
    passed = max_sub_ratio <= gap_ratio + 1e-12

    return FaithfulnessCheck(
        passed=passed,
        condition="F1",
        C_R=C_R if passed else math.inf,
        faithfulness_order=n_order,
        details={
            "abs_c_2": abs_c2,
            "C_Phi": C_phi,
            "max_sub_ratio": max_sub_ratio,
            "gap_ratio_r": gap_ratio,
            "window_range_M": window_range_M,
            "kappa": kappa,
            "shape_norm_bound": shape_norm_bound,
            "window_mellin_bound": window_mellin_bound,
            "sub_ratios": sub_ratios,
        },
    )


# ---------------------------------------------------------------------------
# F2: Spectrally filtered window
# ---------------------------------------------------------------------------

def check_F2(
    *,
    spectrum: Sequence[complex],
    window_mellin_support: Tuple[float, float],
    filter_tolerance: float = 1e-12,
) -> FaithfulnessCheck:
    """
    Check sufficient condition (F2) of paper Prop 5.10.

    The window w is spectrally filtered if its Mellin transform w_tilde
    vanishes outside a band that excludes all sub-leading modes. With
    lambda_1 >= |lambda_2| > |lambda_3| >= ... > 0, set
        u_2 = -log|lambda_2 / lambda_1|       (== T_*/tau)
        u_3 = -log|lambda_3 / lambda_1|
    F2 holds if w_tilde's Mellin support is contained in (a, b) with
        u_2 in (a, b),
        u_3 NOT in (a, b) (i.e. b < u_3).

    Then the modal expansion of O = O[Phi, w] collapses to its lambda_2
    component up to the numerical filter tolerance: R == 0 to tolerance,
    C_R = O(filter_tolerance), n = infinity in the bookkeeping.

    Parameters
    ----------
    spectrum : sequence of complex
        Full spectrum {lambda_1, lambda_2, ...} in decreasing |lambda| order.
    window_mellin_support : (float, float)
        (a, b), the closed support of w_tilde on the Mellin axis (in units
        of -log|lambda/lambda_1|).
    filter_tolerance : float, default 1e-12
        Numerical floor used as the achieved C_R when the check passes.

    Returns
    -------
    FaithfulnessCheck
        Cite: prop:faith-sufficient (F2).
    """
    spec = sorted((abs(complex(x)) for x in spectrum), reverse=True)
    if len(spec) < 3:
        return FaithfulnessCheck(
            passed=False, condition="F2",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "need at least 3 eigenvalues (lambda_1, 2, 3) to check filter band"},
        )
    lam1 = spec[0]
    lam2 = spec[1]
    lam3 = spec[2]
    if lam1 == 0 or lam2 == 0:
        return FaithfulnessCheck(
            passed=False, condition="F2",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "leading eigenvalues zero"},
        )
    u_2 = -math.log(lam2 / lam1)
    u_3 = -math.log(lam3 / lam1) if lam3 > 0 else math.inf
    a, b = float(window_mellin_support[0]), float(window_mellin_support[1])

    contains_u2 = (a <= u_2 <= b)
    excludes_u3 = (b < u_3)

    passed = contains_u2 and excludes_u3

    return FaithfulnessCheck(
        passed=passed,
        condition="F2",
        C_R=float(filter_tolerance) if passed else math.inf,
        faithfulness_order=math.inf if passed else 0,
        details={
            "u_2": u_2,
            "u_3": u_3,
            "window_support": [a, b],
            "contains_u2": contains_u2,
            "excludes_u3": excludes_u3,
        },
    )


# ---------------------------------------------------------------------------
# F3: Doeblin-minorised system with bounded-variation probe
# ---------------------------------------------------------------------------

def check_F3(
    *,
    doeblin_epsilon: float,
    n_0: int,
    probe_tv: float,
    abs_c_2: Optional[float] = None,
    n_modes: int = 8,
    window_mellin_bound: float = 1.0,
    shape_norm_bound: float = 1.0,
) -> FaithfulnessCheck:
    """
    Check sufficient condition (F3) of paper Prop 5.10.

    Under Doeblin minorisation (paper Prop 5.4 case (2)) at level eps^*
    with n_0-step minorisation, |lambda_2| <= (1 - eps^*)^(1/n_0). If the
    probe Phi has bounded variation TV(Phi), then
        |c_j(Phi)| / |c_2(Phi)| <= TV(Phi)/|c_2(Phi)| * (1 - eps^*)^((j-2)/n_0)
    for j >= 3, i.e. modal coefficients decay geometrically.

    The constant C_R is bounded by the explicit summable series:
        C_R <= sum_{j>=3} (1 - eps^*)^((j-2)/n_0) * (TV(Phi)/|c_2|)
             * ||w_tilde||_inf / ||g||_inf

    Parameters
    ----------
    doeblin_epsilon : float
        Doeblin minorization level eps^* in (0, 1). The chain admits a
        coupling with probability eps^* every n_0 steps.
    n_0 : int
        Minorization step count (>= 1).
    probe_tv : float
        Total variation TV(Phi) >= 0 of the probe.
    abs_c_2 : float, optional
        Lower bound on |c_2(Phi)|; required to anchor the modal ratios.
    n_modes : int, default 8
        Number of sub-leading modes to bound; the geometric tail is added
        in closed form beyond this.
    window_mellin_bound, shape_norm_bound : float
        As in F1.

    Returns
    -------
    FaithfulnessCheck
        Cite: prop:faith-sufficient (F3), prop:standard-frameworks (case 2).
    """
    if not (0 < doeblin_epsilon < 1):
        return FaithfulnessCheck(
            passed=False, condition="F3",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "doeblin_epsilon must be in (0, 1)"},
        )
    if n_0 < 1:
        return FaithfulnessCheck(
            passed=False, condition="F3",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "n_0 must be >= 1"},
        )
    if probe_tv < 0:
        return FaithfulnessCheck(
            passed=False, condition="F3",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "probe_tv must be >= 0"},
        )
    if abs_c_2 is None or abs_c_2 <= 0:
        return FaithfulnessCheck(
            passed=False, condition="F3",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "abs_c_2 must be > 0"},
        )

    decay = (1.0 - doeblin_epsilon) ** (1.0 / n_0)
    # geometric tail: sum_{j=3}^infty decay^(j-2) = decay / (1 - decay)
    if decay >= 1.0:
        return FaithfulnessCheck(
            passed=False, condition="F3",
            C_R=math.inf, faithfulness_order=0,
            details={"reason": "geometric decay factor must be < 1"},
        )
    series_sum = decay / (1.0 - decay)
    C_R = (probe_tv / abs_c_2) * series_sum * (window_mellin_bound / shape_norm_bound)

    # The faithfulness order n is the floor of -log(filter_tol)/log(decay) but
    # for the bookkeeping we report the gap-implied bound directly.
    # NOTE (honesty): the 0.1 here is a DECLARED filter-tolerance
    # convention (10%), not a paper-fixed constant; it sets what counts as
    # "negligible" in the order estimate.
    n_order = max(1, int(math.floor(-math.log(0.1) / -math.log(decay))))

    # The check passes when the series converges and gives a meaningful bound.
    passed = C_R < math.inf and series_sum < math.inf

    return FaithfulnessCheck(
        passed=passed,
        condition="F3",
        C_R=C_R if passed else math.inf,
        faithfulness_order=n_order,
        details={
            "doeblin_epsilon": doeblin_epsilon,
            "n_0": n_0,
            "probe_tv": probe_tv,
            "abs_c_2": abs_c_2,
            "decay_factor": decay,
            "geometric_series_sum": series_sum,
            "window_mellin_bound": window_mellin_bound,
            "shape_norm_bound": shape_norm_bound,
        },
    )


# ---------------------------------------------------------------------------
# Utility: KL modal coefficients for a finite-state Markov chain
# ---------------------------------------------------------------------------

def kl_modal_coefficients(
    P: np.ndarray,
    probe: np.ndarray,
    initial: Optional[np.ndarray] = None,
) -> Tuple[np.ndarray, np.ndarray]:
    """
    Compute the KL modal coefficients c_j = <nu, psi_j> <mu_j, Phi> / <mu_j, psi_j>
    for a finite-state Markov chain (paper App B, Lemma B.3).

    Returns
    -------
    eigenvalues : np.ndarray
        Sorted by descending |lambda|.
    coeffs : np.ndarray
        c_j matched to eigenvalues in the same order. c[0] corresponds to
        lambda_1 (typically the invariant mode); c[1] = c_2 etc.

    This is the helper users need to drive `check_F1` from a Markov P matrix
    rather than handing in coefficients by hand. Cite:
    def:faithful, lem:O-admissible (App B.3 in the paper).
    """
    P = np.asarray(P, dtype=float)
    n = P.shape[0]
    if P.shape != (n, n):
        raise ValueError("P must be a square matrix.")
    if initial is None:
        initial = np.ones(n) / n
    initial = np.asarray(initial, dtype=float)
    if initial.shape != (n,):
        raise ValueError("initial distribution shape must match P.")
    probe = np.asarray(probe, dtype=float)
    if probe.shape != (n,):
        raise ValueError("probe shape must match P.")

    # Right eigenvectors psi_j of P, left eigenvectors mu_j of P
    eigvals_r, psi = np.linalg.eig(P)
    eigvals_l, mu = np.linalg.eig(P.T)

    # Sort right eigensystem by |lambda| descending
    idx_r = np.argsort(-np.abs(eigvals_r))
    eigvals_r = eigvals_r[idx_r]
    psi = psi[:, idx_r]

    # Match left eigenvectors to right by eigenvalue (greedy matching)
    eigvals_l_sorted = np.zeros_like(eigvals_r)
    mu_matched = np.zeros_like(mu)
    used = np.zeros(len(eigvals_l), dtype=bool)
    for j, lam in enumerate(eigvals_r):
        best = -1
        best_d = math.inf
        for k in range(len(eigvals_l)):
            if used[k]:
                continue
            d = abs(eigvals_l[k] - lam)
            if d < best_d:
                best_d = d
                best = k
        used[best] = True
        eigvals_l_sorted[j] = eigvals_l[best]
        mu_matched[:, j] = mu[:, best]

    # Compute c_j = <nu_0, psi_j> <mu_j, Phi> / <mu_j, psi_j>
    n_modes = len(eigvals_r)
    coeffs = np.zeros(n_modes, dtype=complex)
    initial_c = initial.astype(complex)
    probe_c = probe.astype(complex)
    for j in range(n_modes):
        psi_j = psi[:, j]
        mu_j = mu_matched[:, j]
        nu_psi = complex(np.dot(initial_c, psi_j))
        mu_phi = complex(np.dot(mu_j, probe_c))
        denom = complex(np.dot(mu_j, psi_j))
        if abs(denom) < 1e-300:
            coeffs[j] = 0.0
            continue
        coeffs[j] = nu_psi * mu_phi / denom

    return eigvals_r, coeffs


# ---------------------------------------------------------------------------
# The Gamma_A feasibility gate --- tail-onset t* and the window condition
# (live-route Thm 2a / Cor. two-probe). Same form as N*delta>1 elsewhere:
# no observation window past the contamination time, no tau.
# ---------------------------------------------------------------------------

def tail_onset_t_star(gamma_A: float, delta: float) -> float:
    """Guaranteed tail-onset time t* = max(0, (2/Delta) * log Gamma_A)
    (autocorrelation form, live-route Thm 2a).

    Gamma_A = ||QA|| / |c_2^A|  = contamination-to-signal of the probe (0 = perfectly
    faithful to the slow mode). Delta = mu_2 - mu_3 > 0 the spectral gap. A well-loaded
    slow mode (Gamma_A <= 1) gives t* = 0 (the tail dominates from the start); a
    contaminated probe (Gamma_A > 1) must wait 2/Delta * log Gamma_A before the mu_2
    term dominates. Delta = +inf or <= 0 (no resolved faster mode) => no contamination
    onset, t* = 0.
    """
    if delta <= 0 or math.isinf(delta):
        return 0.0
    if gamma_A <= 0:
        return 0.0
    return max(0.0, (2.0 / delta) * math.log(gamma_A))


def tail_window_feasible(T: float, gamma_A: float, delta: float) -> bool:
    """The feasibility gate T > t*(Gamma_A, Delta) (live-route Cor. two-probe).

    True iff the observation window T reaches past the guaranteed tail onset, so the
    slow-mode rate mu_2 (hence tau = -1/mu_2) is resolvable. This is the third, and
    the only genuinely computed, condition of the tau OK verdict --- the other two
    (verified isolated gap; faithfulness) are checks that something was supplied.
    """
    return T > tail_onset_t_star(gamma_A, delta)

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Reversible / self-adjoint spectral engine for the certified tau-route.

WHY THIS EXISTS. The live-route theorems
(PACKET_live_route: thm:abscissa, thm:auto) live in the pi-weighted space
L^2(rho_*) with DETAILED BALANCE, not euclidean l^2. Their precondition is
SELF-ADJOINTNESS in <f,g>_pi = sum_i pi_i f_i g_i, i.e. reversibility
pi_i P_ij = pi_j P_ji -- NOT euclidean normality. Two consequences the old
euclidean kappa_2 check got wrong:

  * a reversible chain is generally NOT a symmetric matrix (it is self-adjoint
    only in the pi-weighted inner product) -- checking euclidean symmetry/normality
    would reject exactly the lattice models the programme runs on;
  * a euclidean-normal but non-reversible operator (e.g. a circulant drift walk)
    would pass a normality test yet is OUTSIDE the theorems -- it must be refused.

So this module verifies self-adjointness in the DECLARED inner product, and -- for
a self-adjoint operator -- returns real eigenvalues from `eigh` of the symmetrised
S = D^{1/2} P D^{-1/2} (orthonormal eigenbasis, error ~ eps*||S||; no eigenvector
condition number needed). Jordan blocks / near-defects CANNOT occur on this route
(self-adjoint => always semisimple), which is the whole simplification.

This module returns FACTS (a `ReversibleSpectrum`), never applicability verdicts;
the 4-code verdicts are assembled downstream in codes.py from these facts.

Discipline: pure functions, no analyze() coupling. Detailed balance's pi is derived
from spanning-tree ratios + a Kolmogorov cycle consistency check (NOT the Perron
vector, whose condition ~ 1/(1-lambda_2) is worst exactly in the metastable regime).
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import List, Optional, Sequence, Tuple

import math
import numpy as np

_EPS = np.finfo(float).eps


@dataclass(frozen=True)
class ReversibleSpectrum:
    """Structured facts for the certified tau-route. No verdicts here.

    `selfadjoint` in {"verified","declared","none"}:
      verified -- an operator was supplied and self-adjointness in the declared
                  inner product was CHECKED and holds;
      declared -- only a (real) spectrum was supplied; self-adjointness is asserted
                  by the caller, not checkable here;
      none     -- an operator was supplied and is NOT self-adjoint in the declared
                  inner product (-> NOT_APPLICABLE downstream).
    `embeddable` is False when a negative eigenvalue dominates lambda_2 in modulus
    (no reversible generator embedding; rem:generator). `irreducible` is False when
    the leading eigenvalue 1 is degenerate (pi not unique) or the graph is
    disconnected. `degenerate` is True when lambda_2 has multiplicity > 1 (gap to the
    next DISTINCT level). `gap_below_band` is True when 1-lambda_2 or lambda_2-lambda_next
    sits under the numerical resolution band (cannot certify a strict gap).
    """
    selfadjoint: str
    reason: str = ""
    eigenvalues: Optional[Tuple[float, ...]] = None   # real, descending, when known
    lambda1: Optional[float] = None
    lambda2: Optional[float] = None                   # largest eigenvalue < 1
    lambda_next: Optional[float] = None               # next DISTINCT level below lambda2
    tau_spectral: Optional[float] = None              # -T_star/log(lambda2)
    delta: Optional[float] = None                     # mu2 - mu_next  (> 0)
    mu2: Optional[float] = None
    mu_next: Optional[float] = None
    band_rel: Optional[float] = None                  # dtau/tau ~ eps/(1-lambda2)
    embeddable: bool = True
    irreducible: bool = True
    degenerate: bool = False
    lambda2_multiplicity: int = 1
    gap_below_band: bool = False
    symmetry_residual: Optional[float] = None
    pi: Optional[Tuple[float, ...]] = None
    notes: Tuple[str, ...] = field(default_factory=tuple)


def stationary_by_detailed_balance(
    P: np.ndarray, tol: float
) -> Tuple[Optional[np.ndarray], bool, bool, str]:
    """Derive pi from detailed-balance ratios along a spanning tree, and test
    Kolmogorov consistency on the remaining edges.

    Returns (pi, reversible, irreducible, reason). pi is None when structure is
    incompatible with reversibility. `irreducible` reports graph connectivity
    (the spanning-tree walk covers all states iff connected). This is NOT the
    Perron vector (whose condition ~ 1/(1-lambda_2) degrades in the metastable
    regime); it uses only the ratios pi_j/pi_i = P_ij/P_ji.
    """
    n = P.shape[0]
    # structural reversibility: P_ij > 0 <=> P_ji > 0
    for i in range(n):
        for j in range(i + 1, n):
            if (P[i, j] > 0) != (P[j, i] > 0):
                return None, False, True, (
                    f"structural asymmetry P[{i},{j}]>0 xor P[{j},{i}]>0: no "
                    f"reversible pi exists (an edge is one-directional)."
                )
    log_pi = np.full(n, np.nan)
    log_pi[0] = 0.0
    stack = [0]
    while stack:
        i = stack.pop()
        for j in range(n):
            if j == i or P[i, j] <= 0:
                continue
            # pi_j/pi_i = P_ij/P_ji  =>  log pi_j = log pi_i + log(P_ij/P_ji)
            if math.isnan(log_pi[j]):
                log_pi[j] = log_pi[i] + math.log(P[i, j] / P[j, i])
                stack.append(j)
    if np.any(np.isnan(log_pi)):
        return None, True, False, (
            "the reachability graph is disconnected (reducible): pi is not "
            "unique on a single communicating class."
        )
    log_pi -= log_pi.max()
    pi = np.exp(log_pi)
    pi /= pi.sum()
    # Kolmogorov consistency on ALL edges: pi_i P_ij == pi_j P_ji
    scale = float(np.max(pi[:, None] * P))
    max_resid = 0.0
    for i in range(n):
        for j in range(i + 1, n):
            if P[i, j] > 0:
                resid = abs(pi[i] * P[i, j] - pi[j] * P[j, i])
                max_resid = max(max_resid, resid)
    if scale > 0 and max_resid / scale > tol:
        return pi, False, True, (
            f"detailed balance fails (Kolmogorov cycle residual "
            f"{max_resid/scale:.2e} > tol {tol:.1e}): not reversible."
        )
    return pi, True, True, ""


def _distinct_levels(eigs_desc: np.ndarray, rel_tol: float) -> List[Tuple[float, int]]:
    """Collapse (descending, real) eigenvalues into distinct levels with
    multiplicities, merging values within rel_tol of the spectral scale."""
    scale = max(1.0, float(np.max(np.abs(eigs_desc)))) if eigs_desc.size else 1.0
    levels: List[List[float]] = []
    for v in eigs_desc:
        if levels and abs(v - levels[-1][0]) <= rel_tol * scale:
            levels[-1].append(v)
        else:
            levels.append([v])
    return [(float(np.mean(g)), len(g)) for g in levels]


def analyze_reversible(
    *,
    operator: Optional[Sequence[Sequence[float]]] = None,
    spectrum: Optional[Sequence[float]] = None,
    T_star: float,
    inner_product="euclidean",
    selfadjoint_declared: bool = False,
    band_tol: Optional[float] = None,
) -> ReversibleSpectrum:
    """Verify self-adjointness in the declared inner product and extract the
    certified-route spectral facts. Exactly one of `operator` / `spectrum`.

    inner_product: "euclidean" (S = A, check symmetric), a pi-weight array
    (verify stationarity + detailed balance, S = D^{1/2} A D^{-1/2}), or "auto"
    (derive pi by detailed balance, then symmetrise).
    """
    if (operator is None) == (spectrum is None):
        raise ValueError("supply exactly one of operator= or spectrum=.")
    if T_star <= 0:
        raise ValueError("T_star must be > 0.")

    # ---- spectrum-only tier: self-adjointness is DECLARED, not checkable -------
    if spectrum is not None:
        eigs = np.asarray([complex(x) for x in spectrum])
        if np.max(np.abs(eigs.imag)) > 1e-9:
            return ReversibleSpectrum(
                selfadjoint="none",
                reason=("a bare spectrum with non-real eigenvalues cannot come from "
                        "a self-adjoint generator; declare a real spectrum or supply "
                        "the operator."),
            )
        return _from_real_eigs(
            np.sort(eigs.real)[::-1], T_star,
            selfadjoint=("declared" if selfadjoint_declared else "declared"),
            band_tol=band_tol, pi=None, sym_resid=None,
            note=("self-adjointness is DECLARED (spectrum-only): not verified; "
                  "faithfulness (c_2 != 0) also cannot be checked without a probe."),
        )

    # ---- operator tier: verify self-adjointness in the declared inner product --
    A = np.asarray(operator, dtype=float)
    if A.ndim != 2 or A.shape[0] != A.shape[1]:
        raise ValueError("operator must be a square matrix.")
    n = A.shape[0]
    if band_tol is None:
        band_tol = 100.0 * _EPS * n

    pi = None
    irreducible = True
    if isinstance(inner_product, str) and inner_product == "euclidean":
        S = A
    elif isinstance(inner_product, str) and inner_product == "auto":
        pi, reversible, irreducible, reason = stationary_by_detailed_balance(A, band_tol)
        if pi is None or not reversible:
            return ReversibleSpectrum(selfadjoint="none", reason=reason,
                                      irreducible=irreducible)
        d = np.sqrt(pi)
        S = (d[:, None] * A) / d[None, :]
    else:
        pi = np.asarray(inner_product, dtype=float)
        if pi.shape != (n,) or np.any(pi <= 0):
            raise ValueError("inner_product weight pi must be a positive length-n vector.")
        pi = pi / pi.sum()
        # verify stationarity pi P = pi (only meaningful for a stochastic operator)
        if abs(A.sum(axis=1).mean() - 1.0) < 1e-6:  # looks stochastic
            stat_resid = float(np.max(np.abs(pi @ A - pi)))
            if stat_resid > 1e-6:
                return ReversibleSpectrum(
                    selfadjoint="none",
                    reason=f"declared pi is not stationary (||piP-pi||={stat_resid:.2e}).",
                )
        d = np.sqrt(pi)
        S = (d[:, None] * A) / d[None, :]

    sym_resid = float(np.linalg.norm(S - S.T, "fro"))
    s_norm = float(np.linalg.norm(S, "fro")) or 1.0
    rel_resid = sym_resid / s_norm
    if rel_resid > band_tol:
        return ReversibleSpectrum(
            selfadjoint="none",
            reason=(f"not self-adjoint in the declared inner product "
                    f"(||S-S^T||/||S|| = {rel_resid:.2e} > tol {band_tol:.1e}); "
                    f"the live-route theorems require self-adjointness (reversibility), "
                    f"not mere normality."),
            symmetry_residual=rel_resid,
            pi=tuple(pi.tolist()) if pi is not None else None,
            irreducible=irreducible,
        )
    # self-adjoint: real spectrum via eigh (orthonormal eigenbasis, error ~ eps*||S||)
    Ssym = 0.5 * (S + S.T)
    eigs = np.linalg.eigvalsh(Ssym)[::-1]  # descending
    return _from_real_eigs(
        eigs, T_star, selfadjoint="verified", band_tol=band_tol,
        pi=pi, sym_resid=rel_resid, irreducible=irreducible,
    )


def _from_real_eigs(
    eigs_desc: np.ndarray, T_star: float, *, selfadjoint: str,
    band_tol: Optional[float], pi, sym_resid, irreducible: bool = True,
    note: Optional[str] = None,
) -> ReversibleSpectrum:
    notes: List[str] = [note] if note else []
    if band_tol is None:
        band_tol = 100.0 * _EPS * max(1, eigs_desc.size)
    lam1 = float(eigs_desc[0])
    # leading eigenvalue 1 degenerate -> pi not unique -> reducible
    reducible_lead = eigs_desc.size >= 2 and abs(eigs_desc[1] - lam1) <= band_tol * max(1.0, abs(lam1))
    if reducible_lead:
        return ReversibleSpectrum(
            selfadjoint=selfadjoint, reason="leading eigenvalue is degenerate "
            "(multiplicity > 1): the stationary state is not unique (reducible).",
            eigenvalues=tuple(eigs_desc.tolist()), lambda1=lam1,
            irreducible=False, pi=tuple(pi.tolist()) if pi is not None else None,
            symmetry_residual=sym_resid, notes=tuple(notes),
        )
    # lambda_2 = largest eigenvalue STRICTLY BELOW 1 (== leading), not largest modulus
    below = eigs_desc[eigs_desc < lam1 - band_tol * max(1.0, abs(lam1))]
    if below.size == 0:
        return ReversibleSpectrum(
            selfadjoint=selfadjoint, reason="no subdominant eigenvalue below the "
            "leading one; no relaxation mode.", eigenvalues=tuple(eigs_desc.tolist()),
            lambda1=lam1, irreducible=irreducible,
            pi=tuple(pi.tolist()) if pi is not None else None,
            symmetry_residual=sym_resid, notes=tuple(notes),
        )
    lam2 = float(below[0])
    # negative eigenvalue dominating lambda_2 in modulus -> not generator-embeddable
    most_neg = float(eigs_desc.min())
    embeddable = not (most_neg < 0 and abs(most_neg) > lam2 + band_tol)
    # distinct-level structure around lambda_2 (multiplicity + next distinct level)
    levels = _distinct_levels(eigs_desc, band_tol)
    lam2_level_idx = min(range(len(levels)), key=lambda k: abs(levels[k][0] - lam2))
    lam2_mult = levels[lam2_level_idx][1]
    lam_next = levels[lam2_level_idx + 1][0] if lam2_level_idx + 1 < len(levels) else None
    degenerate = lam2_mult > 1

    tau_spectral = None
    mu2 = mu_next = delta = band_rel = None
    gap_below_band = False
    if 0.0 < lam2 < 1.0:
        tau_spectral = -T_star / math.log(lam2)
        mu2 = math.log(lam2) / T_star
        # numerics band dtau/tau ~ eps/(1-lambda2) (eigh error, diverges lambda2->1)
        band_rel = _EPS / (1.0 - lam2) if lam2 < 1.0 else math.inf
        if (1.0 - lam2) <= band_tol:
            gap_below_band = True
        if lam_next is not None and lam_next > 0:
            mu_next = math.log(lam_next) / T_star
            delta = mu2 - mu_next
            if abs(lam2 - lam_next) <= band_tol * max(1.0, abs(lam2)):
                gap_below_band = True
    return ReversibleSpectrum(
        selfadjoint=selfadjoint,
        eigenvalues=tuple(eigs_desc.tolist()), lambda1=lam1, lambda2=lam2,
        lambda_next=lam_next, tau_spectral=tau_spectral, delta=delta,
        mu2=mu2, mu_next=mu_next, band_rel=band_rel, embeddable=embeddable,
        irreducible=irreducible, degenerate=degenerate,
        lambda2_multiplicity=lam2_mult, gap_below_band=gap_below_band,
        symmetry_residual=sym_resid,
        pi=tuple(pi.tolist()) if pi is not None else None, notes=tuple(notes),
    )

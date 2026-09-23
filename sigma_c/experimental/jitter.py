# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
The jitter channel — EXPERIMENTAL, dial-free.  NOT part of the certified route.

This is the jitter theorem of *The Parrot's Theorems*, Section 5 (`sec:jitter`),
with proofs in Appendices B.3 (`app:B3`) and B.4 (`app:B4`) and the pre-registered
two-channel run in Appendix C (`app:C`).  Zenodo doi:10.5281/zenodo.22066713.

It operates on **clamped-dial repeated readings of several instruments** — a block
of readings taken with every dial held fixed, repeated into a (blocks x instruments)
matrix — NOT on a single time series.  With the dials clamped, a slowly wandering
latent parameter shifts each block mean by a_i*delta, so the block-to-block
covariance is C = Var(delta) * a a^T (App B.3): rank one, leading eigenvector along
the dial susceptibilities a_i = d<O_i>/dgamma, leading eigenvalue Var(delta)|a|^2.
The four parts of the theorem:

  (i)   Detection   — off-diagonal covariance beyond sampling error rules out
                      "deterministic dial response + independent instrument noise".
  (ii)  Counting    — # shared hidden variables = # eigenvalues above the noise
                      floor (rank), read against a shuffle-null.
  (iii) Identification — the leading eigenvector gives the couplings a up to one
                      global scale; compare to the dial susceptibilities via cos(v1,a).
  (iv)  Cross-check — the jitter slope Cov(M1,M2)/Var(M2) equals the dial slope
                      dM1/dM2 from the scan, UNDER detailed balance (Onsager
                      regression); outside it, disagreement is diagnosis, not failure.

The App C two-channel run uses two covariance channels at once: **Channel 1**
(`pair_parity_rank`) = the per-shot covariance of co-measured pair parities (fast,
shot-level), and **Channel 2** (`jitter_rank`) = the drift covariance of block means
(slow); both feed the same detection/counting machinery.

This channel confirms **rank** and **peak location** independently of the dial and
checks detailed balance via the slope identity.  It delivers **no certified tau**:
the route from the shared peak to tau is the fallen sigma_c = rho_star*tau bridge.

Everything here is `validated=False` and emits NO certified `OK` verdict.  The
pre-registered validation (App C: peak agreement within one grid spacing,
cos(v1,a) >= 0.9, shuffle-null n=200 at the 95th percentile, drift-off -> floor) is
in PREREG_jitter.md / VALIDATION_RESULTS.md.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from typing import Optional, Tuple

import numpy as np


# ---------------------------------------------------------------------------
# (i) + (ii) Detection and counting — drift covariance + shuffle-null rank
# ---------------------------------------------------------------------------
def drift_covariance(block_means: np.ndarray, *, n_shots: Optional[int] = None) -> np.ndarray:
    """Block-to-block covariance of clamped-dial block means (App B.3). If `n_shots`
    is given, subtract the per-instrument shot-noise diagonal (1 - m^2)/n_shots, so
    what remains is the drift covariance Var(delta) a a^T."""
    a = np.asarray(block_means, dtype=float)
    C = np.cov(a, rowvar=False)
    if n_shots is not None:
        m = a.mean(axis=0)
        C = C - np.diag((1.0 - m ** 2) / n_shots)
    return C


@dataclass(frozen=True)
class JitterRank:
    """Detection (i) + counting (ii), read against a shuffle-null."""
    eigenvalues: tuple            # descending
    rank: Optional[int]           # # eigenvalues above the family-wise-corrected null
    null_threshold: tuple         # per-order family-wise (Bonferroni) null threshold
    detected: Optional[bool]      # (i): any shared cause above the null
    status: str = "reported"
    note: str = ""
    validated: bool = False

    def summary(self) -> str:
        head = "[experimental] jitter rank — detection (i) + counting (ii), NOT certified"
        if self.rank is not None:
            body = f"  rank above shuffle-null = {self.rank}   (detected shared cause: {self.detected})"
        else:
            body = f"  --  [{self.status}]"
        return f"{head}\n{body}\n  {self.note}"


def jitter_rank(
    block_means: np.ndarray,
    *,
    n_shots: Optional[int] = None,
    n_shuffle: int = 200,
    alpha: float = 0.95,
    seed: int = 0,
) -> JitterRank:
    """Count the shared hidden variables (App C: shuffle-null n=200 at the 95th
    percentile). Each column is permuted independently to destroy cross-instrument
    structure. Detection is a single test on the LEADING eigenvalue at the declared
    level (calibrated ~alpha at every instrument count); the rank is a **step-down**
    count (compare each order from the top to its null, stop at the first that fails),
    so the pure-noise rank cannot inflate with instrument count. (An earlier naive
    per-order sum inflated badly: 9.5% at d=3, 80% at d=30 — replaced.) Reports only
    relative to the null (a finite-N sample covariance is never flat).

    Counting (ii) is CONSERVATIVE by design: the shuffle-null preserves per-column
    variances, so a dominant first factor inflates the order-2 null and absorbs a
    weaker second factor (reported as rank 1). Empirically pinned (App C parameters):
    a DISJOINT second factor below ~30% of the first's power is missed, and a second
    factor under a COMMON dominant factor is missed at ANY strength. For a weak second
    factor use the `tetrad` instead — it stays sensitive well past where this step-down
    fires. This is the price of the false-positive control (it under-counts, not over-
    counts — the safe direction). EXPERIMENTAL.
    """
    a = np.asarray(block_means, dtype=float)
    if a.ndim != 2 or a.shape[1] < 2:
        return JitterRank((), None, (), None, status="not_applicable",
            note="need a (blocks, instruments>=2) matrix of clamped-dial block means.")
    nb, d = a.shape
    if nb <= d or not np.isfinite(a).all():
        return JitterRank((), None, (), None, status="insufficient_data",
            note=f"need finite blocks with N_blocks > instruments ({nb} <= {d}).")

    def _eigs(m):
        C = drift_covariance(m, n_shots=n_shots)
        return np.sort(np.clip(np.linalg.eigvalsh(C), 0.0, None))[::-1]

    obs = _eigs(a)
    rng = np.random.default_rng(seed)
    null = np.empty((n_shuffle, d))
    for i in range(n_shuffle):
        sh = np.empty_like(a)
        for j in range(d):
            sh[:, j] = a[rng.permutation(nb), j]
        null[i] = _eigs(sh)
    # Step-down rank test (parallel-analysis style): compare each eigenvalue order,
    # from the top, to its shuffle-null at the declared level; STOP at the first order
    # that does not exceed its null. Detection is then the single leading test (an
    # independent run confirmed the leading order is calibrated, ~alpha FP at every d),
    # and the step-down stops the pure-noise rank from inflating with instrument count
    # (rank>=2 needs two consecutive exceedances -- rare under the null). This replaces
    # a naive per-order sum, whose family-wise false-positive rate grew with d.
    per_order = np.percentile(null, alpha * 100, axis=0)
    rank = 0
    for _k in range(d):
        if obs[_k] > per_order[_k]:
            rank += 1
        else:
            break
    detected = bool(rank >= 1)
    note = (f"REPORTED (experimental, unvalidated): step-down rank {rank} over "
            f"{d} instruments (leading order at the declared level) -> "
            + ("a shared fluctuating cause is present (i); its count is the rank (ii)."
               if detected else
               "no shared cause above the null (consistent with independent instrument noise).")
            + " App C null; a necessary check, not a validated discriminator (PREREG_jitter.md).")
    return JitterRank(tuple(float(x) for x in obs), rank, tuple(float(x) for x in per_order),
                      detected, status="reported", validated=False, note=note)


# ---------------------------------------------------------------------------
# App C Channel 1 — the per-shot covariance of co-measured pair parities
# ---------------------------------------------------------------------------
def pair_parity_covariance(shot_parities: np.ndarray) -> np.ndarray:
    """App C Channel 1: the per-shot covariance of co-measured pair parities.

    `shot_parities` is (n_shots, n_pairs), each entry a +/-1 parity of a qubit pair
    read in the SAME shot. Unlike Channel 2 (the slow drift of block means), this is
    the FAST, shot-level covariance: a shared per-shot cause (e.g. correlated readout
    error) shows as off-diagonal covariance; independent pairs give none."""
    return drift_covariance(np.asarray(shot_parities, dtype=float))


def pair_parity_rank(
    shot_parities: np.ndarray,
    *,
    n_shuffle: int = 200,
    alpha: float = 0.95,
    seed: int = 0,
) -> JitterRank:
    """App C Channel 1 detection (i) + counting (ii): the rank of the per-shot
    pair-parity covariance against the shuffle-null — the shot-level twin of
    Channel 2's `jitter_rank` (which runs on block means). Same detection/counting
    machinery (step-down leading test + shuffle-null), run on the per-shot matrix with
    no block-mean shot-noise term. Run this ALONGSIDE `jitter_rank` for the full App C
    two-channel acquisition. EXPERIMENTAL — no verdict, validated=False."""
    return jitter_rank(shot_parities, n_shots=None, n_shuffle=n_shuffle,
                       alpha=alpha, seed=seed)


# ---------------------------------------------------------------------------
# (B.4) Tetrad — the k>=3 one-factor blind test
# ---------------------------------------------------------------------------
def tetrad(C: np.ndarray, *, i: int = 0, j: int = 1, k: int = 2,
           instrument_floor: float = 0.0) -> dict:
    """Spearman's tetrad (App B.4): under ONE shared factor, C_ij C_ik / C_jk equals
    a_i^2 sigma_y^2, checkable against C_ii minus the instrument floor. A second
    factor breaks it. Needs >= 3 instruments. Returns the two sides and their ratio
    (~1 => one-factor consistent). EXPERIMENTAL, reported."""
    C = np.asarray(C, dtype=float)
    if C.shape[0] < 3:
        return {"status": "not_applicable", "note": "tetrad needs >= 3 instruments."}
    Cjk = C[j, k]
    if abs(Cjk) < 1e-15:
        return {"status": "insufficient_data", "note": "C_jk ~ 0; tetrad undefined."}
    lhs = C[i, j] * C[i, k] / Cjk          # = a_i^2 sigma_y^2 under one factor
    rhs = C[i, i] - instrument_floor
    ratio = lhs / rhs if abs(rhs) > 1e-15 else float("inf")
    return {"status": "reported", "tetrad_lhs": float(lhs), "diag_minus_floor": float(rhs),
            "ratio": float(ratio), "validated": False,
            "note": ("one-factor consistent (ratio ~ 1)" if 0.7 <= ratio <= 1.4
                     else "one-factor identity broken (ratio far from 1) -> a second factor")}


# ---------------------------------------------------------------------------
# (iii) Identification — leading coupling vs the dial susceptibilities
# ---------------------------------------------------------------------------
def leading_coupling(C: np.ndarray) -> np.ndarray:
    """The leading eigenvector of the drift covariance = the couplings a up to one
    global scale (App B.3). Sign-fixed so its largest component is positive."""
    C = np.asarray(C, dtype=float)
    w, V = np.linalg.eigh(C)
    v1 = V[:, int(np.argmax(w))]
    return v1 * np.sign(v1[int(np.argmax(np.abs(v1)))])


def cos_coupling(v1: np.ndarray, a: np.ndarray) -> float:
    """cos(v1, a): alignment of the jitter's leading eigenvector with the dial
    susceptibilities a_i = d<O_i>/dgamma (App C criterion: >= 0.9). Identification (iii)."""
    v1 = np.asarray(v1, dtype=float); a = np.asarray(a, dtype=float)
    na, nv = np.linalg.norm(a), np.linalg.norm(v1)
    if na < 1e-15 or nv < 1e-15:
        return 0.0
    return float(abs(v1 @ a) / (na * nv))


# ---------------------------------------------------------------------------
# (iv) Cross-check — the jitter slope vs the dial slope (the FDT test)
# ---------------------------------------------------------------------------
def jitter_slope(M1: np.ndarray, M2: np.ndarray) -> float:
    """The jitter slope Cov(M1, M2)/Var(M2) from clamped-dial repeated readings (iv)."""
    M1 = np.asarray(M1, dtype=float); M2 = np.asarray(M2, dtype=float)
    C = np.cov(M1, M2)                       # ddof=1 on both terms, as in the paper
    if C[1, 1] <= 0.0:
        return float("nan")
    return float(C[0, 1] / C[1, 1])


@dataclass(frozen=True)
class TwoChannelReading:
    """(iv) cross-check: the jitter slope vs the dial slope. REPORTED, not pass/fail."""
    jitter_slope: Optional[float]
    dial_slope: Optional[float]
    rel_deviation: Optional[float]
    combined_band: Optional[float]
    consistent_with_db: Optional[bool]
    status: str = "reported"
    note: str = ""
    validated: bool = False


def two_channel_fdt(
    slope_jitter: Optional[float],
    slope_dial: Optional[float],
    *,
    band_jitter: Optional[float] = None,
    band_dial: Optional[float] = None,
) -> TwoChannelReading:
    """Onsager's regression cross-check (iv): the jitter slope must equal the dial
    slope UNDER detailed balance; a deviation beyond the combined band is diagnosis
    (how far from equilibrium), not failure. REPORTED with bands, never a certified
    verdict. EXPERIMENTAL."""
    if slope_jitter is None or slope_dial is None:
        return TwoChannelReading(slope_jitter, slope_dial, None, None, None,
            status="needs_both_channels",
            note="needs BOTH slopes (jitter_slope from clamped readings; dial slope from the scan).")
    if not (np.isfinite(slope_jitter) and np.isfinite(slope_dial)):
        return TwoChannelReading(slope_jitter, slope_dial, None, None, None,
            status="insufficient_data", note="a slope is not finite.")
    denom = max(abs(slope_jitter), abs(slope_dial))
    if denom == 0.0:
        return TwoChannelReading(slope_jitter, slope_dial, None, None, None,
            status="insufficient_data", note="both slopes are zero.")
    rel = abs(slope_jitter - slope_dial) / denom
    if band_jitter is None or band_dial is None:
        return TwoChannelReading(slope_jitter, slope_dial, rel, None, None,
            status="reported",
            note=(f"jitter slope vs dial slope differ by {rel:.2%} REPORTED; consistency "
                  "cannot be judged without both bands (a deviation below the combined "
                  "band is estimation error). Supply band_jitter and band_dial."))
    combined = float(band_jitter + band_dial)
    consistent = rel <= combined
    note = (f"jitter slope = dial slope to {rel:.2%} <= combined band {combined:.2%}: "
            "CONSISTENT with detailed balance (Onsager regression, one curve two ways)."
            if consistent else
            f"jitter slope and dial slope differ by {rel:.2%} > combined band {combined:.2%}: "
            "a measured departure from detailed balance (diagnosis, not failure).")
    return TwoChannelReading(slope_jitter, slope_dial, rel, combined, bool(consistent),
        status="reported", validated=False, note=note)


# ---------------------------------------------------------------------------
# Two-channel peak comparison (App C) — dial chi(gamma) vs jitter lambda_1(gamma)
# ---------------------------------------------------------------------------
@dataclass(frozen=True)
class PeakAgreement:
    dial_peak: Optional[float]
    jitter_peak: Optional[float]
    within_one_cell: Optional[bool]
    cos_v1_a: Optional[float] = None
    status: str = "reported"
    note: str = ""
    validated: bool = False


def peak_agreement(
    grid: np.ndarray,
    dial_chi: np.ndarray,
    jitter_lambda1: np.ndarray,
    *,
    cos_v1_a: Optional[float] = None,
) -> PeakAgreement:
    """App C two-channel signature: the dial channel chi_W(gamma) and the jitter
    channel lambda_1(gamma) must peak at the same grid point (within one grid
    spacing), and cos(v1,a) >= 0.9. REPORTED, not a certified verdict."""
    g = np.asarray(grid, dtype=float)
    cd = np.asarray(dial_chi, dtype=float)
    lj = np.asarray(jitter_lambda1, dtype=float)
    if g.ndim != 1 or g.size < 3 or cd.shape != g.shape or lj.shape != g.shape:
        return PeakAgreement(None, None, None, cos_v1_a, status="insufficient_data",
            note="need matching 1-D arrays grid, dial_chi, jitter_lambda1 of length >= 3.")
    kD, kJ = int(np.argmax(cd)), int(np.argmax(lj))
    within = abs(kD - kJ) <= 1
    cos_ok = "" if cos_v1_a is None else \
        f"; cos(v1,a) = {cos_v1_a:.3f} ({'>=' if cos_v1_a >= 0.9 else '<'} 0.9)"
    note = (f"dial peak gamma={g[kD]:.3g}, jitter peak gamma={g[kJ]:.3g} "
            + ("(within one grid cell)" if within else "(NOT within one grid cell)")
            + cos_ok + ". App C signature; experimental, unvalidated.")
    return PeakAgreement(float(g[kD]), float(g[kJ]), bool(within), cos_v1_a,
                         status="reported", validated=False, note=note)

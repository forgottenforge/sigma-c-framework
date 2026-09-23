# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
The 4-code applicability contract, per output quantity.

Every numeric output carries its OWN verdict --- OK / NOT_APPLICABLE /
NOT_RESOLVABLE / NOT_IDENTIFIED --- with a human reason and, for NOT_IDENTIFIED,
the *remediation*: what the caller could supply to change it. That is the one
code that says "not 'I cannot', but 'from here it is your choice, and here is
what is missing'"; the other three are refusals the caller cannot lift with
these data.

These codes are a FORMALISATION of the verdict analyze() already produces, not a
new threshold. `sigma_c_verdict()` reads classify()'s existing output
(regime.geometric, the A1 declaration, the grid-boundary / operational-floor
flags) and relabels it. On the golden anchors it must reproduce today's verdict
exactly (test_codes_invariance.py); where it would differ, that is a finding
about classify(), recorded in THEOREM_MAP BEFORE any code change --- not a silent
"improvement".

tau's verdicts are the two-field route below (tau_abscissa_verdict /
window_readability_verdict): the certified spectral abscissa, and its
finite-window readability with the T > t*(Gamma_A, Delta) tail-window condition
and the one-sided two-probe note.
"""
from __future__ import annotations

from dataclasses import dataclass
from enum import Enum
from typing import TYPE_CHECKING

if TYPE_CHECKING:  # avoid a circular import at runtime
    from sigma_c.result import Result, TwoProbeResult


class Code(str, Enum):
    """The four applicability codes. Two kinds of 'no':

    - OK: the value is trustworthy under the method's preconditions.
    - NOT_APPLICABLE / NOT_RESOLVABLE: refusals the caller cannot lift with these
      data (precondition of the method violated / below resolution).
    - NOT_IDENTIFIED: a *supplyable* precondition is missing; `remediation` names
      what to bring (a two-probe comparison, a verified spectrum with a gap).
    """

    OK = "OK"
    NOT_APPLICABLE = "NOT_APPLICABLE"
    NOT_RESOLVABLE = "NOT_RESOLVABLE"
    NOT_IDENTIFIED = "NOT_IDENTIFIED"


@dataclass(frozen=True)
class Verdict:
    code: Code
    reason: str
    remediation: str = ""  # only meaningful for NOT_IDENTIFIED

    def __str__(self) -> str:
        s = f"{self.code.value} — {self.reason}"
        if self.remediation:
            s += f"  [to identify, supply: {self.remediation}]"
        return s

    @property
    def ok(self) -> bool:
        return self.code is Code.OK

    def as_dict(self) -> dict:
        """Serialize the verdict. Used by Result.to_dict(), where the verdict is
        DERIVED at serialization time (not stored) so it cannot go stale."""
        return {
            "code": self.code.value,
            "reason": self.reason,
            "remediation": self.remediation,
        }


def sigma_c_verdict(result: "Result") -> Verdict:
    """Per-output code for `result.sigma_c`. A relabel of the existing verdict.

    Mapping (formalisation of what analyze() already does, plus the 6.0
    input-degeneracy / convention-stability guards recorded in THEOREM_MAP):

      - preprocessing NOT scale-equivariant (A1 violated at the observable layer)
        -> NOT_APPLICABLE  (matches today's "exploratory / .falsifiable == False")
      - observable constant to numerical noise (regime.degenerate_O) ->
        NOT_IDENTIFIED  (6.0 guard: sibling of the NaN guard)
      - resolved peak count over the declared ceiling (peak_count_over_ceiling)
        -> NOT_RESOLVABLE  (6.0 guard: many 'modes' off one dial is noise)
      - geometric regime III (chi_O monotone, no interior peak; sigma_c is the
        bottom value / None) -> NOT_RESOLVABLE  (matches today's sigma_c = None)
      - peak set not stable across the declared prominence window r in [r0/2, 2*r0]
        (convention_window_stable() is False) -> NOT_RESOLVABLE  (6.0 guard: an OK
        that holds only under a narrow convention is not a measurement)
      - otherwise regime I or II with an interior peak -> OK. Grid-boundary and
        operational-floor flags are carried as caveats in the reason string.

    sigma_c returns NOT_IDENTIFIED only for the degenerate-observable guard; the
    tau gate is the other NOT_IDENTIFIED path.
    """
    if result.preprocessing_scale_equivariant is False:
        return Verdict(
            Code.NOT_APPLICABLE,
            "A1 scale-equivariance violated at the observable layer: the "
            "preprocessing carries an absolute scale, so the analytic rho_star "
            "anchor does not apply (paper A1); the reading is exploratory, not a "
            "Def 6.7 falsification.",
        )
    # Input-degeneracy guards (siblings of the NaN guard: a silent input
    # degeneracy must not be read as a measurement). Both are recorded in
    # THEOREM_MAP as a deliberate verdict change from the pre-6.0 behaviour.
    if result.regime.degenerate_O:
        return Verdict(
            Code.NOT_IDENTIFIED,
            "the observable has no variation above the numerical-noise floor "
            "(constant to rounding precision): chi_O = |dO/dlog sigma| is 0 up to "
            "float error, so any 'peak' is noise in the log-derivative -- there is "
            "no scale to resolve.",
            "an observable that actually varies with the dial (range(O) above the "
            "numerical-noise floor eps*max|O|)",
        )
    if result.regime.peak_count_over_ceiling:
        return Verdict(
            Code.NOT_RESOLVABLE,
            f"resolved peak count {result.regime.resolved_peak_count} exceeds the "
            f"declared ceiling max_resolved_peaks={result.regime.max_resolved_peaks}: "
            "this many 'modes' read off one dial is noise, not structure (the vector "
            "sigma_c is still reported). If these are noise on a broader feature, "
            "smoothing is YOUR declared preprocessing decision (one of the five "
            "decision points) -- smooth the observable yourself before analysis and "
            "record that you did; the kernel will not smooth silently. Alternatively "
            "declare a coarser resolution (raise min_prominence_ratio / "
            "max_resolved_peaks), which is recorded in the verdict.",
        )
    geom = result.regime.geometric
    if geom == "III_geom":
        return Verdict(
            Code.NOT_RESOLVABLE,
            "no interior susceptibility peak (regime III): chi_O is monotone, so "
            "sigma_c is the bottom value, not a resolved interior scale "
            "(paper Thm 9.1).",
        )
    # Load-bearing convention stability (the last known path to a false OK): a
    # peak SET (scalar OR vector) whose count is not stable across the declared
    # window r in [r0/2, 2*r0] exists only under a narrow convention setting -- a
    # Layer-3 reading that runs away under the convention screw is not a resolved
    # scale, which is itself the finding. Refuse it. White noise (whether it lands
    # in regime I or a low-count regime II for a given r) fails this at every r.
    if result.convention_window_stable() is False:
        r0 = result.regime.min_prominence_ratio
        return Verdict(
            Code.NOT_RESOLVABLE,
            f"the peak set exists only under a narrow prominence convention: the "
            f"resolved peak count is NOT stable across r in [{r0/2:.3g}, "
            f"{min(1.0, r0*2):.3g}] (an octave of min_prominence_ratio={r0:.3g} each "
            f"way). A Layer-3 reading that runs away under the convention screw is "
            f"not a resolved scale. Smoothing is YOUR declared preprocessing decision "
            f"(recorded, not imposed by the kernel): smooth the observable yourself, "
            f"or resolve the convention dependence, before reading sigma_c.",
        )
    caveats = []
    if result.sigma_c_at_grid_boundary:
        caveats.append(
            "peak at the scan-grid boundary --- the true peak may lie outside the "
            "scan range; widen the grid before relying on the value"
        )
    if result.regime.operational_floor_triggered:
        caveats.append(
            "candidate peak below the operational signal/noise floor (paper Def 8.8)"
        )
    reason = f"regime {geom}: interior susceptibility peak resolved"
    if caveats:
        reason += " (caveats: " + "; ".join(caveats) + ")"
    return Verdict(Code.OK, reason)


# Remediation string for tau's NOT_IDENTIFIED: what the caller can SUPPLY to
# lift it (the one code that is a choice-point, not a dead end).
_TAU_REMEDIATION = (
    "a two-probe agreement -- run two_probe_test(r1, r2) on two faithful probes of "
    "the SAME system (the concrete next step when you have only (x, y) and no "
    "spectrum) -- OR a verified spectrum with an isolated real gap (reversible / "
    "self-adjoint, lambda_2 isolated, |lambda_2|>|lambda_3|), passed as "
    "analyze(..., spectrum=..., T_star=..., sigma_axis=...); in either case the "
    "tail-window condition T > t*(Gamma_A, Delta) must hold"
)

# Remediation when the tau-route ran (a spectrum/operator was supplied) but the
# author-pinned sigma_axis declaration is missing: without it there is NO tau
# (the sigma-axis unit and the gate window are undeclared). A major-release break
# of the old OK-on-an-undeclared-assumption behaviour.
_SIGMA_AXIS_REMEDIATION = (
    "declare sigma_axis in {evolution_time, window_time, other} (+ T_obs unless "
    "evolution_time): evolution_time = sigma is time in T_star units (gate uses "
    "max(sigma-grid)); window_time = sigma is time, gate uses T_obs (recording "
    "length); other = sigma is not time, gate uses T_obs, no bridge-vs-spectrum "
    "comparison"
)


def tau_two_probe_verdict(tp: "TwoProbeResult") -> Verdict:
    """tau verdict from a two-probe test. A PASS never yields a bare OK: the test
    is ONE-SIDED (deviation => trouble; passing =/=> faithful, by the common-mode
    failure of Cor. twoprobe). So a pass is OK *with the one-sided caveat carried
    in the reason*, and only for two analytic probes; a fit-based pass cannot
    certify. A deviation is NOT_IDENTIFIED (trouble, tau not pinned); a
    precondition failure of the test itself is NOT_RESOLVABLE / NOT_APPLICABLE.
    """
    if tp.passed:
        if not tp.both_analytic:
            return Verdict(
                Code.NOT_IDENTIFIED,
                "two-probe agreement, but at least one probe is fit-based "
                "(not analytic rho_star): the test is exploratory and cannot "
                "certify tau.",
                "both probes with analytic rho_star (a registry window)",
            )
        return Verdict(
            Code.OK,
            "two analytic faithful probes agree on the sigma_c/rho_star ratio "
            "within the declared threshold. ONE-SIDED: passing does NOT certify "
            "faithfulness --- a common-mode failure (both probes landing on the "
            "same faster mode) also passes (Cor. twoprobe); and T>t*(Gamma_A,Delta) "
            "is not machine-checked. Deviation would have been the firm signal.",
        )
    # failed: the cause decides which 'no' it is
    if tp.cause == "regime_ii":
        return Verdict(
            Code.NOT_RESOLVABLE,
            "at least one probe is multi-mode (regime II): no single tau to agree on.",
        )
    if tp.cause in ("probes_not_distinct", "indeterminate"):
        return Verdict(
            Code.NOT_APPLICABLE,
            f"two-probe precondition failed ({tp.cause}): the test cannot be run as "
            "a faithfulness check on these two results.",
        )
    # faithfulness_break / system_disagreement / any deviation
    return Verdict(
        Code.NOT_IDENTIFIED,
        f"two-probe DEVIATION (cause: {tp.cause}) beyond the threshold: the probes "
        "disagree on tau, so it is not identified (deviation => trouble, the firm "
        "side of the one-sided test).",
        _TAU_REMEDIATION,
    )


# ===========================================================================
# two-field tau verdicts (self-adjoint reversible route)
# ===========================================================================

def tau_abscissa_verdict(result: "Result") -> Verdict:
    """thm:abscissa: is tau_abscissa a certified property of the operator?

    OK needs self-adjointness in the DECLARED inner product (reversibility), a
    discrete isolated lambda_2, and faithfulness (assumed at the observable). NO
    strict gap required (thm:abscissa holds even for a degenerate mu_2). The fault
    line is self-adjointness -> non-self-adjoint / non-embeddable / reducible are
    NOT_APPLICABLE (the theorems do not cover them; a non-normal rate "may survive
    via other estimates" but that is outside what these theorems certify).
    """
    if result.preprocessing_scale_equivariant is False:
        return Verdict(
            Code.NOT_APPLICABLE,
            "A1 scale-equivariance violated at the observable layer; the analytic "
            "anchor does not apply (paper A1). Same precondition as sigma_c.",
        )
    rv = result._rev
    if rv is None:
        return Verdict(
            Code.NOT_IDENTIFIED,
            "no spectrum or operator supplied: the emitted tau_bridge is the "
            "load-dominant read-out (sigma_c/rho_star), not the certified spectral "
            "abscissa; the intrinsic relaxation time is not identified from (x, y) alone.",
            "operator= (with inner_product=) or spectrum= (+ T_star, sigma_axis), OR "
            "run two_probe_test(r1, r2) on two faithful probes of the same system",
        )
    if rv.selfadjoint == "none":
        return Verdict(
            Code.NOT_APPLICABLE,
            "not self-adjoint in the declared inner product: " + rv.reason + " The "
            "live-route theorems (thm:abscissa/thm:auto) live in L^2(rho_*) and require "
            "self-adjointness (reversibility), not mere normality.",
        )
    if not rv.irreducible:
        return Verdict(
            Code.NOT_APPLICABLE,
            "reducible: " + rv.reason + " A unique stationary rho_* (hence a unique "
            "relaxation abscissa) needs an irreducible operator.",
        )
    if not rv.embeddable:
        return Verdict(
            Code.NOT_APPLICABLE,
            "a negative eigenvalue dominates lambda_2 in modulus: the modulus-dominant "
            "mode is sign-alternating and not embeddable in a reversible generator "
            "(rem:generator, all mu_k <= 0), so the spectral-abscissa route does not apply. "
            "(A lazy step (P+I)/2 would change the value -- declare it yourself, the kernel "
            "will not do it silently.)",
        )
    if rv.tau_spectral is None:
        return Verdict(
            Code.NOT_RESOLVABLE,
            "no decaying subdominant mode (no lambda_2 in (0,1)): there is no finite "
            "relaxation time to identify.",
        )
    if rv.selfadjoint == "declared":
        return Verdict(
            Code.OK,
            "certified spectral abscissa tau = -T*/log(lambda_2) (thm:abscissa) -- but on "
            "TWO DECLARED hypotheses (spectrum-only): self-adjointness is asserted, not "
            "verified, and faithfulness (c_2 != 0) cannot be checked without a probe. "
            "Supply the operator + a probe to verify both (see selfadjoint='declared').",
        )
    return Verdict(
        Code.OK,
        "certified spectral abscissa tau = -T*/log(lambda_2) (thm:abscissa): self-adjoint "
        "in the declared inner product (real spectrum, orthonormal eigenbasis), discrete "
        "isolated lambda_2. Faithfulness (c_2 != 0) is assumed at the observable.",
    )


def window_readability_verdict(result: "Result") -> Verdict:
    """thm:auto: is the RATE readable in a finite observation window T?

    Inherits tau_abscissa's refusals (no certified abscissa -> no window read-off).
    Then needs a STRICT gap and T > max(t*, 1/|mu_3|); the 1/|mu_3| floor is a genuine
    hypothesis of thm:auto (the corrected cor:twoprobe gate
    is T > max{t*_A, t*_B, 1/|mu_3|}). NB this gate makes the rate BOUND APPLICABLE, it
    is not itself an accuracy guarantee -- the actual error is result.rate_error_bound,
    which the caller compares to their own tolerance (dominance onset t*, applicability
    floor 1/|mu_3|, and accuracy onset are three distinct thresholds). Reports
    result.rate_error_bound. Spectrum-only (no probe) can never be OK -> NOT_IDENTIFIED.
    """
    aa = tau_abscissa_verdict(result)
    if aa.code is not Code.OK:
        return aa  # NOT_APPLICABLE / NOT_RESOLVABLE / NOT_IDENTIFIED propagate
    rv = result._rev
    if rv.selfadjoint == "declared":
        return Verdict(
            Code.NOT_IDENTIFIED,
            "self-adjointness and faithfulness are only DECLARED (spectrum-only): reading "
            "the RATE in a finite window (thm:auto) needs a VERIFIED operator AND a probe "
            "to compute Gamma_A (contamination-to-signal). It can never be OK from a bare "
            "spectrum.",
            "the operator (matrix) AND a probe vector, so self-adjointness is verified and "
            "Gamma_A can be computed",
        )
    if rv.gap_below_band:
        return Verdict(
            Code.NOT_RESOLVABLE,
            "strict gap not determinable: |lambda_2 - lambda_next| (or 1 - lambda_2) sits "
            "under the numerical resolution band, so an exact degeneracy cannot be "
            "distinguished from a tiny gap.",
        )
    if rv.delta is None or rv.mu_next is None:
        return Verdict(
            Code.NOT_RESOLVABLE,
            "no strict spectral gap below lambda_2 (no next distinct level to set Delta), "
            "so the finite-window rate bound has no gap to decay on.",
        )
    if result.gamma_A is None:
        return Verdict(
            Code.NOT_IDENTIFIED,
            "a strict gap is present but Gamma_A (the probe's contamination-to-signal) was "
            "not supplied, so the tail-window condition cannot be checked.",
            "gamma_A (contamination-to-signal ||QA||/|c_2| of the probe)",
        )
    T = result.tau_window_T
    if T is None:
        return Verdict(
            Code.NOT_IDENTIFIED,
            "a strict gap and gamma_A are present but the observation window T is "
            "undeclared.",
            "sigma_axis (+ T_obs unless evolution_time), which fixes the window T",
        )
    import math
    from sigma_c.core.faithfulness import tail_onset_t_star
    t_star = tail_onset_t_star(result.gamma_A, rv.delta)
    mu3_floor = 1.0 / abs(rv.mu_next) if rv.mu_next != 0 else math.inf
    t_gate = max(t_star, mu3_floor)
    degen_note = ""
    if rv.degenerate:
        degen_note = (
            f" NB lambda_2 is degenerate (multiplicity {rv.lambda2_multiplicity}); "
            "readability rests on rem:gap (the degenerate bound "
            "carries by re-indexing over DISTINCT levels via Pi_2, Delta = mu_2 - mu_next). "
            "For a CROSS-probe in this degenerate case, faithfulness must be <Pi_2 A, Pi_2 B> "
            "!= 0, not merely both non-zero -- the code's single-operator route is unaffected."
        )
    if T > t_gate:
        reb = result.rate_error_bound
        reb_str = f"{reb:.3g}" if reb is not None else "n/a"
        return Verdict(
            Code.OK,
            f"the rate BOUND is applicable in this window: strict gap Delta={rv.delta:.3g}, "
            f"T={T:.3g} > max(t*={t_star:.3g}, 1/|mu_3|={mu3_floor:.3g}) (thm:auto with the "
            f"1/|mu_3| floor; corrected cor:twoprobe gate). The "
            f"relative rate-error bound is {reb_str} -- a REPORTED number, not an accuracy "
            f"verdict (dominance != accuracy): compare it to your own tolerance.{degen_note}",
        )
    return Verdict(
        Code.NOT_IDENTIFIED,
        f"the observation window is too short to read the rate to the thm:auto bound: "
        f"T={T:.3g} <= max(t*={t_star:.3g}, 1/|mu_3|={mu3_floor:.3g}).",
        f"a longer observation window T > {t_gate:.3g}, or a more faithful probe (smaller "
        f"gamma_A)",
    )

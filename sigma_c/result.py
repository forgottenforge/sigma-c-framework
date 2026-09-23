# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Result object — every output carries provenance, regime layers, and
gamma_O stability indicator.

Enforces:
- (Enforcement 2) rho_star_source as a typed field.
- (Enforcement 4) ⊥ as None, never NaN, never default exception.
- (Enforcement 5) two-probe failure with cause branch — see TwoProbeResult.
- (Enforcement 6) detrended flag, off by default.
- (Enforcement 8) smoothing parameters logged.

Cite: paper Def 2.2 (def:sigmac) for the partial-functional convention.
"""
from __future__ import annotations
from dataclasses import dataclass, field
from typing import Any, List, Optional, Tuple, Union

import numpy as np

from sigma_c.framework import Framework


# ---------------------------------------------------------------------------
# Regime verdict — two named convention-operations (peak-count + isolated-gap
# boolean) + a measurement floor. (Was paper §8's "three-layer trichotomy";
# those §8 theorem labels are ARCHIVED (a private math archive, not shipped).)
# ---------------------------------------------------------------------------

GeometricRegime = str  # "I_geom" | "II_geom" | "III_geom"
SpectralRegime = str   # "gap" | "no_gap" | None


@dataclass(frozen=True)
class Trichotomy:
    """
    The regime verdict — two named convention-operations, NOT theorems.

    The paper's §8 "trichotomy" (Sätze 8.3/8.6/8.11/8.16) is a paper-only
    formalism archived privately, not shipped (THEOREM_MAP.md marks those §8
    labels ARCHIVED). The live verdict carries:

      - geometric: a PEAK COUNT under the DECLARED convention
        `min_prominence_ratio` (recorded below — it governs the I/II/III verdict
        and is NOT free);
      - spectral: an ISOLATED-LEADING-GAP boolean ("gap"/"no_gap"/None) on a
        supplied spectrum — spectrum-only, None when no spectrum was provided;
      - operational_floor_triggered: the measurement-side floor bool.
    """
    geometric: GeometricRegime
    """One of 'I_geom' (one peak), 'II_geom' (>=2 peaks), 'III_geom' (none).
    A PEAK COUNT under the declared `min_prominence_ratio` convention — a named
    convention-operation, NOT a theorem (that threshold decides whether a bump
    counts as a peak, hence the regime)."""

    spectral: Optional[SpectralRegime] = None
    """Isolated-leading-gap boolean on a supplied spectrum:
    'gap' iff |lambda_2| < |lambda_1| strictly, 'no_gap' otherwise,
    None when no spectrum was supplied. SPECTRUM-ONLY: it is NOT
    observable-faithful — faithfulness runs through Gamma_A (probe-relative),
    not this gap. A named convention-operation, NOT a theorem."""

    operational_floor_triggered: bool = False
    """True iff the candidate peak amplitude is below eta_O * ||O||.
    Cite: def:noise-floor-diagnostic. Measurement-side diagnostic."""

    degenerate_O: bool = False
    """True iff the observable is constant to a few ULP
    (range(O) <= DEGENERATE_ULP_FACTOR * eps * max|O|, or O == 0). A constant
    observable has chi_O = 0 up to rounding, so any 'peak' is float noise in the
    log-derivative -- NOT a resolved scale. This is the sibling of the NaN guard:
    a silent INPUT DEGENERACY. It makes sigma_c NOT_IDENTIFIED (remediable: supply
    an observable that actually varies with the dial), never a spurious OK. A
    near-constant O whose tiny variation is above this ULP gate but whose chi peaks
    sit under the per-peak noise floor is caught downstream as NOT_RESOLVABLE."""

    resolved_peak_count: int = 0
    """Number of interior peaks that cleared BOTH the relative prominence
    convention AND the absolute numerical-noise floor."""

    max_resolved_peaks: int = 5
    """Declared ceiling on the resolved peak count (a convention, NOT a theorem).
    More than this many distinct 'modes' read off a single dial is noise, not
    structure; above it sigma_c is NOT_RESOLVABLE (the peaks are still listed).
    Tunable via analyze(max_resolved_peaks=...)."""

    peak_count_over_ceiling: bool = False
    """True iff resolved_peak_count > max_resolved_peaks. Makes sigma_c
    NOT_RESOLVABLE: 'peak count beyond the resolution, smooth or declare a coarser
    resolution'. The vector sigma_c is still reported."""

    chi_abs_floor: float = 0.0
    """The per-peak absolute numerical-noise floor used at classify time
    (NOISE_FLOOR_FACTOR * eps * max|O| / min(dlog)). Carried so the derived
    convention-stability recount uses the SAME peak-finding as the verdict."""

    eta_O: float = 0.0       # 0 = pure-mathematics limit
    min_prominence_ratio: float = 0.10
    """Declared prominence convention for the GEOMETRIC layer: a bump is counted
    as a peak only if its prominence >= this fraction of max chi. This governs
    the I/II/III verdict and is therefore a declared convention, NOT a
    paper-free absolute. Default 0.10; no paper clause
    fixes this value — it is a tunable convention, surfaced here so the verdict
    carries the parameter it depends on."""

    def as_dict(self) -> dict:
        return {
            "geometric": self.geometric,
            "spectral": self.spectral,
            "operational_floor_triggered": self.operational_floor_triggered,
            "degenerate_O": self.degenerate_O,
            "resolved_peak_count": self.resolved_peak_count,
            "peak_count_over_ceiling": self.peak_count_over_ceiling,
            "thresholds": {
                "eta_O": self.eta_O,
                "min_prominence_ratio": self.min_prominence_ratio,
                "max_resolved_peaks": self.max_resolved_peaks,
                "chi_abs_floor": self.chi_abs_floor,
            },
        }


# Convenience alias for users who don't want to import Trichotomy directly.
Regime = Trichotomy


# ---------------------------------------------------------------------------
# Blindness map / null cone — the programme's honesty anchor
# ("The Parrot's Theorems", sec:corank + the Layer table in sec:setup)
# ---------------------------------------------------------------------------
#
# WHY THIS EXISTS. The central Parrot contribution is not any single sigma_c
# number; it is the statement of what a given (object, instrument) pairing
# CANNOT see. The paper's sec:corank calls this "the program's honesty anchor":
# every instrument family splits the object into a visible part and a null cone
# (Kalman's observability decomposition; the count of blind directions is the
# corank, Sylvester's inertia). Until now that anchor was NOWHERE in the API.
#
# WHAT THIS IS — AND IS NOT. This map does NOT compute a corank number. The
# corank of a pairing needs the model class (the sensitivity matrix rows
# dG(s_i,.)/dtheta, sec:counting's rank law); this kernel is handed one
# observable scanned over one dial and does not have that matrix. Inventing a
# corank count from a single scan would be exactly the "story-not-science" the
# edition removes. So the map does two honest things instead:
#
#   (1) It classifies each quantity THIS result actually reported by its LAYER,
#       a proven lookup into the paper's Layer table (sec:setup), not a
#       computation. The layer fixes the invariance group, hence how much a
#       relabeling of the dial or a reprocessing of the reading can move the
#       number -- i.e. how trustworthy it is. Layer 1 (counts) is invariant
#       under all dial relabelings AND all invertible post-processing; Layer 2
#       (rates/gaps: tau) under all dial relabelings; Layer 3 (peak LOCATIONS:
#       sigma_c, gamma_c, ...) only under power laws sigma -> a*sigma^b, and so
#       is "spectrum x convention" -- the least trustworthy row.
#
#   (2) It carries the paper's meta-rule (sec:setup / worked in sec:staircase):
#       "layer membership is a property of invariance, not of the defining
#       formula." The geometric verdict is written as a COUNT (Layer-1 formula)
#       but its value depends on the declared prominence convention
#       min_prominence_ratio -- so its READING is Layer 3. The map records both,
#       so a count that secretly rides on a convention cannot pass as intrinsic.
#
#   (3) It states the corank hard edge verbatim in force: a degree of freedom
#       that neither turns with the dial nor rustles into the instrument is
#       invisible to THIS pairing forever, and no formalism repairs it -- only a
#       new instrument does. Whether a second, dial-free channel (the jitter /
#       two-probe channel, sec:jitter) was actually used decides whether the
#       "rustle" direction was even looked at.
#
# Citations render via theorem_map.cite(). The corank / Layer-table anchors live
# in *The Parrot's Theorems* (10.5281/zenodo.22066713), which the
# framework's THEOREM_MAP.md does not yet resolve (that re-anchoring is a
# separate, author-gated step); until then cite() falls back to `paper [label]`,
# which truthfully marks the anchor as not-yet-remapped rather than inventing a
# number.

# Trust ratings attached to each layer (paper sec:setup, "decreasing
# trustworthiness"). These are the paper's words, not a score.
_LAYER_TRUST = {
    1: "intrinsic to the pairing (most trustworthy, coarse)",
    2: "intrinsic to system + probe",
    3: "spectrum x convention (least trustworthy -- a reading, not an invariant)",
}
_LAYER_INVARIANT_UNDER = {
    1: "all dial relabelings AND all invertible post-processing of readings",
    2: "all dial relabelings",
    3: "only power laws sigma -> a*sigma^b",
}


@dataclass(frozen=True)
class LayerEntry:
    """One reported quantity, placed on the Parrot layer ladder (sec:setup).

    `layer_formula` is the layer the quantity's DEFINING FORMULA belongs to
    (a count is Layer 1, a rate Layer 2, a location Layer 3). `layer_reading`
    is the layer its VALUE actually inhabits once you account for the
    representative choices it depends on -- the paper's meta-rule. When the two
    differ (e.g. the geometric count, whose value rides on min_prominence_ratio),
    the reading is the honest one and `convention_deps` names why.
    """
    quantity: str
    layer_formula: int
    layer_reading: int
    invariant_under: str
    trust: str
    convention_deps: Tuple[str, ...] = ()
    note: str = ""

    def as_dict(self) -> dict:
        return {
            "quantity": self.quantity,
            "layer_formula": self.layer_formula,
            "layer_reading": self.layer_reading,
            "invariant_under": self.invariant_under,
            "trust": self.trust,
            "convention_deps": list(self.convention_deps),
            "note": self.note,
        }


@dataclass(frozen=True)
class BlindnessMap:
    """The null cone / blindness map for one Result -- the honesty anchor as
    an output field (paper sec:corank + Layer table sec:setup).

    Built by `Result.blindness()`; do not construct directly. It reports, per
    quantity this result carried, the layer / invariance group / trust, the
    convention dependencies the Layer-3 readouts ride on, and the structural
    null-cone warning. It deliberately reports NO corank number (see the module
    comment: the single-observable kernel does not hold the sensitivity matrix
    that a corank requires).
    """
    layers: Tuple[LayerEntry, ...]
    convention_dependencies: Tuple[str, ...]
    null_cone_note: str
    citations: Tuple[str, ...]

    # NOTE — no `second_channel_used` field on the map. A single dial scan does not
    # run the dial-free jitter channel (the experimental sigma_c.experimental.jitter), so on this path the boolean
    # would always be False and carry zero bits — decoration, the same class as an
    # unresolved citation. The fact lives in null_cone_note as prose, where it reads
    # as a warning, not a queryable capability flag. (The Result itself carries
    # second_channel_used for when a rustle reading is attached.)

    def as_dict(self) -> dict:
        return {
            "layers": [e.as_dict() for e in self.layers],
            "convention_dependencies": list(self.convention_dependencies),
            "null_cone_note": self.null_cone_note,
            "citations": list(self.citations),
        }

    def summary(self, ascii_only: bool = True) -> str:
        from sigma_c.theorem_map import cite

        bullet = "*" if ascii_only else "•"
        arrow = "->" if ascii_only else "→"
        warn = "[!]" if ascii_only else "⚠"

        lines: List[str] = ["Blindness map (null cone -- the honesty anchor):"]
        for e in self.layers:
            if e.layer_formula == e.layer_reading:
                layer_str = f"L{e.layer_reading}"
            else:
                # the meta-rule case: a formula-layer that reads one layer down
                layer_str = f"L{e.layer_formula}-formula/L{e.layer_reading}-reading"
            lines.append(f"  {bullet} {e.quantity:<16} [{layer_str}] {e.trust}")
            lines.append(f"      invariant under: {e.invariant_under}")
            if e.convention_deps:
                lines.append(
                    f"      {warn} value rides on: {', '.join(e.convention_deps)}"
                )
            if e.note:
                lines.append(f"      {e.note}")
        if self.convention_dependencies:
            lines.append(
                "  declared conventions in play: "
                + ", ".join(self.convention_dependencies)
            )
        lines.append(f"  null cone: {self.null_cone_note}")
        for label in self.citations:
            lines.append(f"      backed by {cite(label)}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main Result
# ---------------------------------------------------------------------------

@dataclass
class Result:
    """
    The disciplined-reader output. Cite: paper Def 2.2 (def:sigmac).

    sigma_c is None when chi_O has no interior maximum (paper Def 8.1's bottom
    value) — never NaN, never a thrown exception (Enforcement 4).

    For regime II, sigma_c is a list of peak locations.
    """
    # --- σ_c value(s) ---
    sigma_c: Optional[Union[float, List[float]]]
    """Scalar (regime I), list (regime II), or None (regime III). Cite:
    def:sigmac and def:Onice-all."""

    # --- The certified spectral abscissa (in TIME), when the tau-route ran ---
    tau_abscissa: Optional[float]
    """The certified live-route relaxation time tau = -T_star/log(lambda_2) (in TIME),
    from thm:abscissa, when a self-adjoint spectrum/operator was supplied. None on
    the sigma_c-only path (no spectrum) and when there is no decaying mode. This is
    the FALLEN bridge value
    sigma_c/rho_star now lives ONLY under `tau_bridge` (sigma-axis units). Read
    `tau_abscissa_status` (thm:abscissa) and `window_readability_status` (thm:auto)."""

    # --- Profile constant provenance ---
    rho_star: Optional[float]
    """The profile constant. None when sigma_c is None."""

    rho_star_source: str
    """Provenance of rho_star. Either "analytic:<window_name>" for
    the five canonical windows (Enforcement 2), or "fitted" for the
    legacy fit-based path (cite: def:operational-test)."""

    # --- Regime verdict ---
    regime: Trichotomy
    """The regime verdict: a peak-count + an isolated-gap boolean + a floor flag
    (two named convention-operations, NOT theorems; the §8 'trichotomy' labels are
    ARCHIVED). See the Trichotomy docstring."""

    # --- Stability indicator ---
    gamma_O: Optional[float]
    """Strict-SOC constant gamma_O = -d^2/dsigma^2 chi_O^2 at sigma_c.
    Low gamma_O ⟹ flat peak ⟹ noisy/regime-transition-zone reading.
    None when sigma_c is None. Cite: prop:stability."""

    # --- Reproduction recipe (Enforcement 8) ---
    smoothing: dict = field(default_factory=dict)
    """Smoothing parameters used in computing chi_O from samples.
    Keys: kernel, bandwidth, interpolant_order. Empty when input was
    analytical."""

    window: Optional[str] = None
    """The window name used (bare / gamma2 / … / log_gaussian). PROVENANCE: a
    recorded convention that sets rho_star, so it belongs in to_dict() (adding it
    later would be a schema break). Recoverable from rho_star_source too, but
    stated explicitly so a reproduction is unambiguous."""

    detrended: bool = False
    """Whether power-law detrending was applied (Enforcement 6).
    Off by default — never silently true. If true, the result is marked
    heuristic per paper §C.3."""

    framework: Optional[Framework] = None
    """The declared transfer-operator framework. None when domain is
    continuous-spectrum or unknown (then tau interpretation is
    dominant-scale-probe, not literal spectral gap; Principle 2)."""

    # --- field-found patches ---

    sigma_c_at_grid_boundary: bool = False
    """True iff the chi peak sits at the lowest or highest grid sample,
    so sub-grid quadratic refinement was not possible. Downstream
    cross-tests should treat boundary peaks as soft (the true peak may
    lie outside the scan range entirely). Field-found in
    a downstream-application Phase 1."""

    preprocessing_scale_equivariant: Optional[bool] = None
    """User declaration about the preprocessing operators applied to
    the observable BEFORE analyze() was called.

    - True  : every preprocessing step (filter, smoothing) is scale
              equivariant; rho_star analytic is applicable, A1 holds.
    - False : at least one preprocessing step carries an absolute
              time-scale (e.g. fixed-bandwidth Butterworth bandpass);
              the framework's analytic rho_star is NOT guaranteed valid
              -- rho_star_source becomes "exploratory:..." and
              .falsifiable returns False.
    - None  : not declared (backward-compatible default). The result
              is unchanged, with a single note appended urging
              the user to declare. A1 guard introduced after the
              a downstream-application Phase 1 finding."""

    # --- tau feasibility gate (Gamma_A) inputs, #1 ---
    gamma_A: Optional[float] = None
    """Contamination-to-signal ratio ||QA||/|c_2^A| of the probe (0 = perfectly
    faithful to the slow mode). Supplied by the caller; not inferred from a fit.
    Used with spectral_gap and the observation window T = max(sigma grid) to
    machine-check the tail-window condition T > t*(gamma_A, Delta) in tau_status.
    None when not supplied -> tau cannot be certified OK (NOT_IDENTIFIED)."""

    spectral_gap: Optional[float] = None
    """The spectral gap Delta = mu_2 - mu_3 > 0 (continuous-time rate), computed
    from a supplied gapped spectrum: (log|lambda_2| - log|lambda_3|)/T_star. None
    when no >=3-eigenvalue gapped spectrum was supplied. Feeds t*(gamma_A, Delta).
    For a complex-conjugate leading pair the gap is taken over MODULI after the
    pair (the pair is ONE mode)."""

    # --- live-route spectral tau (the certified value) + its declaration ------
    # These are populated ONLY on the tau-route (a spectrum or operator was
    # supplied). On the sigma_c-only path they stay None / False and result.tau
    # remains the bridge value exactly as before.

    sigma_axis: Optional[str] = None
    """Author-pinned axis declaration for the tau-route: one of 'evolution_time'
    (sigma is time in T_star units; the gate uses T=max(sigma-grid)), 'window_time'
    (sigma is time; the gate uses T_obs; bridge-vs-spectrum comparison is
    instrument-side and OFF), or 'other' (sigma is NOT time; the gate uses T_obs;
    NO bridge-vs-spectrum comparison --- different units). None when not declared:
    then tau is NOT_IDENTIFIED (there is no tau without the declaration)."""

    spectrum_supplied: bool = False
    """True iff a spectrum or operator was supplied (the tau-route ran). Lets the
    tau verdict distinguish 'undeclared sigma_axis on the tau-route' (-> refuse)
    from the plain sigma_c-only bridge path."""

    tau_bridge: Optional[float] = None
    """The FALLEN bridge value sigma_c / rho_star (sigma-axis units). Kept as a
    labelled diagnostic of how far the bridge lies from the certified spectral
    tau; compared to it only per the sigma_axis mode (see tau_bridge_spectral_ratio)."""

    tau_spectral: Optional[float] = None
    """The CERTIFIED live-route relaxation time tau = -T_star/log|lambda_2|
    (= -1/mu_2, in TIME). None when no decaying subdominant mode exists. When
    present it is the value carried in result.tau. Cite: thm:spectral-id (PROVED)
    --- the ONE place the code computes the live-route value, not the bridge."""

    tau_window_T: Optional[float] = None
    """The observation window T fed to the tail-window gate T > t*(Gamma_A,Delta).
    max(sigma-grid) under evolution_time; T_obs under window_time/other. Never
    max(sigma) outside evolution_time (that was the latent gate-leak)."""

    oscillation_period: Optional[float] = None
    """For a complex lambda_2: the oscillation period 2*pi*T_star/|arg lambda_2|,
    reported SEPARATELY, never folded into tau. None for a real lambda_2."""

    spectral_condition_number: Optional[float] = None
    """The eigenvalue condition number kappa_2 = 1/|y^H x| of lambda_2 (a PROXY
    for non-normality, not a theorem). Requires the operator (left+right
    eigenvectors); None when only eigenvalues were supplied."""

    spectral_defect: Optional[str] = None
    """None, or the reason the spectral tau is NOT_RESOLVABLE: 'jordan_block'
    (t^k prefactor, distorted on a finite window) or 'non_normal' (kappa_2 over
    the declared threshold, the transient dominates)."""

    spectrum_provenance: Optional[str] = None
    """'exact' (lattice) or 'estimated' (DMD-like). None on the sigma_c-only path."""

    tau_spectral_band: Optional[float] = None
    """For an ESTIMATED spectrum: the relative band delta_tau/tau =
    delta|lambda|/(|lambda|*|log|lambda||), which DIVERGES as |lambda_2| -> 1.
    None when the spectrum is exact or no eigenvalue error was supplied."""

    tau_bridge_spectral_ratio: Optional[float] = None
    """tau_bridge / tau_spectral, the labelled bridge-vs-spectrum discrepancy.
    Emitted ONLY under sigma_axis='evolution_time' (a literal comparison, same
    units). None under 'window_time' (instrument-side, comparison OFF) and
    'other' (different units --- no comparison)."""

    # --- Phase-2 seam (rides along; the channel itself is a later release) ----
    second_channel_used: bool = False
    """Whether a second, dial-free channel (the jitter channel, sec:jitter) was
    used to look at the 'rustle' direction. A plain analyze() scan does not use it,
    so it stays False; the channel itself is available as the experimental sigma_c.experimental.jitter
    (jitter_rank / pair_parity_rank / two_channel_fdt). The blindness map's
    null-cone note reads THIS field, so attaching a rustle reading updates the note."""

    # --- self-adjoint tau-route facts (from core.reversible) ------------------
    selfadjoint: Optional[str] = None
    """'verified' (operator checked self-adjoint in the declared inner product),
    'declared' (spectrum-only: asserted, not checkable), 'none' (operator supplied,
    NOT self-adjoint -> NOT_APPLICABLE), or None on the sigma_c-only path."""

    inner_product_declared: Optional[Any] = None
    """The declared inner product for the tau-route: 'euclidean', 'auto', or a
    pi-weight vector (stored as a tuple). None on the sigma_c-only path."""

    _rev: Optional[Any] = field(default=None, repr=False)
    """The ReversibleSpectrum facts (core.reversible), when the tau-route ran.
    The two tau verdicts and rate_error_bound are DERIVED from it, never stored."""

    # --- Falsifiability ---
    @property
    def falsifiable(self) -> bool:
        """
        True iff tau was obtained from an analytic rho_star (paper's
        non-circular reading) AND the user has not declared the
        preprocessing as carrying an absolute scale.

        A False preprocessing_scale_equivariant declaration
        downgrades the result to exploratory regardless of rho_star
        source -- the user has explicitly told the framework that A1
        is not satisfied at the observable construction layer.
        """
        if self.preprocessing_scale_equivariant is False:
            return False
        return self.rho_star_source.startswith("analytic:")

    @property
    def sigma_c_status(self) -> "Verdict":
        """Per-output 4-code verdict for sigma_c (OK / NOT_APPLICABLE /
        NOT_RESOLVABLE). A formalisation of the existing regime/flag verdict, not
        a new threshold --- see sigma_c.codes. (tau's verdict is added in a
        second step.)"""
        from sigma_c.codes import sigma_c_verdict
        return sigma_c_verdict(self)

    @property
    def tau_abscissa_status(self) -> "Verdict":
        """Per-output 4-code verdict for tau_abscissa (thm:abscissa, the asymptotic
        spectral abscissa). OK needs self-adjointness in the declared inner product
        + a discrete isolated mu_2 + faithfulness; NO strict gap required.
        Non-self-adjoint / non-embeddable / reducible -> NOT_APPLICABLE."""
        from sigma_c.codes import tau_abscissa_verdict
        return tau_abscissa_verdict(self)

    @property
    def window_readability_status(self) -> "Verdict":
        """Per-output 4-code verdict for reading the RATE in a finite window
        (thm:auto). Adds a strict gap + T > max(t*, 1/|mu_3|) to tau_abscissa's
        preconditions, and reports rate_error_bound. NOT_APPLICABLE inherits from
        tau_abscissa; spectrum-only (no probe) -> NOT_IDENTIFIED."""
        from sigma_c.codes import window_readability_verdict
        return window_readability_verdict(self)

    @property
    def rate_error_bound(self) -> Optional[float]:
        """The proven finite-window rate-error bound (thm:auto), RELATIVE to |mu_2|:
        (|mu_2|+|mu_3|)*Gamma_A^2*e^{-Delta*T} / |mu_2|. A REPORTED number (dominance
        != accuracy), None unless the tau-route produced mu_2, mu_3, Delta and a
        gamma_A + window T are present. Conservative: stays finite as Delta -> 0."""
        import math as _math
        rv = self._rev
        if (rv is None or getattr(rv, "mu2", None) is None
                or getattr(rv, "mu_next", None) is None or rv.delta is None
                or self.gamma_A is None or self.tau_window_T is None):
            return None
        mu2, mu3, delta, T = rv.mu2, rv.mu_next, rv.delta, self.tau_window_T
        if mu2 == 0:
            return None
        return ((abs(mu2) + abs(mu3)) * self.gamma_A ** 2
                * _math.exp(-delta * T)) / abs(mu2)

    @property
    def tau(self):
        """Not part of the API. `result.tau` would conflate the fallen bridge value
        with the certified abscissa. Use `tau_bridge` (sigma-units, load-dominant)
        or `tau_abscissa` (time, certified) instead."""
        raise AttributeError(
            "result.tau is not part of the API (it would conflate the fallen bridge with the "
            "certified abscissa). Use result.tau_bridge (sigma-units, load-dominant) "
            "or result.tau_abscissa (time, certified spectral value)."
        )

    @property
    def tau_status(self):
        """Not part of the API --- split into two honest verdicts."""
        raise AttributeError(
            "result.tau_status is not part of the API. Use result.tau_abscissa_status "
            "(thm:abscissa: is tau a certified operator property) and "
            "result.window_readability_status (thm:auto: is the rate readable in the "
            "observation window)."
        )

    # --- Diagnostics ---
    notes: List[str] = field(default_factory=list)
    """Human-readable diagnostic notes (low gamma_O warnings, threshold-flip
    warnings, etc.). Always present, may be empty."""

    # --- Citation: which theorem this output is anchored on ---
    citations: List[str] = field(default_factory=list)
    """List of paper theorem labels (THEOREM_MAP entries) that backed this
    result. The "every output cites a theorem" half of the perfekt-Definition."""

    # --- Visualization data ---
    _profile_sigma: Optional[Any] = field(default=None, repr=False)
    _profile_chi: Optional[Any] = field(default=None, repr=False)
    title: str = ""
    x_name: str = "sigma (resolution)"
    y_name: str = "chi_O(sigma)"

    # ------------------------------------------------------------------
    # Blindness map / null cone — the honesty anchor as an output field
    # (paper sec:corank + the Layer table in sec:setup)
    # ------------------------------------------------------------------
    def blindness(self) -> "BlindnessMap":
        """Return the null-cone / blindness map for this result.

        Derived purely from fields already on the result, so it can never be
        stale or forgotten in some construction path. See the BlindnessMap
        module comment for what this does and, deliberately, does not compute.
        """
        rg = self.regime
        layers: List[LayerEntry] = []

        # --- the geometric verdict: a COUNT (Layer-1 formula) whose value
        #     rides on the prominence convention -> Layer-3 READING. This is
        #     the paper's meta-rule worked example (sec:staircase).
        layers.append(
            LayerEntry(
                quantity="regime (count)",
                layer_formula=1,
                layer_reading=3,
                invariant_under=_LAYER_INVARIANT_UNDER[1],
                trust=_LAYER_TRUST[1],
                convention_deps=(f"min_prominence_ratio={rg.min_prominence_ratio:g}",),
                note=("a bump below the prominence convention is not counted, so "
                      "the I/II/III value is convention-dependent despite being a count"),
            )
        )

        # --- tau: Layer 2 (a rate/gap; how big), when present.
        if self.tau_abscissa is not None or self.tau_bridge is not None:
            layers.append(
                LayerEntry(
                    quantity="tau",
                    layer_formula=2,
                    layer_reading=2,
                    invariant_under=_LAYER_INVARIANT_UNDER[2],
                    trust=_LAYER_TRUST[2],
                )
            )

        # --- sigma_c: Layer 3 (a peak LOCATION; where), when present.
        if self.sigma_c is not None:
            deps: List[str] = [f"min_prominence_ratio={rg.min_prominence_ratio:g}"]
            if self.detrended:
                deps.append("power-law detrending applied")
            if self.smoothing:
                bw = self.smoothing.get("bandwidth")
                if bw is not None:
                    deps.append(f"smoothing bandwidth={bw}")
            layers.append(
                LayerEntry(
                    quantity="sigma_c (location)",
                    layer_formula=3,
                    layer_reading=3,
                    invariant_under=_LAYER_INVARIANT_UNDER[3],
                    trust=_LAYER_TRUST[3],
                    convention_deps=tuple(deps),
                    note=("a peak location is 'spectrum x convention'; it is a "
                          "reading, not an invariant of the pairing"),
                )
            )

        # --- gamma_O: a stability diagnostic, NOT a clean layer quantity. Its
        #     magnitude is a curvature w.r.t. the dial scale -> parametrization
        #     dependent (Layer-3 reading); only the qualitative flag is used.
        if self.gamma_O is not None:
            layers.append(
                LayerEntry(
                    quantity="gamma_O (SOC)",
                    layer_formula=2,
                    layer_reading=3,
                    invariant_under=_LAYER_INVARIANT_UNDER[3],
                    trust="diagnostic only; magnitude is parametrization-dependent",
                    convention_deps=("dial parametrization of sigma",),
                    note="only the low/high flag is used downstream, not the magnitude",
                )
            )

        # --- the declared conventions this result rides on, collected once.
        #     min_prominence_ratio is the ONE real remaining convention the
        #     verdict rides on. The former epsilon_c/delta_sep entries were
        #     removed: they advertised dead conventions the result
        #     did NOT ride on — epsilon_c was never used, delta_sep only split
        #     I_spec/II_spec, a distinction discarded downstream. A blindness map
        #     naming a dead convention is exactly the honesty defect the map
        #     exists to catch, appearing in the map itself.
        conv: List[str] = [f"min_prominence_ratio={rg.min_prominence_ratio:g}"]
        if rg.eta_O:
            conv.append(f"eta_O={rg.eta_O:g}")
        if self.detrended:
            conv.append("detrended")

        # --- the null-cone hard edge (paper sec:corank), instantiated for this
        #     single-observable pairing. A second, dial-free channel (jitter /
        #     two-probe, sec:jitter) is what looks at the "rustle" direction;
        #     this kernel path does not run one, so the prose says so. (No
        #     boolean flag for it -- see the BlindnessMap note: a constant-False
        #     field carries zero bits.)
        # The 'rustle' clause reads the second_channel_used FIELD (Phase-2 seam),
        # not hard-coded prose: when the jitter/two-probe channel lands and sets
        # the field True, this note updates itself instead of going stale.
        if self.second_channel_used:
            rustle_clause = (
                "A second, dial-free jitter/two-probe channel WAS used, so the "
                "'rustle' direction (the second moment) was looked at."
            )
        else:
            rustle_clause = (
                "This scan used no dial-free jitter channel "
                "(second_channel_used=False), so the 'rustle' direction was not "
                "looked at here (the channel is available separately as the experimental sigma_c.experimental.jitter: "
                "jitter_rank / two_channel_fdt)."
            )
        null_cone_note = (
            "this result came from a SINGLE observable scanned over one dial. "
            "By the corank theorem, any degree of freedom that neither turns "
            "with the dial nor rustles into this observable is invisible to this "
            "pairing -- and no reprocessing recovers it; only a new instrument "
            "(a second observable, or the dial-free jitter/two-probe channel) "
            "does. No corank NUMBER is reported here: counting blind directions "
            "needs the model's sensitivity matrix, which a single scan does not "
            "carry. " + rustle_clause
        )

        return BlindnessMap(
            layers=tuple(layers),
            convention_dependencies=tuple(conv),
            null_cone_note=null_cone_note,
            citations=("sec:corank", "sec:setup"),
        )

    def resolution_band(self) -> Optional[dict]:
        """The resolution/sensitivity band on a scalar sigma_c -- NOT a confidence
        interval (see core.resolution). DERIVED here from the retained chi profile
        and the declared min_prominence_ratio, never stored, so it cannot go stale.

        Returns None unless sigma_c is a finite scalar (regime I) and the profile
        was retained (analyze() retains it automatically). Keys: band_lo/band_hi
        (value interval from the local grid cell), significant_digits +
        sigma_c_display (digits coupled to the band), prominence_stable_range +
        flip_low/flip_high (the per-input range of min_prominence_ratio over which
        the regime is unchanged -- computed from THIS observable, not a hardcoded r).
        """
        if self._profile_sigma is None or self._profile_chi is None:
            return None
        from sigma_c.core.resolution import resolution_band as _band
        return _band(
            self._profile_sigma,
            self._profile_chi,
            self.sigma_c,
            self.regime.min_prominence_ratio,
            chi_abs_floor=self.regime.chi_abs_floor,
        )

    def convention_window_stable(self) -> Optional[bool]:
        """Is the resolved peak COUNT stable across the declared prominence window
        r in [r0/2, 2*r0]? Works for ANY regime (scalar OR vector sigma_c), unlike
        resolution_band() which is scalar-only. The load-bearing OK gate: a peak
        set that holds only under a narrow convention setting is not a measurement.
        None if the profile was not retained."""
        if self._profile_sigma is None or self._profile_chi is None:
            return None
        from sigma_c.core.resolution import peak_count_stable_in_window
        return peak_count_stable_in_window(
            self._profile_sigma,
            self._profile_chi,
            self.regime.min_prominence_ratio,
            chi_abs_floor=self.regime.chi_abs_floor,
        )

    def card(self, save_to: Optional[str] = None):
        """
        Render a visualization card (PNG if save_to is set; matplotlib Figure
        returned regardless).

        Needs matplotlib, which is the opt-in ``[plot]`` extra — the kernel and
        ``analyze()`` run headless without it. Install with
        ``pip install "sigma-c-framework[plot]"``; a clear ImportError names this
        if matplotlib is missing.
        """
        if self._profile_sigma is None or self._profile_chi is None:
            raise RuntimeError(
                "card() requires the chi profile. Use analyze() which "
                "retains it automatically."
            )
        from sigma_c.card import render
        return render(
            self,
            self._profile_sigma,
            self._profile_chi,
            title=self.title,
            x_name=self.x_name,
            y_name=self.y_name,
            save_to=save_to,
        )

    # ------------------------------------------------------------------
    # Structured serialization
    # ------------------------------------------------------------------
    def to_dict(self) -> dict:
        """Machine-readable serialization of this result.

        The blindness map and the two 4-code applicability statuses are DERIVED
        HERE, at serialization time, from the primary fields -- they are NEVER
        stored on the result. A serialized Result therefore cannot carry a
        safeguard frozen at analyze()-time that later went stale: change a
        convention (e.g. regime.min_prominence_ratio) and re-serialize, and both
        the blindness map and the statuses change with it. A blindness map that
        could itself go stale would be worse than none.
        """
        from sigma_c import __version__ as _lib_version
        return {
            # Provenance stamp: a stored result must be reproducible a year later.
            # schema_version is the shape of THIS dict; library_version is the
            # kernel that produced it. schema_version starts at 1 in 6.0 -- adding
            # it now (not later) avoids a silent schema break for stored results.
            "schema_version": 1,
            "library_version": _lib_version,
            "sigma_c": self.sigma_c,
            "tau_bridge": self.tau_bridge,
            "tau_abscissa": self.tau_abscissa,
            "selfadjoint": self.selfadjoint,
            "rate_error_bound": self.rate_error_bound,
            "rho_star": self.rho_star,
            "rho_star_source": self.rho_star_source,
            "window": self.window,
            "falsifiable": self.falsifiable,
            "gamma_O": self.gamma_O,
            "regime": self.regime.as_dict(),
            "framework": self.framework.value if self.framework is not None else None,
            "detrended": self.detrended,
            # --- live-route spectral tau declaration + diagnostics ---
            "sigma_axis": self.sigma_axis,
            "oscillation_period": self.oscillation_period,
            "spectral_condition_number": self.spectral_condition_number,
            "spectral_defect": self.spectral_defect,
            "spectrum_provenance": self.spectrum_provenance,
            "tau_spectral_band": self.tau_spectral_band,
            "tau_bridge_spectral_ratio": self.tau_bridge_spectral_ratio,
            "second_channel_used": self.second_channel_used,
            "citations": list(self.citations),
            "notes": list(self.notes),
            # --- derived safeguards: computed now from the fields, not stored ---
            "blindness": self.blindness().as_dict(),
            "sigma_c_status": self.sigma_c_status.as_dict(),
            "tau_abscissa_status": self.tau_abscissa_status.as_dict(),
            "window_readability_status": self.window_readability_status.as_dict(),
            "resolution_band": self.resolution_band(),
        }

    # ------------------------------------------------------------------
    # Display
    # ------------------------------------------------------------------
    def brief(self, ascii_only: bool = True) -> str:
        """A short summary -- alias for summary(short=True)."""
        return self.summary(ascii_only=ascii_only, short=True)

    def summary(self, ascii_only: bool = True, short: bool = False) -> str:
        """Human-readable summary -- what a Card shows in text form.

        ASCII by default for portability (Windows cp1252 console renders
        Unicode as ?). Set ascii_only=False for the prettier Unicode form.

        short=True gives a compact form (a few lines) that STILL carries the
        4-code status AND the resolution band -- there is no naked-number mode:
        every mode states what the number is worth.
        """
        from sigma_c.theorem_map import cite

        bot = "_|_" if ascii_only else "⊥"  # ⊥
        bullet = "*" if ascii_only else "•"  # •
        warn = "[!]" if ascii_only else "⚠"  # ⚠
        dash = "--" if ascii_only else "—"   # —

        if short:
            from sigma_c.codes import sigma_c_verdict
            v = sigma_c_verdict(self)
            out: List[str] = []
            # value line -- never a naked number: it sits next to the verdict.
            if self.sigma_c is None:
                val = f"sigma_c = {bot}"
            elif isinstance(self.sigma_c, list):
                val = f"sigma_c = [{len(self.sigma_c)} peaks] (regime II, vector)"
            else:
                val = f"sigma_c ~ {self.sigma_c:.4g}"
            out.append(f"VERDICT : {v.code.value}  {dash}  {val}")
            # band -- always shown when a scalar sigma_c has one.
            rb = self.resolution_band()
            if rb is not None:
                out.append(
                    f"band    : {rb['sigma_c_display']} +/- {rb['band_abs']:.2g} "
                    f"(resolution band, NOT a CI; {rb['significant_digits']} sig. digits)"
                )
            elif isinstance(self.sigma_c, list):
                out.append("band    : -- (vector-valued; band is defined for a single peak)")
            else:
                out.append("band    : -- (no resolved scalar to bound)")
            # one-line reason for a refusal, so short mode is still honest.
            if v.code.value != "OK":
                reason = v.reason.split(":")[0] if ":" in v.reason else v.reason
                out.append(f"why     : {reason[:72]}")
            out.append(f"(window={self.window}, min_prominence_ratio="
                       f"{self.regime.min_prominence_ratio:g}; full report: .summary())")
            return "\n".join(out)

        lines: List[str] = []
        lines.append("=" * 60)
        lines.append("sigma_c -- disciplined-reader output")
        lines.append("=" * 60)

        # Status headline: the plain 4-code verdict for sigma_c, up top, so the
        # reader sees the answer (OK / NOT_*) before the explanation.
        from sigma_c.codes import sigma_c_verdict
        _v = sigma_c_verdict(self)
        # Put the answer on the top line: verdict + the value it applies to.
        if self.sigma_c is None:
            _val = f"sigma_c = {bot}"
        elif isinstance(self.sigma_c, list):
            _val = f"sigma_c = [{len(self.sigma_c)} peaks]"
        else:
            _rb = self.resolution_band()
            _val = (f"sigma_c ~ {_rb['sigma_c_display']}" if _rb is not None
                    else f"sigma_c ~ {self.sigma_c:.4g}")
        lines.append(f"VERDICT       : {_v.code.value}  {dash}  {_val}")
        _reason = _v.reason.split(":")[0] if ":" in _v.reason else _v.reason
        if len(_reason) > 60:
            _reason = _reason[:57] + "..."
        lines.append(f"                {_reason}")
        lines.append("")

        # Plain-language header (#4): what sigma_c IS, in one breath, so the number
        # below is readable without the paper.
        lines.append(
            "what this reads: sigma_c is the scale on your dial where the observable"
        )
        lines.append(
            "is most susceptible to a change of scale (the susceptibility peak) -- a"
        )
        lines.append(
            "LOCATION read under a declared convention, not a fitted parameter or bound."
        )
        lines.append("")

        # sigma_c -- with its resolution band (NOT a CI) and band-coupled digits (#1/#2).
        if self.sigma_c is None:
            lines.append(f"sigma_c       : {bot} (no interior maximum -- regime III)")
        elif isinstance(self.sigma_c, list):
            vals = ", ".join(f"{v:.4g}" for v in self.sigma_c)
            lines.append(f"sigma_c       : [{vals}]  (regime II, vector-valued)")
        else:
            rb = self.resolution_band()
            if rb is not None:
                lines.append(
                    f"sigma_c       : {rb['sigma_c_display']} +/- {rb['band_abs']:.2g}"
                    f"  (resolution band, NOT a CI; {rb['significant_digits']} sig. digits)"
                )
                lo, hi = rb["prominence_stable_range"]
                flips = not (rb["prominence_flip_low"] is None
                             and rb["prominence_flip_high"] is None)
                tail = ("  (a bump crosses the counting threshold outside this range)"
                        if flips else "  (stable for every prominence convention)")
                lines.append(
                    f"                regime {self.regime.geometric} holds for prominence "
                    f"r in [{lo:.3g}, {hi:.3g}]{tail}"
                )
            else:
                lines.append(f"sigma_c       : {self.sigma_c:.4g}")

        # tau -- TWO honest fields, each with its verdict INLINE: the fallen
        # bridge (sigma-units, load-dominant) and, on the tau-route, the certified
        # spectral abscissa (time) with BOTH the abscissa (thm:abscissa) and the
        # window-readability (thm:auto) verdicts.
        if self.tau_bridge is not None:
            lines.append(
                f"tau_bridge    : {self.tau_bridge:.4g}  = sigma_c/rho_star "
                f"(load-dominant read-out, NOT the relaxation time)"
            )
        if self.tau_abscissa is not None:
            aa = self.tau_abscissa_status
            aflag = "" if aa.ok else f"  [{aa.code.value}]"
            lines.append(
                f"tau_abscissa  : {self.tau_abscissa:.4g}{aflag}  = -T*/log(lambda_2) "
                f"(certified spectral, in time; selfadjoint={self.selfadjoint})"
            )
            if not aa.ok:
                lines.append(f"                {warn} abscissa {aa.code.value}: {aa.reason}")
            wr = self.window_readability_status
            wflag = "" if wr.ok else f"  [{wr.code.value}]"
            reb = self.rate_error_bound
            reb_str = f"  (rel. rate-error bound {reb:.3g})" if reb is not None else ""
            lines.append(f"window_read.  : {wr.code.value}{wflag}{reb_str}")
            if not wr.ok:
                lines.append(f"                {warn} {wr.reason}")
                if wr.remediation:
                    lines.append(f"                -> to identify, supply: {wr.remediation}")
        elif self.tau_bridge is None:
            lines.append(f"tau           : {dash}")

        # rho_star + provenance
        if self.rho_star is not None:
            kind = "analytic" if self.falsifiable else "FITTED"
            lines.append(
                f"rho_star      : {self.rho_star:.4g}  "
                f"[{self.rho_star_source}, {kind}]"
            )
        else:
            lines.append(f"rho_star      : {dash}")

        # regime
        rg = self.regime
        spec_str = rg.spectral or dash
        floor_str = " (FLOOR triggered)" if rg.operational_floor_triggered else ""
        lines.append(
            f"regime        : geom={rg.geometric}, spec={spec_str}{floor_str}"
        )

        # stability
        if self.gamma_O is not None:
            gamma_warn = "  (LOW -- transition-zone)" if self.gamma_O < 0.1 else ""
            lines.append(f"gamma_O (SOC) : {self.gamma_O:.4g}{gamma_warn}")
        else:
            lines.append(f"gamma_O (SOC) : {dash}")

        # framework
        if self.framework is not None:
            lines.append(
                f"framework     : {self.framework.value} "
                f"({self.framework.reading_kind})"
            )
        else:
            lines.append("framework     : --  (dominant-scale-probe reading only)")

        # falsifiability flag (Enforcement 2) -- only meaningful when rho_star exists
        if self.rho_star is not None and not self.falsifiable:
            lines.append("                " + warn + " rho_star FITTED -- exploratory, not falsifiable")

        # detrending flag (Enforcement 6)
        if self.detrended:
            lines.append("                " + warn + " DETRENDED -- heuristic, see paper §C.3")

        # notes
        if self.notes:
            lines.append("")
            for note in self.notes:
                lines.append(f"  {bullet} {note}")

        # citations (the "every output cites a theorem" half)
        if self.citations:
            lines.append("")
            lines.append("Theorem backing:")
            for label in self.citations:
                lines.append(f"  {bullet} {cite(label)}")

        # blindness map / null cone (the honesty anchor -- paper sec:corank)
        lines.append("")
        lines.append(self.blindness().summary(ascii_only=ascii_only))

        # 4-code applicability legend (#4) -- so the [CODE] tags above can be decoded
        # from the output alone, without the paper.
        lines.append("")
        lines.append("what each verdict means (read it before you trust the number):")
        lines.append("  OK              the number stands.")
        lines.append("  NOT_RESOLVABLE  your data isn't enough for this.")
        lines.append("  NOT_APPLICABLE  wrong tool for this question.")
        lines.append("  NOT_IDENTIFIED  something from you is missing -- it says what to bring.")

        lines.append("=" * 60)
        return "\n".join(lines)

    def __repr__(self) -> str:
        rg = self.regime.geometric
        if self.sigma_c is None:
            return f"Result(sigma_c=None, regime={rg})"
        if isinstance(self.sigma_c, list):
            return f"Result(sigma_c={self.sigma_c}, regime={rg})"
        if self.tau_abscissa is not None:
            aa = self.tau_abscissa_status
            tau_repr = f"tau_abscissa={self.tau_abscissa:.4g}"
            if not aa.ok:
                tau_repr += f" [{aa.code.value}]"
        else:
            tau_repr = f"tau_bridge={self.tau_bridge!s:.6}" if self.tau_bridge is not None else "tau=--"
        return (
            f"Result(sigma_c={self.sigma_c:.4g}, {tau_repr}, "
            f"regime={rg}, source={self.rho_star_source})"
        )


# ---------------------------------------------------------------------------
# Two-probe test result (Enforcement 5)
# ---------------------------------------------------------------------------

@dataclass(frozen=True)
class TwoProbeResult:
    """
    Outcome of the non-circular two-probe test, paper Def 6.4
    (def:operational-test-noncirc).

    A failure carries the *cause* — Enforcement 5 — not just a False.
    """
    passed: bool
    delta: float
    """Measured |log(sigma_c1/rho_star_1) - log(sigma_c2/rho_star_2)|."""
    delta_threshold: float
    """The pre-declared tolerance the user supplied."""

    cause: Optional[str] = None
    """When passed=False, one of:
    - "regime_ii"           : probes resolved different moduli, multi-mode system
    - "faithfulness_break"  : at least one probe is not single-mode faithful
    - "system_disagreement" : probes disagree on what counts as 'same system'
    - "probes_not_distinct" : the two reported sigma_c values are closer than the
                              coarser scan's grid spacing -- a PRECONDITION
                              failure of the two-probe test, not a system
                              property. Field-found in a downstream-application Phase 1:
                              the user's observable construction collapsed two
                              declared windows onto one probe. Fix the
                              observable, not the system.
    - "indeterminate"       : test result inconclusive from data alone
    When passed=True, cause is None.
    """

    tau_1: Optional[float] = None
    tau_2: Optional[float] = None

    rho_star_source_1: str = ""
    rho_star_source_2: str = ""

    @property
    def both_analytic(self) -> bool:
        return (
            self.rho_star_source_1.startswith("analytic:")
            and self.rho_star_source_2.startswith("analytic:")
        )

    @property
    def tau_status(self) -> "Verdict":
        """Per-output 4-code verdict for tau from this two-probe test. A pass is
        OK *with the one-sided caveat* (never a bare OK: passing does not certify
        faithfulness); a deviation is NOT_IDENTIFIED. See sigma_c.codes."""
        from sigma_c.codes import tau_two_probe_verdict
        return tau_two_probe_verdict(self)

    def summary(self) -> str:
        from sigma_c.theorem_map import cite

        head = "PASSED" if self.passed else f"FAILED ({self.cause})"
        cite_str = cite(
            "def:operational-test-noncirc" if self.both_analytic
            else "def:operational-test"
        )
        lines = [
            "Two-probe test:",
            f"  status     : {head}",
            f"  delta      : {self.delta:.4g} (tolerance {self.delta_threshold:.4g})",
            f"  tau_1      : {self.tau_1!r}",
            f"  tau_2      : {self.tau_2!r}",
            f"  provenance : {self.rho_star_source_1} | {self.rho_star_source_2}",
            f"  backed by  : {cite_str}",
        ]
        return "\n".join(lines)

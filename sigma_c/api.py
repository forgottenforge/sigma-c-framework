# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Hero API — the disciplined-reader entry point.

    >>> from sigma_c import analyze
    >>> result = analyze(sigma, O_values)
    >>> print(result.summary())
    >>> result.card("out.png")

Cite (in docstrings throughout): paper labels via THEOREM_MAP.
"""
from __future__ import annotations
from typing import Callable, List, Optional, Sequence, Tuple, Union

import math
import numpy as np

from sigma_c.core.susceptibility import (
    chi_O,
    chi_O_from_callable,
    find_interior_maxima,
    quadratic_peak_in_log_sigma,
    SmoothingSpec,
)
from sigma_c.core.trichotomy import classify
from sigma_c.core.stability import compute_gamma_O
from sigma_c.framework import Framework
from sigma_c.result import Result, TwoProbeResult
from sigma_c.windows import Window, bare, resolve


_SIGMA_AXIS_MODES = ("evolution_time", "window_time", "other")


def analyze(
    sigma: Union[np.ndarray, Sequence[float]],
    O: Union[np.ndarray, Sequence[float], Callable[[np.ndarray], np.ndarray], None] = None,
    *,
    chi: Optional[Union[np.ndarray, Sequence[float]]] = None,
    window: Union[str, Window, None] = "bare",
    framework: Optional[Framework] = None,
    spectrum: Optional[Sequence[complex]] = None,
    operator: Optional[Sequence[Sequence[complex]]] = None,
    inner_product="euclidean",
    T_star: Optional[float] = None,
    sigma_axis: Optional[str] = None,
    T_obs: Optional[float] = None,
    spectrum_provenance: str = "exact",
    lambda2_abs_error: Optional[float] = None,
    non_normal_kappa_max: float = 1e3,
    eta_O: float = 0.0,
    min_prominence_ratio: float = 0.10,
    max_resolved_peaks: int = 5,
    label: str = "",
    preprocessing_scale_equivariant: Optional[bool] = None,
    gamma_A: Optional[float] = None,
) -> Result:
    """
    The hero call. Cite: def:sigmac, prop:structural-reduction. (The regime
    verdict is two named convention-operations, not a cited §8 theorem — the
    §8 "trichotomy" labels are ARCHIVED; see THEOREM_MAP.md.)

    Parameters
    ----------
    sigma : array
        Resolution scale samples, strictly positive.
    O : array or callable
        Observable values O(sigma) — either an array matching sigma or
        a callable.
    window : str or Window
        One of "bare", "gamma2", "gamma3", "exponential", "log_gaussian",
        or a Window instance. Cite: paper §4.2 window-family table.
        The window's rho_star is *analytical*, not fitted (Enforcement 2).
    framework : Framework or None
        The transfer-operator setting per paper Prop 5.4. None means
        continuous-spectrum / unknown-operator — dominant-scale-probe
        reading only (Principle 2).
    spectrum : sequence of complex or None
        Optional full spectrum {lambda_1, ...} for the isolated-leading-gap
        boolean (spectrum-only; not observable-faithful — faithfulness runs
        through gamma_A). Supplying it (or `operator`) invokes the tau-route,
        which then REQUIRES `T_star` and `sigma_axis`.
    operator : 2D array of complex or None
        Optional transfer/generator matrix. When supplied the spectrum is
        derived from it AND the non-normality proxy kappa_2 and the Jordan-block
        check become available (they need left+right eigenvectors). Give either
        `spectrum` or `operator`, not both.
    T_star : float or None
        Per-step time-scale; sets units for tau (tau in TIME). REQUIRED whenever
        the tau-route runs (a spectrum/operator is supplied) — there is no silent
        default 1.0 any more, since a default 1 silently turns per-step into
        "time". Not needed on the sigma_c-only path.
    sigma_axis : {"evolution_time", "window_time", "other"} or None
        MANDATORY declaration for the tau-route (Part B). evolution_time: sigma is
        time in T_star units, the gate uses T=max(sigma-grid) and the
        bridge-vs-spectrum comparison is literal. window_time: sigma is time, the
        gate uses T_obs (mandatory), the comparison is instrument-side and OFF.
        other: sigma is NOT time, the gate uses T_obs (mandatory), NO comparison
        (bridge in sigma-units, spectral tau in time). Without it there is NO tau:
        sigma_c is computed as usual but tau_status is NOT_IDENTIFIED.
    T_obs : float or None
        Observation/recording length. MANDATORY for sigma_axis in
        {window_time, other}; ignored (and max(sigma-grid) used) for
        evolution_time. No model unit-maps: only pure unit conversions may map
        sigma to time — converting a noise amplitude to time is a physical
        hypothesis (Kramers) and never lives here.
    spectrum_provenance : {"exact", "estimated"}
        exact = lattice/known operator; estimated = DMD-like. If estimated and
        lambda2_abs_error is given, the relative band delta_tau/tau is reported.
    lambda2_abs_error : float or None
        The absolute uncertainty delta|lambda_2| for an estimated spectrum; feeds
        the band delta_tau/tau = delta|lambda|/(|lambda|*|log|lambda||).
    non_normal_kappa_max : float
        DECLARED-convention threshold on the eigenvalue condition number kappa_2;
        above it the non-normal transient dominates and tau is NOT_RESOLVABLE.
        Default 1e3 (a convention, not a theorem).
    eta_O : float
        Signal/noise floor for the operational diagnostic (def:noise-floor-
        diagnostic). Default 0.0 = pure-mathematics limit.
    max_resolved_peaks : int
        DECLARED-convention ceiling on the resolved peak count. Above it the peak
        field is noise, not structure, and sigma_c is NOT_RESOLVABLE (the peaks are
        still listed). Default 5 (a convention, not a theorem). Independent of the
        constant-observable guard, which always fires (NOT_IDENTIFIED) when O has
        no variation above the numerical-noise floor.
    label : str
        Human-readable label for cards / reports.

    Returns
    -------
    Result
        See `sigma_c.result.Result`. sigma_c may be None (regime III),
        scalar (regime I), or list (regime II) — never NaN (Enforcement 4).
    """
    win = resolve(window)
    sigma_arr = np.asarray(sigma, dtype=float)

    # --- tau-route declaration validation (Parts A + B) ----------------------
    # The tau-route runs iff a spectrum OR an operator was supplied. On that route
    # T_star is REQUIRED (Part A: no silent default 1) and sigma_axis is a MANDATORY
    # declaration (Part B). Without sigma_axis there is no HARD error --- sigma_c is
    # computed as usual and only tau is refused (NOT_IDENTIFIED, Part B.3).
    tau_route = (spectrum is not None) or (operator is not None)
    if spectrum is not None and operator is not None:
        raise ValueError(
            "supply either spectrum (eigenvalues) or operator (the matrix) for the "
            "tau-route, not both: they would give two different spectra."
        )
    if tau_route:
        if T_star is None:
            raise ValueError(
                "T_star is required on the tau-route (a spectrum/operator was "
                "supplied): tau is in TIME and a silent default 1.0 would turn "
                "per-step into 'time' unannounced. Pass T_star explicitly."
            )
        if T_star <= 0:
            raise ValueError("T_star must be > 0.")
        if sigma_axis is not None and sigma_axis not in _SIGMA_AXIS_MODES:
            raise ValueError(
                f"sigma_axis must be one of {_SIGMA_AXIS_MODES} (or None to refuse "
                f"tau); got {sigma_axis!r}."
            )
        if sigma_axis in ("window_time", "other") and T_obs is None:
            raise ValueError(
                f"T_obs (the observation/recording length) is required for "
                f"sigma_axis={sigma_axis!r}: max(sigma-grid) is the window only "
                f"under evolution_time (the removed gate-leak)."
            )
        if spectrum_provenance not in ("exact", "estimated"):
            raise ValueError(
                "spectrum_provenance must be 'exact' or 'estimated'."
            )
        # P0 finite-input rule, EXTENDED to the tau-route. The kernel-wide guard
        # below covers sigma/O/chi only; a non-finite eigenvalue would sort to
        # lambda_1/lambda_2 and drive a spurious tau / OK verdict. Reject it here,
        # not read it as a spectrum. (np.isfinite on complex checks re AND im.)
        if spectrum is not None:
            _spec_arr = np.asarray([complex(x) for x in spectrum], dtype=complex)
            if not np.all(np.isfinite(_spec_arr)):
                raise ValueError(
                    "spectrum must be finite (no NaN/Inf eigenvalues): a non-finite "
                    "eigenvalue would sort to lambda_1/lambda_2 and drive a spurious "
                    "tau / OK verdict. Rejected, not read as a spectrum."
                )
        if operator is not None:
            _op_arr = np.asarray(operator, dtype=complex)
            if not np.all(np.isfinite(_op_arr)):
                raise ValueError(
                    "operator must be finite (no NaN/Inf entries): rejected, not "
                    "read as an operator (P0 finite-input rule)."
                )

    if chi is not None:
        # Caller supplied a precomputed chi profile (domain-specific
        # normalization, smoothed measurement, etc.). Honest data path:
        # use the chi values verbatim, do not re-derive them.
        chi_arr = np.asarray(chi, dtype=float)
        order = np.argsort(sigma_arr)
        sigma_grid = sigma_arr[order]
        chi = chi_arr[order]
        if O is not None and not callable(O):
            O_vals = np.asarray(O, dtype=float)[order]
        elif callable(O):
            O_vals = O(sigma_grid)
        else:
            # No O supplied; chi-only mode. We still need *some* O for the
            # signal/noise diagnostic; use chi itself as a proxy amplitude.
            O_vals = chi
        from sigma_c.core.susceptibility import SmoothingSpec
        smoothing = SmoothingSpec(kernel="precomputed",
                                  bandwidth=None,
                                  interpolant_order=0)
    elif callable(O):
        sigma_grid, chi, smoothing = chi_O_from_callable(O, sigma_arr)
        O_vals = O(sigma_grid)
    else:
        if O is None:
            from sigma_c.errors import InvalidInputError
            raise InvalidInputError("analyze() needs either O (observable) or chi.")
        O_vals = np.asarray(O, dtype=float)
        sigma_grid, chi, smoothing = chi_O(sigma_arr, O_vals)

    # Finite-input guard for ALL THREE paths (array-O, callable-O, precomputed-chi).
    # chi_O guards the array path at its own layer; this catches the callable-
    # evaluated O and the precomputed chi too, which bypass chi_O. A NaN/Inf anywhere
    # would propagate into a spurious sigma_c/OK verdict -- rejected (ValueError),
    # not read as a measurement. (the first fix guarded only chi_O.)
    if (not np.all(np.isfinite(sigma_grid)) or not np.all(np.isfinite(chi))
            or not np.all(np.isfinite(O_vals))):
        from sigma_c.errors import InvalidInputError
        raise InvalidInputError(
            "sigma, O (or the callable's output), and chi must all be finite "
            "(no NaN/Inf): a non-finite value would propagate into a spurious "
            "sigma_c / OK verdict. Rejected, not read as a measurement."
        )

    citations: List[str] = ["def:sigmac", "prop:structural-reduction"]
    notes: List[str] = []

    # --- live-route spectral analysis (the certified tau) --------------------
    # The ONE place the code computes tau from the spectrum, with the P2.2 gate-
    # hardening branches (complex pair / Jordan / non-normal / provenance band).
    # Delta = mu_2 - mu_3 (over MODULI, after any conjugate pair) comes from here,
    # replacing the old crude descending-sort that returned None for a complex
    # pair (equal moduli).
    # The certified route runs through the REVERSIBLE / self-adjoint engine
    # (core.reversible), NOT the old euclidean spectral_tau. Self-adjointness in the
    # DECLARED inner product (reversibility) is the theorem's precondition; Jordan /
    # near-defects cannot occur on this route (self-adjoint => real spectrum +
    # orthonormal eigenbasis + always semisimple).
    rev = None
    tau_abscissa_val: Optional[float] = None
    selfadjoint_val: Optional[str] = None
    classify_spectrum = spectrum
    spectral_gap_val: Optional[float] = None
    if tau_route:
        from sigma_c.core.reversible import analyze_reversible
        rev = analyze_reversible(
            operator=operator, spectrum=spectrum, T_star=T_star,
            inner_product=inner_product,
            selfadjoint_declared=(spectrum is not None),
        )
        selfadjoint_val = rev.selfadjoint
        tau_abscissa_val = rev.tau_spectral
        spectral_gap_val = rev.delta
        # For the gap boolean, classify() needs the eigenvalue list (real, from eigh
        # on the symmetrised operator). None when the operator is not self-adjoint.
        if rev.eigenvalues is not None:
            classify_spectrum = list(rev.eigenvalues)
        for _n in rev.notes:
            notes.append(_n)

    # Classify: peak-count operation + spectral gap-boolean + floor.
    # (No thm:trichotomy-geometric / thm:trichotomy-spectral citations: those
    # §8 labels are ARCHIVED, not live theorems -- the code runs two named
    # convention-operations, see THEOREM_MAP.md.)
    regime, peak_indices = classify(
        sigma_grid, chi, O_vals,
        spectrum=classify_spectrum,
        eta_O=eta_O,
        min_prominence_ratio=min_prominence_ratio,
        max_resolved_peaks=max_resolved_peaks,
    )
    if eta_O > 0:
        citations.append("def:noise-floor-diagnostic")
    # Framework declaration cites ONLY the framework taxonomy
    # (prop:standard-frameworks). It does NOT cite the PROVED thm:spectral-id:
    # that live-route theorem (tau = spectral abscissa) is not what this code
    # computes. result.tau is the fallen sigma_c=rho_star*tau bridge, cited at the
    # tau emission as thm:spectral-id-B (GAP-KNOWN). A PROVED label renders no
    # warning, so citing it on a value it does not back is a silent misanchoring
    # -- removed. Re-add thm:spectral-id here only when the code
    # actually derives tau via the live route (spectral abscissa from a spectrum).
    if framework is not None:
        if "prop:standard-frameworks" not in citations:
            citations.append("prop:standard-frameworks")

    # The observation window T fed to the tail-window gate: max(sigma-grid) ONLY
    # under evolution_time; T_obs under window_time/other. This is the removed
    # gate-leak (max(sigma) is the window only when sigma IS evolution time).
    tau_window_T: Optional[float] = None
    if tau_route and sigma_axis is not None:
        if sigma_axis == "evolution_time":
            tau_window_T = float(np.max(sigma_grid))
        else:  # window_time / other -> T_obs (validated present above)
            tau_window_T = float(T_obs)

    # Regime-independent tau-route declaration + spectral diagnostics, spread into
    # every return so the tau verdict and the serialized dict see them uniformly.
    spec_fields = dict(
        sigma_axis=sigma_axis if tau_route else None,
        spectrum_supplied=tau_route,
        tau_window_T=tau_window_T,
        spectrum_provenance=(spectrum_provenance if tau_route else None),
        selfadjoint=selfadjoint_val,
        inner_product_declared=(
            (inner_product if isinstance(inner_product, str) else "pi-weight")
            if tau_route else None
        ),
        _rev=rev,
        tau_spectral_band=(rev.band_rel if rev is not None else None),
    )

    # No interior peak → regime III, sigma_c = None (Enforcement 4)
    if len(peak_indices) == 0:
        citations.append("thm:diagnostic")
        if regime.degenerate_O:
            iii_note = (
                "Observable has no variation above the numerical-noise floor "
                "(constant to rounding precision): chi_O = |dO/dlog sigma| is 0 up "
                "to float error, so there is no scale to resolve. sigma_c is "
                "NOT_IDENTIFIED (supply an observable that actually varies with the "
                "dial), NOT a spurious peak read off rounding noise."
            )
        else:
            iii_note = "No interior peak -- regime III (sigma_c = bottom)."
        return Result(
            sigma_c=None, tau_abscissa=tau_abscissa_val,
            rho_star=None, rho_star_source="--",
            window=win.name,
            regime=regime,
            gamma_O=None,
            smoothing=smoothing.as_dict(),
            framework=framework,
            preprocessing_scale_equivariant=preprocessing_scale_equivariant,
            gamma_A=gamma_A, spectral_gap=spectral_gap_val,
            **spec_fields,
            notes=notes + [iii_note],
            citations=citations,
            _profile_sigma=sigma_grid,
            _profile_chi=chi,
            title=label,
        )

    # Stability check on each peak; warn if low gamma_O
    gamma_vals: List[Optional[float]] = [
        compute_gamma_O(sigma_grid, chi, idx) for idx in peak_indices
    ]
    citations.append("prop:stability")

    # Regime II (multi-mode, vector sigma_c)
    if len(peak_indices) >= 2:
        # field-finding 1: sub-grid refinement also for regime II.
        sigma_c_vec_refined: List[float] = []
        any_at_boundary = False
        for idx in peak_indices:
            sc_sub, at_bnd = quadratic_peak_in_log_sigma(sigma_grid, chi, idx)
            sigma_c_vec_refined.append(sc_sub)
            any_at_boundary = any_at_boundary or at_bnd
        gamma_min = min((g for g in gamma_vals if g is not None), default=None)
        notes.append(
            f"Multi-mode (regime II): {len(sigma_c_vec_refined)} resolved peaks "
            "(peak-count operation under min_prominence_ratio). sigma_c is "
            "vector-valued."
        )
        if regime.peak_count_over_ceiling:
            notes.append(
                f"Resolved peak count {regime.resolved_peak_count} exceeds the "
                f"declared ceiling max_resolved_peaks={regime.max_resolved_peaks}: "
                "this many 'modes' off one dial is noise, not structure. sigma_c is "
                "NOT_RESOLVABLE -- smooth the observable or declare a coarser "
                "resolution (raise min_prominence_ratio). The peaks are still listed."
            )
        return Result(
            sigma_c=sigma_c_vec_refined, tau_abscissa=tau_abscissa_val,
            rho_star=win.rho_star, rho_star_source=win.rho_star_source,
            window=win.name,
            regime=regime,
            gamma_O=gamma_min,
            smoothing=smoothing.as_dict(),
            framework=framework,
            sigma_c_at_grid_boundary=any_at_boundary,
            preprocessing_scale_equivariant=preprocessing_scale_equivariant,
            gamma_A=gamma_A, spectral_gap=spectral_gap_val,
            **spec_fields,
            notes=notes,
            citations=citations,
            _profile_sigma=sigma_grid,
            _profile_chi=chi,
            title=label,
        )

    # Regime I (single mode)
    peak_idx = peak_indices[0]
    # field-finding 1: sub-grid quadratic refinement of sigma_c.
    sigma_c_val, at_boundary = quadratic_peak_in_log_sigma(
        sigma_grid, chi, peak_idx,
    )
    # The FALLEN bridge value (sigma-axis units) -- kept ONLY as result.tau_bridge, a
    # labelled diagnostic. The certified spectral abscissa (in TIME) is tau_abscissa_val
    # from the reversible engine; it is a DIFFERENT object, not the bridge.
    tau_bridge_val = sigma_c_val / win.rho_star if win.rho_star > 0 else None
    gamma_val = gamma_vals[0]

    # bridge-vs-abscissa diagnostic: only under evolution_time (both then in T_star-time).
    tau_bridge_spectral_ratio = None
    if (tau_route and tau_abscissa_val is not None and tau_bridge_val is not None
            and sigma_axis == "evolution_time"):
        tau_bridge_spectral_ratio = tau_bridge_val / tau_abscissa_val
        notes.append(
            f"bridge-vs-abscissa (evolution_time, literal): tau_bridge="
            f"{tau_bridge_val:.4g} vs tau_abscissa={tau_abscissa_val:.4g} "
            f"(ratio {tau_bridge_spectral_ratio:.4g}) -- how far the fallen bridge lies "
            f"from the certified abscissa."
        )

    if tau_route:
        # thm:spectral-id (PROVED) backs the abscissa VALUE; the two verdicts
        # (tau_abscissa_status / window_readability_status) gate its trust. The fallen
        # bridge keeps its GAP-KNOWN thm:spectral-id-B flag.
        if "thm:spectral-id" not in citations:
            citations.append("thm:spectral-id")
        citations.append("thm:spectral-id-B")
        notes.append(
            "tau_abscissa = -T_star/log(lambda_2) (in time) is the certified spectral "
            "abscissa (thm:spectral-id, self-adjoint route); tau_bridge = sigma_c/rho_star "
            "is the fallen load-dominant read-out. Read tau_abscissa_status (is tau a "
            "certified operator property) and window_readability_status (is the rate "
            "readable in the observation window)."
        )
    elif tau_bridge_val is not None:
        # sigma_c-only path: only the fallen bridge exists (no spectrum/operator).
        citations.append("thm:spectral-id-B")
        notes.append(
            "result.tau_bridge = sigma_c / rho_star is the LOAD-DOMINANT read-out, not "
            "theorem-backed as the relaxation time: the sigma_c = rho_star*tau "
            "bridge holds only under single-mode faithfulness (THEOREM_MAP "
            "thm:spectral-id-B, GAP-KNOWN). tau_abscissa is NOT_IDENTIFIED here (no "
            "spectrum/operator); supply one, or run two_probe_test(r1, r2)."
        )

    citations.append("thm:cross-obs-concentration")
    if framework is not None and framework.is_experimental:
        notes.append(
            "framework=anisotropic_banach is experimental; spectrum is "
            "functional-space-dependent (paper Prop 5.4 case 6)."
        )

    # NOTE (honesty): the 0.1 cutoff is a DECLARED CONVENTION, not a
    # paper-fixed value. Prop 12.1 motivates "low gamma_O => transition zone";
    # the specific 0.1 has no paper clause and is a tunable threshold.
    if gamma_val is not None and gamma_val < 0.1:
        notes.append(
            f"Low gamma_O = {gamma_val:.3g} -- peak is flat, reading is in "
            f"the regime-transition zone (cf. paper Prop 12.1; the 0.1 cutoff "
            f"is a declared convention, not a paper value)."
        )
        citations.append("prop:transition-zone")

    if regime.operational_floor_triggered:
        notes.append(
            f"Signal/noise floor eta_O={eta_O:.3g} was exceeded; the geometric "
            "verdict is below the operational discrimination threshold "
            "(paper Def 8.8). Operationally classify as regime III."
        )

    if at_boundary:
        notes.append(
            "sigma_c sits at the scan-grid boundary; sub-grid refinement was "
            "not applied. The true peak may lie outside the scan range -- "
            "widen the sigma grid before relying on this value."
        )

    # field-finding 3: A1 preprocessing guard.
    rho_star_source = win.rho_star_source
    if preprocessing_scale_equivariant is False:
        # User declared a fixed-scale preprocessing operator. The framework's
        # analytic rho_star is no longer the right anchor; surface that.
        rho_star_source = (
            f"exploratory:absolute_scale_preprocessing(was {win.rho_star_source})"
        )
        notes.append(
            "preprocessing_scale_equivariant=False: the observable was built "
            "from a chain that carries an absolute time-scale. Analytic "
            "rho_star is NOT applicable (paper A1 axiom violated at the "
            "observable layer). Result is exploratory; .falsifiable returns "
            "False; cross-window/cross-detector two-probe results are "
            "uninterpretable as Def 6.7 falsification tests."
        )
    elif preprocessing_scale_equivariant is None:
        notes.append(
            "preprocessing_scale_equivariant not declared. The framework "
            "assumes scale-equivariant preprocessing (A1 holds) at the "
            "observable layer. Declare True/False explicitly to silence "
            "this note. Introduced after a downstream-application Phase 1 "
            "A1 finding."
        )

    return Result(
        sigma_c=sigma_c_val, tau_abscissa=tau_abscissa_val,
        tau_bridge=tau_bridge_val,
        tau_bridge_spectral_ratio=tau_bridge_spectral_ratio,
        rho_star=win.rho_star, rho_star_source=rho_star_source,
        window=win.name,
        regime=regime,
        gamma_O=gamma_val,
        smoothing=smoothing.as_dict(),
        framework=framework,
        sigma_c_at_grid_boundary=at_boundary,
        preprocessing_scale_equivariant=preprocessing_scale_equivariant,
        gamma_A=gamma_A, spectral_gap=spectral_gap_val,
        **spec_fields,
        notes=notes,
        citations=citations,
        _profile_sigma=sigma_grid,
        _profile_chi=chi,
        title=label,
    )


# ---------------------------------------------------------------------------
# Two-probe test — Enforcement 5 cause-branching
# ---------------------------------------------------------------------------

def two_probe_test(
    result_1: Result,
    result_2: Result,
    *,
    delta_threshold: float = 0.20,
    label_1: str = "probe 1",
    label_2: str = "probe 2",
) -> TwoProbeResult:
    """
    Non-circular operational two-probe test — paper Def 6.4
    (def:operational-test-noncirc).

    Both Result objects must have rho_star_source starting with "analytic:"
    for the test to be falsifiable; otherwise it is a fit-based legacy test
    (def:operational-test) and is marked exploratory.

    Failure carries the cause (Enforcement 5).
    """
    # Coerce to floats / handle None
    if result_1.sigma_c is None or result_2.sigma_c is None:
        # One is in regime III — at least one probe doesn't see a scale.
        s1 = result_1.sigma_c if result_1.sigma_c is not None else float("nan")
        s2 = result_2.sigma_c if result_2.sigma_c is not None else float("nan")
        return TwoProbeResult(
            passed=False,
            delta=float("nan"),
            delta_threshold=delta_threshold,
            cause="indeterminate",
            tau_1=result_1.tau_bridge, tau_2=result_2.tau_bridge,
            rho_star_source_1=result_1.rho_star_source,
            rho_star_source_2=result_2.rho_star_source,
        )

    if isinstance(result_1.sigma_c, list) or isinstance(result_2.sigma_c, list):
        return TwoProbeResult(
            passed=False, delta=float("nan"),
            delta_threshold=delta_threshold,
            cause="regime_ii",
            tau_1=result_1.tau_bridge, tau_2=result_2.tau_bridge,
            rho_star_source_1=result_1.rho_star_source,
            rho_star_source_2=result_2.rho_star_source,
        )

    if (result_1.rho_star is None or result_2.rho_star is None
            or result_1.rho_star <= 0 or result_2.rho_star <= 0):
        return TwoProbeResult(
            passed=False, delta=float("nan"),
            delta_threshold=delta_threshold,
            cause="indeterminate",
            tau_1=result_1.tau_bridge, tau_2=result_2.tau_bridge,
            rho_star_source_1=result_1.rho_star_source,
            rho_star_source_2=result_2.rho_star_source,
        )

    # field-finding 2: probes_not_distinct precheck.
    # When both reported sigma_c values are closer than the coarser scan's
    # grid step, the two probes are not actually distinct -- a precondition
    # failure of the two-probe test, not a system faithfulness break.
    sg1 = getattr(result_1, "_profile_sigma", None)
    sg2 = getattr(result_2, "_profile_sigma", None)
    log_distance_sigma_c = abs(math.log(result_1.sigma_c) - math.log(result_2.sigma_c))
    if sg1 is not None and sg2 is not None and len(sg1) > 1 and len(sg2) > 1:
        log_step_1 = math.log(sg1[1]) - math.log(sg1[0])
        log_step_2 = math.log(sg2[1]) - math.log(sg2[0])
        log_step = max(abs(log_step_1), abs(log_step_2))
        if log_distance_sigma_c < log_step:
            return TwoProbeResult(
                passed=False,
                delta=float(log_distance_sigma_c),
                delta_threshold=delta_threshold,
                cause="probes_not_distinct",
                tau_1=result_1.tau_bridge, tau_2=result_2.tau_bridge,
                rho_star_source_1=result_1.rho_star_source,
                rho_star_source_2=result_2.rho_star_source,
            )

    log_tau_1 = math.log(result_1.sigma_c) - math.log(result_1.rho_star)
    log_tau_2 = math.log(result_2.sigma_c) - math.log(result_2.rho_star)
    delta = abs(log_tau_1 - log_tau_2)

    passed = delta <= delta_threshold
    cause: Optional[str] = None
    if not passed:
        # Failure cause inference (Enforcement 5)
        # If both probes produced regime-I geometric verdicts but disagree
        # on tau by a large amount, it's most likely a faithfulness break or
        # a system-disagreement. With only sigma_c info we cannot perfectly
        # distinguish (a) regime-II coverage from (b) faithfulness, so we
        # mark it as "faithfulness_break" by default and let the user
        # supply more data to refine.
        # NOTE (honesty): the 2.0 boundary is a DECLARED CONVENTION,
        # not a paper-anchored cutoff. It names a cause branch and has no theorem
        # behind the specific value; treat the cause label as heuristic.
        if delta > 2.0:
            cause = "system_disagreement"
        else:
            cause = "faithfulness_break"

    return TwoProbeResult(
        passed=passed,
        delta=delta,
        delta_threshold=delta_threshold,
        cause=cause,
        tau_1=math.exp(log_tau_1),
        tau_2=math.exp(log_tau_2),
        rho_star_source_1=result_1.rho_star_source,
        rho_star_source_2=result_2.rho_star_source,
    )

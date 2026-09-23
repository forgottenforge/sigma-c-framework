# THEOREM_MAP — sigma_c ↔ paper labels

This file is the **single source of truth** binding the kernel's docstrings,
golden tests, and API surfaces to the paper. The code cites **labels**, never
numbers; this file holds the label → number mapping and, separately, a
**proof-status register** so that `cite()` can render a proof's provenance
(not only its location).

When a paper is renumbered (camera-ready, revisions), **update this file only** —
the rest of the code keeps its label citations and follows the new numbers here.

## Paper anchors

The code cites stable LaTeX labels; the numbers below resolve each label against
a fixed numbering scheme, so the code never hard-codes a number. The
results live in *The Parrot's Theorems* (a Zenodo preprint); each result's
status is marked in the register below, and labels will re-anchor to a published
home when the journal publication appears.

- **The Parrot's Theorems (a Zenodo preprint):** Zenodo
  `10.5281/zenodo.22066713`.
- **Maintainer:** ForgottenForge · <nfo@forgottenforge.xyz> · <https://forgottenforge.xyz>

## Mapping table

Stable LaTeX labels are listed verbatim from `paper.aux` of the submission
build. The third column points to the prescription item the label informs
(the third column names the design item each label informs).

### Section 2 — Setup

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `def:Onice` | 2.1 | Admissible observable class | core domain type |
| `def:sigmac` | 2.2 | Susceptibility-peak functional | engine entry point |
| `def:Onice-window` | 2.5 | Windowed admissible class | item G |
| `prop:window-bridge` | 2.6 | Windowed bridge | item G (lift theorem) |
| `def:profile-dec` | 2.8 | Profile decomposition | core profile type |

### Section 3 — Axioms

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `thm:compat` | 3.7 | σ_c satisfies (A1),(A2),(A4),(A5) — FOUR axioms; NOT the locality axiom A3; unconditional | item I (load-bearing direction) |
| `prop:axiomatic-char` | 3.8 | Conditional canonicity under (R†) | item I (out of scope) |

> **Honesty note.** `thm:compat` proves compatibility with **four**
> axioms, (A1),(A2),(A4),(A5); the paper states explicitly that σ_c
> does *not* satisfy the locality axiom (A3). The **five-axiom UNIQUENESS** is a
> *different* result — `thm:characterisation` in a separate, superseded axiomatic manuscript, NOT `thm:compat`. Do not describe
> `thm:compat` as "the five axioms": that conflates a four-axiom compatibility
> theorem with a five-axiom uniqueness theorem in another manuscript, and claims
> a published five-axiom result that does not exist.

## Proof-status register

Provenance applies to the proofs, not only the labels. A label may resolve to a
paper location and still stand on a proof with a known gap; if `cite()` rendered
only the number, a refuted proof would travel as an anchor through every output
until someone read the source. So proofs carry one of the states listed below,
and `cite()` renders anything other than `PROVED`:

- `PROVED` — proof checked by us against the certified core, no known gap. **The
  only status that renders nothing.** Not the default (see below).
- `PUBLISHED` — peer-reviewed and published elsewhere, but NOT re-derived against
  our chain. A real check status, a *different* one than `PROVED`; it renders, so
  the register says which kind of check exists rather than flattening peer-review
  into our own verification.
- `PREPRINT` — stated in a preprint (not yet peer-reviewed, and not
  re-derived against our chain). Renders, so the register never claims a
  peer-review it does not have.
- `GAP-KNOWN` — a specific, located error in the proof as written; the statement
  may still be true, but it is not established.
- `UNDER-REPAIR` — a `GAP-KNOWN` proof with a repair in progress; still not
  established until the repair passes an independent check.
- `SUPERSEDED` — a statement/route replaced by a corrected one elsewhere; must
  not be cited as the current result.
- `REVIEW-CLEAN` — independent cold reads find no remaining gap, but the result
  is NOT yet certified: it is not `PROVED` until a reader of genuinely different
  provenance signs off (the cold readers so far may share a blind spot). Do not
  cite as established.
- `NO-REGISTER-ENTRY` — **the default for any label with no status row.** "Nobody
  looked." Not "probably wrong" — just unlisted/unexamined, which is equally true of a
  definition (it can be inconsistent, vacuous, or double-meaning'd) as of a
  theorem. Renders a marker so the gap is visible. (Renamed from
  `UNVERIFIED`, which an outside user read as "the math is unverified"; the new marker
  says plainly "unlisted, not disproven".)
- `UNVERIFIED` — a valid EXPLICIT status (someone looked, no verdict recorded yet),
  distinct from the NO-REGISTER-ENTRY default. Also renders a marker.
- `ARCHIVED` — a paper-only formalism moved to the math archive
  (kept privately, not shipped with the package) because the live code no longer stands on the labelled
  theorem; it uses named operations instead. Not a proof-quality verdict (not
  `PROVED`/`GAP-KNOWN`): it records that the label is no longer a live theorem
  the code cites. Renders like any other non-`PROVED` status so no output
  silently anchors on an archived label.

`cite()` renders, for any non-`PROVED` label, the status **and** the short reason
in the third column (so a stale, superseded, unverified, not-yet-certified, or
archived anchor travels flagged). Status token must be one of the states above.

| Label | Status | Short reason (rendered by cite) |
|---|---|---|
| `thm:characterisation` | SUPERSEDED | five-axiom/A3 route replaced by `thm:axiomatic-char` (real-analytic, minimal {A1,A2,A4}); do not cite |
| `thm:compat` | PROVED | sigma_c satisfies (A1),(A2),(A4),(A5) and explicitly NOT locality A3 (unconditional chain-rule check) |
| `thm:noninv` | PREPRINT | Non-invariance of peak locations -- stated in The Parrot's Theorems (Zenodo preprint 22066713), not independently re-derived here |
| `prop:profile` | PREPRINT | Profile decomposition sigma_c=rho_star*tau (single-mode) -- stated in The Parrot's Theorems (Zenodo preprint 22066713), not independently re-derived here |
| `thm:spectral-id` | PROVED | tau = spectral abscissa via the two-probe live route (spectral-abscissa + gap-criterion); the sigma_c=rho_star*tau bridge is a SEPARATE claim that does not close (see README, the two tau fields) |
| `thm:spectral-id-A` | GAP-KNOWN | residual-perturbation lemma false in the stated R_{>0} generality; fixable on a compact window |
| `thm:spectral-id-B` | GAP-KNOWN | operator/window mechanism: centering makes chi bimodal (two equal peaks) so sigma_c is ill-defined; holds only on the monotone single-mode branch |
| `thm:cross-obs-concentration` | REVIEW-CLEAN | multi-mode concentration sigma_c=rho_star*tau*(1+O(r^n)) on the live two-probe route; independent cold reads clean, awaits a certifying read |
| `thm:axiomatic-char` | REVIEW-CLEAN | uniqueness from the minimal {A1,A2,A4} for real-analytic scores; cold reads clean, awaits a certifying read; the bare-C3 case is open |
| `thm:trichotomy-geometric` | ARCHIVED | paper-only formalism (moved to the math archive); the live code runs a named peak-count operation under the declared min_prominence_ratio convention |
| `thm:trichotomy-spectral` | ARCHIVED | paper-only formalism (math archive); the live code runs a named isolated-leading-gap boolean (spectrum-only, not observable-faithful) |
| `thm:trichotomy` | ARCHIVED | paper-only umbrella (math archive); the live code rests on the two named convention-operations (peak-count + gap-boolean), not this label |
| `thm:multimode` | ARCHIVED | paper-only formalism (math archive); regime II is decided by the geometric peak-count, never spectrally k-matched |

## Certification model

A proof's register entry has an **author** side and a **checker** side. The
author's own read is *ownership* (the author must know the proof he signs) and
certifies nothing, because it is not independent. Only the checker side moves a
proof toward `PROVED`:

1. **Cold reads** by fresh readers, each of whom wrote or checked no earlier
   part (repairer != checker). Unanimous "no gap" gives `REVIEW-CLEAN`, not
   `PROVED` -- same-provenance readers may share a blind spot.
2. **Certifying read** -- one reader of genuinely different provenance. That read
   is what flips `REVIEW-CLEAN` to `PROVED`.

This is why `cite()` renders every non-`PROVED` status: a result that is only
`REVIEW-CLEAN`, `GAP-KNOWN`, `SUPERSEDED`, or `ARCHIVED` travels flagged, so no
output silently anchors on a proof that is not yet established.

### Section 4 — Structural reduction

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `prop:structural-reduction` | 4.1 | σ_c = ρ_⋆ · τ | central interpretive output |
| `prop:stability` | 4.4 | Stability under C¹ perturbations | item F (γ_O indicator) |
| `rem:interpolation-rho` | (remark, see paper §4) | Interpolation as ρ_⋆ phenomenon | enforcement 8 (smoothing provenance) |
| `obs:probe-rho-star` | (observation, §4) | Probe-dependence is ρ_⋆ under faithfulness | first-stage faithfulness check copy |

### Section 5 — Spectral identification

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `def:transfer-op` | 5.2 | Positive-noise transfer operator | framework input contract |
| `prop:standard-frameworks` | 5.4 | Six standard frameworks satisfying 5.2 | item E (framework taxonomy) |
| `def:faithful` | 5.8 | Single-mode-faithful observable | faithfulness type |
| `prop:faith-sufficient` | 5.10 | Concrete sufficient conditions (F1/F2/F3) | item B |
| `thm:spectral-id-A` | 5.12 | Spectral identification, abstract perturbation | spectral solver |
| `thm:spectral-id-B` | 5.13 | Spectral identification, operator mechanism | spectral solver |
| `thm:spectral-id` | 5.14 | Spectral identification of τ (real positive case) | bridge prerequisite |

### Section 6 — Bridge

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `thm:cross-obs-concentration` | 6.1 | Cross-observable concentration | bridge output |
| `cor:two-probe` | 6.3 | Two-probe verification protocol | item C (entry point) |
| `def:operational-test-noncirc` | 6.4 | Non-circular operational single-mode test | item C (default) |
| `def:operational-test` | 6.6 | Fit-based operational test (legacy variant) | item C (`fit_based=True`) |

### Section 7 — Probe interpretation

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `prop:probe` | 7.1 | Spectral-gap probe in regime (I) | probe interpretation contract |

### Section 8 — Trichotomy

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `def:Onice-all` | 8.1 | Extended admissible class | trichotomy domain type |
| `thm:trichotomy-geometric` | 8.3 | Geometric trichotomy, unconditional | item A (geometric layer) |
| `def:thresholds` | 8.5 | (ε_c, δ_sep) thresholds | item H |
| `thm:trichotomy-spectral` | 8.6 | Spectral attribution under (ε_c, δ_sep) | item A (spectral layer) |
| `def:noise-floor-diagnostic` | 8.8 | Signal/noise diagnostic for regime III | item A (operational layer) + item H (η_O) |
| `thm:trichotomy` | 8.11 | Unified trichotomy reference | item A (umbrella result) |
| `thm:multimode` | 8.16 | Multi-mode diagnostic | regime II solver |

### Section 9 — Spectrally-flat diagnostic

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `thm:diagnostic` | 9.1 | Spectrally-flat diagnostic | regime III output |

### Section 11 — Time domain

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `prop:drift` | 11.1 | σ_c drift in regime (I) near gap-closure | time-domain extension |

### Section 12 — Transition zone

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `prop:transition-zone` | 12.1 | The transition zone | item F (γ_O low-confidence warning) |

### Appendix A — Conditional canonicity proof

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `lem:smooth-argmax-app` | A.1 | Smooth dependence of argmax | (proof infrastructure, not API) |
| `lem:richness-app` | A.2 | Richness of 𝒪_⋆ on local data | (proof infrastructure, not API) |
| `lem:k-zero-conditional` | A.5 | Boundary-killing of scaling degree | (proof infrastructure, not API) |
| `lem:aczel-A1A4` | A.6 | Aczél reduction (A1)+(A4), conditional on (R†) | (out of scope, item I) |
| `lem:aczel-A2` | A.7 | Aczél reduction (A2) | (out of scope, item I) |
| `lem:aczel-A5` | A.10 | Aczél reduction (A5)+(A4) uniqueness | (out of scope, item I) |
| `def:Onice-multi` | A.12 | Multi-mode admissible class | regime II domain type |
| `def:sigmac-multi` | A.13 | Multi-mode susceptibility-peak functional | regime II output type |
| `thm:multi-char` | A.15 | Multi-mode characterisation (sketched) | item A (II layer; honestly sketched) |

### Appendix B — Spectral identification proof

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `def:Phi-O` | B.1 | Canonical resolution-dependent observable | window construction |
| `lem:O-admissible` | B.3 | Admissibility | windowed-spectral solver |
| `lem:modal-sum-general` | B.4 | Modal-sum observables without window | bare-correlator path |
| `def:faithful-formal` | B.6 | Single-mode faithfulness, formal | item B (formal spec) |
| `lem:residual-perturbation` | B.7 | Log-perturbation of σ_c under residuals | bridge error bound |

### Appendix C — Worked examples

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `lem:multimode-perturb-stmt` | C.1 | Multi-mode argmax separation | item A (II local-max bound) |

### Appendix D — Genericity

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `thm:genericity` | D.1 | Genericity of λ_2-coupling | item B (F1 backbone) |
| `cor:generic-two-probe` | D.3 | Generic two-probe faithfulness | item C (genericity flag) |

### Appendix E — Complex λ_2

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `thm:spectral-id-complex` | E.1 | Spectral identification, complex λ_2 | (companion / experimental) |
| `prop:reversibility` | E.2 | Real-positive sufficient condition | item E (framework selection) |

## Golden tests anchored to printed numbers

Per enforcement item 1: the kernel reproduces the paper's printed numerical values as
contract tests. The anchors:

| Paper location | Numerical claim | Tolerance |
|---|---|---|
| Example C.6 (4-state Markov) | spec(P_ε) = {1, 0.6002, 0.3128, 0.1669} | rel 1e-3 |
| Example C.6 | r ≈ 0.521 | abs 1e-3 |
| Example C.6 | τ ≈ 1.959 (T_* = 1) | abs 1e-3 |
| Example C.6 | ρ_⋆ = 2 − √3 ≈ 0.2679 | exact analytic |
| Example C.6 | σ_c[O_A] ≈ 0.5576, σ_c[O_B] ≈ 0.5243 | abs 1e-3 |
| Example C.6 | cross-probe τ-agreement 6.3% | abs 0.5% |
| Example C.7 (1D Ising, βJ=0.5) | τ = −1/log tanh(βJ) ≈ 1.2954 | abs 1e-3 |
| Example C.7 | σ_c[C] = τ (ρ_⋆ = 1) | exact analytic |
| Example C.7 | σ_c[O_W] = (2 − √3)τ ≈ 0.3471 | abs 1e-3 |
| Example C.7 | grid error ≤ 0.08% (paper bound) | bound check |
| §4.2 window table | Gamma-2 ρ_⋆ = 2 − √3 | exact analytic |
| §4.2 window table | Gamma-3 ρ_⋆ = 3 − 2√2 | exact analytic |
| §4.2 window table | exp / log-Gaussian ρ_⋆ = 1 | exact analytic |
| §10 anchor table | E1 ferro σ_c = 8.0, κ = 1.79 | abs 0.01 (from dataset) |
| §10 anchor table | E2 ferro σ_c = 0.36, anti σ_c = 0.91 | abs 0.01 (from dataset) |

These are the contract tests. They are reproduced by
`supplementary/verify_examples.py` shipped with the paper
(41/41 PASS at submission, archived at the Zenodo DOI above); the kernel's
test suite inherits this checking and adds the window-family and
anchor-table cases.

## Revision protocol

1. When the paper is renumbered (camera-ready, arXiv revision, post-acceptance edits): recompile `paper.tex` to get `paper.aux`, extract `\newlabel{...}` lines, update the `Number` columns above in place. Do not change the `Label` columns.
2. Append a new entry to the "Active paper version" block with the new version date and arXiv/journal ID. Keep one history line per revision.
3. The code is unchanged — it cites labels.
4. Golden tests are unchanged unless the paper actually changes a printed value (extremely rare for numerical anchors that have a verification script).

## Out-of-scope clarifications carried from the prescription

- No label of `prop:axiomatic-char` (3.8), `lem:aczel-*` (A.6, A.7, A.10) is referenced by the API code. They live in this map for completeness and to make the out-of-scope-list inspectable.
- The complex-λ_2 spectral identification (`thm:spectral-id-complex`, E.1) is a paper-companion result; the kernel marks it as `experimental`, behind a separate feature flag.
- The `lem:multimode-perturb-stmt` "exactly k vs at least k" subtlety (Rem 8.15 of the paper, the Notion-3 bug 1 fix) is implemented as the conservative "at least k"; the tighter version requires a cross-term suppression check which the kernel exposes as an opt-in audit, not a default verdict.

## Deliberate verdict changes (recorded before the code change, per codes.py)

These are intentional changes to what `classify()` / `sigma_c_verdict()` return, made because the prior behaviour was a false-green on a degenerate input — the sibling class of the NaN guard. The golden anchors (`test_codes_invariance.py`, `test_golden_paper_anchors.py`) are UNCHANGED by these: they cover no constant / >5-peak input, so the formalised verdict still reproduces the anchors exactly. The change touches only inputs the anchors do not exercise. Pinned by `test_input_degeneracy_guards.py`.

- **6.0 — constant / near-constant observable → `NOT_IDENTIFIED` / `NOT_RESOLVABLE`.** A constant observable has `chi_O = |dO/dlog σ| = 0` up to rounding; the peak search, normalising χ to its own max, previously found phantom peaks in that rounding noise and returned `OK`. Two DECOUPLED guards now handle this (decoupling + the constant choices below came out of an adversarial different-provenance review that found a margin false-green and a coupled false-refusal in the first attempt):
  - **Guard 1 (`Trichotomy.degenerate_O`):** `range(O) <= DEGENERATE_ULP_FACTOR·eps·max|O|` (factor 16, **no** `·n` — a max−min does not accumulate over the sample count, and an `·n` term wrongly refused real small signal at large `n`), or `O == 0`. Forces regime III → `NOT_IDENTIFIED` (remediable: supply an observable that varies with the dial). Catches only observables constant to a few ULP.
  - **Guard 2 (per-peak floor):** a χ peak must clear `NOISE_FLOOR_FACTOR·eps·max|O|/min(dlog)` (factor 128, using the **smallest** log-step so a boundary-differencing artifact at the fine end of a linear grid cannot pass). Below it there is no peak → regime III → `NOT_RESOLVABLE`. This is the workhorse: the honest cutoff sits at ≈1e-12 relative amplitude (χ/floor ≈ 20), so a real bump 1000×+ above machine ε resolves to `OK`, while both pure rounding noise and edge artifacts (χ/floor ≈ 2) are refused — consistently, whether they come from a tiny linear ramp or a sub-floor bump.

- **6.0 — convention-window stability gate (the last known path to a false `OK`, closed).** `OK` (regime I or II) now requires the resolved peak COUNT to be unchanged across the declared window `r ∈ [r0/CONVENTION_WINDOW_FACTOR, min(1, CONVENTION_WINDOW_FACTOR·r0)]` (factor 2 — an octave of `min_prominence_ratio` each way; a scale-free declared convention). A peak set that exists only inside a narrow prominence window (e.g. white noise forced to a single peak by a high `r`, or the pre-fix hostile window `r ∈ [0.38, 0.40]`) is a Layer-3 reading that runs away under the convention screw — not a measurement — so it becomes `NOT_RESOLVABLE`. Applies to scalar AND vector `sigma_c` (via `Result.convention_window_stable()`), uses the same `chi_abs_floor` as the verdict, and reuses the sweep already computed for `resolution_band`. Empirically: pure white noise returns `OK` at **no** `r` across 15 seeds × 40 `r`-values; clean single- and multi-mode signals stay `OK`. Pinned by `test_white_noise_never_ok_across_convention_grid` / `test_convention_stable_signal_is_ok`.
- **6.0 — resolved peak count beyond the declared ceiling → `NOT_RESOLVABLE`.** More than `max_resolved_peaks` (default 5, a convention) interior peaks read off one dial is noise, not structure (e.g. a raw decaying time series fed as (σ,O) produced ~100 "peaks" and an `OK`). `Trichotomy.peak_count_over_ceiling` makes the verdict `NOT_RESOLVABLE`; the vector `sigma_c` is still reported. Tunable per call via `analyze(max_resolved_peaks=...)`.

## Known limitations register (carried, not buried; also in the README)

The honest bounds of the shipped kernel, stated where they can be audited:

- **A single time series yields no τ.** σ_c needs a (dial, observable) pairing; one decay curve is not one. A certified τ requires an operator, a verified self-adjoint spectrum with an isolated real gap, or a two-probe agreement (`two_probe_test`). Feeding a raw series as `(σ, O)` correctly returns `NOT_RESOLVABLE`.
- **The jitter / second channel is experimental and uncertified.** `sigma_c.experimental.jitter` emits no 4-code verdict (`validated=False` everywhere), is not imported by the kernel, and must not back a certified claim until its pre-registered validation (`PREREG_jitter.md`) is run and read by a checker of different provenance.
- **`cor:twoprobe` (the finite-window τ gate) is `GAP-KNOWN`.** The corrected tail-window condition `T > max{t*_A, t*_B, 1/|mu_3|}` is applied in `window_readability_verdict`, but the underlying corollary as written carries a located gap; the code renders the gate, never a `PROVED` label on it. See `thm:spectral-id-A/B` (both `GAP-KNOWN`).
- **The Lean proofs are machine-checked, not independently reviewed.** `proofs/lean` (App B.3/B.4 *algebra* only) compiles clean under Lean+Mathlib with no `sorry` and standard axioms, but has had no independent human read; it ships in **neither** wheel nor sdist (git repo only) and stays out of any tag until reviewed.
- **`tau_from_profile` is CONDITIONAL, not certified.** It returns `τ = σ_c/ρ⋆` only when the χ_O profile matches a user-**declared** single-mode template (exponential / Debye) within a declared residual (Prop 2 makes the single-mode bridge exact). The result carries `certified: False`, the residual, and the declared profile; a poor fit is refused (`NOT_APPLICABLE`). Validated in `test_profile_tau.py`: single modes accepted with the right τ, wrong-profile and well-separated (ratio ≥ 3) two-mode curves rejected; barely-separated (ratio 2) modes return an effective τ with a visible residual — the residual is the honesty, not a hidden failure.
- **`bootstrap_sigma_c` is an APPROXIMATE statistical CI.** A nonparametric replicate bootstrap for σ_c (sampling variance — a *different* quantity from the resolution band). The percentile bootstrap of a peak location under-covers; `test_bootstrap.py` pins the measured coverage floor (≈90% at nominal 95%). The result is flagged `is_confidence_interval: True` **and** `approximate: True` with the measured-coverage caveat, so the nominal level is never passed off as exact.

---

# Second anchor set — "The Parrot's Theorems" (Zenodo preprint)

> **Foundation preprint vs the submitted Parrot.** The mapping above binds the
> sigma_c labels to an earlier numbering scheme (used only as a
> stable label-anchor). The results live in *The Parrot's
> Theorems* (`10.5281/zenodo.22066713`). Only the two labels below demonstrably
> survive into the submitted Parrot; the other foundation labels became prose or
> were dropped in the condensation, so they keep resolving to the foundation
> numbers (which is honest -- they ARE foundation-paper theorems) until each is
> re-pointed to a published home.

- **Paper title:** The Parrot's Theorems
- **Maintainer:** ForgottenForge
- **Preprint:** Zenodo `10.5281/zenodo.22066713`. Frozen source:
  `parrots_theorems.tex` — read-only.
- **DOI note:** `22066713` is the **version** DOI (v1), not the concept DOI, so it
  cannot silently change what the package claims.

## Mapping table — Parrot's Theorems

Section numbers are by order of appearance in `parrots_theorems.tex`
(Introduction=§1). Formal-environment numbers are marked pending — they are not
determinable from the source; they need the compiled PDF.

| Label | Number | Title (short) | maps to |
|---|---|---|---|
| `sec:setup` | 2 | Dial, reading, and the group that owns them (the Layer table) | blindness map: layer ladder |
| `sec:corank` | 4 | The corank: what never shows (the null cone) | blindness map: null-cone note (honesty anchor) |
| `sec:jitter` | 5 | Counting what rustles: the jitter theorem | published §5 + App B.3/B.4/C (this DOI); code in `experimental/jitter` (parts i-iv) — dial-free, NOT certified (App C validation pending, see experimental/PREREG_jitter.md) |

**Formal labels that survived the condensation** (the only two): `thm:noninv`
(Non-invariance of peak locations, in §7) and `prop:profile` (Profile
decomposition, §7). Numbers pending the compiled PDF; not added with guessed
numbers, because a wrong "Thm 7.1" is a false locator.

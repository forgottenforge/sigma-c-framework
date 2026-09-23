# Pre-registration — the jitter channel's two-channel run

The pre-registration for this channel is **Appendix C of The Parrot's Theorems**
("The two-channel run (pre-registered)", gauge fixed 2026-08-21, before any result;
Zenodo doi:10.5281/zenodo.22066713). This file restates its fixed criteria so the
validation can be held to them; the paper is the authority.

> **Status:** the code was rebuilt to match Sections 5 + App B.3/B.4/C (the first
> module implemented a time-series `rustle_rate` that is **not in the paper** and
> failed validation — see VALIDATION_RESULTS.md). The rebuilt module HAS since been
> run through this validation by independent different-provenance checkers across
> four rounds — all App C criteria pass (VALIDATION_RESULTS.md); promotion out of
> `experimental/` still awaits the professor's certifying read.

## The channel (Section 5, App B.3/B.4)
Clamp all dials; take repeated block readings of several instruments. The
block-to-block covariance is `C = Var(delta) a a^T` (App B.3). Four parts:
(i) detection (off-diagonal beyond noise), (ii) counting (rank above the noise
floor), (iii) identification (leading eigenvector ∝ the dial susceptibilities
`a_i = d<O_i>/dgamma`), (iv) cross-check (jitter slope `Cov(M1,M2)/Var(M2)` = dial
slope `dM1/dM2`, under detailed balance). Tetrad (App B.4) is the k=3 one-factor
blind test `C12 C13 / C23 = a1^2 sigma_y^2`.

## The run (App C) — fixed criteria
- **Target:** a system with a meaningful dial channel (the AVS decoherence sweep,
  gamma_c ≈ 0.674, kappa = 8.58), not a discovery setup.
- **Acquisition:** a gamma sweep of ~15 points dense around gamma_c; per
  (gamma, basis) at least 10 time-ordered blocks of at least 50 shots; at least two
  bases, one with >= 3 qubit pairs (so Channel 1 = per-shot covariance of pair
  parities runs alongside Channel 2 = drift covariance of block means).
- **Null:** shuffle null, `n_shuffle = 200`, significance at the 95th percentile.
- **Conventions:** linear derivative for the bounded dial; central differences.
- **Pass criteria (declared before any result):**
  1. the dial channel `chi_W(gamma)` and the jitter channel `lambda_1(gamma)` peak at
     the same grid point, **within one grid spacing**;
  2. **`cos(v1, a) >= 0.9`** (the leading eigenvector aligns with the dial
     susceptibilities);
  3. with the drift switched off (dials truly clamped), `lambda_1` returns to the
     shot floor / below the shuffle null (the negative control raises no false positive).

## Reference + negative control
- **Reference generator:** `reproduce_numbers.py` (Zenodo record, fixed seed) —
  Figure 2 signature (both channels peak at the same grid point; drift off ->
  jitter floor) and the five-row logbook cross-check (jitter slope 0.144, analytic
  `2 sech^2(1+2 ybar) = 0.140`, dial slope 0.141). The validation reimplements this
  generator from the App C equations (fixed seed); it does not import across repos.
- **Hardware negative control (documented):** 27 time-ordered measurement-only
  calibration blocks (500 shots each) returned `lambda_1` below both the shot floor
  and the shuffle null's 95th percentile.

## Promotion criteria (all required, unchanged)
1. All three App C criteria met on the reference generator + the logbook cross-check.
2. A different-provenance reader confirms the analysis code was fixed before the run
   and the criteria were not adjusted to fit the outcome.
3. The repair state is frozen and the outside-runs have passed first.

Until then the channel stays in `experimental/`: `validated=False`, no `OK` verdict.

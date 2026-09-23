# Validation results — experimental jitter channel

## Round 1 (first module) — a wrong build, correctly caught

The first module implemented a `rustle_rate` that fitted a decay rate from the
**autocorrelation of a single time series**. An independent, different-provenance
checker ran a pre-registered validation and **V1 (`rustle_rate`) FAILED**: its
reported band was only the log-slope's regression standard error (precision, not
accuracy), undercovering to ~0% even on a noiseless single mode and depending
strongly on the fit window; the pure-noise guard was defeated by a forced window.
V2 (`covariance_structure` shuffle-null) and V3 (a slope comparison) passed.

**The deeper finding (from reading the paper): `rustle_rate` is not in the paper at
all.** The jitter theorem (Section 5, App B.3/B.4/C) operates on **clamped-dial
repeated readings across several instruments**, not a time-series autocorrelation.
Fitting a slow rate from `C(t)` is an inverse-Laplace problem the theory never
poses — which is *why* V1 could not be made honest. Onsager is cited as the lineage
of part (iv), not as a protocol.

## Rebuild — to the paper

`rustle_rate` and the tau-vs-tau comparison were **removed**. The module now
implements the actual theorem, each part citing its source:

- (i) detection + (ii) counting — `jitter_rank` (drift covariance + the shuffle-null
  that PASSED round 1; ~5% false-positive through d=30, retracted-flatness avoided);
- (B.4) the k=3 one-factor `tetrad`;
- (iii) identification — `leading_coupling` + `cos_coupling` (cos(v1, a));
- (iv) cross-check — `jitter_slope` + `two_channel_fdt` (jitter slope vs dial slope,
  the real FDT test, reported with bands, never a verdict);
- App C two-channel signature — `peak_agreement` (dial chi vs jitter lambda_1 peak
  within one grid cell).

Paper-anchored known-answer in the suite: the five logbook rows of Section 5
reproduce jitter slope **0.144**, analytic **0.140**, dial **0.141**.

## Round 2 (App C on the rebuilt module) — one defect, now fixed

An independent, different-provenance checker ran the App C validation on the rebuilt
module (>= 20 seeds; reference generator reimplemented from the paper's equations):

- **Criterion 1 (peak agreement within one grid cell): PASS** (20/20; the checker
  honestly noted the reference grid is coarse relative to the feature, so this is a
  consistency check, not an independent corroboration of the peak).
- **Criterion 2 (cos(v1,a) >= 0.9): PASS** strongly (min 0.988, mean 0.997).
- **Section-5 logbook cross-check: PASS** (jitter slope 0.1435, analytic 0.1402).
- **Tetrad (App B.4): PASS** (one-factor ~1; two-factor broken; a genuine 6-instrument
  two-factor build correctly flagged).
- **Criterion 3 (false-positive control): FAILED, now FIXED.** The leading eigenvalue
  was calibrated (~5%), but the returned flag `detected = rank>=1` summed a per-order
  test across all d eigenvalue orders, so its false-positive rate inflated with
  instrument count: **9.5% at d=3 -> 80% at d=30**. The docstring's "~5% through d=30"
  was refuted for the returned flag.

**Fix applied:** `jitter_rank` now uses a **step-down** test — detection is the single
leading-eigenvalue test at the declared level (calibrated ~alpha at every d), and the
rank counts consecutive orders from the top, stopping at the first that fails its
null, so the pure-noise rank cannot inflate. Re-measured false-positive rate (n_shuffle
=200, 150 seeds): d=3 **6.0%**, d=10 **5.3%**, d=20 **2.0%**, d=30 **6.0%**; rank>=2
false positive <= 1.3%; real drift still detected (rank 1). Docstrings corrected.

**Gap closed:** App C's Channel 1 (per-shot covariance of co-measured pair parities)
is now implemented — `pair_parity_rank` / `pair_parity_covariance`, the shot-level
twin of Channel 2's `jitter_rank` (same step-down detection/counting machinery on the
per-shot matrix). The full App C two-channel acquisition is now in the module. Both
channels should be exercised in the next validation run.

## Round 3 (independent re-check of the fix) — the fix holds

A second, independent different-provenance checker re-ran App C on the fixed module:

- **Criterion 3 (false-positive control): PASS.** detected-FP on iid nulls (n_shuffle
  =200): d=2 5.5%, d=3 5.0%, d=5 4.0%, d=10 5.2% (400 seeds), d=20 5.5%, d=30 5.5%
  (400 seeds) — at the declared ~5% at **every** instrument count, no climb with d
  (was 9.5%→80%); `rank>=2` FP <= 1.5%.
- **Detects real structure: PASS** — drift-ON detected 27/30 (rank 1); a genuine
  comparable two-factor build → rank 2 in 30/30, never over-counting.
- **Criterion 1 (peak agreement): PASS** (30/30). **Criterion 2 (cos(v1,a)>=0.9):
  PASS** (49/50). **Logbook: PASS** (0.1435 / 0.1402 / 0.1413). **Tetrad: PASS**.
  **Verdict hygiene: PASS** (validated=False everywhere, no OK string).

**Honest limitation (not a criterion failure):** counting (ii) is now conservative —
a dominant first factor can absorb a much weaker second one (reported as rank 1);
rank 2+ is reliable only when factors are comparable. This is the price of the FP
control (under-counts rather than over-counts) and is documented in `jitter_rank`.

The round-3 run was launched **before** Channel 1 was added (commit adding
`pair_parity_rank`), so it validated Channel 2 + the fix. A final run should exercise
**both** channels on the complete module before promotion.

## Round 4 (final, complete module — both channels) — PASS across the board

An independent different-provenance checker validated the WHOLE module against App C,
per criterion (Channel 1 on its OWN binomial per-shot data path, not by inheritance).
One dataset (App C three-basis-circuit acquisition) drove both channels; >=20 seeds,
n_shuffle=200.

- **A. Channel 1 (per-shot pair parities) — its own data path:**
  - A1 false-positive control (true null): detected-FP centred ~5% across pairs
    d=2..20 and shots N=400..2000, **no climb with d** (grand mean ~5.3%); rank>=2
    FP <= 2%. The per-shot permutation null on binomial data is calibrated. **PASS.**
  - A2 counting: one shared per-shot cause -> rank 1 (100%); two comparable causes ->
    **rank 2 (100%)**, no over-count. It COUNTS, not merely runs. **PASS.**
  - A3 the binomial shot floor makes no spurious off-diagonal (obs ~= null). **PASS.**
- **B. App C criteria (Channel 2 in the two-channel layout):** B1 shuffle-null n=200
  @95th used by detection; B2 peak agreement 25/25 within one cell; B3 cos(v1,a) min
  0.985, 100% >= 0.9; B4 slope identity consistent (see caveat 2). **all PASS.**
- **C. Logbook (Sec.5):** 0.1435 / 0.1402 / 0.1413 (paper ~0.144/0.140/0.141). **PASS.**
- **D. Negative control, drift OFF (Fig 2):** Channel 2 FP 1.7%, lambda_1 back to the
  floor; Channel 1 on clamped/independent shots does not fire. **PASS.**
- **E. Tetrad touchstone + the counting limitation, EMPIRICALLY PINNED:** one factor
  -> tetrad ratio 1.000; two comparable -> 0.487 (far from 1). The step-down
  `jitter_rank` loses a weaker second factor while the tetrad still catches it:
  - *disjoint blocks:* rank-2 recovered 100% for r >= 0.64, 60% at r=0.36, **0% at
    r <= 0.25** — crossover **r\* ~ 0.30** (r = second/first loading power);
  - *common dominant factor + weaker second:* `jitter_rank` reports rank 1 at **all
    r (even r=1)** — the common mode inflates the order-2 null and hides the second
    factor at every strength — while the tetrad ratio moves 0.487 (r=1) -> 0.796
    (r=0.25) -> 0.917 (r=0.09), statistically != 1 far below where the step-down fires.
  So a genuine weak second factor is missed by counting but caught by the tetrad; the
  module steers users to the tetrad for exactly this. A finding, not a bug. **PASS.**
- **F. Verdict hygiene:** validated=False on every result; no standalone OK verdict
  anywhere. **PASS.**

**Three honest caveats (none a criterion failure):**
1. Peak agreement (B2) is a coarse-grid consistency check (grid spacing 0.053 vs the
   feature at gamma_c=0.674), not an independent localization of the peak.
2. The slope (B4) carries a shot-noise errors-in-variables bias (~6% at 3k shots)
   that vanishes with more shots (0.45% at 30k); consistent within the band, but a
   single acquisition reads a few percent low — carry a band, don't quote the point.
3. Counting (ii) under-counts by design (pinned in E); the tetrad is the instrument
   for a weak second factor.

## Status

**The complete module meets every App C criterion from outside** (rounds 2->3 fixed
and confirmed the false-positive control; round 4 validated both channels per
criterion). It stays `experimental/` — `validated=False`, no `OK`, kernel does not
import it — because promotion out of `experimental/` requires the professor's read
(PREREG criterion 2) and the freeze + outside-runs first (criterion 3). The freeze
itself is not blocked: the whole module has now been seen green from outside once. It delivers **rank** and **peak location**
independent of the dial, and a detailed-balance cross-check; it delivers **no
certified tau** (the peak -> tau step is the fallen sigma_c = rho_star*tau bridge).

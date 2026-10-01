# Changelog

All notable changes to `sigma-c-framework` are documented here. This project
follows [Semantic Versioning](https://semver.org/).

## Migrating from 5.x to 6.0

**`import sigma_c` now gives the disciplined-reader kernel.** In the 5.x line on
PyPI, `import sigma_c` exposed an adapter/utility stack; in 6.0.0 it is the
scale-selector kernel — `analyze()`, the 4-code applicability contract
(`OK` / `NOT_APPLICABLE` / `NOT_RESOLVABLE` / `NOT_IDENTIFIED`), the blindness
map, and the certified/fallen τ split. This is a **breaking change**: code
written against a 5.x `sigma_c` API will not run unchanged against 6.0.

- The PyPI **package name is unchanged** (`sigma-c-framework`); `pip install -U
  sigma-c-framework` moves you from 5.x to 6.0.0 and changes what `import
  sigma_c` means. Pin `sigma-c-framework<6` if you depend on the old surface.
- 6.0 is a **fresh public release** of the kernel; there is no automatic shim
  from the 5.x API. Port callers to `analyze(sigma, O)` and read `result.summary()`.
- Optional quantum-SDK adapters live under `sigma_c.adapters.*` and install via
  extras (`pip install "sigma-c-framework[qiskit]"` etc.); the kernel never
  imports them.

## 6.0.1

### Fixed — one shared dial contract across all entry paths

`analyze()` accepts a dial three ways — an observable array `O`, a callable
`O`, or a precomputed susceptibility `chi`. The array path was validated by
`chi_O`; the callable and `chi` paths routed around it and re-checked only for
`NaN`/`Inf`. As a result a malformed dial was silently processed on those two
paths where the array path cleanly rejected it. The full guard now runs once,
before the paths diverge, so all three reject identically:

- **`sigma <= 0`** — rejected on every path (previously the `chi` path reached
  `log(sigma)` and emitted a `divide by zero` warning; a `NaN` result was
  reachable when a peak sat on a non-positive sample, and otherwise the reading
  came back `OK` with the absolute noise floor silently set to zero).
- **fewer than 5 dial samples** — rejected on every path (previously accepted
  on the `chi` and callable paths).
- **duplicate / non-increasing `sigma`** — rejected on every path (zero
  log-spacing makes the `chi` derivative undefined; previously processed with
  numpy warnings or silently).
- **`chi` (or an `O` array supplied with `chi`) whose length differs from
  `sigma`** — rejected with a clear message (previously an `IndexError` when
  shorter, or a silent truncation/misalignment when longer).

No API or behavioural change for well-formed input; results are unchanged.

### Fixed/hardened — seven items from testing 6.0.0 on measured kinetics

- **Topographic prominence (was a height threshold).** `find_interior_maxima`
  now accepts an interior maximum only if it rises `min_prominence_ratio * max(chi)`
  above the higher of its two flanking saddles — the prominence the parameter name
  always promised. A bare height gate counted a noise wiggle on the flank of the
  dominant peak as a separate peak and fabricated regime II at moderate SNR. The
  convention-stability gate was co-updated so that the degenerate drop to **zero**
  peaks at an extreme convention (no peak on a baseline reaches prominence = max χ)
  is not mistaken for an instability; only a different **non-zero** count is.
- **Regime-I shape note.** A merged two-mode relaxation also yields one χ peak, so
  `OK` now carries a non-fatal note when the single peak's shape residual against a
  one-mode template is large: `OK` means "one resolved peak under the convention",
  not "one physical mode".
- **Profile route propagates the kernel verdict.** `tau_from_profile` no longer
  returns `OK_CONDITIONAL` over a scalar `sigma_c` whose kernel status is
  `NOT_RESOLVABLE`/`NOT_APPLICABLE`.
- **Window weight warning.** `analyze()` warns when a non-bare window is passed on a
  raw-observable path (its weight is not applied to `O`; the window only sets
  `rho_star`).
- **Declarable reversibility tolerance.** `analyze(..., selfadjoint_rtol=...)`
  forwards a declared tolerance to the self-adjointness/reversibility check, so an
  operator built from rate constants rounded to 3–4 digits (symmetry residual ~1e-3,
  far above the ~1e-13 default) can be certified with a reported tolerance instead
  of a brittle refusal. Default unchanged.
- **Summary flags a rate-error bound > 100 %** as DOMINANCE-only.
- **Docs:** `tau_abscissa` is the slowest mode of the whole operator; a faster
  sub-process is never the certified abscissa (isolate it on a restricted operator).

## 6.0.0

### Added
- `analyze()` disciplined-reader entry point; every numeric output carries its
  own 4-code applicability verdict, a blindness map, and provenance.
- `to_dict()` now stamps `schema_version` and `library_version` so a stored
  result is reproducible later.
- `tau_from_profile(sigma, O, profile=…)` — the conditional yes-path for a decay
  curve: declare a single-mode profile (exponential / Debye), the tool checks the
  fit and returns `τ = σ_c/ρ⋆` conditional on the declaration (residual attached),
  or refuses a poor fit. Not certified; carries `certified: False`.
- `bootstrap_sigma_c(sigma, replicates)` — an approximate statistical CI for σ_c
  from replicate curves (sampling variance), distinct from `resolution_band`.
  Honestly flagged `approximate` with its measured-coverage caveat.
- Typed exceptions `SigmaCError` / `InvalidInputError` (the latter subclasses
  `ValueError`, so existing handlers keep working).
- Plateau-aware peak finder: a smooth peak whose maximum falls between two grid
  points (tied top samples) is now resolved instead of missed as regime III.
- Opt-in SDK adapters `sigma_c.adapters.{qiskit,pennylane,braket}` (extras).
- `sigma_c.experimental.jitter` — the dial-free second channel, explicitly
  uncertified (`validated=False`), not imported by the kernel.
- Machine-checked Lean proofs of the jitter algebra (App B.3/B.4) under
  `proofs/lean/` (sdist only).

### Changed — input-degeneracy guards (deliberate verdict change)
- A **constant / near-constant observable** (no variation above the
  numerical-noise floor) now returns `NOT_IDENTIFIED` instead of a spurious `OK`
  read off rounding noise. Remediation: supply an observable that varies with
  the dial.
- A **resolved peak count beyond `max_resolved_peaks`** (default 5, a declared
  convention) now returns `NOT_RESOLVABLE`; the vector `sigma_c` is still listed.
  Tunable via `analyze(max_resolved_peaks=...)`.
- A **peak set that holds only under a narrow prominence convention** now returns
  `NOT_RESOLVABLE`: `OK` requires the resolved peak count to be stable across an
  octave of `min_prominence_ratio` (r ∈ [r0/2, 2·r0]). This closes the last known
  path to a false `OK` — pure white noise returns `OK` at no `min_prominence_ratio`.
- All three guards are recorded in `THEOREM_MAP.md` (§ Deliberate verdict changes)
  and pinned by `test_input_degeneracy_guards.py`. The golden anchors are unchanged.

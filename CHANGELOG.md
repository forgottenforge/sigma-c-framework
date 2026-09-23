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

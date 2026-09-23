<p align="center">
  <img src="sigma_c/c.jpg" width="340" alt="σ_c">
</p>
<p align="center">
  <img src="sigma_c/forgottenforge.jpg" width="220" alt="ForgottenForge">
</p>

# SIGMA C v6 — the disciplined-reader kernel

**Find the scale at which a measurement is most sensitive to a change of scale —
and get told, honestly, how much to trust it, or that it cannot be trusted.**

[![License: AGPL-3.0-or-later OR Commercial](https://img.shields.io/badge/License-AGPL_v3_%7C_Commercial-blue.svg)](#license)
[![PyPI](https://img.shields.io/badge/pip-sigma--c--framework-green.svg)](https://pypi.org/project/sigma-c-framework/)
[![Paper](https://img.shields.io/badge/paper-Zenodo_22066713-blue.svg)](https://doi.org/10.5281/zenodo.22066713)

Package name `sigma-c-framework`, version **6.0** · import path `sigma_c`. The wheel
ships **one** thing: the disciplined-reader kernel, plus opt-in adapters. Code without
tests or a paper anchor is kept out of the distribution on purpose — *only what works
ships.* (Import `sigma_c` changed meaning from the 5.x line — see
[migrating 5.x → 6.0](CHANGELOG.md).)

## What is σ_c?

You have a curve: some measurement `O` plotted against a **dial** you were turning —
a length, a time, a temperature, a probability, whatever you varied. σ_c is the dial
value where the curve reacts **most sharply** to you changing the dial: the scale
where "something happens". Formally it is the peak of the sensitivity `|σ · dO/dσ|`:

```
sigma_c[O] := argmax_sigma | sigma · dO/dsigma |   # argmax = the dial value that makes this biggest
```

(`dO/dσ` is just *how fast the measurement `O` changes as you nudge the dial*, so
`|σ · dO/dσ|` is "how sharply `O` reacts" and σ_c is the dial value where that is biggest.)

That location is a **reading under a convention you declared**, not an intrinsic
constant of the system — and the kernel says so, every time.

Every output carries three things:
- a **4-code verdict** — `OK` (trust it) / `NOT_RESOLVABLE` (your data can't resolve
  it) / `NOT_APPLICABLE` (wrong tool for this) / `NOT_IDENTIFIED` (*you* left out an
  input — it tells you what to bring);
- a **resolution band** — the finest step the tool's ruler can tell apart here (how
  many digits you may honestly write down), **not** a fabricated ± error bar;
- a **blindness map** — what this measurement simply cannot see.

> **New here? Jump to [the learning ladder](#start-here-the-learning-ladder).** It
> opens with "we hid the number 4 in a curve and the tool found 4.00" and needs no
> jargon. The advanced part (the relaxation time τ) has its own section further down.

## What this tool needs — and one thing it will not do

σ_c needs a **dial you varied** and an **observable measured across it** — a curve
`O(σ)`. It reads a *scale on that dial*.

**If you have a single time series and want a relaxation time τ, this is the wrong
tool for that step, and it will refuse rather than guess.** One decaying signal
`x(t)` is a measurement over time, not a (dial, observable) pairing, and there is no
honest way to read an intrinsic τ from it alone: feeding `x(t)` in as `(σ, O)`
returns a `NOT_RESOLVABLE` peak-forest (that is the correct refusal, not a bug). A
certified τ needs one of:
- a **transfer operator / generator** (`analyze(..., operator=M, T_star=..., sigma_axis=...)`), or
- a **verified spectrum** with an isolated real gap (`spectrum=[...]`), or
- **two faithful probes** of the same system, compared with `two_probe_test(r1, r2)`.

If you have none of these, σ_c cannot turn your one curve into τ from first
principles. There is **one conditional yes-path**: if you are willing to *declare*
that your observable follows a single-mode profile (an exponential relaxation or a
Debye peak), `tau_from_profile(sigma, O, profile="exponential")` fits that declared
template, reports the goodness-of-fit residual, and — only if the fit is good (Prop 2:
the bridge is exact for a single mode) — returns `τ = σ_c/ρ⋆` **conditional on your
declaration, with the residual attached**. A poor fit (a multi-mode curve) is
refused, not fudged. This is not the certified spectral abscissa; it is an honest
profile-conditional reading, and it says so.

## What σ_c is *not* for

- **Not a predictor or optimiser.** σ_c is a *diagnostic/observational* instrument:
  it reads where a curve is most scale-sensitive. It does not forecast, and it is not
  a loss to minimise. If you want prediction or optimisation (finance, "predictive
  analytics", hyperparameter tuning as an end in itself), this is the wrong tool.
- **σ_c is not "the optimum".** The peak of `|σ·dO/dσ|` is the scale where the
  observable *changes fastest*, not the scale that is *best* by any objective. For a
  loss-vs-learning-rate sweep it marks the knee (where behaviour turns over), **not**
  the optimal learning rate — reading it as "the optimum" is a category error. σ_c is
  a Layer-3 *location under a convention*, never an extremum of a payoff.
- **Not a statistical confidence interval.** The `resolution_band` is the ruler's
  finest distinguishable step (a discretisation/convention band), **not** a ±CI. If
  you have replicates and want a statistical band, `bootstrap_sigma_c(sigma,
  replicates)` gives one — a nonparametric replicate bootstrap, flagged
  `is_confidence_interval: True` and honestly `approximate` (the percentile bootstrap
  of a peak location under-covers; validated ~90% at nominal 95%). It answers the
  sampling-variance question `resolution_band` does not.

## Install

```bash
pip install sigma-c-framework            # the kernel: numpy only, fully headless
pip install "sigma-c-framework[plot]"    # + matplotlib, for result.card() PNGs
pip install "sigma-c-framework[qiskit]"  # + one quantum-SDK adapter (or [pennylane] / [braket])
```

The kernel depends on **numpy only** and imports no plotting stack — `analyze()`,
the verdict, `to_dict()` all run in a headless pipeline. matplotlib is the `[plot]`
extra (needed only for the PNG card); the SDK adapters are per-SDK extras.

## What's in the box — and where to find it

| You want to… | Use | Docs |
|---|---|---|
| read a scale from a curve `O(σ)` | `analyze(sigma, O)` → `Result` | this file · [`sigma_c/README.md`](sigma_c/README.md) |
| see the verdict / band / blindness map | `result.summary()`, `result.to_dict()` | [`sigma_c/README.md`](sigma_c/README.md) |
| a PNG card | `result.card("out.png")` (needs `[plot]`) | — |
| the certified relaxation time τ | `analyze(..., operator=…, T_star=…, sigma_axis=…)` | [Advanced — τ](#advanced--the-relaxation-time-τ) |
| τ from a **declared** single-mode profile (decay curve) | `tau_from_profile(sigma, O, profile=…)` | conditional yes-path; residual attached |
| a statistical CI from **replicates** | `bootstrap_sigma_c(sigma, replicates)` | approximate, distinct from `resolution_band` |
| feed a quantum-SDK sweep | `sigma_c.adapters.{qiskit,pennylane,braket}` (extras) | adapter docstrings |
| learn by example | the five-step ladder | [`sigma_c/examples/README.md`](sigma_c/examples/README.md) |
| bring your own data | — | [`sigma_c/BYO_GUIDE.md`](sigma_c/BYO_GUIDE.md) |
| know the honest bounds | — | [`sigma_c/LIMITS.md`](sigma_c/LIMITS.md) |
| trace a value to the paper | — | [`sigma_c/THEOREM_MAP.md`](sigma_c/THEOREM_MAP.md) |
| the experimental second (dial-free) channel | `sigma_c.experimental.jitter` — **uncertified** | its `PREREG_jitter.md` |

Not shipped in the PyPI package: the Lean proofs of the jitter algebra
(`proofs/` in the source repo only, pending an independent review).

## Quick start

```python
import numpy as np
from sigma_c import analyze, bare

sigma = np.geomspace(0.1, 100, 400)
O = np.exp(-sigma / 5.0)                 # a smooth fade on scale 5
result = analyze(sigma, O, window=bare())
print(result.summary())                  # sigma_c ~ 5.0, with its band + verdict
```

A curve with no special scale returns `None` honestly, not a made-up number:

```python
result = analyze(sigma, np.sqrt(sigma), window=bare())
assert result.sigma_c is None            # regime III; verdict NOT_RESOLVABLE
```

And if you ask for something you didn't hand it the inputs for, it says so and
**names what to bring** — a different kind of "no", `NOT_IDENTIFIED`:

```python
# ask for the certified relaxation time from a spectrum alone (no operator, no probe)
r = analyze(sigma, O, window=bare(), spectrum=[1.0, 0.6, 0.3],
            T_star=1.0, sigma_axis="evolution_time")
print(r.window_readability_status.code.value)   # NOT_IDENTIFIED
print(r.window_readability_status.remediation)   # "...bring the operator + a probe vector"
```

That is the line the four verdicts draw: **`NOT_APPLICABLE`** = this tool can't answer
that question; **`NOT_IDENTIFIED`** = it *could*, you just haven't handed it the missing
piece — and it names the piece.

## Start here: the learning ladder

New to it? Walk up [`sigma_c/examples/README.md`](sigma_c/examples/README.md) —
five steps from "hello world" to your own data, in plain language:

1. **Hello world** — hide a number in a curve, find it again.
2. **Find a textbook number** — recover the percolation threshold p_c = ½.
3. **The tool says "no"** — why an honest refusal is a good answer (incl. a real
   earthquake-catalogue example: the tool finds your instrument's limit, not physics).
4. **For the advanced** — the certified relaxation time τ (self-adjoint operators).
5. **Your own data** — [`sigma_c/BYO_GUIDE.md`](sigma_c/BYO_GUIDE.md) +
   `examples/TEMPLATE_adapter.py`.

Run one:
```bash
PYTHONPATH=. python sigma_c/examples/demo_percolation.py
```

## Advanced — the relaxation time τ

A `τ` (relaxation time — how long a system takes to settle) is available **only** when
you hand the tool a *self-adjoint* (reversible) operator; then `τ = −T*/log λ₂` is
**certified**. The old shortcut `σ_c = ρ⋆·τ` is **not** an identity: σ_c reads the
*most-susceptible* mode, which is the slowest only under a special condition — so the
kernel keeps two separate fields, `tau_bridge` (a diagnostic) and `tau_abscissa` (the
certified one), and refuses to pass one off as the other. Details in
[`sigma_c/README.md`](sigma_c/README.md) and example 6 of the ladder.

**Plain words for the jargon** (skip all of this if you're new):
- **dial / σ** — whatever you were varying (a length, a time, a probability).
- **argmax** — the input value that makes a quantity biggest.
- **relaxation time τ** — how long a system takes to settle back down.
- **ρ⋆** — a fixed number for the "lens" (window) you chose; divide by it to compare lenses.
- **self-adjoint / reversible** — a system that runs the same forwards and backwards
  (detailed balance); only then is τ certified.
- **L²(π)** — the "space" (weighting) in which a reversible system counts as
  self-adjoint; you never need it unless you're on the certified-τ route.
- **γ_O (SOC)** — a stability indicator of the peak (higher = a healthier, cleaner peak).

## Relation to prior work

σ_c selects a scale as the **argmax of a normalised derivative over a (log-)scale
dial** — the same shape of idea appears across fields, and σ_c is not claiming to
invent it:

- **Scale-space / computer vision (Lindeberg, 1998).** Automatic scale selection
  picks the scale maximising a *γ-normalised* differential response — structurally
  the closest neighbour to σ_c. The difference is not the argmax; it is the
  **discipline around it**: σ_c ships a 4-code applicability verdict, a resolution
  band, a blindness map, and refuses (constant/noise/convention-fragile inputs) —
  Lindeberg's construction is the selection rule, σ_c is the selection rule *plus the
  honesty layer that says when not to trust it*.
- **Loss/dielectric spectroscopy (Debye).** The loss peak against log ω is a σ_c
  reading; the paper cites Debye directly.
- **Systems / control.** A Bode corner frequency `1/RC` is a σ_c on a frequency dial
  (a known-answer check).

The contribution of this kernel is the *reader*, not the peak-finder: same argmax,
but every number carries what it is worth and what it cannot see.

## Known limitations (honest, and carried in the register)

These are stated here and in [`sigma_c/THEOREM_MAP.md`](sigma_c/THEOREM_MAP.md) /
[`sigma_c/LIMITS.md`](sigma_c/LIMITS.md), not buried:

- **A single time series gives no τ.** The most common ask that this tool refuses;
  see the box above. It needs an operator, a verified spectrum, or two probes.
- **The jitter / second channel is experimental and uncertified.**
  `sigma_c.experimental.jitter` emits **no** certified verdict (`validated=False`);
  it is not imported by the kernel and must not back a certified claim.
- **The finite-window τ gate (`cor:twoprobe`) is `GAP-KNOWN`.** The tail-window
  condition is a located, documented gap in the proof as written — the register
  marks it, and the code never renders it as `PROVED`.
- **The Lean proofs are machine-checked but not independently reviewed.** They cover
  only the jitter *algebra* (App B.3/B.4), live in the source repo (not the PyPI
  package), and stay out of any tag until a human other than the author has read them.
- **Dense τ-route only.** No sparse path; time O(n³), memory O(n²) — keep the operator
  dimension in the low thousands ([`LIMITS.md`](sigma_c/LIMITS.md) has measured numbers).

## Documentation

- [`sigma_c/README.md`](sigma_c/README.md) — kernel scope, the 4-code contract,
  the certified-τ route, what is and is not in scope.
- [`sigma_c/examples/README.md`](sigma_c/examples/README.md) — the learning ladder.
- [`sigma_c/BYO_GUIDE.md`](sigma_c/BYO_GUIDE.md) — bring your own data.
- [`sigma_c/THEOREM_MAP.md`](sigma_c/THEOREM_MAP.md) — label → paper-number binding + proof-status.
- [`sigma_c/LIMITS.md`](sigma_c/LIMITS.md) — honest bounds: numerical resolution floor, dense O(n³)/O(n²) τ-route cost (measured), and what the tool refuses to do.

## Papers

> *The Parrot's Theorems* (2026). **Preprint on Zenodo.**
> doi:[10.5281/zenodo.22066713](https://doi.org/10.5281/zenodo.22066713)

The citation is single-sourced in [`CITATION.cff`](CITATION.cff) and
`sigma_c.__paper__` — a journal publication is expected to follow this preprint, and
this kernel is one part of the larger ForgottenForge σ_c / Parrot programme.

Applied validation on real quantum hardware:

> *Operational scale detection in quantum magnetism* (2026).
> AVS Quantum Science **8** (1), 013804. doi:[10.1116/5.0312410](https://doi.org/10.1116/5.0312410)

## License

Dual-licensed: [AGPL-3.0-or-later](LICENSE) for open-source/non-commercial
use, or a commercial license (see [`LICENSE-COMMERCIAL.md`](LICENSE-COMMERCIAL.md) /
[nfo@forgottenforge.xyz](mailto:nfo@forgottenforge.xyz)). See [`LICENSING.md`](LICENSING.md).

## Acknowledgements

Developed at **ForgottenForge** · <nfo@forgottenforge.xyz> · <https://forgottenforge.xyz>

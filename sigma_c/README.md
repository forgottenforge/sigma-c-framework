# sigma_c — the disciplined-reader kernel

The reference implementation of the susceptibility scale-selector σ_c from

> *The Parrot's Theorems* — a Zenodo preprint.
> doi:[10.5281/zenodo.22066713](https://doi.org/10.5281/zenodo.22066713)

> **One rule:** every output either cites a theorem or admits it cannot — and
> both are visible from the outside.

> *Import path: `sigma_c`.*

## What σ_c is (in one breath)

σ_c is the scale on your dial where an observable is **most susceptible to a
change of scale** — the peak of `chi_O(sigma) = |sigma * dO/dsigma|`. It is a
**location read under a declared convention**, not a fitted parameter and not a
predictive bound. The kernel's job is not to hand you a number and stop; it is to
hand you a number **with what that number cannot see**.

## Install / use

```
pip install sigma-c-framework      # >= 6.0.0
```

or, from a clone, editable: `pip install -e .` from the repository root.

```python
import numpy as np
from sigma_c import analyze, gamma_k, Framework

sigma = np.geomspace(0.05, 50.0, 400)
O = np.exp(-sigma / 3.0)

result = analyze(sigma, O, window=gamma_k(2),
                 framework=Framework.REVERSIBLE_MARKOV)

print(result.summary())     # human-readable, self-explaining page
d = result.to_dict()        # machine-readable, statuses + band derived at dump time
result.card("out.png")      # optional PNG card
```

Pass the observable either as an array `O` on the `sigma` grid, as a callable
`O(sigma)`, or as a precomputed susceptibility `chi=`. Windows carry an **analytic**
profile constant `rho_star`: `bare()`, `gamma_k(k)`, `exponential(...)`,
`log_gaussian(...)` (see `WINDOW_REGISTRY`).

## The three things every output carries

**1. A per-output 4-code applicability verdict.** Each quantity (`result.sigma_c_status`,
`result.tau_abscissa_status`, `result.window_readability_status`) is one of:

| code | meaning for you |
|---|---|
| `OK` | trustworthy under the method's preconditions |
| `NOT_APPLICABLE` | a precondition is violated — cannot lift with these data |
| `NOT_RESOLVABLE` | below the method's resolution — cannot lift with these data |
| `NOT_IDENTIFIED` | a **supplyable** precondition is missing — the verdict names what to bring (a two-probe test, a verified spectrum) |

The verdicts are **derived** at read/serialize time from the primary fields, never
stored, so they cannot go stale.

**2. A resolution band on σ_c (NOT a confidence interval).** `result.resolution_band()`
reports how far the reported σ_c can move with the data held fixed: the local
**grid-cell** resolution (which also fixes how many significant digits it is honest
to print) and the **per-input** range of `min_prominence_ratio` over which the
regime verdict is unchanged (the flip points are read off *this* observable, not a
hardcoded threshold). There is no sampling model here, so there is deliberately no
CI — that would be a fabricated error bar.

**3. A blindness map / null cone.** `result.blindness()` places each reported
quantity on the paper's invariance ladder (Layer 1 counts / Layer 2 rates / Layer 3
locations) and states the corank hard edge: a degree of freedom that neither turns
with the dial nor rustles into the observable is invisible to this pairing forever —
no reprocessing recovers it, only a new instrument (a second observable, or the
dial-free jitter / two-probe channel) does. It reports **no corank number** from a
single scan, and says so.

## τ: two honest fields — the load-dominant bridge and the certified abscissa

`result.tau_bridge` = `sigma_c / rho_star` is the **load-dominant** read-out: it reads
the *most-susceptible* mode, which is the slowest one only under single-mode
faithfulness — a diagnostic, **not** a certified relaxation time. `result.tau_abscissa`
is the certified one (below). (Why `sigma_c = rho_star * tau` is a convention factor
rather than an identity is worked out in *The Parrot's Theorems*.)

`result.tau_abscissa` = `-T_star/log(lambda_2)` (in TIME) is the **certified spectral
abscissa** — but only for a **self-adjoint** operator, which is the live-route theorems'
real precondition. Self-adjointness lives in the inner product you **declare**: a
reversible Markov chain is self-adjoint in `L^2(pi)` (detailed balance), *not* symmetric
in euclidean `l^2`.

```python
import numpy as np
from sigma_c import analyze, bare, Framework

sigma = np.geomspace(0.05, 50.0, 400)          # a dial
O = np.exp(-sigma / 3.0)                        # a single-mode relaxation observable
P = np.array([[0.5,  0.5,  0.0],     # a reversible chain (pi = [1,2,1]/4):
              [0.25, 0.5,  0.25],    # self-adjoint in L^2(pi), NOT symmetric in l^2
              [0.0,  0.5,  0.5]])
result = analyze(sigma, O, window=bare(),
                 framework=Framework.REVERSIBLE_MARKOV,
                 operator=P,
                 inner_product="auto",         # verify detailed balance, derive pi, symmetrise
                 gamma_A=0.5,                   # contamination-to-signal (for the window gate)
                 T_star=1.0,                    # required: tau is in TIME
                 sigma_axis="evolution_time")   # declares the sigma-axis (fixes the window T)
# result.selfadjoint == "verified"
# result.tau_abscissa_status.code == OK        (thm:abscissa: certified operator property)
# result.window_readability_status.code        (thm:auto: is the rate readable in the window?)
```

Two verdicts, because *"τ is a property of the operator"* and *"τ is readable in a finite
window"* are different claims:

- **`tau_abscissa_status`** (thm:abscissa): OK iff the operator is self-adjoint in the
  declared inner product, irreducible, and has a decaying subdominant mode — **no strict
  gap required**. A non-self-adjoint / non-reversible operator → `NOT_APPLICABLE` (the
  theorems need self-adjointness, not mere normality — a drift walk is normal but not
  reversible). A **bare spectrum** (no operator) is `selfadjoint="declared"`: OK on the
  *declared* hypotheses, not verified.
- **`window_readability_status`** (thm:auto): OK additionally needs a strict spectral gap
  and a window `T > max(t*, 1/|mu_3|)`, and reports `result.rate_error_bound`. A bare
  spectrum can never be window-OK (no probe → no Gamma_A) → `NOT_IDENTIFIED`.

`sigma_axis` declares what the σ-axis is (it fixes the window `T` for the gate):

| `sigma_axis` | σ-unit | window T for the gate | bridge↔abscissa comparison |
|---|---|---|---|
| `evolution_time` | time (T_star units) | `T = max(sigma-grid)` | literal |
| `window_time` | time | `T = T_obs` (required) | suppressed (instrument-side) |
| `other` | not time | `T = T_obs` (required) | none (different units) |

## Experimental — the jitter channel (dial-free)

`sigma_c.experimental.jitter` implements the **jitter theorem** of *The Parrot's
Theorems* (§5, App B.3/B.4/C). With every dial **clamped**, you take repeated block
readings of several instruments; a slowly wandering latent parameter shifts each
block mean by `a_i·δ`, so the block-to-block covariance is `C = Var(δ)·a aᵀ`. From
that: (i) **detection** and (ii) **counting** (rank above a shuffle-null), (iii)
**identification** (the leading eigenvector ∝ the dial susceptibilities `a`; the
paper's `cos(v₁,a) ≥ 0.9`), and (iv) a **cross-check** (the jitter slope
`Cov(M₁,M₂)/Var(M₂)` equals the dial slope `dM₁/dM₂` under detailed balance — the
real FDT test). It confirms **rank** and **peak location** independently of the dial;
it delivers **no certified τ** (peak→τ is the fallen bridge).

This subpackage is **not part of the certified route**: the kernel does not import
it, it emits no certified `OK` verdict, and everything returns `validated=False`. Its
pre-registered validation (App C — see `experimental/PREREG_jitter.md`) has been run
from outside across four independent rounds, all App C criteria passing
(`experimental/VALIDATION_RESULTS.md`); promotion out of `experimental/` still awaits
the professor's certifying read. Import it explicitly to experiment:

The App C two-channel run is complete: **Channel 1** (`pair_parity_rank`) reads the
per-shot covariance of co-measured pair parities (fast), **Channel 2** (`jitter_rank`)
the drift covariance of block means (slow) — same detection/counting machinery.

```python
from sigma_c.experimental.jitter import (
    jitter_rank, pair_parity_rank, tetrad, leading_coupling, cos_coupling,
    jitter_slope, two_channel_fdt, peak_agreement,
)
```

## Tests

```
python -m pytest sigma_c/tests/ -q
```

Golden values from the paper's `verify_examples.py` plus enforcement invariants,
the 4-code contract, the certified τ-route + its adversarial regressions, and the
resolution band.

## Examples

```
PYTHONPATH=. python sigma_c/examples/01_coffee_cooling.py
PYTHONPATH=. python sigma_c/examples/02_zipf_wealth.py
PYTHONPATH=. python sigma_c/examples/03_two_windows.py
PYTHONPATH=. python sigma_c/examples/04_bimodal_relaxation.py
PYTHONPATH=. python sigma_c/examples/05_ising_chain.py
```

Each writes a PNG card next to the script; see `examples/README.md` for the
regime-coverage table.

## What this kernel does NOT do (out of scope, on purpose)

- **User-facing UI / web app.** Delegated to a separate product layer; the kernel
  ships as a Python library so that every output maps to a paper theorem.
- **A confidence interval on σ_c.** No sampling model → no CI; the resolution band
  is the honest substitute.
- **A corank number from a single scan.** Needs the model's sensitivity matrix.
- **Silent detrending** of power-law backgrounds. Never on by default (paper §C.3).
- **Unconditional Aczél canonicity** (open functional-equation problem, paper §11.3).

## Map

```
sigma_c/
├── api.py            analyze(), two_probe_test()
├── framework.py      the standard frameworks (paper Prop 5.4)
├── windows/          canonical windows with analytic rho_star
├── core/
│   ├── susceptibility.py   chi_O, peak-finding, sub-grid refinement
│   ├── trichotomy.py       regime classifier (peak-count + gap-boolean + floor)
│   ├── reversible.py       certified tau route: self-adjointness in L^2(pi), eigh
│   ├── resolution.py       the resolution/sensitivity band on sigma_c
│   ├── stability.py        gamma_O strict-SOC indicator (Prop 4.4)
│   └── faithfulness.py      F1/F2/F3 checkers + tail-window t*
├── codes.py          the 4-code applicability contract (per output)
├── result.py         Result, Trichotomy, BlindnessMap, TwoProbeResult
├── card.py           card renderer
├── theorem_map.py    cite() — label -> paper number + proof-status
├── examples/         the learning ladder (see examples/README.md)
└── tests/            the golden + contract + regression suite
```

See `THEOREM_MAP.md` (shipped in the package) for the label → paper-number + proof-status binding.

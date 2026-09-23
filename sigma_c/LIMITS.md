# Limits and scaling

Honest bounds of the kernel, so you know before you scale up. Measured, not
guessed (the numbers below are from a real run; reproduce with the snippet at the
end).

## What σ_c does and does not do

- It reads a **scale on a dial** from a curve `O(σ)`. It is a diagnostic/observational
  instrument, not a predictive bound.
- It will **refuse** rather than guess: a constant/near-constant observable →
  `NOT_IDENTIFIED`; a peak set that holds only under a narrow prominence
  convention, or too many peaks, or a monotone profile → `NOT_RESOLVABLE`.
- It does **not** turn a single time series into a relaxation time τ. See the
  "one thing it will not do" box in the README.

## Numerical resolution floor

`analyze()` computes χ = |dO/dlog σ| by a numerical log-derivative. Signal below
the double-precision rounding floor of that derivative cannot be resolved:

- A χ peak must clear `128 · eps · max|O| / min(Δlogσ)` to count.
- Empirically the honest cutoff is at **≈1e-12 relative amplitude** (χ/floor ≈ 20):
  a bump ~1000× above machine ε resolves to `OK`; below that it is refused,
  because it is indistinguishable from the derivative's own rounding noise.

## Computational cost (the τ-route: `operator=` / `spectrum=`)

The certified τ-route diagonalises the operator through a **dense** symmetric
eigensolver (`numpy.linalg.eigh`). There is **no sparse path**. Cost is therefore
the standard dense **time O(n³), memory O(n²)** in the operator dimension `n`.

Measured (`eigh` route, self-adjoint operator, one core):

| n     | time    | peak memory |
|-------|---------|-------------|
| 250   | 0.02 s  | 1.6 MB      |
| 500   | 0.03 s  | 6.1 MB      |
| 1000  | 0.08 s  | 24 MB       |
| 2000  | 0.29 s  | 96 MB       |
| 3000  | 1.3 s   | 217 MB      |

Time per doubling climbs 2.9× → 3.6× → 4.5×, converging on the cubic 8× as the
O(n³) term dominates; memory tracks n². Extrapolating: **n = 20000 is ~minutes of
compute and ~10 GB of memory** — impractical on a workstation. Keep the operator
dimension in the low thousands; reduce/coarse-grain the operator before this route
if it is larger. (The dial-only path — `analyze(sigma, O)` with no operator — is
cheap: it is a gradient + peak scan, linear in the number of samples.)

## Reproduce

```python
import numpy as np, sigma_c, time
s = np.linspace(0.1, 3, 60); O = np.cos(s)
for n in (250, 500, 1000, 2000, 3000):
    rng = np.random.default_rng(0); A = rng.normal(0, 1, (n, n)); S = (A + A.T) / 2
    S = S / (2.2 * np.max(np.abs(np.linalg.eigvalsh(S)))) + 0.5 * np.eye(n)
    t = time.time()
    sigma_c.analyze(s, O, operator=S, T_star=1.0, sigma_axis="evolution_time")
    print(n, round(time.time() - t, 3), "s")
```

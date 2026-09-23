# Bring your own data (Stufe 5)

You have a measurement of your own. Here is how to point the tool at it — and, just
as important, how to read whether you may trust what comes back.

## Step 1 — name your two axes (no code yet)
- **The dial `sigma`.** What did you vary? A length, a temperature, a probability, a
  time-window, a resolution. It must be positive and you sweep it over a range.
- **The observable `O`.** What did you measure at each dial value? One number per
  `sigma`. (An array the same length as `sigma`, or a function `O(sigma)`.)

If you cannot name a single dial and a single measurement, the tool is not for your
question yet — that is fine, it just means there is a modelling step first.

## Step 2 — the decisions (these are yours, not the tool's)
1. **Window (probe shape) — NOT a smoother.** `bare()` or `gamma_k(2)` picks the
   *probe shape* that fixes the analytic profile constant `rho_star`. It does **not**
   smooth your data. If `O` is noisy, the tool may split into many spurious peaks
   (regime II, a vector `sigma_c`) — that is noise, not many real scales. For noisy
   data: **pre-smooth `O` yourself** before calling, and/or raise
   `min_prominence_ratio` (below) so small bumps don't count. A large peak count is
   itself a warning sign.
2. **Noise threshold `min_prominence_ratio`** (default 0.10): a bump counts as a peak
   only if it rises to this fraction of the tallest one. Raise it (e.g. 0.3–0.5) on
   noisy data to keep the noise from being read as structure. This is a declared
   convention — the `resolution_band` reports the range over which your verdict is
   stable to it.
3. **The dial list.** Is `sigma` really the only thing you changed? Anything else
   that moved is a hidden dial and will confuse the reading.
4. **Same quantity across the sweep.** `O` must mean the same thing at every `sigma`
   (not "accuracy" here and "speed" there).
5. **Preprocessing.** If you filtered/normalised `O` with something that carries an
   absolute scale, pass `preprocessing_scale_equivariant=False` — the result is then
   flagged exploratory.
6. **A spectrum?** Only if you actually have a transfer operator / its eigenvalues do
   you touch the τ-route (`operator=`, `T_star=`, `sigma_axis=`) — Stufe 4. Most users
   never need it.

## Step 3 — run it (five lines)
```python
import numpy as np
from sigma_c import analyze, bare, gamma_k

sigma = np.array([...])          # your dial, positive, ascending
O     = np.array([...])          # your measurement, one per sigma
result = analyze(sigma, O, window=gamma_k(2))
print(result.summary())
```

## Step 4 — read the verdict BEFORE the number
`result.summary()` prints `sigma_c` with its **resolution band** (how tightly the
location is pinned — NOT a confidence interval) and a **verdict**:
- **OK** — the number stands.
- **NOT_RESOLVABLE** — your data isn't enough for this.
- **NOT_APPLICABLE** — wrong tool for this question.
- **NOT_IDENTIFIED** — something from you is missing; it says what to bring.

A `NOT_*` verdict is not a failure — it is the tool being honest about what your data
can and cannot say. Trust `sigma_c` only when its status is `OK`, and even then read
the band before you quote digits.

## A ready-to-edit skeleton
Copy `examples/TEMPLATE_adapter.py`, replace the two marked lines with your own
`sigma` and `O`, and run it.

## Where this stops
This gets you a *verified reading* of a scale in your data. Turning that into a new
scientific claim is your work, not the tool's — see "finding something new" in
`examples/README.md`.

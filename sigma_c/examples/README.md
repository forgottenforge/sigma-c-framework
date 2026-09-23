# Learning ladder — start at the bottom

These demos teach the tool **step by step**. Do them in order. The first ones need
no jargon at all; the words like "self-adjoint" only show up near the top, and are
explained when they first appear. If you are curious and fifteen, this is for you.

> **New here?** [`EXPLAINED.md`](EXPLAINED.md) walks through every example card in
> plain words: what each `.py` does, what the card shows, what each number means,
> and what you may conclude — including the certified-τ card (`card_06`).

**What the tool does, in one sentence:** you have a curve (some measurement `O`
plotted against a dial `sigma` — the dial is *whatever you turned*: a length, a
temperature, a probability, a time); the tool finds the dial value where the curve is
*most sensitive* (changes most steeply) — the scale where "something happens" — and
tells you how tightly that location is pinned, or refuses if it cannot honestly say.

Every demo starts with three plain sentences: **what** is measured, **which number**
should come out, and **where that number is known from**. Then it runs.

## The ladder

| Stufe | Demo | What you learn | Known answer |
|---|---|---|---|
| **1 · Hello world** | `demo_raw_shot.py` | we hide a number in a curve; the tool finds it again and says how precisely | we build the curve with 4; the tool returns 4.00 |
| **2 · Find a textbook number** | `demo_percolation.py` | the tool recovers a value someone proved exactly | percolation threshold p_c = 1/2 (marches onto 1/2 as the grid grows) |
| **3 · The tool says "no"** | `demo_no_scale.py` | why *"there is nothing to find here"* is a good, honest answer | a curve with no bump → NOT_RESOLVABLE |
| **3 · A bump is not a discovery** | `demo_seismic.py` | a *real-world* "no": the tool finds a bump, but it is your instrument's edge, not the world | earthquake catalogue: σ_c lands on the completeness Mc≈2.0, **not** the physics (b≈1.0 comes from the standard tool above Mc) |
| **4 · For the advanced** | `demo_reversible_tau.py` | the certified relaxation time τ (and the word *self-adjoint*), with the drift-walk as the refusal for a case that needs this level | τ = −T*/log λ₂ (reversible chain); drift walk → NOT_APPLICABLE |
| **5 · Your own data** | `../BYO_GUIDE.md` + `TEMPLATE_adapter.py` | bring your own measurement and read the verdict before you trust the number | — |

Additional worked examples (older, by roughly the same level): `01_coffee_cooling.py`
(Stufe 1–2), `05_ising_chain.py` (1D correlation length, Stufe 2), `02_zipf_wealth.py`
(another "there is no single scale here" answer, Stufe 3), `03_two_windows.py` (two
probes, Stufe 4), `04_bimodal_relaxation.py` (two scales at once, Stufe 3–4).

## The four verdicts, in plain words

Every number the tool returns carries one of these — read it *before* you trust the number:

- **OK** — the number stands.
- **NOT_RESOLVABLE** — your data isn't enough for this (below what can be resolved).
- **NOT_APPLICABLE** — wrong tool for this question (a precondition is violated).
- **NOT_IDENTIFIED** — something from *you* is missing; the verdict says what to bring.

## Run them

```
PYTHONPATH=. python sigma_c/examples/demo_raw_shot.py        # Stufe 1
PYTHONPATH=. python sigma_c/examples/demo_percolation.py     # Stufe 2
PYTHONPATH=. python sigma_c/examples/demo_no_scale.py        # Stufe 3 (a plain "no")
PYTHONPATH=. python sigma_c/examples/demo_seismic.py         # Stufe 3 (a bump is not a discovery)
PYTHONPATH=. python sigma_c/examples/demo_reversible_tau.py  # Stufe 4
PYTHONPATH=. python sigma_c/examples/TEMPLATE_adapter.py     # Stufe 5 (edit it)
```

## Not on the ladder: finding something *new*

Recovering known values (Stufe 1–4) is "hello world" — it shows the instrument
measures what it should. **Using it to discover something not yet known is a separate
thing, and it is research, not a recipe.** There is no ladder for it here, on purpose:
the tool is a *verification* instrument (it confirms and it refuses), not a discovery
oracle. If you go there, you are doing science, and you own the claims.

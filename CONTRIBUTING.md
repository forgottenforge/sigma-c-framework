# Contributing to sigma-c-framework

Thank you for your interest. The project is dual-licensed
(AGPL-3.0-or-later OR a commercial license — see [`LICENSING.md`](LICENSING.md))
and is maintained by ForgottenForge.

> **Single maintainer.** This project currently has one maintainer. Reviews and
> merges may take time; please be patient.

## Before you start

- Read the **[Code of Conduct](CODE_OF_CONDUCT.md)** — it applies to all
  interactions on this repository.
- Skim the **[README](README.md)** to understand what the kernel does.
- Every public surface must map to a labelled theorem/proposition/definition in
  *The Parrot's Theorems* (Zenodo doi:10.5281/zenodo.22066713); see
  [`THEOREM_MAP.md`](sigma_c/THEOREM_MAP.md).

## How to contribute

### Bug reports

Use the **Bug report** issue template. Include:

- version (`python -c "import sigma_c; print(sigma_c.__version__)"`)
- operating system and Python version
- a minimal reproducible example (≤ 30 lines)
- expected vs observed behaviour

### Feature requests

Use the **Feature request** issue template. Specify which theorem/proposition the
new surface would cite. **If there is no paper anchor for a proposed feature, the
discipline is to not ship it** — open a discussion instead.

### Pull requests

1. Fork the repo and create a feature branch from `main`:
   `git checkout -b feat/short-description`
2. Make your change. Keep it focused — one concern per PR.
3. Run the test suite locally: `python -m pytest sigma_c/tests/ -v`
4. If you add a new public surface, add a `cite("label:name")` reference in its
   docstring pointing at the paper label that backs it, and update
   `THEOREM_MAP.md` if the label is new.
5. Commit with a clear message (see *Commit style*).
6. Open the PR using the pull-request template.

### Commit style

- One concern per commit.
- First line: short imperative ≤ 72 chars, e.g.
  `core/faithfulness: add explicit C_R bound for F3 condition`.
- Body (optional): why, not what.
- **Sign off every commit** (`git commit -s`, adding a `Signed-off-by:` line) —
  this is your Developer Certificate of Origin (DCO) attestation.
- **Do not add `Co-Authored-By:` trailers or tool-attribution links for AI
  assistants.**

### Tests

- New surfaces require a test in `sigma_c/tests/` that either reproduces a
  paper-printed numerical value (golden test) or asserts a discipline invariant
  (provenance, regime, citation). Tests must use tolerances, not exact float
  equality (results vary slightly across BLAS/platforms).
- Do not commit tests that rely on private datasets.

### Code style

- Python: PEP 8, line length 100. Consistency with the surrounding file is
  appreciated.
- Type hints encouraged on public surfaces.
- Docstrings: short, paper-anchored where applicable. Cite via
  `cite("thm:foo")` from `sigma_c.theorem_map`.

## Development setup

```bash
git clone https://github.com/forgottenforge/sigma-c-framework
cd sigma-c-framework
python -m venv .venv
source .venv/bin/activate         # or .venv\Scripts\activate on Windows
pip install -e ".[dev]"
python -m pytest sigma_c/tests/
```

## Reporting security issues

Do not open a public issue for security problems. See
[`SECURITY.md`](SECURITY.md) for responsible-disclosure instructions.

## License of contributions

By submitting a contribution you agree that it is licensed under the project's
dual-license model (see [`LICENSING.md`](LICENSING.md)). Because of the commercial
option, substantial contributions require:

1. a **DCO sign-off** on every commit (`git commit -s`), and
2. a signed **Contributor License Agreement (CLA)** before merging, so the
   project can continue to offer both the AGPL and the commercial license.

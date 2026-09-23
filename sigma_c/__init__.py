# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
sigma_c — The disciplined reader for operational scale selection.

Built on *The Parrot's Theorems* (a preprint on Zenodo
doi:10.5281/zenodo.22066713). THEOREM_MAP.md binds every code citation
to a paper label and its proof-status.

The hero call:
    from sigma_c import analyze, gamma_k, Framework
    result = analyze(sigma, O_values, window=gamma_k(2),
                     framework=Framework.REVERSIBLE_MARKOV)
    print(result.summary())
    result.card("out.png")

Every output either cites a theorem (via THEOREM_MAP) or admits it cannot.
Both are visible from the outside.
"""
from sigma_c.api import analyze, two_probe_test
from sigma_c.result import Result, Regime, Trichotomy
from sigma_c.framework import Framework
from sigma_c.windows import (
    bare,
    gamma_k,
    exponential,
    log_gaussian,
    Window,
    WINDOW_REGISTRY,
)
from sigma_c.core.faithfulness import (
    check_F1,
    check_F2,
    check_F3,
    kl_modal_coefficients,
    FaithfulnessCheck,
)
# The applicability contract, public: the four codes + the per-output verdicts.
# result.sigma_c_status / result.tau_status return these; result.to_dict()
# serializes them (derived, never stored).
from sigma_c.codes import (
    Code,
    Verdict,
    sigma_c_verdict,
    tau_two_probe_verdict,
)
from sigma_c.theorem_map import cite
from sigma_c.bootstrap import bootstrap_sigma_c
from sigma_c.errors import SigmaCError, InvalidInputError
from sigma_c.profile_tau import tau_from_profile, ProfileTauResult
# The second (dial-free) channel is EXPERIMENTAL and deliberately NOT imported
# here: `from sigma_c.experimental.jitter import ...` (see experimental/PREREG_jitter.md).

__version__ = "6.0.0"


def __getattr__(name: str):
    """Helpful migration hint for names that no longer exist at the top level.

    `import sigma_c` changed meaning from the 5.x line (which exposed an
    adapter/utility stack) to the 6.0 kernel. A 5.x name lookup lands here."""
    raise AttributeError(
        f"module 'sigma_c' has no attribute {name!r}. If this is a name from the "
        f"5.x line, the public API changed in 6.0: 'import sigma_c' is now the "
        f"disciplined-reader kernel (analyze, Result, the 4-code Verdict, ...). See "
        f"the 'Migrating 5.x -> 6.0' section of CHANGELOG.md; the SDK adapters now "
        f"live under sigma_c.adapters.* as optional extras (pip install "
        f"'sigma-c-framework[qiskit]' etc.)."
    )

# ---------------------------------------------------------------------------
# The paper — SINGLE SOURCE OF TRUTH for the publication reference.
# Everything else (the aliases below, CITATION.cff, the docs' citations) points
# here. What changes when the journal accepts the paper is `status` and `venue`
# (add the journal volume/issue/pages and, if a journal DOI is minted, `journal_doi`);
# the title, authors and Zenodo DOI do not change. Update this record ONCE on
# publication, then flip the two PREPRINT rows in THEOREM_MAP.md to PUBLISHED and
# refresh CITATION.cff — those are the only places the status is stated.
__paper__ = {
    "title": "The Parrot's Theorems",
    "authors": "ForgottenForge",
    "year": 2026,
    "status": "preprint",              # -> "published" on acceptance
    "venue": "Zenodo",                 # -> journal name/vol/issue/pages on publication
    "zenodo_doi": "10.5281/zenodo.22066713",            # stable across publication
    "journal_doi": None,                               # minted on journal publication
    "url": "https://doi.org/10.5281/zenodo.22066713",
    "note": ("Part of the ForgottenForge sigma_c / Parrot research programme; the "
             "fuller hub (Parrot 2) is in preparation. A journal publication may "
             "follow this preprint."),
}
# Derived (single-sourced) aliases kept for backward compatibility.
__paper_doi__ = __paper__["zenodo_doi"]
__paper_url__ = __paper__["url"]
__paper_version__ = f"{__paper__['title']} ({__paper__['venue']}), {__paper__['year']}"

__all__ = [
    "analyze",
    "two_probe_test",
    "tau_from_profile",
    "ProfileTauResult",
    "bootstrap_sigma_c",
    "Result",
    "Regime",
    "Trichotomy",
    "Framework",
    "Window",
    "WINDOW_REGISTRY",
    "bare",
    "gamma_k",
    "exponential",
    "log_gaussian",
    "check_F1",
    "check_F2",
    "check_F3",
    "kl_modal_coefficients",
    "FaithfulnessCheck",
    "Code",
    "Verdict",
    "sigma_c_verdict",
    "tau_two_probe_verdict",
    "SigmaCError",
    "InvalidInputError",
    "cite",
    "__version__",
    "__paper__",
    "__paper_version__",
    "__paper_doi__",
    "__paper_url__",
]

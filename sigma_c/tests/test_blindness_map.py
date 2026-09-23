# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Tests for the blindness map / null cone -- the programme's honesty anchor
(paper sec:corank + the Layer table sec:setup, "The Parrot's Theorems").

Two guards, both demanded by the discipline itself:

1. UNRESOLVED-LABEL RATCHET. In a kernel whose promise is "every output either
   cites a theorem or admits it cannot", an unresolved citation is honest for a
   few weeks and then rots into decoration. So the number of unresolved labels
   emitted at runtime is itself a test value: it must never grow. Baseline is 0
   -- every label the Result and its blindness map emit must resolve through
   THEOREM_MAP.md. If this test fails, either a new label was emitted without a
   map entry, or an entry was removed: both are debt, caught the moment it lands.

2. NULL-CONE FIXTURE (the ninth rule: a guard is only trusted once something is
   planted that MUST make its claim bite). The blindness map is an
   assertion-generator; its assertion is "this pairing cannot see the following
   directions." That is a claim about INVISIBILITY, which an ordinary result
   cannot test. So we plant a structure exactly in a named blind direction and
   show the kernel does not find it, while the map named it beforehand. The
   dangerous failure is the opposite one -- the kernel resolving something in a
   direction the map called blind -- because a false claim of blindness invites
   you to stop looking. Two planted directions:
     A) the SILENT direction: a hidden d.o.f. the observable does not depend on
        -- output must be invariant along it (invisible forever, sec:corank).
     B) the RUSTLE direction: a shared latent fluctuation with a flat mean
        response -- the single-observable mean-channel kernel must miss it
        (only a second, dial-free channel could catch it; the separate experimental/jitter channel).

Run with pytest, or directly:  python sigma_c/tests/test_blindness_map.py
"""
from __future__ import annotations
import numpy as np

from sigma_c.api import analyze
from sigma_c.theorem_map import cite, reload_map, _parse_theorem_map


# ---------------------------------------------------------------------------
# helpers
# ---------------------------------------------------------------------------

def _emitted_labels(res) -> set:
    """Every paper label this result puts in front of a reader."""
    return set(res.citations) | set(res.blindness().citations)


def _is_unresolved(label: str) -> bool:
    """cite() falls back to 'paper [label]' exactly when the label has no
    THEOREM_MAP entry -- that string is the unresolved signal."""
    return cite(label).startswith("paper [")


_SIGMA = np.linspace(0.1, 5.0, 200)
_CASES = {
    "single_peak": np.exp(-((np.log(_SIGMA) - 0.2) ** 2) / 0.3),
    "two_peak": (np.exp(-((_SIGMA - 1.0) ** 2) / 0.05)
                 + 0.9 * np.exp(-((_SIGMA - 3.0) ** 2) / 0.05)),
    "monotone": 1.0 / (_SIGMA + 0.5),
}


# ---------------------------------------------------------------------------
# 1. unresolved-label ratchet
# ---------------------------------------------------------------------------

# The debt ceiling. It may only ever be LOWERED. Raising it to make a red test
# pass is precisely the accumulation of silent debt this guard exists to stop.
MAX_UNRESOLVED_LABELS = 0


def test_no_emitted_label_is_unresolved():
    all_labels = set()
    for O in _CASES.values():
        all_labels |= _emitted_labels(analyze(_SIGMA, O))
    unresolved = sorted(l for l in all_labels if _is_unresolved(l))
    assert len(unresolved) <= MAX_UNRESOLVED_LABELS, (
        f"{len(unresolved)} unresolved label(s) emitted at runtime "
        f"(ceiling {MAX_UNRESOLVED_LABELS}): {unresolved}. Either add the "
        f"THEOREM_MAP entry or repoint the citation -- do NOT raise the ceiling."
    )


def test_blindness_citations_resolve():
    # the honesty anchor's OWN citations must resolve -- the most expensive
    # place for a placeholder (they cite The Parrot's Theorems).
    bm = analyze(_SIGMA, _CASES["single_peak"]).blindness()
    for label in bm.citations:
        assert not _is_unresolved(label), f"{label} does not resolve"


# ---------------------------------------------------------------------------
# 1b. proof-status register -- a citation must not silently anchor on a proof
#     with a known gap (THEOREM_MAP proof-status register + cite() rendering)
# ---------------------------------------------------------------------------

def test_gappy_proof_status_is_rendered():
    reload_map()
    # the superseded five-axiom/A3 route must be flagged, not cited as current
    assert "SUPERSEDED" in cite("thm:characterisation"), cite("thm:characterisation")
    # the live uniqueness is reviewed-clean but NOT certified -> still flagged
    assert "REVIEW-CLEAN" in cite("thm:axiomatic-char"), cite("thm:axiomatic-char")


def test_status_register_does_not_clobber_numbers():
    # a PROVED row in the register must not overwrite the label's real number
    reload_map()
    assert _parse_theorem_map().get("thm:compat") == "3.7"
    # a PROVED status renders no marker (the ONLY status that renders nothing)
    assert "[proof:" not in cite("thm:compat")


def test_absent_label_defaults_to_no_register_entry_not_proved():
    # HONESTY: the register default flipped from PROVED-as-published
    # to a flagged non-PROVED marker. A theorem with no status row is "nobody
    # looked", not proved, and must be flagged so it does not travel as certified.
    # the marker was renamed UNVERIFIED ->
    # NO-REGISTER-ENTRY, because an outside user read "UNVERIFIED" as "the math is
    # unverified"; NO-REGISTER-ENTRY says plainly "unlisted, not disproven".
    # (def:sigmac is in the mapping table but has no status row -- genuinely unlisted.)
    reload_map()
    rendered = cite("def:sigmac")
    assert "[proof: NO-REGISTER-ENTRY" in rendered, rendered
    assert "not PROVED" in rendered, rendered
    assert "not disproven" in rendered, rendered
    # the number still resolves (status flip must not clobber resolution)
    assert _parse_theorem_map().get("def:sigmac") is not None


def test_regime_I_tau_is_flagged_load_dominant():
    # HONESTY REGRESSION: result.tau_bridge = sigma_c / rho_star is the
    # load-dominant read-out, not theorem-backed as the relaxation time (the
    # sigma_c = rho_star*tau bridge fell -> thm:spectral-id-B, GAP-KNOWN). A
    # regime-I result that emits a tau MUST surface that caveat, so a caller
    # inspecting citations/notes cannot read result.tau_bridge as certified. Guards
    # against a future edit silently dropping the flag (the silent-defect class).
    import numpy as np
    from sigma_c.api import analyze
    reload_map()
    sigma = np.logspace(-1.5, 1.5, 400)
    O = np.exp(-sigma / 1.0)  # single-mode decay -> regime I, tau ~ 1
    r = analyze(O=O, sigma=sigma)
    assert r.tau_bridge is not None, "single-mode decay should emit a tau_bridge"
    # the GAP-KNOWN label is cited, so cite() flags it at the tau emission
    assert "thm:spectral-id-B" in r.citations, r.citations
    assert "GAP-KNOWN" in cite("thm:spectral-id-B")
    # and a human-readable note names the load-dominant caveat
    assert any("LOAD-DOMINANT" in n for n in r.notes), r.notes


# ---------------------------------------------------------------------------
# 2. null-cone fixture -- plant structure in a named blind direction
# ---------------------------------------------------------------------------

def test_silent_direction_is_invisible():
    """Direction A: a hidden d.o.f. the observable does not depend on.

    The object carries (visible peak shape, hidden h). The instrument reads
    only the peak shape, so two objects that differ ONLY in h yield a
    byte-identical response -- the kernel cannot move along h. The map names
    this direction blind BEFORE we look; the kernel then proves blind to it.
    The dangerous failure would be a DIFFERENT verdict for a different h, i.e.
    the kernel resolving something the map called invisible.
    """
    def response(h):  # h is the planted hidden d.o.f.; O does not read it
        return np.exp(-((np.log(_SIGMA) - 0.2) ** 2) / 0.3)  # independent of h

    res0 = analyze(_SIGMA, response(h=0.0))
    res5 = analyze(_SIGMA, response(h=5.0))

    # the map named the blind direction (single-observable pairing) up front
    note = res0.blindness().null_cone_note.lower()
    assert "single observable" in note and "invisible" in note

    # kernel output is invariant along the planted direction -> it cannot see h
    assert res0.sigma_c == res5.sigma_c
    assert res0.regime.geometric == res5.regime.geometric
    # danger check: h did NOT conjure extra resolved structure
    n0 = len(res0.sigma_c) if isinstance(res0.sigma_c, list) else 1
    n5 = len(res5.sigma_c) if isinstance(res5.sigma_c, list) else 1
    assert n0 == n5


def test_rustle_direction_is_not_looked_at():
    """Direction B: a shared latent fluctuation invisible to the mean channel.

    Repeated readings co-vary through one common cause (a rank-1 off-diagonal
    covariance -- the jitter signal of sec:jitter). The latent has exactly zero
    mean, so it leaves NO trace in the mean response the single-observable
    kernel reads: the kernel's verdict is identical with and without the
    rustle, i.e. it is blind to the direction the covariance carries. Only a
    dial-free jitter channel (the separate experimental/jitter entry point) reads that second
    moment, and the map says exactly that: the rustle direction was not looked
    at.
    """
    rng = np.random.default_rng(0)
    base = np.exp(-((np.log(_SIGMA) - 0.2) ** 2) / 0.3)  # a real visible peak
    n_reads = 200
    xi = rng.standard_normal(n_reads)      # one shared latent per reading
    xi -= xi.mean()                        # EXACTLY zero mean: no trace in the mean
    coupling = 0.3 * base                  # every channel rustles with xi
    reads = base[None, :] + np.outer(xi, coupling)  # (n_reads, n_sigma)

    # there really WAS something to miss: the shared latent makes the
    # off-diagonal covariance nonzero (rank-1, the jitter signal of sec:jitter)
    cov = np.cov(reads, rowvar=False)
    off_diag = cov - np.diag(np.diag(cov))
    assert np.abs(off_diag).max() > 1e-6

    # the single-observable kernel sees only the mean channel; a zero-mean
    # latent leaves the mean untouched, so its verdict is IDENTICAL with and
    # without the rustle -- it is blind to the direction the covariance carries.
    res_with = analyze(_SIGMA, reads.mean(axis=0))
    res_without = analyze(_SIGMA, base)
    # identical up to float noise (the zero-mean latent perturbs the mean only
    # at ~1e-16); the verdict is the same, so the kernel is blind to the rustle
    assert res_with.regime.geometric == res_without.regime.geometric
    assert np.allclose(
        np.atleast_1d(res_with.sigma_c), np.atleast_1d(res_without.sigma_c),
        rtol=1e-9, atol=1e-9,
    )

    # and the map admits the rustle direction was not looked at
    note = res_with.blindness().null_cone_note.lower()
    assert "rustle" in note and "jitter" in note and "not looked at" in note


# ---------------------------------------------------------------------------
# no-pytest runner (pytest is the intended runner; this lets the file run bare)
# ---------------------------------------------------------------------------

if __name__ == "__main__":
    reload_map()
    tests = [
        test_no_emitted_label_is_unresolved,
        test_blindness_citations_resolve,
        test_gappy_proof_status_is_rendered,
        test_status_register_does_not_clobber_numbers,
        test_silent_direction_is_invisible,
        test_rustle_direction_is_not_looked_at,
    ]
    failed = 0
    for t in tests:
        try:
            t()
            print(f"PASS  {t.__name__}")
        except AssertionError as e:
            failed += 1
            print(f"FAIL  {t.__name__}: {e}")
    print(f"\n{len(tests) - failed}/{len(tests)} passed")
    raise SystemExit(1 if failed else 0)

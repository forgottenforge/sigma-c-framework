# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Shared dial contract (6.0.1): the three entry paths of ``analyze`` —
array-``O``, callable-``O``, and precomputed ``chi`` — must reject a malformed
dial *identically*. In 6.0.0 the array path was guarded by ``chi_O`` while the
chi and callable paths routed around it, so a non-positive, too-short, tied, or
mismatched dial was silently processed (reaching ``log(sigma<=0)`` or a
zero-spacing derivative) on exactly the paths that skipped ``chi_O``. These
tests pin every path to the same clean ``InvalidInputError``.
"""
import warnings

import numpy as np
import pytest

from sigma_c import analyze
from sigma_c.errors import InvalidInputError

SIGMA = np.linspace(0.5, 3.0, 7)


def _peak(s):
    return 2.5 * np.exp(-((np.log(s) - np.log(1.5)) ** 2) / (2 * 0.25 ** 2)) + 0.1


PEAK = _peak(SIGMA)

# The three entry paths, each as (sigma, values) -> Result. The callable path
# takes its observable from the function, so it ignores ``values`` by design.
PATHS = {
    "array": lambda sig, vals: analyze(sig, np.asarray(vals, dtype=float)),
    "callable": lambda sig, vals: analyze(sig, O=_peak),
    "chi": lambda sig, vals: analyze(sig, chi=np.asarray(vals, dtype=float)),
}

BAD_SIGMA = {
    "zero": np.r_[0.0, SIGMA[1:]],
    "negative": np.r_[-1.0, SIGMA[1:]],
    "too_few_points": SIGMA[:3],
    "nan": np.r_[np.nan, SIGMA[1:]],
    "inf": np.r_[np.inf, SIGMA[1:]],
    "duplicate": np.r_[SIGMA[0], SIGMA[0], SIGMA[2:]],
}


@pytest.mark.parametrize("path", list(PATHS))
@pytest.mark.parametrize("case", list(BAD_SIGMA))
def test_bad_sigma_rejected_on_every_path(path, case):
    """A malformed dial is rejected identically on all three entry paths."""
    sig = BAD_SIGMA[case]
    vals = np.ones(len(sig))
    with pytest.raises(InvalidInputError):
        PATHS[path](sig, vals)


def test_nonpositive_or_tied_sigma_rejected_without_numpy_warning():
    """The pre-check rejects before any log(sigma<=0) / zero-spacing divide, so
    no numpy RuntimeWarning is emitted on the way to the error."""
    with warnings.catch_warnings():
        warnings.simplefilter("error", RuntimeWarning)
        for sig in (np.r_[0.0, SIGMA[1:]], np.r_[SIGMA[0], SIGMA[0], SIGMA[2:]]):
            with pytest.raises(InvalidInputError):
                analyze(sig, chi=np.ones(len(sig)))


@pytest.mark.parametrize("bad_chi", [PEAK[:6], np.r_[PEAK, 0.05]])
def test_chi_length_mismatch_rejected(bad_chi):
    """chi shorter/longer than sigma is rejected (was IndexError / silent truncation)."""
    with pytest.raises(InvalidInputError):
        analyze(SIGMA, chi=bad_chi)


def test_O_array_length_mismatch_with_chi_rejected():
    with pytest.raises(InvalidInputError):
        analyze(SIGMA, O=PEAK[:6], chi=PEAK)


@pytest.mark.parametrize("path", list(PATHS))
def test_valid_dial_accepted_on_every_path(path):
    """The guard does not reject a well-formed dial on any path."""
    result = PATHS[path](SIGMA, PEAK)
    assert result is not None

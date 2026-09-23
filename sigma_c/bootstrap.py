# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
A STATISTICAL confidence interval for sigma_c from REPLICATE curves.

This answers a DIFFERENT question from `resolution_band`. The resolution band is
the ruler's finest distinguishable step (a discretisation/convention band, NOT a
CI). This is sampling variance: given several independent repeats of the same
measurement, how much does sigma_c move? It is a genuine statistical statement, so
it carries its own validation (a coverage check on synthetic replicates in
test_bootstrap.py: a nominal 95% interval covers the true sigma_c about 95% of the
time).

Method: a nonparametric bootstrap over the replicate curves. Each bootstrap draw
resamples the replicates with replacement, averages them into one curve, and reads
sigma_c; the percentile interval of that bootstrap distribution is the CI. Draws
whose averaged curve does not resolve to a scalar sigma_c are recorded and excluded
from the percentiles.
"""
from __future__ import annotations

from typing import Optional, Sequence, Union

import numpy as np


def bootstrap_sigma_c(
    sigma: Union[np.ndarray, Sequence[float]],
    replicates: Union[np.ndarray, Sequence[Sequence[float]]],
    *,
    window: Union[str, object, None] = "bare",
    n_boot: int = 1000,
    ci: float = 0.95,
    min_prominence_ratio: float = 0.10,
    seed: Optional[int] = None,
) -> dict:
    """Statistical CI for sigma_c from replicate curves.

    Parameters
    ----------
    sigma : (M,) array
        The shared dial grid.
    replicates : (R, M) array
        R independent repeat measurements of the observable on `sigma`.
    n_boot : int
        Number of bootstrap resamples.
    ci : float
        Central interval level (e.g. 0.95).

    Returns a dict with `sigma_c` (point estimate: sigma_c of the mean curve),
    `ci_lo`/`ci_hi` (percentile interval), `ci_level`, `n_boot`, `n_resolved`
    (bootstrap draws that yielded a scalar), `is_confidence_interval: True`, and
    `method`. If the mean curve itself does not resolve, `sigma_c` is None.
    """
    from sigma_c import analyze
    sigma = np.asarray(sigma, dtype=float)
    reps = np.asarray(replicates, dtype=float)
    if reps.ndim != 2 or reps.shape[1] != sigma.shape[0]:
        from sigma_c.errors import InvalidInputError
        raise InvalidInputError(
            f"replicates must be (R, M) matching sigma length M={sigma.shape[0]}; "
            f"got {reps.shape}."
        )
    R = reps.shape[0]
    rng = np.random.default_rng(seed)

    def _scalar_sigma_c(curve: np.ndarray):
        r = analyze(sigma, curve, window=window,
                    min_prominence_ratio=min_prominence_ratio)
        sc = r.sigma_c
        return float(sc) if isinstance(sc, (int, float)) else None

    point = _scalar_sigma_c(reps.mean(axis=0))
    draws = []
    for _ in range(int(n_boot)):
        idx = rng.integers(0, R, size=R)
        sc = _scalar_sigma_c(reps[idx].mean(axis=0))
        if sc is not None:
            draws.append(sc)

    lo_q = (1.0 - ci) / 2.0
    hi_q = 1.0 - lo_q
    if draws:
        ci_lo = float(np.quantile(draws, lo_q))
        ci_hi = float(np.quantile(draws, hi_q))
    else:
        ci_lo = ci_hi = None

    return {
        "quantity": "sigma_c",
        "sigma_c": point,
        "ci_lo": ci_lo,
        "ci_hi": ci_hi,
        "ci_level_nominal": ci,
        "n_boot": int(n_boot),
        "n_resolved": len(draws),
        "n_replicates": R,
        "is_confidence_interval": True,
        "approximate": True,
        "method": "nonparametric replicate bootstrap (resample repeats, average, read sigma_c)",
        "note": ("A statistical CI (sampling variance across your replicates) -- "
                 "distinct from result.resolution_band(), which is the ruler step, not "
                 "a CI. APPROXIMATE: the percentile bootstrap of a peak LOCATION "
                 "(an argmax) is known to UNDER-cover; validated coverage is ~90% at "
                 "nominal 95% when the curves resolve (test_bootstrap.py). Treat the "
                 "nominal level as indicative, and note that noisy replicates that do "
                 "not resolve to a scalar are excluded (n_resolved reports how many "
                 "of n_boot draws counted)."),
    }

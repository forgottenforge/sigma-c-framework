# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
#1+#2: the blindness map + 4-code contract as OUTPUT, derived not stored.

Result.to_dict() must DERIVE the blindness map and the applicability statuses at
serialization time from the primary fields — never carry a stored snapshot. The
discriminating test: serialize, change a convention,
serialize again — the map must follow. A blindness map that could itself go
stale would be worse than none. And the contract (Code/Verdict/verdicts) must be
importable from the top-level package.
"""
from __future__ import annotations
import dataclasses

import numpy as np

from sigma_c import (
    analyze, bare,
    Code, Verdict, sigma_c_verdict, tau_two_probe_verdict,
)


def _a_result():
    sigma = np.geomspace(0.05, 100.0, 400)
    O = np.exp(-sigma / 5.0)
    return analyze(sigma, O, window=bare())


def _min_prom_dep(d: dict):
    """The min_prominence_ratio convention-dep string carried in the blindness map."""
    for layer in d["blindness"]["layers"]:
        for dep in layer.get("convention_deps", []):
            if dep.startswith("min_prominence_ratio="):
                return dep
    return None


class TestOutputFieldDerived:

    def test_to_dict_carries_derived_blindness_and_statuses(self):
        d = _a_result().to_dict()
        # primary fields present
        assert "sigma_c" in d and "regime" in d
        # derived safeguards present and shaped
        assert "blindness" in d and "layers" in d["blindness"]
        assert "sigma_c_status" in d and d["sigma_c_status"]["code"] in {c.value for c in Code}
        assert "tau_abscissa_status" in d and d["tau_abscissa_status"]["code"] in {c.value for c in Code}
        assert "window_readability_status" in d and d["window_readability_status"]["code"] in {c.value for c in Code}

    def test_blindness_is_derived_at_serialization_not_stored(self):
        r = _a_result()
        before = _min_prom_dep(r.to_dict())
        assert before is not None and before != "min_prominence_ratio=0.999"

        # change a CONVENTION on the result, then re-serialize
        r.regime = dataclasses.replace(r.regime, min_prominence_ratio=0.999)
        after = _min_prom_dep(r.to_dict())

        assert after == "min_prominence_ratio=0.999", (
            "blindness map did not follow the changed convention: it was stored at "
            "analyze()-time, not derived at serialization"
        )
        assert after != before

    def test_status_is_derived_at_serialization_not_stored(self):
        # regime III => sigma_c_status NOT_RESOLVABLE; flip geometric and the
        # serialized status must change too (derived, not a stored snapshot).
        r = _a_result()
        code_before = r.to_dict()["sigma_c_status"]["code"]
        r.regime = dataclasses.replace(r.regime, geometric="III_geom")
        code_after = r.to_dict()["sigma_c_status"]["code"]
        assert code_after == Code.NOT_RESOLVABLE.value
        assert code_after != code_before

    def test_contract_is_publicly_exported(self):
        assert Code.OK.value == "OK"
        r = _a_result()
        assert isinstance(sigma_c_verdict(r), Verdict)
        assert isinstance(r.sigma_c_status, Verdict)
        assert isinstance(r.tau_abscissa_status, Verdict)
        assert isinstance(r.window_readability_status, Verdict)
        # tau_two_probe_verdict is the two-probe half of the contract
        assert callable(tau_two_probe_verdict)

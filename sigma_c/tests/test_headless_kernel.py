# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
The kernel is headless: `import sigma_c` and analyze() must not pull matplotlib
(it is the opt-in [plot] extra, needed only for result.card()). Integrator
guarantee: JSON in, JSON out, no display stack.
"""
import sys
import numpy as np

from sigma_c import analyze


def test_import_and_analyze_do_not_load_matplotlib():
    # matplotlib may already be imported by another test in this process; the
    # guarantee we can assert unconditionally is that analyze()'s own code path
    # does not require it. Prove it by making an import of matplotlib fail and
    # showing analyze() + to_dict() still work end to end.
    import builtins
    real_import = builtins.__import__

    def _blocked(name, *a, **k):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise ImportError("matplotlib blocked for this test")
        return real_import(name, *a, **k)

    saved = {m: sys.modules[m] for m in list(sys.modules) if m.startswith("matplotlib")}
    for m in saved:
        del sys.modules[m]
    builtins.__import__ = _blocked
    try:
        sigma = np.geomspace(0.1, 100, 400)
        O = np.exp(-sigma / 5.0)
        r = analyze(sigma, O)
        d = r.to_dict()
        assert abs(r.sigma_c - 5.0) < 0.5
        assert d["sigma_c_status"]["code"] == "OK"
        assert d["schema_version"] == 1
        assert r.summary()  # text summary works with no plotting stack
    finally:
        builtins.__import__ = real_import
        sys.modules.update(saved)


def test_card_without_matplotlib_gives_actionable_error():
    import builtins
    real_import = builtins.__import__

    def _blocked(name, *a, **k):
        if name == "matplotlib" or name.startswith("matplotlib."):
            raise ImportError("matplotlib blocked for this test")
        return real_import(name, *a, **k)

    saved = {m: sys.modules[m] for m in list(sys.modules)
             if m.startswith("matplotlib") or m == "sigma_c.card"}
    for m in saved:
        del sys.modules[m]
    builtins.__import__ = _blocked
    try:
        r = analyze(np.geomspace(0.1, 100, 400), np.exp(-np.geomspace(0.1, 100, 400) / 5.0))
        try:
            r.card("/tmp/should_not_be_written.png")
            raised = None
        except ImportError as e:
            raised = str(e)
        assert raised is not None, "card() must raise when matplotlib is absent"
        assert "[plot]" in raised and "matplotlib" in raised
    finally:
        builtins.__import__ = real_import
        sys.modules.update(saved)

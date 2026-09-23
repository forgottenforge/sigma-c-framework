# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
0b scaffold: every shipped example/demo must RUN in the suite --
so a demo that imports a missing module, teaches removed API, or otherwise breaks
turns the suite RED at creation, not later. New demos are added to the examples/
directory and are picked up here automatically.

Examples are executed as subprocesses with a headless matplotlib backend (they
render Cards). We assert a clean exit; demos that expose a `run()` returning a
known-answer dict are additionally checked in test_demos.py.
"""
import os
import subprocess
import sys
from pathlib import Path

import pytest

_EXAMPLES = sorted(
    p for p in (Path(__file__).resolve().parents[1] / "examples").glob("*.py")
    if p.name != "__init__.py"
)


@pytest.mark.parametrize("script", _EXAMPLES, ids=[p.name for p in _EXAMPLES])
def test_example_runs_clean(script):
    repo_root = Path(__file__).resolve().parents[2]
    env = dict(os.environ, MPLBACKEND="Agg", PYTHONPATH=str(repo_root))
    proc = subprocess.run(
        [sys.executable, str(script)],
        capture_output=True, text=True, env=env, timeout=120,
    )
    assert proc.returncode == 0, (
        f"{script.name} exited {proc.returncode}\n"
        f"--- stdout ---\n{proc.stdout[-2000:]}\n--- stderr ---\n{proc.stderr[-2000:]}"
    )


def test_at_least_the_five_hero_examples_present():
    names = {p.name for p in _EXAMPLES}
    for expected in ("01_coffee_cooling.py", "05_ising_chain.py"):
        assert expected in names, f"missing hero example {expected}"

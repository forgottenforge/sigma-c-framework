# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Opt-in adapters that feed the kernel from quantum SDKs.

Each adapter is a thin translator from a platform's parameter-sweep result into the
kernel's `(sigma, O)` input (the dial and the observable), so you can call
`analyze()` on it. They are OPTIONAL: the kernel does not import them, and each needs
its SDK installed as an extra:

    pip install "sigma-c-framework[qiskit]"      # sigma_c.adapters.qiskit
    pip install "sigma-c-framework[pennylane]"   # sigma_c.adapters.pennylane
    pip install "sigma-c-framework[braket]"      # sigma_c.adapters.braket

Import the one you want explicitly, e.g. `from sigma_c.adapters.qiskit import
observable_sweep`. Nothing here is blind: each adapter has a known-answer test
(a swept circuit whose sigma_c is known analytically), run against the SDK's own
local simulator.
"""

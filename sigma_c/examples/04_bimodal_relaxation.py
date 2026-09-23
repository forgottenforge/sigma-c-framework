# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Example 4 -- Bimodal relaxation: regime II in the wild.

A two-state system relaxing via two channels with rates 1/tau_fast and
1/tau_slow:
   O(t) = a * exp(-t / tau_fast) + (1 - a) * exp(-t / tau_slow)

When the two scales are well-separated (tau_slow / tau_fast > 5), chi_O
shows two clear peaks -- regime II. The kernel returns sigma_c as a list of peak
locations, not a single number.

Aha effect:
  A system with two characteristic times has two sigma_c values, not one
  averaged value. The kernel surfaces this honestly instead of silently picking
  the bigger peak.

Real-world analogues:
  - Glass relaxation (alpha and beta processes)
  - Quantum decoherence with two noise channels
  - Server latency: TCP retry + DB query in series
"""
from pathlib import Path
import numpy as np

from sigma_c import analyze, bare


HERE = Path(__file__).resolve().parent


def main() -> None:
    print("=" * 60)
    print("EXAMPLE 4 -- Bimodal relaxation")
    print("  Two relaxation channels -> two sigma_c.")
    print("=" * 60)

    tau_fast = 0.5    # fast channel
    tau_slow = 25.0   # slow channel (50x slower)
    a = 0.6           # weight of the fast channel

    t = np.geomspace(0.02, 500.0, 800)
    O = a * np.exp(-t / tau_fast) + (1 - a) * np.exp(-t / tau_slow)

    result = analyze(
        t, O,
        window=bare(),
        label="Bimodal relaxation (two channels)",
    )
    result.x_name = "time t"
    result.y_name = "chi_O(t)"

    print(result.summary())
    print()
    print(f"Theoretical tau_fast : {tau_fast:.3f}")
    print(f"Theoretical tau_slow : {tau_slow:.3f}")
    if isinstance(result.sigma_c, list):
        print(f"Recovered sigma_c    : {sorted(result.sigma_c)}")
        print("Both channels are resolved -> regime II is the honest answer.")
    elif result.sigma_c is None:
        print("No peak detected (unexpected for this data).")
    else:
        print(f"Single sigma_c = {result.sigma_c:.3f} -- channels collapsed.")
        print("Try increasing the time range or reducing min_prominence_ratio.")

    card_path = HERE / "card_04_bimodal.png"
    from sigma_c.examples import save_card
    save_card(result, card_path)
    print(f"\n--> Card saved: {card_path}")


if __name__ == "__main__":
    main()

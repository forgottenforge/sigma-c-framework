# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Example 1 -- Coffee cooling: when does a hot drink cool by half?

Newton's law of cooling: T(t) = T_room + (T_0 - T_room) * exp(-t / tau_cool).
The temperature difference O(t) := T(t) - T_room is a clean
exponential decay, sigma_c picks out the characteristic cooling time tau_cool
directly: sigma_c = tau_cool (bare window, rho_star = 1).

Aha effect:
  "Where does cooling stop being fast and become slow?" -- exactly at tau_cool.
"""
from pathlib import Path
import math
import numpy as np

from sigma_c import analyze, bare


HERE = Path(__file__).resolve().parent


def main() -> None:
    print("=" * 60)
    print("EXAMPLE 1 -- Coffee cooling")
    print("  When does your coffee cool by half?")
    print("=" * 60)

    # A typical paper-cup coffee in a 21 C room: cools with tau ~ 15 minutes.
    tau_cool = 15.0  # minutes
    t = np.geomspace(1.0, 240.0, 400)  # 1 min .. 4 hours, log-spaced
    O = np.exp(-t / tau_cool)  # normalised temperature difference

    result = analyze(
        t, O,
        window=bare(),                  # no smoothing applied to data
        T_star=1.0,                     # tau is in minutes already
        label="Coffee cooling",
    )
    result.x_name = "time t (minutes)"
    result.y_name = "chi_O(t)  =  |t * dO/dt|"
    print(result.summary())

    # Sanity-check: sigma_c should equal tau_cool to within numerical precision
    print()
    print(f"Theoretical tau_cool : {tau_cool:.4f} min")
    print(f"Recovered sigma_c    : {result.sigma_c:.4f} min")
    print(f"Relative error       : {abs(result.sigma_c - tau_cool) / tau_cool:.2%}")

    card_path = HERE / "card_01_coffee_cooling.png"
    from sigma_c.examples import save_card
    save_card(result, card_path)
    print(f"\n--> Card saved: {card_path}")


if __name__ == "__main__":
    main()

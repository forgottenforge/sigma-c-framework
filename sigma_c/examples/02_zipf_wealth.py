# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Example 2 -- Zipf wealth distribution: regime III as positive output.

The cumulative wealth share W(r) of the bottom-r fraction follows a
Pareto / Zipf law: W(r) ~ r^alpha for some 0 < alpha < 1. This is the
textbook scale-invariant case (paper §C.3 / Thm 9.1).

What the kernel reports: sigma_c = None (bottom value), regime III. The
non-existence is the positive answer -- wealth concentration has no
characteristic scale.

Aha effect:
  Most tools return a number even when there isn't one. The kernel returns
  "no scale exists here" as a typed, falsifiable verdict.
"""
from pathlib import Path
import numpy as np

from sigma_c import analyze


HERE = Path(__file__).resolve().parent


def main() -> None:
    print("=" * 60)
    print("EXAMPLE 2 -- Zipf wealth distribution")
    print("  When the answer is 'no scale exists'.")
    print("=" * 60)

    # Pareto-style cumulative share -- pure power law in fraction r.
    alpha = 0.6  # heavy concentration of wealth at the top
    r = np.geomspace(1e-3, 1.0, 400)
    W = r ** alpha  # share of total wealth held by bottom-r fraction

    result = analyze(
        r, W,
        window="bare",
        label="Zipf wealth concentration",
    )
    result.x_name = "wealth-rank fraction r"
    result.y_name = "chi_W(r)"
    print(result.summary())

    print()
    if result.sigma_c is None:
        print("The kernel correctly reports sigma_c = bottom (no interior scale).")
        print("This is regime III -- the positive answer for scale-invariant data.")
    else:
        print("Unexpected: a scale was detected. Check input data.")

    card_path = HERE / "card_02_zipf_wealth.png"
    from sigma_c.examples import save_card
    save_card(result, card_path)
    print(f"\n--> Card saved: {card_path}")


if __name__ == "__main__":
    main()

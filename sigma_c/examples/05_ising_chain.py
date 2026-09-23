# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Example 5 -- 1D Ising chain: the paper's textbook anchor.

The classical 1D Ising chain at zero field with inverse temperature beta
and coupling J has spin-spin correlation length
   tau = -1 / log tanh(beta * J)
This is the textbook number, derived from the 2x2 transfer matrix.

sigma_c on the bare correlator C(r) = (tanh(beta * J))^r recovers tau exactly
(bare window: rho_star = 1). See paper Example C.7.

Aha effect:
  The framework reproduces a 90-year-old stat-mech result without ever
  diagonalising the transfer matrix -- only by reading the susceptibility
  of the correlator.
"""
from pathlib import Path
import math
import numpy as np

from sigma_c import analyze, bare, Framework


HERE = Path(__file__).resolve().parent


def main() -> None:
    print("=" * 60)
    print("EXAMPLE 5 -- 1D Ising chain (paper anchor C.7)")
    print("  Recover the correlation length from chi_O alone.")
    print("=" * 60)

    beta_J = 0.5
    t = math.tanh(beta_J)
    tau_textbook = -1.0 / math.log(t)
    print(f"Inverse temperature beta * J : {beta_J}")
    print(f"Textbook tau                 : {tau_textbook:.6f}")

    # Bare correlator C(r) = t^r, on a log-spaced lag grid
    sigma = np.geomspace(0.05, 50.0, 600)
    C = t ** sigma

    result = analyze(
        sigma, C,
        window=bare(),
        framework=Framework.TRANSFER_MATRIX_1D,
        label="1D Ising bare correlator (beta*J = 0.5)",
    )
    result.x_name = "lag sigma (lattice units)"
    result.y_name = "chi_C(sigma)"
    print(result.summary())

    print()
    print(f"Theoretical tau     : {tau_textbook:.6f}")
    print(f"Recovered sigma_c   : {result.sigma_c:.6f}")
    err = abs(result.sigma_c - tau_textbook) / tau_textbook
    print(f"Relative error      : {err:.3%}")
    print()
    if err < 0.01:
        print("OK -- reproduces the textbook correlation length.")
    else:
        print("Discrepancy -- check grid resolution.")

    card_path = HERE / "card_05_ising.png"
    from sigma_c.examples import save_card
    save_card(result, card_path)
    print(f"\n--> Card saved: {card_path}")


if __name__ == "__main__":
    main()

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Example 3 (Stufe 4) -- the rho_star principle, and an HONEST two-probe refusal.

A single underlying system (exponential correlator C(r) = exp(-r/tau)) observed
through two different analytical windows:
   bare   : O_1(r) = C(r)                          -> rho_star = 1
   gamma2 : O_2(r) = sigma * tau^2/(sigma+tau)^2   -> rho_star = 2 - sqrt(3)

The rho_star principle: sigma_c[O]/rho_star = tau, so probe-dependent sigma_c is a
CONVENTION factor, not a measurement error -- divide it out and the systems agree.

Two honest wrinkles this demo shows on purpose (it does NOT hide them):
  * The gamma2-windowed observable O_2 is TWO-LOBED, so the tool reports it as
    regime II (two peaks). There is no single sigma_c for it until YOU declare which
    lobe you mean.
  * Therefore `two_probe_test(r_1, r_2)` correctly DECLINES (cause: regime_ii) --
    it refuses to merge two lobes into one tau. That refusal is the tool being
    honest, not a bug. Only after we DECLARE the rising-side (O_*^+) convention does
    a single tau for probe 2 exist, and then it matches probe 1's tau.
"""
from pathlib import Path
import math
import numpy as np

from sigma_c import analyze, two_probe_test, bare, gamma_k


HERE = Path(__file__).resolve().parent


def main() -> None:
    print("=" * 60)
    print("EXAMPLE 3 -- Two windows on one system")
    print("  Different sigma_c values, same tau.")
    print("=" * 60)

    tau_true = 5.0

    # --- Probe 1: bare correlator ---
    sigma = np.geomspace(0.05, 100.0, 600)
    O_1 = np.exp(-sigma / tau_true)
    r_1 = analyze(sigma, O_1, window=bare(),
                  label="Probe 1: bare correlator")
    r_1.x_name = "lag sigma"
    r_1.y_name = "chi_{O_1}(sigma)"

    # --- Probe 2: Gamma-2 windowed correlator ---
    # The analytic windowed observable is sigma * tau^2 / (sigma + tau)^2,
    # which the gamma2 window's known Mellin transform produces.
    O_2 = sigma * tau_true ** 2 / (sigma + tau_true) ** 2
    r_2 = analyze(sigma, O_2, window=gamma_k(2),
                  label="Probe 2: Gamma-2 windowed correlator")
    r_2.x_name = "lag sigma"
    r_2.y_name = "chi_{O_2}(sigma)"

    print("\n--- Probe 1 (bare) ---")
    print(r_1.summary())
    print("\n--- Probe 2 (Gamma-2) ---")
    print(r_2.summary())

    # Non-circular two-probe test. Expect it to DECLINE here (probe 2 is regime II).
    test = two_probe_test(r_1, r_2, delta_threshold=0.05)
    print("\n--- Non-circular two-probe test ---")
    print(test.summary())
    if not test.passed:
        print("  ^ CORRECT refusal: probe 2 is two-lobed (regime II), so there is no")
        print("    single tau to agree on -- the test will not merge two lobes into one.")

    print(f"\nTheoretical tau                 : {tau_true:.4f}")
    # tau_bridge = sigma_c/rho_star is the load-dominant read-out (NOT a certified
    # relaxation time). Probe 1 is single-mode, so it has a single bridge value.
    print(f"probe 1 (bare) tau_bridge       : {r_1.tau_bridge:.4f}")
    # Probe 2 is two-lobed: no single sigma_c until WE declare which lobe. We declare
    # the rising-side (O_*^+) convention EXPLICITLY -- this is a human choice, printed
    # so it is never hidden -- and only then does probe 2 yield a tau, matching probe 1.
    if isinstance(r_2.sigma_c, list):
        rising = min(r_2.sigma_c)
        print(f"probe 2 (gamma2) is two-lobed   : sigma_c = "
              f"[{', '.join(f'{v:.3f}' for v in r_2.sigma_c)}]")
        print(f"  -> DECLARE rising-side (O_*^+)  : sigma_c = {rising:.4f}")
        print(f"  -> tau = sigma_c / rho_star     : {rising / r_2.rho_star:.4f}  "
              f"(matches probe 1 after the convention is declared)")
    elif r_2.tau_bridge is not None:
        print(f"probe 2 (gamma2) tau_bridge     : {r_2.tau_bridge:.4f}")

    card_1 = HERE / "card_03a_probe_bare.png"
    card_2 = HERE / "card_03b_probe_gamma2.png"
    from sigma_c.examples import save_card
    save_card(r_1, card_1)
    save_card(r_2, card_2)
    print(f"\n--> Cards saved: {card_1.name}, {card_2.name}")


if __name__ == "__main__":
    main()

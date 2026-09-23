# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Demo: the certified relaxation time tau -- and the honest refusal.

Raw shot. Two operators, same question ("what is the intrinsic relaxation time?"):

  1. a REVERSIBLE Markov chain -- self-adjoint in L^2(pi) (detailed balance), NOT
     symmetric in euclidean l^2. The instrument VERIFIES self-adjointness (derives
     pi, checks Kolmogorov), computes the certified spectral abscissa
     tau = -T*/log(lambda_2), and -- with a faithful probe and a long enough window
     -- says the rate is READABLE.

  2. a DRIFT walk -- a circulant chain. It is euclidean-NORMAL but NOT reversible
     (Kolmogorov fails around the cycle). The live-route theorems need
     self-adjointness, not mere normality, so the instrument REFUSES:
     NOT_APPLICABLE. This is the point of the whole edition -- it declines the
     case it cannot certify instead of returning a number that looks fine.

Known answer: for the reversible chain the certified tau equals -T*/log(lambda_2)
of its own transfer operator; for the drift walk both tau verdicts are
NOT_APPLICABLE. From here it is your turn: bring your own transfer operator (and a
probe for the window read-out), declare its inner product, and read tau_abscissa +
window_readability off the Result.

Literature anchor: reversible-chain relaxation / spectral gap -- Levin, Peres,
Wilmer, *Markov Chains and Mixing Times* (2nd ed., 2017), ch. 12.
"""
import math

import numpy as np

from sigma_c import analyze, bare, Framework


# a reversible birth-death chain: pi = [1,2,1]/4, self-adjoint in L^2(pi),
# NOT symmetric as a euclidean matrix (P[0,1]=0.5 != P[1,0]=0.25).
P_REVERSIBLE = np.array([[0.5, 0.5, 0.0],
                         [0.25, 0.5, 0.25],
                         [0.0, 0.5, 0.5]])

# a drift walk: circulant (euclidean-normal) but NOT reversible.
P_DRIFT = np.array([[0.0, 0.7, 0.3],
                    [0.3, 0.0, 0.7],
                    [0.7, 0.3, 0.0]])

T_STAR = 1.0


def _observable(tau_slow, sigma_max=60.0, n=1000):
    """A single-mode relaxation observable whose slow scale is tau_slow."""
    sigma = np.geomspace(0.02, sigma_max, n)
    return sigma, np.exp(-sigma / tau_slow)


def run() -> dict:
    """Compute both cases; return the known-answer facts (for the suite)."""
    # analytic lambda_2 of the reversible chain -> the tau it must certify
    lam = np.sort(np.linalg.eigvals(P_REVERSIBLE).real)[::-1]
    lam2 = lam[1]
    tau_expected = -T_STAR / math.log(lam2)

    sigma, O = _observable(tau_expected)
    r_rev = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                    operator=P_REVERSIBLE, inner_product="auto",
                    gamma_A=0.5, T_star=T_STAR, sigma_axis="evolution_time")
    r_drift = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                      operator=P_DRIFT, inner_product="auto",
                      gamma_A=0.5, T_star=T_STAR, sigma_axis="evolution_time")
    return {
        "lambda_2": lam2,
        "tau_expected": tau_expected,
        "reversible": {
            "selfadjoint": r_rev.selfadjoint,
            "tau_abscissa": r_rev.tau_abscissa,
            "abscissa_code": r_rev.tau_abscissa_status.code.value,
            "window_code": r_rev.window_readability_status.code.value,
        },
        "drift": {
            "selfadjoint": r_drift.selfadjoint,
            "abscissa_code": r_drift.tau_abscissa_status.code.value,
            "window_code": r_drift.window_readability_status.code.value,
            "reason": r_drift.tau_abscissa_status.reason,
        },
    }


def main():
    d = run()
    print("=" * 64)
    print("Certified relaxation time tau -- and the honest refusal")
    print("=" * 64)
    print(f"reversible chain  lambda_2 = {d['lambda_2']:.4f} -> "
          f"tau = -T*/log(lambda_2) = {d['tau_expected']:.4f} (expected)")
    rev = d["reversible"]
    print(f"  selfadjoint            : {rev['selfadjoint']}  (verified in L^2(pi))")
    print(f"  tau_abscissa           : {rev['tau_abscissa']:.4f}  [{rev['abscissa_code']}]")
    print(f"  window readable?       : {rev['window_code']}")
    print()
    dr = d["drift"]
    print(f"drift walk (circulant, euclidean-normal but NOT reversible)")
    print(f"  selfadjoint            : {dr['selfadjoint']}")
    print(f"  tau_abscissa           : [{dr['abscissa_code']}]  <- the honest refusal")
    print(f"  why                    : {dr['reason'][:88]}...")
    print()
    print("The instrument certifies tau for the reversible operator and REFUSES the")
    print("drift walk. From here it's your turn: supply your own transfer operator +")
    print("probe, declare its inner product, read tau_abscissa + window_readability.")

    # Save the one example card where tau_abscissa is CERTIFIED (window read? OK).
    # The probe's own dominant scale (3.0) is deliberately DIFFERENT from the
    # operator's certified relaxation time (1.443), so the card shows the contrast:
    # tau_bridge = sigma_c/rho_star = 3 (the load-dominant guess, red) sitting next
    # to the certified tau_abscissa = 1.443 (the number you may actually quote).
    from pathlib import Path
    sigma, O = _observable(3.0)
    r_rev = analyze(sigma, O, window=bare(), framework=Framework.REVERSIBLE_MARKOV,
                    operator=P_REVERSIBLE, inner_product="auto",
                    gamma_A=0.5, T_star=T_STAR, sigma_axis="evolution_time")
    r_rev.title = "Reversible chain — CERTIFIED tau"
    r_rev.x_name = "σ (evolution time)"
    card = Path(__file__).resolve().parent / "card_06_reversible_tau.png"
    from sigma_c.examples import save_card
    save_card(r_rev, card)
    print(f"\n--> Card saved: {card.name}  (tau_abscissa CERTIFIED, window OK)")


if __name__ == "__main__":
    main()

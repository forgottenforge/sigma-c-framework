# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
LADDER STUFE 5 -- bring your own data (a copy-and-edit skeleton).

Replace the two lines marked  # <<< REPLACE  with your own dial `sigma` and your
own measurement `O`, then run. Everything else -- the verdict, the resolution band,
the plain-language read-out -- the tool fills in. See ../BYO_GUIDE.md for the five
decisions behind the choices here.

As shipped it runs on a stand-in dataset so you can see the shape of the output
before you edit; the two REPLACE lines are the only ones you touch.
"""
import numpy as np

from sigma_c import analyze, gamma_k


def load_your_data():
    """Return (sigma, O): your positive, ascending dial and your measurement.

    The stand-in below is a single smooth fade (so this file runs out of the box).
    Delete it and return your own two arrays.
    """
    sigma = np.geomspace(0.1, 50.0, 300)          # <<< REPLACE with your dial
    O = np.exp(-sigma / 7.0)                       # <<< REPLACE with your measurement
    return sigma, O


def run() -> dict:
    sigma, O = load_your_data()
    # Decision 1: window. bare() if O is clean; gamma_k(2) to smooth noisy data.
    result = analyze(sigma, O, window=gamma_k(2))
    return {
        "sigma_c": result.sigma_c,
        "status": result.sigma_c_status.code.value,
        "resolution_band": result.resolution_band(),
        "summary": result.summary(),
    }


def main():
    d = run()
    print(d["summary"])
    print()
    print(">>> Your turn: replace the two marked lines in load_your_data() with your")
    print(">>> own (sigma, O), then trust sigma_c only when the verdict reads OK.")


if __name__ == "__main__":
    main()

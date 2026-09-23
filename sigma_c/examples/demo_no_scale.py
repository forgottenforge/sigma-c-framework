# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
LADDER STUFE 3 -- the tool says "no", and that is a good answer.

Plain words (no jargon):
  Not every curve has a special point. Here is a curve that just keeps rising
  smoothly -- no bump, no bend, nothing that stands out at one place. We ask the tool
  "where is the special scale?" and the honest answer is: there isn't one. The tool
  says so (NOT_RESOLVABLE) instead of inventing a number to please you.

That is the whole lesson of Stufe 3: a tool you can trust is one that refuses when
there is nothing to find. A tool that always returns a confident number is the one
to distrust. (Turn a dial, measure something that only ever goes up -> no tipping
point exists, and the tool declines to name one.)
"""
import numpy as np

from sigma_c import analyze, bare


def run() -> dict:
    sigma = np.geomspace(0.1, 50.0, 400)
    O = np.sqrt(sigma)                    # a smooth, ever-rising curve: no bump anywhere
    r = analyze(sigma, O, window=bare())
    return {
        "sigma_c": r.sigma_c,
        "status": r.sigma_c_status.code.value,
        "reason": r.sigma_c_status.reason,
    }


def main():
    d = run()
    print("=" * 60)
    print("Stufe 3 -- when the honest answer is 'there is no scale here'")
    print("=" * 60)
    print(f"a curve that only ever rises (no bump)")
    print(f"where is the special scale?  : {d['sigma_c']}   <- None, on purpose")
    print(f"the tool's verdict           : {d['status']}")
    print(f"why                          : {d['reason'][:80]}...")
    print()
    print("A tool you can trust refuses when there is nothing to find. That is why the")
    print("four verdicts exist: 'no' is a real answer, not a failure.")


if __name__ == "__main__":
    main()

# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
LADDER STUFE 1 -- hello world.

Plain words (no jargon):
  We make a curve that fades away smoothly, like a hot drink cooling. The left-right
  axis (we call it the "dial") is a scale -- think of it as time. How fast the curve
  fades is set by ONE number -- here we build it with the number 4. The tool finds the
  spot where the curve drops most steeply (its "most sensitive" point) and hands that
  number back: 4.00. It also gives a little band. That's the whole lesson: we hid a
  known number in a curve, and the tool found it again.

Why it must come back as 4: for this curve, O = exp(-dial/4), the steepest / most
sensitive spot sits exactly at 4 -- a closed form, no fitting.

What the band means: it is NOT a +/- error bar from noise. It is the finest step the
tool's ruler can tell apart here -- how fine the dial's own spacing is -- and therefore
how many digits you may honestly write down. Think "ruler markings", not "measurement
scatter".

From here it's your turn: feed your own curve, read the number with its band, and
check the verdict before you trust it. (This is a closed-form self-check; no data.)
"""
import numpy as np

from sigma_c import analyze, bare


TAU_TRUE = 4.0


def run() -> dict:
    sigma = np.geomspace(0.05, 60.0, 1200)
    O = np.exp(-sigma / TAU_TRUE)
    r = analyze(sigma, O, window=bare())
    band = r.resolution_band()
    return {
        "tau_true": TAU_TRUE,
        "sigma_c": r.sigma_c,
        "sigma_c_status": r.sigma_c_status.code.value,
        "band_lo": band["band_lo"],
        "band_hi": band["band_hi"],
        "sigma_c_display": band["sigma_c_display"],
        "significant_digits": band["significant_digits"],
    }


def main():
    d = run()
    print("=" * 60)
    print("Hello world: hide a known number in a curve, find it again")
    print("=" * 60)
    print(f"we built the curve with : {d['tau_true']}")
    print(f"the tool found          : {d['sigma_c_display']}  "
          f"[{d['sigma_c_status']}]  ({d['significant_digits']} sig. digits)")
    print(f"how tightly pinned (band): [{d['band_lo']:.4g}, {d['band_hi']:.4g}]  (NOT a +/- error bar)")
    hit = d["band_lo"] <= d["tau_true"] <= d["band_hi"]
    print(f"the number we hid is     : {'inside the band -- found it' if hit else 'OUTSIDE the band'}")
    print()
    print("From here it's your turn: feed your own curve, read the number with its")
    print("band, and check the verdict before you trust it.")


if __name__ == "__main__":
    main()

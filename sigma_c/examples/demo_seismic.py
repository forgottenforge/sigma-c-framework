# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
LADDER STUFE 3 (real-world companion) -- a bump is not a discovery.

Plain words (no jargon):
  Earthquakes follow a famous law: small ones are common, big ones rare, on a steady
  slope (Gutenberg-Richter; the slope is called "b", and b is about 1). But your
  seismometer misses the smallest quakes -- below some magnitude Mc it just doesn't
  catch them all. So a real catalogue rises to a peak near Mc and then falls off.

  Point the scale-finder at that catalogue and it happily reports a "special scale"
  right at Mc. It looks like a discovery. It is NOT: it is the magnitude where your
  instrument goes deaf -- a fact about your detector, not about the earth. The real
  physics (the slope b) comes from a different, standard calculation above Mc.

The lesson: the tool finding a bump does not mean you found something real. Always ask
"is this bump the world, or is it my measurement's edge?" Here it is the edge.

We use a synthetic catalogue with a KNOWN answer (b = 1.0, Mc = 2.0) so you can check
every step. The real thing is your turn: pull a USGS/ComCat catalogue and declare
your completeness magnitude before you trust any "transition".
Anchor: Gutenberg-Richter law; b-value by Aki & Utsu (1965) maximum likelihood.
"""
import numpy as np

from sigma_c import analyze, gamma_k

B_TRUE, MC_TRUE = 1.0, 2.0


def _synthetic_catalogue(rng, n=200_000):
    """GR magnitudes down to 0.5, thinned by a detection rolloff around Mc."""
    raw = 0.5 + rng.exponential(1.0 / (B_TRUE * np.log(10)), size=n)
    p_detect = 1.0 / (1.0 + np.exp(-(raw - MC_TRUE) / 0.25))
    return raw[rng.random(raw.size) < p_detect]


def run() -> dict:
    rng = np.random.default_rng(11)
    M = _synthetic_catalogue(rng)
    dM = 0.1
    edges = np.arange(0.5, 6.0, dM)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist, _ = np.histogram(M, bins=edges)

    # generator self-check: the GR slope above Mc must be ~ -b (on well-populated bins)
    ok = (hist >= 30) & (centers >= MC_TRUE + 0.2)
    gr_slope = float(np.polyfit(centers[ok], np.log10(hist[ok]), 1)[0])

    # completeness Mc via max-curvature (mode of the non-cumulative FMD + 0.2)
    mc_hat = float(centers[np.argmax(hist)] + 0.2)
    # b via Aki-Utsu MLE, taken safely ABOVE Mc (unbinned mean)
    Mb = M[M >= mc_hat + 0.3]
    b_hat = float(np.log10(np.e) / (Mb.mean() - (mc_hat + 0.3)))

    # the scale-finder on the cumulative catalogue: where does IT say the scale is?
    cum = np.array([(M >= m).sum() for m in centers], float)
    r = analyze(centers, cum, window=gamma_k(2), min_prominence_ratio=0.5)
    return {
        "b_true": B_TRUE, "mc_true": MC_TRUE,
        "gr_slope": gr_slope, "b_hat": b_hat, "mc_hat": mc_hat,
        "sigma_c": r.sigma_c, "sigma_c_status": r.sigma_c_status.code.value,
    }


def main():
    d = run()
    print("=" * 64)
    print("Stufe 3 -- a bump is not a discovery (earthquake catalogue)")
    print("=" * 64)
    print(f"generator check: GR slope above Mc = {d['gr_slope']:.2f}  (should be -b = -1.0)")
    print(f"the real physics (standard tool): b = {d['b_hat']:.2f}  (true {d['b_true']})")
    print(f"completeness magnitude Mc        : {d['mc_hat']:.2f}  (true {d['mc_true']})")
    print(f"the scale-finder's 'special scale': sigma_c = {d['sigma_c']:.2f}  [{d['sigma_c_status']}]")
    print()
    print(f"sigma_c landed on ~{d['mc_true']} -- that is Mc, where the SEISMOMETER goes deaf,")
    print("NOT a law of the earth. The bump is your instrument's edge. The real slope b")
    print("comes from the standard calculation ABOVE Mc. A bump is not a discovery.")


if __name__ == "__main__":
    main()

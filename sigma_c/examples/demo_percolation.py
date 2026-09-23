# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
LADDER STUFE 2 -- recover a textbook number: the percolation threshold.

Plain words (no jargon):
  We randomly connect neighbouring cells on a grid, each connection with
  probability p. We ask: how big is the largest connected blob? As we turn p up,
  there is a sudden point where one blob spans the whole grid -- like asking "how
  wet does ground coffee have to be before water finds a path all the way through?"
  The tool finds that tipping point. For this grid the exact known answer is p = 1/2.

The tool sweeps p, watches the largest-blob fraction, and marks where it is most
sensitive to p. On a 32x32 grid it lands at ~0.51; on bigger grids it marches onto
the exact 1/2 (0.53 at 16, 0.51 at 32, 0.50 at 64). The resolution band says how
tightly the location is pinned. This is a *verification*: the instrument recovers a
value someone else proved exactly (bond percolation on the square lattice, p_c = 1/2;
Kesten 1980).
"""
import numpy as np

from sigma_c import analyze, gamma_k

P_C_EXACT = 0.5


def _largest_blob_fraction(L, p, rng):
    """Fraction of cells in the largest cluster of open bonds on an LxL grid."""
    n = L * L
    parent = np.arange(n)

    def find(x):
        r = x
        while parent[r] != r:
            r = parent[r]
        while parent[x] != r:
            parent[x], x = r, parent[x]
        return r

    for i in range(L):
        for j in range(L):
            s = i * L + j
            if j + 1 < L and rng.random() < p:
                a, b = find(s), find(i * L + j + 1)
                if a != b:
                    parent[a] = b
            if i + 1 < L and rng.random() < p:
                a, b = find(s), find((i + 1) * L + j)
                if a != b:
                    parent[a] = b
    _, counts = np.unique([find(x) for x in range(n)], return_counts=True)
    return counts.max() / n


def run(L=32, seeds=30) -> dict:
    ps = np.linspace(0.30, 0.70, 33)
    O = np.array([
        np.mean([_largest_blob_fraction(L, p, np.random.default_rng(9000 * L + 13 * k))
                 for k in range(seeds)])
        for p in ps
    ])
    r = analyze(ps, O, window=gamma_k(2), min_prominence_ratio=0.5)
    band = r.resolution_band()
    return {
        "p_c_exact": P_C_EXACT,
        "sigma_c": r.sigma_c,
        "status": r.sigma_c_status.code.value,
        "band_lo": band["band_lo"] if band else None,
        "band_hi": band["band_hi"] if band else None,
        "L": L,
    }


def main():
    d = run()
    print("=" * 60)
    print("Stufe 2 -- find the percolation threshold (exact answer 1/2)")
    print("=" * 60)
    print(f"grid                     : {d['L']}x{d['L']}")
    print(f"threshold found (sigma_c): {d['sigma_c']:.4f}  [{d['status']}]")
    print(f"resolution band          : [{d['band_lo']:.3f}, {d['band_hi']:.3f}]")
    print(f"exact known answer        : p_c = {d['p_c_exact']}  (Kesten 1980)")
    print(f"off by                   : {abs(d['sigma_c'] - d['p_c_exact']):.3f} "
          f"(finite-grid; marches onto 1/2 as the grid grows)")
    print()
    print("The tool recovered a value someone proved exactly. That is the point of")
    print("Stufe 2: before trusting it on the unknown, watch it hit the known.")


if __name__ == "__main__":
    main()

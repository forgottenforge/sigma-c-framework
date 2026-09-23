# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Sigma-c Card — publication-quality visualization, light theme.

Sigma-c-specific semantics (not a tipping-point / score card):
  - Regime I  (single mode)        -> blue  (clean, falsifiable reading)
  - Regime II (multi-mode)         -> amber (vector-valued sigma_c)
  - Regime III (no peak / floor)   -> grey  (sigma_c = bottom, a positive output)
  - rho_star fitted (not analytic) -> red flag

Layout is a fixed landscape card (never `bbox_inches="tight"`, which would let a
long footnote stretch the figure into a strip): a title, the chi_O panel, a
collision-free metrics column, a wrapped notes strip, and a dark branded footer.
"""
from __future__ import annotations
import textwrap
from pathlib import Path
from typing import Optional, Union

import numpy as np
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
except ImportError as _exc:  # pragma: no cover - exercised via the friendly message
    raise ImportError(
        "result.card() needs matplotlib, which is an opt-in extra (the kernel and "
        "analyze() run headless without it). Install it with:\n"
        '    pip install "sigma-c-framework[plot]"\n'
        "Everything except the PNG card works without matplotlib."
    ) from _exc

from sigma_c.result import Result


# ---------------------------------------------------------------------------
# Palette — light theme
# ---------------------------------------------------------------------------
C_BG = "#ffffff"
C_CARD = "#f8f9fb"
C_TEXT = "#1a1a2e"
C_DIM = "#6b7280"
C_LIGHT = "#d1d5db"
C_GRID = "#e5e7eb"
C_BRAND = "#17171c"        # footer bar

C_REGIME_I = "#2563eb"     # blue  — single mode, falsifiable
C_REGIME_II = "#d97706"    # amber — multi-mode
C_REGIME_III = "#6b7280"   # grey  — undefined / scale-invariant
C_FLOOR = "#7c3aed"        # purple — operational floor
C_FITTED = "#dc2626"       # red   — fitted rho_star / not certified
C_OK = "#059669"           # green — analytic rho_star / certified

_PKG_DIR = Path(__file__).resolve().parent


def _regime_color(result: Result) -> str:
    if result.regime.operational_floor_triggered:
        return C_FLOOR
    if result.sigma_c is None:
        return C_REGIME_III
    if isinstance(result.sigma_c, list):
        return C_REGIME_II
    return C_REGIME_I


def _regime_text(result: Result) -> str:
    if result.regime.operational_floor_triggered:
        return "REGIME III  ·  operational floor"
    g = result.regime.geometric
    return {"I_geom": "REGIME I  ·  single mode",
            "II_geom": "REGIME II  ·  multi-mode"}.get(g, "REGIME III  ·  no interior scale")


def _short_cite(label: str) -> str:
    """A compact, honest citation for the card: 'Thm 3.7 (thm:compat)' plus a
    one-word status tag when the proof is not PROVED. NEVER the full register
    prose (that is what blew the old card into a wide strip)."""
    from sigma_c.theorem_map import _parse_theorem_map, _parse_proof_status
    num = _parse_theorem_map().get(label)
    kind = {"thm": "Thm", "prop": "Prop", "def": "Def", "lem": "Lem",
            "cor": "Cor", "rem": "Rem", "sec": "§"}.get(label.split(":")[0], "")
    base = f"{kind} {num} ({label})".strip() if num else label
    entry = _parse_proof_status().get(label)
    if entry and entry[0] != "PROVED":
        base += f" [{entry[0]}]"
    return base


def _logo_rgba(name: str, height_px: int = 26):
    """Load a white-on-dark logo, key out the dark background, crop to content.
    Guarded: returns None if Pillow or the file is unavailable (text fallback)."""
    try:
        from PIL import Image
        path = _PKG_DIR / name
        if not path.exists():
            return None
        arr = np.array(Image.open(path).convert("RGBA"))
        bright = arr[:, :, :3].max(axis=2)
        arr[bright < 64, 3] = 0                       # dark bg -> transparent
        alpha = arr[:, :, 3]
        rows, cols = np.any(alpha > 0, axis=1), np.any(alpha > 0, axis=0)
        if not rows.any():
            return None
        r0, r1 = np.where(rows)[0][[0, -1]]
        c0, c1 = np.where(cols)[0][[0, -1]]
        crop = Image.fromarray(arr[r0:r1 + 1, c0:c1 + 1])
        w = max(1, int(height_px * crop.width / crop.height))
        return np.array(crop.resize((w, height_px), Image.LANCZOS))
    except Exception:
        return None


def _draw_footer(fig):
    """Dark branded footer bar: sigma_c mark (left), licence (centre), the
    ForgottenForge mark (right); text branding always, logos when available."""
    from matplotlib.offsetbox import OffsetImage, AnnotationBbox
    bar = fig.add_axes([0, 0, 1, 0.062])
    bar.set_xlim(0, 1); bar.set_ylim(0, 1); bar.axis("off")
    bar.set_facecolor(C_BRAND)
    bar.add_patch(plt.Rectangle((0, 0), 1, 1, color=C_BRAND, zorder=0))

    bar.text(0.085, 0.5, "SIGMA C v6", ha="left", va="center", zorder=3,
             fontsize=10, fontweight="bold", color="#ffffff")
    bar.text(0.5, 0.5,
             "© 2026 ForgottenForge  •  AGPL-3.0 or Commercial  •  forgottenforge.xyz",
             ha="center", va="center", fontsize=7.5, color="#9aa0aa",
             fontfamily="monospace", zorder=3)
    bar.text(0.965, 0.5, "FORGOTTENFORGE", ha="right", va="center", zorder=3,
             fontsize=10, fontweight="bold", color="#ffffff")

    for name, x, box_align in (("c.jpg", 0.04, (0, 0.5)),
                               ("forgottenforge.jpg", 0.80, (1, 0.5))):
        rgba = _logo_rgba(name, height_px=24)
        if rgba is not None:
            ab = AnnotationBbox(OffsetImage(rgba, zoom=0.5), (x, 0.5),
                                xycoords="axes fraction", frameon=False,
                                box_alignment=box_align, zorder=4)
            bar.add_artist(ab)


def render(
    result: Result,
    sigma: np.ndarray,
    chi: np.ndarray,
    *,
    title: str = "",
    x_name: str = "σ (resolution)",
    y_name: str = "χ_O(σ)",
    save_to: Optional[Union[str, Path]] = None,
    figsize=(12.0, 7.0),
):
    """Render a sigma-c card and return the Figure (saving to save_to if given)."""
    fig = plt.figure(figsize=figsize, facecolor=C_BG)
    color = _regime_color(result)

    # ---- Title ----
    fig.text(0.055, 0.945, title or "sigma_c  —  disciplined-reader output",
             fontsize=16, color=C_TEXT, weight="bold", va="top")
    fig.text(0.055, 0.90, _regime_text(result), fontsize=11, color=color,
             weight="bold", va="top")

    # ---- chi_O panel (left) ----
    ax = fig.add_axes([0.055, 0.32, 0.55, 0.50], facecolor=C_CARD)
    ax.plot(sigma, chi, color=color, lw=2.2)
    ax.fill_between(sigma, 0, chi, color=color, alpha=0.10)
    ax.set_xscale("log")
    ax.set_xlabel(x_name, color=C_TEXT, fontsize=11)
    ax.set_ylabel(y_name, color=C_TEXT, fontsize=11)
    ax.tick_params(colors=C_DIM, labelsize=9)
    ax.grid(True, alpha=0.35, color=C_GRID, ls="-", lw=0.5)
    for spine in ax.spines.values():
        spine.set_color(C_LIGHT)

    if isinstance(result.sigma_c, list):
        for s in result.sigma_c:
            ax.axvline(s, color=color, ls="--", lw=1.2, alpha=0.85)
        top = ax.get_ylim()[1]
        for s in result.sigma_c:
            ax.text(s, top * 0.94, f"{s:.3g}", color=color, ha="center", fontsize=8,
                    bbox=dict(boxstyle="round,pad=0.2", facecolor="white",
                              edgecolor=color, lw=0.8))
    elif result.sigma_c is not None:
        ax.axvline(result.sigma_c, color=color, ls="--", lw=1.5)
        ax.text(result.sigma_c, ax.get_ylim()[1] * 0.94, f"σ_c = {result.sigma_c:.3g}",
                color=color, ha="center", fontsize=10, weight="bold",
                bbox=dict(boxstyle="round,pad=0.3", facecolor="white",
                          edgecolor=color, lw=1.0))
    else:
        ax.text(0.5, 0.5, "σ_c = ⊥\n(no interior peak)", transform=ax.transAxes,
                color=C_DIM, fontsize=16, ha="center", va="center", weight="bold", alpha=0.7)

    # ---- metrics column (right) — short labels/values, no collision ----
    mx = fig.add_axes([0.655, 0.32, 0.30, 0.50], facecolor=C_BG)
    mx.axis("off")
    yv = [0.98]

    def row(label, value, vcolor=C_TEXT):
        mx.text(0.0, yv[0], label, fontsize=9, color=C_DIM, va="top", transform=mx.transAxes)
        mx.text(1.0, yv[0], value, fontsize=9.5, color=vcolor, weight="bold", ha="right",
                va="top", transform=mx.transAxes, fontfamily="monospace")
        yv[0] -= 0.082

    if result.sigma_c is None:
        row("σ_c", "⊥", C_REGIME_III)
    elif isinstance(result.sigma_c, list):
        row("σ_c", f"{len(result.sigma_c)}-vec", color)
    else:
        row("σ_c", f"{result.sigma_c:.4g}", color)

    if result.tau_bridge is not None:
        row("τ_bridge", f"{result.tau_bridge:.4g}", C_FITTED)
        mx.text(1.0, yv[0] + 0.02, "load-dominant, not certified", fontsize=7,
                color=C_FITTED, ha="right", va="top", style="italic", transform=mx.transAxes)
        yv[0] -= 0.055
    if result.tau_abscissa is not None:
        aa, wr = result.tau_abscissa_status, result.window_readability_status
        row("τ_abscissa", f"{result.tau_abscissa:.4g}", C_TEXT if aa.ok else C_FITTED)
        row("  window read?", wr.code.value, C_OK if wr.ok else C_DIM)
    elif result.tau_bridge is None:
        row("τ", "—", C_DIM)

    if result.rho_star is not None:
        rc = C_OK if result.falsifiable else C_FITTED
        row("ρ_⋆", f"{result.rho_star:.4g}", rc)
        row("source", str(result.rho_star_source), rc)
    if result.gamma_O is not None:
        row("γ_O (SOC)", f"{result.gamma_O:.3g}", C_FITTED if result.gamma_O < 0.1 else C_OK)
    if result.framework is not None:
        row("framework", result.framework.value, C_TEXT)

    yv[0] -= 0.02
    if result.falsifiable and result.sigma_c is not None:
        badge, bc = "✓ FALSIFIABLE (analytic ρ_⋆)", C_OK
    elif result.sigma_c is not None:
        badge, bc = "⚠ FITTED — exploratory", C_FITTED
    else:
        badge, bc = None, None
    if badge:
        mx.text(0.5, yv[0], badge, transform=mx.transAxes, color=bc, fontsize=9,
                weight="bold", ha="center", va="top",
                bbox=dict(boxstyle="round,pad=0.35", facecolor="white", edgecolor=bc, lw=1.1))

    # ---- notes strip (wrapped, capped) ----
    nx = fig.add_axes([0.055, 0.075, 0.90, 0.145])
    nx.axis("off")
    lines = []
    for note in (result.notes or [])[:2]:
        lines.extend(textwrap.wrap(f"• {note}", width=140)[:2])
    if result.citations:
        cites = "  ·  ".join(_short_cite(c) for c in result.citations[:4])
        lines.extend(textwrap.wrap(f"backed by:  {cites}", width=140)[:2])
    nx.text(0.0, 1.0, "\n".join(lines[:5]), transform=nx.transAxes, color=C_DIM,
            fontsize=8, va="top", linespacing=1.45)

    _draw_footer(fig)

    if save_to is not None:
        fig.savefig(str(save_to), dpi=150, facecolor=C_BG)   # fixed size; NOT tight
    return fig

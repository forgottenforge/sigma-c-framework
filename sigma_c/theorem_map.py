# SPDX-License-Identifier: AGPL-3.0-or-later OR LicenseRef-ForgottenForge-Commercial
# Copyright © 2026 ForgottenForge <nfo@forgottenforge.xyz>
"""
Theorem citation helper.

Cites the paper by LaTeX label (stable under renumbering), resolves the
current number via THEOREM_MAP.md, which ships INSIDE the package
(sigma_c/THEOREM_MAP.md) so a `pip install` keeps the proof-status layer.
A repo-root copy (one level up) is honoured as a dev fallback.
"""
from __future__ import annotations
from pathlib import Path
from typing import Dict, Optional
import re
import warnings


def _resolve_map_path() -> Path:
    """Prefer the shipped in-package map; fall back to the repo-root copy (dev)."""
    in_package = Path(__file__).resolve().parent / "THEOREM_MAP.md"
    if in_package.exists():
        return in_package
    repo_root = Path(__file__).resolve().parent.parent / "THEOREM_MAP.md"
    if repo_root.exists():
        return repo_root
    # Neither present: the proof-status layer would silently vanish. Say so.
    warnings.warn(
        "THEOREM_MAP.md not found (looked in the package and one level up): "
        "citations will render as bare labels with NO paper number and NO "
        "proof-status. This usually means a broken install. Reinstall the package.",
        RuntimeWarning, stacklevel=2,
    )
    return in_package  # nonexistent; .exists() checks downstream handle it


_MAP_PATH = _resolve_map_path()
_CACHE: Optional[Dict[str, str]] = None
_STATUS_CACHE: Optional[Dict[str, str]] = None

_PROOF_STATES = (
    "PROVED", "PUBLISHED", "PREPRINT", "GAP-KNOWN", "UNDER-REPAIR", "SUPERSEDED",
    "REVIEW-CLEAN", "UNVERIFIED", "ARCHIVED",
)
# ARCHIVED = paper-only formalism, moved to the math archive (kept privately,
# not shipped with the package);
# the live code uses named operations instead of the labelled theorem. It is not
# a proof-quality verdict (not PROVED/GAP-KNOWN): it records that the label is no
# longer a live theorem the code stands on. Renders like any other non-PROVED
# status so no output silently anchors on an archived label.
# Default for a label with no status row is NO-REGISTER-ENTRY ("nobody looked"),
# NOT PROVED. Rationale: a checker that finds nothing is
# indistinguishable from one that never ran; treating an un-entered theorem as
# proved silently certifies everything nobody has examined. NO-REGISTER-ENTRY renders
# a marker so the gap is visible. (Renamed from "UNVERIFIED", which an
# outside user read as "the math is unverified"; the new marker says "unlisted, not
# disproven". "UNVERIFIED" remains a valid EXPLICIT register status.) Only PROVED
# renders nothing. PUBLISHED is a
# distinct check status (peer-reviewed elsewhere, not re-derived against our
# chain) and DOES render, so the register says which kind of check exists.
_DEFAULT_STATE = "NO-REGISTER-ENTRY"


def _parse_theorem_map() -> Dict[str, str]:
    """Parse THEOREM_MAP.md tables into {label: number} dict."""
    global _CACHE
    if _CACHE is not None:
        return _CACHE
    mapping: Dict[str, str] = {}
    if not _MAP_PATH.exists():
        _CACHE = mapping
        return mapping
    text = _MAP_PATH.read_text(encoding="utf-8")
    # Match mapping rows: | `label` | number | ... |. The number column must
    # contain a digit -- this excludes the proof-status register rows (whose
    # second column is a status keyword like PROVED / UNDER-REPAIR), so a status
    # row can never clobber a label's real number.
    pattern = re.compile(r"^\|\s*`([\w:.-]+)`\s*\|\s*([\w.]*\d[\w.]*)\s*\|", re.M)
    for match in pattern.finditer(text):
        label, number = match.group(1), match.group(2)
        mapping[label] = number
    _CACHE = mapping
    return mapping


def _parse_proof_status() -> Dict[str, tuple]:
    """Parse the proof-status register into {label: (status, short_reason)}.

    A row `| `label` | STATUS | short reason | note |` where STATUS is one of
    _PROOF_STATES and the short reason (rendered by cite for non-PROVED labels)
    is the third column. Absent labels are treated as PROVED-as-published."""
    global _STATUS_CACHE
    if _STATUS_CACHE is not None:
        return _STATUS_CACHE
    status: Dict[str, tuple] = {}
    if not _MAP_PATH.exists():
        _STATUS_CACHE = status
        return status
    text = _MAP_PATH.read_text(encoding="utf-8")
    pattern = re.compile(
        r"^\|\s*`([\w:.-]+)`\s*\|\s*"
        r"(PROVED|PUBLISHED|PREPRINT|GAP-KNOWN|UNDER-REPAIR|SUPERSEDED|REVIEW-CLEAN|UNVERIFIED|ARCHIVED)\s*\|"
        r"\s*([^|]*?)\s*\|",
        re.M,
    )
    for match in pattern.finditer(text):
        short = match.group(3).strip()
        if short in ("", "—", "-"):
            short = ""
        status[match.group(1)] = (match.group(2), short)
    _STATUS_CACHE = status
    return status


def cite(label: str, note: str = "") -> str:
    """
    Render a paper citation by label.

    Example:
        >>> cite("thm:compat")
        'paper Thm 3.7 (thm:compat)'

    If the label is missing from THEOREM_MAP.md, falls back to the label alone.
    """
    mapping = _parse_theorem_map()
    number = mapping.get(label)
    if number is None:
        base = f"paper [{label}]"
    else:
        # Heuristic: deduce kind from label prefix
        kind = {
            "thm": "Thm",
            "prop": "Prop",
            "def": "Def",
            "lem": "Lem",
            "cor": "Cor",
            "rem": "Rem",
            "obs": "Obs",
            "sec": "§",
            "app": "App",
        }.get(label.split(":")[0], "")
        base = f"paper {kind} {number} ({label})"
    # Provenance of the PROOF, not only the label: render any non-PROVED status
    # so no output silently anchors on a proof with a known gap (THEOREM_MAP
    # proof-status register). Only PROVED renders nothing; a label with NO row
    # defaults to NO-REGISTER-ENTRY ("nobody looked") and is flagged, not treated
    # as proved.
    entry = _parse_proof_status().get(label)
    if entry is None:
        base = (
            f"{base} [proof: {_DEFAULT_STATE} — not in the proof-status register "
            f"(means UNLISTED, not disproven; the register default is not PROVED)]"
        )
    elif entry[0] != "PROVED":
        state, short = entry
        # Cap the rendered reason: the register's reason column can carry a long
        # rev-history (e.g. thm:cross-obs-concentration), and user-facing output
        # must not dump it. The full text stays in the register.
        if short and len(short) > 100:
            short = short[:97].rstrip() + "…"
        base = f"{base} [proof: {state}" + (f" — {short}" if short else "") + "]"
    if note:
        return f"{base} — {note}"
    return base


def reload_map() -> None:
    """Force re-read of THEOREM_MAP.md (after paper renumbering)."""
    global _CACHE, _STATUS_CACHE
    _CACHE = None
    _STATUS_CACHE = None
    _parse_theorem_map()
    _parse_proof_status()

"""IMGT-scheme numbering for VH sequences.

Uses ANARCI when available; falls back to simple sequential numbering.
"""

from __future__ import annotations

import logging
from typing import Any

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# IMGT region boundaries (inclusive).  Positions are integers except for
# insertion codes (e.g. "111.1").  For boundary checks we use the integer
# part only.
# ---------------------------------------------------------------------------
IMGT_REGIONS: dict[str, tuple[int, int]] = {
    "FR1": (1, 26),
    "CDR1": (27, 38),
    "FR2": (39, 55),
    "CDR2": (56, 65),
    "FR3": (66, 104),
    "CDR3": (105, 117),
    "FR4": (118, 128),
}


def classify_region(imgt_position: str) -> str:
    """Return the IMGT region name for a given position string.

    Parameters
    ----------
    imgt_position:
        An IMGT position such as ``"27"``, ``"82A"``, or ``"111.1"``.

    Returns
    -------
    str
        One of ``FR1``, ``CDR1``, ``FR2``, ``CDR2``, ``FR3``, ``CDR3``,
        ``FR4``, or ``"unknown"``.
    """
    # Extract the integer portion of the position.
    num_str = ""
    for ch in imgt_position:
        if ch.isdigit():
            num_str += ch
        else:
            break
    if not num_str:
        return "unknown"
    pos_int = int(num_str)

    for region, (start, end) in IMGT_REGIONS.items():
        if start <= pos_int <= end:
            return region
    return "unknown"


def _number_with_anarci(sequence: str) -> dict[str, str] | None:
    """Attempt to number *sequence* using ANARCI with the IMGT scheme.

    Returns *None* if ANARCI is unavailable or numbering fails.
    """
    try:
        from anarci import anarci as run_anarci  # type: ignore[import-untyped]
    except ImportError:
        logger.debug("ANARCI not installed – will use fallback numbering")
        return None

    try:
        results: Any = run_anarci([("H", sequence)], scheme="imgt", output=False)
        numbering_results, _ = results
        if not numbering_results or numbering_results[0] is None:
            logger.warning("ANARCI returned no numbering for the input sequence")
            return None

        domain_numbering = numbering_results[0]
        # domain_numbering is a list of (numbering_list, chain_type) tuples.
        # We take the first domain.
        if not domain_numbering:
            return None
        positions_list, _chain = domain_numbering[0]

        numbered: dict[str, str] = {}
        for (pos_int, insertion), aa in positions_list:
            if aa == "-":
                continue
            if insertion and insertion.strip():
                key = f"{pos_int}{insertion.strip()}"
            else:
                key = str(pos_int)
            numbered[key] = aa
        return numbered

    except Exception:
        logger.warning("ANARCI numbering failed – falling back", exc_info=True)
        return None


def _fallback_numbering(sequence: str) -> dict[str, str]:
    """Assign sequential IMGT-like positions 1..N to each residue."""
    logger.info("Using sequential fallback numbering (no ANARCI)")
    return {str(i + 1): aa for i, aa in enumerate(sequence)}


def number_sequence(sequence: str) -> dict[str, str]:
    """Number a VH amino-acid sequence using the IMGT scheme.

    Tries ANARCI first; if it is not installed or fails, falls back to
    simple sequential numbering.

    Parameters
    ----------
    sequence:
        Single-letter amino-acid VH sequence.

    Returns
    -------
    dict[str, str]
        Mapping of IMGT position string → single-letter amino acid.
    """
    if not sequence:
        raise ValueError("Empty sequence provided to number_sequence")

    result = _number_with_anarci(sequence)
    if result is not None:
        return result
    return _fallback_numbering(sequence)

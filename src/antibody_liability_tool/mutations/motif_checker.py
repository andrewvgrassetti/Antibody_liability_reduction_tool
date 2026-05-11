"""Check whether a mutation introduces problematic sequence motifs."""

from __future__ import annotations

import logging
import re
from dataclasses import dataclass

logger = logging.getLogger(__name__)


@dataclass(frozen=True)
class MotifHit:
    """A detected problematic motif at a specific location."""

    motif_name: str
    """Short identifier (e.g. ``'N-glycosylation'``)."""
    pattern: str
    """The regex or literal pattern that matched."""
    start_position: int
    """0-based start index in the linear sequence context."""
    matched_text: str
    """The actual substring that matched."""


# Compiled patterns for efficiency.
_PATTERNS: list[tuple[str, re.Pattern[str]]] = [
    ("N-glycosylation", re.compile(r"N[^P][ST]")),
    ("deamidation", re.compile(r"N[GS]")),
    ("isomerization", re.compile(r"D[GS]")),
]


def check_motifs(
    sequence_context: str,
    *,
    check_free_cys: bool = True,
    check_oxidation_met: bool = True,
    surface_positions: set[int] | None = None,
) -> list[MotifHit]:
    """Scan *sequence_context* for problematic motifs.

    Parameters
    ----------
    sequence_context:
        A short amino-acid string representing the local sequence
        around the mutation site (typically ±3 residues).
    check_free_cys:
        If ``True``, flag any cysteine (``C``) as a potential
        unpaired-cysteine liability.
    check_oxidation_met:
        If ``True``, flag any methionine (``M``) as a potential
        surface-oxidation liability.  When *surface_positions* is
        provided, only Met residues at those 0-based offsets are flagged.
    surface_positions:
        Optional set of 0-based indices within *sequence_context* that
        are surface-exposed.  Used to restrict Met-oxidation checks.

    Returns
    -------
    list[MotifHit]
        All detected motif hits.
    """
    hits: list[MotifHit] = []

    for name, pattern in _PATTERNS:
        for m in pattern.finditer(sequence_context):
            hits.append(
                MotifHit(
                    motif_name=name,
                    pattern=pattern.pattern,
                    start_position=m.start(),
                    matched_text=m.group(),
                )
            )

    if check_free_cys:
        for i, aa in enumerate(sequence_context):
            if aa == "C":
                hits.append(
                    MotifHit(
                        motif_name="free_cysteine",
                        pattern="C",
                        start_position=i,
                        matched_text="C",
                    )
                )

    if check_oxidation_met:
        for i, aa in enumerate(sequence_context):
            if aa == "M":
                if surface_positions is None or i in surface_positions:
                    hits.append(
                        MotifHit(
                            motif_name="oxidation_met",
                            pattern="M",
                            start_position=i,
                            matched_text="M",
                        )
                    )

    return hits

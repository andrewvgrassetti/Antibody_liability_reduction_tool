"""Mutation candidate generation for identified liabilities."""

from __future__ import annotations

import json
import logging
from dataclasses import dataclass
from pathlib import Path
from typing import Any

from antibody_liability_tool.liabilities.detector import Liability
from antibody_liability_tool.mutations.motif_checker import check_motifs

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data"
_FREQ_FILE = _DATA_DIR / "abysis_vh_frequencies.json"

# Amino-acid property sets for liability-reducing logic.
_POLAR_AA = set("STNQDEHKR")
_NEUTRAL_AA = set("STAGQNDE")
_HYDROPHOBIC_AA = set("FILMVWY")


@dataclass
class CandidateMutation:
    """A proposed point mutation to reduce a liability."""

    position: str
    """IMGT position string."""
    original_aa: str
    """Wild-type single-letter amino acid."""
    proposed_aa: str
    """Suggested replacement amino acid."""
    human_frequency: float
    """Frequency of *proposed_aa* at this position in the Abysis human database."""
    rationale: str
    """Human-readable explanation for why this substitution is suggested."""


_freq_data: dict[str, dict[str, float]] | None = None


def _load_frequency_data() -> dict[str, dict[str, float]]:
    """Load and cache the Abysis VH frequency JSON."""
    global _freq_data  # noqa: PLW0603
    if _freq_data is not None:
        return _freq_data

    if not _FREQ_FILE.exists():
        raise FileNotFoundError(
            f"Abysis frequency data not found: {_FREQ_FILE}. "
            "Ensure data/abysis_vh_frequencies.json exists."
        )
    with _FREQ_FILE.open() as fh:
        raw = json.load(fh)
    _freq_data = raw.get("positions", raw)
    logger.debug("Loaded Abysis frequencies for %d positions", len(_freq_data))
    return _freq_data


def _build_context(
    numbered: dict[str, str],
    target_pos: str,
    proposed_aa: str,
    window: int = 3,
) -> str:
    """Build a linear sequence context around *target_pos* with the proposed substitution.

    Positions are sorted numerically and a window of ±*window* positions
    around *target_pos* is extracted.
    """

    def _sort_key(p: str) -> tuple[int, str]:
        num = ""
        suffix = ""
        for ch in p:
            if ch.isdigit():
                num += ch
            else:
                suffix += ch
        return (int(num) if num else 0, suffix)

    sorted_positions = sorted(numbered.keys(), key=_sort_key)
    try:
        idx = sorted_positions.index(target_pos)
    except ValueError:
        return proposed_aa

    start = max(0, idx - window)
    end = min(len(sorted_positions), idx + window + 1)
    context_chars: list[str] = []
    for i in range(start, end):
        p = sorted_positions[i]
        if p == target_pos:
            context_chars.append(proposed_aa)
        else:
            context_chars.append(numbered.get(p, "X"))
    return "".join(context_chars)


def _reduces_hydrophobic_liability(original: str, proposed: str) -> bool:
    """Return True if the substitution reduces hydrophobic character."""
    return original in _HYDROPHOBIC_AA and proposed in _POLAR_AA


def _reduces_charge_liability(original: str, proposed: str) -> bool:
    """Return True if the substitution neutralises positive charge."""
    return original in {"R", "K"} and proposed in _NEUTRAL_AA


def generate_mutations(
    liabilities: list[Liability],
    numbered_sequence: dict[str, str],
    config: dict[str, Any] | None = None,
) -> list[CandidateMutation]:
    """Propose human-frequent substitutions that reduce each liability.

    Parameters
    ----------
    liabilities:
        Output of :func:`~antibody_liability_tool.liabilities.detector.detect_liabilities`.
    numbered_sequence:
        IMGT position → amino acid mapping for the VH domain.
    config:
        Full tool configuration dict.

    Returns
    -------
    list[CandidateMutation]
        Candidate mutations sorted by human frequency (descending).
    """
    if config is None:
        config = {}
    mcfg = config.get("mutations", {})
    min_freq: float = mcfg.get("min_human_frequency", 0.05)
    avoid = mcfg.get("avoid_motifs", {})
    check_cys: bool = bool(avoid.get("free_cys", True))
    check_met: bool = bool(avoid.get("oxidation_met", True))

    freq_data = _load_frequency_data()
    candidates: list[CandidateMutation] = []

    for liab in liabilities:
        pos = liab.position
        original = liab.residue
        pos_freqs = freq_data.get(pos)
        if pos_freqs is None:
            logger.debug("No Abysis data for position %s – skipping", pos)
            continue

        for aa, freq in pos_freqs.items():
            if aa == "other" or aa == original:
                continue
            if freq < min_freq:
                continue

            # Check that the substitution actually reduces the liability.
            is_hydrophobic_fix = _reduces_hydrophobic_liability(original, aa)
            is_charge_fix = _reduces_charge_liability(original, aa)
            if not (is_hydrophobic_fix or is_charge_fix):
                continue

            # Check for introduced motif liabilities.
            context = _build_context(numbered_sequence, pos, aa)
            motif_hits = check_motifs(
                context,
                check_free_cys=check_cys,
                check_oxidation_met=check_met,
            )
            if motif_hits:
                motif_names = ", ".join(h.motif_name for h in motif_hits)
                logger.debug(
                    "Skipping %s%s%s – introduces motif(s): %s",
                    original, pos, aa, motif_names,
                )
                continue

            if is_hydrophobic_fix:
                rationale = (
                    f"Replace hydrophobic {original} with polar {aa} "
                    f"(human frequency {freq:.0%} at IMGT {pos})"
                )
            else:
                rationale = (
                    f"Neutralise positive charge {original} with {aa} "
                    f"(human frequency {freq:.0%} at IMGT {pos})"
                )

            candidates.append(
                CandidateMutation(
                    position=pos,
                    original_aa=original,
                    proposed_aa=aa,
                    human_frequency=freq,
                    rationale=rationale,
                )
            )

    candidates.sort(key=lambda c: -c.human_frequency)
    logger.info(
        "Generated %d candidate mutations from %d liabilities",
        len(candidates), len(liabilities),
    )
    return candidates

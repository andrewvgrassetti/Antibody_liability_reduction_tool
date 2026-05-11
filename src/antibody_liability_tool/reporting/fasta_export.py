"""FASTA export of top candidate mutant sequences.

Exports parent and top-ranked mutant sequences in FASTA format with
descriptive headers containing mutation labels, composite scores, and
key metrics.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Sequence

from antibody_liability_tool.optimization.scorer import ScoredCandidate

logger = logging.getLogger(__name__)

_LINE_WIDTH = 80


def _wrap_sequence(sequence: str, width: int = _LINE_WIDTH) -> str:
    """Wrap a sequence string to fixed-width lines."""
    return "\n".join(sequence[i : i + width] for i in range(0, len(sequence), width))


def _build_header(
    seq_id: str,
    description: str,
    score: float | None = None,
    metrics: dict[str, float] | None = None,
) -> str:
    """Build a FASTA header line with optional score and metrics.

    Parameters
    ----------
    seq_id : str
        Sequence identifier.
    description : str
        Human-readable description.
    score : float, optional
        Composite score to include.
    metrics : dict, optional
        Key metrics to append.

    Returns
    -------
    str
        FASTA header (without leading ``>``).
    """
    parts = [seq_id, description]
    if score is not None:
        parts.append(f"score={score:.4f}")
    if metrics:
        metric_str = " ".join(f"{k}={v:.4f}" for k, v in metrics.items())
        parts.append(metric_str)
    return " ".join(parts)


def export_fasta(
    parent_sequence: str,
    ranked_candidates: Sequence[ScoredCandidate],
    output_path: str | Path,
    top_n: int = 10,
    parent_id: str = "parent_VH",
) -> Path:
    """Export parent and top-ranked mutant sequences as a FASTA file.

    Parameters
    ----------
    parent_sequence : str
        The parent (wild-type) amino-acid sequence.
    ranked_candidates : Sequence[ScoredCandidate]
        Scored and ranked mutation candidates (must have ``sequence`` set).
    output_path : str or Path
        Path to write the FASTA file.
    top_n : int
        Maximum number of candidates to include (in addition to parent).
    parent_id : str
        Identifier for the parent sequence.

    Returns
    -------
    Path
        Absolute path to the generated FASTA file.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    lines: list[str] = []

    # Parent sequence
    parent_header = _build_header(parent_id, "wild-type parent sequence")
    lines.append(f">{parent_header}")
    lines.append(_wrap_sequence(parent_sequence))

    # Top N candidates
    for i, cand in enumerate(ranked_candidates[:top_n]):
        if not cand.sequence:
            logger.warning("Candidate %s has no sequence – skipping FASTA entry", cand.label)
            continue

        seq_id = f"mutant_{i + 1}_{cand.label.replace('+', '_')}"
        key_metrics: dict[str, float] = {}
        for metric_name in ("PSH", "PPC", "PNC", "oasis_humanness", "SAP_score"):
            if metric_name in cand.mutant_metrics:
                key_metrics[metric_name] = cand.mutant_metrics[metric_name]

        header = _build_header(
            seq_id=seq_id,
            description=cand.label,
            score=cand.composite_score,
            metrics=key_metrics,
        )
        lines.append(f">{header}")
        lines.append(_wrap_sequence(cand.sequence))

    content = "\n".join(lines) + "\n"
    output.write_text(content, encoding="utf-8")
    n_written = min(top_n, len(ranked_candidates)) + 1
    logger.info("FASTA export: %d sequences written to %s", n_written, output)
    return output.resolve()

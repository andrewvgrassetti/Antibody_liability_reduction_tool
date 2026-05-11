"""Surface-exposure classification for IMGT-numbered VH positions."""

from __future__ import annotations

import json
import logging
from pathlib import Path
from typing import Any

logger = logging.getLogger(__name__)

_DATA_DIR = Path(__file__).resolve().parent.parent.parent.parent / "data"
_EXPOSURE_FILE = _DATA_DIR / "imgt_surface_exposure.json"

_exposure_data: dict[str, Any] | None = None


def _load_exposure_data() -> dict[str, Any]:
    """Load and cache the exposure JSON on first call."""
    global _exposure_data  # noqa: PLW0603
    if _exposure_data is not None:
        return _exposure_data

    if not _EXPOSURE_FILE.exists():
        raise FileNotFoundError(
            f"Surface-exposure data file not found: {_EXPOSURE_FILE}. "
            "Ensure the repository data/ directory is present."
        )
    with _EXPOSURE_FILE.open() as fh:
        raw = json.load(fh)
    _exposure_data = raw.get("positions", raw)
    logger.debug("Loaded surface-exposure data for %d positions", len(_exposure_data))
    return _exposure_data


def classify_surface_exposure(
    numbered_sequence: dict[str, str],
) -> dict[str, dict[str, str]]:
    """Classify each position in a numbered VH sequence by surface exposure.

    Parameters
    ----------
    numbered_sequence:
        Mapping of IMGT position string → amino acid (as returned by
        :func:`~antibody_liability_tool.numbering.imgt.number_sequence`).

    Returns
    -------
    dict[str, dict[str, str]]
        For every position present in the numbered sequence **and** in the
        exposure data file, a dict with keys ``"exposure"`` and ``"region"``.
        Positions without exposure data are omitted.
    """
    data = _load_exposure_data()
    result: dict[str, dict[str, str]] = {}
    for pos, aa in numbered_sequence.items():
        if pos in data:
            result[pos] = {
                "exposure": data[pos]["exposure"],
                "region": data[pos]["region"],
                "residue": aa,
            }
    return result


def get_exposed_positions(
    exposure_map: dict[str, dict[str, str]],
    *,
    include_partial: bool = True,
) -> dict[str, dict[str, str]]:
    """Filter *exposure_map* to only exposed (and optionally partially-exposed) positions.

    Parameters
    ----------
    exposure_map:
        Output of :func:`classify_surface_exposure`.
    include_partial:
        If ``True`` (default), partially-exposed positions are included.
    """
    allowed = {"exposed"}
    if include_partial:
        allowed.add("partially_exposed")
    return {pos: info for pos, info in exposure_map.items() if info["exposure"] in allowed}


def classify_surface_exposure_structural(
    numbered_sequence: dict[str, str],
    *,
    sasa_threshold_exposed: float = 40.0,
    sasa_threshold_partial: float = 20.0,
) -> dict[str, dict[str, str]]:
    """Placeholder for FreeSASA-based structural surface-exposure classification.

    This function is intended to be implemented when the ``structural``
    optional dependency (FreeSASA + ABodyBuilder) is available.

    Raises
    ------
    NotImplementedError
        Always; this is a future integration point.
    """
    raise NotImplementedError(
        "Structural surface-exposure classification requires FreeSASA and "
        "ABodyBuilder. Install the 'structural' extra and implement this method."
    )

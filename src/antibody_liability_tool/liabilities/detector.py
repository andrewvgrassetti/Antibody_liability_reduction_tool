"""Liability detection for surface-exposed VH positions."""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

from antibody_liability_tool.numbering.imgt import classify_region

logger = logging.getLogger(__name__)


@dataclass(order=True)
class Liability:
    """A single detected liability at an IMGT position.

    Instances are ordered by *severity* (descending) then *position*.
    """

    # Sort key: higher severity first, then ascending position.
    _sort_key: tuple[int, str] = field(init=False, repr=False, compare=True)

    position: str
    """The IMGT position string (e.g. ``"56"``)."""
    residue: str
    """Single-letter amino acid at this position."""
    imgt_number: str
    """Same as *position*; kept for backwards compatibility with downstream modules."""
    region: str
    """IMGT region (FR1, CDR1, …)."""
    reason: str
    """Human-readable description of the liability."""
    severity: int
    """Integer severity: 1 = low, 2 = medium, 3 = high."""

    def __post_init__(self) -> None:
        self._sort_key = (-self.severity, self.position)


_SEVERITY_HIGH = 3
_SEVERITY_MEDIUM = 2
_SEVERITY_LOW = 1


def _extract_position_int(pos: str) -> int:
    """Return the integer portion of an IMGT position for distance calculations."""
    num = ""
    for ch in pos:
        if ch.isdigit():
            num += ch
        else:
            break
    return int(num) if num else 0


def _detect_hydrophobic(
    numbered_sequence: dict[str, str],
    exposure_map: dict[str, dict[str, str]],
    hydrophobic_set: set[str],
) -> list[Liability]:
    """Flag hydrophobic residues at surface-exposed positions."""
    liabilities: list[Liability] = []
    for pos, info in exposure_map.items():
        if info["exposure"] == "buried":
            continue
        aa = numbered_sequence.get(pos, "")
        if aa in hydrophobic_set:
            region = classify_region(pos)
            sev = _SEVERITY_HIGH if info["exposure"] == "exposed" else _SEVERITY_MEDIUM
            liabilities.append(
                Liability(
                    position=pos,
                    residue=aa,
                    imgt_number=pos,
                    region=region,
                    reason=f"Hydrophobic residue {aa} at surface-{info['exposure']} position",
                    severity=sev,
                )
            )
    return liabilities


def _detect_positive_charge_clusters(
    numbered_sequence: dict[str, str],
    exposure_map: dict[str, dict[str, str]],
    charge_set: set[str],
    cluster_size: int,
    max_distance: int,
) -> list[Liability]:
    """Flag clusters of positively-charged residues at surface positions."""
    # Collect surface charge positions.
    charge_positions: list[tuple[str, str, int]] = []
    for pos, info in exposure_map.items():
        if info["exposure"] == "buried":
            continue
        aa = numbered_sequence.get(pos, "")
        if aa in charge_set:
            charge_positions.append((pos, aa, _extract_position_int(pos)))

    # Sort by integer position for clustering.
    charge_positions.sort(key=lambda x: x[2])

    # Simple single-linkage clustering.
    clusters: list[list[tuple[str, str, int]]] = []
    for item in charge_positions:
        added = False
        for cluster in clusters:
            if any(abs(item[2] - member[2]) <= max_distance for member in cluster):
                cluster.append(item)
                added = True
                break
        if not added:
            clusters.append([item])

    liabilities: list[Liability] = []
    for cluster in clusters:
        if len(cluster) < cluster_size:
            continue
        positions_str = ", ".join(c[0] for c in cluster)
        for pos, aa, _ in cluster:
            region = classify_region(pos)
            liabilities.append(
                Liability(
                    position=pos,
                    residue=aa,
                    imgt_number=pos,
                    region=region,
                    reason=(
                        f"Positive-charge cluster ({aa}) with {len(cluster)} "
                        f"charged residues in proximity (positions {positions_str})"
                    ),
                    severity=_SEVERITY_MEDIUM,
                )
            )
    return liabilities


def detect_liabilities(
    numbered_sequence: dict[str, str],
    exposure_map: dict[str, dict[str, str]],
    config: dict[str, Any] | None = None,
) -> list[Liability]:
    """Detect surface-exposed liabilities in a numbered VH sequence.

    Parameters
    ----------
    numbered_sequence:
        IMGT position → amino acid mapping.
    exposure_map:
        Output of :func:`~antibody_liability_tool.surface.exposure.classify_surface_exposure`.
    config:
        Full tool configuration dict. Falls back to built-in defaults for
        any missing keys.

    Returns
    -------
    list[Liability]
        Liabilities sorted by severity (highest first), then position.
    """
    if config is None:
        config = {}
    lcfg = config.get("liabilities", {})

    hydrophobic_set: set[str] = set(
        lcfg.get("hydrophobic_residues", ["F", "I", "L", "M", "V", "W", "Y"])
    )
    charge_set: set[str] = set(lcfg.get("positive_charge_residues", ["R", "K"]))
    cluster_size: int = lcfg.get("positive_charge_cluster_size", 2)
    max_distance: int = lcfg.get("charge_cluster_distance", 5)

    liabilities: list[Liability] = []
    liabilities.extend(
        _detect_hydrophobic(numbered_sequence, exposure_map, hydrophobic_set)
    )
    liabilities.extend(
        _detect_positive_charge_clusters(
            numbered_sequence, exposure_map, charge_set, cluster_size, max_distance
        )
    )

    liabilities.sort()
    logger.info("Detected %d total liabilities", len(liabilities))
    return liabilities

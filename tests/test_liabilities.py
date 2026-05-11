"""Tests for the liability detection module."""

from __future__ import annotations

import pytest

from antibody_liability_tool.liabilities.detector import Liability, detect_liabilities


def _make_exposure(pos: str, exposure: str, region: str = "FR1") -> dict:
    return {pos: {"exposure": exposure, "region": region}}


class TestDetectLiabilities:
    """Tests for detect_liabilities()."""

    def test_detect_hydrophobic_surface_liabilities(self) -> None:
        """Hydrophobic residues at exposed positions should be flagged."""
        numbered = {"56": "W", "57": "A"}
        exposure = {
            "56": {"exposure": "exposed", "region": "CDR2"},
            "57": {"exposure": "exposed", "region": "CDR2"},
        }
        liabilities = detect_liabilities(numbered, exposure)
        hydrophobic_liabs = [l for l in liabilities if "Hydrophobic" in l.reason]
        assert len(hydrophobic_liabs) >= 1
        assert any(l.residue == "W" for l in hydrophobic_liabs)

    def test_detect_positive_charge_clusters(self) -> None:
        """Adjacent positively charged residues at surface should trigger cluster detection."""
        numbered = {"60": "R", "62": "K", "63": "R"}
        exposure = {
            "60": {"exposure": "exposed", "region": "CDR2"},
            "62": {"exposure": "exposed", "region": "CDR2"},
            "63": {"exposure": "exposed", "region": "CDR2"},
        }
        liabilities = detect_liabilities(numbered, exposure)
        charge_liabs = [l for l in liabilities if "charge cluster" in l.reason.lower()]
        assert len(charge_liabs) >= 2

    def test_no_liabilities_for_ideal_sequence(self) -> None:
        """A sequence with only non-hydrophobic, non-charged surface residues."""
        numbered = {"10": "S", "11": "T", "12": "A", "13": "G"}
        exposure = {
            "10": {"exposure": "exposed", "region": "FR1"},
            "11": {"exposure": "exposed", "region": "FR1"},
            "12": {"exposure": "exposed", "region": "FR1"},
            "13": {"exposure": "exposed", "region": "FR1"},
        }
        liabilities = detect_liabilities(numbered, exposure)
        assert len(liabilities) == 0

    def test_liability_severity_ordering(self) -> None:
        """Liabilities should be sorted: highest severity first."""
        numbered = {"10": "W", "50": "F"}
        exposure = {
            "10": {"exposure": "exposed", "region": "FR1"},
            "50": {"exposure": "partially_exposed", "region": "FR2"},
        }
        liabilities = detect_liabilities(numbered, exposure)
        assert len(liabilities) == 2
        # Exposed = severity 3, partially_exposed = severity 2
        assert liabilities[0].severity >= liabilities[1].severity

    def test_detect_liabilities_returns_sorted(self) -> None:
        """Output list should be sorted by the Liability dataclass ordering."""
        numbered = {"10": "W", "50": "F", "60": "L"}
        exposure = {
            "10": {"exposure": "exposed", "region": "FR1"},
            "50": {"exposure": "partially_exposed", "region": "FR2"},
            "60": {"exposure": "exposed", "region": "CDR2"},
        }
        liabilities = detect_liabilities(numbered, exposure)
        assert liabilities == sorted(liabilities)

    def test_buried_positions_not_flagged(self) -> None:
        """Hydrophobic residues at buried positions should not be flagged."""
        numbered = {"10": "W"}
        exposure = {"10": {"exposure": "buried", "region": "FR1"}}
        liabilities = detect_liabilities(numbered, exposure)
        assert len(liabilities) == 0

    def test_liability_dataclass_fields(self) -> None:
        """Verify that the Liability dataclass has the expected fields."""
        liab = Liability(
            position="56",
            residue="W",
            imgt_number="56",
            region="CDR2",
            reason="test reason",
            severity=3,
        )
        assert liab.position == "56"
        assert liab.residue == "W"
        assert liab.imgt_number == "56"
        assert liab.region == "CDR2"
        assert liab.severity == 3

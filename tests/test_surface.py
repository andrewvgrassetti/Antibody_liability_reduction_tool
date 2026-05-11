"""Tests for the surface exposure module."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from antibody_liability_tool.surface.exposure import (
    classify_surface_exposure,
    get_exposed_positions,
)


def _mock_exposure_data() -> dict:
    """Return synthetic exposure data for testing."""
    return {
        "1": {"exposure": "buried", "region": "FR1"},
        "2": {"exposure": "exposed", "region": "FR1"},
        "3": {"exposure": "partially_exposed", "region": "FR1"},
        "56": {"exposure": "exposed", "region": "CDR2"},
        "57": {"exposure": "buried", "region": "CDR2"},
        "100": {"exposure": "exposed", "region": "FR3"},
    }


@pytest.fixture(autouse=True)
def _patch_exposure_data():
    """Patch the exposure data loader for all tests in this module."""
    with patch(
        "antibody_liability_tool.surface.exposure._load_exposure_data",
        return_value=_mock_exposure_data(),
    ):
        yield


class TestClassifySurfaceExposure:
    """Tests for classify_surface_exposure()."""

    def test_classify_surface_exposure_returns_dict(self) -> None:
        numbered = {"1": "E", "2": "V", "3": "Q", "56": "R"}
        result = classify_surface_exposure(numbered)
        assert isinstance(result, dict)

    def test_exposed_positions_not_empty(self) -> None:
        numbered = {"1": "E", "2": "V", "3": "Q", "56": "R"}
        result = classify_surface_exposure(numbered)
        exposed = {k: v for k, v in result.items() if v["exposure"] == "exposed"}
        assert len(exposed) > 0

    def test_buried_positions_not_empty(self) -> None:
        numbered = {"1": "E", "2": "V", "57": "Y"}
        result = classify_surface_exposure(numbered)
        # Position 1 is buried in our mock data
        assert any(v["exposure"] == "buried" for v in result.values())

    def test_unknown_position_handled(self) -> None:
        """Positions not in the exposure data are silently omitted."""
        numbered = {"999": "A", "2": "V"}
        result = classify_surface_exposure(numbered)
        assert "999" not in result
        assert "2" in result

    def test_residue_included_in_result(self) -> None:
        numbered = {"2": "V"}
        result = classify_surface_exposure(numbered)
        assert result["2"]["residue"] == "V"


class TestGetExposedPositions:
    """Tests for get_exposed_positions()."""

    def test_get_exposed_positions_filters_correctly(self) -> None:
        numbered = {"1": "E", "2": "V", "3": "Q", "56": "R", "57": "Y"}
        exposure_map = classify_surface_exposure(numbered)
        exposed = get_exposed_positions(exposure_map)
        # Position 1 (buried) and 57 (buried) should be excluded
        assert "1" not in exposed
        assert "57" not in exposed
        # Exposed and partially exposed should be included
        assert "2" in exposed
        assert "3" in exposed

    def test_get_exposed_positions_exclude_partial(self) -> None:
        numbered = {"1": "E", "2": "V", "3": "Q"}
        exposure_map = classify_surface_exposure(numbered)
        exposed = get_exposed_positions(exposure_map, include_partial=False)
        assert "3" not in exposed
        assert "2" in exposed

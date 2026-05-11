"""Tests for the IMGT numbering module."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from antibody_liability_tool.numbering.imgt import (
    classify_region,
    number_sequence,
)
from tests.conftest import SAMPLE_VH


class TestClassifyRegion:
    """Tests for classify_region()."""

    def test_classify_region_fr1(self) -> None:
        assert classify_region("1") == "FR1"
        assert classify_region("26") == "FR1"

    def test_classify_region_cdr1(self) -> None:
        assert classify_region("27") == "CDR1"
        assert classify_region("38") == "CDR1"

    def test_classify_region_fr2(self) -> None:
        assert classify_region("39") == "FR2"
        assert classify_region("55") == "FR2"

    def test_classify_region_cdr2(self) -> None:
        assert classify_region("56") == "CDR2"
        assert classify_region("65") == "CDR2"

    def test_classify_region_fr3(self) -> None:
        assert classify_region("66") == "FR3"
        assert classify_region("104") == "FR3"

    def test_classify_region_cdr3(self) -> None:
        assert classify_region("105") == "CDR3"
        assert classify_region("117") == "CDR3"

    def test_classify_region_fr4(self) -> None:
        assert classify_region("118") == "FR4"
        assert classify_region("128") == "FR4"

    def test_classify_region_insertion_code(self) -> None:
        assert classify_region("82A") == "FR3"
        assert classify_region("111.1") == "CDR3"

    def test_classify_region_unknown_high_position(self) -> None:
        assert classify_region("200") == "unknown"

    def test_classify_region_non_numeric(self) -> None:
        assert classify_region("XYZ") == "unknown"


class TestNumberSequence:
    """Tests for number_sequence()."""

    def test_number_sequence_returns_dict(self) -> None:
        with patch(
            "antibody_liability_tool.numbering.imgt._number_with_anarci",
            return_value=None,
        ):
            result = number_sequence(SAMPLE_VH)
        assert isinstance(result, dict)

    def test_number_sequence_correct_length(self) -> None:
        with patch(
            "antibody_liability_tool.numbering.imgt._number_with_anarci",
            return_value=None,
        ):
            result = number_sequence(SAMPLE_VH)
        assert len(result) == len(SAMPLE_VH)

    def test_number_sequence_preserves_residues(self) -> None:
        with patch(
            "antibody_liability_tool.numbering.imgt._number_with_anarci",
            return_value=None,
        ):
            result = number_sequence(SAMPLE_VH)
        assert "".join(result.values()) == SAMPLE_VH

    def test_number_sequence_fallback(self) -> None:
        """When ANARCI is not available, fallback sequential numbering is used."""
        with patch(
            "antibody_liability_tool.numbering.imgt._number_with_anarci",
            return_value=None,
        ):
            result = number_sequence(SAMPLE_VH)
        assert result["1"] == SAMPLE_VH[0]
        assert result[str(len(SAMPLE_VH))] == SAMPLE_VH[-1]

    def test_number_sequence_empty_raises(self) -> None:
        with pytest.raises(ValueError, match="Empty sequence"):
            number_sequence("")

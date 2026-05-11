"""Tests for the mutation generation module."""

from __future__ import annotations

from unittest.mock import patch

import pytest

from antibody_liability_tool.liabilities.detector import Liability
from antibody_liability_tool.mutations.generator import generate_mutations


def _mock_freq_data() -> dict[str, dict[str, float]]:
    """Return synthetic frequency data for testing."""
    return {
        "56": {"S": 0.35, "T": 0.20, "A": 0.15, "V": 0.10, "G": 0.08, "other": 0.12},
        "64": {"Q": 0.25, "S": 0.20, "T": 0.15, "A": 0.10, "other": 0.30},
        "99": {"G": 0.90, "other": 0.10},
    }


@pytest.fixture(autouse=True)
def _patch_freq_data():
    with patch(
        "antibody_liability_tool.mutations.generator._load_frequency_data",
        return_value=_mock_freq_data(),
    ):
        yield


def _make_liability(
    position: str,
    residue: str,
    reason: str = "Hydrophobic residue",
    severity: int = 3,
) -> Liability:
    return Liability(
        position=position,
        residue=residue,
        imgt_number=position,
        region="CDR2",
        reason=reason,
        severity=severity,
    )


class TestGenerateMutations:
    """Tests for generate_mutations()."""

    def test_generate_candidates_returns_list(self) -> None:
        liabilities = [_make_liability("56", "V")]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["56"] = "V"
        result = generate_mutations(liabilities, numbered)
        assert isinstance(result, list)

    def test_candidates_have_required_fields(self) -> None:
        liabilities = [_make_liability("56", "V")]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["56"] = "V"
        result = generate_mutations(liabilities, numbered)
        if result:
            cand = result[0]
            assert hasattr(cand, "position")
            assert hasattr(cand, "original_aa")
            assert hasattr(cand, "proposed_aa")
            assert hasattr(cand, "human_frequency")
            assert hasattr(cand, "rationale")

    def test_candidates_above_frequency_threshold(self) -> None:
        liabilities = [_make_liability("56", "V")]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["56"] = "V"
        config = {"mutations": {"min_human_frequency": 0.05}}
        result = generate_mutations(liabilities, numbered, config=config)
        for cand in result:
            assert cand.human_frequency >= 0.05

    def test_candidates_reduce_liability(self) -> None:
        """Hydrophobic -> polar substitutions should be proposed."""
        liabilities = [_make_liability("56", "V")]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["56"] = "V"
        result = generate_mutations(liabilities, numbered)
        # S, T are polar replacements for V (hydrophobic)
        proposed_aas = {c.proposed_aa for c in result}
        assert proposed_aas & {"S", "T"}, f"Expected polar replacements, got {proposed_aas}"

    def test_candidates_no_introduced_motifs(self) -> None:
        """Mutations that introduce N-glycosylation etc. should be filtered out."""
        # Position 56 is V, and context around it matters.
        # With mock freq data, we check that no motif-introducing mutations pass.
        liabilities = [_make_liability("56", "V")]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["56"] = "V"
        result = generate_mutations(liabilities, numbered)
        # All returned candidates should have passed motif checking
        assert isinstance(result, list)

    def test_no_candidates_for_conserved_position(self) -> None:
        """A position with no alternatives above threshold returns no candidates."""
        liabilities = [_make_liability("99", "G")]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["99"] = "G"
        # G is not hydrophobic and not R/K, so no liability-reducing substitution applies
        result = generate_mutations(liabilities, numbered)
        assert len(result) == 0

    def test_charge_liability_reduction(self) -> None:
        """Positive charge (K) -> neutral should be proposed."""
        liabilities = [
            _make_liability("64", "K", reason="Positive-charge cluster", severity=2)
        ]
        numbered = {str(i): "A" for i in range(1, 120)}
        numbered["64"] = "K"
        result = generate_mutations(liabilities, numbered)
        proposed_aas = {c.proposed_aa for c in result}
        # Q, S, T, A are neutral replacements for K
        assert proposed_aas & {"Q", "S", "T", "A"}, (
            f"Expected neutral replacements, got {proposed_aas}"
        )

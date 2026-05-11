"""Tests for the combinatorial expansion module."""

from __future__ import annotations

from antibody_liability_tool.mutations.generator import CandidateMutation
from antibody_liability_tool.optimization.combinatorial import (
    MutationCombination,
    generate_combinations,
)


def _make_candidate(position: str, original: str, proposed: str) -> CandidateMutation:
    return CandidateMutation(
        position=position,
        original_aa=original,
        proposed_aa=proposed,
        human_frequency=0.3,
        rationale=f"Replace {original} with {proposed}",
    )


class TestGenerateCombinations:
    """Tests for generate_combinations()."""

    def test_generate_doubles(self) -> None:
        mutations = [
            _make_candidate("56", "V", "S"),
            _make_candidate("64", "K", "Q"),
        ]
        combos = generate_combinations(mutations, max_order=2)
        doubles = [c for c in combos if c.order == 2]
        assert len(doubles) >= 1
        # The double should combine positions 56 and 64
        for d in doubles:
            positions = {m.position for m in d.mutations}
            assert len(positions) == 2

    def test_generate_triples(self) -> None:
        mutations = [
            _make_candidate("56", "V", "S"),
            _make_candidate("64", "K", "Q"),
            _make_candidate("80", "F", "T"),
        ]
        combos = generate_combinations(mutations, max_order=3)
        triples = [c for c in combos if c.order == 3]
        assert len(triples) >= 1

    def test_no_same_position_combinations(self) -> None:
        """Multiple mutations at the same position should not be combined together."""
        mutations = [
            _make_candidate("56", "V", "S"),
            _make_candidate("56", "V", "T"),
            _make_candidate("64", "K", "Q"),
        ]
        combos = generate_combinations(mutations, max_order=2)
        for combo in combos:
            if combo.order > 1:
                positions = [m.position for m in combo.mutations]
                assert len(positions) == len(set(positions)), (
                    f"Duplicate position in combination: {positions}"
                )

    def test_empty_input(self) -> None:
        combos = generate_combinations([], max_order=3)
        assert len(combos) == 0

    def test_single_position_returns_singles_only(self) -> None:
        """With only one position, no multi-mutants can be generated."""
        mutations = [
            _make_candidate("56", "V", "S"),
            _make_candidate("56", "V", "T"),
        ]
        combos = generate_combinations(mutations, max_order=3)
        for combo in combos:
            assert combo.order == 1

    def test_combination_label(self) -> None:
        mutations = [
            _make_candidate("56", "V", "S"),
            _make_candidate("64", "K", "Q"),
        ]
        combos = generate_combinations(mutations, max_order=2)
        doubles = [c for c in combos if c.order == 2]
        assert len(doubles) >= 1
        # Label should contain both mutations
        assert "V56S" in doubles[0].label
        assert "K64Q" in doubles[0].label

    def test_apply_mutation(self) -> None:
        """MutationCombination.apply() should produce a mutated sequence."""
        m1 = _make_candidate("2", "V", "S")
        combo = MutationCombination(mutations=(m1,))
        numbered = {"1": "E", "2": "V", "3": "Q"}
        result = combo.apply("EVQ", numbered)
        assert result == "ESQ"

"""Combinatorial expansion of validated single mutations into multi-mutant candidates.

Generates all double and triple mutation combinations from validated single
mutations, ensuring that only mutations at *different* IMGT positions are
combined.  When the combinatorial space exceeds a configurable threshold the
module delegates to :mod:`antibody_liability_tool.optimization.bayesian` for
efficient exploration.
"""

from __future__ import annotations

import itertools
import logging
from collections import defaultdict
from dataclasses import dataclass, field
from typing import Any, Callable, Sequence

from antibody_liability_tool.mutations.generator import CandidateMutation

logger = logging.getLogger(__name__)


@dataclass
class MutationCombination:
    """A combination of one or more single-point mutations.

    Attributes
    ----------
    mutations : tuple[CandidateMutation, ...]
        The individual mutations in this combination.
    order : int
        Number of simultaneous mutations (1 = single, 2 = double, etc.).
    label : str
        Human-readable label such as ``"V56S+K64Q"``.
    """

    mutations: tuple[CandidateMutation, ...] = ()
    order: int = 0
    label: str = ""
    _positions: frozenset[str] = field(init=False, repr=False, compare=False)

    def __post_init__(self) -> None:
        self._positions = frozenset(m.position for m in self.mutations)
        if not self.order:
            self.order = len(self.mutations)
        if not self.label:
            self.label = "+".join(
                f"{m.original_aa}{m.position}{m.proposed_aa}" for m in self.mutations
            )

    @property
    def positions(self) -> frozenset[str]:
        """Set of IMGT positions affected by this combination."""
        return self._positions

    def apply(self, sequence: str, numbered: dict[str, str]) -> str:
        """Apply this combination to a numbered sequence and return the mutated string.

        Parameters
        ----------
        sequence : str
            The original linear amino-acid sequence.
        numbered : dict[str, str]
            IMGT position → amino acid mapping for the parent sequence.

        Returns
        -------
        str
            Mutated sequence.
        """
        mutated = dict(numbered)
        for m in self.mutations:
            mutated[m.position] = m.proposed_aa

        def _sort_key(p: str) -> tuple[int, str]:
            num, suffix = "", ""
            for ch in p:
                if ch.isdigit():
                    num += ch
                else:
                    suffix += ch
            return (int(num) if num else 0, suffix)

        return "".join(mutated[p] for p in sorted(mutated, key=_sort_key))


def _group_by_position(
    mutations: Sequence[CandidateMutation],
) -> dict[str, list[CandidateMutation]]:
    """Group mutations by their IMGT position."""
    groups: dict[str, list[CandidateMutation]] = defaultdict(list)
    for m in mutations:
        groups[m.position].append(m)
    return dict(groups)


def generate_combinations(
    validated_mutations: Sequence[CandidateMutation],
    max_order: int = 3,
    bayesian_threshold: int = 50,
    evaluator_func: Callable[..., Any] | None = None,
    config: dict[str, Any] | None = None,
) -> list[MutationCombination]:
    """Expand validated single mutations into multi-mutant combinations.

    Parameters
    ----------
    validated_mutations : Sequence[CandidateMutation]
        Single-point mutations that have passed all validation filters.
    max_order : int
        Maximum combination order (2 = doubles only, 3 = doubles + triples).
    bayesian_threshold : int
        If the total number of combinations exceeds this threshold,
        Bayesian optimisation is used instead of exhaustive enumeration.
    evaluator_func : callable, optional
        Evaluation function for Bayesian optimisation (see
        :class:`~antibody_liability_tool.optimization.bayesian.BayesianMutationOptimizer`).
    config : dict, optional
        Full project configuration dict.

    Returns
    -------
    list[MutationCombination]
        All generated mutation combinations (singles included).
    """
    if config is None:
        config = {}

    position_groups = _group_by_position(validated_mutations)
    positions = sorted(position_groups.keys())
    n_positions = len(positions)

    logger.info(
        "Combinatorial expansion: %d mutations across %d positions (max_order=%d)",
        len(validated_mutations),
        n_positions,
        max_order,
    )

    # --- Singles ---
    singles: list[MutationCombination] = [
        MutationCombination(mutations=(m,)) for m in validated_mutations
    ]

    if n_positions < 2 or max_order < 2:
        logger.info("Returning %d single-mutation candidates only", len(singles))
        return singles

    # --- Estimate total combinations ---
    total = _estimate_combinations(position_groups, positions, max_order)
    logger.info("Estimated %d total multi-mutant combinations", total)

    if total > bayesian_threshold and evaluator_func is not None:
        logger.info(
            "Combination count (%d) exceeds threshold (%d) – using Bayesian optimisation",
            total,
            bayesian_threshold,
        )
        return _bayesian_expansion(
            singles, position_groups, positions, max_order, evaluator_func, config
        )

    # --- Exhaustive enumeration ---
    combos = list(singles)
    combos.extend(_enumerate_combos(position_groups, positions, max_order))

    logger.info("Generated %d total candidates (singles + combinations)", len(combos))
    return combos


def _estimate_combinations(
    position_groups: dict[str, list[CandidateMutation]],
    positions: list[str],
    max_order: int,
) -> int:
    """Estimate the number of combinations without generating them."""
    total = 0
    sizes = [len(position_groups[p]) for p in positions]
    for order in range(2, min(max_order, len(positions)) + 1):
        for idx_combo in itertools.combinations(range(len(positions)), order):
            product = 1
            for idx in idx_combo:
                product *= sizes[idx]
            total += product
    return total


def _enumerate_combos(
    position_groups: dict[str, list[CandidateMutation]],
    positions: list[str],
    max_order: int,
) -> list[MutationCombination]:
    """Exhaustively enumerate all multi-mutant combinations."""
    results: list[MutationCombination] = []

    for order in range(2, min(max_order, len(positions)) + 1):
        for pos_combo in itertools.combinations(positions, order):
            per_position = [position_groups[p] for p in pos_combo]
            for mut_combo in itertools.product(*per_position):
                results.append(MutationCombination(mutations=tuple(mut_combo)))

    return results


def _bayesian_expansion(
    singles: list[MutationCombination],
    position_groups: dict[str, list[CandidateMutation]],
    positions: list[str],
    max_order: int,
    evaluator_func: Callable[..., Any],
    config: dict[str, Any],
) -> list[MutationCombination]:
    """Use Bayesian optimisation to explore the combinatorial space."""
    from antibody_liability_tool.optimization.bayesian import BayesianMutationOptimizer

    bcfg = config.get("combinatorial", {}).get("bayesian", {})
    ref_point = bcfg.get("reference_point", [-1.0, -1.0, 0.0])

    optimizer = BayesianMutationOptimizer(
        evaluator_func=evaluator_func,
        reference_point=ref_point,
    )

    # Build search space: list of all possible combos (as indices)
    search_space: list[MutationCombination] = list(singles)
    search_space.extend(_enumerate_combos(position_groups, positions, max_order))

    results = optimizer.optimize(
        search_space=search_space,
        n_initial=bcfg.get("n_initial_samples", 20),
        n_iterations=bcfg.get("n_iterations", 30),
        batch_size=bcfg.get("batch_size", 5),
    )

    logger.info(
        "Bayesian optimisation selected %d candidates from %d possibilities",
        len(results),
        len(search_space),
    )
    return results

"""Bayesian optimisation for multi-objective mutation space exploration.

Uses BoTorch with a Gaussian Process surrogate and Expected Hypervolume
Improvement (EHVI) when available.  Falls back to random sampling with
scoring when BoTorch / PyTorch are not installed.
"""

from __future__ import annotations

import logging
import random
from dataclasses import dataclass
from typing import TYPE_CHECKING, Any, Callable, Sequence

import numpy as np

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Optional BoTorch / PyTorch imports
# ---------------------------------------------------------------------------
_BOTORCH_AVAILABLE = False
try:
    import torch
    from botorch.acquisition.multi_objective import (  # type: ignore[import-untyped]
        qExpectedHypervolumeImprovement,
    )
    from botorch.fit import fit_gpytorch_mll  # type: ignore[import-untyped]
    from botorch.models.gp_regression import SingleTaskGP  # type: ignore[import-untyped]
    from botorch.models.model_list_gp_regression import ModelListGP  # type: ignore[import-untyped]
    from botorch.models.transforms.outcome import Standardize  # type: ignore[import-untyped]
    from botorch.utils.multi_objective.box_decompositions.non_dominated import (  # type: ignore[import-untyped]
        FastNondominatedPartitioning,
    )
    from gpytorch.mlls.sum_marginal_log_likelihood import (  # type: ignore[import-untyped]
        SumMarginalLogLikelihood,
    )

    _BOTORCH_AVAILABLE = True
    logger.debug("BoTorch available – using GP-based Bayesian optimisation.")
except ImportError:
    logger.debug("BoTorch/PyTorch not installed – falling back to random sampling.")

if TYPE_CHECKING:
    from antibody_liability_tool.optimization.combinatorial import MutationCombination


# ---------------------------------------------------------------------------
# Feature encoding
# ---------------------------------------------------------------------------

# Standard amino-acid physicochemical properties for encoding.
_AA_PROPERTIES: dict[str, tuple[float, float, float, float]] = {
    # (hydrophobicity, charge, size, polarity) – normalised 0-1
    "A": (0.62, 0.0, 0.17, 0.0),
    "R": (0.0, 1.0, 0.60, 1.0),
    "N": (0.14, 0.0, 0.37, 1.0),
    "D": (0.15, -1.0, 0.36, 1.0),
    "C": (0.35, 0.0, 0.30, 0.0),
    "E": (0.16, -1.0, 0.47, 1.0),
    "Q": (0.17, 0.0, 0.48, 1.0),
    "G": (0.50, 0.0, 0.0, 0.0),
    "H": (0.23, 0.5, 0.50, 1.0),
    "I": (1.0, 0.0, 0.46, 0.0),
    "L": (0.94, 0.0, 0.46, 0.0),
    "K": (0.07, 1.0, 0.53, 1.0),
    "M": (0.74, 0.0, 0.50, 0.0),
    "F": (1.0, 0.0, 0.63, 0.0),
    "P": (0.32, 0.0, 0.30, 0.0),
    "S": (0.25, 0.0, 0.17, 1.0),
    "T": (0.27, 0.0, 0.30, 1.0),
    "W": (0.88, 0.0, 0.80, 0.0),
    "Y": (0.63, 0.0, 0.67, 1.0),
    "V": (0.86, 0.0, 0.33, 0.0),
}
_DEFAULT_PROPS = (0.5, 0.0, 0.4, 0.5)


def encode_mutations(
    combinations: Sequence[MutationCombination],
    max_mutations: int = 3,
) -> np.ndarray:
    """Convert mutation combinations to numeric feature vectors.

    Each mutation contributes four physicochemical features (original AA)
    plus four features (proposed AA) plus a normalised position index.
    Unused mutation slots are zero-padded.

    Parameters
    ----------
    combinations : Sequence[MutationCombination]
        Mutation combinations to encode.
    max_mutations : int
        Maximum number of mutation slots per combination (determines
        the feature vector width).

    Returns
    -------
    np.ndarray
        Feature matrix of shape ``(len(combinations), max_mutations * 9)``.
    """
    n_features = max_mutations * 9  # 4 orig + 4 proposed + 1 position per slot
    X = np.zeros((len(combinations), n_features), dtype=np.float64)

    for i, combo in enumerate(combinations):
        for j, mut in enumerate(combo.mutations):
            if j >= max_mutations:
                break
            offset = j * 9
            orig_props = _AA_PROPERTIES.get(mut.original_aa, _DEFAULT_PROPS)
            prop_props = _AA_PROPERTIES.get(mut.proposed_aa, _DEFAULT_PROPS)
            X[i, offset : offset + 4] = orig_props
            X[i, offset + 4 : offset + 8] = prop_props
            # Normalised IMGT position (extract integer part, divide by 128)
            pos_num = "".join(ch for ch in mut.position if ch.isdigit())
            X[i, offset + 8] = int(pos_num) / 128.0 if pos_num else 0.0

    return X


# ---------------------------------------------------------------------------
# Bayesian Mutation Optimizer
# ---------------------------------------------------------------------------


@dataclass
class OptimizationResult:
    """Result container for a single evaluated candidate."""

    combination: Any  # MutationCombination
    objectives: np.ndarray  # shape (n_objectives,)
    score: float = 0.0


class BayesianMutationOptimizer:
    """Multi-objective Bayesian optimisation over the mutation space.

    Minimises PSH, minimises PPC, and maximises OASis humanness using
    Expected Hypervolume Improvement (EHVI).

    Parameters
    ----------
    evaluator_func : callable
        Function that takes a :class:`MutationCombination` and returns a
        dict with keys ``"PSH"``, ``"PPC"``, ``"oasis_humanness"`` (floats).
    reference_point : list[float]
        Reference point for hypervolume computation.  Objectives are
        transformed so that all are *minimised*: ``[-PSH, -PPC, +OASis]``.
    n_objectives : int
        Number of objectives (default 3).
    """

    def __init__(
        self,
        evaluator_func: Callable[..., dict[str, float]],
        reference_point: list[float] | None = None,
        n_objectives: int = 3,
    ) -> None:
        self.evaluator_func = evaluator_func
        self.n_objectives = n_objectives
        self.reference_point = reference_point or [-1.0, -1.0, 0.0]
        self._use_botorch = _BOTORCH_AVAILABLE

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def optimize(
        self,
        search_space: Sequence[MutationCombination],
        n_initial: int = 20,
        n_iterations: int = 30,
        batch_size: int = 5,
    ) -> list[MutationCombination]:
        """Run multi-objective optimisation over the search space.

        Parameters
        ----------
        search_space : Sequence[MutationCombination]
            All candidate mutation combinations.
        n_initial : int
            Number of initial random evaluations.
        n_iterations : int
            Number of Bayesian optimisation iterations.
        batch_size : int
            Candidates evaluated per iteration.

        Returns
        -------
        list[MutationCombination]
            Ranked list of the best candidates found.
        """
        if len(search_space) == 0:
            logger.warning("Empty search space – returning empty results.")
            return []

        if self._use_botorch:
            return self._optimize_botorch(search_space, n_initial, n_iterations, batch_size)
        return self._optimize_random(search_space, n_initial, n_iterations, batch_size)

    # ------------------------------------------------------------------
    # BoTorch implementation
    # ------------------------------------------------------------------

    def _optimize_botorch(
        self,
        search_space: Sequence[MutationCombination],
        n_initial: int,
        n_iterations: int,
        batch_size: int,
    ) -> list[MutationCombination]:
        """GP + EHVI multi-objective optimisation."""
        space = list(search_space)
        n_initial = min(n_initial, len(space))

        # Encode full search space
        X_all = encode_mutations(space)
        X_all_tensor = torch.tensor(X_all, dtype=torch.double)  # type: ignore[possibly-undefined]

        # Initial random evaluations
        initial_indices = random.sample(range(len(space)), n_initial)
        evaluated_indices: list[int] = list(initial_indices)
        objectives_list: list[list[float]] = []

        for idx in initial_indices:
            obj = self._evaluate_candidate(space[idx])
            objectives_list.append(obj)

        logger.info("Completed %d initial evaluations", n_initial)

        # Iterative BO loop
        for iteration in range(n_iterations):
            X_train = X_all_tensor[evaluated_indices]
            Y_train = torch.tensor(objectives_list, dtype=torch.double)

            # Build independent GP per objective
            models = []
            for obj_idx in range(self.n_objectives):
                y = Y_train[:, obj_idx].unsqueeze(-1)
                gp = SingleTaskGP(X_train, y, outcome_transform=Standardize(m=1))
                models.append(gp)

            model = ModelListGP(*models)
            mll = SumMarginalLogLikelihood(model.likelihood, model)

            try:
                fit_gpytorch_mll(mll)
            except Exception as exc:  # noqa: BLE001
                logger.warning("GP fitting failed at iteration %d: %s", iteration, exc)
                continue

            # EHVI acquisition
            ref_point_tensor = torch.tensor(self.reference_point, dtype=torch.double)
            partitioning = FastNondominatedPartitioning(
                ref_point=ref_point_tensor,
                Y=Y_train,
            )

            acq = qExpectedHypervolumeImprovement(
                model=model,
                ref_point=ref_point_tensor,
                partitioning=partitioning,
                sampler=None,
            )

            # Evaluate acquisition function on unevaluated candidates
            unevaluated = [i for i in range(len(space)) if i not in set(evaluated_indices)]
            if not unevaluated:
                logger.info("All candidates evaluated at iteration %d", iteration)
                break

            # Score unevaluated candidates and pick top batch_size
            acq_values: list[tuple[int, float]] = []
            for idx in unevaluated:
                x = X_all_tensor[idx].unsqueeze(0).unsqueeze(0)
                with torch.no_grad():
                    val = acq(x).item()
                acq_values.append((idx, val))

            acq_values.sort(key=lambda t: -t[1])
            selected = [t[0] for t in acq_values[:batch_size]]

            for idx in selected:
                obj = self._evaluate_candidate(space[idx])
                evaluated_indices.append(idx)
                objectives_list.append(obj)

            logger.debug(
                "Iteration %d/%d: evaluated %d candidates",
                iteration + 1,
                n_iterations,
                len(selected),
            )

        # Rank all evaluated candidates by Pareto dominance
        return self._rank_evaluated(space, evaluated_indices, objectives_list)

    # ------------------------------------------------------------------
    # Random sampling fallback
    # ------------------------------------------------------------------

    def _optimize_random(
        self,
        search_space: Sequence[MutationCombination],
        n_initial: int,
        n_iterations: int,
        batch_size: int,
    ) -> list[MutationCombination]:
        """Random sampling fallback when BoTorch is unavailable."""
        logger.info("Using random sampling fallback (BoTorch not available)")
        space = list(search_space)
        total_budget = min(n_initial + n_iterations * batch_size, len(space))
        sample_indices = random.sample(range(len(space)), total_budget)

        evaluated_indices: list[int] = []
        objectives_list: list[list[float]] = []

        for idx in sample_indices:
            obj = self._evaluate_candidate(space[idx])
            evaluated_indices.append(idx)
            objectives_list.append(obj)

        logger.info("Random sampling evaluated %d candidates", len(evaluated_indices))
        return self._rank_evaluated(space, evaluated_indices, objectives_list)

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _evaluate_candidate(self, combo: MutationCombination) -> list[float]:
        """Evaluate a single candidate and return transformed objectives.

        Objectives are transformed so all should be *minimised*:
        - PSH → negative (lower is better, so negate to minimise)
        - PPC → negative
        - OASis → positive (higher is better, so keep positive for min)

        Wait – for EHVI we actually want to *maximise* the hypervolume.
        Convention: negate metrics that should be minimised so all objectives
        are "higher is better" for the GP / EHVI framework.
        """
        metrics = self.evaluator_func(combo)
        psh = metrics.get("PSH", 0.0)
        ppc = metrics.get("PPC", 0.0)
        oasis = metrics.get("oasis_humanness", 0.0)
        # Transform: negate PSH and PPC so "higher is better" for all
        return [-psh, -ppc, oasis]

    def _rank_evaluated(
        self,
        space: list[MutationCombination],
        indices: list[int],
        objectives: list[list[float]],
    ) -> list[MutationCombination]:
        """Rank candidates by a simple weighted sum of transformed objectives."""
        scored: list[tuple[float, int]] = []
        for i, obj in zip(indices, objectives):
            # Equal weighting across normalised objectives
            score = sum(obj) / len(obj)
            scored.append((score, i))

        scored.sort(key=lambda t: -t[0])  # highest combined score first
        return [space[idx] for _, idx in scored]

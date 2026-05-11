"""Abstract base classes for antibody sequence evaluators."""

from __future__ import annotations

import asyncio
import logging
from abc import ABC, abstractmethod
from dataclasses import dataclass, field
from typing import Any, Sequence

logger = logging.getLogger(__name__)


@dataclass
class EvaluationResult:
    """Container for a single evaluation output.

    Attributes
    ----------
    sequence : str
        The input amino-acid sequence that was evaluated.
    evaluator_name : str
        Name of the evaluator that produced this result.
    metrics : dict[str, float]
        Metric name → numeric value mapping.
    raw : dict[str, Any]
        Raw / unprocessed data returned by the evaluator.
    success : bool
        Whether the evaluation completed without error.
    error_message : str
        Human-readable error description when *success* is ``False``.
    """

    sequence: str
    evaluator_name: str
    metrics: dict[str, float] = field(default_factory=dict)
    raw: dict[str, Any] = field(default_factory=dict)
    success: bool = True
    error_message: str = ""


@dataclass
class ComparisonResult:
    """Result of comparing a mutant evaluation to its parent.

    Attributes
    ----------
    parent : EvaluationResult
        Evaluation of the parent (wild-type) sequence.
    mutant : EvaluationResult
        Evaluation of the mutant sequence.
    deltas : dict[str, float]
        Metric name → (mutant − parent) value.
    improved_metrics : list[str]
        Metrics where the mutant is better (lower liability).
    worsened_metrics : list[str]
        Metrics where the mutant is worse.
    passes_filter : bool
        Whether the mutant passes the evaluator's acceptance criteria.
    """

    parent: EvaluationResult
    mutant: EvaluationResult
    deltas: dict[str, float] = field(default_factory=dict)
    improved_metrics: list[str] = field(default_factory=list)
    worsened_metrics: list[str] = field(default_factory=list)
    passes_filter: bool = True


class BaseEvaluator(ABC):
    """Abstract base class that every evaluator must implement.

    Parameters
    ----------
    config : dict[str, Any]
        Evaluator-specific configuration section from the project config.
    name : str | None
        Human-readable evaluator name (defaults to the class name).
    """

    def __init__(self, config: dict[str, Any] | None = None, name: str | None = None) -> None:
        self.config: dict[str, Any] = config or {}
        self.name: str = name or self.__class__.__name__
        self._logger = logging.getLogger(f"{__name__}.{self.name}")

    # ------------------------------------------------------------------
    # Core abstract interface
    # ------------------------------------------------------------------

    @abstractmethod
    def evaluate(self, sequence: str) -> EvaluationResult:
        """Evaluate a single antibody sequence.

        Parameters
        ----------
        sequence : str
            Amino-acid sequence (VH or VL).

        Returns
        -------
        EvaluationResult
        """

    def evaluate_batch(self, sequences: Sequence[str]) -> list[EvaluationResult]:
        """Evaluate multiple sequences.

        The default implementation iterates sequentially.  Subclasses may
        override for true batch / parallel support.

        Parameters
        ----------
        sequences : Sequence[str]
            Iterable of amino-acid sequences.

        Returns
        -------
        list[EvaluationResult]
        """
        results: list[EvaluationResult] = []
        for seq in sequences:
            try:
                results.append(self.evaluate(seq))
            except Exception as exc:  # noqa: BLE001
                self._logger.error("Evaluation failed for sequence %.20s…: %s", seq, exc)
                results.append(
                    EvaluationResult(
                        sequence=seq,
                        evaluator_name=self.name,
                        success=False,
                        error_message=str(exc),
                    )
                )
        return results

    def compare(
        self,
        parent_result: EvaluationResult,
        mutant_result: EvaluationResult,
        *,
        lower_is_better: set[str] | None = None,
        higher_is_better: set[str] | None = None,
    ) -> ComparisonResult:
        """Compare a mutant evaluation against its parent.

        Parameters
        ----------
        parent_result, mutant_result : EvaluationResult
            Results from :meth:`evaluate`.
        lower_is_better : set[str] | None
            Metric names where a decrease is desirable.
        higher_is_better : set[str] | None
            Metric names where an increase is desirable.

        Returns
        -------
        ComparisonResult
        """
        lower = lower_is_better or set()
        higher = higher_is_better or set()

        deltas: dict[str, float] = {}
        improved: list[str] = []
        worsened: list[str] = []

        all_keys = set(parent_result.metrics) | set(mutant_result.metrics)
        for key in sorted(all_keys):
            p_val = parent_result.metrics.get(key, 0.0)
            m_val = mutant_result.metrics.get(key, 0.0)
            delta = m_val - p_val
            deltas[key] = delta

            if key in lower:
                if delta < 0:
                    improved.append(key)
                elif delta > 0:
                    worsened.append(key)
            elif key in higher:
                if delta > 0:
                    improved.append(key)
                elif delta < 0:
                    worsened.append(key)

        return ComparisonResult(
            parent=parent_result,
            mutant=mutant_result,
            deltas=deltas,
            improved_metrics=improved,
            worsened_metrics=worsened,
            passes_filter=len(worsened) == 0,
        )

    # ------------------------------------------------------------------
    # Async helpers
    # ------------------------------------------------------------------

    async def evaluate_async(self, sequence: str) -> EvaluationResult:
        """Evaluate a sequence asynchronously (runs in a thread executor).

        Subclasses that have native async support should override.
        """
        loop = asyncio.get_running_loop()
        return await loop.run_in_executor(None, self.evaluate, sequence)

    async def evaluate_batch_async(self, sequences: Sequence[str]) -> list[EvaluationResult]:
        """Evaluate a batch asynchronously, running evaluations concurrently."""
        tasks = [self.evaluate_async(seq) for seq in sequences]
        return list(await asyncio.gather(*tasks, return_exceptions=False))

    # ------------------------------------------------------------------
    # Utility
    # ------------------------------------------------------------------

    def __repr__(self) -> str:
        return f"{self.__class__.__name__}(name={self.name!r})"

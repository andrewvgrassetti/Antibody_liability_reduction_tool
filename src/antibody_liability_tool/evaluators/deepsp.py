"""DeepSP evaluator for antibody spatial property prediction.

Wraps the DeepSP package (optional dependency) to produce 30 spatial
property descriptors and derive aggregation/developability metrics.
When DeepSP is not installed a lightweight mock implementation is used so
the rest of the pipeline can still run.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any, Sequence

from antibody_liability_tool.caching import cached
from antibody_liability_tool.evaluators.base import (
    BaseEvaluator,
    ComparisonResult,
    EvaluationResult,
)

logger = logging.getLogger(__name__)

# ──────────────────────────────────────────────────────────────────────
# Optional import
# ──────────────────────────────────────────────────────────────────────

_DEEPSP_AVAILABLE = False
try:
    import DeepSP as _deepsp_mod  # type: ignore[import-untyped]

    _DEEPSP_AVAILABLE = True
    logger.debug("DeepSP package found.")
except ImportError:
    _deepsp_mod = None
    logger.debug("DeepSP not installed – mock mode will be used.")


# ──────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────

# The 30 canonical DeepSP descriptor names (grouped logically).
DEEPSP_DESCRIPTORS: list[str] = [
    # Hydrophobicity
    "hyd_D1",
    "hyd_D2",
    "hyd_D3",
    "hyd_D4",
    "hyd_D5",
    "hyd_D6",
    # Charge
    "chg_D1",
    "chg_D2",
    "chg_D3",
    "chg_D4",
    "chg_D5",
    "chg_D6",
    # Flexibility
    "flx_D1",
    "flx_D2",
    "flx_D3",
    "flx_D4",
    "flx_D5",
    "flx_D6",
    # Surface accessibility
    "sac_D1",
    "sac_D2",
    "sac_D3",
    "sac_D4",
    "sac_D5",
    "sac_D6",
    # Aggregation-related
    "agg_D1",
    "agg_D2",
    "agg_D3",
    "agg_D4",
    "agg_D5",
    "agg_D6",
]

# Descriptors for which *lower* is better (aggregation-prone).
AGGREGATION_DESCRIPTORS: set[str] = {
    "agg_D1",
    "agg_D2",
    "agg_D3",
    "agg_D4",
    "agg_D5",
    "agg_D6",
    "hyd_D1",
    "hyd_D2",
    "hyd_D3",
    "hyd_D4",
    "hyd_D5",
    "hyd_D6",
}

# Focus metrics used for quick summarization.
FOCUS_METRICS: list[str] = ["SAP_score", "SCM_score", "developability_index"]


@dataclass
class DeepSPResult:
    """Container for DeepSP prediction output.

    Attributes
    ----------
    descriptors : dict[str, float]
        All 30 spatial property descriptors.
    SAP_score : float
        Spatial Aggregation Propensity score.
    SCM_score : float
        Spatial Charge Map score.
    developability_index : float
        Combined developability metric.
    is_mock : bool
        ``True`` when results come from the fallback mock.
    """

    descriptors: dict[str, float] = field(default_factory=dict)
    SAP_score: float = 0.0
    SCM_score: float = 0.0
    developability_index: float = 0.0
    is_mock: bool = False

    def to_metrics(self) -> dict[str, float]:
        metrics: dict[str, float] = dict(self.descriptors)
        metrics["SAP_score"] = self.SAP_score
        metrics["SCM_score"] = self.SCM_score
        metrics["developability_index"] = self.developability_index
        return metrics


# ──────────────────────────────────────────────────────────────────────
# Mock prediction (used when DeepSP is unavailable)
# ──────────────────────────────────────────────────────────────────────


def _mock_predict(sequence: str) -> DeepSPResult:
    """Deterministic mock that derives pseudo-scores from the sequence.

    Values are length-normalised so that tests remain reproducible.
    """
    n = len(sequence)
    descriptors = {
        name: round(((i + 1) * n) % 100 / 100.0, 4) for i, name in enumerate(DEEPSP_DESCRIPTORS)
    }
    sap = round(sum(descriptors.get(d, 0) for d in list(AGGREGATION_DESCRIPTORS)[:6]) / 6, 4)
    scm = round(sum(descriptors.get(f"chg_D{i}", 0) for i in range(1, 7)) / 6, 4)
    dev_idx = round((sap + scm) / 2, 4)
    return DeepSPResult(
        descriptors=descriptors,
        SAP_score=sap,
        SCM_score=scm,
        developability_index=dev_idx,
        is_mock=True,
    )


# ──────────────────────────────────────────────────────────────────────
# Real prediction wrapper
# ──────────────────────────────────────────────────────────────────────


@cached(ttl=604800)
def _predict_deepsp(sequence: str) -> DeepSPResult:
    """Run DeepSP on *sequence* and return structured results (cached)."""
    if not _DEEPSP_AVAILABLE or _deepsp_mod is None:
        return _mock_predict(sequence)

    raw_output = _deepsp_mod.predict(sequence)

    # DeepSP may return a dict or a list of dicts depending on version.
    if isinstance(raw_output, list):
        raw_output = raw_output[0]

    descriptors: dict[str, float] = {}
    for name in DEEPSP_DESCRIPTORS:
        val = raw_output.get(name, 0.0)
        try:
            descriptors[name] = float(val)
        except (TypeError, ValueError):
            descriptors[name] = 0.0

    sap = float(raw_output.get("SAP_score", 0.0))
    scm = float(raw_output.get("SCM_score", 0.0))
    dev_idx = float(raw_output.get("developability_index", (sap + scm) / 2))

    return DeepSPResult(
        descriptors=descriptors,
        SAP_score=sap,
        SCM_score=scm,
        developability_index=dev_idx,
        is_mock=False,
    )


# ──────────────────────────────────────────────────────────────────────
# Evaluator
# ──────────────────────────────────────────────────────────────────────


class DeepSPEvaluator(BaseEvaluator):
    """Evaluate antibody sequences using DeepSP spatial property predictions.

    Parameters
    ----------
    config : dict[str, Any] | None
        ``deepsp`` section from the project config.
    """

    def __init__(self, config: dict[str, Any] | None = None) -> None:
        super().__init__(config=config, name="DeepSP")
        self.max_stability_delta: float = float(self.config.get("max_stability_delta", 0.1))
        self.focus_metrics: list[str] = list(self.config.get("focus_metrics", FOCUS_METRICS))

    @property
    def is_available(self) -> bool:
        """Return whether the real DeepSP package is importable."""
        return _DEEPSP_AVAILABLE

    # ------------------------------------------------------------------
    # Core prediction
    # ------------------------------------------------------------------

    def run(self, sequence: str) -> DeepSPResult:
        """Run DeepSP on a single sequence.

        Returns
        -------
        DeepSPResult
            Contains all 30 descriptors plus focus metrics.
        """
        return _predict_deepsp(sequence)

    # ------------------------------------------------------------------
    # BaseEvaluator interface
    # ------------------------------------------------------------------

    def evaluate(self, sequence: str) -> EvaluationResult:
        try:
            dsp = self.run(sequence)
            return EvaluationResult(
                sequence=sequence,
                evaluator_name=self.name,
                metrics=dsp.to_metrics(),
                raw={"deepsp_result": dsp},
                success=True,
            )
        except Exception as exc:  # noqa: BLE001
            self._logger.error("DeepSP evaluation failed: %s", exc)
            return EvaluationResult(
                sequence=sequence,
                evaluator_name=self.name,
                success=False,
                error_message=str(exc),
            )

    # ------------------------------------------------------------------
    # Comparison helpers
    # ------------------------------------------------------------------

    def compare_to_parent(
        self,
        parent_result: EvaluationResult,
        mutant_result: EvaluationResult,
    ) -> ComparisonResult:
        """Compare mutant vs parent on stability / aggregation metrics.

        A mutant *worsens* stability when any aggregation-prone descriptor
        increases beyond ``max_stability_delta``.
        """
        comp = self.compare(
            parent_result,
            mutant_result,
            lower_is_better=AGGREGATION_DESCRIPTORS | {"SAP_score", "developability_index"},
            higher_is_better={"SCM_score"},
        )
        # Apply delta tolerance
        tolerance_worsened = [
            m
            for m in comp.worsened_metrics
            if abs(comp.deltas.get(m, 0.0)) > self.max_stability_delta
        ]
        comp.worsened_metrics = tolerance_worsened
        comp.passes_filter = len(tolerance_worsened) == 0
        return comp

    def filter_mutants(
        self,
        parent_sequence: str,
        mutant_sequences: Sequence[str],
    ) -> list[str]:
        """Return only mutants whose aggregation descriptors ≤ parent values.

        Parameters
        ----------
        parent_sequence : str
            Wild-type sequence.
        mutant_sequences : Sequence[str]
            Candidate mutant sequences.

        Returns
        -------
        list[str]
            Sequences that pass the aggregation filter.
        """
        parent_result = self.evaluate(parent_sequence)
        if not parent_result.success:
            self._logger.warning("Parent evaluation failed; returning all mutants unfiltered.")
            return list(mutant_sequences)

        accepted: list[str] = []
        for seq in mutant_sequences:
            mutant_result = self.evaluate(seq)
            if not mutant_result.success:
                self._logger.warning("Skipping failed mutant %.20s…", seq)
                continue
            comp = self.compare_to_parent(parent_result, mutant_result)
            if comp.passes_filter:
                accepted.append(seq)
            else:
                self._logger.debug("Filtered out mutant (worsened: %s)", comp.worsened_metrics)
        self._logger.info(
            "DeepSP filter: %d / %d mutants accepted",
            len(accepted),
            len(mutant_sequences),
        )
        return accepted

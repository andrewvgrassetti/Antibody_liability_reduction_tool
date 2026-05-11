"""Composite scoring and ranking of mutation candidates.

Combines multiple evaluation metrics (PSH, PPC, OASis, stability) into a
single weighted score for ranking and selection.
"""

from __future__ import annotations

import logging
from dataclasses import dataclass, field
from typing import Any

logger = logging.getLogger(__name__)


@dataclass
class ScoredCandidate:
    """A mutation candidate annotated with its composite score and metrics.

    Attributes
    ----------
    label : str
        Human-readable mutation label (e.g. ``"V56S+K64Q"``).
    sequence : str
        Mutated amino-acid sequence.
    parent_metrics : dict[str, float]
        Evaluation metrics for the parent (wild-type) sequence.
    mutant_metrics : dict[str, float]
        Evaluation metrics for this mutant.
    composite_score : float
        Weighted composite score (higher is better).
    deltas : dict[str, float]
        Per-metric delta (mutant − parent).
    metadata : dict[str, Any]
        Arbitrary additional data (mutation combination, etc.).
    """

    label: str = ""
    sequence: str = ""
    parent_metrics: dict[str, float] = field(default_factory=dict)
    mutant_metrics: dict[str, float] = field(default_factory=dict)
    composite_score: float = 0.0
    deltas: dict[str, float] = field(default_factory=dict)
    metadata: dict[str, Any] = field(default_factory=dict)


class CompositeScorer:
    """Configurable composite scorer for mutation candidates.

    Weights control the relative importance of each metric improvement.
    All metric contributions are normalised to ``[0, 1]`` before weighting.

    Parameters
    ----------
    psh_reduction : float
        Weight for PSH (Patches of Surface Hydrophobicity) reduction.
    ppc_reduction : float
        Weight for PPC (Patches of Positive Charge) reduction.
    oasis_delta : float
        Weight for OASis humanness improvement.
    stability_penalty : float
        Weight for stability penalty (penalises worsened aggregation).
    oasis_threshold : float
        Minimum OASis score a candidate must achieve.
    max_stability_delta : float
        Maximum acceptable worsening in stability metric.
    """

    def __init__(
        self,
        psh_reduction: float = 0.30,
        ppc_reduction: float = 0.25,
        oasis_delta: float = 0.25,
        stability_penalty: float = 0.20,
        oasis_threshold: float = 0.0,
        max_stability_delta: float = 0.5,
    ) -> None:
        self.weights = {
            "psh_reduction": psh_reduction,
            "ppc_reduction": ppc_reduction,
            "oasis_delta": oasis_delta,
            "stability_penalty": stability_penalty,
        }
        self.oasis_threshold = oasis_threshold
        self.max_stability_delta = max_stability_delta

        total = sum(self.weights.values())
        if total > 0:
            self.weights = {k: v / total for k, v in self.weights.items()}

        logger.debug("CompositeScorer initialised with weights: %s", self.weights)

    @classmethod
    def from_config(cls, config: dict[str, Any]) -> CompositeScorer:
        """Create a scorer from the project configuration dict.

        Parameters
        ----------
        config : dict
            Full project config (expects ``scoring.weights`` and
            ``scoring.oasis_threshold`` keys).

        Returns
        -------
        CompositeScorer
        """
        scfg = config.get("scoring", {})
        weights = scfg.get("weights", {})
        return cls(
            psh_reduction=float(weights.get("psh_reduction", 0.30)),
            ppc_reduction=float(weights.get("ppc_reduction", 0.25)),
            oasis_delta=float(weights.get("oasis_delta", 0.25)),
            stability_penalty=float(weights.get("stability_penalty", 0.20)),
            oasis_threshold=float(scfg.get("oasis_threshold", 0.0)),
            max_stability_delta=float(
                config.get("deepsp", {}).get("max_stability_delta", 0.5)
            ),
        )

    # ------------------------------------------------------------------
    # Scoring
    # ------------------------------------------------------------------

    def score(
        self,
        parent_metrics: dict[str, float],
        mutant_metrics: dict[str, float],
    ) -> float:
        """Compute a composite score for a single mutant relative to its parent.

        Parameters
        ----------
        parent_metrics : dict[str, float]
            Metrics from evaluating the parent sequence.
        mutant_metrics : dict[str, float]
            Metrics from evaluating the mutant sequence.

        Returns
        -------
        float
            Composite score in approximately ``[0, 1]`` (higher is better).
        """
        # PSH reduction: positive if mutant PSH < parent PSH
        psh_parent = parent_metrics.get("PSH", 0.0)
        psh_mutant = mutant_metrics.get("PSH", 0.0)
        psh_reduction = self._safe_reduction(psh_parent, psh_mutant)

        # PPC reduction
        ppc_parent = parent_metrics.get("PPC", 0.0)
        ppc_mutant = mutant_metrics.get("PPC", 0.0)
        ppc_reduction = self._safe_reduction(ppc_parent, ppc_mutant)

        # OASis delta: improvement = mutant − parent (higher is better)
        oasis_parent = parent_metrics.get("oasis_humanness", 0.0)
        oasis_mutant = mutant_metrics.get("oasis_humanness", 0.0)
        oasis_delta = oasis_mutant - oasis_parent
        # Clip to [−1, 1] and shift to [0, 1]
        oasis_contrib = max(0.0, min(1.0, (oasis_delta + 1.0) / 2.0))

        # Stability penalty: penalise SAP_score increase
        sap_parent = parent_metrics.get("SAP_score", 0.0)
        sap_mutant = mutant_metrics.get("SAP_score", 0.0)
        sap_delta = sap_mutant - sap_parent
        # Convert to penalty: 1.0 = no worsening, 0.0 = worst
        divisor = max(self.max_stability_delta, 1e-6)
        stability_contrib = max(0.0, 1.0 - max(0.0, sap_delta) / divisor)

        composite = (
            self.weights["psh_reduction"] * psh_reduction
            + self.weights["ppc_reduction"] * ppc_reduction
            + self.weights["oasis_delta"] * oasis_contrib
            + self.weights["stability_penalty"] * stability_contrib
        )

        return round(composite, 6)

    # ------------------------------------------------------------------
    # Ranking and filtering
    # ------------------------------------------------------------------

    def rank(
        self,
        candidates: list[dict[str, Any]],
    ) -> list[ScoredCandidate]:
        """Score and rank a list of candidates.

        Parameters
        ----------
        candidates : list[dict]
            Each dict must contain at minimum:
            - ``"label"`` : str
            - ``"sequence"`` : str
            - ``"parent_metrics"`` : dict[str, float]
            - ``"mutant_metrics"`` : dict[str, float]
            Additional keys are preserved in ``metadata``.

        Returns
        -------
        list[ScoredCandidate]
            Candidates sorted by composite score (descending).
        """
        scored: list[ScoredCandidate] = []

        for cand in candidates:
            parent_m = cand.get("parent_metrics", {})
            mutant_m = cand.get("mutant_metrics", {})

            composite = self.score(parent_m, mutant_m)

            deltas: dict[str, float] = {}
            for key in set(parent_m) | set(mutant_m):
                deltas[key] = mutant_m.get(key, 0.0) - parent_m.get(key, 0.0)

            scored.append(
                ScoredCandidate(
                    label=cand.get("label", ""),
                    sequence=cand.get("sequence", ""),
                    parent_metrics=parent_m,
                    mutant_metrics=mutant_m,
                    composite_score=composite,
                    deltas=deltas,
                    metadata={
                        k: v
                        for k, v in cand.items()
                        if k not in {"label", "sequence", "parent_metrics", "mutant_metrics"}
                    },
                )
            )

        scored.sort(key=lambda s: -s.composite_score)
        logger.info("Ranked %d candidates", len(scored))
        return scored

    def filter_candidates(
        self,
        scored: list[ScoredCandidate],
        min_oasis: float | None = None,
        max_stability_delta: float | None = None,
    ) -> list[ScoredCandidate]:
        """Filter scored candidates by constraint thresholds.

        Parameters
        ----------
        scored : list[ScoredCandidate]
            Pre-scored candidates.
        min_oasis : float, optional
            Minimum OASis humanness score. Defaults to ``self.oasis_threshold``.
        max_stability_delta : float, optional
            Maximum acceptable SAP_score increase. Defaults to ``self.max_stability_delta``.

        Returns
        -------
        list[ScoredCandidate]
            Candidates passing all constraints.
        """
        oasis_min = min_oasis if min_oasis is not None else self.oasis_threshold
        stab_max = (
            max_stability_delta if max_stability_delta is not None
            else self.max_stability_delta
        )

        accepted: list[ScoredCandidate] = []
        for cand in scored:
            oasis_score = cand.mutant_metrics.get("oasis_humanness", 1.0)
            if oasis_score < oasis_min:
                logger.debug(
                    "Filtered %s: OASis %.4f < threshold %.4f",
                    cand.label,
                    oasis_score,
                    oasis_min,
                )
                continue

            sap_delta = cand.deltas.get("SAP_score", 0.0)
            if sap_delta > stab_max:
                logger.debug(
                    "Filtered %s: SAP delta %.4f > max %.4f",
                    cand.label,
                    sap_delta,
                    stab_max,
                )
                continue

            accepted.append(cand)

        logger.info(
            "Filter: %d / %d candidates passed constraints",
            len(accepted),
            len(scored),
        )
        return accepted

    # ------------------------------------------------------------------
    # Private helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _safe_reduction(parent_val: float, mutant_val: float) -> float:
        """Compute fractional reduction clipped to [0, 1].

        Returns 1.0 for complete elimination and 0.0 for no improvement or worse.
        """
        if parent_val <= 0.0:
            return 0.5  # neutral if parent has no liability
        reduction = (parent_val - mutant_val) / parent_val
        return max(0.0, min(1.0, reduction))

"""Tests for the composite scorer module."""

from __future__ import annotations

import pytest

from antibody_liability_tool.optimization.scorer import CompositeScorer, ScoredCandidate


class TestCompositeScorer:
    """Tests for CompositeScorer."""

    def test_score_perfect_improvement(self) -> None:
        """Mutant eliminates all PSH/PPC liabilities and improves OASis."""
        scorer = CompositeScorer()
        parent = {"PSH": 1.0, "PPC": 1.0, "oasis_humanness": 0.5, "SAP_score": 0.3}
        mutant = {"PSH": 0.0, "PPC": 0.0, "oasis_humanness": 0.9, "SAP_score": 0.3}
        score = scorer.score(parent, mutant)
        assert score > 0.7

    def test_score_no_change(self) -> None:
        """Identical metrics should produce a moderate score (no improvement)."""
        scorer = CompositeScorer()
        metrics = {"PSH": 0.5, "PPC": 0.5, "oasis_humanness": 0.7, "SAP_score": 0.3}
        score = scorer.score(metrics, dict(metrics))
        # No PSH/PPC reduction → 0, oasis delta 0 → 0.5 contrib, stability neutral → 1.0
        assert 0.0 <= score <= 1.0

    def test_score_worsened_stability(self) -> None:
        """Worsened SAP_score should reduce the composite score."""
        scorer = CompositeScorer()
        parent = {"PSH": 1.0, "PPC": 1.0, "oasis_humanness": 0.5, "SAP_score": 0.1}
        mutant_good = {"PSH": 0.5, "PPC": 0.5, "oasis_humanness": 0.5, "SAP_score": 0.1}
        mutant_bad = {"PSH": 0.5, "PPC": 0.5, "oasis_humanness": 0.5, "SAP_score": 0.6}
        score_good = scorer.score(parent, mutant_good)
        score_bad = scorer.score(parent, mutant_bad)
        assert score_good > score_bad

    def test_rank_orders_by_score(self) -> None:
        scorer = CompositeScorer()
        parent = {"PSH": 1.0, "PPC": 1.0, "oasis_humanness": 0.5, "SAP_score": 0.3}
        candidates = [
            {
                "label": "bad",
                "sequence": "AAA",
                "parent_metrics": parent,
                "mutant_metrics": {"PSH": 1.0, "PPC": 1.0, "oasis_humanness": 0.5, "SAP_score": 0.3},
            },
            {
                "label": "good",
                "sequence": "BBB",
                "parent_metrics": parent,
                "mutant_metrics": {"PSH": 0.0, "PPC": 0.0, "oasis_humanness": 0.9, "SAP_score": 0.1},
            },
        ]
        ranked = scorer.rank(candidates)
        assert ranked[0].label == "good"
        assert ranked[1].label == "bad"
        assert ranked[0].composite_score > ranked[1].composite_score

    def test_filter_by_oasis_threshold(self) -> None:
        scorer = CompositeScorer(oasis_threshold=0.6)
        parent = {"PSH": 1.0, "PPC": 1.0, "oasis_humanness": 0.7, "SAP_score": 0.3}
        candidates = [
            {
                "label": "passes",
                "sequence": "AAA",
                "parent_metrics": parent,
                "mutant_metrics": {"PSH": 0.5, "PPC": 0.5, "oasis_humanness": 0.7, "SAP_score": 0.3},
            },
            {
                "label": "fails",
                "sequence": "BBB",
                "parent_metrics": parent,
                "mutant_metrics": {"PSH": 0.5, "PPC": 0.5, "oasis_humanness": 0.4, "SAP_score": 0.3},
            },
        ]
        scored = scorer.rank(candidates)
        filtered = scorer.filter_candidates(scored, min_oasis=0.6)
        labels = [c.label for c in filtered]
        assert "passes" in labels
        assert "fails" not in labels

    def test_custom_weights(self) -> None:
        """Custom weights should change relative scoring."""
        scorer_psh = CompositeScorer(psh_reduction=1.0, ppc_reduction=0.0, oasis_delta=0.0, stability_penalty=0.0)
        scorer_ppc = CompositeScorer(psh_reduction=0.0, ppc_reduction=1.0, oasis_delta=0.0, stability_penalty=0.0)
        parent = {"PSH": 1.0, "PPC": 0.0, "oasis_humanness": 0.5, "SAP_score": 0.3}
        mutant = {"PSH": 0.0, "PPC": 0.0, "oasis_humanness": 0.5, "SAP_score": 0.3}
        score_psh = scorer_psh.score(parent, mutant)
        score_ppc = scorer_ppc.score(parent, mutant)
        # PSH-focused scorer should give high score (PSH reduced from 1->0)
        assert score_psh > 0.9
        # PPC-focused scorer: PPC was 0->0, neutral
        assert score_ppc < score_psh

    def test_from_config(self) -> None:
        config = {
            "scoring": {
                "weights": {
                    "psh_reduction": 0.4,
                    "ppc_reduction": 0.3,
                    "oasis_delta": 0.2,
                    "stability_penalty": 0.1,
                },
                "oasis_threshold": 0.5,
            },
            "deepsp": {"max_stability_delta": 0.3},
        }
        scorer = CompositeScorer.from_config(config)
        assert scorer.oasis_threshold == 0.5
        assert scorer.max_stability_delta == 0.3

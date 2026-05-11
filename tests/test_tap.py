"""Tests for the TAP evaluator module."""

from __future__ import annotations

from pathlib import Path

from antibody_liability_tool.evaluators.tap import (
    TAPEvaluator,
    TAPResult,
    _parse_tap_csv,
    _parse_tap_html,
)


class TestTAPResult:
    """Tests for the TAPResult dataclass."""

    def test_tap_result_dataclass(self) -> None:
        result = TAPResult(
            patches_positive_charge=2,
            patches_negative_charge=1,
            patches_hydrophobic=3,
            PSH=0.45,
            PPC=0.30,
            PNC=0.15,
            SFvCSP=0.22,
        )
        assert result.patches_positive_charge == 2
        assert result.PSH == 0.45
        assert result.SFvCSP == 0.22

    def test_tap_result_to_metrics(self) -> None:
        result = TAPResult(PSH=0.5, PPC=0.3, PNC=0.1, SFvCSP=0.2)
        metrics = result.to_metrics()
        assert metrics["PSH"] == 0.5
        assert metrics["PPC"] == 0.3
        assert "patches_hydrophobic" in metrics


class TestTAPEvaluator:
    """Tests for TAPEvaluator methods."""

    def test_manual_mode_write_batch(self, tmp_path: Path) -> None:
        """write_batch() should create a FASTA file."""
        evaluator = TAPEvaluator(config={"mode": "manual", "batch_dir": str(tmp_path)})
        sequences = {"parent": "EVQLVES", "mutant1": "EVQLVQS"}
        fasta_path = evaluator.write_batch(sequences)
        assert fasta_path.exists()
        content = fasta_path.read_text()
        assert ">parent" in content
        assert ">mutant1" in content
        assert "EVQLVES" in content

    def test_parse_tap_results_csv(self) -> None:
        csv_text = (
            "PSH,PPC,PNC,SFvCSP,patches_positive_charge,patches_negative_charge,patches_hydrophobic\n"
            "0.45,0.30,0.15,0.22,2,1,3\n"
            "0.50,0.35,0.20,0.25,3,2,4\n"
        )
        results = _parse_tap_csv(csv_text)
        assert len(results) == 2
        assert results[0].PSH == 0.45
        assert results[0].patches_positive_charge == 2
        assert results[1].PSH == 0.50

    def test_parse_tap_html(self) -> None:
        html = """
        <html><body>
        PSH = 0.45
        PPC = 0.30
        PNC = 0.15
        SFvCSP = 0.22
        positive charge patches = 2
        negative charge patches = 1
        hydrophobic patches = 3
        </body></html>
        """
        result = _parse_tap_html(html)
        assert result.PSH == 0.45
        assert result.PPC == 0.30
        assert result.patches_positive_charge == 2

    def test_compare_parent_mutant(self) -> None:
        """TAP comparison should detect when mutant is better."""
        evaluator = TAPEvaluator(config={"mode": "manual"})

        from antibody_liability_tool.evaluators.base import EvaluationResult

        parent_result = EvaluationResult(
            sequence="PARENT",
            evaluator_name="TAP",
            metrics={"PSH": 0.5, "PPC": 0.4, "PNC": 0.3, "SFvCSP": 0.2},
            success=True,
        )
        mutant_result = EvaluationResult(
            sequence="MUTANT",
            evaluator_name="TAP",
            metrics={"PSH": 0.3, "PPC": 0.2, "PNC": 0.2, "SFvCSP": 0.1},
            success=True,
        )
        comp = evaluator.compare_to_parent(parent_result, mutant_result)
        assert comp.passes_filter is True
        assert len(comp.improved_metrics) > 0

    def test_manual_mode_evaluate_returns_success(self) -> None:
        evaluator = TAPEvaluator(config={"mode": "manual"})
        result = evaluator.evaluate("EVQLVES")
        assert result.success is True
        assert result.raw.get("mode") == "manual"

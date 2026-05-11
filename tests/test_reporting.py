"""Tests for the reporting / FASTA export module."""

from __future__ import annotations

from pathlib import Path

import pytest

from antibody_liability_tool.optimization.scorer import ScoredCandidate
from antibody_liability_tool.reporting.fasta_export import export_fasta
from tests.conftest import SAMPLE_VH


class TestFastaExport:
    """Tests for export_fasta()."""

    def test_fasta_export_creates_file(self, tmp_path: Path) -> None:
        output = tmp_path / "test_output.fasta"
        candidates = [
            ScoredCandidate(
                label="V56S",
                sequence="EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWSARIYPTNGYTRYADSVKGRFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS",
                parent_metrics={"PSH": 0.5},
                mutant_metrics={"PSH": 0.3},
                composite_score=0.85,
            ),
        ]
        result_path = export_fasta(SAMPLE_VH, candidates, output)
        assert result_path.exists()
        content = result_path.read_text()
        assert len(content) > 0

    def test_fasta_export_contains_parent(self, tmp_path: Path) -> None:
        output = tmp_path / "test_output.fasta"
        candidates = [
            ScoredCandidate(
                label="V56S",
                sequence="MUTANTSEQ",
                composite_score=0.75,
            ),
        ]
        export_fasta(SAMPLE_VH, candidates, output, parent_id="trastuzumab_VH")
        content = output.read_text()
        assert ">trastuzumab_VH" in content
        assert SAMPLE_VH[:40] in content

    def test_fasta_export_limits_top_n(self, tmp_path: Path) -> None:
        output = tmp_path / "test_output.fasta"
        candidates = [
            ScoredCandidate(label=f"mut{i}", sequence=f"SEQ{i}", composite_score=float(i))
            for i in range(20)
        ]
        export_fasta(SAMPLE_VH, candidates, output, top_n=5)
        content = output.read_text()
        # Should have parent + 5 candidates = 6 headers
        headers = [line for line in content.splitlines() if line.startswith(">")]
        assert len(headers) == 6

    def test_fasta_export_skips_empty_sequence(self, tmp_path: Path) -> None:
        output = tmp_path / "test_output.fasta"
        candidates = [
            ScoredCandidate(label="empty", sequence="", composite_score=0.5),
            ScoredCandidate(label="valid", sequence="ACDEFG", composite_score=0.5),
        ]
        export_fasta(SAMPLE_VH, candidates, output)
        content = output.read_text()
        assert "empty" not in content
        assert "valid" in content

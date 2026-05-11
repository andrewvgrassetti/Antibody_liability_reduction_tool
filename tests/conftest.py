"""Shared fixtures for the Antibody Liability Reduction Tool test suite."""

from __future__ import annotations

import os
import shutil
from pathlib import Path
from unittest.mock import patch

import pytest

# Trastuzumab VH sequence (well-characterised therapeutic antibody).
SAMPLE_VH = (
    "EVQLVESGGGLVQPGGSLRLSCAASGFNIKDTYIHWVRQAPGKGLEWVARIYPTNGYTRYADSVKG"
    "RFTISADTSKNTAYLQMNSLRAEDTAVYYCSRWGGDGFYAMDYWGQGTLVTVSS"
)


@pytest.fixture()
def sample_numbered_sequence() -> dict[str, str]:
    """Return a numbered Trastuzumab VH sequence using fallback numbering."""
    from antibody_liability_tool.numbering.imgt import number_sequence

    with patch(
        "antibody_liability_tool.numbering.imgt._number_with_anarci",
        return_value=None,
    ):
        return number_sequence(SAMPLE_VH)


@pytest.fixture()
def sample_config() -> dict[str, str]:
    """Return a minimal test configuration dict."""
    return {
        "liabilities": {
            "hydrophobic_residues": ["F", "I", "L", "M", "V", "W", "Y"],
            "positive_charge_residues": ["R", "K"],
            "positive_charge_cluster_size": 2,
            "charge_cluster_distance": 5,
        },
        "mutations": {
            "min_human_frequency": 0.05,
            "avoid_motifs": {
                "free_cys": True,
                "oxidation_met": True,
            },
        },
        "cache": {
            "enabled": True,
            "directory": ".cache/test_cache",
            "ttl_seconds": 60,
        },
        "output": {
            "directory": "output",
            "formats": ["fasta"],
            "top_n_report": 5,
        },
    }


@pytest.fixture()
def tmp_cache_dir(tmp_path: Path) -> Path:
    """Return a temporary cache directory that is cleaned up after the test."""
    cache_dir = tmp_path / "test_cache"
    cache_dir.mkdir()
    yield cache_dir
    if cache_dir.exists():
        shutil.rmtree(cache_dir)

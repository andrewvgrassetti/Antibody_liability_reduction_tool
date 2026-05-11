"""Configuration loading and validation for the Antibody Liability Reduction Tool."""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any

import yaml

logger = logging.getLogger(__name__)

# Minimal defaults that guarantee every required key exists.
_DEFAULTS: dict[str, Any] = {
    "input": {
        "numbering_scheme": "imgt",
        "chain_type": "H",
    },
    "surface": {
        "method": "rule_based",
        "sasa_threshold_exposed": 40.0,
        "sasa_threshold_partial": 20.0,
    },
    "liabilities": {
        "hydrophobic_residues": ["F", "I", "L", "M", "V", "W", "Y"],
        "positive_charge_residues": ["R", "K"],
        "positive_charge_cluster_size": 2,
        "charge_cluster_distance": 5,
    },
    "mutations": {
        "min_human_frequency": 0.05,
        "avoid_motifs": {
            "n_glycosylation": "N[^P][ST]",
            "deamidation_ng": "NG",
            "deamidation_ns": "NS",
            "isomerization_dg": "DG",
            "isomerization_ds": "DS",
            "oxidation_met": True,
            "free_cys": True,
        },
    },
    "cache": {
        "enabled": True,
        "directory": ".cache/antibody_tool",
        "ttl_seconds": 604800,
    },
    "logging": {
        "level": "INFO",
        "format": "%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    },
    "output": {
        "directory": "output",
        "formats": ["html", "fasta", "csv"],
        "top_n_report": 10,
    },
}

_REQUIRED_SECTIONS = ("liabilities", "mutations")


def _deep_merge(base: dict[str, Any], override: dict[str, Any]) -> dict[str, Any]:
    """Recursively merge *override* into *base*, returning a new dict."""
    merged: dict[str, Any] = dict(base)
    for key, value in override.items():
        if key in merged and isinstance(merged[key], dict) and isinstance(value, dict):
            merged[key] = _deep_merge(merged[key], value)
        else:
            merged[key] = value
    return merged


def _validate(cfg: dict[str, Any]) -> None:
    """Raise ``ValueError`` if required sections/fields are missing."""
    for section in _REQUIRED_SECTIONS:
        if section not in cfg:
            raise ValueError(f"Missing required config section: '{section}'")
    liabilities = cfg["liabilities"]
    if not liabilities.get("hydrophobic_residues"):
        raise ValueError("liabilities.hydrophobic_residues must be a non-empty list")
    if not liabilities.get("positive_charge_residues"):
        raise ValueError("liabilities.positive_charge_residues must be a non-empty list")
    mutations = cfg["mutations"]
    freq = mutations.get("min_human_frequency")
    if freq is None or freq <= 0.0 or freq > 1.0:
        raise ValueError(
            "mutations.min_human_frequency must be between 0 and 1 (exclusive/inclusive)"
        )


def load_config(path: str | Path | None = None) -> dict[str, Any]:
    """Load a YAML configuration file and merge it with built-in defaults.

    Parameters
    ----------
    path:
        Path to a YAML config file. If *None* or the file does not exist,
        the built-in defaults are returned.

    Returns
    -------
    dict[str, Any]
        Merged configuration dictionary.
    """
    if path is not None:
        config_path = Path(path)
        if config_path.exists():
            logger.debug("Loading config from %s", config_path)
            with config_path.open() as fh:
                user_cfg = yaml.safe_load(fh) or {}
            merged = _deep_merge(_DEFAULTS, user_cfg)
        else:
            logger.warning("Config file %s not found – using defaults", config_path)
            merged = dict(_DEFAULTS)
    else:
        merged = dict(_DEFAULTS)

    _validate(merged)
    return merged

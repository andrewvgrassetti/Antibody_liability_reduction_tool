"""Tests for the configuration module."""

from __future__ import annotations

from pathlib import Path

import pytest
import yaml

from antibody_liability_tool.config import _DEFAULTS, _deep_merge, load_config


class TestLoadConfig:
    """Tests for load_config()."""

    def test_load_default_config(self) -> None:
        """Loading with no path should return defaults."""
        cfg = load_config(None)
        assert "liabilities" in cfg
        assert "mutations" in cfg
        assert cfg["mutations"]["min_human_frequency"] == 0.05

    def test_merge_configs(self, tmp_path: Path) -> None:
        """User config should override defaults while preserving unset keys."""
        user_cfg = {
            "mutations": {"min_human_frequency": 0.10},
            "custom_key": "custom_value",
        }
        cfg_file = tmp_path / "config.yaml"
        cfg_file.write_text(yaml.dump(user_cfg))

        cfg = load_config(cfg_file)
        assert cfg["mutations"]["min_human_frequency"] == 0.10
        assert cfg["custom_key"] == "custom_value"
        # Defaults should still be present
        assert "liabilities" in cfg
        assert cfg["liabilities"]["hydrophobic_residues"] == ["F", "I", "L", "M", "V", "W", "Y"]

    def test_missing_config_file(self) -> None:
        """A nonexistent path should fall back to defaults."""
        cfg = load_config("/nonexistent/path/config.yaml")
        assert "liabilities" in cfg
        assert "mutations" in cfg


class TestDeepMerge:
    """Tests for _deep_merge()."""

    def test_deep_merge_nested(self) -> None:
        base = {"a": {"b": 1, "c": 2}, "d": 3}
        override = {"a": {"b": 10, "e": 5}}
        merged = _deep_merge(base, override)
        assert merged["a"]["b"] == 10
        assert merged["a"]["c"] == 2
        assert merged["a"]["e"] == 5
        assert merged["d"] == 3

    def test_deep_merge_does_not_mutate_base(self) -> None:
        base = {"a": {"b": 1}}
        override = {"a": {"b": 2}}
        _deep_merge(base, override)
        assert base["a"]["b"] == 1

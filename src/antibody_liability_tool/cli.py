"""Click-based command-line interface for the Antibody Liability Reduction Tool."""

from __future__ import annotations

import logging
import sys
from pathlib import Path

import click
import yaml

from antibody_liability_tool.config import load_config

logger = logging.getLogger(__name__)


def _setup_logging(level: str, fmt: str | None = None) -> None:
    """Configure root logging based on config values."""
    numeric_level = getattr(logging, level.upper(), logging.INFO)
    log_format = fmt or "%(asctime)s - %(name)s - %(levelname)s - %(message)s"
    logging.basicConfig(level=numeric_level, format=log_format, stream=sys.stderr)


def _read_sequence(sequence: str | None, fasta: str | None) -> str:
    """Return an amino-acid sequence from either a raw string or a FASTA file."""
    if sequence:
        return sequence.strip()
    if fasta:
        fasta_path = Path(fasta)
        if not fasta_path.exists():
            raise click.BadParameter(f"FASTA file not found: {fasta}")
        lines = fasta_path.read_text().splitlines()
        seq_lines = [line.strip() for line in lines if not line.startswith(">")]
        return "".join(seq_lines)
    raise click.UsageError("Provide --sequence or --fasta")


@click.group()
@click.version_option(package_name="antibody-liability-tool")
def main() -> None:
    """Antibody Liability Reduction Tool – identify and fix VH surface liabilities."""


@main.command()
@click.option("--sequence", "-s", default=None, help="Raw VH amino-acid sequence.")
@click.option("--fasta", "-f", default=None, help="Path to a FASTA file with the VH sequence.")
@click.option(
    "--config",
    "-c",
    "config_path",
    default="config/default.yaml",
    show_default=True,
    help="Path to YAML configuration file.",
)
@click.option(
    "--output-dir",
    "-o",
    default="output",
    show_default=True,
    help="Directory for results.",
)
def run(
    sequence: str | None,
    fasta: str | None,
    config_path: str,
    output_dir: str,
) -> None:
    """Run the full liability-reduction pipeline."""
    cfg = load_config(config_path)
    _setup_logging(
        cfg.get("logging", {}).get("level", "INFO"),
        cfg.get("logging", {}).get("format"),
    )
    logger.info("Loading configuration from %s", config_path)

    seq = _read_sequence(sequence, fasta)
    logger.info("Input sequence length: %d", len(seq))

    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)

    # Stage 1-3: Analyze liabilities
    from antibody_liability_tool.liabilities.detector import detect_liabilities
    from antibody_liability_tool.mutations.generator import generate_mutations
    from antibody_liability_tool.numbering.imgt import number_sequence
    from antibody_liability_tool.surface.exposure import classify_surface_exposure

    numbered = number_sequence(seq)
    exposure = classify_surface_exposure(numbered)
    liabilities = detect_liabilities(numbered, exposure, cfg)
    logger.info("Detected %d liabilities", len(liabilities))

    # Stage 4: Generate mutations
    mutations = generate_mutations(liabilities, numbered, cfg)
    logger.info("Generated %d candidate mutations", len(mutations))

    # Write summary
    summary_path = out / "summary.yaml"
    summary = {
        "sequence_length": len(seq),
        "liabilities_found": len(liabilities),
        "mutations_proposed": len(mutations),
        "liabilities": [
            {
                "position": li.position,
                "residue": li.residue,
                "imgt_number": li.imgt_number,
                "region": li.region,
                "reason": li.reason,
                "severity": li.severity,
            }
            for li in liabilities
        ],
        "mutations": [
            {
                "position": m.position,
                "original": m.original_aa,
                "proposed": m.proposed_aa,
                "frequency": m.human_frequency,
                "rationale": m.rationale,
            }
            for m in mutations
        ],
    }
    summary_path.write_text(yaml.dump(summary, default_flow_style=False, sort_keys=False))
    click.echo(f"Pipeline complete. Results written to {out}")


@main.command()
@click.option("--sequence", "-s", default=None, help="Raw VH amino-acid sequence.")
@click.option("--fasta", "-f", default=None, help="Path to a FASTA file with the VH sequence.")
@click.option(
    "--config",
    "-c",
    "config_path",
    default="config/default.yaml",
    show_default=True,
    help="Path to YAML configuration file.",
)
def analyze(
    sequence: str | None,
    fasta: str | None,
    config_path: str,
) -> None:
    """Analyze liabilities only (stages 1-3: number, expose, detect)."""
    cfg = load_config(config_path)
    _setup_logging(
        cfg.get("logging", {}).get("level", "INFO"),
        cfg.get("logging", {}).get("format"),
    )

    seq = _read_sequence(sequence, fasta)

    from antibody_liability_tool.liabilities.detector import detect_liabilities
    from antibody_liability_tool.numbering.imgt import number_sequence
    from antibody_liability_tool.surface.exposure import classify_surface_exposure

    numbered = number_sequence(seq)
    exposure = classify_surface_exposure(numbered)
    liabilities = detect_liabilities(numbered, exposure, cfg)

    click.echo(f"Detected {len(liabilities)} liabilities:")
    for li in liabilities:
        click.echo(f"  [{li.severity}] {li.imgt_number} {li.residue} ({li.region}): {li.reason}")


@main.command("refresh-abysis")
@click.option(
    "--output",
    "-o",
    default="data/abysis_vh_frequencies.json",
    show_default=True,
    help="Output path for the Abysis frequency cache.",
)
def refresh_abysis(output: str) -> None:
    """Refresh the Abysis human VH frequency cache (placeholder)."""
    click.echo(
        f"Abysis cache refresh is not yet automated. "
        f"Please update {output} manually from http://www.abysis.org"
    )

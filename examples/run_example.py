#!/usr/bin/env python3
"""Example script demonstrating the Antibody Liability Reduction Tool pipeline.

This script runs the analyze-only pipeline (stages 1-3) on the trastuzumab VH
sequence and prints the detected liabilities and candidate mutations.
"""

from __future__ import annotations

from pathlib import Path

from antibody_liability_tool.config import load_config
from antibody_liability_tool.liabilities.detector import detect_liabilities
from antibody_liability_tool.mutations.generator import generate_mutations
from antibody_liability_tool.numbering.imgt import number_sequence
from antibody_liability_tool.surface.exposure import classify_surface_exposure


def read_fasta(path: str | Path) -> str:
    """Read the first sequence from a FASTA file."""
    lines = Path(path).read_text().splitlines()
    return "".join(line.strip() for line in lines if not line.startswith(">"))


def main() -> None:
    # Load configuration
    config = load_config()

    # Read example sequence
    example_fasta = Path(__file__).parent / "example_input.fasta"
    sequence = read_fasta(example_fasta)
    print(f"Input sequence ({len(sequence)} aa):")
    print(f"  {sequence[:60]}...")

    # Stage 1: IMGT numbering
    print("\n--- Stage 1: IMGT Numbering ---")
    numbered = number_sequence(sequence)
    print(f"  Numbered {len(numbered)} positions")

    # Stage 2: Surface exposure classification
    print("\n--- Stage 2: Surface Exposure ---")
    exposure = classify_surface_exposure(numbered)
    surface_count = sum(
        1 for v in exposure.values() if v.get("exposure") == "surface"
    )
    print(f"  {surface_count} surface-exposed positions out of {len(exposure)}")

    # Stage 3: Liability detection
    print("\n--- Stage 3: Liability Detection ---")
    liabilities = detect_liabilities(numbered, exposure, config)
    print(f"  Found {len(liabilities)} liabilities:")
    for li in liabilities:
        print(f"    [{li.severity}] {li.imgt_number} {li.residue} ({li.region}): {li.reason}")

    # Stage 4: Mutation generation
    print("\n--- Stage 4: Mutation Generation ---")
    mutations = generate_mutations(liabilities, numbered, config)
    print(f"  Generated {len(mutations)} candidate mutations:")
    for m in mutations[:10]:
        print(
            f"    {m.position} {m.original_aa} -> {m.proposed_aa} "
            f"(freq={m.human_frequency:.2f}, {m.rationale})"
        )
    if len(mutations) > 10:
        print(f"    ... and {len(mutations) - 10} more")

    # --- Stages 5-10 require external services ---
    # To run the full pipeline including TAP, DeepSP, and OASis evaluation:
    #
    # from antibody_liability_tool.pipeline import LiabilityReductionPipeline
    # pipeline = LiabilityReductionPipeline(output_dir="output")
    # result = pipeline.run(sequence)
    # print(f"Ranked {len(result.ranked_candidates)} candidates")

    print("\nDone! To run the full pipeline (stages 5-10), see the CLI:")
    print("  antibody-liability-tool run --fasta examples/example_input.fasta")


if __name__ == "__main__":
    main()

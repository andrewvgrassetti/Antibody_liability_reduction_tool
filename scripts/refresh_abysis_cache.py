#!/usr/bin/env python3
"""Refresh the Abysis human VH germline frequency cache.

This script fetches amino-acid frequency data from the Abysis database
(http://www.abysis.org) for human VH sequences at each IMGT-numbered
position and writes the results to data/abysis_vh_frequencies.json.

IMPORTANT:
- Respect Abysis Terms of Service and rate limits.
- Add appropriate delays between requests (>=1 second recommended).
- This script is provided as a skeleton; the actual scraping endpoints
  may change and should be verified before use.

Usage:
    python scripts/refresh_abysis_cache.py [--output data/abysis_vh_frequencies.json]
"""

from __future__ import annotations

import argparse
import json
import logging
import time
from pathlib import Path

import requests

logger = logging.getLogger(__name__)

ABYSIS_BASE_URL = "http://www.abysis.org/abysis/sequence_input/key_annotation"

# IMGT positions for VH (1-128, with insertion positions)
VH_IMGT_POSITIONS = [str(i) for i in range(1, 129)]

REQUEST_DELAY_SECONDS = 1.5


def fetch_position_frequencies(position: str) -> dict[str, float]:
    """Fetch amino-acid frequencies for a single IMGT position.

    Parameters
    ----------
    position : str
        IMGT position number (e.g., "1", "27", "104").

    Returns
    -------
    dict[str, float]
        Mapping of amino acid single-letter codes to their observed
        frequency (0.0 – 1.0) at this position in human VH sequences.

    Raises
    ------
    NotImplementedError
        This is a skeleton implementation. The actual HTTP request and
        HTML parsing logic must be implemented based on the current
        Abysis web interface.
    """
    raise NotImplementedError(
        "Abysis scraping is not yet implemented. "
        "Please implement the HTTP request and response parsing for "
        f"position {position} based on the current Abysis web interface."
    )


def refresh_cache(output_path: Path) -> None:
    """Fetch frequencies for all VH positions and write to JSON.

    Parameters
    ----------
    output_path : Path
        File path for the output JSON cache.
    """
    frequencies: dict[str, dict[str, float]] = {}

    logger.info("Fetching Abysis frequencies for %d positions", len(VH_IMGT_POSITIONS))
    for position in VH_IMGT_POSITIONS:
        try:
            freq = fetch_position_frequencies(position)
            frequencies[position] = freq
            logger.info("Position %s: %d amino acids", position, len(freq))
        except NotImplementedError:
            logger.warning(
                "Scraping not implemented – using empty data for position %s",
                position,
            )
            frequencies[position] = {}
        except requests.RequestException as exc:
            logger.error("Failed to fetch position %s: %s", position, exc)
            frequencies[position] = {}

        time.sleep(REQUEST_DELAY_SECONDS)

    output_path.parent.mkdir(parents=True, exist_ok=True)
    output_path.write_text(json.dumps(frequencies, indent=2, sort_keys=True))
    logger.info("Wrote frequency cache to %s", output_path)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Refresh Abysis human VH frequency cache"
    )
    parser.add_argument(
        "--output",
        "-o",
        default="data/abysis_vh_frequencies.json",
        help="Output path for the frequency cache JSON (default: data/abysis_vh_frequencies.json)",
    )
    parser.add_argument(
        "--verbose",
        "-v",
        action="store_true",
        help="Enable verbose logging",
    )
    args = parser.parse_args()

    logging.basicConfig(
        level=logging.DEBUG if args.verbose else logging.INFO,
        format="%(asctime)s - %(name)s - %(levelname)s - %(message)s",
    )

    output_path = Path(args.output)
    logger.info("Abysis cache will be written to: %s", output_path)
    refresh_cache(output_path)


if __name__ == "__main__":
    main()

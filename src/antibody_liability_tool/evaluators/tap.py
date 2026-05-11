"""TAP (Therapeutic Antibody Profiler) evaluator via SAbPred.

Supports two modes:
- **automation**: POST sequences directly to the TAP web service and parse
  the HTML response.
- **manual**: Write sequences to FASTA files for manual submission, then
  parse downloaded results.
"""

from __future__ import annotations

import csv
import io
import logging
import re
import time
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Sequence

import requests

from antibody_liability_tool.caching import cached
from antibody_liability_tool.evaluators.base import (
    BaseEvaluator,
    ComparisonResult,
    EvaluationResult,
)

logger = logging.getLogger(__name__)

_DEFAULT_TAP_URL = "https://opig.stats.ox.ac.uk/webapps/sabdab-sabpred/sabpred/tap"

# Metrics where a *lower* value means fewer liabilities.
_LOWER_IS_BETTER: set[str] = {"PSH", "PPC", "PNC", "SFvCSP"}


# ──────────────────────────────────────────────────────────────────────
# Data classes
# ──────────────────────────────────────────────────────────────────────


@dataclass
class TAPResult:
    """Structured result from a TAP analysis.

    Attributes
    ----------
    patches_positive_charge : int
        Number of positive charge patches.
    patches_negative_charge : int
        Number of negative charge patches.
    patches_hydrophobic : int
        Number of hydrophobic patches.
    PSH : float
        Patches of Surface Hydrophobicity score.
    PPC : float
        Patches of Positive Charge score.
    PNC : float
        Patches of Negative Charge score.
    SFvCSP : float
        Structural Fv Charge Symmetry Parameter.
    raw_html : str
        Raw HTML response (kept for debugging).
    """

    patches_positive_charge: int = 0
    patches_negative_charge: int = 0
    patches_hydrophobic: int = 0
    PSH: float = 0.0
    PPC: float = 0.0
    PNC: float = 0.0
    SFvCSP: float = 0.0
    raw_html: str = ""

    def to_metrics(self) -> dict[str, float]:
        """Return a flat metrics dict suitable for :class:`EvaluationResult`."""
        return {
            "patches_positive_charge": float(self.patches_positive_charge),
            "patches_negative_charge": float(self.patches_negative_charge),
            "patches_hydrophobic": float(self.patches_hydrophobic),
            "PSH": self.PSH,
            "PPC": self.PPC,
            "PNC": self.PNC,
            "SFvCSP": self.SFvCSP,
        }


# ──────────────────────────────────────────────────────────────────────
# HTML parsing helpers
# ──────────────────────────────────────────────────────────────────────

_METRIC_PATTERNS: dict[str, re.Pattern[str]] = {
    "PSH": re.compile(r"PSH\s*[=:]\s*([\d.eE+-]+)", re.IGNORECASE),
    "PPC": re.compile(r"PPC\s*[=:]\s*([\d.eE+-]+)", re.IGNORECASE),
    "PNC": re.compile(r"PNC\s*[=:]\s*([\d.eE+-]+)", re.IGNORECASE),
    "SFvCSP": re.compile(r"SFvCSP\s*[=:]\s*([\d.eE+-]+)", re.IGNORECASE),
}

_PATCH_PATTERNS: dict[str, re.Pattern[str]] = {
    "patches_positive_charge": re.compile(
        r"positive\s+charge\s+patch(?:es)?\s*[=:]\s*(\d+)", re.IGNORECASE
    ),
    "patches_negative_charge": re.compile(
        r"negative\s+charge\s+patch(?:es)?\s*[=:]\s*(\d+)", re.IGNORECASE
    ),
    "patches_hydrophobic": re.compile(r"hydrophobic\s+patch(?:es)?\s*[=:]\s*(\d+)", re.IGNORECASE),
}

# Fallback: look for values in HTML table cells (<td>)
_TD_PATTERN = re.compile(r"<td[^>]*>\s*([\d.eE+-]+)\s*</td>", re.IGNORECASE)


def _parse_tap_html(html: str) -> TAPResult:
    """Extract TAP metrics from the SAbPred HTML response."""
    result = TAPResult(raw_html=html)

    for name, pattern in _METRIC_PATTERNS.items():
        match = pattern.search(html)
        if match:
            setattr(result, name, float(match.group(1)))

    for name, pattern in _PATCH_PATTERNS.items():
        match = pattern.search(html)
        if match:
            setattr(result, name, int(match.group(1)))

    # Fallback: try to pull numbers from table cells if structured patterns
    # didn't match (the TAP page sometimes returns a simple HTML table).
    if result.PSH == 0.0 and result.PPC == 0.0:
        values = [float(m.group(1)) for m in _TD_PATTERN.finditer(html)]
        if len(values) >= 7:
            result.patches_positive_charge = int(values[0])
            result.patches_negative_charge = int(values[1])
            result.patches_hydrophobic = int(values[2])
            result.PSH = values[3]
            result.PPC = values[4]
            result.PNC = values[5]
            result.SFvCSP = values[6]

    return result


def _parse_tap_csv(text: str) -> list[TAPResult]:
    """Parse a CSV file downloaded from the TAP results page."""
    results: list[TAPResult] = []
    reader = csv.DictReader(io.StringIO(text))
    for row in reader:
        result = TAPResult()
        for attr in ("PSH", "PPC", "PNC", "SFvCSP"):
            if attr in row:
                try:
                    setattr(result, attr, float(row[attr]))
                except (ValueError, TypeError):
                    pass
        for attr in ("patches_positive_charge", "patches_negative_charge", "patches_hydrophobic"):
            if attr in row:
                try:
                    setattr(result, attr, int(row[attr]))
                except (ValueError, TypeError):
                    pass
        results.append(result)
    return results


# ──────────────────────────────────────────────────────────────────────
# Cached submission helper
# ──────────────────────────────────────────────────────────────────────


@cached(ttl=604800)
def _submit_sequence_to_tap(
    sequence: str,
    url: str,
    timeout: float,
) -> TAPResult:
    """POST *sequence* to TAP and return parsed metrics (cached on disk)."""
    payload = {"vh_sequence": sequence, "submit": "Submit"}
    resp = requests.post(url, data=payload, timeout=timeout)
    resp.raise_for_status()
    return _parse_tap_html(resp.text)


# ──────────────────────────────────────────────────────────────────────
# Evaluator
# ──────────────────────────────────────────────────────────────────────


class TAPEvaluator(BaseEvaluator):
    """Evaluate antibody sequences with the TAP web service.

    Parameters
    ----------
    config : dict[str, Any] | None
        TAP-specific configuration (``tap`` section of the project config).
    """

    def __init__(self, config: dict[str, Any] | None = None) -> None:
        super().__init__(config=config, name="TAP")
        self.url: str = self.config.get("url", _DEFAULT_TAP_URL)
        self.mode: str = self.config.get("mode", "manual")
        self.rate_limit: float = float(self.config.get("rate_limit_seconds", 5.0))
        self.max_retries: int = int(self.config.get("max_retries", 3))
        self.retry_delay: float = float(self.config.get("retry_delay_seconds", 10.0))
        self.timeout: float = float(self.config.get("timeout_seconds", 120))
        self.batch_dir: Path = Path(self.config.get("batch_dir", "output/tap_batches"))
        self._last_request_time: float = 0.0

    # ------------------------------------------------------------------
    # Rate-limiting
    # ------------------------------------------------------------------

    def _rate_limit_wait(self) -> None:
        """Block until the rate-limit interval has elapsed."""
        elapsed = time.monotonic() - self._last_request_time
        remaining = self.rate_limit - elapsed
        if remaining > 0:
            self._logger.debug("Rate-limiting: sleeping %.1fs", remaining)
            time.sleep(remaining)
        self._last_request_time = time.monotonic()

    # ------------------------------------------------------------------
    # Automation mode
    # ------------------------------------------------------------------

    def _submit_with_retry(self, sequence: str) -> TAPResult:
        """Submit with exponential-backoff retry logic."""
        last_exc: Exception | None = None
        for attempt in range(1, self.max_retries + 1):
            try:
                self._rate_limit_wait()
                return _submit_sequence_to_tap(sequence, self.url, self.timeout)
            except Exception as exc:  # noqa: BLE001
                last_exc = exc
                delay = self.retry_delay * (2 ** (attempt - 1))
                self._logger.warning(
                    "TAP attempt %d/%d failed (%s). Retrying in %.0fs…",
                    attempt,
                    self.max_retries,
                    exc,
                    delay,
                )
                time.sleep(delay)
        raise RuntimeError(f"TAP submission failed after {self.max_retries} retries") from last_exc

    # ------------------------------------------------------------------
    # Manual / batch mode
    # ------------------------------------------------------------------

    def write_batch(
        self,
        sequences: dict[str, str],
        filename: str = "tap_batch.fasta",
    ) -> Path:
        """Write sequences to a FASTA file for manual TAP submission.

        Parameters
        ----------
        sequences : dict[str, str]
            Mapping of sequence identifiers to amino-acid strings.
        filename : str
            Name of the output FASTA file.

        Returns
        -------
        Path
            Path to the written FASTA file.
        """
        self.batch_dir.mkdir(parents=True, exist_ok=True)
        fasta_path = self.batch_dir / filename
        with fasta_path.open("w") as fh:
            for seq_id, seq in sequences.items():
                fh.write(f">{seq_id}\n{seq}\n")
        self._logger.info("Wrote %d sequences to %s", len(sequences), fasta_path)
        return fasta_path

    def parse_results(self, path: str | Path) -> list[TAPResult]:
        """Parse downloaded TAP results from a CSV or HTML file.

        Parameters
        ----------
        path : str | Path
            Path to the results file.

        Returns
        -------
        list[TAPResult]
        """
        file_path = Path(path)
        text = file_path.read_text()
        if file_path.suffix.lower() == ".csv":
            return _parse_tap_csv(text)
        return [_parse_tap_html(text)]

    # ------------------------------------------------------------------
    # BaseEvaluator interface
    # ------------------------------------------------------------------

    def evaluate(self, sequence: str) -> EvaluationResult:
        """Evaluate a VH sequence using TAP.

        In *automation* mode the sequence is submitted to the web service.
        In *manual* mode an error-free stub result is returned so the
        pipeline can proceed; real data should be loaded via
        :meth:`parse_results`.
        """
        if self.mode == "automation":
            try:
                tap = self._submit_with_retry(sequence)
                return EvaluationResult(
                    sequence=sequence,
                    evaluator_name=self.name,
                    metrics=tap.to_metrics(),
                    raw={"tap_result": tap},
                    success=True,
                )
            except Exception as exc:  # noqa: BLE001
                self._logger.error("TAP evaluation failed: %s", exc)
                return EvaluationResult(
                    sequence=sequence,
                    evaluator_name=self.name,
                    success=False,
                    error_message=str(exc),
                )

        # Manual mode – return placeholder
        self._logger.info("TAP running in manual mode. Use write_batch() and parse_results().")
        return EvaluationResult(
            sequence=sequence,
            evaluator_name=self.name,
            success=True,
            raw={"mode": "manual"},
        )

    def evaluate_batch(self, sequences: Sequence[str]) -> list[EvaluationResult]:
        """Batch evaluation; in manual mode writes a FASTA file as well."""
        if self.mode == "manual":
            seq_dict = {f"seq_{i}": s for i, s in enumerate(sequences)}
            self.write_batch(seq_dict)
        return super().evaluate_batch(sequences)

    def compare_to_parent(
        self,
        parent_result: EvaluationResult,
        mutant_result: EvaluationResult,
    ) -> ComparisonResult:
        """Compare mutant vs parent using TAP-specific criteria.

        A mutant passes if PSH, PPC, and PNC do not increase.
        """
        return self.compare(
            parent_result,
            mutant_result,
            lower_is_better=_LOWER_IS_BETTER,
        )

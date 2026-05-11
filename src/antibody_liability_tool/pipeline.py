"""Main pipeline orchestrator for the Antibody Liability Reduction Tool.

Coordinates the full 10-stage workflow from sequence input through
liability detection, mutation generation, evaluation, optimisation,
and reporting.
"""

from __future__ import annotations

import logging
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

from antibody_liability_tool.config import load_config
from antibody_liability_tool.evaluators.base import EvaluationResult
from antibody_liability_tool.evaluators.deepsp import DeepSPEvaluator
from antibody_liability_tool.evaluators.oasis import OASisEvaluator
from antibody_liability_tool.evaluators.tap import TAPEvaluator
from antibody_liability_tool.liabilities.detector import Liability, detect_liabilities
from antibody_liability_tool.mutations.generator import CandidateMutation, generate_mutations
from antibody_liability_tool.numbering.imgt import number_sequence
from antibody_liability_tool.optimization.combinatorial import (
    MutationCombination,
    generate_combinations,
)
from antibody_liability_tool.optimization.scorer import CompositeScorer, ScoredCandidate
from antibody_liability_tool.reporting.fasta_export import export_fasta
from antibody_liability_tool.reporting.html_report import generate_html_report
from antibody_liability_tool.reporting.visualization import create_liability_map, create_radar_plot
from antibody_liability_tool.surface.exposure import classify_surface_exposure

logger = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Pipeline result container
# ---------------------------------------------------------------------------


@dataclass
class PipelineResult:
    """Container for all outputs produced by the pipeline.

    Attributes
    ----------
    sequence : str
        Input parent amino-acid sequence.
    numbered_sequence : dict[str, str]
        IMGT-numbered sequence.
    exposure_map : dict[str, dict[str, str]]
        Surface exposure classification per position.
    liabilities : list[Liability]
        Detected surface-exposed liabilities.
    candidate_mutations : list[CandidateMutation]
        Single-point mutation candidates.
    combinations : list[MutationCombination]
        Combinatorial expansion of mutations.
    parent_evaluations : dict[str, EvaluationResult]
        Parent sequence evaluations keyed by evaluator name.
    mutant_evaluations : dict[str, list[EvaluationResult]]
        Mutant evaluations keyed by evaluator name.
    ranked_candidates : list[ScoredCandidate]
        Final scored and ranked candidates.
    report_paths : dict[str, Path]
        Paths to generated output files.
    stages_completed : list[str]
        Names of stages that have been completed.
    elapsed_seconds : float
        Total pipeline runtime in seconds.
    errors : list[str]
        Non-fatal error messages collected during the run.
    """

    sequence: str = ""
    numbered_sequence: dict[str, str] = field(default_factory=dict)
    exposure_map: dict[str, dict[str, str]] = field(default_factory=dict)
    liabilities: list[Liability] = field(default_factory=list)
    candidate_mutations: list[CandidateMutation] = field(default_factory=list)
    combinations: list[MutationCombination] = field(default_factory=list)
    parent_evaluations: dict[str, EvaluationResult] = field(default_factory=dict)
    mutant_evaluations: dict[str, list[EvaluationResult]] = field(default_factory=dict)
    ranked_candidates: list[ScoredCandidate] = field(default_factory=list)
    report_paths: dict[str, Path] = field(default_factory=dict)
    stages_completed: list[str] = field(default_factory=list)
    elapsed_seconds: float = 0.0
    errors: list[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Pipeline class
# ---------------------------------------------------------------------------


class LiabilityReductionPipeline:
    """Orchestrates the full antibody liability reduction workflow.

    Parameters
    ----------
    config : dict[str, Any] | None
        Project configuration dict.  Loaded from defaults if not provided.
    config_path : str or Path, optional
        Path to YAML config file (alternative to passing ``config``).
    output_dir : str or Path
        Directory for pipeline outputs.
    """

    # Stage names for ordering and resume support
    STAGES = [
        "numbering",
        "surface_exposure",
        "liability_detection",
        "mutation_generation",
        "combinatorial_expansion",
        "parent_evaluation",
        "mutant_evaluation",
        "scoring",
        "filtering",
        "reporting",
    ]

    def __init__(
        self,
        config: dict[str, Any] | None = None,
        config_path: str | Path | None = None,
        output_dir: str | Path = "output",
    ) -> None:
        if config is not None:
            self.config = config
        elif config_path is not None:
            self.config = load_config(config_path)
        else:
            self.config = load_config()

        self.output_dir = Path(output_dir)
        self.output_dir.mkdir(parents=True, exist_ok=True)

        # Evaluators
        self.tap_evaluator = TAPEvaluator(self.config.get("tap"))
        self.deepsp_evaluator = DeepSPEvaluator(self.config.get("deepsp"))
        self.oasis_evaluator = OASisEvaluator(self.config.get("oasis"))
        self.scorer = CompositeScorer.from_config(self.config)

        self._result = PipelineResult()

    # ------------------------------------------------------------------
    # Full pipeline run
    # ------------------------------------------------------------------

    def run(
        self,
        sequence: str,
        resume_from: str | None = None,
    ) -> PipelineResult:
        """Execute the full pipeline.

        Parameters
        ----------
        sequence : str
            Input VH amino-acid sequence.
        resume_from : str, optional
            Stage name to resume from (skips prior stages, requires that
            ``self._result`` contains data from previous runs).

        Returns
        -------
        PipelineResult
            Complete pipeline output.
        """
        start_time = time.monotonic()
        self._result.sequence = sequence

        stages = list(self.STAGES)
        if resume_from:
            if resume_from not in stages:
                raise ValueError(
                    f"Unknown stage '{resume_from}'. Must be one of: {stages}"
                )
            start_idx = stages.index(resume_from)
            stages = stages[start_idx:]
            logger.info("Resuming pipeline from stage '%s'", resume_from)
        else:
            logger.info("Starting full pipeline for sequence of length %d", len(sequence))

        stage_methods = {
            "numbering": self._stage_numbering,
            "surface_exposure": self._stage_surface_exposure,
            "liability_detection": self._stage_liability_detection,
            "mutation_generation": self._stage_mutation_generation,
            "combinatorial_expansion": self._stage_combinatorial_expansion,
            "parent_evaluation": self._stage_parent_evaluation,
            "mutant_evaluation": self._stage_mutant_evaluation,
            "scoring": self._stage_scoring,
            "filtering": self._stage_filtering,
            "reporting": self._stage_reporting,
        }

        for stage_name in stages:
            logger.info("=== Stage: %s ===", stage_name)
            stage_start = time.monotonic()
            try:
                stage_methods[stage_name]()
                self._result.stages_completed.append(stage_name)
                elapsed = time.monotonic() - stage_start
                logger.info(
                    "Stage '%s' completed in %.2fs", stage_name, elapsed
                )
            except Exception as exc:
                error_msg = f"Stage '{stage_name}' failed: {exc}"
                logger.error(error_msg, exc_info=True)
                self._result.errors.append(error_msg)
                # Try to continue to reporting if possible
                if stage_name not in ("numbering", "surface_exposure", "liability_detection"):
                    logger.warning("Attempting to continue pipeline despite error")
                    continue
                break

        self._result.elapsed_seconds = time.monotonic() - start_time
        logger.info(
            "Pipeline finished in %.2fs (%d stages, %d errors)",
            self._result.elapsed_seconds,
            len(self._result.stages_completed),
            len(self._result.errors),
        )
        return self._result

    # ------------------------------------------------------------------
    # Individual stages
    # ------------------------------------------------------------------

    def _stage_numbering(self) -> None:
        """Stage 1: IMGT numbering."""
        self._result.numbered_sequence = number_sequence(self._result.sequence)
        logger.info(
            "Numbered %d positions", len(self._result.numbered_sequence)
        )

    def _stage_surface_exposure(self) -> None:
        """Stage 2: Surface exposure classification."""
        self._result.exposure_map = classify_surface_exposure(
            self._result.numbered_sequence
        )
        logger.info(
            "Classified %d positions for surface exposure",
            len(self._result.exposure_map),
        )

    def _stage_liability_detection(self) -> None:
        """Stage 3: Liability detection."""
        self._result.liabilities = detect_liabilities(
            self._result.numbered_sequence,
            self._result.exposure_map,
            self.config,
        )
        logger.info("Detected %d liabilities", len(self._result.liabilities))

    def _stage_mutation_generation(self) -> None:
        """Stage 4: Candidate mutation generation."""
        self._result.candidate_mutations = generate_mutations(
            self._result.liabilities,
            self._result.numbered_sequence,
            self.config,
        )
        logger.info(
            "Generated %d candidate mutations",
            len(self._result.candidate_mutations),
        )

    def _stage_combinatorial_expansion(self) -> None:
        """Stage 5: Combinatorial expansion of mutations."""
        ccfg = self.config.get("combinatorial", {})
        max_order = ccfg.get("max_order", 3)
        threshold = ccfg.get("bayesian_threshold", 50)

        self._result.combinations = generate_combinations(
            self._result.candidate_mutations,
            max_order=max_order,
            bayesian_threshold=threshold,
            evaluator_func=self._quick_evaluate,
            config=self.config,
        )
        logger.info(
            "Expanded to %d mutation combinations",
            len(self._result.combinations),
        )

    def _stage_parent_evaluation(self) -> None:
        """Stage 6: Evaluate parent sequence with all evaluators."""
        seq = self._result.sequence

        tap_result = self.tap_evaluator.evaluate(seq)
        self._result.parent_evaluations["TAP"] = tap_result

        deepsp_result = self.deepsp_evaluator.evaluate(seq)
        self._result.parent_evaluations["DeepSP"] = deepsp_result

        oasis_result = self.oasis_evaluator.evaluate(seq)
        self._result.parent_evaluations["OASis"] = oasis_result

        logger.info("Parent evaluation complete for %d evaluators", 3)

    def _stage_mutant_evaluation(self) -> None:
        """Stage 7: Evaluate mutant sequences."""
        numbered = self._result.numbered_sequence
        combos = self._result.combinations

        if not combos:
            logger.warning("No mutation combinations to evaluate")
            return

        self._result.mutant_evaluations = {"TAP": [], "DeepSP": [], "OASis": []}

        for combo in combos:
            seq = combo.apply(self._result.sequence, numbered)

            tap_result = self.tap_evaluator.evaluate(seq)
            self._result.mutant_evaluations["TAP"].append(tap_result)

            deepsp_result = self.deepsp_evaluator.evaluate(seq)
            self._result.mutant_evaluations["DeepSP"].append(deepsp_result)

            oasis_result = self.oasis_evaluator.evaluate(seq)
            self._result.mutant_evaluations["OASis"].append(oasis_result)

        logger.info("Evaluated %d mutant sequences", len(combos))

    def _stage_scoring(self) -> None:
        """Stage 8: Score and rank candidates."""
        parent_metrics = self._merge_parent_metrics()
        combos = self._result.combinations
        numbered = self._result.numbered_sequence

        candidates: list[dict[str, Any]] = []
        for i, combo in enumerate(combos):
            mutant_metrics = self._merge_mutant_metrics(i)
            mutations_data = [
                {
                    "position": m.position,
                    "original_aa": m.original_aa,
                    "proposed_aa": m.proposed_aa,
                }
                for m in combo.mutations
            ]
            candidates.append(
                {
                    "label": combo.label,
                    "sequence": combo.apply(self._result.sequence, numbered),
                    "parent_metrics": parent_metrics,
                    "mutant_metrics": mutant_metrics,
                    "combination": combo,
                    "mutations": mutations_data,
                }
            )

        self._result.ranked_candidates = self.scorer.rank(candidates)
        logger.info("Scored and ranked %d candidates", len(self._result.ranked_candidates))

    def _stage_filtering(self) -> None:
        """Stage 9: Filter candidates by constraints."""
        scfg = self.config.get("scoring", {})
        oasis_threshold = scfg.get("oasis_threshold", 0.0)

        pre_count = len(self._result.ranked_candidates)
        self._result.ranked_candidates = self.scorer.filter_candidates(
            self._result.ranked_candidates,
            min_oasis=oasis_threshold,
        )
        logger.info(
            "Filtering: %d → %d candidates",
            pre_count,
            len(self._result.ranked_candidates),
        )

    def _stage_reporting(self) -> None:
        """Stage 10: Generate reports."""
        top_n = self.config.get("output", {}).get("top_n_report", 10)
        formats = self.config.get("output", {}).get("formats", ["html", "fasta"])

        if "html" in formats:
            html_path = self.output_dir / "report.html"
            generate_html_report(
                parent_sequence=self._result.sequence,
                numbered_sequence=self._result.numbered_sequence,
                liabilities=self._result.liabilities,
                ranked_candidates=self._result.ranked_candidates[:top_n],
                output_path=html_path,
            )
            self._result.report_paths["html"] = html_path
            logger.info("HTML report: %s", html_path)

        if "fasta" in formats:
            fasta_path = self.output_dir / "candidates.fasta"
            export_fasta(
                parent_sequence=self._result.sequence,
                ranked_candidates=self._result.ranked_candidates,
                output_path=fasta_path,
                top_n=top_n,
            )
            self._result.report_paths["fasta"] = fasta_path
            logger.info("FASTA export: %s", fasta_path)

        # Visualisations
        try:
            parent_metrics = self._merge_parent_metrics()

            radar_path = self.output_dir / "radar_plot.html"
            create_radar_plot(
                parent_metrics=parent_metrics,
                candidates=self._result.ranked_candidates[:5],
                output_path=radar_path,
            )
            self._result.report_paths["radar"] = radar_path

            map_path = self.output_dir / "liability_map.html"
            create_liability_map(
                numbered_sequence=self._result.numbered_sequence,
                liabilities=self._result.liabilities,
                output_path=map_path,
            )
            self._result.report_paths["liability_map"] = map_path
        except Exception as exc:  # noqa: BLE001
            msg = f"Visualisation generation failed: {exc}"
            logger.warning(msg)
            self._result.errors.append(msg)

    # ------------------------------------------------------------------
    # Helper methods
    # ------------------------------------------------------------------

    def _merge_parent_metrics(self) -> dict[str, float]:
        """Merge metrics from all parent evaluations into a single dict."""
        merged: dict[str, float] = {}
        for eval_result in self._result.parent_evaluations.values():
            if eval_result.success:
                merged.update(eval_result.metrics)
        return merged

    def _merge_mutant_metrics(self, index: int) -> dict[str, float]:
        """Merge metrics from all evaluators for the i-th mutant."""
        merged: dict[str, float] = {}
        for evaluator_name, results_list in self._result.mutant_evaluations.items():
            if index < len(results_list) and results_list[index].success:
                merged.update(results_list[index].metrics)
        return merged

    def _quick_evaluate(self, combo: MutationCombination) -> dict[str, float]:
        """Quick evaluation for Bayesian optimisation loop.

        Uses DeepSP and OASis (fast) evaluators only.
        """
        seq = combo.apply(self._result.sequence, self._result.numbered_sequence)
        metrics: dict[str, float] = {}

        deepsp_result = self.deepsp_evaluator.evaluate(seq)
        if deepsp_result.success:
            metrics.update(deepsp_result.metrics)

        oasis_result = self.oasis_evaluator.evaluate(seq)
        if oasis_result.success:
            metrics.update(oasis_result.metrics)

        return metrics

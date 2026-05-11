"""Evaluator sub-package for the Antibody Liability Reduction Tool.

Re-exports all evaluator classes and result dataclasses for convenient access::

    from antibody_liability_tool.evaluators import TAPEvaluator, DeepSPEvaluator, OASisEvaluator
"""

from antibody_liability_tool.evaluators.base import (
    BaseEvaluator,
    ComparisonResult,
    EvaluationResult,
)
from antibody_liability_tool.evaluators.deepsp import (
    AGGREGATION_DESCRIPTORS,
    DEEPSP_DESCRIPTORS,
    DeepSPEvaluator,
    DeepSPResult,
)
from antibody_liability_tool.evaluators.oasis import (
    OASisEvaluator,
    OASisResult,
)
from antibody_liability_tool.evaluators.tap import (
    TAPEvaluator,
    TAPResult,
)

__all__ = [
    # Base
    "BaseEvaluator",
    "ComparisonResult",
    "EvaluationResult",
    # TAP
    "TAPEvaluator",
    "TAPResult",
    # DeepSP
    "DeepSPEvaluator",
    "DeepSPResult",
    "DEEPSP_DESCRIPTORS",
    "AGGREGATION_DESCRIPTORS",
    # OASis
    "OASisEvaluator",
    "OASisResult",
]

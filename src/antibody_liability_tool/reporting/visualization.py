"""Visualisation module for antibody liability analysis.

Creates interactive Plotly-based radar/spider plots comparing mutation
candidates against the parent sequence, and a liability map showing
liability positions along the sequence.
"""

from __future__ import annotations

import logging
from pathlib import Path
from typing import Any, Sequence

from antibody_liability_tool.optimization.scorer import ScoredCandidate

logger = logging.getLogger(__name__)

# Metrics displayed on the radar plot
_RADAR_METRICS = ["PSH", "PPC", "PNC", "oasis_humanness", "SAP_score", "developability_index"]
_RADAR_LABELS = ["PSH", "PPC", "PNC", "OASis", "SAP", "Developability"]


def _normalise_for_radar(
    parent_metrics: dict[str, float],
    candidates: Sequence[ScoredCandidate],
    metrics: list[str],
) -> tuple[list[float], list[list[float]], float]:
    """Normalise metrics to [0, 1] range for radar display.

    For metrics where *lower* is better (PSH, PPC, PNC, SAP,
    developability), the display is inverted so that *outer* = better.

    Returns
    -------
    tuple
        (parent_values, list of candidate_values, max_val used for normalisation)
    """
    lower_is_better = {"PSH", "PPC", "PNC", "SAP_score", "developability_index"}

    # Collect all values for normalisation
    all_vals: dict[str, list[float]] = {m: [] for m in metrics}
    for m in metrics:
        all_vals[m].append(parent_metrics.get(m, 0.0))
        for cand in candidates:
            all_vals[m].append(cand.mutant_metrics.get(m, 0.0))

    # Normalise
    ranges: dict[str, tuple[float, float]] = {}
    for m in metrics:
        vals = all_vals[m]
        mn, mx = min(vals), max(vals)
        if mx - mn < 1e-9:
            ranges[m] = (mn - 1.0, mx + 1.0)
        else:
            ranges[m] = (mn, mx)

    def _norm(metric: str, value: float) -> float:
        mn, mx = ranges[metric]
        normed = (value - mn) / (mx - mn)
        if metric in lower_is_better:
            normed = 1.0 - normed
        return max(0.0, min(1.0, normed))

    parent_vals = [_norm(m, parent_metrics.get(m, 0.0)) for m in metrics]
    cand_vals = [
        [_norm(m, c.mutant_metrics.get(m, 0.0)) for m in metrics] for c in candidates
    ]
    return parent_vals, cand_vals, 1.0


def create_radar_plot(
    parent_metrics: dict[str, float],
    candidates: Sequence[ScoredCandidate],
    max_candidates: int = 5,
    output_path: str | Path | None = None,
) -> str:
    """Create an interactive radar/spider plot comparing candidates vs parent.

    Parameters
    ----------
    parent_metrics : dict[str, float]
        Metrics from the parent sequence evaluation.
    candidates : Sequence[ScoredCandidate]
        Ranked mutation candidates.
    max_candidates : int
        Maximum number of candidates to display.
    output_path : str or Path, optional
        If provided, save as a standalone HTML file.

    Returns
    -------
    str
        HTML ``<div>`` element containing the Plotly chart.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        logger.warning("Plotly not installed – returning placeholder HTML")
        return "<div><p>Plotly is required for radar plots.</p></div>"

    top = list(candidates[:max_candidates])
    metrics = _RADAR_METRICS
    labels = _RADAR_LABELS

    parent_vals, cand_vals, _ = _normalise_for_radar(parent_metrics, top, metrics)

    fig = go.Figure()

    # Parent trace
    fig.add_trace(
        go.Scatterpolar(
            r=parent_vals + [parent_vals[0]],
            theta=labels + [labels[0]],
            fill="toself",
            name="Parent",
            line=dict(color="#2c3e50", width=2, dash="dash"),
            opacity=0.6,
        )
    )

    # Candidate traces
    colours = ["#3498db", "#2ecc71", "#e74c3c", "#f39c12", "#9b59b6"]
    for i, (cand, vals) in enumerate(zip(top, cand_vals)):
        fig.add_trace(
            go.Scatterpolar(
                r=vals + [vals[0]],
                theta=labels + [labels[0]],
                fill="toself",
                name=cand.label[:30],
                line=dict(color=colours[i % len(colours)], width=2),
                opacity=0.5,
            )
        )

    fig.update_layout(
        polar=dict(
            radialaxis=dict(visible=True, range=[0, 1], showticklabels=False),
        ),
        showlegend=True,
        title="Candidate Comparison (outer = better)",
        font=dict(family="Segoe UI, sans-serif"),
        width=700,
        height=500,
    )

    html_div = fig.to_html(full_html=False, include_plotlyjs="cdn")

    if output_path:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        full_html = fig.to_html(full_html=True, include_plotlyjs="cdn")
        out.write_text(full_html, encoding="utf-8")
        logger.info("Radar plot saved to %s", out)

    return html_div


def create_liability_map(
    numbered_sequence: dict[str, str],
    liabilities: list[Any],
    output_path: str | Path | None = None,
) -> str:
    """Create a visual liability map showing positions along the sequence.

    Parameters
    ----------
    numbered_sequence : dict[str, str]
        IMGT numbered sequence mapping.
    liabilities : list
        Detected liabilities with ``position``, ``residue``, ``severity``.
    output_path : str or Path, optional
        If provided, save as a standalone HTML file.

    Returns
    -------
    str
        HTML ``<div>`` containing the liability map.
    """
    try:
        import plotly.graph_objects as go
    except ImportError:
        logger.warning("Plotly not installed – returning placeholder HTML")
        return "<div><p>Plotly is required for liability maps.</p></div>"

    def _sort_key(p: str) -> tuple[int, str]:
        num, suffix = "", ""
        for ch in p:
            if ch.isdigit():
                num += ch
            else:
                suffix += ch
        return (int(num) if num else 0, suffix)

    sorted_positions = sorted(numbered_sequence.keys(), key=_sort_key)
    liability_map: dict[str, int] = {}
    for li in liabilities:
        liability_map[li.position] = li.severity

    x_vals: list[int] = []
    y_vals: list[float] = []
    colours: list[str] = []
    hover_text: list[str] = []
    severity_colours = {3: "#e74c3c", 2: "#f39c12", 1: "#f1c40f"}

    for i, pos in enumerate(sorted_positions):
        aa = numbered_sequence[pos]
        x_vals.append(i)
        if pos in liability_map:
            sev = liability_map[pos]
            y_vals.append(float(sev))
            colours.append(severity_colours.get(sev, "#95a5a6"))
            hover_text.append(f"IMGT {pos}: {aa} (severity {sev})")
        else:
            y_vals.append(0.0)
            colours.append("#bdc3c7")
            hover_text.append(f"IMGT {pos}: {aa}")

    fig = go.Figure()
    fig.add_trace(
        go.Bar(
            x=x_vals,
            y=y_vals,
            marker_color=colours,
            text=hover_text,
            hoverinfo="text",
        )
    )

    # Region annotations
    from antibody_liability_tool.numbering.imgt import IMGT_REGIONS

    region_boundaries: list[tuple[str, int, int]] = []
    for region_name, (start, end) in IMGT_REGIONS.items():
        # Find index range in sorted positions
        start_idx = None
        end_idx = None
        for idx, pos in enumerate(sorted_positions):
            pos_int = _sort_key(pos)[0]
            if start_idx is None and pos_int >= start:
                start_idx = idx
            if pos_int <= end:
                end_idx = idx
        if start_idx is not None and end_idx is not None:
            region_boundaries.append((region_name, start_idx, end_idx))

    for region_name, start_idx, end_idx in region_boundaries:
        mid = (start_idx + end_idx) / 2
        is_cdr = "CDR" in region_name
        fig.add_shape(
            type="rect",
            x0=start_idx - 0.5,
            x1=end_idx + 0.5,
            y0=-0.3,
            y1=-0.1,
            fillcolor="#3498db" if is_cdr else "#ecf0f1",
            line=dict(width=0),
        )
        fig.add_annotation(
            x=mid,
            y=-0.5,
            text=region_name,
            showarrow=False,
            font=dict(size=9),
        )

    fig.update_layout(
        title="Liability Map (bar height = severity)",
        xaxis_title="Sequence Position",
        yaxis_title="Severity",
        yaxis=dict(range=[-0.8, 4]),
        showlegend=False,
        width=1000,
        height=350,
        font=dict(family="Segoe UI, sans-serif"),
    )

    html_div = fig.to_html(full_html=False, include_plotlyjs="cdn")

    if output_path:
        out = Path(output_path)
        out.parent.mkdir(parents=True, exist_ok=True)
        full_html = fig.to_html(full_html=True, include_plotlyjs="cdn")
        out.write_text(full_html, encoding="utf-8")
        logger.info("Liability map saved to %s", out)

    return html_div

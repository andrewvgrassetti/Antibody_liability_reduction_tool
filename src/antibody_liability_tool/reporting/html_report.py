"""HTML report generation using Jinja2 templates.

Generates a self-contained HTML report with inline CSS covering:
- Parent sequence analysis with highlighted liabilities
- Ranked table of evaluated mutants with metrics
- Sequence alignment of top candidates
- Summary statistics
"""

from __future__ import annotations

import html
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from jinja2 import BaseLoader, Environment

from antibody_liability_tool.optimization.scorer import ScoredCandidate

logger = logging.getLogger(__name__)

# ---------------------------------------------------------------------------
# Inline CSS
# ---------------------------------------------------------------------------
_CSS = """\
* { box-sizing: border-box; margin: 0; padding: 0; }
body {
    font-family: 'Segoe UI', Tahoma, Geneva, Verdana, sans-serif;
    line-height: 1.6; color: #333; background: #f8f9fa; padding: 20px;
}
.container { max-width: 1200px; margin: 0 auto; }
h1 { color: #2c3e50; margin-bottom: 5px; }
h2 { color: #34495e; margin: 25px 0 10px; border-bottom: 2px solid #3498db; padding-bottom: 5px; }
h3 { color: #555; margin: 15px 0 8px; }
.subtitle { color: #7f8c8d; font-size: 0.9em; margin-bottom: 20px; }
.summary-grid {
    display: grid; grid-template-columns: repeat(auto-fit, minmax(200px, 1fr));
    gap: 15px; margin-bottom: 20px;
}
.summary-card {
    background: white; border-radius: 8px; padding: 15px;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1); text-align: center;
}
.summary-card .value { font-size: 2em; font-weight: bold; color: #2c3e50; }
.summary-card .label { font-size: 0.85em; color: #7f8c8d; }
.seq-display {
    font-family: 'Courier New', monospace; font-size: 13px;
    background: white; border: 1px solid #ddd; border-radius: 5px;
    padding: 12px; word-break: break-all; line-height: 1.8;
    overflow-x: auto; margin-bottom: 15px;
}
.liability { background: #e74c3c; color: white; padding: 1px 3px; border-radius: 3px; }
.liability-medium { background: #f39c12; color: white; padding: 1px 3px; border-radius: 3px; }
.liability-low { background: #f1c40f; color: #333; padding: 1px 3px; border-radius: 3px; }
table {
    width: 100%; border-collapse: collapse; background: white;
    border-radius: 8px; overflow: hidden;
    box-shadow: 0 2px 4px rgba(0,0,0,0.1); margin-bottom: 20px;
}
th {
    background: #2c3e50; color: white; padding: 10px 12px;
    text-align: left; font-size: 0.85em; text-transform: uppercase;
}
td { padding: 8px 12px; border-bottom: 1px solid #eee; font-size: 0.9em; }
tr:hover { background: #f5f6fa; }
tr:nth-child(even) { background: #fafbfc; }
.score-bar {
    height: 20px; border-radius: 10px; background: #ecf0f1;
    position: relative; min-width: 100px;
}
.score-fill {
    height: 100%; border-radius: 10px;
    background: linear-gradient(90deg, #e74c3c, #f39c12, #2ecc71);
    transition: width 0.3s;
}
.improved { color: #27ae60; font-weight: 600; }
.worsened { color: #e74c3c; font-weight: 600; }
.neutral { color: #95a5a6; }
.alignment {
    font-family: 'Courier New', monospace; font-size: 12px;
    background: white; border: 1px solid #ddd; border-radius: 5px;
    padding: 10px; overflow-x: auto; margin-bottom: 15px; white-space: pre;
}
.mutation-pos { background: #2ecc71; color: white; padding: 1px 2px; border-radius: 2px; }
.footer { margin-top: 30px; text-align: center; color: #95a5a6; font-size: 0.8em; }
"""

# ---------------------------------------------------------------------------
# Jinja2 Template
# ---------------------------------------------------------------------------
_TEMPLATE = """\
<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>Antibody Liability Reduction Report</title>
    <style>{{ css }}</style>
</head>
<body>
<div class="container">
    <h1>Antibody Liability Reduction Report</h1>
    <p class="subtitle">Generated {{ timestamp }}</p>

    <h2>Summary</h2>
    <div class="summary-grid">
        <div class="summary-card">
            <div class="value">{{ sequence_length }}</div>
            <div class="label">Sequence Length</div>
        </div>
        <div class="summary-card">
            <div class="value">{{ n_liabilities }}</div>
            <div class="label">Liabilities Detected</div>
        </div>
        <div class="summary-card">
            <div class="value">{{ n_mutations }}</div>
            <div class="label">Candidate Mutations</div>
        </div>
        <div class="summary-card">
            <div class="value">{{ n_evaluated }}</div>
            <div class="label">Candidates Evaluated</div>
        </div>
        {% if top_score is not none %}
        <div class="summary-card">
            <div class="value">{{ "%.3f"|format(top_score) }}</div>
            <div class="label">Top Composite Score</div>
        </div>
        {% endif %}
    </div>

    <h2>Parent Sequence</h2>
    <div class="seq-display">{{ parent_seq_html }}</div>

    {% if liabilities %}
    <h2>Detected Liabilities</h2>
    <table>
        <thead>
            <tr>
                <th>Position</th><th>Residue</th><th>Region</th>
                <th>Reason</th><th>Severity</th>
            </tr>
        </thead>
        <tbody>
        {% for li in liabilities %}
            <tr>
                <td>{{ li.position }}</td>
                <td>{{ li.residue }}</td>
                <td>{{ li.region }}</td>
                <td>{{ li.reason }}</td>
                <td>{{ li.severity_label }}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    {% endif %}

    {% if ranked_candidates %}
    <h2>Ranked Candidates</h2>
    <table>
        <thead>
            <tr>
                <th>Rank</th><th>Mutation</th><th>Score</th>
                <th>ΔPSH</th><th>ΔPPC</th><th>ΔPNC</th>
                <th>ΔOASis</th><th>ΔSAP</th>
            </tr>
        </thead>
        <tbody>
        {% for cand in ranked_candidates %}
            <tr>
                <td>{{ loop.index }}</td>
                <td>{{ cand.label }}</td>
                <td>
                    <div class="score-bar">
                        <div class="score-fill"
                             style="width: {{ (cand.composite_score * 100)|int }}%"></div>
                    </div>
                    {{ "%.4f"|format(cand.composite_score) }}
                </td>
                <td class="{{ cand.psh_class }}">{{ cand.psh_delta }}</td>
                <td class="{{ cand.ppc_class }}">{{ cand.ppc_delta }}</td>
                <td class="{{ cand.pnc_class }}">{{ cand.pnc_delta }}</td>
                <td class="{{ cand.oasis_class }}">{{ cand.oasis_delta }}</td>
                <td class="{{ cand.sap_class }}">{{ cand.sap_delta }}</td>
            </tr>
        {% endfor %}
        </tbody>
    </table>
    {% endif %}

    {% if alignment_html %}
    <h2>Sequence Alignment (Top Candidates)</h2>
    <div class="alignment">{{ alignment_html }}</div>
    {% endif %}

    <div class="footer">
        <p>Antibody Liability Reduction Tool v{{ version }} &mdash; {{ timestamp }}</p>
    </div>
</div>
</body>
</html>
"""


# ---------------------------------------------------------------------------
# Helper functions
# ---------------------------------------------------------------------------


def _severity_label(severity: int) -> str:
    """Convert integer severity to a label."""
    return {3: "High", 2: "Medium", 1: "Low"}.get(severity, "Unknown")


def _delta_class(delta: float, lower_is_better: bool = True) -> str:
    """Return CSS class for a metric delta."""
    if abs(delta) < 1e-6:
        return "neutral"
    if lower_is_better:
        return "improved" if delta < 0 else "worsened"
    return "improved" if delta > 0 else "worsened"


def _format_delta(value: float) -> str:
    """Format a delta value with sign."""
    if abs(value) < 1e-6:
        return "0.000"
    return f"{value:+.3f}"


def _build_parent_seq_html(
    sequence: str,
    numbered: dict[str, str],
    liability_positions: set[str],
    severity_map: dict[str, int],
) -> str:
    """Build HTML of the parent sequence with liability positions highlighted."""

    def _sort_key(p: str) -> tuple[int, str]:
        num, suffix = "", ""
        for ch in p:
            if ch.isdigit():
                num += ch
            else:
                suffix += ch
        return (int(num) if num else 0, suffix)

    sorted_positions = sorted(numbered.keys(), key=_sort_key)
    parts: list[str] = []
    for pos in sorted_positions:
        aa = html.escape(numbered[pos])
        if pos in liability_positions:
            sev = severity_map.get(pos, 1)
            css_class = {3: "liability", 2: "liability-medium"}.get(sev, "liability-low")
            parts.append(f'<span class="{css_class}" title="IMGT {pos}">{aa}</span>')
        else:
            parts.append(aa)
    return "".join(parts)


def _build_alignment_html(
    parent_numbered: dict[str, str],
    candidates: list[ScoredCandidate],
    max_candidates: int = 5,
) -> str:
    """Build a simple text alignment of top candidates vs parent."""

    def _sort_key(p: str) -> tuple[int, str]:
        num, suffix = "", ""
        for ch in p:
            if ch.isdigit():
                num += ch
            else:
                suffix += ch
        return (int(num) if num else 0, suffix)

    sorted_positions = sorted(parent_numbered.keys(), key=_sort_key)
    parent_seq = "".join(parent_numbered[p] for p in sorted_positions)

    lines: list[str] = []
    lines.append(f"Parent:  {parent_seq}")

    for i, cand in enumerate(candidates[:max_candidates]):
        mutant_seq = list(parent_seq)
        mutations = cand.metadata.get("mutations", [])
        mutation_indices: set[int] = set()

        for mut in mutations:
            pos = mut.get("position", "")
            proposed = mut.get("proposed_aa", "")
            if pos in sorted_positions and proposed:
                idx = sorted_positions.index(pos)
                mutant_seq[idx] = proposed
                mutation_indices.add(idx)

        # Build display with dots for unchanged positions
        display: list[str] = []
        for idx, (parent_aa, mut_aa) in enumerate(zip(parent_seq, mutant_seq)):
            if parent_aa == mut_aa:
                display.append(".")
            else:
                display.append(mut_aa)

        label = cand.label[:8].ljust(8)
        display_str = "".join(display)
        lines.append(f"{label} {display_str}")

    return html.escape("\n".join(lines))


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------


def generate_html_report(
    parent_sequence: str,
    numbered_sequence: dict[str, str],
    liabilities: list[Any],
    ranked_candidates: list[ScoredCandidate],
    output_path: str | Path,
    version: str = "0.1.0",
) -> Path:
    """Generate a comprehensive HTML report.

    Parameters
    ----------
    parent_sequence : str
        The parent (wild-type) amino-acid sequence.
    numbered_sequence : dict[str, str]
        IMGT numbered sequence mapping.
    liabilities : list
        Detected liabilities (``Liability`` dataclass instances).
    ranked_candidates : list[ScoredCandidate]
        Scored and ranked mutation candidates.
    output_path : str or Path
        Path to write the HTML file.
    version : str
        Tool version string for the footer.

    Returns
    -------
    Path
        Absolute path to the generated report.
    """
    output = Path(output_path)
    output.parent.mkdir(parents=True, exist_ok=True)

    # Build liability data for template
    liability_positions: set[str] = set()
    severity_map: dict[str, int] = {}
    liability_rows: list[dict[str, str]] = []
    for li in liabilities:
        liability_positions.add(li.position)
        severity_map[li.position] = li.severity
        liability_rows.append(
            {
                "position": li.imgt_number,
                "residue": li.residue,
                "region": li.region,
                "reason": li.reason,
                "severity_label": _severity_label(li.severity),
            }
        )

    # Build parent sequence HTML
    parent_seq_html = _build_parent_seq_html(
        parent_sequence, numbered_sequence, liability_positions, severity_map
    )

    # Build candidate rows
    candidate_rows: list[dict[str, Any]] = []
    for cand in ranked_candidates:
        d = cand.deltas
        candidate_rows.append(
            {
                "label": cand.label,
                "composite_score": cand.composite_score,
                "psh_delta": _format_delta(d.get("PSH", 0.0)),
                "psh_class": _delta_class(d.get("PSH", 0.0), lower_is_better=True),
                "ppc_delta": _format_delta(d.get("PPC", 0.0)),
                "ppc_class": _delta_class(d.get("PPC", 0.0), lower_is_better=True),
                "pnc_delta": _format_delta(d.get("PNC", 0.0)),
                "pnc_class": _delta_class(d.get("PNC", 0.0), lower_is_better=True),
                "oasis_delta": _format_delta(d.get("oasis_humanness", 0.0)),
                "oasis_class": _delta_class(d.get("oasis_humanness", 0.0), lower_is_better=False),
                "sap_delta": _format_delta(d.get("SAP_score", 0.0)),
                "sap_class": _delta_class(d.get("SAP_score", 0.0), lower_is_better=True),
            }
        )

    # Build alignment
    alignment_html = _build_alignment_html(numbered_sequence, ranked_candidates)

    # Render template
    env = Environment(loader=BaseLoader(), autoescape=False)
    template = env.from_string(_TEMPLATE)

    rendered = template.render(
        css=_CSS,
        timestamp=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
        sequence_length=len(parent_sequence),
        n_liabilities=len(liabilities),
        n_mutations=len(ranked_candidates),
        n_evaluated=len(ranked_candidates),
        top_score=ranked_candidates[0].composite_score if ranked_candidates else None,
        parent_seq_html=parent_seq_html,
        liabilities=liability_rows,
        ranked_candidates=candidate_rows,
        alignment_html=alignment_html,
        version=version,
    )

    output.write_text(rendered, encoding="utf-8")
    logger.info("HTML report written to %s", output)
    return output.resolve()

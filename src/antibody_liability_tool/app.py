"""Streamlit-based web UI for the Antibody Liability Reduction Tool."""

from __future__ import annotations

import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

import streamlit as st
import streamlit.components.v1 as components

from antibody_liability_tool.config import _deep_merge, load_config
from antibody_liability_tool.pipeline import LiabilityReductionPipeline, PipelineResult

# ---------------------------------------------------------------------------
# Page config
# ---------------------------------------------------------------------------

st.set_page_config(
    page_title="Antibody Liability Reduction Tool",
    page_icon="🧬",
    layout="wide",
)

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

_STAGE_LABELS: dict[str, str] = {
    "numbering": "Stage 1 — IMGT Numbering",
    "surface_exposure": "Stage 2 — Surface Exposure Classification",
    "liability_detection": "Stage 3 — Liability Detection",
    "mutation_generation": "Stage 4 — Mutation Generation",
    "combinatorial_expansion": "Stage 5 — Combinatorial Expansion",
    "parent_evaluation": "Stage 6 — Parent Evaluation",
    "mutant_evaluation": "Stage 7 — Mutant Evaluation",
    "scoring": "Stage 8 — Scoring & Ranking",
    "filtering": "Stage 9 — Filtering",
    "reporting": "Stage 10 — Report Generation",
}


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------


def _read_sequence_from_fasta(content: str) -> str:
    """Parse a FASTA string and return the amino-acid sequence."""
    lines = content.splitlines()
    seq_lines = [line.strip() for line in lines if not line.startswith(">")]
    return "".join(seq_lines)


def _build_config_overrides(
    min_human_freq: float,
    max_order: int,
    oasis_threshold: float,
    top_n: int,
) -> dict:
    """Return a nested dict of config overrides from sidebar widgets."""
    return {
        "mutations": {"min_human_frequency": min_human_freq},
        "combinatorial": {"max_order": max_order},
        "scoring": {"oasis_threshold": oasis_threshold},
        "output": {"top_n_report": top_n},
    }


# ---------------------------------------------------------------------------
# Sidebar — Input & Configuration
# ---------------------------------------------------------------------------

st.sidebar.title("🧬 Input & Configuration")

st.sidebar.subheader("Sequence Input")
fasta_file = st.sidebar.file_uploader(
    "Upload FASTA file", type=["fasta", "fa"], key="fasta_upload"
)
raw_sequence = st.sidebar.text_area(
    "Or paste a raw VH sequence", height=120, key="raw_seq"
)

st.sidebar.subheader("Pipeline Settings")
min_human_freq = st.sidebar.slider(
    "Min human frequency",
    min_value=0.01,
    max_value=1.0,
    value=0.05,
    step=0.01,
    help="mutations.min_human_frequency",
)
max_order = st.sidebar.slider(
    "Max combination order",
    min_value=1,
    max_value=3,
    value=3,
    help="combinatorial.max_order",
)
oasis_threshold = st.sidebar.slider(
    "OASis humanness threshold",
    min_value=0.0,
    max_value=1.0,
    value=0.5,
    step=0.01,
    help="scoring.oasis_threshold",
)
top_n = st.sidebar.number_input(
    "Top-N candidates to report",
    min_value=1,
    max_value=100,
    value=10,
    help="output.top_n_report",
)

st.sidebar.divider()
run_full = st.sidebar.button("🚀 Run Pipeline", type="primary", use_container_width=True)
run_analyze = st.sidebar.button(
    "🔍 Analyze Only (Stages 1-3)", use_container_width=True
)


# ---------------------------------------------------------------------------
# Main area — Title
# ---------------------------------------------------------------------------

st.title("🧬 Antibody Liability Reduction Tool")
st.markdown(
    "Identify and reduce surface-exposed liabilities in antibody VH sequences. "
    "Upload a FASTA file or paste a sequence, configure parameters, then run the pipeline."
)

# ---------------------------------------------------------------------------
# Sequence resolution
# ---------------------------------------------------------------------------


def _resolve_sequence() -> str | None:
    """Return the sequence or show an error and return None."""
    has_fasta = fasta_file is not None
    has_raw = bool(raw_sequence and raw_sequence.strip())

    if has_fasta and has_raw:
        st.error("Please provide **either** a FASTA file **or** a raw sequence — not both.")
        return None
    if not has_fasta and not has_raw:
        return None  # nothing submitted yet

    if has_fasta:
        content = fasta_file.read().decode("utf-8")  # type: ignore[union-attr]
        return _read_sequence_from_fasta(content)
    return raw_sequence.strip()  # type: ignore[union-attr]


# ---------------------------------------------------------------------------
# Display helpers
# ---------------------------------------------------------------------------


def _show_summary(result: PipelineResult, analyze_only: bool = False) -> None:
    """Render the Summary tab."""
    cols = st.columns(4 if not analyze_only else 3)
    cols[0].metric("Sequence Length", len(result.sequence))
    cols[1].metric("Liabilities Found", len(result.liabilities))
    if not analyze_only:
        cols[2].metric("Mutations Proposed", len(result.candidate_mutations))
        cols[3].metric(
            "Candidates After Filtering", len(result.ranked_candidates)
        )
    else:
        cols[2].metric("Mutations Proposed", len(result.candidate_mutations))

    st.metric("Pipeline Runtime", f"{result.elapsed_seconds:.2f} s")


def _show_liabilities(result: PipelineResult) -> None:
    """Render the Liabilities tab."""
    import pandas as pd

    if not result.liabilities:
        st.info("No liabilities detected.")
        return

    rows = [
        {
            "IMGT Number": li.imgt_number,
            "Residue": li.residue,
            "Region": li.region,
            "Reason": li.reason,
            "Severity": li.severity,
        }
        for li in result.liabilities
    ]
    st.dataframe(pd.DataFrame(rows), use_container_width=True)


def _show_candidates(result: PipelineResult) -> None:
    """Render the Candidates tab."""
    import pandas as pd

    if not result.ranked_candidates:
        st.info("No candidates to display.")
        return

    rows = []
    for rank, cand in enumerate(result.ranked_candidates, 1):
        rows.append(
            {
                "Rank": rank,
                "Label": cand.label,
                "Composite Score": round(cand.composite_score, 4),
                "PSH Reduction": round(cand.deltas.get("PSH", 0.0), 4),
                "PPC Reduction": round(cand.deltas.get("PPC", 0.0), 4),
                "OASis Delta": round(cand.deltas.get("oasis_humanness", 0.0), 4),
                "Stability Penalty": round(cand.deltas.get("stability_penalty", 0.0), 4),
            }
        )

    df = pd.DataFrame(rows)
    st.dataframe(df, use_container_width=True)

    csv = df.to_csv(index=False).encode("utf-8")
    st.download_button(
        "⬇️ Download as CSV",
        data=csv,
        file_name="candidates.csv",
        mime="text/csv",
    )


def _show_visualizations(result: PipelineResult) -> None:
    """Render the Visualizations tab using Plotly figures inline."""
    from antibody_liability_tool.reporting.visualization import (
        create_liability_map,
        create_radar_plot,
    )

    if result.ranked_candidates and result.parent_evaluations:
        # Merge parent metrics
        parent_metrics: dict[str, float] = {}
        for ev in result.parent_evaluations.values():
            if ev.success:
                parent_metrics.update(ev.metrics)

        st.subheader("Radar Plot — Candidate Comparison")
        radar_html = create_radar_plot(
            parent_metrics=parent_metrics,
            candidates=result.ranked_candidates[:5],
        )
        components.html(radar_html, height=550, scrolling=True)
    else:
        st.info("Radar plot requires scored candidates and parent evaluations.")

    if result.numbered_sequence and result.liabilities:
        st.subheader("Liability Map")
        map_html = create_liability_map(
            numbered_sequence=result.numbered_sequence,
            liabilities=result.liabilities,
        )
        components.html(map_html, height=400, scrolling=True)
    else:
        st.info("Liability map requires numbered sequence and liabilities data.")


def _show_report(result: PipelineResult) -> None:
    """Render the Report tab."""
    html_path = result.report_paths.get("html")
    if html_path and Path(html_path).exists():
        html_content = Path(html_path).read_text(encoding="utf-8")
        components.html(html_content, height=800, scrolling=True)
        st.download_button(
            "⬇️ Download HTML Report",
            data=html_content,
            file_name="report.html",
            mime="text/html",
        )
    else:
        st.info("HTML report not available.")


def _show_fasta(result: PipelineResult) -> None:
    """Render the FASTA Export tab."""
    fasta_path = result.report_paths.get("fasta")
    if fasta_path and Path(fasta_path).exists():
        fasta_content = Path(fasta_path).read_text(encoding="utf-8")
        st.code(fasta_content, language="text")
        st.download_button(
            "⬇️ Download FASTA",
            data=fasta_content,
            file_name="candidates.fasta",
            mime="text/plain",
        )
    else:
        st.info("FASTA export not available.")


def _show_errors(result: PipelineResult) -> None:
    """Display pipeline warnings/errors."""
    for err in result.errors:
        st.warning(err)


# ---------------------------------------------------------------------------
# Pipeline execution
# ---------------------------------------------------------------------------

if run_full or run_analyze:
    seq = _resolve_sequence()
    if seq is None and not (run_full or run_analyze):
        pass  # user hasn't provided input yet
    elif seq is None:
        st.error("Please provide a VH sequence using one of the input methods in the sidebar.")
    else:
        # Build config
        base_cfg = load_config()
        overrides = _build_config_overrides(min_human_freq, max_order, oasis_threshold, top_n)
        cfg = _deep_merge(base_cfg, overrides)

        try:
            with tempfile.TemporaryDirectory() as tmp_dir:
                pipeline = LiabilityReductionPipeline(
                    config=cfg, output_dir=tmp_dir
                )

                if run_analyze:
                    # --- Analyze-only mode (stages 1-3) ---
                    import time as _time

                    start_time = _time.monotonic()
                    pipeline._result.sequence = seq

                    with st.status("Running analysis (Stages 1-3)…", expanded=True) as status:
                        analyze_stages = ["numbering", "surface_exposure", "liability_detection"]
                        failed = False
                        for stage_name in analyze_stages:
                            st.write(f"⏳ {_STAGE_LABELS[stage_name]}")
                            try:
                                {
                                    "numbering": pipeline._stage_numbering,
                                    "surface_exposure": pipeline._stage_surface_exposure,
                                    "liability_detection": pipeline._stage_liability_detection,
                                }[stage_name]()
                                pipeline._result.stages_completed.append(stage_name)
                                st.write(f"✅ {_STAGE_LABELS[stage_name]}")
                            except Exception as exc:
                                error_msg = f"Stage '{stage_name}' failed: {exc}"
                                pipeline._result.errors.append(error_msg)
                                st.error(f"Critical stage failed: {error_msg}")
                                status.update(label="Analysis failed", state="error")
                                failed = True
                                break

                        if not failed:
                            status.update(label="Analysis complete!", state="complete")

                    pipeline._result.elapsed_seconds = _time.monotonic() - start_time

                    result = pipeline._result

                    # Show warnings
                    _show_errors(result)

                    # Show Summary and Liabilities tabs only
                    tab_summary, tab_liabilities = st.tabs(["📊 Summary", "⚠️ Liabilities"])
                    with tab_summary:
                        _show_summary(result, analyze_only=True)
                    with tab_liabilities:
                        _show_liabilities(result)

                else:
                    # --- Full pipeline ---
                    pipeline._result.sequence = seq

                    with st.status("Running full pipeline…", expanded=True) as status:
                        import time as _time

                        start_time = _time.monotonic()

                        for stage_name in LiabilityReductionPipeline.STAGES:
                            st.write(f"⏳ {_STAGE_LABELS[stage_name]}")
                            stage_method = {
                                "numbering": pipeline._stage_numbering,
                                "surface_exposure": pipeline._stage_surface_exposure,
                                "liability_detection": pipeline._stage_liability_detection,
                                "mutation_generation": pipeline._stage_mutation_generation,
                                "combinatorial_expansion": pipeline._stage_combinatorial_expansion,
                                "parent_evaluation": pipeline._stage_parent_evaluation,
                                "mutant_evaluation": pipeline._stage_mutant_evaluation,
                                "scoring": pipeline._stage_scoring,
                                "filtering": pipeline._stage_filtering,
                                "reporting": pipeline._stage_reporting,
                            }[stage_name]

                            try:
                                stage_method()
                                pipeline._result.stages_completed.append(stage_name)
                                st.write(f"✅ {_STAGE_LABELS[stage_name]}")
                            except Exception as exc:
                                error_msg = f"Stage '{stage_name}' failed: {exc}"
                                pipeline._result.errors.append(error_msg)
                                if stage_name in (
                                    "numbering",
                                    "surface_exposure",
                                    "liability_detection",
                                ):
                                    st.error(f"Critical stage failed: {error_msg}")
                                    status.update(
                                        label="Pipeline failed", state="error"
                                    )
                                    break
                                else:
                                    st.warning(f"Non-critical failure: {error_msg}")
                        else:
                            status.update(
                                label="Pipeline complete!", state="complete"
                            )

                        pipeline._result.elapsed_seconds = (
                            _time.monotonic() - start_time
                        )

                    result = pipeline._result

                    # Show warnings
                    _show_errors(result)

                    # Tabs
                    tab_names = [
                        "📊 Summary",
                        "⚠️ Liabilities",
                        "🏆 Candidates",
                        "📈 Visualizations",
                        "📄 Report",
                        "🧬 FASTA Export",
                    ]
                    tabs = st.tabs(tab_names)

                    with tabs[0]:
                        _show_summary(result)
                    with tabs[1]:
                        _show_liabilities(result)
                    with tabs[2]:
                        _show_candidates(result)
                    with tabs[3]:
                        _show_visualizations(result)
                    with tabs[4]:
                        _show_report(result)
                    with tabs[5]:
                        _show_fasta(result)

        except Exception:
            st.error("An unexpected error occurred:")
            st.code(traceback.format_exc(), language="python")

elif not run_full and not run_analyze:
    st.info("👈 Provide a VH sequence in the sidebar and click a run button to begin.")


# ---------------------------------------------------------------------------
# Programmatic launch
# ---------------------------------------------------------------------------


def main() -> None:
    """Programmatically launch the Streamlit app."""
    app_path = Path(__file__).resolve()
    subprocess.run([sys.executable, "-m", "streamlit", "run", str(app_path)], check=True)

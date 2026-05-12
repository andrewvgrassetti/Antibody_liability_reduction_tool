# Antibody Liability Reduction Tool

A Python pipeline for systematically identifying and reducing surface-exposed liabilities (hydrophobic patches, positive charge clusters) in antibody VH sequences while maintaining humanness and structural stability.

## Pipeline Overview

The tool runs a 10-stage pipeline:

1. **IMGT Numbering** – Number residues via ANARCI with IMGT scheme
2. **Surface Exposure** – Classify each position as surface or buried
3. **Liability Detection** – Identify surface-exposed hydrophobic/charged liabilities
4. **Mutation Generation** – Propose humanness-constrained mutations using Abysis frequencies
5. **Combinatorial Expansion** – Build double/triple mutation combinations (Bayesian optimization)
6. **Parent Evaluation** – Score the parent sequence (TAP, DeepSP, OASis)
7. **Mutant Evaluation** – Score all candidate sequences
8. **Scoring** – Composite ranking of candidates
9. **Filtering** – Apply humanness and developability thresholds
10. **Reporting** – Generate HTML reports, FASTA exports, and visualizations

## Installation

```bash
# Clone the repository
git clone https://github.com/andrewvgrassetti/Antibody_liability_reduction_tool.git
cd Antibody_liability_reduction_tool

# Install in development mode
pip install -e ".[dev]"
```

> **Note:** Some dependencies (ANARCI, BioPhi, DeepSP) have complex installation requirements. See their documentation for details.

## Quick Start

### CLI Usage

```bash
# Analyze liabilities only (stages 1-3)
antibody-liability-tool analyze --fasta examples/example_input.fasta

# Run the full pipeline
antibody-liability-tool run --fasta examples/example_input.fasta --output-dir output/

# With a custom config
antibody-liability-tool run -s "EVQLVES..." -c config/default.yaml -o output/
```

### Streamlit Web UI

The tool includes an interactive web interface built with [Streamlit](https://streamlit.io/). Launch it with:

```bash
# Using the installed entry point
antibody-liability-ui

# Or directly with Streamlit
streamlit run src/antibody_liability_tool/app.py
```

The web UI provides:

- **Sequence input** — Upload a FASTA file or paste a raw VH sequence
- **Interactive configuration** — Adjust pipeline parameters (human frequency threshold, combination order, OASis threshold, top-N) via sidebar controls
- **Two run modes** — "Analyze Only" (stages 1–3) for quick liability detection, or "Run Pipeline" for the full 10-stage workflow
- **Live progress** — Real-time stage-by-stage status updates
- **Results tabs** — Summary metrics, liability table, ranked candidates (with CSV download), interactive Plotly visualizations (radar plot, liability map), full HTML report, and FASTA export

### Python API

```python
from antibody_liability_tool.config import load_config
from antibody_liability_tool.numbering.imgt import number_sequence
from antibody_liability_tool.surface.exposure import classify_surface_exposure
from antibody_liability_tool.liabilities.detector import detect_liabilities

config = load_config()
numbered = number_sequence(sequence)
exposure = classify_surface_exposure(numbered)
liabilities = detect_liabilities(numbered, exposure, config)
```

See `examples/run_example.py` for a complete walkthrough.

## Configuration

All pipeline parameters are controlled via `config/default.yaml`. Key sections:

- **liabilities** – Severity thresholds and residue classifications
- **mutations** – Abysis frequency thresholds, allowed substitutions
- **combinatorial** – Max mutation order, Bayesian optimization settings
- **scoring** – Evaluator weights and filtering thresholds
- **output** – Report formats and top-N settings

## Project Structure

```
├── config/
│   └── default.yaml            # Default pipeline configuration
├── data/
│   ├── abysis_vh_frequencies.json  # Human VH germline frequencies
│   └── imgt_surface_exposure.json  # IMGT surface exposure reference
├── examples/
│   ├── example_input.fasta     # Example trastuzumab VH sequence
│   └── run_example.py          # Example pipeline script
├── scripts/
│   └── refresh_abysis_cache.py # Refresh Abysis frequency data
├── src/antibody_liability_tool/
│   ├── numbering/              # IMGT numbering (ANARCI wrapper)
│   ├── surface/                # Surface exposure classification
│   ├── liabilities/            # Liability detection
│   ├── mutations/              # Mutation generation & motif checking
│   ├── evaluators/             # TAP, DeepSP, OASis evaluators
│   ├── optimization/           # Combinatorial expansion & Bayesian opt
│   ├── reporting/              # HTML reports, FASTA export, plots
│   ├── app.py                  # Streamlit web UI
│   ├── cli.py                  # Click CLI
│   ├── config.py               # Configuration loader
│   ├── pipeline.py             # Main pipeline orchestrator
│   └── caching.py              # Disk cache utilities
├── tests/                      # Unit and integration tests
├── pyproject.toml              # Project metadata and dependencies
└── .github/workflows/ci.yml   # CI pipeline
```

## Data Files

- **`data/abysis_vh_frequencies.json`** – Amino acid frequencies at each IMGT position from human VH germline sequences (Abysis database). Used to constrain mutations to human-like substitutions.
- **`data/imgt_surface_exposure.json`** – Reference data for IMGT-position surface exposure classification based on known antibody structures.

## External Services

| Service | Purpose | Required? |
|---------|---------|-----------|
| **TAP** | Therapeutic Antibody Profiler – surface patch analysis | Optional (stages 6-7) |
| **DeepSP** | Deep learning spatial property predictor | Optional (stages 6-7) |
| **BioPhi/OASis** | Humanness scoring via Observed Antibody Space | Optional (stages 6-7) |

Stages 1-4 run fully offline. Stages 5-10 benefit from external evaluators but can run with subsets.

## Development

```bash
# Install dev dependencies
pip install -e ".[dev]"

# Run linter
ruff check src/ tests/

# Run tests (excluding integration/slow)
pytest tests/ -v --tb=short -m "not integration and not slow"

# Run all tests
pytest tests/ -v
```

## Limitations & Future Work

- Abysis frequency scraping requires manual cache refresh (`scripts/refresh_abysis_cache.py`)
- TAP evaluation requires browser automation (Playwright)
- Currently supports VH sequences only; VL/scFv support planned
- Bayesian optimization requires PyTorch/BoTorch for full functionality
- Structural validation (FreeSASA) is an optional extension

## License

MIT

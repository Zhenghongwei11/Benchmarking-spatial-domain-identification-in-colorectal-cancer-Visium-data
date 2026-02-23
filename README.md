# Benchmarking spatial domain identification in CRC Visium data

This repository contains code and derived tables for a benchmark of spatial domain identification methods in colorectal cancer (CRC) 10x Genomics Visium datasets. The benchmark emphasizes pre-specified statistical decision criteria, stability auditing, and transparent reporting.

## Why this study matters
Spatial transcriptomics makes it possible to see how tumor cells, stroma, and immune compartments are organized in situ, but many downstream analyses depend on an upstream “spatial domain” map that is often chosen by eye. In colorectal cancer, domain boundaries can be gradual and mixed, so small analytic choices can change the apparent tissue structure. This project benchmarks a commonly used Bayesian spatial clustering method (BayesSpace) against simple baselines under fixed settings and reports quantitative evidence for domain quality and stability, with a focus on transparent, reproducible decision-making.

Zenodo DOIs:
- Archived release used for reporting (v1.0.4): https://doi.org/10.5281/zenodo.18733963
- Concept DOI (all versions; resolves to latest): https://doi.org/10.5281/zenodo.18733930

## Quick start (reproduce key tables)
Prerequisites: Python (3.x) and R (with `Rscript`) available on PATH.

Run the minimal stages (these scripts will download public GEO data and create a local `.venv` automatically):
- `bash scripts/run_crc_stage2_local.sh`
- `BAYES_INSTALL=1 bash scripts/run_crc_stage3_full_replication.sh`
- `bash scripts/run_crc_stage3d_spatial_ward_baseline.sh`
- `BAYES_INSTALL=1 bash scripts/run_crc_stage5_bayesspace_robustness.sh`
- `BAYES_INSTALL=1 bash scripts/run_crc_stage6_bayesspace_multiseed.sh`
- `BAYES_INSTALL=1 bash scripts/run_crc_stage7_bayesspace_crosscohort.sh`
- `BAYES_INSTALL=1 bash scripts/run_crc_stage8_bayesspace_crosscohort_multiseed.sh`
- `BAYES_INSTALL=1 bash scripts/run_crc_stage9_bayesspace_crosscohort_sample3.sh`
- `bash scripts/run_crc_stage10_bayesspace_rigor_backfill.sh`

Then rebuild the claim-gate table and derived artifacts:
- `Rscript scripts/build_statistical_gate_summary.R`
- `python3 scripts/build_required_artifacts.py`

See `docs/review_bundle/REPRODUCE.md` for details and expected outputs.

### One-click (end-to-end)
To run the full pipeline (tables + figures) with a single command:
- `bash scripts/reproduce_one_click.sh`

## Regenerating figures (optional)
If you want to regenerate publication figures locally, install the Python dependencies and rerun figure scripts:
- `python3 -m venv .venv && source .venv/bin/activate`
- `python -m pip install -r requirements.txt`
- `python scripts/make_publication_figures_v2.py`
- `python scripts/make_supplementary_figures.py`

## Data sources
Public GEO accessions used in this benchmark:
- GSE267401
- GSE311294
- GSE285505

The download URLs and file sizes are recorded in `docs/DATA_MANIFEST.tsv`.

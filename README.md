# Benchmarking spatial domain identification in CRC Visium data

This repository contains code and derived tables for a benchmark of spatial domain identification methods in colorectal cancer (CRC) 10x Genomics Visium datasets. The benchmark emphasizes pre-specified statistical decision criteria, stability auditing, and transparent reporting.

Zenodo DOI (concept): https://doi.org/10.5281/zenodo.18733931

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

See `docs/review_bundle/REVIEWER_GUIDE.md` for details and expected outputs.

### One-click (end-to-end)
To run the full pipeline (tables + figures) with a single command:
- `bash scripts/reproduce_one_click.sh`

## Manuscript and figures
- Main figures: `plots/publication/png/figure1.png` to `figure4.png`
- Supplementary figures: `plots/publication/png/figureS1.png`, `figureS2.png`
Submission materials (manuscript files, cover letters, checklists) are maintained outside the public repository.

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

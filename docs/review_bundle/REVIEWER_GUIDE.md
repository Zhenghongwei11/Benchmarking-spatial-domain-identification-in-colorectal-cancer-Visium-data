# Reviewer guide (how to reproduce key numbers)

## Repository (recommended)
- Repository (HTTPS): https://github.com/Zhenghongwei11/Benchmarking-spatial-domain-identification-in-colorectal-cancer-Visium-data
- Repository (SSH): `git@github.com:Zhenghongwei11/Benchmarking-spatial-domain-identification-in-colorectal-cancer-Visium-data.git`

Clone (choose one):
- `git clone https://github.com/Zhenghongwei11/Benchmarking-spatial-domain-identification-in-colorectal-cancer-Visium-data`
- `git clone git@github.com:Zhenghongwei11/Benchmarking-spatial-domain-identification-in-colorectal-cancer-Visium-data.git`

## Canonical reproducibility bundle (alternative)
If you received `crc_spatial_benchmark_review_bundle.zip`, verify integrity using the provided `CHECKSUMS.sha256` manifest before running any commands.

## Minimal reproduction (tables)
### One-click (recommended)
Run the full pipeline end-to-end:
- `bash scripts/reproduce_one_click.sh`

Environment controls:
- `BAYES_INSTALL=1` (default): allow BayesSpace installation if missing.
- `RUN_FIGURES=0`: tables-only reproduction.

### Step-by-step (tables only)
1. Ensure Python and R are available.
2. Run the stages (baseline + BayesSpace where applicable):
   - `bash scripts/run_crc_stage2_local.sh`
   - `BAYES_INSTALL=1 bash scripts/run_crc_stage3_full_replication.sh`
   - `bash scripts/run_crc_stage3d_spatial_ward_baseline.sh`
   - `BAYES_INSTALL=1 bash scripts/run_crc_stage5_bayesspace_robustness.sh`
   - `BAYES_INSTALL=1 bash scripts/run_crc_stage6_bayesspace_multiseed.sh`
   - `BAYES_INSTALL=1 bash scripts/run_crc_stage7_bayesspace_crosscohort.sh`
   - `BAYES_INSTALL=1 bash scripts/run_crc_stage8_bayesspace_crosscohort_multiseed.sh`
   - `BAYES_INSTALL=1 bash scripts/run_crc_stage9_bayesspace_crosscohort_sample3.sh`
   - `bash scripts/run_crc_stage10_bayesspace_rigor_backfill.sh`
3. Rebuild the claim-gate table:
   - `Rscript scripts/build_statistical_gate_summary.R`
4. Rebuild derived artifacts:
    - `python3 scripts/build_required_artifacts.py`

## Regenerate figures (optional)
To regenerate the publication figures after running the stages above:
1. Install Python plotting dependencies:
   - `python3 -m venv .venv && source .venv/bin/activate`
   - `python -m pip install -r requirements.txt`
2. Recreate figures:
   - `python scripts/make_publication_figures_v2.py`
   - `python scripts/make_supplementary_figures.py`

## What to check
- Claim gates: `results/benchmarks/statistical_gate_summary.tsv`
- Effect sizes with CIs: `results/effect_sizes/claim_effects.tsv`
- Dataset coverage: `results/dataset_summary.tsv`
- Replication tables: `results/replication/*.tsv`

## Notes
- Bayesian MCMC methods can be slow; long runtimes are expected and should be recorded in the benchmark tables.
- The bundle intentionally excludes raw data; all datasets are public and can be downloaded from GEO.
- If your environment pins SSH host keys, validate GitHub host keys via official GitHub documentation before trusting any fingerprint. (Record provided by the authors: `SHA256:qk4DtsWl/9BBnPK5csvpuHeooFvh/WTLJMBmFdpoBoA`.)

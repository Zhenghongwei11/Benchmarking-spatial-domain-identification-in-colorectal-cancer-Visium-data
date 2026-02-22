#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
cd "${ROOT_DIR}"

# Controls
# - Set BAYES_INSTALL=1 to allow the BayesSpace installer to run if needed.
# - Set RUN_FIGURES=0 to skip regenerating figures (tables-only reproduction).
# - This script intentionally does not build submission documents (DOCX/PDF).
BAYES_INSTALL="${BAYES_INSTALL:-1}"
RUN_FIGURES="${RUN_FIGURES:-1}"

echo "[1/6] Running baseline stages (downloads public GEO data)"
bash scripts/run_crc_stage2_local.sh

echo "[2/6] Running full replication baselines"
BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage3_full_replication.sh

echo "[3/6] Running spatial Ward baseline (M2) sensitivity baseline"
bash scripts/run_crc_stage3d_spatial_ward_baseline.sh

echo "[4/6] Running BayesSpace stages (robustness + multiseed + cross-cohort samples)"
BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage5_bayesspace_robustness.sh
BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage6_bayesspace_multiseed.sh
BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage7_bayesspace_crosscohort.sh
BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage8_bayesspace_crosscohort_multiseed.sh
BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage9_bayesspace_crosscohort_sample3.sh

echo "[5/6] Backfilling rigor rows and rebuilding locked tables"
bash scripts/run_crc_stage10_bayesspace_rigor_backfill.sh
Rscript scripts/build_statistical_gate_summary.R
python3 scripts/build_required_artifacts.py

if [[ "${RUN_FIGURES}" == "1" ]]; then
  echo "[6/6] Regenerating figures"
  if [[ -d ".venv" ]]; then
    # Stage scripts create the venv; activate it for figure generation.
    # shellcheck disable=SC1091
    source ".venv/bin/activate"
  fi
  python scripts/make_publication_figures_v2.py
  python scripts/make_supplementary_figures.py
else
  echo "[6/6] Skipping figure regeneration (RUN_FIGURES=0)"
fi

echo "Done."

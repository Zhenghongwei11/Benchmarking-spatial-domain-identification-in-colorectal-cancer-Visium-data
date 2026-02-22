#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SEEDS="${SEEDS:-11,23}"
K_GRID="${K_GRID:-4,6}"
NOTE_PREFIX="${NOTE_PREFIX:-stage10-rigor-backfill}"
STAGE_ID="${STAGE_ID:-stage10-bayesspace-rigor-backfill}"
DATASET_FILTER="${DATASET_FILTER:-}"
SAMPLE_FILTER="${SAMPLE_FILTER:-GSM9322957_TR11_206,GSM9322958_TR11_16184,GSM9322959_TR11_18105,GSM8265211_CTC21P,GSM8265212_CTC21M,GSM8265213_CTC17P}"
BAYES_INSTALL="${BAYES_INSTALL:-0}"

cd "${ROOT_DIR}"

if [[ "${BAYES_INSTALL}" == "1" ]]; then
  Rscript scripts/check_or_install_bayesspace.R >/dev/null
else
  Rscript scripts/check_or_install_bayesspace.R --no-install >/dev/null
fi

SEEDS="${SEEDS}" \
K_GRID="${K_GRID}" \
NOTE_PREFIX="${NOTE_PREFIX}" \
STAGE_ID="${STAGE_ID}" \
DATASET_FILTER="${DATASET_FILTER}" \
SAMPLE_FILTER="${SAMPLE_FILTER}" \
bash scripts/run_bayesspace_rigor_backfill.sh

echo "Stage-10 completed with BayesSpace rigor backfill."

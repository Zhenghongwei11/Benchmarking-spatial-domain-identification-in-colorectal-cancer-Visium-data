#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SUMMARY_TSV="${ROOT_DIR}/results/benchmarks/stage4_run_summary.tsv"
BAYES_INSTALL="${BAYES_INSTALL:-1}"

cd "${ROOT_DIR}"

count_bayespace_rows() {
  if [[ ! -f results/benchmarks/method_benchmark.tsv ]]; then
    echo "0"
    return
  fi
  tail -n +2 results/benchmarks/method_benchmark.tsv | awk -F'\t' '$3=="BayesSpace"{n++} END{print n+0}'
}

before_count="$(count_bayespace_rows)"
started_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
stage3c_status="success"

if ! STAGES=stage3c BAYES_INSTALL="${BAYES_INSTALL}" bash scripts/run_crc_stage3_full_replication.sh; then
  stage3c_status="failed"
fi

after_count="$(count_bayespace_rows)"
added_count=$((after_count - before_count))
finished_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"

printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
  "stage4-bayesspace" "${started_utc}" "${finished_utc}" "${BAYES_INSTALL}" \
  "${stage3c_status}" "${before_count}" "${after_count}" "${added_count}" \
  "stage4-bayesspace-only-run" >> "${SUMMARY_TSV}"

if [[ "${stage3c_status}" != "success" || "${added_count}" -le 0 ]]; then
  echo "Stage-4 completed with BayesSpace not added. Check failure_log.tsv and stage3_run_summary.tsv"
  exit 1
fi

echo "Stage-4 completed with BayesSpace rows added: ${added_count}"

#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SUMMARY_TSV="${ROOT_DIR}/results/benchmarks/stage5_run_summary.tsv"
DATASET_ID="${DATASET_ID:-GSE311294}"
DATASET_ROOT="${DATASET_ROOT:-data/raw/${DATASET_ID}/extracted}"
SAMPLE_ID="${SAMPLE_ID:-GSM9322958_TR11_16184}"
K_GRID="${K_GRID:-4,6}"
SEED="${SEED:-11}"
BAYES_INSTALL="${BAYES_INSTALL:-0}"
NOTE_PREFIX="${NOTE_PREFIX:-stage5-bayesspace-robustness}"

cd "${ROOT_DIR}"

count_bayespace_rows() {
  if [[ ! -f results/benchmarks/method_benchmark.tsv ]]; then
    echo "0"
    return
  fi
  tail -n +2 results/benchmarks/method_benchmark.tsv | awk -F'\t' '$3=="BayesSpace"{n++} END{print n+0}'
}

append_failure_log() {
  local dataset_id="$1"
  local sample_id="$2"
  local err_msg="$3"
  local now
  now="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  err_msg="$(echo "${err_msg}" | tr '\t' ' ' | tr '\n' ' ' | cut -c1-500)"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${dataset_id}" "${sample_id}" "BayesSpace" "bayesspace_default" "K${K_GRID}" "${K_GRID}" "${SEED}" \
    "bayesspace_error" "${err_msg}" "" "" "" "" "${now}" "${NOTE_PREFIX}" \
    >> "${ROOT_DIR}/results/benchmarks/failure_log.tsv"
}

append_stage5_summary() {
  local stage_id="$1"
  local started="$2"
  local finished="$3"
  local status="$4"
  local before_count="$5"
  local after_count="$6"
  local added_count="$7"
  local note="$8"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${stage_id}" "${DATASET_ID}" "${SAMPLE_ID}" "${K_GRID}" "${SEED}" "${BAYES_INSTALL}" \
    "${started}" "${finished}" "${status}" "${before_count}" "${after_count}" "${added_count}" "${note}" \
    >> "${SUMMARY_TSV}"
}

before_count="$(count_bayespace_rows)"
started_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
status="success"
tmp_tsv="$(mktemp)"
tmp_log="$(mktemp)"

if [[ "${BAYES_INSTALL}" == "1" ]]; then
  R_CHECK_CMD=(Rscript scripts/check_or_install_bayesspace.R)
else
  R_CHECK_CMD=(Rscript scripts/check_or_install_bayesspace.R --no-install)
fi

if ! "${R_CHECK_CMD[@]}" >"${tmp_log}" 2>&1; then
  status="failed"
  append_failure_log "${DATASET_ID}" "${SAMPLE_ID}" "$(tail -n 40 "${tmp_log}")"
else
  if ! Rscript scripts/run_bayesspace_baseline.R \
    --dataset-id "${DATASET_ID}" \
    --dataset-root "${DATASET_ROOT}" \
    --sample-id "${SAMPLE_ID}" \
    --k-grid "${K_GRID}" \
    --seed "${SEED}" \
    --note "${NOTE_PREFIX}" \
    --output-tsv "${tmp_tsv}" >"${tmp_log}" 2>&1; then
    status="failed"
    append_failure_log "${DATASET_ID}" "${SAMPLE_ID}" "$(tail -n 40 "${tmp_log}")"
  else
    tail -n +2 "${tmp_tsv}" >> "${ROOT_DIR}/results/benchmarks/method_benchmark.tsv"
    tail -n +2 "${tmp_tsv}" | sed 's/^/Fig2B\t/' >> "${ROOT_DIR}/results/figures/fig2_benchmark_summary.tsv"
    while IFS=$'\t' read -r ds sid mid mf prep pset kval seed_cnt s_ari s_iqr sp sp_iqr mk mk_iqr wt mem fail notes; do
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "Fig4B" "${ds}" "${sid}" "${mid}" "${kval}" "success" "${wt}" "${mem}" "1" "0" "local-first" "${notes}" \
        >> "${ROOT_DIR}/results/figures/fig4_compute_and_guidance.tsv"
    done < <(tail -n +2 "${tmp_tsv}")
  fi
fi

after_count="$(count_bayespace_rows)"
added_count=$((after_count - before_count))
finished_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
append_stage5_summary "stage5-bayesspace-robustness" "${started_utc}" "${finished_utc}" "${status}" "${before_count}" "${after_count}" "${added_count}" "${NOTE_PREFIX}"

if [[ "${status}" != "success" || "${added_count}" -le 0 ]]; then
  echo "Stage-5 finished without BayesSpace additions. Check stage5_run_summary.tsv and failure_log.tsv"
  exit 1
fi

echo "Stage-5 completed with BayesSpace rows added: ${added_count}"

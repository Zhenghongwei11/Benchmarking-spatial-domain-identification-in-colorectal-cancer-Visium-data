#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
SUMMARY_TSV="${ROOT_DIR}/results/benchmarks/stage6_run_summary.tsv"
DATASET_ID="${DATASET_ID:-GSE311294}"
DATASET_ROOT="${DATASET_ROOT:-data/raw/${DATASET_ID}/extracted}"
SAMPLE_ID="${SAMPLE_ID:-GSM9322959_TR11_18105}"
SEEDS="${SEEDS:-23,37}"
K_GRID="${K_GRID:-4,6}"
BAYES_INSTALL="${BAYES_INSTALL:-0}"
NOTE_PREFIX="${NOTE_PREFIX:-stage6-bayesspace-multiseed}"

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
  local seed="$3"
  local err_msg="$4"
  local now
  now="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  err_msg="$(echo "${err_msg}" | tr '\t' ' ' | tr '\n' ' ' | cut -c1-500)"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${dataset_id}" "${sample_id}" "BayesSpace" "bayesspace_default" "K${K_GRID}" "${K_GRID}" "${seed}" \
    "bayesspace_error" "${err_msg}" "" "" "" "" "${now}" "${NOTE_PREFIX}" \
    >> "${ROOT_DIR}/results/benchmarks/failure_log.tsv"
}

append_stage6_summary() {
  local stage_id="$1"
  local dataset_id="$2"
  local sample_id="$3"
  local seed="$4"
  local started="$5"
  local finished="$6"
  local status="$7"
  local before_count="$8"
  local after_count="$9"
  local added_count="${10}"
  local note="${11}"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${stage_id}" "${dataset_id}" "${sample_id}" "${K_GRID}" "${seed}" "${BAYES_INSTALL}" \
    "${started}" "${finished}" "${status}" "${before_count}" "${after_count}" "${added_count}" "${SEEDS}" "${note}" \
    >> "${SUMMARY_TSV}"
}

if [[ "${BAYES_INSTALL}" == "1" ]]; then
  R_CHECK_CMD=(Rscript scripts/check_or_install_bayesspace.R)
else
  R_CHECK_CMD=(Rscript scripts/check_or_install_bayesspace.R --no-install)
fi
"${R_CHECK_CMD[@]}" >/dev/null

before_total="$(count_bayespace_rows)"
overall_status="success"

IFS=',' read -r -a seed_list <<< "${SEEDS}"
for seed in "${seed_list[@]}"; do
  seed="$(echo "${seed}" | xargs)"
  if [[ -z "${seed}" ]]; then
    continue
  fi

  before_count="$(count_bayespace_rows)"
  started_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  status="success"
  note="${NOTE_PREFIX}-seed${seed}"
  tmp_tsv="$(mktemp)"
  tmp_log="$(mktemp)"

  if ! Rscript scripts/run_bayesspace_baseline.R \
    --dataset-id "${DATASET_ID}" \
    --dataset-root "${DATASET_ROOT}" \
    --sample-id "${SAMPLE_ID}" \
    --k-grid "${K_GRID}" \
    --seed "${seed}" \
    --note "${note}" \
    --output-tsv "${tmp_tsv}" >"${tmp_log}" 2>&1; then
    status="failed"
    overall_status="failed"
    append_failure_log "${DATASET_ID}" "${SAMPLE_ID}" "${seed}" "$(tail -n 40 "${tmp_log}")"
  else
    tail -n +2 "${tmp_tsv}" >> "${ROOT_DIR}/results/benchmarks/method_benchmark.tsv"
    tail -n +2 "${tmp_tsv}" | sed 's/^/Fig2C\t/' >> "${ROOT_DIR}/results/figures/fig2_benchmark_summary.tsv"
    while IFS=$'\t' read -r ds sid mid mf prep pset kval seed_cnt s_ari s_iqr sp sp_iqr mk mk_iqr wt mem fail notes; do
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "Fig4C" "${ds}" "${sid}" "${mid}" "${kval}" "success" "${wt}" "${mem}" "1" "0" "local-first" "${notes}" \
        >> "${ROOT_DIR}/results/figures/fig4_compute_and_guidance.tsv"
    done < <(tail -n +2 "${tmp_tsv}")
  fi

  after_count="$(count_bayespace_rows)"
  added_count=$((after_count - before_count))
  finished_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  append_stage6_summary "stage6-bayesspace-multiseed" "${DATASET_ID}" "${SAMPLE_ID}" "${seed}" "${started_utc}" "${finished_utc}" "${status}" "${before_count}" "${after_count}" "${added_count}" "${note}"
done

after_total="$(count_bayespace_rows)"
added_total=$((after_total - before_total))

if [[ "${overall_status}" != "success" || "${added_total}" -le 0 ]]; then
  echo "Stage-6 finished without full BayesSpace multiseed success. Check stage6_run_summary.tsv and failure_log.tsv"
  exit 1
fi

echo "Stage-6 completed with BayesSpace rows added: ${added_total}"

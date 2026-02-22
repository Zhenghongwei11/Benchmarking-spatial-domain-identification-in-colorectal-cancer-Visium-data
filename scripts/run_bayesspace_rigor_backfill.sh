#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
METHOD_TSV="${ROOT_DIR}/results/benchmarks/method_benchmark.tsv"
BACKFILL_TSV="${ROOT_DIR}/results/benchmarks/bayesspace_rigor_backfill.tsv"
RUNLOG_TSV="${ROOT_DIR}/results/benchmarks/bayesspace_rigor_backfill_runs.tsv"
STAGE_SUMMARY_TSV="${ROOT_DIR}/results/benchmarks/stage10_run_summary.tsv"
FAILURE_TSV="${ROOT_DIR}/results/benchmarks/failure_log.tsv"
FIG2_TSV="${ROOT_DIR}/results/figures/fig2_benchmark_summary.tsv"
FIG4_TSV="${ROOT_DIR}/results/figures/fig4_compute_and_guidance.tsv"

SEEDS="${SEEDS:-11,23}"
K_GRID="${K_GRID:-4,6}"
NOTE_PREFIX="${NOTE_PREFIX:-rigor-backfill}"
STAGE_ID="${STAGE_ID:-stage10-bayesspace-rigor-backfill}"
DATASET_FILTER="${DATASET_FILTER:-}"
SAMPLE_FILTER="${SAMPLE_FILTER:-}"

cd "${ROOT_DIR}"

if [[ ! -f "${METHOD_TSV}" ]]; then
  echo "Missing ${METHOD_TSV}"
  exit 1
fi

if [[ ! -f "${BACKFILL_TSV}" ]]; then
  head -n 1 "${METHOD_TSV}" > "${BACKFILL_TSV}"
fi

if [[ ! -f "${RUNLOG_TSV}" ]]; then
  printf "dataset_id\tsample_id\tk_grid\tseeds\tstarted_utc\tfinished_utc\tstatus\trows_added\tnote\n" > "${RUNLOG_TSV}"
fi

if [[ ! -f "${STAGE_SUMMARY_TSV}" ]]; then
  printf "stage_id\tdataset_id\tsample_id\tk_grid\tseeds\tstarted_utc\tfinished_utc\tstatus\trows_added\tnote\n" > "${STAGE_SUMMARY_TSV}"
fi

count_bayes_rows() {
  tail -n +2 "${METHOD_TSV}" | awk -F'\t' '$3=="BayesSpace"{n++} END{print n+0}'
}

append_failure_log() {
  local dataset_id="$1"
  local sample_id="$2"
  local err_msg="$3"
  local now
  now="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  err_msg="$(echo "${err_msg}" | tr '\t' ' ' | tr '\n' ' ' | cut -c1-500)"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${dataset_id}" "${sample_id}" "BayesSpace" "bayesspace_default" "K${K_GRID}" "${K_GRID}" "${SEEDS}" \
    "bayesspace_error" "${err_msg}" "" "" "" "" "${now}" "${NOTE_PREFIX}" \
    >> "${FAILURE_TSV}"
}

sample_allowed() {
  local sample_id="$1"
  if [[ -z "${SAMPLE_FILTER}" ]]; then
    return 0
  fi
  IFS=',' read -r -a sample_list <<< "${SAMPLE_FILTER}"
  for target in "${sample_list[@]}"; do
    target="$(echo "${target}" | xargs)"
    if [[ "${sample_id}" == "${target}" ]]; then
      return 0
    fi
  done
  return 1
}

target_count="$(
  tail -n +2 "${METHOD_TSV}" \
    | awk -F'\t' '$3=="BayesSpace" && $1!="" && $2!=""{print $1 "\t" $2}' \
    | sort -u \
    | wc -l \
    | tr -d ' '
)"

if [[ "${target_count}" == "0" ]]; then
  echo "No BayesSpace rows found in ${METHOD_TSV}"
  exit 1
fi

echo "Backfill targets discovered: ${target_count}"
overall_status="success"
executed_count=0

while IFS=$'\t' read -r dataset_id sample_id; do

  if [[ -n "${DATASET_FILTER}" && "${dataset_id}" != "${DATASET_FILTER}" ]]; then
    continue
  fi
  if ! sample_allowed "${sample_id}"; then
    continue
  fi
  executed_count=$((executed_count + 1))

  dataset_root="data/raw/${dataset_id}/extracted"
  note="${NOTE_PREFIX}-${dataset_id}-${sample_id}"
  started_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  before_count="$(count_bayes_rows)"
  status="success"
  tmp_tsv="$(mktemp)"
  tmp_log="$(mktemp)"

  echo "[${started_utc}] Running ${dataset_id} / ${sample_id} (K=${K_GRID}, seeds=${SEEDS})"
  if ! Rscript scripts/run_bayesspace_baseline.R \
    --dataset-id "${dataset_id}" \
    --dataset-root "${dataset_root}" \
    --sample-id "${sample_id}" \
    --k-grid "${K_GRID}" \
    --seeds "${SEEDS}" \
    --note "${note}" \
    --output-tsv "${tmp_tsv}" >"${tmp_log}" 2>&1; then
    status="failed"
    append_failure_log "${dataset_id}" "${sample_id}" "$(tail -n 40 "${tmp_log}")"
  else
    tail -n +2 "${tmp_tsv}" >> "${METHOD_TSV}"
    tail -n +2 "${tmp_tsv}" >> "${BACKFILL_TSV}"
    tail -n +2 "${tmp_tsv}" | sed 's/^/Fig2R\t/' >> "${FIG2_TSV}"
    while IFS=$'\t' read -r ds sid mid mf prep pset kval seed_cnt s_ari s_iqr sp sp_iqr mk mk_iqr wt mem fail notes; do
      printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
        "Fig4R" "${ds}" "${sid}" "${mid}" "${kval}" "success" "${wt}" "${mem}" "1" "0" "local-first" "${notes}" \
        >> "${FIG4_TSV}"
    done < <(tail -n +2 "${tmp_tsv}")
  fi

  after_count="$(count_bayes_rows)"
  rows_added=$((after_count - before_count))
  finished_utc="$(date -u +"%Y-%m-%dT%H:%M:%SZ")"
  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${dataset_id}" "${sample_id}" "${K_GRID}" "${SEEDS}" "${started_utc}" "${finished_utc}" "${status}" "${rows_added}" "${note}" \
    >> "${RUNLOG_TSV}"

  if [[ "${status}" == "failed" ]]; then
    overall_status="failed"
    echo "[${finished_utc}] FAILED ${dataset_id} / ${sample_id}"
  else
    echo "[${finished_utc}] OK ${dataset_id} / ${sample_id}, rows_added=${rows_added}"
  fi

  printf "%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\t%s\n" \
    "${STAGE_ID}" "${dataset_id}" "${sample_id}" "${K_GRID}" "${SEEDS}" "${started_utc}" "${finished_utc}" "${status}" "${rows_added}" "${note}" \
    >> "${STAGE_SUMMARY_TSV}"
done < <(
  tail -n +2 "${METHOD_TSV}" \
    | awk -F'\t' '$3=="BayesSpace" && $1!="" && $2!=""{print $1 "\t" $2}' \
    | sort -u
)

if [[ "${executed_count}" -eq 0 ]]; then
  echo "No targets matched filters."
  exit 1
fi

echo "Backfill run finished. executed=${executed_count} status=${overall_status}"
if [[ "${overall_status}" != "success" ]]; then
  exit 1
fi

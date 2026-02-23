#!/usr/bin/env bash
set -euo pipefail

ROOT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
OUT_DIR="${ROOT_DIR}/docs/review_bundle"
ZIP_NAME="${ZIP_NAME:-crc_spatial_benchmark_review_bundle.zip}"

mkdir -p "${OUT_DIR}"

cat > "${OUT_DIR}/POLICY.md" <<'MD'
# Public review bundle policy (canonical zip)

## Purpose
Provide a single, clean reproducibility bundle suitable for peer review and public release.

## Canonical bundle rule
- There is exactly one canonical ZIP in this folder.
- The journal ZIP and GitHub release ZIP should be identical (verified by checksums).

## Included (high level)
- `scripts/`: analysis entrypoints and helper scripts
- `results/`: benchmark tables and derived summary tables
- `docs/`: protocol and reporting artifacts needed to interpret/reproduce results (excluding submission manuscripts)
- `docs/audit_runs/`: lightweight run provenance (environment + checksums), when safe

## Excluded (by default)
- Submission-only materials: `docs/submissions/`, `docs/manuscript/`, cover letters, checklists tied to a specific submission UI
- Internal scaffolding not needed for reproduction (e.g., local spec tooling and agent/config directories)
- Raw data and large intermediates: `data/` (reviewers can download public data separately)
- Local environments/caches: `.venv/`, `__pycache__/`, OS/editor metadata

## Rationale
Reviewers should see the science and the reproducibility artifacts, not internal spec tooling or submission packaging.
MD

cat > "${OUT_DIR}/REPRODUCE.md" <<'MD'
# Reproduce key numbers (minimal instructions)

## Minimal reproduction (tables)
### One-click (recommended)
Run the full pipeline end-to-end:
- `bash scripts/reproduce_one_click.sh`

Environment controls:
- `BAYES_INSTALL=1` (default): allow BayesSpace installation if missing.
- `RUN_FIGURES=0`: tables-only reproduction.

### Step-by-step (tables only)
1. Ensure Python and R are available.
2. Run the local stages (baseline + BayesSpace where applicable):
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

## What to check
- Claim gates: `results/benchmarks/statistical_gate_summary.tsv`
- Effect sizes with CIs: `results/effect_sizes/claim_effects.tsv`
- Dataset coverage: `results/dataset_summary.tsv`
- Replication tables: `results/replication/*.tsv`

## Notes
- Bayesian MCMC methods can be slow; long runtimes are expected and should be recorded in the benchmark tables.
- The bundle intentionally excludes raw data; all datasets are public and can be downloaded from GEO.
MD

python3 - <<'PY' "${ROOT_DIR}" "${OUT_DIR}"
import os
import sys
from pathlib import Path

root = Path(sys.argv[1]).resolve()
out_dir = Path(sys.argv[2]).resolve()

include_roots = [
    root / "scripts",
    root / "results",
]

include_docs_files = [
    root / "docs" / "STATISTICAL_DECISION_RULES.md",
    root / "docs" / "SOURCE_DATA_MAP.tsv",
    root / "docs" / "FIGURE_PROVENANCE.tsv",
    root / "docs" / "DATA_MANIFEST.tsv",
]

optional_dirs = []

exclude_dir_prefixes = [
    root / "data",
    root / ".venv",
    root / "docs" / "submissions",
    root / "docs" / "manuscript",
    root / "docs" / "review_bundle",
]

def is_excluded(path: Path) -> bool:
    for prefix in exclude_dir_prefixes:
        try:
            path.resolve().relative_to(prefix.resolve())
            return True
        except Exception:
            continue
    return False

def collect_files() -> list[Path]:
    files: list[Path] = []
    for base in include_roots:
        if not base.exists():
            continue
        for p in base.rglob("*"):
            if not p.is_file():
                continue
            if is_excluded(p):
                continue
            if p.name in [".DS_Store"]:
                continue
            if p.suffix in [".zip"]:
                continue
            files.append(p)
    for p in include_docs_files:
        if p.exists() and p.is_file():
            if not is_excluded(p):
                files.append(p)
    for d in optional_dirs:
        if not d.exists():
            continue
        for p in d.rglob("*"):
            if p.is_file() and not is_excluded(p):
                files.append(p)
    # De-dup and sort by repo-relative path
    uniq = {}
    for p in files:
        rel = p.resolve().relative_to(root)
        uniq[str(rel)] = rel
    return [root / uniq[k] for k in sorted(uniq.keys())]

files = collect_files()
filelist_path = out_dir / "FILELIST.txt"
with filelist_path.open("w", encoding="utf-8") as h:
    for p in files:
        h.write(str(p.resolve().relative_to(root)) + "\n")

print(f"file_count={len(files)}")
PY

FILELIST="${OUT_DIR}/FILELIST.txt"
ZIP_PATH="${OUT_DIR}/${ZIP_NAME}"

rm -f "${ZIP_PATH}"

(cd "${ROOT_DIR}" && zip -q -@ "${ZIP_PATH}" < "${FILELIST}")

(
  cd "${OUT_DIR}"
  shasum -a 256 "${ZIP_NAME}" > CHECKSUMS.sha256
)

echo "Wrote ${ZIP_PATH}"

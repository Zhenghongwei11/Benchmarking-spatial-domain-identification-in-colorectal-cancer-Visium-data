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

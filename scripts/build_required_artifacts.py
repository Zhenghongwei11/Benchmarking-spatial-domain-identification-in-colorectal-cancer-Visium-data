#!/usr/bin/env python3

from __future__ import annotations

import csv
import os
import re
from collections import defaultdict
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Iterable


RE_INT = re.compile(r"(-?\d+)")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with path.open("r", encoding="utf-8", newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return [dict(row) for row in reader]


def write_tsv(path: Path, rows: Iterable[dict[str, Any]], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in fieldnames})


def parse_int_maybe(value: str) -> int | None:
    if value is None:
        return None
    match = RE_INT.search(str(value))
    if not match:
        return None
    try:
        return int(match.group(1))
    except Exception:
        return None


def count_extracted_samples(dataset_root: Path) -> int:
    if not dataset_root.exists():
        return 0
    # Flat mtx format
    prefixes = set()
    for path in dataset_root.glob("*_matrix.mtx.gz"):
        prefixes.add(path.name.replace("_matrix.mtx.gz", ""))
    if prefixes:
        return len(prefixes)
    # 10x h5 format
    h5s = list(dataset_root.glob("*_filtered_feature_bc_matrix.h5"))
    if h5s:
        return len(h5s)
    return 0


@dataclass(frozen=True)
class MethodKey:
    dataset_id: str
    sample_id: str
    method_id: str
    K: str


def main() -> int:
    repo_root = Path(__file__).resolve().parent.parent
    docs_dir = repo_root / "docs"
    results_dir = repo_root / "results"

    dataset_landscape_path = docs_dir / "DATASET_LANDSCAPE.tsv"
    method_benchmark_path = results_dir / "benchmarks" / "method_benchmark.tsv"
    gate_summary_path = results_dir / "benchmarks" / "statistical_gate_summary.tsv"

    if not dataset_landscape_path.exists():
        raise FileNotFoundError(dataset_landscape_path)
    if not method_benchmark_path.exists():
        raise FileNotFoundError(method_benchmark_path)
    if not gate_summary_path.exists():
        raise FileNotFoundError(gate_summary_path)

    dataset_landscape = read_tsv(dataset_landscape_path)
    method_benchmark = read_tsv(method_benchmark_path)
    gate_summary = read_tsv(gate_summary_path)

    # -------------------------
    # results/dataset_summary.tsv
    # -------------------------
    mb_by_dataset_samples: dict[str, set[str]] = defaultdict(set)
    mb_by_dataset_methods: dict[str, set[str]] = defaultdict(set)
    mb_by_dataset_bayespace_samples: dict[str, set[str]] = defaultdict(set)
    mb_by_dataset_bayespace_rows: dict[str, int] = defaultdict(int)

    for row in method_benchmark:
        ds = row.get("dataset_id", "").strip()
        sid = row.get("sample_id", "").strip()
        mid = row.get("method_id", "").strip()
        if not ds or not sid or not mid:
            continue
        mb_by_dataset_samples[ds].add(sid)
        mb_by_dataset_methods[ds].add(mid)
        if mid == "BayesSpace":
            mb_by_dataset_bayespace_samples[ds].add(sid)
            mb_by_dataset_bayespace_rows[ds] += 1

    ds_rows: list[dict[str, Any]] = []
    for row in dataset_landscape:
        ds = (row.get("dataset_id_or_name") or "").strip()
        if not ds.startswith("GSE"):
            continue
        extracted_root = repo_root / "data" / "raw" / ds / "extracted"
        samples_on_disk = count_extracted_samples(extracted_root)

        role = (row.get("intended_role") or "").strip()
        modality = (row.get("modality") or "").strip()
        sample_size_declared = parse_int_maybe(row.get("sample_size") or "")
        methods_covered = sorted(mb_by_dataset_methods.get(ds, set()))
        bayespace_samples = sorted(mb_by_dataset_bayespace_samples.get(ds, set()))

        has_flat_mtx = bool(list(extracted_root.glob("*_matrix.mtx.gz")))
        has_h5 = bool(list(extracted_root.glob("*_filtered_feature_bc_matrix.h5")))
        format_hint = "mtx" if has_flat_mtx else ("h5" if has_h5 else "unknown")

        ds_rows.append(
            {
                "dataset_id": ds,
                "intended_role": role,
                "modality": modality,
                "sample_size_declared": sample_size_declared if sample_size_declared is not None else "",
                "samples_on_disk": samples_on_disk,
                "format_hint": format_hint,
                "samples_in_method_benchmark": len(mb_by_dataset_samples.get(ds, set())),
                "methods_covered": ",".join(methods_covered),
                "bayesspace_samples_covered": len(bayespace_samples),
                "bayesspace_rows": mb_by_dataset_bayespace_rows.get(ds, 0),
                "notes": (row.get("notes") or "").strip(),
            }
        )

    write_tsv(
        results_dir / "dataset_summary.tsv",
        ds_rows,
        fieldnames=[
            "dataset_id",
            "intended_role",
            "modality",
            "sample_size_declared",
            "samples_on_disk",
            "format_hint",
            "samples_in_method_benchmark",
            "methods_covered",
            "bayesspace_samples_covered",
            "bayesspace_rows",
            "notes",
        ],
    )

    # -------------------------
    # results/effect_sizes/claim_effects.tsv
    # -------------------------
    effect_rows: list[dict[str, Any]] = []
    for row in gate_summary:
        note = (row.get("notes") or "").strip()
        n_val = None
        for key in ["n_samples", "n_sample_k", "n"]:
            match = re.search(rf"{key}=(\d+)", note)
            if match:
                n_val = int(match.group(1))
                break
        effect_rows.append(
            {
                "claim_id": row.get("claim_id", ""),
                "dataset_id": "meta",
                "outcome": row.get("metric_id", ""),
                "model": row.get("comparison_id", ""),
                "analysis_unit": row.get("analysis_unit", ""),
                "test_name": row.get("test_name", ""),
                "effect_type": row.get("effect_size_name", ""),
                "effect": row.get("effect_size_value", ""),
                "ci_lower": row.get("effect_size_ci_lower", ""),
                "ci_upper": row.get("effect_size_ci_upper", ""),
                "pvalue": row.get("pvalue", ""),
                "fdr": row.get("fdr", ""),
                "n": n_val if n_val is not None else "",
                "overall_gate_status": row.get("overall_gate_status", ""),
                "support_tier": row.get("support_tier", ""),
                "notes": note,
            }
        )

    write_tsv(
        results_dir / "effect_sizes" / "claim_effects.tsv",
        effect_rows,
        fieldnames=[
            "claim_id",
            "dataset_id",
            "outcome",
            "model",
            "analysis_unit",
            "test_name",
            "effect_type",
            "effect",
            "ci_lower",
            "ci_upper",
            "pvalue",
            "fdr",
            "n",
            "overall_gate_status",
            "support_tier",
            "notes",
        ],
    )

    # -------------------------
    # results/replication/* (lightweight tables)
    # -------------------------
    index: dict[MethodKey, dict[str, str]] = {}
    for row in method_benchmark:
        key = MethodKey(
            dataset_id=row.get("dataset_id", ""),
            sample_id=row.get("sample_id", ""),
            method_id=row.get("method_id", ""),
            K=row.get("K", ""),
        )
        if not key.dataset_id or not key.sample_id or not key.method_id or not key.K:
            continue
        index[key] = row

    bayespace_keys = [k for k in index.keys() if k.method_id == "BayesSpace"]
    deltas: list[dict[str, Any]] = []
    for k in bayespace_keys:
        bs = index[k]
        for baseline in ["M0_expr_kmeans", "M1_spatial_concat_kmeans", "M2_spatial_ward"]:
            base_key = MethodKey(k.dataset_id, k.sample_id, baseline, k.K)
            if base_key not in index:
                continue
            base = index[base_key]
            def fnum(val: str) -> float | None:
                try:
                    return float(val)
                except Exception:
                    return None

            bs_sp = fnum(bs.get("spatial_coherence_median", ""))
            bs_mk = fnum(bs.get("marker_coherence_median", ""))
            base_sp = fnum(base.get("spatial_coherence_median", ""))
            base_mk = fnum(base.get("marker_coherence_median", ""))
            if bs_sp is None or bs_mk is None or base_sp is None or base_mk is None:
                continue
            deltas.append(
                {
                    "dataset_id": k.dataset_id,
                    "sample_id": k.sample_id,
                    "K": k.K,
                    "baseline_method_id": baseline,
                    "delta_spatial_coherence": bs_sp - base_sp,
                    "delta_marker_coherence": bs_mk - base_mk,
                    "bayesspace_note": bs.get("notes", ""),
                    "baseline_note": base.get("notes", ""),
                }
            )

    write_tsv(
        results_dir / "replication" / "domain_quality_deltas_by_sample.tsv",
        deltas,
        fieldnames=[
            "dataset_id",
            "sample_id",
            "K",
            "baseline_method_id",
            "delta_spatial_coherence",
            "delta_marker_coherence",
            "bayesspace_note",
            "baseline_note",
        ],
    )

    # BayesSpace stability summaries per sample/K (for replication reporting)
    stability_rows: list[dict[str, Any]] = []
    for k in bayespace_keys:
        bs = index[k]
        stability_rows.append(
            {
                "dataset_id": k.dataset_id,
                "sample_id": k.sample_id,
                "K": k.K,
                "seed_count": bs.get("seed_count", ""),
                "stability_ari_median": bs.get("stability_ari_median", ""),
                "stability_ari_iqr": bs.get("stability_ari_iqr", ""),
                "notes": bs.get("notes", ""),
            }
        )

    write_tsv(
        results_dir / "replication" / "bayesspace_stability_by_sample.tsv",
        stability_rows,
        fieldnames=[
            "dataset_id",
            "sample_id",
            "K",
            "seed_count",
            "stability_ari_median",
            "stability_ari_iqr",
            "notes",
        ],
    )

    # Simple per-dataset/method/K median summaries (for replication overview)
    grouped: dict[tuple[str, str, str], list[dict[str, str]]] = defaultdict(list)
    for row in method_benchmark:
        ds = row.get("dataset_id", "")
        mid = row.get("method_id", "")
        kval = row.get("K", "")
        if not ds or not mid or not kval:
            continue
        grouped[(ds, mid, kval)].append(row)

    def median(values: list[float]) -> float | None:
        if not values:
            return None
        values = sorted(values)
        m = len(values) // 2
        if len(values) % 2 == 1:
            return values[m]
        return 0.5 * (values[m - 1] + values[m])

    summary_rows: list[dict[str, Any]] = []
    for (ds, mid, kval), rows in sorted(grouped.items()):
        def collect(col: str) -> list[float]:
            out = []
            for r in rows:
                try:
                    out.append(float(r.get(col, "")))
                except Exception:
                    continue
            return out

        summary_rows.append(
            {
                "dataset_id": ds,
                "method_id": mid,
                "K": kval,
                "n_samples": len({r.get("sample_id", "") for r in rows if r.get("sample_id", "")}),
                "spatial_coherence_median_of_samples": median(collect("spatial_coherence_median")) or "",
                "marker_coherence_median_of_samples": median(collect("marker_coherence_median")) or "",
                "stability_ari_median_of_samples": median(collect("stability_ari_median")) or "",
                "wall_time_sec_median_of_samples": median(collect("wall_time_sec_median")) or "",
                "failure_rate_median_of_samples": median(collect("failure_rate")) or "",
            }
        )

    write_tsv(
        results_dir / "replication" / "method_benchmark_dataset_level_summary.tsv",
        summary_rows,
        fieldnames=[
            "dataset_id",
            "method_id",
            "K",
            "n_samples",
            "spatial_coherence_median_of_samples",
            "marker_coherence_median_of_samples",
            "stability_ari_median_of_samples",
            "wall_time_sec_median_of_samples",
            "failure_rate_median_of_samples",
        ],
    )

    return 0


if __name__ == "__main__":
    raise SystemExit(main())

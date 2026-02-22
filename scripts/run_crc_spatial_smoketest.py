#!/usr/bin/env python3
"""Run a local smoke test for CRC spatial-domain benchmarking."""

from __future__ import annotations

import argparse
import csv
import gzip
import itertools
import pathlib
import time
from datetime import datetime, timezone
from typing import TypedDict

import h5py
import numpy as np
import pandas as pd
import psutil
from scipy import sparse
from scipy.io import mmread
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.metrics import adjusted_rand_score
from sklearn.neighbors import NearestNeighbors


def utc_now() -> str:
    return datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")


def _read_lines(path: pathlib.Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return [line.strip() for line in handle]
    return path.read_text(encoding="utf-8").splitlines()


class SampleEntry(TypedDict):
    sample_id: str
    matrix_path: str
    barcodes_path: str
    coords_path: str


def _discover_mtx_samples(dataset_extract_dir: pathlib.Path) -> list[SampleEntry]:
    entries: list[SampleEntry] = []
    for matrix_path in dataset_extract_dir.rglob("matrix.mtx*"):
        matrix_dir = matrix_path.parent
        if not (matrix_dir / "barcodes.tsv.gz").exists() and not (matrix_dir / "barcodes.tsv").exists():
            continue
        sample_root = matrix_dir.parent
        coords_candidates = list(sample_root.rglob("tissue_positions*.csv*"))
        if not coords_candidates:
            continue
        barcodes_path = resolve_file(matrix_dir, ["barcodes.tsv.gz", "barcodes.tsv"])
        entries.append(
            {
                "sample_id": sample_root.name,
                "matrix_path": str(matrix_path),
                "barcodes_path": str(barcodes_path),
                "coords_path": str(coords_candidates[0]),
            }
        )
    return entries


def _discover_flat_mtx_samples(dataset_extract_dir: pathlib.Path) -> list[SampleEntry]:
    entries: list[SampleEntry] = []
    for matrix_path in dataset_extract_dir.glob("*_matrix.mtx*"):
        matrix_name = matrix_path.name
        prefix = matrix_name.replace("_matrix.mtx.gz", "").replace("_matrix.mtx", "")
        barcodes_candidates = [
            dataset_extract_dir / f"{prefix}_barcodes.tsv.gz",
            dataset_extract_dir / f"{prefix}_barcodes.tsv",
        ]
        coords_candidates = [
            dataset_extract_dir / f"{prefix}_tissue_positions.csv.gz",
            dataset_extract_dir / f"{prefix}_tissue_positions.csv",
            dataset_extract_dir / f"{prefix}_tissue_positions_list.csv.gz",
            dataset_extract_dir / f"{prefix}_tissue_positions_list.csv",
        ]
        barcodes_path = None
        coords_path = None
        for candidate in barcodes_candidates:
            if candidate.exists():
                barcodes_path = candidate
                break
        for candidate in coords_candidates:
            if candidate.exists():
                coords_path = candidate
                break
        if barcodes_path is None or coords_path is None:
            continue
        entries.append(
            {
                "sample_id": prefix,
                "matrix_path": str(matrix_path),
                "barcodes_path": str(barcodes_path),
                "coords_path": str(coords_path),
            }
        )
    return entries


def _discover_h5_samples(dataset_extract_dir: pathlib.Path) -> list[SampleEntry]:
    entries: list[SampleEntry] = []
    for h5_path in dataset_extract_dir.rglob("*_filtered_feature_bc_matrix.h5"):
        sample_id = h5_path.name.replace("_filtered_feature_bc_matrix.h5", "")
        coords_candidates = [
            dataset_extract_dir / f"{sample_id}_tissue_positions_list.csv",
            dataset_extract_dir / f"{sample_id}_tissue_positions_list.csv.gz",
            dataset_extract_dir / f"{sample_id}_tissue_positions.csv",
            dataset_extract_dir / f"{sample_id}_tissue_positions.csv.gz",
        ]
        coords_path = None
        for candidate in coords_candidates:
            if candidate.exists():
                coords_path = candidate
                break
        if coords_path is None:
            continue
        entries.append(
            {
                "sample_id": sample_id,
                "matrix_path": str(h5_path),
                "barcodes_path": "",
                "coords_path": str(coords_path),
            }
        )
    return entries


def find_sample_entries(dataset_extract_dir: pathlib.Path) -> list[SampleEntry]:
    entries = _discover_mtx_samples(dataset_extract_dir)
    entries.extend(_discover_flat_mtx_samples(dataset_extract_dir))
    entries.extend(_discover_h5_samples(dataset_extract_dir))
    dedup: dict[str, SampleEntry] = {}
    for entry in entries:
        dedup[entry["sample_id"]] = entry
    return [dedup[key] for key in sorted(dedup)]


def resolve_file(parent: pathlib.Path, names: list[str]) -> pathlib.Path:
    for name in names:
        candidate = parent / name
        if candidate.exists():
            return candidate
    raise FileNotFoundError(f"Missing one of {names} under {parent}")


def load_sample(sample_entry: SampleEntry) -> tuple[str, sparse.csr_matrix, np.ndarray]:
    matrix_path = pathlib.Path(sample_entry["matrix_path"])
    coords_path = pathlib.Path(sample_entry["coords_path"])

    if matrix_path.suffix == ".h5":
        with h5py.File(matrix_path, "r") as handle:
            group = handle["matrix"]
            data = np.array(group["data"])
            indices = np.array(group["indices"])
            indptr = np.array(group["indptr"])
            shape = tuple(np.array(group["shape"]).tolist())
            matrix = sparse.csc_matrix((data, indices, indptr), shape=shape).tocsr()
            barcodes = [
                value.decode("utf-8") if isinstance(value, (bytes, bytearray)) else str(value)
                for value in np.array(group["barcodes"])
            ]
    else:
        barcodes_path = pathlib.Path(sample_entry["barcodes_path"])
        matrix = mmread(matrix_path).tocsr()
        barcodes = _read_lines(barcodes_path)

    if matrix.shape[1] != len(barcodes):
        raise ValueError(
            f"Barcode mismatch for {sample_entry['sample_id']}: matrix spots={matrix.shape[1]}, "
            f"barcodes={len(barcodes)}"
        )

    coords = pd.read_csv(coords_path, header=None, compression="infer")
    if isinstance(coords.iloc[0, 0], str) and coords.iloc[0, 0] == "barcode":
        coords = pd.read_csv(coords_path, compression="infer")
    else:
        coords.columns = [
            "barcode",
            "in_tissue",
            "array_row",
            "array_col",
            "pxl_row_in_fullres",
            "pxl_col_in_fullres",
        ]

    coords = coords[coords["barcode"].isin(barcodes)].copy()
    coords["barcode"] = pd.Categorical(coords["barcode"], categories=barcodes, ordered=True)
    coords = coords.sort_values("barcode")
    if "in_tissue" in coords.columns:
        in_tissue = coords["in_tissue"].to_numpy().astype(int) == 1
    else:
        in_tissue = np.ones(len(coords), dtype=bool)

    # 10x matrix orientation is genes x spots; transpose to spots x genes.
    matrix = matrix.transpose().tocsr()
    matrix = matrix[in_tissue, :]
    coord_array = coords.loc[in_tissue, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()

    sample_id = sample_entry["sample_id"]
    return sample_id, matrix, coord_array


def normalize_and_select(matrix: sparse.csr_matrix, max_genes: int) -> np.ndarray:
    counts_per_spot = np.asarray(matrix.sum(axis=1)).ravel()
    counts_per_spot[counts_per_spot == 0] = 1.0
    scale = 1e4 / counts_per_spot
    normalized = matrix.multiply(scale[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)

    means = np.asarray(normalized.mean(axis=0)).ravel()
    sq_means = np.asarray(normalized.power(2).mean(axis=0)).ravel()
    variances = np.maximum(sq_means - means**2, 0.0)
    top_idx = np.argsort(variances)[::-1][:max_genes]
    dense = normalized[:, top_idx].toarray().astype(np.float32)
    return dense


def build_pcs(x: np.ndarray, n_components: int = 20) -> np.ndarray:
    n_components = max(2, min(n_components, x.shape[0] - 1, x.shape[1] - 1))
    model = PCA(n_components=n_components, random_state=0)
    return model.fit_transform(x)


def spatial_coherence(labels: np.ndarray, coords: np.ndarray, neighbors: int = 6) -> float:
    nn = NearestNeighbors(n_neighbors=min(neighbors + 1, len(coords)), algorithm="auto")
    nn.fit(coords)
    idx = nn.kneighbors(return_distance=False)
    neighbor_idx = idx[:, 1:]
    matches = (labels[:, None] == labels[neighbor_idx]).mean()
    return float(matches)


def spatial_connectivity_graph(coords: np.ndarray, neighbors: int = 6) -> sparse.csr_matrix:
    nn = NearestNeighbors(n_neighbors=min(neighbors + 1, len(coords)), algorithm="auto")
    nn.fit(coords)
    graph = nn.kneighbors_graph(coords, mode="connectivity")
    graph = graph.maximum(graph.T)
    return graph.tocsr()


def marker_separation_score(x: np.ndarray, labels: np.ndarray) -> float:
    scores: list[float] = []
    for cluster in np.unique(labels):
        in_mask = labels == cluster
        out_mask = ~in_mask
        if in_mask.sum() < 3 or out_mask.sum() < 3:
            continue
        diff = x[in_mask].mean(axis=0) - x[out_mask].mean(axis=0)
        top = np.sort(diff)[-20:]
        scores.append(float(np.mean(top)))
    return float(np.median(scores)) if scores else float("nan")


def append_rows(path: pathlib.Path, rows: list[dict[str, object]]) -> None:
    if not rows:
        return
    with path.open("r", encoding="utf-8", newline="") as existing:
        header = existing.readline().rstrip("\n").split("\t")
    with path.open("a", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=header, delimiter="\t", extrasaction="ignore")
        for row in rows:
            writer.writerow({k: row.get(k, "") for k in header})


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-id", default="GSE285505")
    parser.add_argument("--dataset-root", default="data/raw/GSE285505/extracted")
    parser.add_argument("--max-samples", type=int, default=1)
    parser.add_argument("--max-genes", type=int, default=2000)
    parser.add_argument("--k-grid", default="4,6")
    parser.add_argument("--seeds", default="11,23,37")
    parser.add_argument(
        "--methods",
        default="M0_expr_kmeans,M1_spatial_concat_kmeans",
        help="Comma-separated method IDs to run (baselines only).",
    )
    parser.add_argument("--results-bench", default="results/benchmarks")
    parser.add_argument("--results-fig", default="results/figures")
    parser.add_argument("--hardware-id", default="local-macbook-air-8gb-smoketest")
    parser.add_argument("--note", default="smoke-test")
    args = parser.parse_args()

    dataset_root = pathlib.Path(args.dataset_root)
    bench_dir = pathlib.Path(args.results_bench)
    fig_dir = pathlib.Path(args.results_fig)
    k_grid = [int(x) for x in args.k_grid.split(",") if x.strip()]
    seeds = [int(x) for x in args.seeds.split(",") if x.strip()]
    methods_requested = [m.strip() for m in args.methods.split(",") if m.strip()]

    sample_entries = find_sample_entries(dataset_root)[: args.max_samples]
    if not sample_entries:
        raise FileNotFoundError(f"No sample entries found under {dataset_root}")

    process = psutil.Process()
    method_rows: list[dict[str, object]] = []
    stability_rows: list[dict[str, object]] = []
    spatial_rows: list[dict[str, object]] = []
    marker_rows: list[dict[str, object]] = []
    runtime_rows: list[dict[str, object]] = []
    failure_rows: list[dict[str, object]] = []
    control_rows: list[dict[str, object]] = []
    sensitivity_rows: list[dict[str, object]] = []

    fig1_rows: list[dict[str, object]] = []
    fig2_rows: list[dict[str, object]] = []
    fig3_rows: list[dict[str, object]] = []
    fig4_rows: list[dict[str, object]] = []

    for sample_entry in sample_entries:
        sample_id, matrix, coords = load_sample(sample_entry)
        x = normalize_and_select(matrix, max_genes=args.max_genes)
        pcs = build_pcs(x)
        coords_scaled = (coords - coords.mean(axis=0)) / (coords.std(axis=0) + 1e-6)
        spatial_features = np.concatenate([pcs, 0.5 * coords_scaled], axis=1)
        spatial_graph = spatial_connectivity_graph(coords, neighbors=6)

        method_configs: list[tuple[str, str, np.ndarray, list[int]]] = []
        if "M0_expr_kmeans" in methods_requested:
            method_configs.append(("M0_expr_kmeans", "kmeans", pcs, seeds))
        if "M1_spatial_concat_kmeans" in methods_requested:
            method_configs.append(("M1_spatial_concat_kmeans", "kmeans", spatial_features, seeds))
        if "M2_spatial_ward" in methods_requested:
            # Spatially constrained Ward clustering (expression PCs with coordinate connectivity constraint).
            method_configs.append(("M2_spatial_ward", "ward", pcs, [seeds[0]]))

        fig1_rows.append(
            {
                "panel_id": "Fig1A",
                "record_type": "sample_overview",
                "dataset_id": args.dataset_id,
                "sample_id": sample_id,
                "role": "smoketest",
                "modality": "Visium",
                "platform": "10x",
                "organism": "Homo sapiens",
                "tissue": "CRC",
                "n_spots": int(matrix.shape[0]),
                "n_genes": int(matrix.shape[1]),
                "notes": f"{args.note} sample",
            }
        )

        method_control_seed_offsets = {
            "M0_expr_kmeans": 0,
            "M1_spatial_concat_kmeans": 100,
            "M2_spatial_ward": 200,
        }

        for method_id, method_kind, features, method_seeds in method_configs:
            for k in k_grid:
                labels_by_seed: dict[int, np.ndarray] = {}
                run_spatial_scores: list[float] = []
                run_marker_scores: list[float] = []
                run_times: list[float] = []
                run_mem: list[float] = []

                for seed in method_seeds:
                    started = time.perf_counter()
                    started_utc = utc_now()
                    status = "success"
                    error_message = ""
                    labels = None
                    try:
                        if method_kind == "kmeans":
                            model = KMeans(
                                n_clusters=k,
                                random_state=seed,
                                n_init=10,
                                max_iter=300,
                            )
                            labels = model.fit_predict(features)
                        elif method_kind == "ward":
                            model = AgglomerativeClustering(
                                n_clusters=k,
                                linkage="ward",
                                connectivity=spatial_graph,
                            )
                            labels = model.fit_predict(features)
                        else:
                            raise ValueError(f"Unknown method_kind={method_kind}")
                    except Exception as exc:  # pragma: no cover - defensive logging
                        status = "failed"
                        error_message = str(exc)
                    elapsed = time.perf_counter() - started
                    finished_utc = utc_now()
                    peak_rss_mb = process.memory_info().rss / (1024 * 1024)

                    runtime_rows.append(
                        {
                            "dataset_id": args.dataset_id,
                            "sample_id": sample_id,
                            "method_id": method_id,
                            "preprocessing_id": "log1p_hvg",
                            "param_set_id": f"K{k}",
                            "K": k,
                            "seed": seed,
                            "status": status,
                            "wall_time_sec": round(elapsed, 6),
                            "cpu_time_sec": round(elapsed, 6),
                            "peak_rss_mb": round(peak_rss_mb, 3),
                            "disk_tmp_bytes": "",
                            "started_utc": started_utc,
                            "finished_utc": finished_utc,
                            "hardware_id": args.hardware_id,
                            "software_versions": "python-smoketest",
                            "notes": "",
                        }
                    )

                    if labels is None:
                        failure_rows.append(
                            {
                                "dataset_id": args.dataset_id,
                                "sample_id": sample_id,
                                "method_id": method_id,
                                "preprocessing_id": "log1p_hvg",
                                "param_set_id": f"K{k}",
                                "K": k,
                                "seed": seed,
                                "failure_type": "runtime_error",
                                "error_message": error_message,
                                "log_path": "",
                                "stacktrace_path": "",
                                "wall_time_sec": round(elapsed, 6),
                                "peak_rss_mb": round(peak_rss_mb, 3),
                                "finished_utc": finished_utc,
                                "notes": "",
                            }
                        )
                        continue

                    labels_by_seed[seed] = labels
                    sp_score = spatial_coherence(labels, coords)
                    mk_score = marker_separation_score(x, labels)
                    run_spatial_scores.append(sp_score)
                    run_marker_scores.append(mk_score)
                    run_times.append(elapsed)
                    run_mem.append(peak_rss_mb)

                if not labels_by_seed:
                    continue

                ari_values: list[float] = []
                for a, b in itertools.combinations(labels_by_seed, 2):
                    ari = adjusted_rand_score(labels_by_seed[a], labels_by_seed[b])
                    ari_values.append(float(ari))
                    stability_rows.append(
                        {
                            "dataset_id": args.dataset_id,
                            "sample_id": sample_id,
                            "method_id": method_id,
                            "preprocessing_id": "log1p_hvg",
                            "param_set_id": f"K{k}",
                            "K": k,
                            "stability_type": "seed_pair_ari",
                            "replicate_id": f"{a}_vs_{b}",
                            "seed": a,
                            "seed2": b,
                            "metric_name": "ari",
                            "metric_value": round(float(ari), 6),
                            "n_spots": int(matrix.shape[0]),
                            "notes": "",
                        }
                    )

                for seed, labels in labels_by_seed.items():
                    sp_score = spatial_coherence(labels, coords)
                    spatial_rows.append(
                        {
                            "dataset_id": args.dataset_id,
                            "sample_id": sample_id,
                            "method_id": method_id,
                            "preprocessing_id": "log1p_hvg",
                            "param_set_id": f"K{k}",
                            "K": k,
                            "metric_name": "neighbor_agreement",
                            "metric_value": round(sp_score, 6),
                            "bootstrap_ci_lower": "",
                            "bootstrap_ci_upper": "",
                            "n_spots": int(matrix.shape[0]),
                            "notes": "",
                        }
                    )

                reference_seed = method_seeds[0]
                ref_labels = labels_by_seed[reference_seed]
                for cluster in np.unique(ref_labels):
                    in_mask = ref_labels == cluster
                    out_mask = ~in_mask
                    if in_mask.sum() < 3 or out_mask.sum() < 3:
                        continue
                    diff = x[in_mask].mean(axis=0) - x[out_mask].mean(axis=0)
                    marker_rows.append(
                        {
                            "dataset_id": args.dataset_id,
                            "sample_id": sample_id,
                            "method_id": method_id,
                            "preprocessing_id": "log1p_hvg",
                            "param_set_id": f"K{k}",
                            "K": k,
                            "domain_id": int(cluster),
                            "top_marker_n": 20,
                            "marker_reproducibility_within_dataset": round(
                                float(np.mean(np.sort(diff)[-20:])), 6
                            ),
                            "marker_reproducibility_across_datasets": "",
                            "enrichment_coherence_score": round(float(np.max(diff)), 6),
                            "notes": "",
                        }
                    )

                offset = method_control_seed_offsets.get(method_id, 0)
                random_rng = np.random.default_rng(20260211 + k + offset)
                rand_labels = random_rng.integers(0, k, size=ref_labels.shape[0])
                rand_sp = spatial_coherence(rand_labels, coords)
                real_sp = spatial_coherence(ref_labels, coords)
                control_rows.append(
                    {
                        "dataset_id": args.dataset_id,
                        "sample_id": sample_id,
                        "method_id": method_id,
                        "preprocessing_id": "log1p_hvg",
                        "param_set_id": f"K{k}",
                        "K": k,
                        "seed": reference_seed,
                        "control_type": "random_labels",
                        "metric_name": "neighbor_agreement",
                        "metric_value": round(rand_sp, 6),
                        "metric_value_real": round(real_sp, 6),
                        "delta_vs_real": round(rand_sp - real_sp, 6),
                        "notes": "",
                    }
                )

                shuffled_coords = coords.copy()
                random_rng.shuffle(shuffled_coords)
                shuf_sp = spatial_coherence(ref_labels, shuffled_coords)
                control_rows.append(
                    {
                        "dataset_id": args.dataset_id,
                        "sample_id": sample_id,
                        "method_id": method_id,
                        "preprocessing_id": "log1p_hvg",
                        "param_set_id": f"K{k}",
                        "K": k,
                        "seed": reference_seed,
                        "control_type": "shuffled_coordinates",
                        "metric_name": "neighbor_agreement",
                        "metric_value": round(shuf_sp, 6),
                        "metric_value_real": round(real_sp, 6),
                        "delta_vs_real": round(shuf_sp - real_sp, 6),
                        "notes": "",
                    }
                )

                median_ari = float(np.median(ari_values)) if ari_values else float("nan")
                iqr_ari = (
                    float(np.percentile(ari_values, 75) - np.percentile(ari_values, 25))
                    if ari_values
                    else float("nan")
                )
                median_sp = float(np.median(run_spatial_scores)) if run_spatial_scores else float("nan")
                iqr_sp = (
                    float(np.percentile(run_spatial_scores, 75) - np.percentile(run_spatial_scores, 25))
                    if run_spatial_scores
                    else float("nan")
                )
                median_mk = float(np.median(run_marker_scores)) if run_marker_scores else float("nan")
                iqr_mk = (
                    float(np.percentile(run_marker_scores, 75) - np.percentile(run_marker_scores, 25))
                    if run_marker_scores
                    else float("nan")
                )
                failure_rate = 1.0 - (len(labels_by_seed) / len(method_seeds))

                summary_row = {
                    "dataset_id": args.dataset_id,
                    "sample_id": sample_id,
                    "method_id": method_id,
                    "method_family": "baseline",
                    "preprocessing_id": "log1p_hvg",
                    "param_set_id": f"K{k}",
                    "K": k,
                    "seed_count": len(method_seeds),
                    "stability_ari_median": round(median_ari, 6),
                    "stability_ari_iqr": round(iqr_ari, 6),
                    "spatial_coherence_median": round(median_sp, 6),
                    "spatial_coherence_iqr": round(iqr_sp, 6),
                    "marker_coherence_median": round(median_mk, 6),
                    "marker_coherence_iqr": round(iqr_mk, 6),
                    "wall_time_sec_median": round(float(np.median(run_times)), 6),
                    "peak_rss_mb_median": round(float(np.median(run_mem)), 3),
                    "failure_rate": round(failure_rate, 6),
                    "notes": args.note,
                }
                method_rows.append(summary_row)
                fig2_rows.append({"panel_id": "Fig2A", **summary_row})
                fig3_rows.append(
                    {
                        "panel_id": "Fig3A",
                        "dataset_id": args.dataset_id,
                        "sample_id": sample_id,
                        "method_id": method_id,
                        "sensitivity_factor": "K",
                        "factor_value": k,
                        "K": k,
                        "control_type": "none",
                        "metric_name": "spatial_coherence_median",
                        "metric_value": round(median_sp, 6),
                        "metric_ci_lower": "",
                        "metric_ci_upper": "",
                        "notes": args.note,
                    }
                )
                fig4_rows.append(
                    {
                        "panel_id": "Fig4A",
                        "dataset_id": args.dataset_id,
                        "sample_id": sample_id,
                        "method_id": method_id,
                        "K": k,
                        "status": "success",
                        "wall_time_sec": round(float(np.median(run_times)), 6),
                        "peak_rss_mb": round(float(np.median(run_mem)), 3),
                        "is_laptop_feasible": int(np.median(run_mem) < 6000),
                        "is_cloud_required": int(np.median(run_mem) >= 6000),
                        "recommendation": "local-first",
                        "notes": args.note,
                    }
                )
                sensitivity_rows.append(
                    {
                        "dataset_id": args.dataset_id,
                        "sample_id": sample_id,
                        "method_id": method_id,
                        "preprocessing_id": "log1p_hvg",
                        "param_set_id": f"K{k}",
                        "sensitivity_factor": "K",
                        "factor_value": k,
                        "K": k,
                        "metric_name": "spatial_coherence_median",
                        "metric_value": round(median_sp, 6),
                        "metric_ci_lower": "",
                        "metric_ci_upper": "",
                        "notes": args.note,
                    }
                )

    if not failure_rows:
        failure_rows.append(
            {
                "dataset_id": args.dataset_id,
                "sample_id": sample_entries[0]["sample_id"],
                "method_id": "none",
                "preprocessing_id": "log1p_hvg",
                "param_set_id": "NA",
                "K": "",
                "seed": "",
                "failure_type": "none",
                "error_message": "",
                "log_path": "",
                "stacktrace_path": "",
                "wall_time_sec": "",
                "peak_rss_mb": "",
                "finished_utc": utc_now(),
                "notes": f"{args.note} completed without runtime failures",
            }
        )

    append_rows(bench_dir / "method_benchmark.tsv", method_rows)
    append_rows(bench_dir / "stability.tsv", stability_rows)
    append_rows(bench_dir / "spatial_coherence.tsv", spatial_rows)
    append_rows(bench_dir / "marker_coherence.tsv", marker_rows)
    append_rows(bench_dir / "runtime_memory.tsv", runtime_rows)
    append_rows(bench_dir / "failure_log.tsv", failure_rows)
    append_rows(bench_dir / "negative_controls.tsv", control_rows)
    append_rows(bench_dir / "sensitivity_summary.tsv", sensitivity_rows)

    append_rows(fig_dir / "fig1_overview.tsv", fig1_rows)
    append_rows(fig_dir / "fig2_benchmark_summary.tsv", fig2_rows)
    append_rows(fig_dir / "fig3_sensitivity.tsv", fig3_rows)
    append_rows(fig_dir / "fig4_compute_and_guidance.tsv", fig4_rows)

    print(
        f"[done] dataset={args.dataset_id} samples={len(sample_entries)} "
        f"method_rows={len(method_rows)} runtime_rows={len(runtime_rows)}"
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

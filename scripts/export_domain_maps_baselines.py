#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import gzip
import pathlib

import numpy as np
import pandas as pd
from scipy import sparse
from scipy.io import mmread
from sklearn.cluster import AgglomerativeClustering, KMeans
from sklearn.decomposition import PCA
from sklearn.neighbors import NearestNeighbors


def _read_lines(path: pathlib.Path) -> list[str]:
    if path.suffix == ".gz":
        with gzip.open(path, "rt", encoding="utf-8") as handle:
            return [line.strip() for line in handle]
    return path.read_text(encoding="utf-8").splitlines()


def load_flat_sample(dataset_root: pathlib.Path, sample_id: str) -> tuple[sparse.csr_matrix, np.ndarray]:
    matrix_path = dataset_root / f"{sample_id}_matrix.mtx.gz"
    barcodes_path = dataset_root / f"{sample_id}_barcodes.tsv.gz"
    coords_path = dataset_root / f"{sample_id}_tissue_positions_list.csv.gz"
    if not coords_path.exists():
        coords_path = dataset_root / f"{sample_id}_tissue_positions.csv.gz"

    if not matrix_path.exists():
        raise FileNotFoundError(matrix_path)
    if not barcodes_path.exists():
        raise FileNotFoundError(barcodes_path)
    if not coords_path.exists():
        raise FileNotFoundError(coords_path)

    counts = mmread(matrix_path).tocsr()
    barcodes = _read_lines(barcodes_path)
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

    if counts.shape[1] != len(barcodes):
        raise ValueError(f"Barcode mismatch: matrix spots={counts.shape[1]} barcodes={len(barcodes)}")

    coords = coords[coords["barcode"].isin(barcodes)].copy()
    coords["barcode"] = pd.Categorical(coords["barcode"], categories=barcodes, ordered=True)
    coords = coords.sort_values("barcode")
    in_tissue = coords["in_tissue"].to_numpy().astype(int) == 1
    if in_tissue.sum() < 50:
        raise ValueError("Too few in-tissue spots")

    # counts is genes x spots; transpose to spots x genes
    counts = counts.transpose().tocsr()
    counts = counts[in_tissue, :]
    xy = coords.loc[in_tissue, ["pxl_col_in_fullres", "pxl_row_in_fullres"]].to_numpy()
    return counts, xy


def normalize_and_hvg(counts: sparse.csr_matrix, n_hvg: int = 2000) -> np.ndarray:
    counts_per_spot = np.asarray(counts.sum(axis=1)).ravel()
    counts_per_spot[counts_per_spot == 0] = 1.0
    scale = 1e4 / counts_per_spot
    normalized = counts.multiply(scale[:, None]).tocsr()
    normalized.data = np.log1p(normalized.data)

    means = np.asarray(normalized.mean(axis=0)).ravel()
    sq_means = np.asarray(normalized.power(2).mean(axis=0)).ravel()
    variances = np.maximum(sq_means - means**2, 0.0)
    top_idx = np.argsort(variances)[::-1][: min(n_hvg, normalized.shape[1])]
    dense = normalized[:, top_idx].toarray().astype(np.float32)
    return dense


def build_pcs(x: np.ndarray, n_components: int = 20) -> np.ndarray:
    n_components = max(2, min(n_components, x.shape[0] - 1, x.shape[1] - 1))
    model = PCA(n_components=n_components, random_state=0)
    return model.fit_transform(x)


def main() -> int:
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset-id", required=True)
    parser.add_argument("--dataset-root", required=True)
    parser.add_argument("--sample-id", required=True)
    parser.add_argument("--k", type=int, default=6)
    parser.add_argument("--seed", type=int, default=11)
    parser.add_argument("--output-tsv", required=True)
    args = parser.parse_args()

    dataset_root = pathlib.Path(args.dataset_root)
    counts, xy = load_flat_sample(dataset_root, args.sample_id)
    x = normalize_and_hvg(counts, n_hvg=2000)
    pcs = build_pcs(x, n_components=20)
    coords_scaled = (xy - xy.mean(axis=0)) / (xy.std(axis=0) + 1e-6)
    spatial_features = np.concatenate([pcs, 0.5 * coords_scaled], axis=1)

    nn = NearestNeighbors(n_neighbors=min(7, len(xy)), algorithm="auto")
    nn.fit(xy)
    spatial_graph = nn.kneighbors_graph(xy, mode="connectivity")
    spatial_graph = spatial_graph.maximum(spatial_graph.T)

    rows: list[dict[str, object]] = []
    for method_id, method_kind, features in [
        ("M0_expr_kmeans", "kmeans", pcs),
        ("M1_spatial_concat_kmeans", "kmeans", spatial_features),
        ("M2_spatial_ward", "ward", pcs),
    ]:
        if method_kind == "kmeans":
            model = KMeans(n_clusters=args.k, random_state=args.seed, n_init=10, max_iter=300)
            labels = model.fit_predict(features)
        elif method_kind == "ward":
            model = AgglomerativeClustering(n_clusters=args.k, linkage="ward", connectivity=spatial_graph)
            labels = model.fit_predict(features)
        else:
            raise ValueError(f"Unknown method_kind={method_kind}")

        for (xv, yv), lab in zip(xy, labels, strict=True):
            rows.append(
                {
                    "dataset_id": args.dataset_id,
                    "sample_id": args.sample_id,
                    "method_id": method_id,
                    "K": args.k,
                    "x": float(xv),
                    "y": float(yv),
                    "domain_label": int(lab) + 1,
                    "notes": "baseline-domain-map",
                }
            )

    out_path = pathlib.Path(args.output_tsv)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["dataset_id", "sample_id", "method_id", "K", "x", "y", "domain_label", "notes"]
    with out_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        for row in rows:
            writer.writerow(row)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())


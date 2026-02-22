#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  library(Matrix)
  library(SingleCellExperiment)
  library(BayesSpace)
})

parse_args <- function(input) {
  result <- list(
    dataset_id = "",
    dataset_root = "",
    sample_id = "",
    k = "6",
    seed = "11",
    output_tsv = "",
    note = "bayesspace-domain-map"
  )
  i <- 1L
  while (i <= length(input)) {
    key <- input[[i]]
    if (startsWith(key, "--") && i < length(input)) {
      value <- input[[i + 1L]]
      name <- gsub("^--", "", key)
      name <- gsub("-", "_", name)
      result[[name]] <- value
      i <- i + 2L
    } else {
      i <- i + 1L
    }
  }
  result
}

args <- parse_args(commandArgs(trailingOnly = TRUE))
if (!nzchar(args$dataset_id) || !nzchar(args$dataset_root) || !nzchar(args$sample_id) || !nzchar(args$output_tsv)) {
  stop("Missing required args: --dataset-id --dataset-root --sample-id --output-tsv")
}

dataset_root <- args$dataset_root
dataset_id <- args$dataset_id
sample_id <- args$sample_id
q <- as.integer(args$k)
seed <- as.integer(args$seed)

matrix_file <- file.path(dataset_root, paste0(sample_id, "_matrix.mtx.gz"))
barcodes_file <- file.path(dataset_root, paste0(sample_id, "_barcodes.tsv.gz"))
features_file <- file.path(dataset_root, paste0(sample_id, "_features.tsv.gz"))
coords_file <- file.path(dataset_root, paste0(sample_id, "_tissue_positions_list.csv.gz"))
if (!file.exists(coords_file)) {
  coords_file <- file.path(dataset_root, paste0(sample_id, "_tissue_positions.csv.gz"))
}

if (!file.exists(matrix_file) || !file.exists(barcodes_file) || !file.exists(features_file) || !file.exists(coords_file)) {
  stop("Missing flat-matrix Visium files for requested sample_id under dataset_root")
}

counts <- readMM(gzfile(matrix_file))
barcodes <- read.delim(gzfile(barcodes_file), header = FALSE, stringsAsFactors = FALSE)
features <- read.delim(gzfile(features_file), header = FALSE, stringsAsFactors = FALSE)
coords <- read.csv(gzfile(coords_file), stringsAsFactors = FALSE, header = FALSE)

if (is.character(coords[1, 1]) && coords[1, 1] == "barcode") {
  coords <- read.csv(gzfile(coords_file), stringsAsFactors = FALSE)
} else {
  colnames(coords) <- c("barcode", "in_tissue", "array_row", "array_col", "pxl_row_in_fullres", "pxl_col_in_fullres")
}

rownames(counts) <- make.unique(features[[2]])
colnames(counts) <- barcodes[[1]]

coords <- coords[coords$barcode %in% colnames(counts), , drop = FALSE]
coords <- coords[match(colnames(counts), coords$barcode), , drop = FALSE]
in_tissue <- coords$in_tissue == 1
if (sum(in_tissue) < 50) {
  stop("Too few in-tissue spots for BayesSpace domain-map export")
}

counts <- counts[, in_tissue, drop = FALSE]
coords <- coords[in_tissue, , drop = FALSE]

sce <- SingleCellExperiment(assays = list(counts = counts))
colData(sce)$array_row <- coords$array_row
colData(sce)$array_col <- coords$array_col
colData(sce)$row <- coords$array_row
colData(sce)$col <- coords$array_col
colData(sce)$imagerow <- coords$pxl_row_in_fullres
colData(sce)$imagecol <- coords$pxl_col_in_fullres

set.seed(seed)
sce_pre <- spatialPreprocess(
  sce,
  platform = "Visium",
  n.HVGs = min(2000, nrow(sce)),
  n.PCs = min(15, ncol(sce) - 1L),
  log.normalize = TRUE
)

cluster_args <- list(
  sce = sce_pre,
  q = q,
  platform = "Visium",
  d = min(15, {
    dim_names <- reducedDimNames(sce_pre)
    if ("PCA" %in% dim_names) {
      ncol(reducedDim(sce_pre, "PCA"))
    } else if (length(dim_names) > 0) {
      ncol(reducedDim(sce_pre, dim_names[[1]]))
    } else {
      15L
    }
  }),
  init.method = "mclust",
  model = "t",
  gamma = 3,
  nrep = 200,
  save.chain = FALSE
)

if ("burn.in" %in% names(formals(BayesSpace::spatialCluster))) {
  cluster_args[["burn.in"]] <- min(100L, as.integer(cluster_args$nrep) - 1L)
}
if ("verbose" %in% names(formals(BayesSpace::spatialCluster))) {
  cluster_args$verbose <- FALSE
}

sce_q <- do.call(BayesSpace::spatialCluster, cluster_args)
labels <- as.integer(colData(sce_q)$spatial.cluster)
imagecol <- as.numeric(colData(sce_q)$imagecol)
imagerow <- as.numeric(colData(sce_q)$imagerow)

out <- data.frame(
  dataset_id = dataset_id,
  sample_id = sample_id,
  method_id = "BayesSpace",
  K = q,
  x = imagecol,
  y = imagerow,
  domain_label = labels,
  notes = args$note,
  stringsAsFactors = FALSE
)

write.table(out, file = args$output_tsv, sep = "\t", quote = FALSE, row.names = FALSE)
cat(sprintf("Wrote %d domain-map rows to %s\n", nrow(out), args$output_tsv))

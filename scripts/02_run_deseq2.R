args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(DESeq2)
  library(ggplot2)
})

config <- yaml::read_yaml(file.path(project_dir, "config", "analysis_config.yml"))
metadata <- readr::read_csv(file.path(project_dir, config$inputs$metadata_file), show_col_types = FALSE)
count_path <- file.path(project_dir, config$inputs$count_matrix)

if (!file.exists(count_path)) {
  stop("Raw count matrix not found at: ", count_path)
}

count_tbl <- readr::read_tsv(count_path, show_col_types = FALSE)
sample_id_col <- config$inputs$sample_id_column
batch_col <- config$design$batch_column
condition_col <- config$design$condition_column
annotation_cols <- unique(unlist(config$inputs$annotation_columns))
annotation_cols <- annotation_cols[annotation_cols %in% colnames(count_tbl)]
sample_cols <- setdiff(colnames(count_tbl), annotation_cols)

missing_samples <- setdiff(metadata[[sample_id_col]], sample_cols)
if (length(missing_samples) > 0) {
  stop("Count matrix is missing samples from metadata: ", paste(missing_samples, collapse = ", "))
}

metadata <- metadata[match(sample_cols, metadata[[sample_id_col]]), , drop = FALSE]
rownames(metadata) <- metadata[[sample_id_col]]

count_mat <- as.matrix(count_tbl[, sample_cols, drop = FALSE])
storage.mode(count_mat) <- "integer"
rownames(count_mat) <- count_tbl[[config$inputs$gene_id_column]]

if (!batch_col %in% colnames(metadata) || length(unique(stats::na.omit(metadata[[batch_col]]))) <= 1) {
  design_formula <- as.formula("~ condition")
} else {
  metadata[[batch_col]] <- factor(metadata[[batch_col]])
  design_formula <- as.formula(config$design$design_formula)
}

metadata[[condition_col]] <- factor(metadata[[condition_col]])

dds <- DESeqDataSetFromMatrix(
  countData = count_mat,
  colData = metadata,
  design = design_formula
)
dds <- dds[rowSums(counts(dds)) > 10, ]
dds <- DESeq(dds)

norm_counts <- counts(dds, normalized = TRUE)
vst_obj <- vst(dds, blind = FALSE)
vst_mat <- assay(vst_obj)

readr::write_csv(
  cbind(data.frame(gene_id = rownames(norm_counts), stringsAsFactors = FALSE), as.data.frame(norm_counts, check.names = FALSE)),
  file.path(project_dir, "results", "deseq2", "normalized_counts.csv")
)
readr::write_csv(
  cbind(data.frame(gene_id = rownames(vst_mat), stringsAsFactors = FALSE), as.data.frame(vst_mat, check.names = FALSE)),
  file.path(project_dir, "results", "deseq2", "vst_counts.csv")
)
writeLines(resultsNames(dds), con = file.path(project_dir, "results", "deseq2", "deseq2_resultsNames.txt"))

pca <- prcomp(t(vst_mat), center = TRUE, scale. = TRUE)
var_explained <- (pca$sdev^2) / sum(pca$sdev^2)
pca_df <- data.frame(
  sample_id = rownames(pca$x),
  PC1 = pca$x[, 1],
  PC2 = pca$x[, 2],
  metadata[rownames(pca$x), , drop = FALSE],
  stringsAsFactors = FALSE
)
readr::write_csv(pca_df, file.path(project_dir, "results", "deseq2", "pca_coordinates.csv"))

p <- ggplot(pca_df, aes(x = PC1, y = PC2, color = .data[[condition_col]])) +
  geom_point(size = 3) +
  labs(
    title = "RNA-seq PCA (VST counts)",
    x = sprintf("PC1 (%.1f%%)", 100 * var_explained[1]),
    y = sprintf("PC2 (%.1f%%)", 100 * var_explained[2]),
    color = condition_col
  ) +
  theme_bw(base_size = 12)
ggsave(file.path(project_dir, "results", "figures", "01_pca_vst.png"), p, width = 8, height = 6, dpi = 300)

for (contrast_def in config$contrasts) {
  res <- results(
    dds,
    contrast = c(contrast_def$factor, contrast_def$numerator, contrast_def$denominator),
    alpha = config$thresholds$padj
  )
  res_df <- as.data.frame(res)
  res_df$gene_id <- rownames(res_df)
  out_name <- paste0(contrast_def$name, "_deseq2_results.csv")
  readr::write_csv(res_df, file.path(project_dir, "results", "deseq2", out_name))
}

message("DESeq2 workflow complete.")

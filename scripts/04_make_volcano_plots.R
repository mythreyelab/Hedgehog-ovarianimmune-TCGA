args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(ggplot2)
})

config <- yaml::read_yaml(file.path(project_dir, "config", "analysis_config.yml"))
padj_cutoff <- config$thresholds$padj
lfc_cutoff <- config$thresholds$abs_log2fc

for (contrast_def in config$contrasts) {
  res_path <- file.path(project_dir, "results", "deseq2", paste0(contrast_def$name, "_deseq2_results.csv"))
  if (!file.exists(res_path)) {
    next
  }

  res_df <- readr::read_csv(res_path, show_col_types = FALSE)
  res_df$neglog10_padj <- -log10(res_df$padj)
  res_df$significance <- "Not significant"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < padj_cutoff & !is.na(res_df$log2FoldChange) & res_df$log2FoldChange >= lfc_cutoff] <- "Up"
  res_df$significance[!is.na(res_df$padj) & res_df$padj < padj_cutoff & !is.na(res_df$log2FoldChange) & res_df$log2FoldChange <= -lfc_cutoff] <- "Down"

  p <- ggplot(res_df, aes(x = log2FoldChange, y = neglog10_padj, color = significance)) +
    geom_point(alpha = 0.6, size = 1.2) +
    geom_vline(xintercept = c(-lfc_cutoff, lfc_cutoff), linetype = 2, color = "gray50") +
    geom_hline(xintercept = NULL, yintercept = -log10(padj_cutoff), linetype = 2, color = "gray50") +
    scale_color_manual(values = c("Up" = "#B2182B", "Down" = "#2166AC", "Not significant" = "gray75")) +
    labs(
      title = paste("Volcano plot:", contrast_def$name),
      x = "log2 fold change",
      y = "-log10 adjusted p-value",
      color = NULL
    ) +
    theme_bw(base_size = 12)

  ggsave(
    file.path(project_dir, "results", "figures", paste0("02_volcano_", contrast_def$name, ".png")),
    p,
    width = 7,
    height = 6,
    dpi = 300
  )
}

message("Volcano plot workflow complete.")

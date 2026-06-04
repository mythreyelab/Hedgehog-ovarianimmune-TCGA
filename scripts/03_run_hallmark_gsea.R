args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
  library(dplyr)
  library(msigdbr)
  library(fgsea)
  library(ggplot2)
})

config <- yaml::read_yaml(file.path(project_dir, "config", "analysis_config.yml"))
hallmark_tbl <- msigdbr(species = config$gsea$species, category = config$gsea$category)
hallmark_sets <- split(hallmark_tbl$gene_symbol, hallmark_tbl$gs_name)

for (contrast_def in config$contrasts) {
  res_path <- file.path(project_dir, "results", "deseq2", paste0(contrast_def$name, "_deseq2_results.csv"))
  if (!file.exists(res_path)) {
    next
  }

  res_df <- readr::read_csv(res_path, show_col_types = FALSE) %>%
    filter(!is.na(stat))

  ranks <- res_df$stat
  names(ranks) <- res_df$gene_id
  ranks <- sort(ranks, decreasing = TRUE)

  gsea_res <- fgsea(
    pathways = hallmark_sets,
    stats = ranks,
    minSize = config$gsea$min_size,
    maxSize = config$gsea$max_size,
    nperm = config$gsea$nperm
  ) %>%
    arrange(padj, desc(abs(NES)))

  out_csv <- file.path(project_dir, "results", "gsea", paste0(contrast_def$name, "_hallmark_fgsea.csv"))
  readr::write_csv(gsea_res, out_csv)

  plot_df <- gsea_res %>%
    filter(!is.na(padj)) %>%
    slice_head(n = 10) %>%
    mutate(gs_name = factor(gs_name, levels = rev(gs_name)))

  if (nrow(plot_df) > 0) {
    p <- ggplot(plot_df, aes(x = NES, y = gs_name, color = padj, size = size)) +
      geom_point() +
      labs(
        title = paste("Hallmark GSEA:", contrast_def$name),
        x = "Normalized enrichment score",
        y = NULL,
        color = "FDR",
        size = "Geneset size"
      ) +
      theme_bw(base_size = 12)
    ggsave(
      file.path(project_dir, "results", "figures", paste0("03_hallmark_gsea_", contrast_def$name, ".png")),
      p,
      width = 8,
      height = 6,
      dpi = 300
    )
  }
}

message("Hallmark GSEA workflow complete.")

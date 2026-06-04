args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
})

config <- yaml::read_yaml(file.path(project_dir, "config", "analysis_config.yml"))

manifest <- data.frame(
  figure_id = c("FIG_PCA", "FIG_BARCODE", "FIG_METHODS"),
  output_file = c(
    file.path("results", "figures", "01_pca_vst.png"),
    file.path("results", "manifests", "figure_traceability.csv"),
    "RNAseq_methods_paper_ready.md"
  ),
  generating_script = c(
    file.path("scripts", "02_run_deseq2.R"),
    file.path("scripts", "05_build_figure_traceability.R"),
    "manual_methods_file"
  ),
  primary_inputs = c(
    paste(config$inputs$count_matrix, config$inputs$metadata_file, sep = "; "),
    paste(file.path("scripts", "02_run_deseq2.R"), file.path("scripts", "03_run_hallmark_gsea.R"), file.path("scripts", "04_make_volcano_plots.R"), sep = "; "),
    paste(file.path("config", "analysis_config.yml"), file.path("data", "metadata", "sample_metadata_confirmed_local.csv"), sep = "; ")
  ),
  status = c("planned", "generated", "generated"),
  stringsAsFactors = FALSE
)

for (contrast_def in config$contrasts) {
  manifest <- rbind(
    manifest,
    data.frame(
      figure_id = paste0("VOLCANO_", contrast_def$name),
      output_file = file.path("results", "figures", paste0("02_volcano_", contrast_def$name, ".png")),
      generating_script = file.path("scripts", "04_make_volcano_plots.R"),
      primary_inputs = file.path("results", "deseq2", paste0(contrast_def$name, "_deseq2_results.csv")),
      status = ifelse(file.exists(file.path(project_dir, "results", "figures", paste0("02_volcano_", contrast_def$name, ".png"))), "generated", "planned"),
      stringsAsFactors = FALSE
    ),
    data.frame(
      figure_id = paste0("GSEA_", contrast_def$name),
      output_file = file.path("results", "figures", paste0("03_hallmark_gsea_", contrast_def$name, ".png")),
      generating_script = file.path("scripts", "03_run_hallmark_gsea.R"),
      primary_inputs = file.path("results", "deseq2", paste0(contrast_def$name, "_deseq2_results.csv")),
      status = ifelse(file.exists(file.path(project_dir, "results", "figures", paste0("03_hallmark_gsea_", contrast_def$name, ".png"))), "generated", "planned"),
      stringsAsFactors = FALSE
    )
  )
}

readr::write_csv(manifest, file.path(project_dir, "results", "manifests", "figure_traceability.csv"))
message("Figure traceability manifest written.")

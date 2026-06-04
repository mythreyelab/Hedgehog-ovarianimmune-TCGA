args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
script_dir <- dirname(script_path)
project_dir <- normalizePath(file.path(script_dir, ".."))

steps <- c(
  "01_prepare_inputs.R",
  "02_run_deseq2.R",
  "03_run_hallmark_gsea.R",
  "04_make_volcano_plots.R",
  "05_build_figure_traceability.R"
)

count_matrix <- file.path(project_dir, "data", "raw_counts", "featurecounts_gene_counts.tsv")

for (step in steps) {
  if (!file.exists(count_matrix) && step %in% c("02_run_deseq2.R", "03_run_hallmark_gsea.R", "04_make_volcano_plots.R")) {
    message("Skipping ", step, " because the raw count matrix is not present.")
    next
  }

  step_path <- file.path(script_dir, step)
  status <- system2(file.path(R.home("bin"), "Rscript"), args = shQuote(step_path))
  if (!identical(status, 0L)) {
    stop("Workflow failed at step: ", step)
  }
}

message("Workflow completed.")

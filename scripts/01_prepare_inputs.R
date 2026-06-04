args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

suppressPackageStartupMessages({
  library(yaml)
  library(readr)
})

config <- yaml::read_yaml(file.path(project_dir, "config", "analysis_config.yml"))
metadata_path <- file.path(project_dir, config$inputs$metadata_file)
count_path <- file.path(project_dir, config$inputs$count_matrix)
audit_path <- file.path(project_dir, "results", "manifests", "input_audit.csv")

metadata <- readr::read_csv(metadata_path, show_col_types = FALSE)

required_cols <- c(
  config$inputs$sample_id_column,
  config$design$condition_column
)
missing_cols <- setdiff(required_cols, colnames(metadata))
if (length(missing_cols) > 0) {
  stop("Metadata is missing required columns: ", paste(missing_cols, collapse = ", "))
}

audit <- data.frame(
  item = c(
    "metadata_file_present",
    "count_matrix_present",
    "n_samples_in_metadata",
    "unique_conditions",
    "alignment_version_verified",
    "featurecounts_version_verified",
    "default_design_formula",
    "local_confirmed_factorial_design"
  ),
  value = c(
    file.exists(metadata_path),
    file.exists(count_path),
    nrow(metadata),
    length(unique(metadata[[config$design$condition_column]])),
    config$methods$alignment_version != "VERIFY_VERSION",
    config$methods$counting_version != "VERIFY_VERSION",
    config$design$design_formula,
    config$design$local_confirmed_factorial_design
  ),
  stringsAsFactors = FALSE
)

readr::write_csv(metadata, file.path(project_dir, "results", "manifests", "sample_metadata_validated.csv"))
readr::write_csv(audit, audit_path)

message("Prepared inputs and wrote audit to: ", audit_path)

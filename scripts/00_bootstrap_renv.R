args <- commandArgs(trailingOnly = FALSE)
file_arg <- grep("^--file=", args, value = TRUE)
script_path <- normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)))
project_dir <- normalizePath(file.path(dirname(script_path), ".."))

if (!requireNamespace("renv", quietly = TRUE)) {
  install.packages("renv", repos = "https://cloud.r-project.org")
}

required_packages <- c(
  "yaml",
  "readr",
  "dplyr",
  "tibble",
  "ggplot2",
  "DESeq2",
  "apeglm",
  "ashr",
  "msigdbr",
  "fgsea",
  "EnhancedVolcano"
)

if (!file.exists(file.path(project_dir, "renv.lock"))) {
  renv::init(project = project_dir, bare = TRUE)
}

renv::settings$snapshot.type("implicit", project = project_dir)
renv::install(required_packages, project = project_dir)
renv::snapshot(project = project_dir, prompt = FALSE)

message("renv bootstrap complete for: ", project_dir)

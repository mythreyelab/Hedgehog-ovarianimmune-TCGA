file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(file_arg) == 1) {
  normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("scripts/08_make_ratio_method_figures.R", winslash = "/", mustWork = FALSE)
}
project_dir <- normalizePath(file.path(dirname(script_path), ".."), winslash = "/", mustWork = FALSE)

args <- commandArgs(trailingOnly = TRUE)

input_file <- if (length(args) >= 1) {
  args[[1]]
} else {
  file.path(project_dir, "data", "tcga_expression", "data_mrna_seq_v2_rsem.txt")
}

output_base_dir <- if (length(args) >= 2) {
  args[[2]]
} else {
  file.path(project_dir, "results")
}

figure_dir <- file.path(output_base_dir, "figures")
table_dir <- file.path(output_base_dir, "tables")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

activation_genes <- c("GLI1", "GLI2", "PTCH1", "PTCH2", "HHIP")
ligand_genes <- c("IHH", "SHH")
direct_markers <- c("ARG1", "MRC1", "CD80", "FOXP3", "CD8A", "CD274", "PDCD1", "CCL22", "CCL17", "CCL5", "CXCL10")
suppressive_markers <- c("ARG1", "MRC1", "CD80", "FOXP3", "CD274", "PDCD1", "CCL22", "CCL17")
effector_markers <- c("CD8A", "CCL5", "CXCL10")
ratio_metrics <- c("arg1_cxcl10_ratio", "mrc1_cd8a_ratio")
analysis_metrics <- c(direct_markers, ratio_metrics)

expr_raw <- read.delim(
  input_file,
  sep = "\t",
  header = TRUE,
  check.names = FALSE,
  stringsAsFactors = FALSE,
  na.strings = c("NA", "", "NaN")
)

expr_raw <- expr_raw[!is.na(expr_raw$Hugo_Symbol) & expr_raw$Hugo_Symbol != "", ]
expr_raw <- expr_raw[!duplicated(expr_raw$Hugo_Symbol), ]

sample_ids <- setdiff(colnames(expr_raw), c("Hugo_Symbol", "Entrez_Gene_Id"))
needed_genes <- unique(c(activation_genes, ligand_genes, direct_markers))
missing_genes <- setdiff(needed_genes, expr_raw$Hugo_Symbol)

if (length(missing_genes) > 0) {
  stop(sprintf("Missing genes in expression matrix: %s", paste(missing_genes, collapse = ", ")))
}

subset_raw <- expr_raw[match(needed_genes, expr_raw$Hugo_Symbol), c("Hugo_Symbol", sample_ids)]
expr_matrix <- as.matrix(subset_raw[, sample_ids, drop = FALSE])
rownames(expr_matrix) <- subset_raw$Hugo_Symbol
storage.mode(expr_matrix) <- "numeric"

sample_df <- data.frame(SAMPLE_ID = sample_ids, stringsAsFactors = FALSE)
for (gene in needed_genes) {
  sample_df[[gene]] <- as.numeric(expr_matrix[gene, ])
}

floor_raw <- function(x, floor_value = 1) {
  x <- as.numeric(x)
  x[x < floor_value] <- floor_value
  x
}

floor_and_log2p1 <- function(x, floor_value = 1) {
  log2(floor_raw(x, floor_value = floor_value) + 1)
}

safe_zscore <- function(x) {
  s <- sd(x, na.rm = TRUE)
  if (is.na(s) || s == 0) {
    return(rep(NA_real_, length(x)))
  }
  (x - mean(x, na.rm = TRUE)) / s
}

compute_spearman <- function(x, y) {
  keep <- complete.cases(x, y)
  n <- sum(keep)
  if (n < 3) {
    return(c(rho = NA_real_, p_value = NA_real_, n = n))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p_value = test$p.value, n = n)
}

transformed <- sample_df
for (gene in needed_genes) {
  transformed[[paste0(gene, "_floored")]] <- floor_raw(sample_df[[gene]])
  transformed[[paste0(gene, "_log")]] <- floor_and_log2p1(sample_df[[gene]])
}

transformed$hh_activation_score <- rowMeans(
  sapply(activation_genes, function(gene) safe_zscore(transformed[[paste0(gene, "_log")]])),
  na.rm = FALSE
)
transformed$hh_ligand_score <- rowMeans(
  sapply(ligand_genes, function(gene) safe_zscore(transformed[[paste0(gene, "_log")]])),
  na.rm = FALSE
)

transformed$arg1_cxcl10_ratio_old <- log2((transformed$ARG1_log + 1) / (transformed$CXCL10_log + 1))
transformed$mrc1_cd8a_ratio_old <- log2((transformed$MRC1_log + 1) / (transformed$CD8A_log + 1))
transformed$arg1_cxcl10_ratio_corrected <- log2((transformed$ARG1_floored + 1) / (transformed$CXCL10_floored + 1))
transformed$mrc1_cd8a_ratio_corrected <- log2((transformed$MRC1_floored + 1) / (transformed$CD8A_floored + 1))

metric_category <- function(metric) {
  if (metric %in% suppressive_markers || metric %in% ratio_metrics) {
    return("suppressive")
  }
  "effector"
}

run_analysis <- function(ratio_suffix) {
  out <- data.frame(
    ratio_method = character(),
    score_type = character(),
    metric = character(),
    category = character(),
    rho = numeric(),
    p_value = numeric(),
    fdr = numeric(),
    n = integer(),
    stringsAsFactors = FALSE
  )

  for (score_name in c("hh_activation_score", "hh_ligand_score")) {
    for (metric in direct_markers) {
      stats <- compute_spearman(transformed[[score_name]], transformed[[paste0(metric, "_log")]])
      out <- rbind(
        out,
        data.frame(
          ratio_method = ratio_suffix,
          score_type = score_name,
          metric = metric,
          category = metric_category(metric),
          rho = as.numeric(stats[["rho"]]),
          p_value = as.numeric(stats[["p_value"]]),
          fdr = NA_real_,
          n = as.integer(stats[["n"]]),
          stringsAsFactors = FALSE
        )
      )
    }

    for (metric in ratio_metrics) {
      stats <- compute_spearman(transformed[[score_name]], transformed[[paste0(metric, "_", ratio_suffix)]])
      out <- rbind(
        out,
        data.frame(
          ratio_method = ratio_suffix,
          score_type = score_name,
          metric = metric,
          category = metric_category(metric),
          rho = as.numeric(stats[["rho"]]),
          p_value = as.numeric(stats[["p_value"]]),
          fdr = NA_real_,
          n = as.integer(stats[["n"]]),
          stringsAsFactors = FALSE
        )
      )
    }

    keep <- out$ratio_method == ratio_suffix & out$score_type == score_name
    out$fdr[keep] <- p.adjust(out$p_value[keep], method = "BH")
  }

  out
}

results_old <- run_analysis("old")
results_corrected <- run_analysis("corrected")
all_results <- rbind(results_old, results_corrected)
all_results <- all_results[order(all_results$ratio_method, all_results$score_type, all_results$fdr, -abs(all_results$rho)), ]

comparison <- merge(
  results_old[, c("score_type", "metric", "category", "rho", "fdr")],
  results_corrected[, c("score_type", "metric", "rho", "fdr")],
  by = c("score_type", "metric"),
  suffixes = c("_old", "_corrected")
)
comparison$abs_rho_delta <- abs(comparison$rho_corrected) - abs(comparison$rho_old)
comparison$direction_changed <- sign(comparison$rho_old) != sign(comparison$rho_corrected)
comparison <- comparison[order(comparison$score_type, -abs(comparison$abs_rho_delta)), ]

write.csv(
  transformed[, c(
    "SAMPLE_ID",
    "hh_activation_score",
    "hh_ligand_score",
    paste0(direct_markers, "_log"),
    "arg1_cxcl10_ratio_old",
    "mrc1_cd8a_ratio_old",
    "arg1_cxcl10_ratio_corrected",
    "mrc1_cd8a_ratio_corrected"
  )],
  file = file.path(table_dir, "sample_level_scores_ratio_method_comparison.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  all_results,
  file = file.path(table_dir, "hh_scores_vs_immune_correlations_ratio_method_comparison.csv"),
  row.names = FALSE,
  quote = TRUE
)

write.csv(
  comparison,
  file = file.path(table_dir, "ratio_method_outcome_comparison.csv"),
  row.names = FALSE,
  quote = TRUE
)

summary_lines <- c(
  "TCGA All floor-aware score analysis with ratio-method comparison",
  "",
  paste("Input file:", normalizePath(input_file, winslash = "/", mustWork = FALSE)),
  paste("Samples:", nrow(transformed)),
  paste("Activation genes:", paste(activation_genes, collapse = ", ")),
  paste("Ligand genes:", paste(ligand_genes, collapse = ", ")),
  paste("Direct markers:", paste(direct_markers, collapse = ", ")),
  "",
  "Transformation for direct genes: values < 1 set to 1, then log2(x + 1).",
  "Old ratio method: log2((log-transformed numerator + 1) / (log-transformed denominator + 1))",
  "Corrected ratio method: log2((floored raw numerator + 1) / (floored raw denominator + 1))",
  ""
)

for (ratio_method in c("old", "corrected")) {
  summary_lines <- c(summary_lines, paste("Top activation-score correlations using", ratio_method, "ratios:"))
  top_activation <- head(all_results[all_results$ratio_method == ratio_method & all_results$score_type == "hh_activation_score", ], 8)
  summary_lines <- c(summary_lines, apply(top_activation, 1, function(row) {
    sprintf(
      "%s (%s): rho=%.3f, FDR=%.3g, n=%s",
      row[["metric"]],
      row[["category"]],
      as.numeric(row[["rho"]]),
      as.numeric(row[["fdr"]]),
      row[["n"]]
    )
  }))
  summary_lines <- c(summary_lines, "", paste("Top ligand-score correlations using", ratio_method, "ratios:"))
  top_ligand <- head(all_results[all_results$ratio_method == ratio_method & all_results$score_type == "hh_ligand_score", ], 8)
  summary_lines <- c(summary_lines, apply(top_ligand, 1, function(row) {
    sprintf(
      "%s (%s): rho=%.3f, FDR=%.3g, n=%s",
      row[["metric"]],
      row[["category"]],
      as.numeric(row[["rho"]]),
      as.numeric(row[["fdr"]]),
      row[["n"]]
    )
  }), "")
}

summary_lines <- c(summary_lines, "Metrics most affected by changing the ratio derivation:")
summary_lines <- c(summary_lines, apply(head(comparison, 12), 1, function(row) {
  sprintf(
    "%s / %s: old rho=%.3f (FDR=%.3g), corrected rho=%.3f (FDR=%.3g), abs_delta=%.3f",
    row[["score_type"]],
    row[["metric"]],
    as.numeric(row[["rho_old"]]),
    as.numeric(row[["fdr_old"]]),
    as.numeric(row[["rho_corrected"]]),
    as.numeric(row[["fdr_corrected"]]),
    as.numeric(row[["abs_rho_delta"]])
  )
}))

writeLines(summary_lines, con = file.path(table_dir, "ratio_method_summary.txt"))

metric_order <- c(direct_markers, ratio_metrics)
suppressive_color <- "#D55E00"
effector_color <- "#0072B2"
plot_colors <- ifelse(metric_order %in% c(suppressive_markers, ratio_metrics), suppressive_color, effector_color)

plot_panel <- function(results_df, main_title) {
  for (score_name in c("hh_activation_score", "hh_ligand_score")) {
    plot_df <- results_df[results_df$score_type == score_name, ]
    plot_df <- plot_df[match(metric_order, plot_df$metric), ]
    bp <- barplot(
      plot_df$rho,
      names.arg = plot_df$metric,
      las = 2,
      col = plot_colors,
      ylim = c(-0.7, 0.7),
      ylab = "Spearman rho",
      main = ifelse(score_name == "hh_activation_score", paste(main_title, "Activation"), paste(main_title, "Ligand"))
    )
    abline(h = 0, lty = 2, col = "gray40")
    sig_labels <- ifelse(plot_df$fdr < 0.001, "***", ifelse(plot_df$fdr < 0.01, "**", ifelse(plot_df$fdr < 0.05, "*", "")))
    text(bp, plot_df$rho, labels = sig_labels, pos = ifelse(plot_df$rho >= 0, 3, 1), cex = 0.9)
    legend(
      "topleft",
      legend = c("Suppressive / checkpoint", "Effector / inflammatory"),
      fill = c(suppressive_color, effector_color),
      bty = "n",
      cex = 0.75
    )
  }
}

png(
  filename = file.path(figure_dir, "ratio_method_comparison_figure.png"),
  width = 1800,
  height = 1600,
  res = 180
)
par(mfrow = c(2, 2), mar = c(10, 5, 4, 2))
plot_panel(results_old, "Old Ratio")
plot_panel(results_corrected, "Corrected Ratio")
dev.off()

png(
  filename = file.path(figure_dir, "corrected_ratio_activation_vs_ligand_correlations.png"),
  width = 1800,
  height = 900,
  res = 180
)
par(mfrow = c(1, 2), mar = c(10, 5, 4, 2))
plot_panel(results_corrected, "Corrected Ratio")
dev.off()

pdf(
  file = file.path(figure_dir, "corrected_ratio_activation_vs_ligand_correlations.pdf"),
  width = 12,
  height = 6
)
par(mfrow = c(1, 2), mar = c(10, 5, 4, 2))
plot_panel(results_corrected, "Corrected Ratio")
dev.off()

get_corrected_fdr <- function(metric_name) {
  hit <- results_corrected[
    results_corrected$score_type == "hh_activation_score" &
      results_corrected$metric == metric_name,
  ]
  if (nrow(hit) == 0) {
    return(NA_real_)
  }
  as.numeric(hit$fdr[[1]])
}

plot_scatter_with_lm <- function(x, y, xlab, ylab, main_title, metric_name, point_color = "#0072B280") {
  keep <- complete.cases(x, y)
  x <- x[keep]
  y <- y[keep]
  stats <- compute_spearman(x, y)
  fdr_value <- get_corrected_fdr(metric_name)

  plot(
    x,
    y,
    pch = 16,
    col = point_color,
    xlab = xlab,
    ylab = ylab,
    main = main_title
  )
  abline(lm(y ~ x), col = "#333333", lwd = 2)
  legend(
    "topleft",
    legend = c(
      sprintf("Spearman rho = %.3f", as.numeric(stats[["rho"]])),
      sprintf("FDR = %.3g", fdr_value)
    ),
    bty = "n",
    cex = 0.8
  )
}

png(
  filename = file.path(figure_dir, "hh_activation_selected_scatterplots.png"),
  width = 1800,
  height = 600,
  res = 180
)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
plot_scatter_with_lm(
  transformed$hh_activation_score,
  transformed$arg1_cxcl10_ratio_corrected,
  "HH Activation Score",
  "ARG1/CXCL10 ratio",
  "HH Activation vs ARG1/CXCL10",
  "arg1_cxcl10_ratio"
)
plot_scatter_with_lm(
  transformed$hh_activation_score,
  transformed$mrc1_cd8a_ratio_corrected,
  "HH Activation Score",
  "MRC1/CD8A ratio",
  "HH Activation vs MRC1/CD8A",
  "mrc1_cd8a_ratio"
)
plot_scatter_with_lm(
  transformed$hh_activation_score,
  transformed$FOXP3_log,
  "HH Activation Score",
  "FOXP3 expression",
  "HH Activation vs FOXP3",
  "FOXP3"
)
dev.off()

pdf(
  file = file.path(figure_dir, "hh_activation_selected_scatterplots.pdf"),
  width = 14,
  height = 5
)
par(mfrow = c(1, 3), mar = c(5, 5, 4, 1))
plot_scatter_with_lm(
  transformed$hh_activation_score,
  transformed$arg1_cxcl10_ratio_corrected,
  "HH Activation Score",
  "ARG1/CXCL10 ratio",
  "HH Activation vs ARG1/CXCL10",
  "arg1_cxcl10_ratio"
)
plot_scatter_with_lm(
  transformed$hh_activation_score,
  transformed$mrc1_cd8a_ratio_corrected,
  "HH Activation Score",
  "MRC1/CD8A ratio",
  "HH Activation vs MRC1/CD8A",
  "mrc1_cd8a_ratio"
)
plot_scatter_with_lm(
  transformed$hh_activation_score,
  transformed$FOXP3_log,
  "HH Activation Score",
  "FOXP3 expression",
  "HH Activation vs FOXP3",
  "FOXP3"
)
dev.off()

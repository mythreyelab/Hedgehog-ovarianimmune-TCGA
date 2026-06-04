file_arg <- grep("^--file=", commandArgs(), value = TRUE)
script_path <- if (length(file_arg) == 1) {
  normalizePath(gsub("~\\+~", " ", sub("^--file=", "", file_arg)), winslash = "/", mustWork = FALSE)
} else {
  normalizePath("scripts/07_make_hedgehog_immune_heatmap.R", winslash = "/", mustWork = FALSE)
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

hedgehog_genes <- c("GLI1", "GLI2", "PTCH1", "PTCH2", "HHIP", "IHH", "SHH")
immune_genes <- c("ARG1", "MRC1", "CD80", "FOXP3", "CD8A", "CD274", "PDCD1", "CXCL10", "CCL5", "CCL17", "CCL22")

figure_dir <- file.path(output_base_dir, "figures")
table_dir <- file.path(output_base_dir, "tables")
dir.create(figure_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(table_dir, showWarnings = FALSE, recursive = TRUE)

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

sample_cols <- setdiff(colnames(expr_raw), c("Hugo_Symbol", "Entrez_Gene_Id"))
expr_matrix <- as.matrix(expr_raw[, sample_cols, drop = FALSE])
rownames(expr_matrix) <- expr_raw$Hugo_Symbol
storage.mode(expr_matrix) <- "numeric"

missing_hedgehog <- setdiff(hedgehog_genes, rownames(expr_matrix))
missing_immune <- setdiff(immune_genes, rownames(expr_matrix))

hedgehog_genes <- intersect(hedgehog_genes, rownames(expr_matrix))
immune_genes <- intersect(immune_genes, rownames(expr_matrix))

if (length(hedgehog_genes) == 0 || length(immune_genes) == 0) {
  stop("No overlapping Hedgehog or immune genes were found in the expression matrix.")
}

compute_pair <- function(x, y) {
  keep <- complete.cases(x, y)
  n <- sum(keep)
  if (n < 3) {
    return(c(rho = NA_real_, p_value = NA_real_, n = n))
  }
  test <- suppressWarnings(cor.test(x[keep], y[keep], method = "spearman", exact = FALSE))
  c(rho = unname(test$estimate), p_value = test$p.value, n = n)
}

results <- data.frame(
  hedgehog_gene = character(),
  immune_gene = character(),
  rho = numeric(),
  p_value = numeric(),
  fdr = numeric(),
  n = integer(),
  stringsAsFactors = FALSE
)

for (hh in hedgehog_genes) {
  for (imm in immune_genes) {
    stats <- compute_pair(expr_matrix[hh, ], expr_matrix[imm, ])
    results <- rbind(
      results,
      data.frame(
        hedgehog_gene = hh,
        immune_gene = imm,
        rho = as.numeric(stats[["rho"]]),
        p_value = as.numeric(stats[["p_value"]]),
        fdr = NA_real_,
        n = as.integer(stats[["n"]]),
        stringsAsFactors = FALSE
      )
    )
  }
}

results$fdr <- p.adjust(results$p_value, method = "BH")
results <- results[order(results$fdr, -abs(results$rho), results$hedgehog_gene, results$immune_gene), ]

write.csv(
  results,
  file = file.path(table_dir, "hedgehog_immune_spearman_correlations.csv"),
  row.names = FALSE,
  quote = TRUE
)

rho_matrix <- xtabs(rho ~ hedgehog_gene + immune_gene, data = results)
fdr_matrix <- xtabs(fdr ~ hedgehog_gene + immune_gene, data = results)

write.csv(
  as.data.frame.matrix(rho_matrix),
  file = file.path(table_dir, "hedgehog_immune_rho_matrix.csv"),
  quote = TRUE
)

write.csv(
  as.data.frame.matrix(fdr_matrix),
  file = file.path(table_dir, "hedgehog_immune_fdr_matrix.csv"),
  quote = TRUE
)

sig_count <- sum(!is.na(results$fdr) & results$fdr < 0.05)
top_hits <- head(results[!is.na(results$fdr), ], 15)

summary_lines <- c(
  "Hedgehog vs immune gene correlation analysis",
  "",
  paste("Input file:", normalizePath(input_file, winslash = "/", mustWork = FALSE)),
  paste("Samples:", length(sample_cols)),
  paste("Hedgehog genes analyzed:", paste(hedgehog_genes, collapse = ", ")),
  paste("Immune genes analyzed:", paste(immune_genes, collapse = ", ")),
  paste("Missing Hedgehog genes from requested panel:", if (length(missing_hedgehog)) paste(missing_hedgehog, collapse = ", ") else "None"),
  paste("Missing immune genes from requested panel:", if (length(missing_immune)) paste(missing_immune, collapse = ", ") else "None"),
  paste("Total gene pairs tested:", nrow(results)),
  paste("Significant pairs at FDR < 0.05:", sig_count),
  "",
  "Top associations by FDR:"
)

if (nrow(top_hits) == 0) {
  summary_lines <- c(summary_lines, "No evaluable gene pairs.")
} else {
  top_text <- apply(top_hits, 1, function(row) {
    sprintf(
      "%s vs %s: rho=%.3f, p=%.3g, FDR=%.3g, n=%s",
      row[["hedgehog_gene"]],
      row[["immune_gene"]],
      as.numeric(row[["rho"]]),
      as.numeric(row[["p_value"]]),
      as.numeric(row[["fdr"]]),
      row[["n"]]
    )
  })
  summary_lines <- c(summary_lines, top_text)
}

writeLines(summary_lines, con = file.path(table_dir, "hedgehog_immune_heatmap_summary.txt"))

plot_matrix <- rho_matrix[hedgehog_genes, immune_genes, drop = FALSE]
plot_fdr <- fdr_matrix[hedgehog_genes, immune_genes, drop = FALSE]
label_matrix <- matrix("", nrow = nrow(plot_fdr), ncol = ncol(plot_fdr))
label_matrix[!is.na(plot_fdr) & plot_fdr < 0.05] <- "*"
label_matrix[!is.na(plot_fdr) & plot_fdr < 0.01] <- "**"
label_matrix[!is.na(plot_fdr) & plot_fdr < 0.001] <- "***"

png(
  filename = file.path(figure_dir, "hedgehog_immune_correlation_heatmap.png"),
  width = 1700,
  height = 950,
  res = 180
)

par(mar = c(12, 8, 4, 2))
image(
  x = seq_len(ncol(plot_matrix)),
  y = seq_len(nrow(plot_matrix)),
  z = t(plot_matrix[nrow(plot_matrix):1, , drop = FALSE]),
  col = colorRampPalette(c("#2166AC", "#F7F7F7", "#B2182B"))(200),
  zlim = c(-1, 1),
  xaxt = "n",
  yaxt = "n",
  xlab = "",
  ylab = "",
  main = "Spearman Correlation: Hedgehog vs Immune Genes"
)

axis(1, at = seq_len(ncol(plot_matrix)), labels = colnames(plot_matrix), las = 2, cex.axis = 0.9)
axis(2, at = seq_len(nrow(plot_matrix)), labels = rev(rownames(plot_matrix)), las = 2, cex.axis = 0.9)

for (i in seq_len(nrow(plot_matrix))) {
  for (j in seq_len(ncol(plot_matrix))) {
    text(
      x = j,
      y = nrow(plot_matrix) - i + 1,
      labels = label_matrix[i, j],
      cex = 0.9
    )
  }
}

dev.off()

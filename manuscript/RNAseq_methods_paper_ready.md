# RNA-seq Methods

RNA sequencing reads were aligned to the reference genome using STAR (`vX.X`, pending verification from the original alignment workflow), and gene-level count matrices were generated using featureCounts (`vX.X`, pending verification). The modular workflow in this folder expects a featureCounts-derived raw count matrix as input.

The locally confirmed sample structure comprised 16 RNA-seq libraries representing four experimental groups: `PCEFL_2D` (`n = 4`), `PCEFL_Suspension_24hr` (`n = 4`), `MAPK15OE_2D` (`n = 4`), and `MAPK15OE_Suspension_24hr` (`n = 4`). In the consolidated workflow, sample metadata are recorded in `data/metadata/sample_metadata_confirmed_local.csv`.

Differential expression analysis is organized around a template design formula of `~ batch + condition`, where `condition` encodes the experimental group and `batch` can be included when applicable. Because the locally archived sample sheets also support a two-factor design (`mapk15_status` and `attachment_state`), the underlying experiment may alternatively be modeled as `~ mapk15_status * attachment_state` when interaction testing is required. The current workflow documents both formulations so that the manuscript can align with the final analytical choice.

For differential expression, genes are intended to be tested with DESeq2 using the Benjamini-Hochberg procedure for multiple-testing correction. Genes meeting an adjusted `p < 0.05` and an absolute `log2` fold change greater than 1 are designated as differentially expressed in the current workflow specification.

Hallmark pathway analysis is implemented downstream using ranked differential expression results and Hallmark gene sets. The modular workflow uses `msigdbr` to retrieve Hallmark gene sets and `fgsea` to calculate enrichment statistics. Volcano plots are generated from DESeq2 result tables using the same adjusted `p`-value and absolute `log2` fold-change thresholds applied in the differential-expression summaries.

## Accuracy note

The sample grouping described above is confirmed locally from archived sample sheets and result folders. However, the original STAR command, featureCounts command, and exact upstream DESeq2 execution script were not present in the current workspace at the time this folder was assembled. Those fields therefore remain explicitly marked as verification items and should be updated before final manuscript submission if the original pipeline records become available.

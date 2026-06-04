# Raw Count Matrix

Place the gene-level raw count matrix generated from `featureCounts` in this directory.

Expected file name:

- `featurecounts_gene_counts.tsv`

Expected columns:

- annotation columns such as `gene_id` and `gene_name`
- one column per RNA-seq library matching the `sample_id` values in `data/metadata/sample_metadata_confirmed_local.csv`

This file was not present in the local workspace at the time the modular workflow was created, so the DESeq2, GSEA, and volcano scripts are scaffolded but were not executed end-to-end here.

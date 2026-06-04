# Hedgehog TCGA and RNA-seq Analysis Repo

## Purpose

This repository packages the final paper-facing Hedgehog and immune-correlation figures together with a cleaned RNA-seq workflow scaffold. It is structured so the included final figures remain traceable to code, while larger upstream inputs that are not ready for public distribution stay out of version control.

## Included

- final paper figures in `results/figures/`
- figure-generating R scripts in `scripts/07_make_hedgehog_immune_heatmap.R` and `scripts/08_make_ratio_method_figures.R`
- modular RNA-seq workflow scripts in `scripts/00_bootstrap_renv.R` through `scripts/06_run_all.R`
- configuration and metadata templates in `config/` and `data/metadata/`
- figure traceability manifests in `results/manifests/`
- supporting result tables used for the final figure set in `results/tables/`

## Not Included

- the TCGA expression matrix `data_mrna_seq_v2_rsem.txt`
- the RNA-seq raw `featureCounts` count matrix
- large intermediate outputs from DESeq2 or GSEA
- any protected or non-public source data

## Repository Layout

- `config/`: workflow configuration
- `data/metadata/`: sample metadata included with the workflow
- `data/raw_counts/`: placeholder for raw RNA-seq count matrix
- `data/tcga_expression/`: placeholder for the TCGA expression matrix used by the figure scripts
- `renv/`: renv bootstrap notes
- `results/figures/`: final figure files included in the paper package
- `results/manifests/`: code-to-figure traceability manifests
- `results/tables/`: supporting correlation and scoring tables
- `scripts/`: numbered workflow and figure-generation R scripts

## Reproducing the Included Figure Set

1. Place `data_mrna_seq_v2_rsem.txt` into `data/tcga_expression/`.
2. Run `Rscript scripts/07_make_hedgehog_immune_heatmap.R`.
3. Run `Rscript scripts/08_make_ratio_method_figures.R`.

These scripts write figures to `results/figures/` and supporting tables to `results/tables/`.

## Running the RNA-seq Workflow

1. Review `config/analysis_config.yml`.
2. Place `featurecounts_gene_counts.tsv` into `data/raw_counts/`.
3. Review `data/metadata/sample_metadata_confirmed_local.csv`.
4. Run `Rscript scripts/00_bootstrap_renv.R`.
5. Run `Rscript scripts/06_run_all.R`.

## Notes on Confirmed Versus Placeholder Metadata

The local workspace supports the downstream analysis structure and confirms the 16-sample MAPK15 experiment metadata. However, the original `STAR` version, `featureCounts` version, and some upstream alignment details were not recoverable from local scripts alone and remain marked as placeholders in the config and methods files.

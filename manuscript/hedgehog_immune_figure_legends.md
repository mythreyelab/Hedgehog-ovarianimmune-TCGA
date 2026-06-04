# Final Figure Legends

## Figure 1. Hedgehog-Immune Correlation Heatmap

### Description
Across TCGA ovarian tumors, hedgehog pathway genes showed selective rather than uniform relationships with immune-associated transcripts. Correlation patterns differed across individual pathway components, but a recurring theme was that hedgehog signaling aligned more strongly with specific immune regulatory programs than with a broad, nonspecific immune transcriptional state. Together, these data support the idea that hedgehog pathway activity is linked to structured immune remodeling within the tumor microenvironment.

### Legend
**Figure 1. Spearman correlation heatmap of hedgehog pathway genes and immune-related genes in TCGA ovarian tumors.** Heatmap showing Spearman rank correlation coefficients between hedgehog pathway genes (`GLI1`, `GLI2`, `PTCH1`, `PTCH2`, `HHIP`, `IHH`, and `SHH`) and selected immune-related genes across TCGA ovarian tumor samples. Gene expression values were derived from normalized mRNA expression data and analyzed after log2 transformation (`log2[RSEM + 1]`). Color intensity indicates the direction and magnitude of the correlation, with warm colors representing positive correlations and cool colors representing negative correlations. Statistical significance is indicated within the heatmap (`*P < 0.05`, `**P < 0.01`, `***P < 0.001`).

## Figure 2. Corrected Ratio Activation vs Ligand Correlations

### Description
Using the floor-aware corrected ratio approach, both hedgehog activation and hedgehog ligand scores were associated with a pattern consistent with altered immune balance rather than uniform association with all immune genes. Higher hedgehog signaling tracked with stronger suppressive-to-effector skewing when immune features were summarized as composite ratios such as `ARG1/CXCL10` and `MRC1/CD8A`, while several effector and inflammatory markers showed inverse relationships. These findings support a model in which hedgehog pathway activity is associated with a more suppressive immune expression context in TCGA ovarian tumors.

### Legend
**Figure 2. Correlation of corrected hedgehog activation and ligand scores with immune-related genes and ratio features.** Bar plots showing Spearman correlation coefficients between the corrected hedgehog activation score or corrected hedgehog ligand score and selected immune-related features across TCGA ovarian tumor samples. Features include individual suppressive/checkpoint-associated genes, effector/inflammatory genes, and the composite ratios `ARG1/CXCL10` and `MRC1/CD8A`. Orange bars indicate suppressive/checkpoint-associated features, and blue bars indicate effector/inflammatory features. Ratio features were calculated using the floor-aware corrected approach, defined as `log2((floored numerator + 1) / (floored denominator + 1))`, where higher values indicate relatively greater `ARG1` expression compared with `CXCL10` or relatively greater `MRC1` expression compared with `CD8A`. Asterisks denote FDR-adjusted significance (`*FDR < 0.05`, `**FDR < 0.01`, `***FDR < 0.001`).

## Figure 3. HH Activation Selected Scatterplots

### Description
The selected scatterplots provide a sample-level view of the immune-skewing relationships associated with hedgehog activation. Higher hedgehog activation was associated with increased values of the `ARG1/CXCL10` and `MRC1/CD8A` ratios, consistent with relative enrichment of suppressive myeloid features over effector or inflammatory features. In contrast, the relationship with `FOXP3` indicates that not all immunoregulatory genes follow the same pattern, emphasizing that hedgehog signaling is linked to a specific configuration of immune remodeling rather than a uniform increase in all suppressive markers.

### Legend
**Figure 3. Selected scatterplots of hedgehog activation versus immune ratio and immune gene features.** Scatterplots showing the relationship between the hedgehog activation score and three selected immune features across TCGA ovarian tumor samples. The left panel shows the corrected `ARG1/CXCL10` ratio, the middle panel shows the corrected `MRC1/CD8A` ratio, and the right panel shows `FOXP3` expression. Corrected ratio features were calculated using the floor-aware method, defined as `log2((floored numerator + 1) / (floored denominator + 1))`, such that higher ratio values indicate relatively greater `ARG1` over `CXCL10` or `MRC1` over `CD8A`. Each point represents one tumor sample, the solid line indicates the fitted linear trend, and statistical association was assessed using Spearman correlation.

# Final Figure Methods

## Dataset and expression matrix

Gene expression analyses were performed using the TCGA ovarian cancer mRNA expression matrix (`data_mrna_seq_v2_rsem.txt`). The input table contained Hugo gene symbols, Entrez gene identifiers, and tumor sample-level RSEM-normalized expression values. Rows with missing gene symbols were excluded, and duplicate Hugo symbols were removed by retaining the first occurrence in the analysis scripts used for the final figures. All analyses were performed across 300 tumor samples.

## Figure 1. Hedgehog-immune correlation heatmap

### Methods
For the hedgehog-immune correlation heatmap, a predefined hedgehog gene panel (`GLI1`, `GLI2`, `PTCH1`, `PTCH2`, `HHIP`, `IHH`, and `SHH`) was compared with a predefined immune-associated gene panel (`ARG1`, `MRC1`, `CD80`, `FOXP3`, `CD8A`, `CD274`, `PDCD1`, `CXCL10`, `CCL5`, `CCL17`, and `CCL22`). Immune and chemokine markers for the human validation analysis were selected on the basis of established roles as immunosuppressive or immune-regulatory mediators in human high-grade serous ovarian cancer, together with chemokine features highlighted by spatial transcriptomic inflammatory-hallmark analyses in the murine model. Pairwise associations between each hedgehog gene and each immune gene were assessed using Spearman rank correlation across all tumor samples. P values from all tested gene pairs were adjusted for multiple comparisons using the Benjamini-Hochberg false discovery rate (FDR) method. Correlation coefficients and adjusted significance values were exported as matrices and visualized as a heatmap. Heatmap cells were colored according to Spearman rho values, and significance markers were overlaid using thresholds of `P < 0.05`, `P < 0.01`, and `P < 0.001` in the plotting workflow.

## Figure 2. Corrected ratio activation vs ligand correlations

### Methods
To evaluate pathway-level associations with immune features, hedgehog activation and ligand scores were generated at the sample level. The hedgehog activation score was defined as the mean z score of log-transformed expression values for `GLI1`, `GLI2`, `PTCH1`, `PTCH2`, and `HHIP`. The hedgehog ligand score was defined as the mean z score of log-transformed expression values for `IHH` and `SHH`. For direct immune markers, expression values below 1 were first floored to 1 and then transformed as `log2(x + 1)`.

Two composite immune-balance metrics were evaluated: `ARG1/CXCL10` and `MRC1/CD8A`. In the final corrected approach, these ratios were calculated from floored raw expression values using `log2((floored numerator + 1) / (floored denominator + 1))`. This approach was selected to reduce distortion from very low denominator values while preserving interpretable suppressive-to-effector contrasts. The analysis also retained a comparison to the earlier ratio derivation, defined as `log2((log-transformed numerator + 1) / (log-transformed denominator + 1))`, in order to quantify the impact of ratio construction on downstream correlations.

Spearman correlations were then calculated between each hedgehog score and a panel of immune-related features, including direct immune markers (`ARG1`, `MRC1`, `CD80`, `FOXP3`, `CD8A`, `CD274`, `PDCD1`, `CCL22`, `CCL17`, `CCL5`, and `CXCL10`) as well as the two ratio metrics. P values were adjusted within each score-specific analysis using the Benjamini-Hochberg FDR method. The corrected ratio figure displays the final floor-aware results for the activation and ligand scores.

## Figure 3. HH activation selected scatterplots

### Methods
Selected scatterplots were generated to visualize sample-level relationships between the hedgehog activation score and representative immune features derived from the corrected ratio workflow. The plotted outcomes were the corrected `ARG1/CXCL10` ratio, the corrected `MRC1/CD8A` ratio, and `FOXP3` expression after flooring values below 1 and applying `log2(x + 1)`. Each scatterplot displays one point per tumor sample and includes a fitted linear trend line for visualization. Statistical associations shown in the analysis output were based on Spearman rank correlation.

## Statistical analysis and software

All analyses were performed in R. Nonparametric associations were assessed using `cor.test(..., method = "spearman", exact = FALSE)`, and multiple-testing correction was performed using the Benjamini-Hochberg procedure (`p.adjust`, method = `"BH"`). Figures were generated using base R graphics and exported as PNG and PDF files.

## Note for manuscript consistency

The current heatmap rerun script computes hedgehog-immune correlations directly from the processed expression matrix after filtering missing and duplicate gene symbols. If the manuscript text or legend is intended to state that these correlations were calculated from additional log2-transformed values, that wording should be checked against the exact script used to generate the final heatmap so the methods and legends match precisely.

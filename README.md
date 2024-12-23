## scRNA-seq and bulk RNA-seq analysis of pseudomyxoma peritonei (PMP)

### Abstract
Backgrounds: Psuedomyxoma peritonei (PMP) is a rare colon disease, whose symptoms and treatments have been numerously investigated, However, omics profiling of this disease remains significantly underexplored. Here, we present single-cell transcriptomic profiling of five PMP cases and bulk RNA-Seq profiling of 19 fresh frozen tissues and 25 FFPE samples to identify gene characteristics specific to each cell type associated with PMP pathogenesis.

Experimental Design: We conducted single-cell transcriptomic profiling on five PMP cases and three normal peritoneum cases to identify cell type-specific gene signatures associated with PMP. Validation was performed using bulk RNA-seq datasets from two independent cohorts: 19 fresh frozen tissue samples (12 PMPs) and 34 formalin-fixed paraffin-embedded samples (25 PMPs). 

Results: Single-cell transcriptomic analysis revealed the cellular diversity of PMP, contrasting the coexistence of epithelial and mesenchymal characteristics within PMP cells against the primarily mesenchymal composition of normal peritoneum. These findings were confirmed via bulk RNA-seq by highlighting increased expressions of genes associated with cellular functions such as proliferation, movement, and immune responses in PMP tissues. 

Conclusion: This study investigates the intricate PMP molecular landscape, demonstrating the coexistence of distinct epithelial and mesenchymal cells. 

### Descriptions of files in the repository
R Code

<b>Bulk RNA-seq differential expression gene (DEG) analysis</b> : The R code for identifying DEGs in bulk RNA-seq requires the data from GSE228376 (19 fresh frozen samples) and GSE228375 (25 FFPE samples).

<b>PMP_scRNA_analysis</b> : The R code for scRNA-seq analysis requires data from GSE228377 (5 PMP single-cell samples) and GSE130888 (GSM3755693, GSM3755694, GSM3755695; normal peritoneum samples).

### How to use
To use the code provided in this repository, please first install R and run RStudio.

The analyses were conducted using R version 4.1.1 and implemented with Seurat version 4.3.

The rds files for the five PMP samples are registered under GSE228377.

The log2CPM of the bulk RNA-seq data can be referenced in the log2CPM.xlsx file provided within the GSE228375 and GSE228376 datasets.

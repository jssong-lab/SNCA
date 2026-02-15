* Overview

The Seurat pipelines for processing scRNA-seq data from 10x genomics
are included in the SNCA_Seurat_pipeline.R and IGH_Seurat_pipeline.R
files for SNCA-overexpression in mouse and IGH analysis in human,
respectively.  Each of these generates all UMAP plots, gene expression
scatter and violin plots, and dotplots. Heat maps and gene ontology
(GO) related figures are produced separately. Furthermore, the Seurat
pipelines perform differential expression analysis which is used
downstream.


Additional analysis scripts for SNCA include:

  - 'SNCA_GO_category_scatter_plotter.py' produces scatter plots of GO
    terms (Extended Figure 7H).  It takes as input a csv file of
    labeled GO terms and a text file containing all adjusted p-values
    and enrichment scores used for unlabeled points (background).  All
    data was obtained from the ouptut of The Database for Annotation,
    Visualization, and Integrated Discovery (DAVID) using the methods
    described in Methods.
  - 'SNCA_heatmap.R' produces the clustered heatmaps in Extended
    Figure 8. It requires the Seurat object processed by the SNCA
    pipeline as well as gene lists for each gene ontology category,
    and a corresponding table of adjusted p-values for those
    genes. P-values are obtained from the differential expression (DE)
    analysis files output by SNCA pipeline.


Additional analysis scripts for IGH include:

  - 'IGH_metap_combiner.py' performs p-value combination and
    adjustment for all patients.  It takes as input the DE analysis
    output from the IGH pipeline and constructs a single table with
    row-wise genes and column-wise patient log2-FCs, p-values, and
    their respective combined values.

  - 'IGH_heatmap.R' which produced extended figure 9E.  It takes as
    input all gene expressions for the desired GO categories (see
    Methods) for all patients, which is obtained from the IGH
    pipeline.

  - 'IGH_neurodegeneration_analysis.py' which produces Extended Data
    Figure 9F.  It requires a folder of GO categories to be used in
    the heatmap and the data structure (constructed in
    'IGH_metap_combiner.py') of all gene expressions and log2-FC
    (lesional vs. perilesional) across all patients.




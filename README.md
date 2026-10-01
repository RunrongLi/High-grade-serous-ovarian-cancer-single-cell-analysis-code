# High-grade-serous-ovarian-cancer-single-cell-analysis-code.
This repository contains the single-cell RNA-seq analysis code for the study of senescent cancer-associated fibroblasts (CAFs) in high-grade serous ovarian cancer (HGSOC).

## Analysis example

The scripts in this repository use **GSE165897** as an example dataset to demonstrate the complete single-cell analysis workflow, including quality control, integration, clustering, and downstream analysis.

Other datasets listed below are provided for data availability and can be analyzed by adapting the input paths and sample metadata.

## Data availability

The scRNA-seq datasets used in this study are publicly available.

| Dataset | Description | Source |
|---|---|---|
| GSE165897 | Human HGSOC scRNA-seq | [GEO GSE165897](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE165897) |
| GSE144735 | Colorectal cancer (CRC) scRNA-seq | [GEO GSE144735](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE144735) |
| Mendeley dataset | Additional HGSOC dataset | [Mendeley Data rc47y6m9mp](https://data.mendeley.com/datasets/rc47y6m9mp/1) |


## Acknowledgement
Seurat was used for single-cell analysis. Please cite:
Hao Y, et al. Integrated analysis of multimodal single-cell data. Cell. 2021;184(13):3573-87.e29.

The `StackedVlnPlot` function in `03_gene_visualization.R` is adapted from Ming Tang's blog post:
https://divingintogeneticsandgenomics.rbind.io/post/stacked-violin-plot-for-visualizing-single-cell-data-in-seurat/

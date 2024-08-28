# YMC Multi-omics Analysis

## Introduction

This repository contains the core code of our analysis. It includes:
1. Singular value decomposition (SVD) of each omics data matrix. 
2. The polarization sample- and molecule-eigenvectors by sorting their loadings.
3. The enrichment analysis on the polarized gene-eigenvectors of the transcriptome.
4. The alignment of oxygen curves by DDTW.
5. The computation of time differences between different omics datasets.
6. Other results.


### Code language

All code in this repository is based on `R`.

### Usage

```bash
cd code
bash run.sh
```
All computational results can be found in the `results` directory generated.

## Directory Description


- ### `code` directory

* **`run.sh`**  
This script is used for running other scripts all in sequence.

* **`config.yaml`**  
This is the configuration file that contains basic information about datasets and script parameters.

* **`data_process.R`**
This script reads the data matrices of the transcriptome, metabolome, and epigenome, and prepares them for the following analysis.

* **`matrix_decomposition.R` & `matrix_decomposition_metabolic.R`**  
This script performs singular value decomposition(SVD) on data matrices, and obtains singular values and the polarized molecule- and sample-eigenvectors. These are saved in `results` directory.

* **`gene_enrichment_analysis.R` & `format_enrichment_result.R`**  
This script performs the gene set enrichment analysis by the Wilcoxon scoring method on the polarized gene eigenvector. The enrichment result is saved in `results` directory.

* **`enrichment_utils.R`**  
This script contains functions used to perform enrichment analysis.




### data folder

`pathway` directory contains the pathway data of [GO](https://geneontology.org/) (GO.BP, GO.CC, and GO.MF), [KEGG](https://www.kegg.jp/), and [Reactome](https://reactome.org/).

**sweeden_total_result_3.csv** and **reads_polished_reg_dp_fixed.csv** are the gene expression profile of human and rat. The rat RNA-seq sequencing reads of rats were downloaded from GEO under accession number GSE81811 and normalized by MUREN. The human Affymetrix microarray data were downloaded from GEO under accession number GSE38642 and normalized by sub-sub normalization.

**HuGene-1_0-st-v1.na36.hg19.probeset-single-name.csv** is the annotation file of the human microarray.

**ensembl_convertID_102.tsv** is the annotation file of the conversion among Ensembl id, Entrez id, and gene symbol.

### figures folder

This folder contains most of the figures in our analysis. All figures in the folder are plotted by `R`.

